#!/usr/bin/env python3
"""Rigorous Arb certificate for the four-cycle candidate.

The candidate is the continuous piecewise-linear Lienard system

    xdot = F(x) - y,   ydot = x,

with the decimal parameters in ``four-cycle-001.json`` interpreted as exact
rational numbers.  In a linearity zone put ``p = F(x) - y``.  Then

    xdot = p,   pdot = m p - x,

so every flight has an exact elementary flow.  For each proposed cycle this
script applies a Krawczyk test to the section amplitude and all flight times.
It then proves the prescribed branch by directed Arb ball arithmetic on a
normalized-time subdivision, proves every switch transverse, and proves that
the divergence integral is nonzero.

The first cycle uses a five-variable chart (section amplitude plus four
flights); the other cycles use nine-variable charts (amplitude plus eight
flights).  Floating-point or mpmath arithmetic is used only to construct a
preconditioner.  Every acceptance decision is made with Arb enclosures.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Iterable

import flint
from flint import arb, ctx, fmpq
import mpmath as mp


DEFAULT_PRECISION = 480
DEFAULT_RADIUS = "1e-25"
DEFAULT_SLICES = 512

# Numerically computed centers.  They are only seed data: the Krawczyk test
# below is what proves that an exact zero exists in each surrounding box.
# Each row is (section amplitude s, consecutive flight times).
CENTERS: tuple[tuple[str, ...], ...] = (
    (
        "18.543511823474473787728678940457730779622096964584535124563812372634893141430204",
        "3.3367132325803590379398256194866607050468155245524258188781025009593877175744979",
        "0.17701113250010402098101992534131276281696505545928125959987665779030489982631031",
        "3.2237638910571879065894351133295232964551437166265336715066464488677272501530714",
        "0.05299479879943667859017000963946931067964553826513493147821517418556808155399764",
    ),
    (
        "24.731127033873463213572344095983047275442411821170972961884280243029151792247137",
        "3.3367132325803590379398256194866607050468155245524258188781025009593877175744979",
        "0.13037775000695843464149406619350170984553674838662267639569639433903458286608847",
        "1.2945522370676825942193679295640784917557822730063796146054110461896785664797733",
        "0.046164026896590735175199037018695139543800555901082957313375478300747462092477258",
        "1.1692723099017226687196343022100509184313709769947552941655754332013531617608129",
        "0.038048226323509399306740453654646971876757392002660094713593387951656220048761535",
        "0.67798178811555500663661182070523464959016977651236976115634091292657300281925736",
        "0.039904225946609159012091519799819750062471090082322414375084631804341783325022028",
    ),
    (
        "35.794748987988437365913112052737345194179591365532136371389882621828829014989302",
        "3.3367132325803590379398256194866607050468155245524258188781025009593877175744979",
        "0.08870406191219561473263852620412598047481599122835384712210867785783336470212617",
        "0.86540806821641899117481094765619391035391673659454250886850631339072095327929332",
        "0.025231479582373257354555810780713787456015133242807806055376766413972435697527191",
        "1.9427157019496539907939774695115639745261957480129331427311155181467689720737195",
        "0.016048878232996533747504397553784113551309589801901386820800595009499030847076863",
        "0.39763275421520468064278525064791833490350253152095398278385237334981681797491396",
        "0.027680834231504168094814597838551762525955012735908002642625956855860370351797766",
    ),
    (
        "234.31929081370307740910995934351290012228737591390531658384819892691352360134609",
        "3.3367132325803590379398256194866607050468155245524258188781025009593877175744979",
        "0.013195059996341769037958827237999361267664738351676792710180986688488192081951357",
        "0.14841710359788900260856152121386534445411600157968354524199230243931216695217126",
        "0.0044716200491140515918946417511719837610042748258469738583529803303126010440182949",
        "3.1134499352939853219751167462775846353446754092852511704437214814618638259788336",
        "0.0016780636434162142581615173897086757040831021688250500756729335009716485061949525",
        "0.051750403910575525112525642203595639101909712787610597264392735702469834018521845",
        "0.0042615688089118733823348461578931335070856986417231045092084379025123045141171203",
    ),
)


@dataclass(frozen=True)
class ExactModel:
    break_text: tuple[str, str, str]
    slope_text: tuple[str, str, str, str]
    breaks_q: tuple[fmpq, fmpq, fmpq]
    slopes_q: tuple[fmpq, fmpq, fmpq, fmpq]
    intercepts_fraction: tuple[Fraction, Fraction, Fraction, Fraction]

    @property
    def breaks(self) -> list[arb]:
        return [arb(value) for value in self.breaks_q]

    @property
    def slopes(self) -> list[arb]:
        return [arb(value) for value in self.slopes_q]


@dataclass(frozen=True)
class Chart:
    slopes: list[arb]
    nodes: list[arb]
    turns: dict[int, int]
    itinerary: str


class CertificateFailure(RuntimeError):
    """Raised when a strict interval acceptance test cannot be proved."""


def exact_q(text: str) -> fmpq:
    """Parse a finite decimal as a rational without binary floating point."""

    value = Fraction(text)
    return fmpq(value.numerator, value.denominator)


def exact_fraction(text: str) -> Fraction:
    return Fraction(text)


def around(text: str, radius: str) -> arb:
    """Outward Arb enclosure of the exact rational interval c +/- r."""

    return arb(exact_q(text), exact_q(radius))


def ball_text(value: arb, digits: int = 36) -> str:
    return value.str(digits)


def minimum_ball(values: Iterable[arb]) -> arb:
    iterator = iter(values)
    try:
        result = next(iterator)
    except StopIteration as error:
        raise ValueError("minimum of an empty collection") from error
    for value in iterator:
        result = result.min(value)
    return result


def maximum_ball(values: Iterable[arb]) -> arb:
    iterator = iter(values)
    try:
        result = next(iterator)
    except StopIteration as error:
        raise ValueError("maximum of an empty collection") from error
    for value in iterator:
        result = result.max(value)
    return result


def require_positive(value: arb, label: str) -> arb:
    if not value > 0:
        raise CertificateFailure(
            f"could not prove {label} positive: {ball_text(value)}"
        )
    return value


def focus_flow(m, time, x, p):
    """Exact affine-focus flow in (x,p), valid for Arb or Dual values."""

    alpha = m / 2
    omega = (1 - alpha * alpha).sqrt()
    ratio = alpha / omega
    phase = omega * time
    cosine = phase.cos()
    sine = phase.sin()
    scale = (alpha * time).exp()
    x_out = scale * (
        (cosine - ratio * sine) * x + (sine / omega) * p
    )
    p_out = scale * (
        (-sine / omega) * x + (cosine + ratio * sine) * p
    )
    return x_out, p_out


class Dual:
    """First-order forward automatic differentiation over Arb balls."""

    def __init__(self, value: arb, derivative: list[arb]):
        self.value = value
        self.derivative = derivative

    @classmethod
    def constant(cls, value, dimension: int) -> "Dual":
        value = value if isinstance(value, arb) else arb(value)
        return cls(value, [arb(0) for _ in range(dimension)])

    def _coerce(self, other) -> "Dual":
        if isinstance(other, Dual):
            return other
        return Dual.constant(other, len(self.derivative))

    def __add__(self, other) -> "Dual":
        other = self._coerce(other)
        return Dual(
            self.value + other.value,
            [a + b for a, b in zip(self.derivative, other.derivative)],
        )

    __radd__ = __add__

    def __neg__(self) -> "Dual":
        return Dual(-self.value, [-value for value in self.derivative])

    def __sub__(self, other) -> "Dual":
        return self + (-self._coerce(other))

    def __rsub__(self, other) -> "Dual":
        return self._coerce(other) - self

    def __mul__(self, other) -> "Dual":
        other = self._coerce(other)
        return Dual(
            self.value * other.value,
            [
                left * other.value + self.value * right
                for left, right in zip(
                    self.derivative, other.derivative
                )
            ],
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "Dual":
        denominator = self.value * self.value
        return Dual(
            1 / self.value,
            [-value / denominator for value in self.derivative],
        )

    def __truediv__(self, other) -> "Dual":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other) -> "Dual":
        return self._coerce(other) * self.reciprocal()

    def sin(self) -> "Dual":
        return Dual(
            self.value.sin(),
            [self.value.cos() * value for value in self.derivative],
        )

    def cos(self) -> "Dual":
        return Dual(
            self.value.cos(),
            [-self.value.sin() * value for value in self.derivative],
        )

    def exp(self) -> "Dual":
        value = self.value.exp()
        return Dual(
            value,
            [value * derivative for derivative in self.derivative],
        )


def build_chart(short: bool, model: ExactModel) -> Chart:
    slopes = model.slopes
    breaks = model.breaks
    if short:
        return Chart(
            slopes=[slopes[0], slopes[0], slopes[1], slopes[0]],
            nodes=[arb(0), arb(0), breaks[0], breaks[0], arb(0)],
            turns={0: -1, 2: 1},
            itinerary="xi1+,xi1-",
        )
    return Chart(
        slopes=[
            slopes[0],
            slopes[0],
            slopes[1],
            slopes[2],
            slopes[3],
            slopes[2],
            slopes[1],
            slopes[0],
        ],
        nodes=[
            arb(0),
            arb(0),
            breaks[0],
            breaks[1],
            breaks[2],
            breaks[2],
            breaks[1],
            breaks[0],
            arb(0),
        ],
        turns={0: -1, 4: 1},
        itinerary="xi1+,xi2+,xi3+,xi3-,xi2-,xi1-",
    )


def residual_dual(
    short: bool, values: list[arb], model: ExactModel
) -> list[Dual]:
    chart = build_chart(short, model)
    dimension = len(values)
    variables: list[Dual] = []
    for index, value in enumerate(values):
        derivative = [arb(0) for _ in range(dimension)]
        derivative[index] = arb(1)
        variables.append(Dual(value, derivative))

    # On x=0, F(0)=0, hence p=-y=-s.
    p = -variables[0]
    residuals: list[Dual] = []
    for slope, time, x_in, x_target in zip(
        chart.slopes, variables[1:], chart.nodes[:-1], chart.nodes[1:]
    ):
        x_out, p = focus_flow(slope, time, x_in, p)
        residuals.append(x_out - x_target)
    # Returning to x=0 with p=-s closes the full state (x,y).
    residuals.append(p + variables[0])
    return residuals


def mp_chart(short: bool, model: ExactModel):
    slopes = [mp.mpf(text) for text in model.slope_text]
    breaks = [mp.mpf(text) for text in model.break_text]
    if short:
        return (
            [slopes[0], slopes[0], slopes[1], slopes[0]],
            [mp.mpf(0), mp.mpf(0), breaks[0], breaks[0], mp.mpf(0)],
        )
    return (
        [
            slopes[0],
            slopes[0],
            slopes[1],
            slopes[2],
            slopes[3],
            slopes[2],
            slopes[1],
            slopes[0],
        ],
        [
            mp.mpf(0),
            mp.mpf(0),
            breaks[0],
            breaks[1],
            breaks[2],
            breaks[2],
            breaks[1],
            breaks[0],
            mp.mpf(0),
        ],
    )


def mp_flow(m, time, x, p):
    alpha = m / 2
    omega = mp.sqrt(1 - alpha * alpha)
    ratio = alpha / omega
    phase = omega * time
    cosine = mp.cos(phase)
    sine = mp.sin(phase)
    scale = mp.exp(alpha * time)
    return (
        scale
        * ((cosine - ratio * sine) * x + (sine / omega) * p),
        scale
        * ((-sine / omega) * x + (cosine + ratio * sine) * p),
    )


def mp_residual(short: bool, values, model: ExactModel):
    slopes, nodes = mp_chart(short, model)
    p = -values[0]
    residuals = []
    for slope, time, x_in, x_target in zip(
        slopes, values[1:], nodes[:-1], nodes[1:]
    ):
        x_out, p = mp_flow(slope, time, x_in, p)
        residuals.append(x_out - x_target)
    residuals.append(p + values[0])
    return residuals


def arb_dyadic_from_mp(value) -> arb:
    """Return a fixed dyadic preconditioner entry near an mpmath value."""

    rational_enclosure = arb(exact_q(mp.nstr(value, 150)))
    return rational_enclosure.mid()


def matmul(left, right):
    rows = len(left)
    shared = len(right)
    columns = len(right[0])
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(shared)),
                arb(0),
            )
            for j in range(columns)
        ]
        for i in range(rows)
    ]


def matvec(matrix, vector):
    return [
        sum(
            (matrix[i][j] * vector[j] for j in range(len(vector))),
            arb(0),
        )
        for i in range(len(matrix))
    ]


def exact_box_data(center_text: str, radius_text: str):
    center = exact_fraction(center_text)
    radius = exact_fraction(radius_text)
    return center - radius, center + radius


def prove_box_positive(center_text: str, radius_text: str, label: str) -> None:
    lower, _upper = exact_box_data(center_text, radius_text)
    if lower <= 0:
        raise CertificateFailure(f"{label} box is not strictly positive")


def krawczyk_check(
    short: bool,
    centers: tuple[str, ...],
    radius: str,
    model: ExactModel,
) -> tuple[dict[str, object], list[arb]]:
    dimension = len(centers)
    expected_dimension = 5 if short else 9
    if dimension != expected_dimension:
        raise CertificateFailure(
            f"wrong chart dimension {dimension}, expected {expected_dimension}"
        )

    for index, text in enumerate(centers):
        name = "section amplitude" if index == 0 else f"time {index}"
        prove_box_positive(text, radius, name)

    mp_centers = [mp.mpf(text) for text in centers]
    jacobian = mp.matrix(dimension)
    for row in range(dimension):
        for column in range(dimension):

            def component(value, row=row, column=column):
                variables = mp_centers.copy()
                variables[column] = value
                return mp_residual(short, variables, model)[row]

            jacobian[row, column] = mp.diff(
                component, mp_centers[column]
            )
    inverse = jacobian**-1
    preconditioner = [
        [arb_dyadic_from_mp(inverse[i, j]) for j in range(dimension)]
        for i in range(dimension)
    ]

    full_box = [around(text, radius) for text in centers]
    center_points = [box.mid() for box in full_box]

    # Verify that each outward Arb box contains the intended exact rational
    # box.  The latter is used below for a human-readable strict inclusion.
    exact_bounds: list[tuple[arb, arb]] = []
    for box, text in zip(full_box, centers):
        lower_fraction, upper_fraction = exact_box_data(text, radius)
        lower = arb(
            fmpq(lower_fraction.numerator, lower_fraction.denominator)
        )
        upper = arb(
            fmpq(upper_fraction.numerator, upper_fraction.denominator)
        )
        if not box.contains(lower) or not box.contains(upper):
            raise CertificateFailure("Arb box did not enclose rational box")
        exact_bounds.append((lower, upper))

    point_residual = [
        item.value for item in residual_dual(short, center_points, model)
    ]
    interval_residual = residual_dual(short, full_box, model)
    interval_jacobian = [
        [interval_residual[i].derivative[j] for j in range(dimension)]
        for i in range(dimension)
    ]
    product = matmul(preconditioner, interval_jacobian)
    defect = [
        [
            (arb(1) if i == j else arb(0)) - product[i][j]
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    correction = matvec(preconditioner, point_residual)
    centered_box = [
        full_box[i] - center_points[i] for i in range(dimension)
    ]
    defect_box = matvec(defect, centered_box)
    krawczyk = [
        center_points[i] - correction[i] + defect_box[i]
        for i in range(dimension)
    ]

    inclusion_margins: list[arb] = []
    for enclosure, (lower, upper) in zip(krawczyk, exact_bounds):
        left_margin = enclosure - lower
        right_margin = upper - enclosure
        inclusion_margins.extend((left_margin, right_margin))
        require_positive(left_margin, "left Krawczyk inclusion margin")
        require_positive(right_margin, "right Krawczyk inclusion margin")

    row_sums = [
        sum((abs(entry).upper() for entry in row), arb(0))
        for row in defect
    ]
    if not all(row_sum < 1 for row_sum in row_sums):
        raise CertificateFailure("could not prove Krawczyk contraction")

    q_upper = maximum_ball(row_sums)
    max_residual = maximum_ball(
        abs(value).upper() for value in point_residual
    )
    min_inclusion = minimum_ball(inclusion_margins)
    return (
        {
            "dimension": dimension,
            "flight_count": dimension - 1,
            "strict_inclusion": True,
            "contraction_norm_strictly_below_one": True,
            "defect_infinity_norm_upper_enclosure": ball_text(q_upper),
            "point_residual_max_abs_enclosure": ball_text(max_residual),
            "strict_inclusion_margin_enclosure": ball_text(min_inclusion),
            "strict_inclusion_margin_lower_bound": ball_text(
                min_inclusion.lower()
            ),
        },
        full_box,
    )


def normalized_cell(index: int, slices: int) -> arb:
    return arb(
        fmpq(2 * index + 1, 2 * slices),
        fmpq(1, 2 * slices),
    )


def branch_and_multiplier_check(
    short: bool,
    centers: tuple[str, ...],
    radius: str,
    slices: int,
    model: ExactModel,
) -> dict[str, object]:
    chart = build_chart(short, model)
    variables = [around(text, radius) for text in centers]
    section = variables[0]
    times = variables[1:]
    p = -section
    endpoints: list[arb] = []
    velocity_guards: list[arb] = []
    position_guards: list[arb] = []
    segment_reports: list[dict[str, object]] = []

    for segment, (slope, total_time, x_in, x_target) in enumerate(
        zip(chart.slopes, times, chart.nodes[:-1], chart.nodes[1:])
    ):
        p_in = p
        segment_velocity: list[arb] = []
        segment_position: list[arb] = []
        if segment not in chart.turns:
            direction = 1 if x_target > x_in else -1
            for cell_index in range(slices):
                fraction = normalized_cell(cell_index, slices)
                _x_cell, p_cell = focus_flow(
                    slope, fraction * total_time, x_in, p_in
                )
                signed_velocity = direction * p_cell
                require_positive(
                    signed_velocity,
                    f"segment {segment} monotonicity guard",
                )
                segment_velocity.append(signed_velocity)
        else:
            side = chart.turns[segment]
            first_fraction = normalized_cell(0, slices)
            last_fraction = normalized_cell(slices - 1, slices)
            _x_first, p_first = focus_flow(
                slope, first_fraction * total_time, x_in, p_in
            )
            _x_last, p_last = focus_flow(
                slope, last_fraction * total_time, x_in, p_in
            )
            departure = side * p_first
            arrival = -side * p_last
            require_positive(departure, f"segment {segment} departure")
            require_positive(arrival, f"segment {segment} arrival")
            segment_velocity.extend((departure, arrival))

            # The first and last cells are controlled by the velocity signs
            # and the exact endpoint equations.  On every intervening cell,
            # prove that the excursion remains on the prescribed side.
            for cell_index in range(1, slices - 1):
                fraction = normalized_cell(cell_index, slices)
                x_cell, _p_cell = focus_flow(
                    slope, fraction * total_time, x_in, p_in
                )
                side_clearance = side * (x_cell - x_in)
                require_positive(
                    side_clearance,
                    f"segment {segment} turn-side clearance",
                )
                segment_position.append(side_clearance)

                # The short orbit turns in zone 1 and must not reach xi_2.
                if short and segment == 2:
                    upper_clearance = model.breaks[1] - x_cell
                    require_positive(
                        upper_clearance,
                        "short-orbit xi_2 clearance",
                    )
                    segment_position.append(upper_clearance)

        _x_out, p = focus_flow(slope, total_time, x_in, p_in)
        endpoints.append(p)
        velocity_guards.extend(segment_velocity)
        position_guards.extend(segment_position)
        segment_report: dict[str, object] = {
            "index": segment,
            "turn": segment in chart.turns,
            "time_enclosure": ball_text(total_time),
        }
        if segment_velocity:
            velocity_min = minimum_ball(segment_velocity)
            segment_report["velocity_guard_lower_bound"] = ball_text(
                velocity_min.lower()
            )
        if segment_position:
            position_min = minimum_ball(segment_position)
            segment_report["position_guard_lower_bound"] = ball_text(
                position_min.lower()
            )
        segment_reports.append(segment_report)

    crossing_margins: list[arb] = []
    for index, value in enumerate(endpoints):
        if value.contains(0):
            raise CertificateFailure(
                f"endpoint {index} transversality enclosure contains zero"
            )
        crossing_margins.append(abs(value))

    divergence = sum(
        (slope * time for slope, time in zip(chart.slopes, times)),
        arb(0),
    )
    if divergence.contains(0):
        raise CertificateFailure("divergence integral enclosure contains zero")
    multiplier = divergence.exp()
    if (multiplier - 1).contains(0):
        raise CertificateFailure("Floquet multiplier enclosure contains one")

    velocity_min = minimum_ball(velocity_guards)
    position_min = minimum_ball(position_guards)
    crossing_min = minimum_ball(crossing_margins)
    require_positive(velocity_min, "global velocity guard")
    require_positive(position_min, "global position guard")
    require_positive(crossing_min, "global crossing guard")

    return {
        "normalized_time_slices": slices,
        "first_event_branch_certified": True,
        "all_switches_transverse": True,
        "minimum_velocity_guard_lower_bound": ball_text(
            velocity_min.lower()
        ),
        "minimum_position_guard_lower_bound": ball_text(
            position_min.lower()
        ),
        "minimum_crossing_speed_lower_bound": ball_text(
            crossing_min.lower()
        ),
        "endpoint_p_enclosures": [
            ball_text(value) for value in endpoints
        ],
        "characteristic_exponent_enclosure": ball_text(divergence),
        "characteristic_exponent_excludes_zero": True,
        "floquet_multiplier_enclosure": ball_text(multiplier),
        "floquet_multiplier_excludes_one": True,
        "segments": segment_reports,
    }


def load_model(document: dict[str, object]) -> ExactModel:
    raw_breaks = document["breakpoints"]
    raw_slopes = document["slopes"]
    if not isinstance(raw_breaks, list) or not all(
        isinstance(value, str) for value in raw_breaks
    ):
        raise CertificateFailure("breakpoints must be decimal strings")
    if not isinstance(raw_slopes, list) or not all(
        isinstance(value, str) for value in raw_slopes
    ):
        raise CertificateFailure("slopes must be decimal strings")
    break_text = tuple(raw_breaks)
    slope_text = tuple(raw_slopes)
    if len(break_text) != 3 or len(slope_text) != 4:
        raise CertificateFailure("candidate must have three folds and four slopes")

    break_fraction = tuple(exact_fraction(value) for value in break_text)
    slope_fraction = tuple(exact_fraction(value) for value in slope_text)
    if not break_fraction[0] < break_fraction[1] < break_fraction[2]:
        raise CertificateFailure("breakpoints are not strictly increasing")
    if not all(
        slope_fraction[index] != slope_fraction[index + 1]
        for index in range(3)
    ):
        raise CertificateFailure("one of the advertised folds is absent")
    if not all(4 - slope * slope > 0 for slope in slope_fraction):
        raise CertificateFailure("not every zone has focus type")
    if not all(point > 0 for point in break_fraction):
        raise CertificateFailure("certificate assumes all switches are right of zero")

    intercepts = [Fraction(0)] * 4
    for index in range(3):
        intercepts[index + 1] = (
            intercepts[index]
            + (slope_fraction[index] - slope_fraction[index + 1])
            * break_fraction[index]
        )

    return ExactModel(
        break_text=break_text,  # type: ignore[arg-type]
        slope_text=slope_text,  # type: ignore[arg-type]
        breaks_q=tuple(exact_q(value) for value in break_text),  # type: ignore[arg-type]
        slopes_q=tuple(exact_q(value) for value in slope_text),  # type: ignore[arg-type]
        intercepts_fraction=tuple(intercepts),  # type: ignore[arg-type]
    )


def validate_candidate_metadata(document: dict[str, object]) -> None:
    if document.get("system") != "xdot=F(x)-y, ydot=x":
        raise CertificateFailure("unexpected system metadata")
    if document.get("normalization") != "F(0)=0":
        raise CertificateFailure("certificate requires F(0)=0")

    brackets = document.get("root_brackets")
    itineraries = document.get("expected_itineraries")
    if not isinstance(brackets, list) or len(brackets) != 4:
        raise CertificateFailure("candidate must contain four root brackets")
    if not isinstance(itineraries, list) or len(itineraries) != 4:
        raise CertificateFailure("candidate must contain four itineraries")

    expected_short = [[0, 1], [0, -1]]
    expected_full = [
        [0, 1],
        [1, 1],
        [2, 1],
        [2, -1],
        [1, -1],
        [0, -1],
    ]
    if itineraries[0] != expected_short:
        raise CertificateFailure("unexpected short itinerary")
    if any(itinerary != expected_full for itinerary in itineraries[1:]):
        raise CertificateFailure("unexpected full itinerary")

    for index, (center, bracket) in enumerate(
        zip(CENTERS, brackets)
    ):
        if not isinstance(bracket, list) or len(bracket) != 2:
            raise CertificateFailure(f"malformed root bracket {index}")
        left = exact_fraction(str(bracket[0]))
        right = exact_fraction(str(bracket[1]))
        amplitude = exact_fraction(center[0])
        # The numerical sign brackets can be much narrower than the
        # independently chosen Krawczyk box.  They locate the seed only and
        # play no role in the interval acceptance test.
        if not left < amplitude < right:
            raise CertificateFailure(f"cycle {index + 1} seed misses bracket")


def certify(
    document: dict[str, object],
    *,
    precision: int,
    radius: str,
    slices: int,
    candidate_sha256: str,
) -> dict[str, object]:
    ctx.prec = precision
    # mpmath only builds C; these digits do not enter an acceptance test.
    mp.mp.dps = max(170, precision // 3 + 20)

    validate_candidate_metadata(document)
    model = load_model(document)
    root_reports: list[dict[str, object]] = []
    section_boxes: list[arb] = []

    for index, centers in enumerate(CENTERS):
        short = index == 0
        chart = build_chart(short, model)
        krawczyk, boxes = krawczyk_check(
            short, centers, radius, model
        )
        branch = branch_and_multiplier_check(
            short, centers, radius, slices, model
        )
        section_boxes.append(boxes[0])
        root_reports.append(
            {
                "cycle": index + 1,
                "chart": "short" if short else "full",
                "itinerary": chart.itinerary,
                "section_amplitude_enclosure": ball_text(boxes[0]),
                "positive_section_amplitude": True,
                "all_flight_times_positive": True,
                "krawczyk": krawczyk,
                "branch": branch,
            }
        )

    for left in range(len(section_boxes)):
        for right in range(left + 1, len(section_boxes)):
            if section_boxes[left].overlaps(section_boxes[right]):
                raise CertificateFailure(
                    f"section boxes {left + 1} and {right + 1} overlap"
                )

    return {
        "status": "rigorous-directed-ball-certificate-pass",
        "claim": (
            "At least four distinct hyperbolic limit cycles exist for "
            "the exact rational candidate parameters."
        ),
        "scope": (
            "Existence, distinctness, prescribed first-event itineraries, "
            "switch transversality, and hyperbolicity are certified; this "
            "does not assert that there are exactly four cycles."
        ),
        "arithmetic": {
            "backend": "python-flint Arb directed real ball arithmetic",
            "python_flint_version": flint.__version__,
            "flint_version": flint.__FLINT_VERSION__,
            "precision_bits": precision,
            "decimal_parameters": "interpreted as exact rationals",
            "box_radius": radius,
            "normalized_time_slices": slices,
            "acceptance_uses_binary_float": False,
        },
        "candidate_sha256": candidate_sha256,
        "exact_parameters": {
            "breakpoints": list(model.break_text),
            "slopes": list(model.slope_text),
            "intercepts_as_exact_fractions": [
                str(value) for value in model.intercepts_fraction
            ],
            "all_four_zones_focus_type": True,
            "all_three_folds_genuine": True,
            "continuous_at_switches_by_intercept_construction": True,
        },
        "pairwise_disjoint_section_boxes": True,
        "cycles": root_reports,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "candidate",
        nargs="?",
        type=Path,
        default=Path("experiments/candidates/four-cycle-001.json"),
    )
    parser.add_argument("--precision", type=int, default=DEFAULT_PRECISION)
    parser.add_argument("--radius", default=DEFAULT_RADIUS)
    parser.add_argument("--slices", type=int, default=DEFAULT_SLICES)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/generated/four-cycle-001-arb-certificate.json"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.precision < 256:
        raise SystemExit("use at least 256 bits of Arb precision")
    if args.slices < 64:
        raise SystemExit("use at least 64 normalized-time slices")
    radius = exact_fraction(args.radius)
    if radius <= 0:
        raise SystemExit("box radius must be positive")

    candidate_bytes = args.candidate.read_bytes()
    document = json.loads(candidate_bytes.decode("utf-8"))
    try:
        report = certify(
            document,
            precision=args.precision,
            radius=args.radius,
            slices=args.slices,
            candidate_sha256=hashlib.sha256(candidate_bytes).hexdigest(),
        )
    except CertificateFailure as error:
        print(f"CERTIFICATE FAIL: {error}", file=sys.stderr)
        return 1

    report["checker_sha256"] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()

    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(args.output)

    print("CERTIFICATE PASS")
    print(report["claim"])
    print(f"report: {args.output}")
    for cycle in report["cycles"]:
        multiplier = cycle["branch"]["floquet_multiplier_enclosure"]
        print(
            f"cycle {cycle['cycle']}: "
            f"s={cycle['section_amplitude_enclosure']}; "
            f"multiplier={multiplier}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
