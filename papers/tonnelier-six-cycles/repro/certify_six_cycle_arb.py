#!/usr/bin/env python3
"""Rigorous directed-Arb certificate for the q=2 six-cycle candidate.

The section x=0 lies in zone 2.  One proposed orbit turns in zone 1 and
five proposed orbits visit all four zones.  Decimal parameters and box
centres are interpreted as exact rationals.  Floating-point-style mpmath
values are used only to construct fixed dyadic Krawczyk preconditioners;
every acceptance decision is made with directed Arb ball arithmetic.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import flint
from flint import arb
import mpmath as mp

import certify_candidate_arb as base


SHORT_ITINERARY = [[1, -1], [1, 1]]
FULL_ITINERARY = [
    [1, -1],
    [0, -1],
    [0, 1],
    [1, 1],
    [2, 1],
    [2, -1],
]


@dataclass(frozen=True)
class Chart:
    name: str
    slope_indices: tuple[int, ...]
    node_indices: tuple[int | None, ...]
    turns: dict[int, int]
    itinerary: list[list[int]]


SHORT_CHART = Chart(
    name="short",
    slope_indices=(2, 1, 2),
    node_indices=(None, 1, 1, None),
    # Segment 2 leaves xi_2 to the right, turns inside zone 2, and returns
    # to the oriented section x=0 with p<0.
    turns={1: -1, 2: 1},
    itinerary=SHORT_ITINERARY,
)
FULL_CHART = Chart(
    name="full",
    slope_indices=(2, 1, 0, 1, 2, 3, 2),
    node_indices=(None, 1, 0, 0, 1, 2, 2, None),
    turns={2: -1, 5: 1},
    itinerary=FULL_ITINERARY,
)


def chart_values(
    chart: Chart, model: base.ExactModel
) -> tuple[list[arb], list[arb]]:
    slopes = model.slopes
    breaks = model.breaks
    return (
        [slopes[index] for index in chart.slope_indices],
        [arb(0) if index is None else breaks[index] for index in chart.node_indices],
    )


def mp_chart_values(
    chart: Chart, model: base.ExactModel
) -> tuple[list[mp.mpf], list[mp.mpf]]:
    slopes = [mp.mpf(value) for value in model.slope_text]
    breaks = [mp.mpf(value) for value in model.break_text]
    return (
        [slopes[index] for index in chart.slope_indices],
        [mp.mpf(0) if index is None else breaks[index] for index in chart.node_indices],
    )


def residual_dual(
    chart: Chart, values: list[arb], model: base.ExactModel
) -> list[base.Dual]:
    slopes, nodes = chart_values(chart, model)
    dimension = len(values)
    variables: list[base.Dual] = []
    for index, value in enumerate(values):
        derivative = [arb(0) for _ in range(dimension)]
        derivative[index] = arb(1)
        variables.append(base.Dual(value, derivative))

    # F(0)=0, y=s>0 on the section, hence p=F-y=-s.
    p = -variables[0]
    residuals: list[base.Dual] = []
    for slope, time, x_in, x_target in zip(
        slopes, variables[1:], nodes[:-1], nodes[1:]
    ):
        x_out, p = base.focus_flow(slope, time, x_in, p)
        residuals.append(x_out - x_target)
    residuals.append(p + variables[0])
    return residuals


def mp_residual(chart: Chart, values, model: base.ExactModel):
    slopes, nodes = mp_chart_values(chart, model)
    p = -values[0]
    residuals = []
    for slope, time, x_in, x_target in zip(
        slopes, values[1:], nodes[:-1], nodes[1:]
    ):
        x_out, p = base.mp_flow(slope, time, x_in, p)
        residuals.append(x_out - x_target)
    residuals.append(p + values[0])
    return residuals


def krawczyk_check(
    chart: Chart,
    centers: tuple[str, ...],
    radius: str,
    model: base.ExactModel,
) -> tuple[dict[str, object], list[arb]]:
    dimension = len(centers)
    expected_dimension = len(chart.slope_indices) + 1
    if dimension != expected_dimension:
        raise base.CertificateFailure(
            f"wrong {chart.name} chart dimension {dimension}, "
            f"expected {expected_dimension}"
        )
    for index, value in enumerate(centers):
        label = "section amplitude" if index == 0 else f"time {index}"
        base.prove_box_positive(value, radius, label)

    mp_centers = [mp.mpf(value) for value in centers]
    jacobian = mp.matrix(dimension)
    for row in range(dimension):
        for column in range(dimension):

            def component(value, row=row, column=column):
                variables = mp_centers.copy()
                variables[column] = value
                return mp_residual(chart, variables, model)[row]

            jacobian[row, column] = mp.diff(component, mp_centers[column])
    inverse = jacobian**-1
    preconditioner = [
        [base.arb_dyadic_from_mp(inverse[i, j]) for j in range(dimension)]
        for i in range(dimension)
    ]

    full_box = [base.around(value, radius) for value in centers]
    center_points = [box.mid() for box in full_box]
    exact_bounds: list[tuple[arb, arb]] = []
    for box, value in zip(full_box, centers):
        lower_fraction, upper_fraction = base.exact_box_data(value, radius)
        lower = arb(base.fmpq(lower_fraction.numerator, lower_fraction.denominator))
        upper = arb(base.fmpq(upper_fraction.numerator, upper_fraction.denominator))
        if not box.contains(lower) or not box.contains(upper):
            raise base.CertificateFailure("Arb box did not enclose rational box")
        exact_bounds.append((lower, upper))

    point_residual = [
        item.value for item in residual_dual(chart, center_points, model)
    ]
    interval_residual = residual_dual(chart, full_box, model)
    interval_jacobian = [
        [interval_residual[i].derivative[j] for j in range(dimension)]
        for i in range(dimension)
    ]
    product = base.matmul(preconditioner, interval_jacobian)
    defect = [
        [
            (arb(1) if i == j else arb(0)) - product[i][j]
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    correction = base.matvec(preconditioner, point_residual)
    centered_box = [
        full_box[i] - center_points[i] for i in range(dimension)
    ]
    defect_box = base.matvec(defect, centered_box)
    krawczyk = [
        center_points[i] - correction[i] + defect_box[i]
        for i in range(dimension)
    ]

    inclusion_margins: list[arb] = []
    for enclosure, (lower, upper) in zip(krawczyk, exact_bounds):
        left_margin = enclosure - lower
        right_margin = upper - enclosure
        base.require_positive(left_margin, "left Krawczyk inclusion margin")
        base.require_positive(right_margin, "right Krawczyk inclusion margin")
        inclusion_margins.extend((left_margin, right_margin))

    row_sums = [
        sum((abs(entry).upper() for entry in row), arb(0)) for row in defect
    ]
    if not all(row_sum < 1 for row_sum in row_sums):
        raise base.CertificateFailure("could not prove Krawczyk contraction")

    q_upper = base.maximum_ball(row_sums)
    max_residual = base.maximum_ball(
        abs(value).upper() for value in point_residual
    )
    min_inclusion = base.minimum_ball(inclusion_margins)
    return (
        {
            "dimension": dimension,
            "flight_count": dimension - 1,
            "strict_inclusion": True,
            "contraction_norm_strictly_below_one": True,
            "defect_infinity_norm_upper_enclosure": base.ball_text(q_upper),
            "point_residual_max_abs_enclosure": base.ball_text(max_residual),
            "strict_inclusion_margin_enclosure": base.ball_text(min_inclusion),
            "strict_inclusion_margin_lower_bound": base.ball_text(
                min_inclusion.lower()
            ),
        },
        full_box,
    )


def branch_and_multiplier_check(
    chart: Chart,
    centers: tuple[str, ...],
    radius: str,
    slices: int,
    model: base.ExactModel,
) -> dict[str, object]:
    slopes, nodes = chart_values(chart, model)
    variables = [base.around(value, radius) for value in centers]
    times = variables[1:]
    p = -variables[0]
    endpoints: list[arb] = []
    velocity_guards: list[arb] = []
    position_guards: list[arb] = []
    segment_reports: list[dict[str, object]] = []

    for segment, (slope, total_time, x_in, x_target) in enumerate(
        zip(slopes, times, nodes[:-1], nodes[1:])
    ):
        p_in = p
        segment_velocity: list[arb] = []
        segment_position: list[arb] = []
        if segment not in chart.turns:
            direction = 1 if x_target > x_in else -1
            for cell_index in range(slices):
                fraction = base.normalized_cell(cell_index, slices)
                _x_cell, p_cell = base.focus_flow(
                    slope, fraction * total_time, x_in, p_in
                )
                signed_velocity = direction * p_cell
                base.require_positive(
                    signed_velocity, f"segment {segment} monotonicity guard"
                )
                segment_velocity.append(signed_velocity)
        else:
            side = chart.turns[segment]
            cells = [
                base.focus_flow(
                    slope,
                    base.normalized_cell(cell_index, slices) * total_time,
                    x_in,
                    p_in,
                )
                for cell_index in range(slices)
            ]
            _x_first, p_first = cells[0]
            _x_last, p_last = cells[-1]
            departure = side * p_first
            arrival = -side * p_last
            base.require_positive(departure, f"segment {segment} departure")
            base.require_positive(arrival, f"segment {segment} arrival")
            segment_velocity.extend((departure, arrival))

            # Near a boundary endpoint direct evaluation of x(t)-x(0)
            # suffers avoidable interval cancellation.  Cover the maximal
            # initial block by the rigorously signed velocity instead.  For
            # a return to the same switch, do the analogous argument
            # backwards from the exact terminal endpoint.  Direct position
            # enclosures cover every remaining cell.
            prefix = 0
            for _x_cell, p_cell in cells:
                signed_velocity = side * p_cell
                if not signed_velocity > 0:
                    break
                segment_velocity.append(signed_velocity)
                prefix += 1
            same_endpoint = (
                chart.node_indices[segment]
                == chart.node_indices[segment + 1]
            )
            suffix = 0
            if same_endpoint:
                for _x_cell, p_cell in reversed(cells):
                    signed_velocity = -side * p_cell
                    if not signed_velocity > 0:
                        break
                    segment_velocity.append(signed_velocity)
                    suffix += 1

            for cell_index, (x_cell, _p_cell) in enumerate(cells):
                velocity_covered = (
                    cell_index < prefix
                    or (
                        same_endpoint
                        and cell_index >= slices - suffix
                    )
                )
                if not velocity_covered:
                    side_clearance = side * (x_cell - x_in)
                    base.require_positive(
                        side_clearance,
                        f"segment {segment} turn-side clearance",
                    )
                    segment_position.append(side_clearance)

                # The short orbit turns inside zone 1: it must remain to the
                # right of xi_1, rather than silently entering zone 0.
                if chart.name == "short" and segment == 1:
                    competitor_clearance = x_cell - model.breaks[0]
                    base.require_positive(
                        competitor_clearance,
                        "short-orbit xi_1 clearance",
                    )
                    segment_position.append(competitor_clearance)
                if chart.name == "short" and segment == 2:
                    competitor_clearance = model.breaks[2] - x_cell
                    base.require_positive(
                        competitor_clearance,
                        "short-orbit xi_3 clearance",
                    )
                    segment_position.append(competitor_clearance)

        _x_out, p = base.focus_flow(slope, total_time, x_in, p_in)
        endpoints.append(p)
        velocity_guards.extend(segment_velocity)
        position_guards.extend(segment_position)
        report: dict[str, object] = {
            "index": segment,
            "turn": segment in chart.turns,
            "time_enclosure": base.ball_text(total_time),
        }
        if segment_velocity:
            report["velocity_guard_lower_bound"] = base.ball_text(
                base.minimum_ball(segment_velocity).lower()
            )
        if segment_position:
            report["position_guard_lower_bound"] = base.ball_text(
                base.minimum_ball(segment_position).lower()
            )
        segment_reports.append(report)

    crossing_margins: list[arb] = []
    for index, value in enumerate(endpoints):
        if value.contains(0):
            raise base.CertificateFailure(
                f"endpoint {index} transversality enclosure contains zero"
            )
        crossing_margins.append(abs(value))

    divergence = sum(
        (slope * time for slope, time in zip(slopes, times)), arb(0)
    )
    if divergence.contains(0):
        raise base.CertificateFailure("divergence integral enclosure contains zero")
    multiplier = divergence.exp()
    if (multiplier - 1).contains(0):
        raise base.CertificateFailure("Floquet multiplier enclosure contains one")

    velocity_min = base.minimum_ball(velocity_guards)
    position_min = base.minimum_ball(position_guards)
    crossing_min = base.minimum_ball(crossing_margins)
    base.require_positive(velocity_min, "global velocity guard")
    base.require_positive(position_min, "global position guard")
    base.require_positive(crossing_min, "global crossing guard")
    return {
        "normalized_time_slices": slices,
        "first_event_branch_certified": True,
        "all_switches_transverse": True,
        "minimum_velocity_guard_lower_bound": base.ball_text(
            velocity_min.lower()
        ),
        "minimum_position_guard_lower_bound": base.ball_text(
            position_min.lower()
        ),
        "minimum_crossing_speed_lower_bound": base.ball_text(
            crossing_min.lower()
        ),
        "endpoint_p_enclosures": [base.ball_text(value) for value in endpoints],
        "characteristic_exponent_enclosure": base.ball_text(divergence),
        "characteristic_exponent_excludes_zero": True,
        "floquet_multiplier_enclosure": base.ball_text(multiplier),
        "floquet_multiplier_excludes_one": True,
        "segments": segment_reports,
    }


def load_model(document: dict[str, object]) -> base.ExactModel:
    raw_breaks = document.get("breakpoints")
    raw_slopes = document.get("slopes")
    if not isinstance(raw_breaks, list) or not all(
        isinstance(value, str) for value in raw_breaks
    ):
        raise base.CertificateFailure("breakpoints must be decimal strings")
    if not isinstance(raw_slopes, list) or not all(
        isinstance(value, str) for value in raw_slopes
    ):
        raise base.CertificateFailure("slopes must be decimal strings")
    if len(raw_breaks) != 3 or len(raw_slopes) != 4:
        raise base.CertificateFailure("candidate must have three folds and four slopes")

    break_text = tuple(raw_breaks)
    slope_text = tuple(raw_slopes)
    breaks = tuple(base.exact_fraction(value) for value in break_text)
    slopes = tuple(base.exact_fraction(value) for value in slope_text)
    if not breaks[0] < breaks[1] < 0 < breaks[2]:
        raise base.CertificateFailure("certificate requires xi_1 < xi_2 < 0 < xi_3")
    if not all(slopes[i] != slopes[i + 1] for i in range(3)):
        raise base.CertificateFailure("one of the advertised folds is absent")
    if not all(4 - slope * slope > 0 for slope in slopes):
        raise base.CertificateFailure("not every zone has focus type")

    # Anchor F(0)=0 in zone 2 and propagate continuity in both directions.
    intercepts = [Fraction(0) for _ in range(4)]
    anchor = 2
    for index in range(anchor, 3):
        intercepts[index + 1] = (
            intercepts[index]
            + (slopes[index] - slopes[index + 1]) * breaks[index]
        )
    for index in range(anchor - 1, -1, -1):
        intercepts[index] = (
            intercepts[index + 1]
            + (slopes[index + 1] - slopes[index]) * breaks[index]
        )

    return base.ExactModel(
        break_text=break_text,  # type: ignore[arg-type]
        slope_text=slope_text,  # type: ignore[arg-type]
        breaks_q=tuple(base.exact_q(value) for value in break_text),  # type: ignore[arg-type]
        slopes_q=tuple(base.exact_q(value) for value in slope_text),  # type: ignore[arg-type]
        intercepts_fraction=tuple(intercepts),  # type: ignore[arg-type]
    )


def load_centers(document: dict[str, object]) -> tuple[tuple[str, ...], ...]:
    raw = document.get("krawczyk_centers")
    brackets = document.get("root_brackets")
    itineraries = document.get("expected_itineraries")
    if not isinstance(raw, list) or len(raw) != 6:
        raise base.CertificateFailure("six-cycle candidate must contain six centers")
    if not isinstance(brackets, list) or len(brackets) != len(raw):
        raise base.CertificateFailure("root bracket count does not match centers")
    if not isinstance(itineraries, list) or len(itineraries) != len(raw):
        raise base.CertificateFailure("itinerary count does not match centers")
    if itineraries[0] != SHORT_ITINERARY:
        raise base.CertificateFailure("first cycle is not the short itinerary")
    if any(value != FULL_ITINERARY for value in itineraries[1:]):
        raise base.CertificateFailure("a proposed outer cycle is not full")

    centers: list[tuple[str, ...]] = []
    for index, (row, bracket) in enumerate(zip(raw, brackets)):
        expected = 4 if index == 0 else 8
        if (
            not isinstance(row, list)
            or len(row) != expected
            or not all(isinstance(value, str) for value in row)
        ):
            raise base.CertificateFailure(f"malformed Krawczyk center {index + 1}")
        if not isinstance(bracket, list) or len(bracket) != 2:
            raise base.CertificateFailure(f"malformed root bracket {index + 1}")
        left, right = (Fraction(str(value)) for value in bracket)
        if not left < Fraction(row[0]) < right:
            raise base.CertificateFailure(
                f"cycle {index + 1} seed misses its root bracket"
            )
        centers.append(tuple(row))
    return tuple(centers)


def validate_metadata(document: dict[str, object]) -> None:
    if document.get("system") != "xdot=F(x)-y, ydot=x":
        raise base.CertificateFailure("unexpected system metadata")
    if document.get("normalization") != "F(0)=0":
        raise base.CertificateFailure("certificate requires F(0)=0")
    if document.get("origin_zone") != 2:
        raise base.CertificateFailure("certificate requires origin_zone=2")


def certify(
    document: dict[str, object],
    *,
    precision: int,
    radius: str,
    slices: int,
    candidate_sha256: str,
) -> dict[str, object]:
    base.ctx.prec = precision
    mp.mp.dps = max(170, precision // 3 + 20)
    validate_metadata(document)
    centers = load_centers(document)
    model = load_model(document)

    root_reports: list[dict[str, object]] = []
    section_boxes: list[arb] = []
    for index, center in enumerate(centers):
        chart = SHORT_CHART if index == 0 else FULL_CHART
        krawczyk, boxes = krawczyk_check(chart, center, radius, model)
        branch = branch_and_multiplier_check(
            chart, center, radius, slices, model
        )
        section_boxes.append(boxes[0])
        root_reports.append(
            {
                "cycle": index + 1,
                "chart": chart.name,
                "itinerary": chart.itinerary,
                "section_amplitude_enclosure": base.ball_text(boxes[0]),
                "positive_section_amplitude": True,
                "all_flight_times_positive": True,
                "krawczyk": krawczyk,
                "branch": branch,
            }
        )

    for left in range(len(section_boxes)):
        for right in range(left + 1, len(section_boxes)):
            if section_boxes[left].overlaps(section_boxes[right]):
                raise base.CertificateFailure(
                    f"section boxes {left + 1} and {right + 1} overlap"
                )

    return {
        "status": "rigorous-directed-ball-certificate-pass",
        "claim": (
            "At least six distinct hyperbolic limit cycles exist for "
            "the exact rational candidate parameters."
        ),
        "scope": (
            "Existence, distinctness, prescribed first-event itineraries, "
            "switch transversality, and hyperbolicity are certified; this "
            "does not assert that there are exactly six cycles."
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
            "origin_zone": 2,
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
        default=Path("experiments/candidates/six-cycle-001.json"),
    )
    parser.add_argument("--precision", type=int, default=480)
    parser.add_argument("--radius", default="1e-25")
    parser.add_argument("--slices", type=int, default=512)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/generated/six-cycle-001-arb-certificate.json"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.precision < 256:
        raise SystemExit("use at least 256 bits of Arb precision")
    if args.slices < 64:
        raise SystemExit("use at least 64 normalized-time slices")
    if Fraction(args.radius) <= 0:
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
    except base.CertificateFailure as error:
        print(f"CERTIFICATE FAIL: {error}", file=sys.stderr)
        return 1

    report["checker_sha256"] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    report["core_checker_sha256"] = hashlib.sha256(
        Path(base.__file__).read_bytes()
    ).hexdigest()
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(output, encoding="utf-8")
    temporary.replace(args.output)

    print("CERTIFICATE PASS")
    print(report["claim"])
    print(f"report: {args.output}")
    for cycle in report["cycles"]:
        print(
            f"cycle {cycle['cycle']}: "
            f"s={cycle['section_amplitude_enclosure']}; "
            f"multiplier={cycle['branch']['floquet_multiplier_enclosure']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
