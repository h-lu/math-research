#!/usr/bin/env python3
"""Fresh floating-point harvesting-fold diagnostic; JSON goes to stdout.

Usage (from the repository root):
    python3 archival-verification/scripts/diagnose_harvesting.py
    python3 archival-verification/scripts/diagnose_harvesting.py --csv orbit.csv

This is a numerical diagnostic, not a validated existence or uniqueness proof.
No stored orbit or precomputed solution is read. See notes/HARVESTING-DIAGNOSTIC.md.
"""

if not __debug__:
    raise SystemExit("Optimized Python is unsupported; run without -O/-OO.")

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


U = 21.0 / 100.0
EFFORT = 1.0e-6
STARTUP_TIME = 0.01
END_TIME_LIMIT = 100.0
SAMPLE_COUNT = 1001
RESIDUAL_SCALES = np.array([1.0e-4, 1.0e-2, 1.0e-2])


@dataclass(frozen=True)
class Accuracy:
    rtol: float
    atol: float
    max_step: float


STANDARD = Accuracy(1.0e-12, 1.0e-14, 0.1)
TIGHTER = Accuracy(3.0e-14, 3.0e-16, 0.05)


def physical_parameters(parameters):
    beta, nu, lower = parameters
    y_equilibrium = 4.0 / 5.0 - EFFORT * (1.0 - 5.0 * lower)
    d = 8.0 * beta + 1.0 / 5.0 + nu
    h = y_equilibrium * (d - 5.0 * beta * y_equilibrium)
    return beta, d, h, lower, y_equilibrium


def field_and_jacobian(state, parameters):
    """Derivatives of the original-time smooth middle field (x > 0)."""
    x, y = state[:2]
    beta, d, h, lower, _ = physical_parameters(parameters)
    vector = np.array([
        x * (1.0 - x - y) - EFFORT * (x - lower),
        y * (d - beta * y / x) - h,
    ])
    jacobian = np.array([
        [1.0 - 2.0 * x - y - EFFORT, -x],
        [beta * y * y / (x * x), d - 2.0 * beta * y / x],
    ])
    return vector, jacobian


def initial_state(parameters):
    lower = parameters[2]
    y = 1.0 - U - EFFORT * (U - lower) / U
    return np.array([U, y, 0.0])


def integrate_one_turn(parameters, accuracy):
    def rhs(_time, state):
        vector, jacobian = field_and_jacobian(state, parameters)
        return np.array([vector[0], vector[1], np.trace(jacobian)])

    def maximum(_time, state):
        return field_and_jacobian(state, parameters)[0][0]

    maximum.direction = -1
    maximum.terminal = True

    def minimum(_time, state):
        return field_and_jacobian(state, parameters)[0][0]

    minimum.direction = 1
    minimum.terminal = False
    options = dict(method="DOP853", rtol=accuracy.rtol,
                   atol=accuracy.atol, max_step=accuracy.max_step,
                   dense_output=True)
    # The first maximum is the initial condition. Move off it before enabling
    # the decreasing-F event, to avoid terminating on the t = 0 event.
    startup = solve_ivp(rhs, (0.0, STARTUP_TIME), initial_state(parameters),
                        **options)
    if not startup.success:
        raise RuntimeError("Startup integration failed")
    turn = solve_ivp(rhs, (STARTUP_TIME, END_TIME_LIMIT), startup.y[:, -1],
                     events=(maximum, minimum), **options)
    if not turn.success:
        raise RuntimeError("One-turn integration failed")
    if len(turn.t_events[0]) != 1 or len(turn.t_events[1]) != 1:
        raise RuntimeError("Expected one detected minimum followed by a maximum")
    if not STARTUP_TIME < turn.t_events[1][0] < turn.t_events[0][0]:
        raise RuntimeError("Unexpected extremum event order")
    endpoint = turn.y_events[0][0]
    bottom = turn.y_events[1][0]
    residual = np.array([endpoint[0] - U, endpoint[2],
                         bottom[0] - parameters[2]])
    return residual, startup, turn


def samples(startup, turn):
    period = turn.t_events[0][0]
    times = np.linspace(0.0, period, SAMPLE_COUNT)
    data = np.empty((SAMPLE_COUNT, 3))
    before = times <= STARTUP_TIME
    data[before] = startup.sol(times[before]).T
    data[~before] = turn.sol(times[~before]).T
    return times, data


def variational_check(parameters, period):
    """Propagate the full fundamental matrix independently, at fixed time."""
    def rhs(_time, state):
        vector, jacobian = field_and_jacobian(state, parameters)
        matrix = state[3:].reshape(2, 2)
        return np.concatenate((vector, [np.trace(jacobian)],
                               (jacobian @ matrix).ravel()))

    start = np.concatenate((initial_state(parameters), np.eye(2).ravel()))
    result = solve_ivp(rhs, (0.0, period), start, method="DOP853",
                       rtol=TIGHTER.rtol, atol=TIGHTER.atol,
                       max_step=TIGHTER.max_step)
    if not result.success:
        raise RuntimeError("Variational integration failed")
    endpoint = result.y[:, -1]
    matrix = endpoint[3:].reshape(2, 2)
    determinant = np.linalg.det(matrix)
    liouville_error = determinant - np.exp(endpoint[2])
    # The initial curve F = 0 has y'(x) = -1-uL/x^2. At the terminal
    # maximum, F = 0, so the event-time correction has zero x component.
    tangent = np.array([1.0, -1.0 - EFFORT * parameters[2] / (U * U)])
    section_derivative = (matrix @ tangent)[0]
    assert abs(liouville_error) < 1.0e-10, "Liouville mismatch"
    assert abs(section_derivative - 1.0) < 1.0e-9, "Return derivative mismatch"
    assert np.max(np.abs(endpoint[:2] - start[:2])) < 1.0e-10
    return {
        "fundamental_matrix": matrix.tolist(),
        "fundamental_matrix_determinant": determinant,
        "liouville_difference": liouville_error,
        "maximum_section_return_derivative": section_derivative,
        "integrated_divergence": endpoint[2],
    }


def rounded(value):
    """Stable compact output; underlying acceptance checks use full precision."""
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [rounded(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(format(value, ".12g"))
    return value


def run_diagnostic(csv_path=None):
    # Initial guesses are only the rounded caption parameters; the trajectory,
    # return conditions and solution are recomputed on every invocation.
    beta_guess = 0.031443254
    d_guess = 0.451546877
    guess = np.array([beta_guess, d_guess - 8.0 * beta_guess - 0.2,
                      0.191154841])

    def residual(parameters):
        return integrate_one_turn(parameters, STANDARD)[0] / RESIDUAL_SCALES

    solution = root(residual, guess, method="hybr", tol=1.0e-8,
                    options={"eps": 1.0e-7})
    if not solution.success:
        raise RuntimeError("Nonlinear solve did not converge")
    parameters = solution.x
    normal_residual, _, normal_turn = integrate_one_turn(parameters, STANDARD)
    tighter_residual, startup, turn = integrate_one_turn(parameters, TIGHTER)
    assert np.max(np.abs(normal_residual)) < 1.0e-10
    assert np.max(np.abs(tighter_residual)) < 1.0e-10
    period = turn.t_events[0][0]
    assert abs(period - normal_turn.t_events[0][0]) < 1.0e-8
    beta, d, h, lower, y_equilibrium = physical_parameters(parameters)
    assert beta > 0.0 and d > 0.0 and h > 0.0 and 0.0 < lower < U < 1.0
    assert np.max(np.abs(field_and_jacobian([0.2, y_equilibrium], parameters)[0])) < 1.0e-14
    times, data = samples(startup, turn)
    assert np.all(data[:, :2] > 0.0)
    assert np.min(data[:, 0]) > lower - 1.0e-10
    assert np.max(data[:, 0]) < U + 1.0e-10
    initial = initial_state(parameters)
    endpoint = turn.y_events[0][0]
    assert np.max(np.abs(endpoint[:2] - initial[:2])) < 1.0e-10

    def contact(state):
        vector, _ = field_and_jacobian(state, parameters)
        return {"x": state[0], "y": state[1], "G": vector[1],
                "x_second_time_derivative": -state[0] * vector[1]}

    bottom = contact(turn.y_events[1][0])
    top = contact(initial)
    assert bottom["x_second_time_derivative"] > 0.0
    assert top["x_second_time_derivative"] < 0.0
    # A comparison with the caption's displayed precision, not fitted digits.
    assert abs(beta - beta_guess) < 0.5e-9
    assert abs(d - d_guess) < 0.5e-9
    assert abs(h - 0.260619080) < 0.5e-9
    assert abs(lower - 0.191154841) < 0.5e-9

    report = {
        "status": "PASS: nonvalidated floating-point diagnostic",
        "source": "source/sections/08-realization.tex: smooth middle field and perturbed chart",
        "method": "SciPy DOP853; scipy.optimize.root(method=hybr)",
        "unknowns": ["beta", "nu", "L"],
        "parameters": {"u": EFFORT, "U": U, "beta": beta,
                       "nu": parameters[1], "L": lower, "d": d, "h": h},
        "equilibrium": {"x": 0.2, "y": y_equilibrium},
        "initial_maximum": top,
        "minimum": bottom,
        "minimum_time": turn.t_events[1][0],
        "period": period,
        "standard_accuracy": vars(STANDARD),
        "tighter_accuracy": vars(TIGHTER),
        "residual_order": ["maximum_return_x-U", "integrated_divergence", "minimum_x-L"],
        "standard_residual": normal_residual.tolist(),
        "tighter_residual": tighter_residual.tolist(),
        "period_tightening_difference": period - normal_turn.t_events[0][0],
        "multiplier_from_divergence": np.exp(tighter_residual[1]),
        "variational": variational_check(parameters, period),
        "sample_count": SAMPLE_COUNT,
        "sampled_positive_quadrant_margin": float(np.min(data[:, :2])),
        "sampled_y_range": [float(np.min(data[:, 1])), float(np.max(data[:, 1]))],
        "caption_rounding_comparison": "agrees for beta,d,h,L at 9 digits after the decimal point",
        "acceptance_bounds": {"absolute_return_residual": 1.0e-10,
                              "absolute_integrated_divergence": 1.0e-10,
                              "absolute_period_tightening_change": 1.0e-8,
                              "absolute_variational_return_derivative_error": 1.0e-9},
        "limitations": [
            "No interval arithmetic or rigorous integration error enclosure.",
            "No proof of existence, uniqueness, nondegeneracy, or exact neutrality of a cyclic fold.",
            "No certified observation rank, exact six-cycle configuration, or all-parameter estimate.",
            "Event detection and finite sampling do not prove global extremum counts or domain bounds.",
            "The upper effort v is unused because this orbit is integrated in the smooth middle extension.",
        ],
    }
    if csv_path is not None:
        with Path(csv_path).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["phase", "time", "x", "y", "integral_divergence"])
            for time, state in zip(times, data):
                writer.writerow(format(item, ".12g")
                                for item in [time / period, time, *state])
    return rounded(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", help="Optionally write 1001 fresh orbit samples")
    arguments = parser.parse_args()
    print(json.dumps(run_diagnostic(arguments.csv), indent=2, sort_keys=True,
                     allow_nan=False))


if __name__ == "__main__":
    main()
