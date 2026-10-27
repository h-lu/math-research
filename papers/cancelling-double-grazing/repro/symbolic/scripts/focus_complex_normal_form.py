#!/usr/bin/env python3
"""Independent physical-time complex normal form of the harvesting core.

This uses polynomial Lie transforms of the vector field, not the energy
homological equations in source/scripts/check-focus.py.  Run from any folder.
Only exact arithmetic is used.  See notes/FOCUS-COMPLEX.md for conventions.
"""

from __future__ import annotations

if not __debug__:
    raise SystemExit("Optimized Python is unsupported; run without -O/-OO.")

import json

import sympy as s
from sympy.polys.rings import ring


def require(condition, message):
    """Raise on a failed comparison without relying on assert."""
    if not condition:
        raise ValueError(message)


def compute(max_degree=5):
    w = s.Symbol("w", positive=True)
    field = s.QQ_I.frac_field(w)
    polynomials, z, zb = ring("z,zb", field)
    W = field.from_sympy(w)
    I = field.from_sympy(s.I)

    def scalar(p, q=1):
        return field(p) / field(q)

    def truncate(p, degree=max_degree):
        return polynomials.from_dict({m: c for m, c in p.items() if sum(m) <= degree})

    def homogeneous(p, degree):
        return polynomials.from_dict({m: c for m, c in p.items() if sum(m) == degree})

    def conjugate(p):
        # w is real; complex conjugation also interchanges z and zb.
        return polynomials.from_dict({
            (m[1], m[0]): field.from_sympy(s.conjugate(field.to_sympy(c)))
            for m, c in p.items()
        })

    def bracket(h, f):
        # [h,f] = Df h - Dh f, so exp(ad_h) f is the pullback by
        # the time-one map of h.  This is a state change, not a clock change.
        return tuple(truncate(
            f[j].diff(z) * h[0] + f[j].diff(zb) * h[1]
            - h[j].diff(z) * f[0] - h[j].diff(zb) * f[1]
        ) for j in range(2))

    def show(c):
        return s.sstr(s.factor(field.to_sympy(c)))

    def terms(p):
        return {f"z^{j}*zb^{k}": show(c) for (j, k), c in sorted(p.items())}

    # Inputs are the original model equations, with u=nu=0.  The rational
    # reciprocal is expanded about x=1/5, not a precomputed rotated field.
    beta = (25 * W**2 + 1) / 80
    d = 8 * beta + scalar(1, 5)
    harvest = scalar(16, 5) * beta + scalar(4, 25)
    X, Z = (z + zb) / 2, (z - zb) / (2 * I)
    x, y = X + scalar(1, 5), 5 * W * Z - X + scalar(4, 5)
    inverse_x = 5 * sum((-5 * X)**j for j in range(max_degree + 1))
    require(truncate(x * inverse_x) == polynomials.one, "truncated reciprocal of x failed")
    f = truncate(x * (1 - x - y))
    g = truncate(y * d - beta * y**2 * inverse_x - harvest)
    original = truncate(f + I * (f + g) / (5 * W))
    require(homogeneous(original, 0) == 0, "original-model equilibrium failed")
    require(homogeneous(original, 1) == I * W * z, "positive rotation failed")
    vector = (original, conjugate(original))
    report = {
        "method": "complex Poincare-Dulac vector-field Lie transforms",
        "arithmetic": "exact Gaussian rational-function field QQ_I(w)",
        "sympy_version": s.__version__,
        "clock": "original physical t; no time rescaling",
        "coordinates": "z=X+i*Z; X=x-1/5; Z=(x+y-1)/(5*w)",
        "input_beta": show(beta),
        "input_d": show(d),
        "input_h": show(harvest),
        "original_complex_field": terms(original),
        "normalization_steps": {},
    }
    for n in range(2, max_degree + 1):
        before = homogeneous(vector[0], n)
        h = polynomials.from_dict({
            (j, k): c / (I * W * (j - k - 1))
            for (j, k), c in before.items() if j - k != 1
        })
        generator = (h, conjugate(h))
        transformed = vector
        current = vector
        factorial = 1
        for power in range(1, (max_degree - 1) // (n - 1) + 1):
            current = bracket(generator, current)
            factorial *= power
            transformed = tuple(
                truncate(transformed[j] + current[j] / factorial)
                for j in range(2)
            )
        vector = transformed
        after = homogeneous(vector[0], n)
        require(all(j - k == 1 for j, k in after), f"degree {n} removal failed")
        require(vector[1] == conjugate(vector[0]), f"degree {n} reality failed")
        require(
            all(j - k == 1 for j, k in vector[0] if 1 <= j + k <= n),
            f"degree {n} transformation disturbed an earlier normal-form degree",
        )
        report["normalization_steps"][str(n)] = {
            "generator_z_component": terms(h),
            "retained_resonances": terms(after),
            "removed_nonresonant_coefficients": len(h),
        }
    report["normal_form_z_component"] = terms(vector[0])
    return w, field, vector[0], report


def main():
    w, field, normal, report = compute()
    c3 = field.to_sympy(normal.get((2, 1), field.zero))
    cubic_real = s.factor(s.re(c3))
    beta = s.Symbol("beta", real=True)
    k4 = s.factor(cubic_real.subs(w**2, (80 * beta - 1) / 25))
    expected_k4 = 25 * beta * (2000 * beta**2 + 160 * beta - 7) / (4 * (80 * beta - 1))
    require(s.cancel(k4 - expected_k4) == 0, "cubic focus formula mismatch")
    roots = s.solve(s.together(k4).as_numer_denom()[0], beta)
    beta0 = next(root for root in roots if root.is_positive and (80 * root - 1).is_positive)
    w0_squared = s.simplify((80 * beta0 - 1) / 25)
    derivative = s.simplify(s.diff(k4, beta).subs(beta, beta0))
    require(s.simplify(beta0 - (s.sqrt(51) - 4) / 100) == 0, "focus parameter mismatch")
    require(
        s.simplify(derivative - s.Rational(17, 2) - 4 * s.sqrt(51)) == 0,
        "cubic transversality mismatch",
    )
    require(w0_squared.is_positive and derivative.is_positive, "focus sign mismatch")
    report["results"] = {
        "cubic_normal_form_complex": s.sstr(s.factor(c3)),
        "cubic_normal_form_real": s.sstr(cubic_real),
        "k4_beta": s.sstr(k4),
        "beta0_derived_from_k4": s.sstr(beta0),
        "w0_squared": s.sstr(w0_squared),
        "k4_prime_beta0": s.sstr(derivative),
    }
    c5 = field.to_sympy(normal.get((3, 2), field.zero))
    quintic_real = s.factor(s.re(c5))
    k6 = s.radsimp(s.simplify(quintic_real.subs(w**2, w0_squared)))
    require(
        s.simplify(k6 + s.Rational(155, 48) + 235 * s.sqrt(51) / 576) == 0,
        "quintic focus formula mismatch",
    )
    require(k6.is_negative, "quintic focus sign mismatch")
    require(
        not any(value.has(s.Float) for value in (c3, c5, k4, beta0, derivative, k6)),
        "floating-point input encountered",
    )
    report["results"].update({
        "quintic_normal_form_complex": s.sstr(s.factor(c5)),
        "quintic_normal_form_real": s.sstr(quintic_real),
        "quintic_real_at_focus": s.sstr(k6),
        "one_turn_angle_clock_quintic_at_focus": s.sstr(2 * s.pi * k6 / s.sqrt(w0_squared)),
        "angular_frequency_at_focus": s.sstr(s.sqrt(w0_squared)),
        "comparisons": "PASS: k4(beta), beta0, k4_prime(beta0), k6(beta0), exact signs",
    })
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
