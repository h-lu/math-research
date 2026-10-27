#!/usr/bin/env python3
"""Exact finite focus calculation for the smooth harvesting core.

Run from the source root: python3 scripts/check-focus.py [--details]
Dependency: SymPy. No floating-point arithmetic or trajectory integration.
The assertions check finite algebra; the analytic return estimates remain
the arguments of the manuscript.
"""

import argparse
import json
import sympy as s


def exact(expr):
    return s.factor(s.cancel(expr))


def coefficients(expr, degree):
    p = s.Poly(s.expand(expr), X, Z)
    return s.Matrix([
        exact(p.coeff_monomial(X ** (degree - j) * Z ** j))
        for j in range(degree + 1)
    ])


def polynomial(vector):
    n = len(vector) - 1
    return sum(vector[j] * X ** (n - j) * Z ** j for j in range(n + 1))


def circle_moment(n, j):
    """Mean of cos(theta)^(n-j) sin(theta)^j, with no integration oracle."""
    if n % 2 or j % 2:
        return s.S.Zero
    a, b = (n - j) // 2, j // 2
    return R(
        s.factorial(2 * a) * s.factorial(2 * b),
        4 ** (a + b) * s.factorial(a) * s.factorial(b) * s.factorial(a + b),
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--details", action="store_true", help="also print all V_n coefficients as JSON")
args = parser.parse_args()
R = s.Rational
X, Z, beta = s.symbols("X Z beta")
w = s.symbols("w", nonzero=True)

# Original smooth model, u=nu=0, in the equilibrium-preserving chart.
beta_w = (25 * w ** 2 + 1) / 80
x, y = X + R(1, 5), 5 * w * Z - X + R(4, 5)
d, harvest = 8 * beta_w + R(1, 5), R(16, 5) * beta_w + R(4, 25)
f = x * (1 - x - y)
g = y * (d - beta_w * y / x) - harvest
phi = -5 * w * X * Z
psi = -(
    80 * w * X ** 2 * Z + (125 * w ** 2 + 5) * X ** 2
    + (6 * w - 250 * w ** 3) * X * Z
    + (125 * w ** 4 + 5 * w ** 2) * Z ** 2
) / (16 * w * (5 * X + 1))
rotation_residuals = [exact(f + w * Z - phi), exact((f + g) / (5 * w) - w * X - psi)]
assert rotation_residuals == [0, 0]
print("rotation-coordinate identities: PASS (2 exact zero residuals)", flush=True)

# Expand directly from the original model. Only nonlinear degrees 2,...,5
# enter dV/dt through degree 6. The finite geometric sum is exact to that order.
numerator = exact(((f + g) / (5 * w) - w * X) * (1 + 5 * X))
psi_series = s.expand(numerator * sum((-5 * X) ** j for j in range(6)))
P = {2: phi}
Q = {n: polynomial(coefficients(psi_series, n)) for n in range(2, 6)}
V = {2: (X ** 2 + Z ** 2) / 2}
K, details = {}, {}

for n in range(3, 7):
    H = coefficients(sum(
        s.diff(V[j], X) * P.get(n - j + 1, 0)
        + s.diff(V[j], Z) * Q.get(n - j + 1, 0)
        for j in range(2, n)
    ), n)
    K[n] = exact(sum(H[j] * circle_moment(n, j) for j in range(n + 1)))
    radial = (X ** 2 + Z ** 2) ** (n // 2) if n % 2 == 0 else 0
    rhs = coefficients(K[n] * radial, n) - H
    operator = s.zeros(n + 1)
    for j in range(n + 1):
        if j:
            operator[j - 1, j] = w * j
        if j < n:
            operator[j + 1, j] = -w * (n - j)
    system, target = operator, rhs
    if n % 2 == 0:
        # Fix the radial kernel by [X^n] V_n = 0, as in the appendix.
        system = system.col_join(s.Matrix([[1] + [0] * n]))
        target = target.col_join(s.Matrix([0]))
    vector, free = system.gauss_jordan_solve(target)
    assert free.rows == 0
    vector = vector.applyfunc(exact)
    residuals = (operator * vector - rhs).applyfunc(exact)
    assert all(value == 0 for value in residuals)
    assert n % 2 or vector[0] == 0
    V[n] = polynomial(vector)
    details[n] = {"k": s.sstr(K[n]), "V_coefficients": [s.sstr(v) for v in vector]}
    print(f"degree {n} homological equation: PASS ({n + 1} exact zero residuals)", flush=True)

assert K[3] == K[5] == 0
k4_beta = exact(K[4].subs(w ** 2, (80 * beta - 1) / 25))
expected_k4 = 25 * beta * (2000 * beta ** 2 + 160 * beta - 7) / (4 * (80 * beta - 1))
assert exact(k4_beta - expected_k4) == 0
root = s.sqrt(51)
beta0 = (root - 4) / 100
w0_squared = (-21 + 4 * root) / 125
k4_at_focus = s.simplify(k4_beta.subs(beta, beta0))
k4_derivative = s.simplify(s.diff(k4_beta, beta).subs(beta, beta0))
k6_at_focus = s.radsimp(s.simplify(K[6].subs(w ** 2, w0_squared)))
assert k4_at_focus == 0
assert s.simplify(k4_derivative - (R(17, 2) + 4 * root)) == 0
assert s.simplify(k6_at_focus + R(155, 48) + 235 * root / 576) == 0
assert w0_squared > 0 and k4_derivative > 0 and k6_at_focus < 0
assert all(not expression.has(s.Float) for expression in [*K.values(), k4_beta, k4_derivative, k6_at_focus])

print("k4(beta) =", k4_beta)
print("k4(beta0) =", k4_at_focus)
print("k4_prime(beta0) =", k4_derivative)
print("k6(beta0) =", k6_at_focus)
print("focus frequency squared > 0, k4_prime > 0, k6 < 0: PASS (exact signs)")
print("all finite-algebra checks: PASS; no floating-point arithmetic")
if args.details:
    print(json.dumps(details, indent=2))
