#!/usr/bin/env python3
"""Independent finite algebra for the quartic and curvature-graph reductions.

Print deterministic JSON; never write files. All arithmetic is exact. The
bundled original-source inputs are checked before the calculations. Ordinary
state derivatives hold the physical controls fixed. See notes/QUARTIC-CURVATURE.md
for the distinction between these identities and the analytic return-map claims.
"""

if not __debug__:
    raise SystemExit("Optimized Python is unsupported; run without -O/-OO.")


from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sp


CHECKS = []


def require(condition, label):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def equal(actual, expected, label):
    require(sp.simplify(actual - expected) == 0, label)


def elementary(values, degree):
    return sum((sp.prod(c) for c in combinations(values, degree)), sp.S.Zero)


def text(expr):
    return sp.sstr(sp.factor(expr))


def source_provenance():
    root = Path(__file__).resolve().parents[1] / "inputs" / "original-source"
    inputs = {
        "sections/03-sharp-cyclicity.tex": (
            "2cabfdc6c56df854cd6ed453f82533935f9dfbe04b305bd6138a2b5faa2eb76f",
            ["eq:J-exact", "eq:V-boundary", "eq:W-boundary",
             "eq:cubic-critical-value", "eq:quartic-Vieta"],
        ),
        "sections/04-saturation-necessity.tex": (
            "bcce62e5d19a36f130bb4df1f69ebf31dc67c8e5da09905b0535426b2fdd2ecb",
            ["sec:main-saturation", "prop:k4-necessary"],
        ),
        "sections/06-actual-bifurcations.tex": (
            "719f4b95ffcd9233b5c5649c4d12f2df09e06d25bf3bf9045f97ae8c8e27ee04",
            ["quartic:centre", "quartic:joint-quotient", "quartic:preparation",
             "quartic:double-surface"],
        ),
        "scripts/check-local-diagrams.py": (
            "053d03b09431fbca135bf9cd45562c1ad41c0846c8935c8d35da854e440f2fa6",
            ["sp.discriminant", "sample_labels"],
        ),
    }
    result = {}
    for name, (digest, anchors) in inputs.items():
        raw = (root / name).read_bytes()
        require(sha256(raw).hexdigest() == digest, "pinned input hash: " + name)
        for anchor in anchors:
            require(anchor in raw.decode("utf-8"), "source anchor: " + anchor)
        result[name] = {"sha256": digest, "anchors": anchors}
    return result


def quartic_checks():
    x, y, t, z = sp.symbols("x y t z", real=True)
    A, S, T = sp.symbols("A S T", real=True)
    C, B = sp.symbols("C B", positive=True)
    alpha, beta, d = sp.symbols("alpha beta d", real=True)
    roots = sp.symbols("r1:5", nonnegative=True)
    e = {j: elementary(roots, j) for j in range(1, 5)}

    # Start with the product of four unspecified roots; solve the physical
    # coefficient matching equations, rather than insert the Vieta controls.
    physical = alpha + beta*x**2 + x**3 - C*x**4 + 3*B*d*x/2
    product = -C*sp.prod(x-r for r in roots)
    difference = sp.Poly(physical-product, x)
    controls = sp.solve([difference.coeff_monomial(x**j) for j in (0, 1, 2)],
                        (alpha, beta, d), dict=True)[0]
    equal(controls[alpha], -C*e[4], "Vieta offset from coefficient matching")
    equal(controls[beta], -C*e[2], "Vieta multiplier from coefficient matching")
    equal(controls[d], 2*C*e[3]/(3*B), "Vieta separation from coefficient matching")
    equal(difference.coeff_monomial(x**3), 1-C*e[1], "Vieta root sum constraint")
    constrained = {roots[-1]: 1/C-sum(roots[:-1])}
    equal((physical.subs(controls)-product).subs(constrained), 0,
          "full physical quartic after root sum constraint")

    phi = A-S*x**2+x**3-x**4+T*x
    normalized = sp.Poly(phi+sp.prod(x-r for r in roots), x)
    normalized_map = sp.solve([normalized.coeff_monomial(x**j) for j in (0, 1, 2)],
                              (A, S, T), dict=True)[0]
    for variable, expected in ((A, -e[4]), (S, e[2]), (T, e[3])):
        equal(normalized_map[variable], expected, "normalized coefficient: " + str(variable))
    equal(normalized.coeff_monomial(x**3), 1-e[1], "normalized root sum")

    # The four equations themselves determine the quadruple configuration.
    H = sp.Matrix([sp.diff(phi, x, j) for j in range(4)])
    solutions = sp.solve(list(H), (A, S, T, x), dict=True)
    require(len(solutions) == 1, "unique algebraic quadruple configuration")
    centre = solutions[0]
    for variable, expected in ((A, -sp.Rational(1, 256)), (S, sp.Rational(3, 8)),
                               (T, sp.Rational(1, 16)), (x, sp.Rational(1, 4))):
        equal(centre[variable], expected, "quadruple centre: " + str(variable))
    unfolding = H.jacobian((A, S, T, x)).subs(centre)
    equal(unfolding.det(), -48, "quadruple four-equation determinant")
    equal(H[:3, :].jacobian((A, S, T)).subs(centre).det(), 2,
          "three-observation unfolding determinant")
    equal(sp.diff(phi, x, 4), -24, "quartic fourth derivative")

    # Depression is derived by removing the cubic coefficient of -Phi.
    shift = sp.solve(sp.Poly(-phi.subs(x, y+t), y).coeff_monomial(y**3), t)[0]
    depressed = sp.Poly(-phi.subs(x, y+shift), y)
    u, v, w = sp.symbols("u v w", real=True)
    coordinate_map = {u: depressed.coeff_monomial(y**2),
                      v: depressed.coeff_monomial(y), w: depressed.coeff_monomial(1)}
    equal(shift, sp.Rational(1, 4), "depressing translation")
    for variable, expected in ((u, S-sp.Rational(3, 8)),
                               (v, S/2-T-sp.Rational(1, 8)),
                               (w, -A+S/16-T/4-sp.Rational(3, 256))):
        equal(coordinate_map[variable], expected, "exact limiting preparation map: " + str(variable))
    coordinate_det = sp.Matrix(list(coordinate_map.values())).jacobian((A, S, T)).det()
    equal(coordinate_det, 1, "prepared control determinant")

    # Derive strata by solving the defining derivative equations.
    p = y**4+u*y**2+v*y+w
    fold = sp.solve([p.subs(y, t), sp.diff(p, y).subs(y, t)], (v, w), dict=True)[0]
    equal(fold[v], -4*t**3-2*u*t, "derived double-root surface v")
    equal(fold[w], 3*t**4+u*t**2, "derived double-root surface w")
    equal(sp.discriminant(p, y).subs(fold), 0, "discriminant vanishes on double surface")
    triple_u = sp.solve(sp.diff(p, y, 2).subs(fold).subs(y, t), u)[0]
    triple = {u: triple_u, **{key: value.subs(u, triple_u) for key, value in fold.items()}}
    equal(p.subs(triple), (y-t)**3*(y+3*t), "derived triple-root factorization")
    equal(sp.diff(p, y, 3).subs(triple).subs(y, t), 24*t,
          "triple curve becomes quadruple only at zero")
    s = sp.symbols("s", real=True)
    two_double_product = sp.Poly((y-t)**2*(y-s)**2, y)
    opposite = sp.solve(two_double_product.coeff_monomial(y**3), s)[0]
    equal(opposite, -t, "two double roots are opposite after depression")
    two_double = sp.Poly(two_double_product.as_expr().subs(s, opposite), y)
    for power, expected in ((2, -2*t**2), (1, 0), (0, t**4)):
        equal(two_double.coeff_monomial(y**power), expected,
              "derived two-double coefficient degree " + str(power))

    # The repeated-root part of the K4 boundary has a nonnegative root t. The remaining
    # two roots have sum 1-2t and product z; expand this independent factor.
    boundary = sp.Poly(-(x-t)**2*(x**2-(1-2*t)*x+z), x)
    boundary_map = {A: boundary.coeff_monomial(1), S: -boundary.coeff_monomial(x**2),
                    T: boundary.coeff_monomial(x)}
    for variable, expected in ((A, -t**2*z), (S, z+2*t-3*t**2),
                               (T, 2*t*z+t**2-2*t**3)):
        equal(boundary_map[variable], expected, "K4 repeated-root boundary: " + str(variable))
    equal(phi.subs(boundary_map), boundary.as_expr(), "K4 repeated-root full factorization")
    remaining_discriminant = sp.discriminant(x**2-(1-2*t)*x+z, x)
    equal(remaining_discriminant, (1-2*t)**2-4*z, "remaining-root admissibility discriminant")
    k4_triple = {A: -t**3*(1-3*t), S: 3*t-6*t**2, T: 3*t**2-8*t**3}
    equal(phi.subs(k4_triple), -(x-t)**3*(x-(1-3*t)), "K4 triple-root curve")
    partner = sp.Rational(1, 2)-t
    k4_two_double = {A: -t**2*partner**2, S: sp.Rational(1, 4)+t-2*t**2,
                     T: t*partner}
    equal(phi.subs(k4_two_double), -(x-t)**2*(x-partner)**2, "K4 two-double curve")

    # Critical-value difference is an integral of the factored derivative.
    # This recovers the switch between the two lower offset arcs.
    cu, cv, cw = sp.symbols("cu cv cw", real=True)
    derivative = -4*(x-cu)*(x-cv)*(x-cw)
    value_difference = sp.integrate(derivative, (x, cu, cw))
    equal(value_difference, (cw-cu)**3*(cu-2*cv+cw)/3,
          "equal-maxima critical-value identity")
    middle = sp.solve([cu+cv+cw-sp.Rational(3, 4), cu-2*cv+cw], (cv, cw), dict=True)[0][cv]
    equal(middle, sp.Rational(1, 4), "equal maxima middle critical point")
    equal(sp.diff(phi, x).subs(x, middle), T-S/2+sp.Rational(1, 8),
          "equal-maxima coefficient plane")

    # Explicit rational interior configuration; no floating-point root finding.
    sample_roots = [sp.Rational(j, 10) for j in (1, 2, 3, 4)]
    sample = {A: -elementary(sample_roots, 4), S: elementary(sample_roots, 2),
              T: elementary(sample_roots, 3)}
    equal(phi.subs(sample), -sp.prod(x-r for r in sample_roots), "rational K4 interior example")
    sample_negative = -(y**2+sample[S]*y-sample[A])
    sample_companions = sp.solve(sample_negative, y)
    require(len(sample_companions) == 2 and all(r.is_negative is True for r in sample_companions),
            "rational interior example has two distinct negative companions")
    negative_at_centre = (A-S*y-y**2).subs(centre)
    companion_centres = sp.solve(negative_at_centre, y)
    for actual, expected in zip(companion_centres,
                                [(-3-2*sp.sqrt(2))/16, (-3+2*sp.sqrt(2))/16]):
        equal(actual, expected, "quadruple negative companion centre")

    return {
        "physical_vieta_controls": {"alpha4": "-C*e4", "beta4": "-C*e2",
                                    "d4": "2*C*e3/(3*B)", "root_sum": "1/C"},
        "normalized_map": {"A": "-e4", "S": "e2", "T": "e3", "root_sum": "1"},
        "quadruple_centre": {str(key): text(value) for key, value in centre.items()},
        "quadruple_matrix_columns": ["A", "S", "T", "x"],
        "quadruple_matrix": [[text(vv) for vv in row] for row in unfolding.tolist()],
        "quadruple_determinant": text(unfolding.det()),
        "prepared_coordinates": {str(key): text(value) for key, value in coordinate_map.items()},
        "prepared_coordinate_determinant": text(coordinate_det),
        "depressed_discriminant": text(sp.discriminant(p, y)),
        "double_surface": {str(key): text(value) for key, value in fold.items()},
        "triple_curve": {str(key): text(value) for key, value in triple.items()},
        "two_double_curve": {"u": "-2*t**2", "v": "0", "w": "t**4"},
        "K4_repeated_root_boundary": {str(key): text(value) for key, value in boundary_map.items()},
        "K4_repeated_root_admissibility": "0 <= t <= 1/2; 0 <= z <= (1-2*t)**2/4",
        "K4_triple_curve": {str(key): text(value) for key, value in k4_triple.items()},
        "K4_triple_root_range": "0 <= t <= 1/3",
        "K4_two_double_curve": {str(key): text(value) for key, value in k4_two_double.items()},
        "K4_two_double_root_range": "0 <= t <= 1/2",
        "equal_maxima_value_difference": text(value_difference),
        "equal_maxima_plane": "T = S/2 - 1/8; three distinct positive critical points assumed",
        "rational_interior_example": {
            "positive_roots": [text(rr) for rr in sample_roots],
            "coefficients": {str(key): text(value) for key, value in sample.items()},
            "negative_companions": [text(rr) for rr in sample_companions],
        },
        "quadruple_negative_companions": [text(rr) for rr in companion_centres],
    }


def simplex_checks():
    roots = sp.symbols("r1:5", nonnegative=True)
    e = {j: elementary(roots, j) for j in range(1, 5)}
    pair_products = [a*b for a, b in combinations(roots, 2)]
    equal(sp.prod(pair_products), e[4]**3, "six-pair-product AM-GM witness")
    equal(3*e[1]**2-8*e[2], sum((a-b)**2 for a, b in combinations(roots, 2)),
          "S upper-bound sum-of-squares witness")
    # At any face maximum with at least j positive coordinates, these
    # derivative differences force all positive coordinates to agree.
    candidates, maxima = {}, {}
    for degree in (2, 3, 4):
        for i, j in combinations(range(4), 2):
            others = [roots[k] for k in range(4) if k not in (i, j)]
            equal(sp.diff(e[degree], roots[i])-sp.diff(e[degree], roots[j]),
                  (roots[j]-roots[i])*elementary(others, degree-2),
                  "simplex stationary difference e%d pair %d%d" % (degree, i+1, j+1))
        values = {n: sp.binomial(n, degree)/sp.Integer(n)**degree if n >= degree else sp.S.Zero
                  for n in range(1, 5)}
        maxima[degree] = max(values.values())
        candidates["e"+str(degree)] = {str(n): text(value) for n, value in values.items()}
    for degree, expected in ((2, sp.Rational(3, 8)), (3, sp.Rational(1, 16)),
                             (4, sp.Rational(1, 256))):
        equal(maxima[degree], expected, "simplex exact maximum e"+str(degree))
    a, b, c, d = roots
    invariant_sos = ((a-b)**2*(c-d)**2+(a-c)**2*(b-d)**2+(a-d)**2*(b-c)**2)/2
    equal(e[2]**2-3*e[1]*e[3]+12*e[4], invariant_sos,
          "additional K4 invariant sum-of-squares")
    t, x, S = sp.symbols("t x S", real=True)
    zero_face_roots = [t, 1-t, 0, 0]
    equal(elementary(zero_face_roots, 3), 0, "T=0 face root representation")
    equal(elementary(zero_face_roots, 4), 0, "T=0 face offset")
    equal(elementary(zero_face_roots, 2), t*(1-t), "T=0 face S parameter")
    equal(sp.Rational(1, 4)-t*(1-t), (t-sp.Rational(1, 2))**2,
          "T=0 face S upper bound")
    equal(-x**4+x**3-S*x**2, x**2*(-x**2+x-S), "T=A=0 face factorization")
    return {
        "assumptions": "r_i >= 0; sum(r_i) = 1",
        "face_stationary_values_by_support_size": candidates,
        "extrema": {"A_min": "-1/256", "A_max": "0", "S_max": "3/8", "T_max": "1/16"},
        "maximizing_root_vector": ["1/4"]*4,
        "T_zero_face": "A=0; 0<=S<=1/4; roots=(t,1-t,0,0), 0<=t<=1",
        "pair_product_identity": "product_{i<j}(r_i*r_j) = e4**3",
        "AM_GM_consequence": "e2 >= 6*sqrt(e4); e2**2-4*e4 >= 32*e4 > 0 when all r_i>0",
        "additional_consistency_inequality": "S**2-3*T-12*A >= 0",
        "additional_inequality_status": "exact SOS consequence, not a claimed characterization of K4",
        "inequality_reasoning": "nonnegative roots, compact simplex, face stationarity, and six-term AM-GM; see note",
    }


def curvature_checks():
    u, q, delta, r, v, m = sp.symbols("u q delta r v m", real=True)
    A, k, c = sp.symbols("A k c", nonzero=True, real=True)
    A1, A2, ell1, ell2, H1, H2 = sp.symbols("A1 A2 ell1 ell2 H1 H2", real=True)
    Rp, W = sp.symbols("Rp W", positive=True)
    Rpp = sp.symbols("Rpp", real=True)

    # Differentiate the signed-root definition, not the manuscript's
    # already differentiated formulas.
    profile = sp.Function("T")(u)
    g = u*profile
    ell = profile**2+u*sp.diff(profile**2, u)/2
    H = sp.diff(ell, u)/2
    equal(sp.diff(g**2, u)/(2*u), ell, "first passage derivative from signed-root identity")
    equal(sp.diff(ell, u)/(2*u), H/u, "second passage derivative from signed-root identity")
    b, d = sp.symbols("b d", real=True)
    jet = 1+A1*u/2+b*u**2+d*u**3
    ell_jet = sp.expand(jet**2+u*sp.diff(jet**2, u)/2)
    H_jet = sp.diff(ell_jet, u)/2
    equal(ell_jet.subs(u, 0), 1, "passage derivative boundary")
    equal(H_jet.subs(u, 0), 3*A1/4, "passage curvature boundary from strength coefficient")

    G1, G2, R = sp.Function("Gamma1"), sp.Function("Gamma2"), sp.Function("R")
    z = sp.symbols("z", real=True)
    downstream = G1(q)-delta
    Q = delta+G2(downstream)
    G2p = sp.Subs(sp.diff(G2(z), z), z, downstream)
    G2pp = sp.Subs(sp.diff(G2(z), z, 2), z, downstream)
    Qp = G2p*sp.diff(G1(q), q)
    Qpp = G2pp*sp.diff(G1(q), q)**2+G2p*sp.diff(G1(q), q, 2)
    equal(sp.diff(Q, q), Qp.doit(), "actual nested-passage first chain rule at fixed delta")
    equal(sp.diff(Q, q, 2), Qpp.doit(), "actual nested-passage second chain rule at fixed delta")
    Ppp = (sp.Subs(sp.diff(R(z), z, 2), z, Q)*sp.diff(Q, q)**2
           +sp.Subs(sp.diff(R(z), z), z, Q)*sp.diff(Q, q, 2))
    equal(sp.diff(R(Q), q, 2), Ppp.doit(), "closing-leg curvature chain rule at fixed delta")
    root_Qp = ell2*ell1
    root_Qpp = ell2*H1/r+H2*ell1**2/(r*sp.sqrt(W))
    J = sp.expand(r*(Rpp*root_Qp**2+Rp*root_Qpp))
    equal(J, r*Rpp*(ell2*ell1)**2+Rp*(ell2*H1+H2*ell1**2/sp.sqrt(W)),
          "root curvature exact formula after fixed-control differentiation")
    J0 = J.subs({r: 0, ell1: 1, ell2: 1, H1: 3*A1/4, H2: 3*A2/4,
                 Rp: k, W: 1-v})
    equal(J0, 3*k*(A1+A2/sp.sqrt(1-v))/4, "root curvature boundary")
    equal(sp.diff(J0, v).subs({v: 0, A1: A, A2: -A}), -3*k*A/8,
          "root curvature v transversality sign and coefficient")

    # Retain both strengths' independent m dependence, and that of R'.
    a1, a2, km = sp.symbols("a1 a2 km", real=True)
    strength1 = A+a1*m+a2*m**2
    strength2 = m-strength1
    boundary_with_m = J0.subs({A1: strength1, A2: strength2, k: k+km*m})
    Jm = sp.diff(boundary_with_m, m).subs({m: 0, v: 0})
    Jv = sp.diff(boundary_with_m, v).subs({m: 0, v: 0})
    equal(Jm, 3*k/4, "root curvature independent-strength derivative")
    # Solve first in z=sqrt(1-v)>0, preserving the positive-root condition.
    boundary_z = sp.solve(A1*z+A2, z)[0]
    Vboundary = 1-boundary_z**2
    equal(Vboundary, 1-(-A2/A1)**2, "curvature-zero boundary graph")
    Vm_boundary = sp.diff(Vboundary.subs({A1: strength1, A2: strength2}), m).subs(m, 0)
    equal(Vm_boundary, 2/A, "boundary graph strength derivative with moving strengths")

    # Aligned finite cancellation Q=q+q^(5/2)L(sqrt(q)) supplies Jr(0).
    # Only the indicated finite jets can contribute to that coefficient.
    l0, l1, a0, r3 = sp.symbols("l0 l1 a0 r3", real=True)
    Qaligned = r**2+r**5*(l0+l1*r)
    Paligned = a0+k*Qaligned+c*Qaligned**2+r3*Qaligned**3
    state_derivative = lambda expression: sp.diff(expression, r)/(2*r)
    Jaligned = sp.cancel(r*state_derivative(state_derivative(Paligned)))
    equal(Jaligned.subs(r, 0), 0, "aligned curvature boundary cancellation")
    Jr = sp.diff(Jaligned, r).subs(r, 0)
    equal(Jr, 2*c, "aligned curvature linear coefficient")
    Vm, Vr = sp.simplify(-Jm/Jv), sp.simplify(-Jr/Jv)
    equal(Vm, 2/A, "implicit graph Vm from independent curvature derivatives")
    equal(Vr, 16*c/(3*k*A), "implicit graph Vr from independent curvature derivatives")
    Delta_leading = r**2*(Vm*m+Vr*r)
    turn_factor = sp.cancel(sp.diff(Delta_leading, r)/r)
    equal(sp.diff(turn_factor, m), 4/A, "turn factor strength derivative")
    equal(sp.diff(turn_factor, r), 16*c/(k*A), "turn factor root derivative")
    rho = sp.solve(turn_factor, r)[0]
    gamma = sp.cancel(rho/m)
    kappa = sp.cancel(Delta_leading.subs(r, rho)/m**3)
    equal(gamma, -k/(4*c), "signed turn location leading coefficient")
    equal(kappa, k**2/(24*A*c**2), "cubic separation critical-value coefficient")

    x, C, B = sp.symbols("x C B", positive=True)
    normalized_graph = sp.cancel((3*B*C**2/(2*m**3)*Delta_leading)
                                 .subs({r: m*x/C, c: -C, k: 1, A: B}))
    equal(normalized_graph, 3*x**2-8*x**3, "physical graph main-scale normalization")
    # Derive the normalized state-derivative operator from q=m^2*x^2/C^2.
    F = sp.Function("F")(x)
    D = m**4*F/C**3
    dqdx = sp.diff(m**2*x**2/C**2, x)
    Dq = sp.diff(D, x)/dqdx
    Dqq = sp.diff(Dq, x)/dqdx
    equal(Dq, m**2*sp.diff(F, x)/(2*C*x), "main-scale first physical state derivative")
    equal(Dqq, C*(x*sp.diff(F, x, 2)-sp.diff(F, x))/(4*x**3),
          "main-scale second physical state derivative")
    offset, slope, separation = sp.symbols("A4 S4 T4", real=True)
    phi = offset-slope*x**2+x**3-x**4+separation*x
    curvature_numerator = sp.expand(x*sp.diff(phi, x, 2)-sp.diff(phi, x))
    equal(curvature_numerator, normalized_graph-separation, "independent quartic and graph curvature agreement")
    positive_turns = [rr for rr in sp.solve(sp.diff(normalized_graph, x), x) if rr.is_positive]
    require(positive_turns == [sp.Rational(1, 4)], "unique positive turn of limiting cubic")
    turn = positive_turns[0]
    equal(normalized_graph.subs(x, turn), sp.Rational(1, 16), "limiting turning value")
    equal(sp.diff(normalized_graph, x, 2).subs(x, turn), -6, "limiting turn nondegeneracy")
    equal(normalized_graph.subs(x, sp.Rational(3, 8)), 0, "zero-separation nonzero curvature root")
    return {
        "normalization": "q=r**2; delta=r**2*v only AFTER physical state differentiation",
        "passage_boundary": {"ell(0)": "1", "H(0)": "3*A_nu/4"},
        "J_exact": text(J),
        "J_boundary": text(J0),
        "J_first_derivatives_at_cancellation": {"r": text(Jr), "m": text(Jm), "v": text(Jv)},
        "V_boundary": text(Vboundary), "V_m": text(Vm), "V_r": text(Vr),
        "turn_factor_derivatives": {"m": text(sp.diff(turn_factor, m)),
                                    "r": text(sp.diff(turn_factor, r))},
        "rho_over_m_at_zero": text(gamma), "delta_c_over_m_cubed_at_zero": text(kappa),
        "reference_specialization": "k=1, c=-C, A=B; B>0, C>0",
        "normalized_curvature_graph": text(normalized_graph),
        "quartic_curvature_numerator": text(curvature_numerator),
        "positive_turn": text(turn), "turning_value": text(normalized_graph.subs(x, turn)),
        "zero_separation_nonzero_root": "3/8",
    }


def main():
    provenance = source_provenance()
    quartic = quartic_checks()
    simplex = simplex_checks()
    curvature = curvature_checks()
    print(json.dumps({
        "status": "PASS", "arithmetic": "exact symbolic and rational; no numerical approximation",
        "python_version": ".".join(map(str, sys.version_info[:3])), "sympy_version": sp.__version__,
        "source_commit": "31515b525a12ff95a1e84875e42ad9b7fe28d76c",
        "original_source_inputs": provenance,
        "existing_check_inventory": {
            "script": "scripts/check-local-diagrams.py",
            "scope": "exact prepared-quartic strata/discriminant and five exact rational root-count labels",
            "extension": "coefficient reconstruction, K4 witnesses, quadruple unfolding, and curvature normalization",
        },
        "quartic": quartic, "simplex": simplex, "curvature": curvature,
        "checks_passed": len(CHECKS), "checks": CHECKS,
        "not_certified": [
            "actual first-return domain, positivity margins and field/model admissibility",
            "uniform common boxes, event-layer estimates and all-parameter remainder bounds",
            "implicit-function existence/uniqueness, Hadamard division and smooth preparation hypotheses",
            "weighted confinement, actual six-cycle existence or all-rate counting",
            "finite-scale saturation boundary corrections and Hausdorff convergence",
        ],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
