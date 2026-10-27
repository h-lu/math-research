#!/usr/bin/env python3
"""Exact finite passage and control-rank algebra; print deterministic JSON.

Usage: python3 archival-verification/scripts/check_passage_rank.py
Only the standard library and SymPy are used. This is a finite symbolic
calculation, not certification of common domains or an actual orbit.
"""

if not __debug__:
    raise SystemExit("Optimized Python is unsupported; run without -O/-OO.")

import json
import platform

import sympy as sp


SOURCE_COMMIT = "31515b525a12ff95a1e84875e42ad9b7fe28d76c"
checks = []


def equal(name, actual, expected):
    residual = sp.factor(sp.cancel(sp.expand(actual - expected)))
    if residual != 0:
        raise RuntimeError(f"{name}: nonzero exact residual {residual}")
    checks.append(name)


def poly(expr, variable, degree):
    return sp.series(expr, variable, 0, degree + 1).removeO().expand()


def string(expr):
    return sp.sstr(sp.factor(expr))


# Start with dz/dt = (z-t^2)*beta(t,z)/(1+(z-t^2)*alpha(t,z)).
# Put t=r*s, z=r^2*w(s,r), entry s=-1, w(-1,r)=1.
# Taylor coefficients b_tt and b_z below multiply t^2 and z directly;
# b_tt is half the usual second t derivative.
r, s, t = sp.symbols("r s t")
b, b_t, b_z, b_tt, a = sp.symbols("b b_t b_z b_tt a")
w = sp.Integer(1)
terms = []
for degree in range(1, 4):
    beta_jet = b + b_t*r*s + r**2*(b_z*w + b_tt*s**2)
    # Terms omitted from alpha, beta and the inverse denominator enter
    # dw/ds only at r^4 or above, so cannot change the computed coefficients.
    rhs = r*(w-s**2)*beta_jet*(1-a*r**2*(w-s**2))
    forcing = poly(rhs, r, degree).coeff(r, degree)
    term = sp.integrate(forcing, (s, -1, s)).expand()
    terms.append(term)
    w += r**degree*term
    equal(f"entry condition w_{degree}(-1)=0", term.subs(s, -1), 0)

rhs_exact_jet = (r*(w-s**2)*(b+b_t*r*s+r**2*(b_z*w+b_tt*s**2))
                 / (1+a*r**2*(w-s**2)))
for degree in range(4):
    equal(f"scaled ODE coefficient r^{degree}",
          poly(sp.diff(w, s)-rhs_exact_jet, r, 3).coeff(r, degree), 0)

exit_coefficients = []
for degree in range(1, 4):
    unknown = sp.Symbol(f"exit_{degree}")
    exit_s = (1+sum(exit_coefficients[j-1]*r**j
                  for j in range(1, degree))+unknown*r**degree)
    equation = poly(w.subs(s, exit_s)-exit_s**2, r, degree).coeff(r, degree)
    exit_coefficients.append(sp.factor(sp.solve(equation, unknown)[0]))

e1, e2, e3 = exit_coefficients
equal("ODE-derived root quadratic coefficient", e1, sp.Rational(2, 3)*b)
equal("ODE-derived root cubic coefficient", e2, e1**2)
equal("ODE-derived root quartic coefficient", e3,
      -sp.Rational(8, 15)*a*b+sp.Rational(44, 135)*b**3
      +sp.Rational(2, 15)*b*b_t+sp.Rational(2, 15)*b_tt
      +sp.Rational(2, 3)*b_z)
exit_s = 1+e1*r+e2*r**2+e3*r**3
for degree in range(4):
    equal(f"moving-exit equation coefficient r^{degree}",
          poly(w.subs(s, exit_s)-exit_s**2, r, 3).coeff(r, degree), 0)
g = r*exit_s
gamma = poly(g**2, r, 5)
C1 = gamma.coeff(r, 3)
equal("strength from the ODE", C1, sp.Rational(4, 3)*b)
equal("forced ordinary quadratic passage coefficient",
      gamma.coeff(r, 4), sp.Rational(3, 4)*C1**2)

# Independently solve the formal inverse and compose two such ODE-derived
# jets. The free quartic root coefficients retain all higher local data.
A, B, C, E = sp.symbols("A B C E")
g1 = t+A*t**2+A**2*t**3+C*t**4
g2 = t+B*t**2+B**2*t**3+E*t**4
inverse_coefficients = []
for degree in range(2, 5):
    unknown = sp.Symbol(f"inverse_{degree}")
    candidate = (t+sum(inverse_coefficients[j-2]*t**j
                       for j in range(2, degree))+unknown*t**degree)
    equation = poly(g1.subs(t, candidate)-t, t, degree).coeff(t, degree)
    inverse_coefficients.append(sp.factor(sp.solve(equation, unknown)[0]))
inverse = t+sum(inverse_coefficients[j-2]*t**j for j in range(2, 5))
equal("inverse signed companion agrees with reversal through degree four",
      inverse, -g1.subs(t, -t))
inverse_physical = poly(inverse**2, t, 5)
equal("activation inverse leading correction", inverse_physical.coeff(t, 3), -2*A)
equal("activation inverse quadratic correction", inverse_physical.coeff(t, 4), 3*A**2)
composition = poly(g2.subs(t, g1), t, 4)
equal("composed root quadratic coefficient", composition.coeff(t, 2), A+B)
equal("composed root cubic coefficient", composition.coeff(t, 3), (A+B)**2)
composed_physical = poly(composition**2, t, 5)
equal("aligned strength sum", composed_physical.coeff(t, 3), 2*(A+B))
equal("aligned forced quadratic coefficient", composed_physical.coeff(t, 4), 3*(A+B)**2)
equal("opposite strengths leave a free fifth-half coefficient",
      composed_physical.coeff(t, 5).subs(B, -A), 2*(C+E))

kappa = sp.Symbol("kappa")
first = {a: 0, b_t: 0, b_z: 0, b_tt: kappa}
second = {a: 0, b_t: 0, b_z: 0, b_tt: 0, b: -b}
generic_remainder = sp.factor(2*(e3.subs(first)+e3.subs(second)))
equal("explicit opposite-strength noninverse example", generic_remainder,
      sp.Rational(4, 15)*kappa)
z = sp.Symbol("z")
f_constant = b*(z-t**2)
equal("special constant-beta exact time-reversal symmetry",
      f_constant.subs(b, -b), -f_constant.subs(t, -t))

# The supported control construction: integrate dz/dtau=eta0+eta1*z.
# tau is accumulated cutoff time, and its endpoint is one.
eta0, eta1, tau, z0 = sp.symbols("eta0 eta1 tau z0")
flow = sp.exp(eta1*tau)*z0+eta0*(sp.exp(eta1*tau)-1)/eta1
equal("closing-box affine ODE", sp.diff(flow, tau), eta0+eta1*flow)
equal("closing-box entrance", flow.subs(tau, 0), z0)
phi = (sp.exp(eta1)-1)/eta1
equal("closing-box phi(0)", sp.limit(phi, eta1, 0), 1)
c0, beta_q, k, d_control, beta0 = sp.symbols("c0 beta_q k d_control beta0")
q_contact, e0, clock0 = sp.symbols("q_contact e0 clock0")
beta_contact = beta0+beta_q*q_contact
h0 = e0*(q_contact-t**2)
hd = h0-d_control*e0
base_inactive = sp.Matrix([clock0, 0])
base_difference = sp.Matrix([clock0*a/e0, clock0*beta_contact/e0])
base_active = base_inactive+h0*base_difference
controlled_active = (base_active+(hd-h0)*base_difference
                     +sp.Matrix([0, k*hd*clock0/e0]))
controlled_difference = controlled_active-base_inactive
equal("controlled seam matching in original time t component",
      controlled_difference[0], hd*clock0*a/e0)
equal("controlled seam matching in original time q component",
      controlled_difference[1], hd*clock0*(beta_contact+k)/e0)
equal("controlled normalized active transverse component",
      (controlled_difference[1]/clock0).subs(q_contact, z+d_control),
      (z-t**2)*(beta0+beta_q*(z+d_control)+k))
controlled_beta_at_contact = (beta_contact+k).subs(q_contact, d_control)
H0 = eta0*(1+eta1/2+eta1**2/6)
R0 = lambda x: x+c0*x**2
a_observed = R0(H0)
b_observed = (1+2*c0*H0)*sp.exp(eta1)-1
m_observed = C1.subs(b, controlled_beta_at_contact)-sp.Rational(4, 3)*beta0
local_observations = sp.Matrix([a_observed, b_observed, m_observed, d_control])
local_controls = [eta0, eta1, k, d_control]
local_jacobian = local_observations.jacobian(local_controls).subs(
    dict.fromkeys(local_controls, 0))
expected_local = sp.Matrix([[1, 0, 0, 0], [2*c0, 1, 0, 0],
                            [0, 0, sp.Rational(4, 3), sp.Rational(4, 3)*beta_q],
                            [0, 0, 0, 1]])
for row in range(4):
    for column in range(4):
        equal(f"local control Jacobian ({row+1},{column+1})",
              local_jacobian[row, column], expected_local[row, column])
equal("local control determinant", local_jacobian.det(), sp.Rational(4, 3))

# Reconstruct middle/exterior branches and the physical chart, rather than
# defining an observation matrix as if its unknown sensitivities were known.
x, y, beta, nu, u, v, L, U = sp.symbols("x y beta nu u v L U")
d_physical, h_physical = sp.symbols("d_physical h_physical")
Fell = x*(1-x-y)
Fc = Fell-u*(x-L)
Fr = Fc-v*(x-U)
G = y*(d_physical-beta*y/x)-h_physical
equal("lower branch matched difference", Fell-Fc, (L-x)*(-u))
equal("upper branch matched difference", Fr-Fc, (x-U)*(-v))
equal("lower traces agree", (Fell-Fc).subs(x, L), 0)
equal("upper traces agree", (Fr-Fc).subs(x, U), 0)
ye = sp.Rational(4, 5)-u*(1-5*L)
d_chart = 8*beta+sp.Rational(1, 5)+nu
h_chart = ye*(d_chart-5*beta*ye)
equilibrium = {x: sp.Rational(1, 5), y: ye,
               d_physical: d_chart, h_physical: h_chart}
equal("perturbed chart F equilibrium", Fc.subs(equilibrium), 0)
equal("perturbed chart G equilibrium", G.subs(equilibrium), 0)
equal("u=0 core harvest chart", h_chart.subs(u, 0),
      sp.Rational(16, 5)*beta+sp.Rational(4, 25)+sp.Rational(4, 5)*nu)
J = sp.Matrix([Fc, G]).jacobian([x, y]).subs(equilibrium)
equal("u=0 core linear trace", J.subs(u, 0).trace(), nu)
equal("u=0 core linear determinant", J.subs(u, 0).det(),
      (80*beta-1-5*nu)/25)
frequency_squared = (80*beta-1-5*nu)/25-nu**2/4
linear_rate_one_turn = sp.pi*nu/sp.sqrt(frequency_squared)
linear_rate_nu_at_hopf = sp.diff(linear_rate_one_turn, nu).subs(nu, 0)
equal("canonical linear one-turn radial rate derivative",
      linear_rate_nu_at_hopf, 5*sp.pi/sp.sqrt(80*beta-1))
physical_chart = sp.Matrix([beta, d_chart, h_chart, u, v, L, U])
chart_jacobian = physical_chart.jacobian([nu, L, v, U])
embedding_minor = chart_jacobian.extract([1, 5, 4, 6], [0, 1, 2, 3])
equal("physical chart embedding minor (d,L,v,U)", embedding_minor.det(), 1)

# The contact-clock conversion is a chain rule identity. Strict signs and
# the existence of I are hypotheses; neither is supplied by the symbols.
sigma, sx, Bx, omega, Fcy, Gp = sp.symbols("sigma sx Bx omega Fcy Gp")
e_contact = sigma/sx
omega_squared = -sx*Fcy*Gp/2
equal("contact clock second-derivative identity",
      -2*e_contact*omega_squared, sigma*Fcy*Gp)
normalized_beta = sp.cancel(e_contact*sx*Bx/omega)
equal("physical-to-clock transverse coefficient", normalized_beta, sigma*Bx/omega)
omega1, omega2 = sp.symbols("omega1 omega2", positive=True)
strength1 = C1.subs(b, normalized_beta.subs({sigma: -1, Bx: -u, omega: omega1}))
strength2 = C1.subs(b, normalized_beta.subs({sigma: 1, Bx: -v, omega: omega2}))
K2 = sp.Rational(4, 3)/omega2
equal("lower effort strength", strength1, sp.Rational(4, 3)*u/omega1)
equal("upper effort strength", strength2, -K2*v)
effort_ratio = sp.solve(strength1+strength2, v)[0]
equal("positive effort cancellation ratio", effort_ratio, u*omega2/omega1)
equal("upper effort sensitivity", sp.diff(strength1+strength2, v), -K2)
upper_y = sp.solve(Fc.subs(x, U), y)[0]
s1, s2, Iy = sp.symbols("s1 s2 Iy")
delta_U = (s2+Iy*sp.diff(upper_y, U)).subs(Iy, 0)
equal("upper contact chain-rule sensitivity", delta_U, s2)
lower_y_u0 = sp.solve(Fc.subs({u: 0, x: L}), y)[0]
zeta_L = (s1+Iy*sp.diff(lower_y_u0, L)).subs(Iy, 0)
equal("lower contact chain-rule sensitivity", zeta_L, s1)

# Moving-coordinate cancellation at an actual smooth fold. P is a fixed
# return germ, T_l(x)=k(l)*(x-zeta(l))+mu(l)*(x-zeta(l))^2.
# No finite jet is used here to assert an actual fold exists.
l, zeta_prime, scale_prime, mu0, mu_prime, cubic = sp.symbols(
    "l zeta_prime scale_prime mu0 mu_prime cubic")
q = sp.Symbol("q")
P = q+c0*q**2+cubic*q**3
center = zeta_prime*l
scale = 1+scale_prime*l
mu = mu0+mu_prime*l
Pcenter = P.subs(q, center)
T_of_Pcenter = scale*(Pcenter-center)+mu*(Pcenter-center)**2
R_derivative = ((scale+2*mu*(Pcenter-center))*sp.diff(P, q).subs(q, center)/scale)
equal("moving lower coordinate offset sensitivity", sp.diff(T_of_Pcenter, l).subs(l, 0), 0)
equal("moving lower coordinate multiplier sensitivity",
      sp.diff(R_derivative, l).subs(l, 0), 2*c0*zeta_prime)
a_nu, a_L, b_nu, b_L = sp.symbols("a_nu a_L b_nu b_L")
m_nu, m_L, m_U, delta_nu, delta_L = sp.symbols("m_nu m_L m_U delta_nu delta_L")
natural = sp.Matrix([[a_nu, a_L, 0, 0], [b_nu, b_L, 0, 0],
                     [m_nu, m_L, -K2, m_U], [delta_nu, delta_L, 0, delta_U]])
natural_determinant = sp.factor(natural.det())
equal("natural structural observation determinant", natural_determinant,
      -K2*delta_U*(a_nu*b_L-a_L*b_nu))
fold_determinant = natural_determinant.subs({a_L: 0, b_L: 2*c0*zeta_L})
equal("u=0 natural observation determinant", fold_determinant,
      -2*K2*delta_U*a_nu*c0*zeta_L)

# The retained radial flow is used only for its finite fold algebra.
rho, radial, eta_one, eta_two, ell = sp.symbols("rho radial eta_one eta_two ell")
drift = eta_one*radial+eta_two*radial**3+ell*radial**5
fold_equations = sp.Matrix([drift, sp.diff(drift, radial)])
radial_parameter_jacobian = fold_equations.jacobian([eta_one, eta_two])
equal("retained radial fold parameter determinant",
      radial_parameter_jacobian.det(), 2*radial**3)
fold_parameters = sp.solve(list(fold_equations.subs(radial, rho)), [eta_one, eta_two])
equal("retained radial fold eta1", fold_parameters[eta_one], ell*rho**4)
equal("retained radial fold eta2", fold_parameters[eta_two], -2*ell*rho**2)
equal("retained radial fold curvature",
      sp.diff(drift, radial, 2).subs(fold_parameters).subs(radial, rho), 8*ell*rho**3)
weighted_z, s_one, s_two = sp.symbols("weighted_z s_one s_two")
weighted_drift = sp.expand(drift.subs({radial: rho*weighted_z,
                                      eta_one: rho**4*s_one,
                                      eta_two: rho**2*s_two})/rho**5)
weighted_equations = sp.Matrix([weighted_drift, sp.diff(weighted_drift, weighted_z)])
equal("scaled limiting fold determinant",
      weighted_equations.jacobian([s_one, s_two]).det().subs(weighted_z, 1), 2)

result = {
    "status": "PASS",
    "method": "Exact symbolic polynomial, rational, differential and chain-rule identities; no numerical sampling",
    "source_commit": SOURCE_COMMIT,
    "software": {"python": platform.python_version(), "sympy": sp.__version__},
    "source_anchors": {
        "source/sections/02-actual-returns.tex": [
            "eq:passage:local-geometry", "eq:passage:quadratic-jet",
            "eq:passage:strength-jet", "lem:passage:finite-inverse",
            "eq:passage:threshold-formula"],
        "source/sections/08-realization.tex": [
            "eq:specified-contact-controls", "eq:specified-four-jacobian",
            "eq:harvesting:model", "eq:harvesting:perturbed-chart",
            "eq:harvesting:clock", "eq:harvesting:strengths",
            "eq:harvesting:effort-relation", "eq:harvesting:rank-matrix",
            "eq:harvesting:lower-rank"],
        "source/sections/09-focus-appendix.tex": [
            "eq:harvesting:core-chart", "eq:harvesting:weighted-flow",
            "eq:harvesting:scaled-return", "eq:harvesting:fold-data"]},
    "passage": {
        "input": "dz/dt=(z-t^2)*beta(t,z)/(1+(z-t^2)*alpha(t,z))",
        "scaled_initial_value_problem": "dw/ds=r*(w-s^2)*beta(r*s,r^2*w)/(1+r^2*(w-s^2)*alpha(r*s,r^2*w)); w(-1,r)=1",
        "jet_convention": "beta=b+b_t*t+b_z*z+b_tt*t^2+terms of higher weighted degree; alpha=a+terms of positive weighted degree; weight(t)=1, weight(z)=2",
        "w_coefficients_r_1_to_3": [string(term) for term in terms],
        "exit_s_coefficients_r_1_to_3": [string(e) for e in exit_coefficients],
        "physical_Gamma_r_squared_through_r_5": string(gamma),
        "C1": string(C1),
        "C2": string(gamma.coeff(r, 4)),
        "aligned_composition_root_through_degree_4": string(composition),
        "inverse_root_through_degree_4": string(inverse),
        "activation_inverse_in_root_delta_through_degree_5": string(inverse_physical),
        "noninverse_example": {
            "first_field": "dt/dtime=1; dz/dtime=(z-t^2)*(b+kappa*t^2)",
            "second_field": "dt/dtime=1; dz/dtime=(z-t^2)*(-b)",
            "opposite_leading_strengths": ["4*b/3", "-4*b/3"],
            "composition_q_5_over_2_coefficient": string(generic_remainder)},
        "special_symmetry": "For alpha=0 and constant beta=b, f_{-b}(t,z)=-f_b(-t,z); time reversal pairs full local passages. This extra symmetry is not assumed by finite cancellation.",
        "activation_scope": "The second physical threshold is Gamma1^{-1}(delta), not generally delta. For delta<=0 it is delta; positive-delta inverse jet above uses t=sqrt(delta), A=C1/2."},
    "local_control_rank": {
        "control_order": ["eta0", "eta1", "kappa", "d"],
        "observation_order": ["a", "b", "m", "delta"],
        "controlled_contact_beta_at_touch": string(controlled_beta_at_contact),
        "jacobian": [[string(e) for e in row] for row in local_jacobian.tolist()],
        "determinant": string(local_jacobian.det())},
    "harvesting_rank_algebra": {
        "equilibrium_y": string(ye),
        "physical_chart_d": string(d_chart),
        "physical_chart_h": string(h_chart),
        "canonical_linear_one_turn_radial_rate_nu_derivative_at_hopf":
            string(linear_rate_nu_at_hopf),
        "linear_rate_scope": "A canonical linear angular normalization; nonzero for beta>1/80. It does not bound the actual return sensitivity at a positive fold radius.",
        "physical_chart_embedding_jacobian_rows_beta_d_h_u_v_L_U":
            [[string(e) for e in row] for row in chart_jacobian.tolist()],
        "embedding_minor_rows_d_L_v_U": string(embedding_minor.det()),
        "strengths": [string(strength1), string(strength2)],
        "canceling_effort_v": string(effort_ratio),
        "upper_contact_y": string(upper_y),
        "natural_control_order": ["nu", "L", "v", "U"],
        "natural_observation_jacobian": [[string(e) for e in row] for row in natural.tolist()],
        "structural_determinant": string(natural_determinant),
        "u_0_fold_determinant": string(fold_determinant),
        "required_analytic_nonzero_factors": ["a_nu", "c0", "s1=zeta_L", "s2=delta_U", "K2=4/(3*omega2)"],
        "actual_rank_status": "Conditional algebra verified. Actual-return sensitivities, existence, signs and persistence are analytic obligations, not values computed by this script."},
    "retained_radial_fold": {
        "drift": string(drift),
        "parameter_determinant": string(radial_parameter_jacobian.det()),
        "eta1": string(fold_parameters[eta_one]),
        "eta2": string(fold_parameters[eta_two]),
        "curvature": string(sp.diff(drift, radial, 2).subs(fold_parameters).subs(radial, rho)),
        "scaled_limiting_parameter_determinant": "2",
        "scope": "Finite polynomial algebra, not the actual-return remainder estimate or actual fold existence"},
    "checks_passed": len(checks),
    "checks": checks,
    "limitations": [
        "No existence or uniform bounds for the active ODE, moving-exit branch, signed companion or inverse are certified.",
        "No common first-return tube, first entry/exit, strict activity, no extra seam event or port admissibility is certified.",
        "Finite Taylor identities do not certify equality or inversion of complete generic passage germs.",
        "The existence and orientation of transported inactive first integrals and positive contact clocks remain analytic assumptions.",
        "No actual numerical fold or actual sensitivities a_nu, c0, I_x(p1), I_x(p2) are computed.",
        "Rank-four persistence from u=0 to fixed sufficiently small positive u uses smooth dependence and a nonzero-margin argument not certified here.",
        "Positive ecological parameters, threshold ordering, compact positive-quadrant containment and fitted threshold existence are not certified.",
        "Normal-form validity, derivative-controlled remainder estimates, parameter smoothness and the IFT step for an actual cyclic fold are not certified.",
        "No uniformity as rho or u tends to zero is claimed; the manuscript fixes a small positive rho, then a small positive u before the unfolding scale."]}

print(json.dumps(result, indent=2, sort_keys=True))
