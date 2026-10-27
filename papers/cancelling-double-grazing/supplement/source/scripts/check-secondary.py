#!/usr/bin/env python3
"""Finite exact checks of the extracted cubic geometry and figure labels.

This checks polynomial identities, not the analytic actual-return estimates
or existence of physical cycles. Run from source/: python3 scripts/check-secondary.py.
"""
import json
from pathlib import Path
import sympy as s
z, t, alpha, d, B, u, v = s.symbols('z t alpha d B u v', real=True)
Bpos = s.symbols('B', positive=True)
Lam=s.symbols('Lambda', real=True)
F=z**3-z**2+Lam*z+alpha
assert s.expand(s.discriminant(F,z)) == Lam**2-4*Lam**3+4*alpha-18*Lam*alpha-27*alpha**2
boundary={alpha:-t**2*(1-2*t),Lam:2*t-3*t**2}
assert s.expand(F.subs(boundary)) == s.expand((z-t)**2*(z-(1-2*t)))
assert s.expand(F.subs({alpha:-s.Rational(1,27),Lam:s.Rational(1,3)})) == s.expand((z-s.Rational(1,3))**3)
# Finite first-only minimum and simultaneous endpoint constants.
eta=s.symbols('eta',real=True)
fine=eta-z**2+Bpos*z**3
assert s.simplify(fine.subs(z,2/(3*Bpos))) == eta-4/(27*Bpos**2)
assert s.factor(F.subs(alpha,0)) == z*(Lam+z**2-z)
assert s.discriminant(z**2-z+Lam,z) == 1-4*Lam
# Ordinary triple-root continuation and its two-control unfolding.
G=z**3-z**2+s.Rational(3,2)*B*d*z+alpha
J=s.Matrix([G,s.diff(G,z),s.diff(G,z,2)]).jacobian([alpha,d,z])
assert s.simplify(J.det().subs({z:s.Rational(1,3)})) == 9*B
curve=s.Matrix([-t**2+2*t**3,2*t*(2-3*t)/(3*B)])
assert s.simplify(s.det(s.Matrix.hstack(curve.diff(t,2),curve.diff(t,3))).subs(t,s.Rational(1,3))) == 48/B
# Prepared cubic, with exactly three simple companions outside its chart.
P=z**3+u*z+v
fold={u:-3*t**2,v:2*t**3}
assert s.expand(P.subs(fold).subs(z,t)) == 0
assert s.expand(s.diff(P,z).subs(fold).subs(z,t)) == 0
assert s.expand((4*u**3+27*v**2).subs(fold)) == 0
figure=(Path(__file__).resolve().parent.parent/'figures/cusp-local-diagram.tex').read_text()
assert r'\t^' not in figure
samples=[]
for U,V,total in [('-0.88','0.015',6),('0.26','0.25',4),('-0.72','0.44',4)]:
    poly=s.Poly(P.subs({u:s.Rational(U),v:s.Rational(V)}),z)
    intervals=s.polys.polytools.intervals(poly,eps=s.Rational(1,10**10))
    assert all(m==1 for _,m in intervals)
    assert len(intervals)+3==total
    samples.append({'u':U,'v':V,'real_roots':len(intervals),'total_simple_cycles':total})
# Check multiplicities on a nonzero fold and at the cusp, over exact rationals.
for controls,expected in [({u:-3,v:2},[1,2]),({u:0,v:0},[3])]:
    mult=sorted(m for _,m in s.polys.polytools.intervals(s.Poly(P.subs(controls),z)))
    assert mult==expected
print(json.dumps({
    'status':'all finite exact checks passed',
    'scope':'secondary cubic discriminant and boundary; positive-displacement constants; cusp Jacobians and prepared figure counts',
    'secondary_domain_discriminant':str(s.discriminant(F,z)),
    'triple_root_control':['-1/27','2/(9B)','1/3'],
    'positive_displacement_endpoints':['4/(27B^2)','1/(6B)'],
    'cusp_branch_determinant':'9B',
    'fold_curve_determinant':'48/B',
    'prepared_cubic_samples':samples,
    'fold_total_distinct_cycles':5,
    'cusp_total_distinct_cycles':4,
    'sympy_version':s.__version__,
    'limitations':'does not verify analytic estimates, root localization, actual periodic orbits, or theorem correctness'
},indent=2))
