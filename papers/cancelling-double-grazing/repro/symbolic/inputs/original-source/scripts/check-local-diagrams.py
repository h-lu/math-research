#!/usr/bin/env python3
"""Reproduce the finite algebra checks for the quartic local diagram.

Run from the project root: python3 scripts/check-local-diagrams.py
Only JSON is printed; no source or output files are changed.
Scope: exact quartic discriminant identities and sample region labels for
quartic-local-diagram.tex.
Geometry is produced directly by the retained TikZ source using an exact
polynomial discriminant, not numerical continuation or simulated ecological
trajectories.  Exact rational real-root isolation checks the five printed
region labels; floating-point midpoints below are only readable summaries.
These are finite production checks, not existence proofs for actual cycles.
"""
import json
from pathlib import Path
import re
import sympy as sp

y, t, u, v, w = sp.symbols('y t u v w', real=True)
quartic = y**4 + u*y**2 + v*y + w
quartic_fold = {u: -1, v: 2*t-4*t**3, w: 3*t**4-t**2}

# The retained PGF formulas need parentheses around a possibly negative
# parameter before exponentiation.  Otherwise PGF interprets "-.5^2"
# as a negative number, distorting the negative-t half of the diagram.
figure_dir = Path(__file__).resolve().parent.parent/'figures'
source = (figure_dir/'quartic-local-diagram.tex').read_text()
assert r'\t^' not in source

def zero(expr):
    return sp.simplify(expr) == 0

assert zero(quartic.subs(quartic_fold).subs(y, t))
assert zero(sp.diff(quartic, y).subs(quartic_fold).subs(y, t))
assert zero(sp.discriminant(quartic, y).subs(quartic_fold))
for root in [-sp.sqrt(sp.Rational(1, 6)), sp.sqrt(sp.Rational(1, 6))]:
    vals = {key: sp.simplify(sp.sympify(value).subs(t, root))
            for key, value in quartic_fold.items()}
    assert zero(sp.diff(quartic, y, 2).subs(vals).subs(y, root))
    assert not zero(sp.diff(quartic, y, 3).subs(vals).subs(y, root))
assert sp.factor(quartic.subs({u: -1, v: 0, w: sp.Rational(1, 4)})) == (2*y**2-1)**2/4

def sample(coefficients, companions):
    poly = sp.Poly(sum(sp.Rational(str(a))*y**(len(coefficients)-i-1)
                       for i, a in enumerate(coefficients)), y)
    exact_intervals = sp.polys.polytools.intervals(poly, eps=sp.Rational(1, 10**10))
    # intervals returns only real root intervals, with exact multiplicity.
    roots = [float(sum(interval)/2) for interval, _ in exact_intervals]
    multiplicities = [multiplicity for _, multiplicity in exact_intervals]
    assert all(m == 1 for m in multiplicities)
    return {'different_real_roots': len(roots),
            'total_simple_cycles': len(roots)+companions,
            'isolated_roots': roots, 'multiplicities': multiplicities}

quartic_samples = {}
sample_labels = [('0', '0.46', 2), ('0', '0.12', 6),
                 ('-0.76', '0.13', 4), ('0.76', '0.13', 4), ('0', '-0.105', 4)]
printed_labels = re.findall(
    r'\\node\[region[^\]]*\] at \(([^,]+),([^\)]+)\) \{\$(\d+)\$\};', source)
assert [(sp.Rational(V), sp.Rational(W), int(count))
        for V, W, count in printed_labels] == [
    (sp.Rational(V), sp.Rational(W), count) for V, W, count in sample_labels]
for V, W, expected in sample_labels:
    result = sample([1, 0, -1, V, W], 2)
    assert result['total_simple_cycles'] == expected
    quartic_samples[f'V={V},W={W}'] = result

print(json.dumps({
    'status': 'all retained quartic exact identities and root-count labels passed',
    'scope': 'Exact quartic discriminant identities and sample region labels for quartic-local-diagram.tex',
    'geometry_status': 'EXACT quartic discriminant in actual prepared controls',
    'quartic_slice': {'u': '-epsilon^2', 'similarity': 'y=epsilon Y, V=v/epsilon^3, W=w/epsilon^4',
                      'fold': ['2t-4t^3', '3t^4-t^2'],
                      'triple_points': ['(-4/(3sqrt(6)),-1/12)', '(4/(3sqrt(6)),-1/12)'],
                      'two_double_point': '(0,1/4)',
                      'printed_region_labels_checked': len(printed_labels),
                      'samples': quartic_samples},
    'root_isolation': 'SymPy exact rational polynomial real intervals, epsilon=1e-10',
    'sympy_version': sp.__version__,
}, indent=2))
