# Exact finite focus calculation

This standalone script checks the finite algebra used for the smooth harvesting core. It derives the nonlinear field from the original model with `u = nu = 0`, rather than loading stored focus coefficients. It uses Python 3 and SymPy; it has been run with Python 3.14.4 and SymPy 1.14.0.

From the source root:

```sh
python3 scripts/check-focus.py
python3 scripts/check-focus.py --details
```

The script writes to standard output only and exits with a nonzero status if an assertion fails. The optional second command also prints the coefficients of `V_3,...,V_6` in the order `X^(n-j) Z^j`, `j = 0,...,n`. No data file, network access, ODE solver, parameter scan, interval arithmetic, or numerical approximation is required.

## Relation to the manuscript

1. The smooth core is equation `eq:harvesting:model` with `u = 0` and chart `eq:harvesting:core-chart` with `nu = 0`. The script substitutes `x = X + 1/5`, `y = 5*w*Z - X + 4/5`, and `beta = (25*w^2 + 1)/80`. Both coordinate identities in `eq:harvesting:rotation` are checked as exact rational identities.
2. It expands `(1 + 5*X)^(-1)` by a finite geometric sum, retaining only nonlinear degrees 2 through 5. These are precisely the field degrees that enter the derivative of `V_2+...+V_6` through degree 6.
3. For `n = 3,4,5,6`, it constructs `H_n` and solves `eq:harvesting:homological` in the monomial basis. It computes the radial obstruction using the exact circle moments
   
   `(2a)!(2b)! / [4^(a+b) a! b! (a+b)!]`.
   
   At even degrees it fixes the radial kernel by setting the `X^n` coefficient of `V_n` to zero. Every coefficient residual of the full homological equation is checked as identically zero.
4. The resulting `k_4(beta)`, `k_4(beta_0)`, `k_4'(beta_0)`, and `k_6(beta_0)` are compared exactly with `eq:harvesting:focus-data` and the immediately following derivative formula. The assertions also check the exact signs `omega_0^2 > 0`, `k_4'(beta_0) > 0`, and `k_6(beta_0) < 0`.

Here `k_4` is the energy focus coefficient denoted `ell_1(beta)` in the manuscript. The script does not identify every convention for a first or second Lyapunov coefficient. At `k_4 = 0`, the identity `V = R^2/2 + O(R^3)` and the sign of `k_6` give the nonzero negative quintic radial drift by the analytic argument in the manuscript.

The script checks only finite rational/symbolic algebra. It does not replace the analytic arguments for a smooth parameter-dependent finite normal form, its remainders, the derivative-controlled actual return, the cyclic fold, or the four-parameter realization.

## Provenance

This portable version was extracted from the existing exact homogeneous-coefficient calculation and checked against the current manuscript's core and appendix equations. It omits prior review metadata and filesystem assumptions, adds assertions for every stated identity, and makes no numerical claim.
