# Fixed-parameter interval witnesses

The checker programs, exact parameter inputs and retained reports are
unchanged. Each of two fixed harvesting-model tuples has a retained
directed-interval certificate for at least six distinct simple periodic orbits.
They do not certify an exact global count or proximity to a neutral reference.

These entry points require Python 3 and python-flint (reference version 0.9.0,
native FLINT 3.6.0). From this directory run:

```sh
python3 certify_candidate.py --output main-rerun.json
python3 certify_secondary.py --output secondary-rerun.json
```

The main exact tuple is embedded in certify_candidate.py; the secondary
tuple is read from secondary-candidate.json. Retained expected reports are
interval-certificate.json and secondary-interval-certificate.json. The broader
historical environment is in requirements.txt; these two certificate entry
points only need the standard library and python-flint. Elapsed time and
runtime-version metadata can vary.
No interval certificate or high-precision integration was rerun in this
packaging pass. The analytic theorems are proved in the papers.
