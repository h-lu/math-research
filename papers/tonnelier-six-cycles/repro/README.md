# Directed-interval six-cycle verification

These files derive from Supplement S1 v1.0.0. The candidate retains all
mathematical inputs; one machine-name provenance field was removed for
portable distribution. The checker sources and reference output are unchanged. A fresh output will
record a different `candidate_sha256` because that field hashes the entire
JSON file, including provenance; this packaging change leaves every
mathematical input unchanged.
The original frozen attachment remains available at https://doi.org/10.5281/zenodo.23113416.

Reference environment: Python 3.14.4, python-flint 0.9.0, FLINT 3.6.0,
mpmath 1.3.0. Dependencies are listed in `requirements.txt`.

From this directory, run:

```sh
python3 certify_six_cycle_arb.py six-cycle-001.json --precision 480 --radius 1e-25 --slices 512 --output certificate.json
```

This certifies the fixed six-cycle witness with directed Arb arithmetic.
It does not formally verify the general analytic finiteness or lower-bound proofs.
A fresh certificate replay was not run during this collection preparation.
