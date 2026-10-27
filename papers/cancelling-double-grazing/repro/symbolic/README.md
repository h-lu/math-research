# Symbolic calculations and diagnostic

The supplied programs, frozen mathematical inputs and expected outputs
are unchanged. Python dependencies are recorded in requirements.txt.
From this directory the following entry points reproduce exact checks:

```sh
python3 scripts/check_original.py focus
python3 scripts/check_original.py diagrams
python3 scripts/check_passage_rank.py
python3 scripts/check_quartic_curvature.py
python3 scripts/focus_complex_normal_form.py
```

The stored outputs include exact focus, passage-rank and quartic calculations.
The separate scripts/diagnose_harvesting.py and its recorded diagnostic/trajectory
outputs use floating arithmetic and are not interval certificates.
These checks supplement the analytic proofs; they do not formally verify them.
The additional independent calculations and floating diagnostic were not rerun
during this public-edition packaging. Main and supplement author checkers may
be run separately using their source guides.
