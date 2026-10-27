# Maintaining this collection

The single source for the paper index is `papers/<slug>/paper.json`.
Edit that metadata and run `python3 scripts/catalog.py`; do not hand-edit
the generated root README or catalog.json. Python 3 with its standard
library is sufficient. No service account, hosting framework, or external
publishing credential is required.

## Update a paper

1. Keep its slug stable. Edit author source in its `source/` directory.
2. Build using that paper's documented requirements. Inspect the resulting
   PDF, citations, figures, and author AI-use disclosures.
3. Replace `paper.pdf` with the intended author manuscript. Update the
   version and artifact SHA-256 values in paper.json; record the revision
   in that paper's CHANGELOG.md.
4. Record exactly what any calculation or formal verification checks. Keep
   metadata limited to mathematical identity, version, artifacts and scope.
5. Run `python3 scripts/catalog.py`, then `python3 scripts/catalog.py --check`.
   Review the exact changed files before committing.

Use Git history for earlier revisions. Important frozen versions may receive
descriptive tags; do not copy multiple full version directories or move an
existing frozen tag. Update the DOI only when its corresponding deposit exists.

## Add a paper

Copy `templates/paper/` to `papers/<stable-lowercase-slug>/`. Replace every
placeholder, add `paper.pdf`, author sources and only necessary reproducibility
materials. List every public artifact with its role and SHA-256 in paper.json.
Keep supporting artifacts self-contained; preserve relative build paths.
Do not add internal reviews, correspondence, incomplete research notes,
credentials, or machine-specific paths. Verify rights for each paper and
attachment; an absent license is not permission to invent one.

Record the original authors, exact title, manuscript version, citation, license evidence if available, and validation scope.
The template is deliberately incomplete and will not pass validation until
real artifacts and their checksums are supplied.
