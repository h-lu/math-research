# Maintaining this collection

The single source for the paper index is `papers/<slug>/paper.json`.
Edit that metadata and run `python3 scripts/catalog.py`; do not hand-edit
the generated root README or catalog.json. Python 3 with its standard
library is sufficient. No service account, hosting framework, or external
publishing credential is required.

## Update a paper

1. Edit author source, keep the slug stable, and update the version and CHANGELOG.
2. Rebuild using the paper's source guide, copy the intended PDF to `paper.pdf`,
   and inspect text, figures, citations and the author's AI-use disclosure.
3. Run the paper's documented checks. Remove unlisted generated PDF copies or
   explicitly list any new intended artifact in paper.json before refreshing.
4. Refresh only the chosen paper, regenerate the index, and check the result:

```sh
python3 scripts/catalog.py --refresh-hashes your-paper-slug
python3 scripts/catalog.py
python3 scripts/catalog.py --check
python3 -m unittest discover -s tests
```

`--refresh-hashes` only updates checksums for artifacts already explicitly
listed in that paper.json. It refuses missing, unsafe, duplicate, or unlisted
files and validates the other papers before writing. It does not add artifacts,
regenerate indexes, commit, push, or publish. Review the diff before any Git
commit or push. Normal generation and `--check` always enforce the saved hashes.
Build directories, Python caches, and the LaTeX/Python intermediate suffixes
listed in scripts/catalog.py are excluded from inventory discovery; they are
never added to the public artifact list by refresh. An explicitly listed
supporting `.aux` file remains a checked source dependency.

Use Git history for earlier revisions. Important frozen versions may receive
descriptive tags; do not copy multiple full version directories or move an
existing frozen tag. Update the DOI only when its corresponding deposit exists.

## Add a paper

Copy `templates/paper/` to `papers/<stable-lowercase-slug>/`. Replace every
placeholder, add `paper.pdf`, author sources and only necessary reproducibility
materials. List every public artifact with its role and SHA-256 in paper.json.
For new artifacts, use 64 zeroes as the initial hash and then run the explicit
refresh command after supplying the files. The required paper guide is README.md;
list paper.pdf, source/README.md and the real author TeX source as artifacts.
Keep supporting artifacts self-contained; preserve relative build paths.
Do not add internal reviews, correspondence, incomplete research notes,
credentials, or machine-specific paths. Verify rights for each paper and
attachment; an absent license is not permission to invent one.

Record the original authors, exact title, manuscript version, citation, license evidence if available, and validation scope.
The template is deliberately incomplete and will not pass validation until
real artifacts and their checksums are supplied.
