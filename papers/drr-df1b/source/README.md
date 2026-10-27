# Author source

Run `python3 build.py` from this directory. It writes build intermediates
to build/ and the sanitized PDF to DF1b-exact2.pdf. Preserve build.py,
tools/sanitize-pdf-metadata.py and figures/df1b-graphic.pdf. The build requires
Python 3 with pypdf, plus latexmk, pdfLaTeX and BibTeX. No mathematical
computation is run by the build entry point.

Requires the LaTeX packages named in the supplied source. The frozen author
PDF is ../paper.pdf. Build output bytes can depend on the TeX environment.
