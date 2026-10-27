# Author source

Run `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex`
from this directory. The output is build/main.pdf. The direct-sections/
directory includes all figure and retained generated-coordinate TeX inputs.
The bibliography is references.bib. Requires latexmk, pdfLaTeX and BibTeX.

Requires the LaTeX packages named in the supplied source. The frozen author
PDF is ../paper.pdf. Build output bytes can depend on the TeX environment.
