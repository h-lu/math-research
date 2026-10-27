# Author source

Run `python3 build.py` from this directory. It invokes pdfLaTeX with shell
escape disabled for three passes. Requires pdfLaTeX, Poppler pdfinfo and the
LaTeX packages named in paper/preamble.tex. The output is paper/main.pdf.
The frozen author PDF is ../paper.pdf. The bibliography and four TikZ figures
are embedded in the supplied author sources.
