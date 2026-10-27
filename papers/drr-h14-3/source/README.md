# Author source

From this directory run pdfLaTeX three times:

```sh
mkdir -p build
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
```

The output is build/main.pdf. All figures and the bibliography are TeX inputs.

Requires the LaTeX packages named in the supplied source. The frozen author
PDF is ../paper.pdf. Build output bytes can depend on the TeX environment.
