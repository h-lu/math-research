# Author source

Build the main paper first:

```sh
cd manuscript
make check
cd ../supplements
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build ghost-paths.tex
```

The mathematical supplement imports the main build’s auxiliary labels.
`make check` also uses ripgrep (`rg`) for log checks. For a TeX-only build,
run the latexmk command in manuscript/ directly. Keep the manuscript/ and
supplements/ relative layout intact. Outputs are manuscript/build/main.pdf
and supplements/build/ghost-paths.pdf.

Requires the LaTeX packages named in the supplied source. The frozen author
PDF is ../paper.pdf. Build output bytes can depend on the TeX environment.
