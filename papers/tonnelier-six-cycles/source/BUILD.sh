#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
latexmk -g -pdf -interaction=nonstopmode -halt-on-error \
  -usepretex='\def\RewriteRelease{1}' -outdir=build paper/manuscript.tex
