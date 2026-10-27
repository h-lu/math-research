#!/bin/sh
set -eu
export SOURCE_DATE_EPOCH=1792454400
export FORCE_SOURCE_DATE=1
export TZ=UTC
exec latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=build main.tex
