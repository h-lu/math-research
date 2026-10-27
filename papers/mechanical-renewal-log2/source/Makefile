SHELL := /bin/bash
export SOURCE_DATE_EPOCH := 1792713600
export FORCE_SOURCE_DATE := 1
TEX := main.tex prelim.tex mechanical-proof.tex counting.tex replication.tex renewal-conjecture.tex

.PHONY: all check clean
all: main.pdf

main.pdf: $(TEX) references.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
	cp build/main.pdf main.pdf

check:
	python3 repro/verify.py

clean:
	latexmk -C -outdir=build main.tex
