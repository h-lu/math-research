.PHONY: all check clean

all: check

check:
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
	! rg 'LaTeX Warning: (Reference|Citation).*undefined|There were undefined references|multiply defined|Overfull' build/main.log

clean:
	latexmk -c -outdir=build main.tex
