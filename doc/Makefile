.PHONY : clean cleanall

formulation.pdf : formulation.tex
	latexmk -pdf formulation

clean :
	rm -f formulation.aux formulation.log formulation.out formulation.fdb_latexmk

cleanall : clean
	rm -f formulation.pdf

