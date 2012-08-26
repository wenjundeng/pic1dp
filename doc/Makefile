.PHONY : clean cleanall

formulation.pdf : formulation.tex
	pdflatex formulation
	pdflatex formulation

clean :
	rm -f formulation.aux formulation.log formulation.out

cleanall : clean
	rm -f formulation.pdf

