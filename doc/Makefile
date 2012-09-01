# Copyright 2012 Wenjun Deng <wdeng@wdeng.info>
#
# This file is part of PIC1D-PETSc
#
# PIC1D-PETSc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIC1D-PETSc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PIC1D-PETSc.  If not, see <http://www.gnu.org/licenses/>.


# Makefile for compiling PIC1D-PETSc documentation

.PHONY : clean cleanall

formulation.pdf : formulation.tex
	latexmk -pdf formulation

clean :
	rm -f formulation.aux formulation.log formulation.out formulation.fdb_latexmk

cleanall : clean
	rm -f formulation.pdf

