# ============================================================================
# Name        : BIZEP - Plastic BEM for COPPE class BEM-II
# Author      : J. C. DE FARIA TELLES 
# Version     :
# Copyright   : J. C. DE FARIA TELLES
# Description : Makefile for debug Fortran source
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran

all: src/BIZEP.F
	$(FORTRAN_COMPILER) -Og -g -Wall -Wextra -fcheck=all -std=legacy -static --fixed-form \
		-o bin/BIZEP.exe \
		src/BIZEP.F

clean:
	rm -f bin/BIZEP *.mod
