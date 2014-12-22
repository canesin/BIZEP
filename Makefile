# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran

all: src/BIZEP.F
	$(FORTRAN_COMPILER) -O2 -std=legacy --fixed-form -g \
		-o bin/BIZEP \
		src/BIZEP.F

clean:
	rm -f bin/BIZEP *.mod
