#
# Automagically generated by Approximatrix Simply Fortran 2.15
#
FC="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\gfortran.exe"
CC="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\gcc.exe"
AR="C:\Program Files (x86)\Simply Fortran 2\mingw-w64\bin\ar.exe"
WRC="C:\PROGRA~2\SIMPLY~1\MINGW-~1\bin\windres.exe"
RM=rm -f

IDIR=
# -I error: Directory C:\Users\fcanesin\AppData\Local\\sfpm\64\include does not exist

LDIR=
# -L error: Directory C:\Users\fcanesin\AppData\Local\\sfpm\64\lib does not exist


OPTFLAGS= -g

SPECIALFLAGS=$(IDIR)

RCFLAGS=-O coff

PRJ_FFLAGS=

PRJ_CFLAGS=

PRJ_LFLAGS=

FFLAGS=$(SPECIALFLAGS) $(OPTFLAGS) $(PRJ_FFLAGS) -Jmodules 

CFLAGS=$(SPECIALFLAGS) $(OPTFLAGS) $(PRJ_CFLAGS)

"build\BIZEPB.o": ".\BIZEPB.FOR"
	@echo Compiling .\BIZEPB.FOR
	@$(FC) -c -o "build\BIZEPB.o" $(FFLAGS) ".\BIZEPB.FOR"

clean: .SYMBOLIC
	@echo Deleting build\BIZEPB.o and related files
	@$(RM) "build\BIZEPB.o"
	@echo Deleting target.exe
	@$(RM) "target.exe"

"target.exe":  "build\BIZEPB.o"
	@echo Generating target.exe
	@$(FC) -o "target.exe" -static "build\BIZEPB.o" $(LDIR) $(PRJ_LFLAGS)

all: "target.exe" .SYMBOLIC
