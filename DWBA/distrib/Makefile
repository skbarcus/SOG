# This is a make file for fitr.  The following commands are defined:
# make			Create program on source directory (don't install
#                       in destination)
# make all			"
# make program			"
# make clean		Remove object files from the source directory
# make install		Move the program from the source directory to the
#                       destination
# make update		Combines program creation with installation in
#                       destination directory

# Final result of this make file
PROGRAM=fitr

# Required object files
OBJS=	fitr.o \
	ravpro.o \
	chgden-fb.o \
	indeq2.o \
	deq2.o \
	invrt2.o \
	prompt.o \
	hellom.o \
	tod.o \
	nospac.o \
	nsp1.o \
	nsp2.o \
	ccout.o \
	cpubgn.o

# Source file names (syntax is the same as SRCS=$(patsubst %.o,%.F,$(OBJS)))
SRCS=$(OBJS:%.o=%.F)

# Place to store the final result (no last /)
DEST=$(HOME)/Tritium/Analysis/SOG/DWBA/distrib

# Name of this file
MAKEFILE=Makefile

# Library directories (include the -L) -- note that HP-UX uses a
# different syntax for library path names, but that is handled
# later in the HP specific code.
LIBDIR=-L$(HOME)/lib/$(OSNAME)

# Libraries to link to (other operating system specific libraries are defined
# by OSLIBS).
LIBS=

ifeq ($(OSNAME),HP-UX)
# Fortran command
FC=f77

# Fortran flags for HP-UX (not all are absolutely necessary)
# +E1	Permit calls to "industry standard" (though not Fortran standard)
#	subroutines and functions (DATE, TIME, RAN, ...)
# +E6	Permit variable formats
# -K	Static variables in subroutines (ie they keep their values accross
#	calls
# +T	Generate traceback information on an error
# -O	Optimization at level 2
# +ppu	Put a _ after all external references (necessary to link to CERNLIB
#	and therefore our standard).
# -V	To get a listing in .l file.
# +R	Include symbol table cross reference in the .l file
# -C    Enable range checking
FFLAGS=+E1 +E6 -K +T -O +ppu

# Loader options for HP-UX
# +FPD			Enable so-called "fast" mode of dealing with
#			underflows, which means set the result to 0.
# +U77			Permit linking to the Fortran run-time library
# -Wl,-L <directory>	Look in <directory> for libraries
LDFLAGS=+FPD +U77 $(LIBDIR:-L%=-Wl,-L%)

# Libraries required on this operating system
OSLIBS=
endif

ifeq ($(OSNAME),SunOS)
# Fortran command
FC=f77

# Fortran flags for SunOS
# -C	Enable range checking
# -O	Optimization at level 3
# -vax=misalign Selects all VMS Fortran extensions
# -Xlist Generate a listing in .lst file
FFLAGS=-O -vax=misalign

# Loader options for SunOS
# -L <directory>	Look in <directory> for libraries
LDFLAGS=$(LIBDIR)

# Libraries required on this operating system
OSLIBS=
endif

ifeq ($(OSNAME),OSF)
# Fortran command
FC=f77

# Fortran flags for OSF (eg axherx at hermes)
# -static Static variables in subroutines (ie they keep their values across
#         calls).
FFLAGS=-static

# Loader options for OSF
# -L <directory>        Look in <directory> for libraries
LDFLAGS=$(LIBDIR)

# Libraries required on this operating system
OSLIBS=
endif

ifeq ($(OSNAME),IRIX)
# Fortran command
FC=f77

# Fortran flags for SGI running IRIX
# -static Static variables in subroutines (ie they keep their values across
#         calls).
# -G0     No data is accessed from the global pointer
# -mips2  Generate code using the MIPS 2 instruction set
# -O      Invoke the global ucode optimizer
FFLAGS=-static -G0 -mips2 -O

# Loader options for SGI running IRIX
# -L <directory>        Look in <directory> for libraries
LDFLAGS=$(LIBDIR)

# Libraries required on this operating system
OSLIBS=
endif

# The following code creates a variable that is "yes" if the operating system
# is Linux or Cygwin.
ifeq ($(OSNAME),Linux)
LINCYG=yes
else
ifeq ($(OSNAME),Cygwin)
LINCYG=yes
else
LINCYG=no
endif
endif

ifeq ($(LINCYG),yes)
ifdef ABSOFT
# If the ABSOFT Fortran 77 compiler is available, use it.

# Fortran command
FC=f77

# Fortran flags for the Absoft compiler on Linux
# -g77			Fold all symbols to lower case, static storage of
#			local variables, character args passed as g77,
#			append underscore to subroutine names. (equivalent
#			to -f -s -N90 -B108 for the f77 compiler)
# -V			Accept VAX tab source code (tab to start) (f77 compiler)
# -N3			Include record length information for SEQUENTIAL
#			UNFORMATTED files
# -N51			Interpret the RECL specifier as the number of
#			32-bit words in a record
# -N34			Align common variables to avoid alignment errors.
# -g			Debugging enabled
# -O1			Optimize the code (could be O2 for even
#			beter optimization)
# -Dlinux		Define "linux" for the preprocessor
# -Dunix		Define "unix" for the preprocessor
# -Dabsoft		Define "absoft" for the preprocessor
FFLAGS=-Dlinux -Dunix -Dabsoft -g77 -N3 -N51 -N34 -V -O1

# Loader options for Linux with Absoft.  Note that the environment
# variable ABSOFT should point to the root directory of the ABSOFT
# compilers.  The environment variable ABSOFT_F2C should contain
# the -L and -l options required to load the appropriate f2c library.
# -L <directory>        Look in <directory> for libraries
# -X-static		Tell ld to use only static libraries (so that
#			executables can be copied to and run on other systems).
LDFLAGS=-X-static -L $(ABSOFT)/lib $(LIBDIR) 

# Libraries required on this operating system
OSLIBS=-lV77 -lU77 $(ABSOFT_F2C)

else
# There is no Absoft compiler.
ifeq ($(OSNAME),Linux)
# On systems where the "modern" gnu Fortran (gfortran) is available,
# we use that.  Note though that gfortran does not understand "tab"
# formatted (VAX Fortran).  So we use a script to run "notab" to
# convert what IS accepted.
FC=gfortran-run.bash

# Fortran (gfortran) flags for Linux
# -I.			Make sure the current directory is searched for
#			INCLUDE files.
# -fbounds-check	Enable generation of run-time checks
#                       for array subscripts and substring start
#                       and end points against the (locally) declared
#                       minimum and maximum values.
# -fno-backslash	Treat \ as a normal character
# -fno-automatic	Statically allocate variables in subroutines
# -ffloat-store 	Causes gfortran to try to remove the extra precision
#			produced in double precision calculations, when
#			the floating-point unit is set to work in IEEE 854
#			"extended precision" mode, by "spilling"
#			data from floating-point registers into memory.  This
#			typically involves a big performance hit, but can
#			correct the operation of programs that will NEVER
#			find two double precision numbers to be equal (say
#			when one is in a register and the other is in memory).
#			(gcc option)
# -ffixed-line-length-0 The entire code line is meaningful and continued
#			character constants never have implicit spaces
#			appendFixed form was traditionally used in older Fortran programsed to them to fill out the line
# -ffixed-form		Use traditional fixed format command lines (statement
#			number in col 1-5, continuation in 6, command in 7...)
# -g			Debugging enabled
# -O			Optimize the code (equivalent to -O1, this
#			could be -O2 or even -O3 for even beter
#			optimization)
# -Wall			Warn when variables are unused (except
#			for declaration) and when uninitialized
#			variables are used.
# -Dlinux		Define "linux" for the preprocessor
# -Dunix		Define "unix" for the preprocessor
# -Dgfortran		Define "gfortran" for the preprocessor
# -save-temps		Save temporary files (like .cpp)
FFLAGS=-I. -Wall -Dlinux -Dunix -Dgfortran -fno-backslash -ffixed-line-length-0 -fno-automatic -O

# Loader options for Linux
# -L <directory>        Look in <directory> for libraries
LDFLAGS=$(LIBDIR)

# Libraries required on this operating system
OSLIBS=
# Also include trig routines for degrees on Linux.  One of the OSLIBS+=
# lines below should be un-commented depending on whether single or
# double precision routines are needed.
#OSLIBS+= -ldtrigd
#OSLIBS+= -ltrigd

else
# On Cygwin compile with g77
# Fortran command
FC=g77

# Fortran (g77) flags for Linux
# -I.			Make sure the current directory is searched for
#			INCLUDE files.
# -fvxt			Enable VAX Fortran (rather than G77)
#			interpretation of various constructs.
# -ffortran-bounds-check Enable generation of run-time checks
#                       for array subscripts and substring start
#                       and end points against the (locally) declared
#                       minimum and maximum values.
# -fno-backslash	Treat \ as a normal character
# -fcase-lower		Map to lower case (case insensitive)
#			(note that most libraries are compiled
#			with lower case symbols).
# -fno-automatic	Statically allocate variables in
#			subroutines
# -finit-local-zero	Set un-initialized variables to 0
#			(though normally this should be done
#			explicitly in the code).
# -ffloat-store 	Causes g77 to try to remove the extra precision
#			produced in double precision calculations, when
#			the floating-point unit is set to work in IEEE 854
#			"extended precision" mode, by "spilling"
#			data from floating-point registers into memory.  This
#			typically involves a big performance hit, but can
#			correct the operation of programs that will NEVER
#			find two double precision numbers to be equal (say
#			when one is in a register and the other is in memory).
#			(gcc option)
# -g			Debugging enabled
# -O			Optimize the code (equivalent to -O1, this
#			could be -O2 or even -O3 for even beter
#			optimization)
# -Wall			Warn when variables are unused (except
#			for declaration) and when uninitialized
#			variables are used.
# -Dlinux		Define "linux" for the preprocessor
# -Dunix		Define "unix" for the preprocessor
# -save-temps		Save temporary files (like .cpp)
FFLAGS=-I. -Wall -Dlinux -Dunix -fvxt -fno-backslash -fcase-lower -fno-automatic -O

# Loader options for Linux
# -L <directory>        Look in <directory> for libraries
LDFLAGS=$(LIBDIR)

# Libraries required on this operating system
OSLIBS=
# Also include trig routines for degrees on Linux.  One of the OSLIBS+=
# lines below should be un-commented depending on whether single or
# double precision routines are needed.
#OSLIBS+= -ldtrigd
#OSLIBS+= -ltrigd
endif
endif
endif

#Attempt to ignore OSNAME stuff and just use gfortran (copied from above).
# On systems where the "modern" gnu Fortran (gfortran) is available,
# we use that.  Note though that gfortran does not understand "tab"
# formatted (VAX Fortran).  So we use a script to run "notab" to
# convert what IS accepted.
FC=gfortran

# Fortran (gfortran) flags for Linux
# -I.			Make sure the current directory is searched for
#			INCLUDE files.
# -fbounds-check	Enable generation of run-time checks
#                       for array subscripts and substring start
#                       and end points against the (locally) declared
#                       minimum and maximum values.
# -fno-backslash	Treat \ as a normal character
# -fno-automatic	Statically allocate variables in subroutines
# -ffloat-store 	Causes gfortran to try to remove the extra precision
#			produced in double precision calculations, when
#			the floating-point unit is set to work in IEEE 854
#			"extended precision" mode, by "spilling"
#			data from floating-point registers into memory.  This
#			typically involves a big performance hit, but can
#			correct the operation of programs that will NEVER
#			find two double precision numbers to be equal (say
#			when one is in a register and the other is in memory).
#			(gcc option)
# -ffixed-line-length-0 The entire code line is meaningful and continued
#			character constants never have implicit spaces
#			appendFixed form was traditionally used in older Fortran programsed to them to fill out the line
# -ffixed-form		Use traditional fixed format command lines (statement
#			number in col 1-5, continuation in 6, command in 7...)
# -g			Debugging enabled
# -O			Optimize the code (equivalent to -O1, this
#			could be -O2 or even -O3 for even beter
#			optimization)
# -Wall			Warn when variables are unused (except
#			for declaration) and when uninitialized
#			variables are used.
# -Dlinux		Define "linux" for the preprocessor
# -Dunix		Define "unix" for the preprocessor
# -Dgfortran		Define "gfortran" for the preprocessor
# -save-temps		Save temporary files (like .cpp)
FFLAGS=-I. -Wall -Dlinux -Dunix -Dgfortran -fno-backslash -ffixed-line-length-0 -fno-automatic -O

# Loader options for Linux
# -L <directory>        Look in <directory> for libraries
LDFLAGS=$(LIBDIR)

# Libraries required on this operating system
OSLIBS=
# Also include trig routines for degrees on Linux.  One of the OSLIBS+=
# lines below should be un-commented depending on whether single or
# double precision routines are needed.
#OSLIBS+= -ldtrigd
#OSLIBS+= -ltrigd
#End attempt to avoid OSNAME stuff.

# "make all" and "make program" are the same as "make"
all:		$(PROGRAM)
program:        $(PROGRAM)
$(PROGRAM):     $(OBJS)
		$(FC) $(LDFLAGS) $(OBJS) $(LIBS) $(OSLIBS) -o $(PROGRAM)

# "make install" transfers the program to the destination directory
install:	$(PROGRAM)
		mv $(PROGRAM) $(DEST)

# "make update" combines the creation of the program with its installation.
update:		$(DEST)/$(PROGRAM)
$(DEST)/$(PROGRAM): $(SRCS)
		@make -f $(MAKEFILE) DEST=$(DEST) install

# "make clean" removes all old object files.  This uses a "phony target"
# (ie target without dependencies)
.PHONY:		clean
clean:
		rm -f $(OBJS)
