# This Makefile.in created by simplemake.  Do not edit

# f90/misc/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/f90/misc
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/f90/misc
LIBS            = -lpthread 
MAKE            = make
FC               = /scratch/ptc/prefix-O2/mpich2/1.4.1p1/bin/mpif90
FCFLAGS         =  
FC_COMPILE      = $(FC) $(FCFLAGS) $(FCINCLUDES)
FC_LINK        = $(FC)
FC_SHL          = @FC_SHL@
FC_COMPILE_SHL = $(FC_SHL) $(FCFLAGS) $(FCINCLUDES)
LDFLAGS     =    

VPATH = .

.SUFFIXES:
.SUFFIXES: .o .lo .f90 
.f90.o:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  FC              $<" ; \
	else \
	  echo $(FC_COMPILE) -c $< ; \
	fi
	@$(FC_COMPILE) -c $<
.f90.lo:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  FC              $<" ; \
	else \
	  echo $(FC_COMPILE_SHL) -c $< -o _s$*.o ; \
	echo mv -f _s$*.o $*.lo ; \
	fi
	@$(FC_COMPILE_SHL) -c $< -o _s$*.o
	@mv -f _s$*.o $*.lo
.f90:
	$(FC_COMPILE) -o $* $< $(LDFLAGS) $(LIBS)
apply: 
	$(ACTION) $(ACTION_INPUT)
clean: 
	-@rm -f *.o ${srcdir}/*.o sizeof2
	-@rm -f sizeof2.exe
	-@rm -f core core.[0-9]*
	-@rm -f *.lo ${srcdir}/*.lo
	-@rm -f *.irlog*
distclean:  clean distclean-xxx remove-makefile
distclean-xxx:
	-@rm -f Makefile
	-@rm -rf autom4te*.cache
	-@rm -rf $(DEPS_DIR)
	-@rm -f TAGS
remove-makefile:
	rm -f Makefile
remove-genmakefiles:
	rm -f  Makefile
all-redirect: sizeof2

sizeof2: sizeof2.o   
	$(FC_LINK) $(FCFLAGS) $(LDFLAGS) -o sizeof2 sizeof2.o  ${LIBS}
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
SOURCES         =  sizeof2.f90
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------

.POSIX:

FORCE_TARGET:

