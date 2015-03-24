# This Makefile.in created by simplemake.  Do not edit

# f77/comm/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/f77/comm
bindir          = ${exec_prefix}/bin
exec_prefix     = ${prefix}
prefix          = /scratch/ptc/prefix-O2/mpich2/1.4.1p1
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/f77/comm
LIBS            = -lpthread 
MAKE            = make
F77             = /scratch/ptc/prefix-O2/mpich2/1.4.1p1/bin/mpif77
FFLAGS          =  
F77_COMPILE     = $(F77) $(FFLAGS) $(F77INCLUDES)
F77_LINK       = $(F77)
LDFLAGS     =    ../util/mtestf.o 

VPATH = .

.SUFFIXES:
.SUFFIXES: .o .f 
.f.o:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  F77             $<" ; \
	else \
	  echo $(F77_COMPILE) -c $< ; \
	fi
	@$(F77_COMPILE) -c $<
.f:
	$(F77_COMPILE) -o $* $< $(LDFLAGS) $(LIBS)
apply: 
	$(ACTION) $(ACTION_INPUT)
clean: 
	-@rm -f *.o ${srcdir}/*.o commerrf commnamef
	-@rm -f commerrf.exe commnamef.exe
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
all-redirect: commnamef commerrf

commerrf: commerrf.o    ../util/mtestf.o
	$(F77_LINK) $(FFLAGS) $(LDFLAGS) -o commerrf commerrf.o  ${LIBS}
commnamef: commnamef.o    ../util/mtestf.o
	$(F77_LINK) $(FFLAGS) $(LDFLAGS) -o commnamef commnamef.o  ${LIBS}
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
SOURCES         =  commnamef.f commerrf.f
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------


testing:
	../../runtests -srcdir=$(srcdir) -tests=testlist \
			-mpiexec=$(bindir)/mpiexec \
		   	-xmlfile=summary.xml

../util/mtestf.o: 
	(cd ../util && $(MAKE) mtestf.o)

.POSIX:

FORCE_TARGET:

