# This Makefile.in created by simplemake.  Do not edit

# errors/attr/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/errors/attr
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/errors/attr
DEFS            = -DHAVE_CONFIG_H -I. -I${srcdir}
CPPFLAGS        = 
LIBS            = -lpthread 
MAKE            = make
CC              = /scratch/ptc/prefix-O2/mpich2/1.4.1p1/bin/mpicc
CFLAGS          =   $(MPICH2_MAKE_CFLAGS)
C_COMPILE       = $(CC) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)
C_LINK         = $(CC)
LDFLAGS     =    ../../util/mtest.o

VPATH = .
INCLUDES = -I../../include -I${srcdir}/../../include

.SUFFIXES:
.SUFFIXES: .o .c 
.c.o:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CC              $<" ; \
	else \
	  echo $(C_COMPILE) -c $< ; \
	fi
	@$(C_COMPILE) -c $<
.c:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CC              $<" ; \
	else \
	  echo $(C_COMPILE) -o $* $< $(LDFLAGS) $(LIBS) ; \
	fi
	@$(C_COMPILE) -o $* $< $(LDFLAGS) $(LIBS)
apply: 
	$(ACTION) $(ACTION_INPUT)
clean: 
	-@rm -f *.o ${srcdir}/*.o keyvalmis
	-@rm -f keyvalmis.exe
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
all-redirect: keyvalmis

keyvalmis: keyvalmis.o    ../../util/mtest.o
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o keyvalmis keyvalmis.o  ${LIBS}
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
SOURCES         =  keyvalmis.c
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------

# Use extra programs to list programs that should not be built by default.


../../util/mtest.o:
	(cd ../../util && $(MAKE) mtest.o)

.POSIX:

FORCE_TARGET:

