# This Makefile.in created by simplemake.  Do not edit

# cxx/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/cxx
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/cxx
MAKE            = make

VPATH = .

.SUFFIXES:
# Some make programs complain if no suffixes are set
.SUFFIXES: .c .o
apply: 
	$(ACTION) $(ACTION_INPUT)
	for dir in util attr pt2pt init info comm coll datatype errhan io rma spawn ; do \
		(cd $$dir && ${MAKE} apply ) ; done
clean: 
	-@rm -f *.o ${srcdir}/*.o 
	-@rm -f *.lo ${srcdir}/*.lo
	for dir in util attr pt2pt init info comm coll datatype errhan io rma spawn ; do \
		(cd $$dir && ${MAKE} clean ) ; done
distclean:  clean distclean-xxx remove-makefile
distclean-xxx:
	for dir in util attr pt2pt init info comm coll datatype errhan io rma spawn ; do \
		if [ ! -s $$dir/Makefile ] ; then continue ; fi ;\
		(cd $$dir && ${MAKE} distclean ) ; done
	-@rm -f Makefile testlist
	-@rm -rf autom4te*.cache
	-@rm -rf $(DEPS_DIR)
	-@rm -f TAGS
remove-makefile:
	rm -f Makefile
remove-genmakefiles:
	rm -f  Makefile
all-redirect:
	@(cd util && ${MAKE} all )
	@(cd attr && ${MAKE} all )
	@(cd pt2pt && ${MAKE} all )
	@(cd init && ${MAKE} all )
	@(cd info && ${MAKE} all )
	@(cd comm && ${MAKE} all )
	@(cd coll && ${MAKE} all )
	@(cd datatype && ${MAKE} all )
	@(cd errhan && ${MAKE} all )
	@(cd io && ${MAKE} all )
	@(cd rma && ${MAKE} all )
	@(cd spawn && ${MAKE} all )

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------

# A special testing target to simplify testing just the C++ code
testing:
	../runtests -srcdir=$(srcdir) -tests=testlist \
		   -mpiexec=/scratch/ptc/prefix-O2/mpich2/1.4.1p1/bin/mpiexec \
		   -xmlfile=summary.xml


.POSIX:

FORCE_TARGET:

