# This Makefile.in created by simplemake.  Do not edit

# threads/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/threads
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/threads
MAKE            = make

VPATH = .

.SUFFIXES:
# Some make programs complain if no suffixes are set
.SUFFIXES: .c .o
apply: 
	$(ACTION) $(ACTION_INPUT)
	for dir in util pt2pt comm init spawn ; do \
		(cd $$dir && ${MAKE} apply ) ; done
clean:  clean-local
	-@rm -f *.o ${srcdir}/*.o 
	-@rm -f *.lo ${srcdir}/*.lo
	for dir in util pt2pt comm init spawn ; do \
		(cd $$dir && ${MAKE} clean ) ; done
distclean:  clean distclean-xxx remove-makefile
distclean-xxx:
	for dir in util pt2pt comm init spawn ; do \
		if [ ! -s $$dir/Makefile ] ; then continue ; fi ;\
		(cd $$dir && ${MAKE} distclean ) ; done
	-@rm -f Makefile testlist
	-@rm -rf autom4te*.cache
	if [ -s f77/Makefile ] ; then (cd f77 && $(MAKE) distclean ) ; fi
	if [ -s cxx/Makefile ] ; then (cd cxx && $(MAKE) distclean ) ; fi
	-@rm -rf $(DEPS_DIR)
	-@rm -f TAGS
remove-makefile:
	rm -f Makefile
remove-genmakefiles:
	rm -f  Makefile
all-redirect:
	@(cd util && ${MAKE} all )
	@(cd pt2pt && ${MAKE} all )
	@(cd comm && ${MAKE} all )
	@(cd init && ${MAKE} all )
	@(cd spawn && ${MAKE} all )

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------
#SUBDIRS = util attr basic datatype coll comm errhan group info init pt2pt rma spawn topo io  f77 f90 cxx .

# A special testing target to simplify running just the thread tests
testing:
	../runtests -srcdir=$(srcdir) -tests=testlist \
		   -mpiexec=/scratch/ptc/prefix-O2/mpich2/1.4.1p1/bin/mpiexec \
		   -xmlfile=summary.xml
clean-local:
	-rm -f summary.xml

.POSIX:

FORCE_TARGET:

