# This Makefile.in created by simplemake.  Do not edit

# src/mpi/timer/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
AR              = ar
AR_FLAGS        = cr
RANLIB          = ranlib
top_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpi/timer
master_top_srcdir  = /scratch/ptc/ppm/work/mpich2-1.4.1p1
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpi/timer
DEFS            = -DHAVE_CONFIG_H -I. -I${srcdir}
CPPFLAGS        =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers
MAKE            = make
MPILIBNAME      = mpich
PMPILIBNAME     = pmpich

CC              = icc
CFLAGS          = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2   -O2 $(MPICH2_MAKE_CFLAGS)
C_COMPILE       = $(CC) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)
LIBTOOL         = 
CC_SHL          = true
C_COMPILE_SHL   = $(CC_SHL) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)

VPATH = .
mpi_sources = wtime.c wtick.c
INCLUDES = -I../../include -I${master_top_srcdir}/src/include

.SUFFIXES:
.SUFFIXES: .o .lo .c 
.c.o:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CC              $<" ; \
	else \
	  echo $(C_COMPILE) -c $< ; \
	fi
	@$(C_COMPILE) -c $<
.c.lo:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CC              $<" ; \
	else \
	  echo $(C_COMPILE_SHL) -c $< -o _s$*.o ; \
	echo mv -f _s$*.o $*.lo ; \
	fi
	@$(C_COMPILE_SHL) -c $< -o _s$*.o
	@mv -f _s$*.o $*.lo
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
	-@rm -f *.o ${srcdir}/*.o 
	-@rm -f .libstamp0 .libstamp1 .libstamp2 .libstamp3 
	-@rm -f *.lo ${srcdir}/*.lo
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
all-redirect: .libstamp0 ../../../lib/lib${MPILIBNAME}.a
	@if [ -n "" ] ; then \
	    ${MAKE} .libstamp2 || exit 1; \
	fi
	@if [ "none" != "none" ] ; then \
	    ${MAKE} .libstamp1 || exit 1; \
	    if [ -n "" ] ; then \
	        ${MAKE} .libstamp3 || exit 1; \
	    fi ; \
	fi

.libstamp0 ../../../lib/lib${MPILIBNAME}.a: wtime.o wtick.o mpidtime.o 
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  AR ${AR_FLAGS} ../../../lib/lib${MPILIBNAME}.a $?" ; \
	else \
	  echo ${AR} ${AR_FLAGS} ../../../lib/lib${MPILIBNAME}.a $? ; \
	fi
	@${AR} ${AR_FLAGS} ../../../lib/lib${MPILIBNAME}.a $?
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  RANLIB          lib${MPILIBNAME}.a" ; \
	else \
	  echo ${RANLIB} ../../../lib/lib${MPILIBNAME}.a ; \
	fi
	@${RANLIB} ../../../lib/lib${MPILIBNAME}.a
	date > .libstamp0
.libstamp1 ../../../lib/lib${MPILIBNAME}.la: wtime.lo wtick.lo mpidtime.lo 
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  AR ${AR_FLAGS} ../../../lib/lib${MPILIBNAME}.la $?" ; \
	else \
	  echo ${AR} ${AR_FLAGS} ../../../lib/lib${MPILIBNAME}.la $? ; \
	fi
	@${AR} ${AR_FLAGS} ../../../lib/lib${MPILIBNAME}.la $?
	date > .libstamp1

# --------------------------------------------------------------------------
.libstamp2: _wtime.o _wtick.o 
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  AR ${AR_FLAGS} ../../../lib/lib${PMPILIBNAME}.a $?" ; \
	else \
	  echo ${AR} ${AR_FLAGS} ../../../lib/lib${PMPILIBNAME}.a $? ; \
	fi
	@${AR} ${AR_FLAGS} ../../../lib/lib${PMPILIBNAME}.a $?
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  RANLIB          lib${PMPILIBNAME}.a" ; \
	else \
	  echo ${RANLIB} ../../../lib/lib${PMPILIBNAME}.a ; \
	fi
	@${RANLIB} ../../../lib/lib${PMPILIBNAME}.a
	date > .libstamp2
.libstamp3: _wtime.lo _wtick.lo 
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  AR ${AR_FLAGS} ../../../lib/lib${PMPILIBNAME}.la $?" ; \
	else \
	  echo ${AR} ${AR_FLAGS} ../../../lib/lib${PMPILIBNAME}.la $? ; \
	fi
	@${AR} ${AR_FLAGS} ../../../lib/lib${PMPILIBNAME}.la $?
	date > .libstamp3
_wtime.o: wtime.c
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CC             -D<profiling> ${srcdir}/wtime.c" ; \
	else \
	  echo $(C_COMPILE) -c  $(srcdir)/wtime.c -o _wtime.o ; \
	fi
	@$(C_COMPILE) -c  $(srcdir)/wtime.c -o _wtime.o
_wtime.lo: wtime.c
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo " CC           -D<profiling> $(srcdir)/wtime.c" ; \
	else \
	  echo $(C_COMPILE_SHL) -c  $(srcdir)/wtime.c -o _wtime.lo ; \
	fi
	@$(C_COMPILE_SHL) -c  $(srcdir)/wtime.c -o _wtime.lo
_wtick.o: wtick.c
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CC             -D<profiling> ${srcdir}/wtick.c" ; \
	else \
	  echo $(C_COMPILE) -c  $(srcdir)/wtick.c -o _wtick.o ; \
	fi
	@$(C_COMPILE) -c  $(srcdir)/wtick.c -o _wtick.o
_wtick.lo: wtick.c
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo " CC           -D<profiling> $(srcdir)/wtick.c" ; \
	else \
	  echo $(C_COMPILE_SHL) -c  $(srcdir)/wtick.c -o _wtick.lo ; \
	fi
	@$(C_COMPILE_SHL) -c  $(srcdir)/wtick.c -o _wtick.lo
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
SOURCES         =  wtime.c wtick.c mpidtime.c
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------

#

.POSIX:

FORCE_TARGET:

