# This Makefile.in created by simplemake.  Do not edit

# src/mpid/common/datatype/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
AR              = ar
AR_FLAGS        = cr
RANLIB          = ranlib
top_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype
master_top_srcdir  = /scratch/ptc/ppm/work/mpich2-1.4.1p1
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype
DEFS            = -DHAVE_CONFIG_H -I. -I${srcdir}
CPPFLAGS        =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers
MAKE            = make
MPILIBNAME      = mpich
CC              = icc
CFLAGS          = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2   -O2 $(MPICH2_MAKE_CFLAGS)
C_COMPILE       = $(CC) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)
LIBTOOL         = 
CC_SHL          = true
C_COMPILE_SHL   = $(CC_SHL) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)

VPATH = .
HEADERS = dataloop/dataloop_parts.h dataloop/dataloop_create.h mpid_dataloop.h mpid_datatype.h \
	mpid_ext32_segment.h segment_states.h
INCLUDES = -I${srcdir} -I${top_builddir}/src/include -I${master_top_srcdir}/src/include

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
	(cd dataloop && ${MAKE} apply )
clean: 
	-@rm -f *.o ${srcdir}/*.o 
	-@rm -f .libstamp0 .libstamp1 
	-@rm -f *.lo ${srcdir}/*.lo
	(cd dataloop && ${MAKE} clean )
distclean:  clean distclean-xxx remove-makefile
distclean-xxx:
	if [ -s dataloop/Makefile ] ; then (cd dataloop && ${MAKE} distclean ) ; fi
	-@rm -f Makefile
	-@rm -rf autom4te*.cache
	-@rm -rf $(DEPS_DIR)
	-@rm -f TAGS
remove-makefile:
	rm -f Makefile
remove-genmakefiles:
	rm -f  Makefile
all-redirect:
	@(cd dataloop && ${MAKE} all )
	${MAKE} all-local
all-local: .libstamp0 ../../../../lib/lib${MPILIBNAME}.a
	@if [ "none" != "none" ] ; then \
	    ${MAKE} .libstamp1 || exit 1; \
	fi

.libstamp0 ../../../../lib/lib${MPILIBNAME}.a: mpid_type_debug.o mpid_segment.o \
	mpid_ext32_segment.o mpid_ext32_datatype.o mpid_type_dup.o mpir_type_get_contig_blocks.o \
	mpir_type_flatten.o mpid_type_commit.o mpid_type_zerolen.o mpid_type_vector.o \
	mpid_type_contiguous.o mpid_type_create_resized.o mpid_type_blockindexed.o \
	mpid_type_indexed.o mpid_type_struct.o mpid_type_get_envelope.o mpid_datatype_free.o \
	mpid_datatype_contents.o mpid_type_get_contents.o mpid_type_create_pairtype.o \
	mpid_contents_support.o 
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  AR ${AR_FLAGS} ../../../../lib/lib${MPILIBNAME}.a $?" ; \
	else \
	  echo ${AR} ${AR_FLAGS} ../../../../lib/lib${MPILIBNAME}.a $? ; \
	fi
	@${AR} ${AR_FLAGS} ../../../../lib/lib${MPILIBNAME}.a $?
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  RANLIB          lib${MPILIBNAME}.a" ; \
	else \
	  echo ${RANLIB} ../../../../lib/lib${MPILIBNAME}.a ; \
	fi
	@${RANLIB} ../../../../lib/lib${MPILIBNAME}.a
	date > .libstamp0
.libstamp1 ../../../../lib/lib${MPILIBNAME}.la: mpid_type_debug.lo mpid_segment.lo \
	mpid_ext32_segment.lo mpid_ext32_datatype.lo mpid_type_dup.lo \
	mpir_type_get_contig_blocks.lo mpir_type_flatten.lo mpid_type_commit.lo \
	mpid_type_zerolen.lo mpid_type_vector.lo mpid_type_contiguous.lo \
	mpid_type_create_resized.lo mpid_type_blockindexed.lo mpid_type_indexed.lo \
	mpid_type_struct.lo mpid_type_get_envelope.lo mpid_datatype_free.lo \
	mpid_datatype_contents.lo mpid_type_get_contents.lo mpid_type_create_pairtype.lo \
	mpid_contents_support.lo 
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  AR ${AR_FLAGS} ../../../../lib/lib${MPILIBNAME}.la $?" ; \
	else \
	  echo ${AR} ${AR_FLAGS} ../../../../lib/lib${MPILIBNAME}.la $? ; \
	fi
	@${AR} ${AR_FLAGS} ../../../../lib/lib${MPILIBNAME}.la $?
	date > .libstamp1
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
SOURCES         =  mpid_type_debug.c mpid_segment.c mpid_ext32_segment.c mpid_ext32_datatype.c \
	mpid_type_dup.c mpir_type_get_contig_blocks.c mpir_type_flatten.c mpid_type_commit.c \
	mpid_type_zerolen.c mpid_type_vector.c mpid_type_contiguous.c mpid_type_create_resized.c \
	mpid_type_blockindexed.c mpid_type_indexed.c mpid_type_struct.c mpid_type_get_envelope.c \
	mpid_datatype_free.c mpid_datatype_contents.c mpid_type_get_contents.c \
	mpid_type_create_pairtype.c mpid_contents_support.c 
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------

.POSIX:

FORCE_TARGET:

