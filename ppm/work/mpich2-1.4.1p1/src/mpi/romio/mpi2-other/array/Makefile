CC          = icc
AR          = ar cr
LIBNAME     = /scratch/ptc/ppm/work/mpich2-1.4.1p1/lib/libmpich.a
RANLIB      = ranlib
srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpi/romio/mpi2-other/array
CC_SHL      = true
SHLIBNAME   = /scratch/ptc/ppm/work/mpich2-1.4.1p1/lib/libmpich

INCLUDE_DIR = -I. -I../../include -I${srcdir}/../../mpi-io -I${srcdir}/../../adio/include -I../../adio/include
CFLAGS      =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -O2  -DFORTRANUNDERSCORE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_ROMIOCONF_H -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 $(MPIOPROFILE) $(INCLUDE_DIR)

top_builddir  = /scratch/ptc/ppm/work/mpich2-1.4.1p1
LIBTOOL       = 
C_COMPILE_SHL = $(CC_SHL)

VPATH=.:${srcdir}

MPIO_OBJECTS = subarray.o darray.o

MPIO_TMP_POBJECTS = subarray.p darray.p

MPIO_REAL_POBJECTS = _subarray.o _darray.o

all: $(LIBNAME)
	@if [ "none" != "none" ] ; then \
	    $(MAKE) $(SHLIBNAME).la ;\
	fi

.SUFFIXES: $(SUFFIXES) .p .sp .lo

.c.o:
	$(CC) $(CFLAGS) -c $<
.c.lo:
	$(C_COMPILE_SHL) $(CFLAGS) -c $< -o _s$*.o
	@mv -f _s$*.o $*.lo

.c.p:
	@cp $(srcdir)/$*.c _$*.c
	$(CC) $(CFLAGS) -c _$*.c
	@rm -f _$*.c

.c.sp:
	$(C_COMPILE_SHL) $(CFLAGS) -c $< -o _$*.lo

profile: $(MPIO_REAL_POBJECTS)
	$(AR) $(LIBNAME) $(MPIO_REAL_POBJECTS)
	$(RANLIB) $(LIBNAME) 
	@if [ "none" != "none" ] ; then \
	    $(MAKE) P$(SHLIBNAME).la ;\
	fi
	@rm -f _*.o

$(LIBNAME): $(MPIO_OBJECTS)
	$(AR) $(LIBNAME) $(MPIO_OBJECTS)
	$(RANLIB) $(LIBNAME) 

MPIO_LOOBJECTS = $(MPIO_OBJECTS:.o=.lo)
$(SHLIBNAME).la: $(MPIO_LOOBJECTS) 
	$(AR) $(SHLIBNAME).la $(MPIO_LOOBJECTS)

# --------------------------------------------------------------------------
# We use P$(SHLBNAME) simply to distinguish the regular lib build from the
# profiling library build
MPIO_LOPOBJECTS = $(MPIO_REAL_POBJECTS:.o=.lo)
MPIO_TMP_LOPOBJECTS = ${MPIO_TMP_POBJECTS:.p=.sp}
P$(SHLIBNAME).la: $(MPIO_TMP_LOPOBJECTS)
	$(AR) $(SHLIBNAME).la $(MPIO_LOPOBJECTS)
# --------------------------------------------------------------------------

clean:
	@rm -f *.o *.lo *.gcno *.gcda *.bb *.bbg
	@rm -f ${srcdir}/*.gcno ${srcdir}/*.gcda 
	@rm -f ${srcdir}/*.bb ${srcdir}/*.bbg

# Rules for the profiling objects
_subarray.o: subarray.c
	$(CC) $(CFLAGS) -c $(srcdir)/subarray.c -o _subarray.o
_darray.o: darray.c
	$(CC) $(CFLAGS) -c $(srcdir)/darray.c -o _darray.o

