CC          = icc
AR          = ar cr
RANLIB      = ranlib
LIBNAME     = /scratch/ptc/ppm/work/mpich2-1.4.1p1/lib/libmpich.a
srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpi/romio/adio/ad_zoidfs
CC_SHL      = true
SHLIBNAME   = /scratch/ptc/ppm/work/mpich2-1.4.1p1/lib/libmpich

INCLUDES    = -I. -I${srcdir}/../include -I../include -I../../include  -I${srcdir}/../../../../include -I../../../../include
CFLAGS      =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -O2  -DFORTRANUNDERSCORE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_ROMIOCONF_H -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 $(INCLUDES)

top_builddir  = /scratch/ptc/ppm/work/mpich2-1.4.1p1
LIBTOOL       = 
C_COMPILE_SHL = $(CC_SHL)

VPATH=.:${srcdir}

AD_ZOIDFS_OBJECTS = ad_zoidfs.o ad_zoidfs_close.o ad_zoidfs_common.o \
	ad_zoidfs_delete.o ad_zoidfs_fcntl.o ad_zoidfs_flush.o \
	ad_zoidfs_io.o ad_zoidfs_open.o ad_zoidfs_resize.o ad_zoidfs_features.o\
	ad_zoidfs_read_list.o ad_zoidfs_write_list.o

default: $(LIBNAME)
	@if [ "none" != "none" ] ; then \
	    $(MAKE) $(SHLIBNAME).la ;\
	fi

.SUFFIXES: $(SUFFIXES) .p .lo

.c.o:
	$(CC) $(CFLAGS) -c $<
.c.lo:
	$(C_COMPILE_SHL) $(CFLAGS) -c $< -o _s$*.o
	@mv -f _s$*.o $*.lo


$(LIBNAME): $(AD_ZOIDFS_OBJECTS)
	$(AR) $(LIBNAME) $(AD_ZOIDFS_OBJECTS)
	$(RANLIB) $(LIBNAME) 

AD_ZOIDFS_LOOBJECTS=$(AD_ZOIDFS_OBJECTS:.o=.lo)
$(SHLIBNAME).la: $(AD_ZOIDFS_LOOBJECTS)
	$(AR) $(SHLIBNAME).la $(AD_ZOIDFS_LOOBJECTS)

coverage:
	-@for file in  ${AD_ZOIDFS_OBJECTS:.o=.c} ; do \
		gcov -b -f $$file ; done

clean:
	@rm -f *.o *.lo *.gcno *.gcda *.bb *.bbg
	@rm -f ${srcdir}/*.gcno ${srcdir}/*.gcda 
	@rm -f ${srcdir}/*.bb ${srcdir}/*.bbg
