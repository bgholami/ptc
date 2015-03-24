CC          = icc
AR          = ar cr
RANLIB      = ranlib
LIBNAME     = /scratch/ptc/ppm/work/mpich2-1.4.1p1/lib/libmpich.a
srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpi/romio/adio/ad_testfs
CC_SHL      = true
SHLIBNAME   = /scratch/ptc/ppm/work/mpich2-1.4.1p1/lib/libmpich

#DEFS            = -DHAVE_CONFIG_H -I. -I${srcdir}
CPPFLAGS        =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers
INCLUDES = -I. -I${srcdir}/../include -I../include -I../../include  -I${srcdir}/../../../../include -I../../../../include
#CFLAGS      =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -O2  -DFORTRANUNDERSCORE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_ROMIOCONF_H -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 $(INCLUDES) 
CFLAGS          = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -O2  -DFORTRANUNDERSCORE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_ROMIOCONF_H -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

top_builddir  = /scratch/ptc/ppm/work/mpich2-1.4.1p1
LIBTOOL       = 
C_COMPILE       = $(CC) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)
C_COMPILE_SHL   = $(CC_SHL) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)
#C_COMPILE_SHL = $(CC_SHL) 

VPATH=.:${srcdir}

AD_TESTFS_OBJECTS = ad_testfs_close.o ad_testfs_read.o ad_testfs_rdcoll.o \
	ad_testfs_wrcoll.o ad_testfs_open.o ad_testfs_write.o \
	ad_testfs_done.o ad_testfs_fcntl.o ad_testfs_iread.o \
	ad_testfs_iwrite.o ad_testfs_wait.o ad_testfs_flush.o \
	ad_testfs_seek.o ad_testfs_resize.o ad_testfs_hints.o \
	ad_testfs_delete.o ad_testfs.o

default: $(LIBNAME)
	@if [ "none" != "none" ] ; then \
	    $(MAKE) $(SHLIBNAME).la ;\
	fi

.SUFFIXES: $(SUFFIXES) .p .lo

#	$(CC) $(CFLAGS) -c $<
.c.o:
	@if [ "x$(VERBOSE)" != "x1" ] ; then \
	  echo "  CC              $<" ; \
	else \
	  echo $(C_COMPILE) -c $< ; \
	fi
	@$(C_COMPILE) -c $<
#	$(C_COMPILE_SHL) $(CFLAGS) -c $< -o _s$*.o
#	@mv -f _s$*.o $*.lo
.c.lo:
	@if [ "x$(VERBOSE)" != "x1" ] ; then \
	  echo "  CC              $<" ; \
	else \
	  echo $(C_COMPILE_SHL) -c $< -o _s$*.o ; \
	  echo mv -f _s$*.o $*.lo ; \
	fi
	@$(C_COMPILE_SHL) -c $< -o _s$*.o
	@mv -f _s$*.o $*.lo

$(LIBNAME): $(AD_TESTFS_OBJECTS)
	$(AR) $(LIBNAME) $(AD_TESTFS_OBJECTS)
	$(RANLIB) $(LIBNAME) 

AD_TESTFS_LOOBJECTS = $(AD_TESTFS_OBJECTS:.o=.lo)

$(SHLIBNAME).la: $(AD_TESTFS_LOOBJECTS)
	$(AR) $(SHLIBNAME).la $(AD_TESTFS_LOOBJECTS)

clean:
	@rm -f *.o *.lo *.gcno *.gcda *.bb *.bbg
	@rm -f ${srcdir}/*.gcno ${srcdir}/*.gcda 
	@rm -f ${srcdir}/*.bb ${srcdir}/*.bbg

tags: TAGS
SOURCES = ${AD_TESTFS_OBJECTS:.o=.c}
HEADERS = 
TAGS:${HEADERS} ${SOURCES}
	here=`cd ../../../../.. && pwd` ; cd ${srcdir} && etags -o $$here/TAGS --append ${HEADERS} ${SOURCES}
