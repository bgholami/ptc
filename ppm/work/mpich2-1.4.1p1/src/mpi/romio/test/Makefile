CC = /scratch/ptc/ppm/work/mpich2-1.4.1p1/bin/mpicc
F77 = /scratch/ptc/ppm/work/mpich2-1.4.1p1/bin/mpif77
INCLUDE_DIR = 
# because := is not universally avalible, we have to play games to use the
# user-specified LDFLAGS and OUR_LIBS env. variables (if set)
OUR_LIBS =   ${LDFLAGS} ${LIBS}
USER_CFLAGS =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers  $(INCLUDE_DIR)
USER_FFLAGS =  -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpl/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/openpa/src -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/datatype -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/common/locks -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/include -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpid/ch3/channels/nemesis/nemesis/utils/monitor -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers -I/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/util/wrappers  $(INCLUDE_DIR)
CTESTS = simple perf async coll_test coll_perf misc file_info excl large_array \
     atomicity noncontig i_noncontig noncontig_coll split_coll shared_fp \
     large_file psimple error status noncontig_coll2 aggregation1 aggregation2 \
     async-multiple ordered_fp hindexed
FTESTS = fcoll_test fperf fmisc pfcoll_test 
srcdir=/scratch/ptc/ppm/work/mpich2-1.4.1p1/src/mpi/romio/test

VPATH=.:${srcdir}


all: c-test fperf fcoll_test fmisc pfcoll_test

c-test: $(CTESTS)

default: all

.c.o:
	$(CC) $(USER_CFLAGS) -c $<

.f.o:
	$(F77) $(USER_FFLAGS) -c $<

.o:
	$(CC) $(USER_CFLAGS) -o  $@ $< $(OUR_LIBS)

clean: 
	rm -f .P* PI* *.o 
	rm -f $(CTESTS) $(FTESTS)
	@-rm -f work.pc work.pcl

#
# The Fortran tests must be built with the Fortran linker
fperf: fperf.o
	$(F77) $(USER_FFLAGS) -o  fperf fperf.o $(OUR_LIBS)

fcoll_test: fcoll_test.o
	$(F77) $(USER_FFLAGS) -o  fcoll_test fcoll_test.o $(OUR_LIBS)

fmisc: fmisc.o
	$(F77) $(USER_FFLAGS) -o  fmisc fmisc.o $(OUR_LIBS)

pfcoll_test: pfcoll_test.o	
	$(F77) $(USER_FFLAGS) -o  pfcoll_test pfcoll_test.o $(OUR_LIBS)

#testing: 
#	-./runtests $(TESTARGS)

cleanall:
	-@rm -f .P* PI* misc.c large_file.c fcoll_test.f fmisc.f fperf.f \
	mpif.h $(CTESTS) $(FTESTS) *.out runtests pfcoll_test.f
