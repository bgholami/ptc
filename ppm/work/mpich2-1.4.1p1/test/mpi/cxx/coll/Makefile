# This Makefile.in created by simplemake.  Do not edit

# cxx/coll/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/cxx/coll
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/test/mpi/cxx/coll
DEFS            = -DHAVE_CONFIG_H -I. -I${srcdir}
CPPFLAGS        = 
LIBS            = -lpthread 
MAKE            = make
CXX             = /scratch/ptc/prefix-O2/mpich2/1.4.1p1/bin/mpicxx
CXXFLAGS        = 
CXX_COMPILE     = $(CXX) $(DEFS) $(INCLUDES) $(CXXFLAGS) $(CPPFLAGS)
CXX_LINK       = $(CXX)
LDFLAGS     =    ../util/mtest.o 

VPATH = .
INCLUDES = -I../../include -I${srcdir}/../../include

.SUFFIXES:
.SUFFIXES: .o .cxx 
.cxx.o:
	@if [ "x$(VERBOSE)" != "x1" -a "x$(V)" != "x1" ] ; then \
	  echo "  CXX             $<" ; \
	else \
	  echo $(CXX_COMPILE) -c $< ; \
	fi
	@$(CXX_COMPILE) -c $<
.cxx:
	$(CXX_COMPILE) -o $* $< $(LDFLAGS) $(LIBS)
apply: 
	$(ACTION) $(ACTION_INPUT)
clean: 
	-@rm -f *.o ${srcdir}/*.o uallredx icallreducex reduceboolx alltoallw2x icscattervx arcomplex icbcastx uallreduce uscan icgathervx icgatherx uexscan icscatterx ureduce icbarrierx icalltoallx icallgatherx icreducex icallgathervx
	-@rm -f uallredx.exe icallreducex.exe reduceboolx.exe alltoallw2x.exe icscattervx.exe arcomplex.exe icbcastx.exe uallreduce.exe uscan.exe icgathervx.exe icgatherx.exe uexscan.exe icscatterx.exe ureduce.exe icbarrierx.exe icalltoallx.exe icallgatherx.exe icreducex.exe icallgathervx.exe
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
all-redirect: arcomplex uallredx uallreduce ureduce uscan uexscan icbcastx icbarrierx icgatherx\
	 icgathervx icscatterx icscattervx icallgatherx icallgathervx icallreducex icreducex\
	 icalltoallx alltoallw2x reduceboolx

uallredx: uallredx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o uallredx uallredx.o  ${LIBS}
icallreducex: icallreducex.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icallreducex icallreducex.o  ${LIBS}
reduceboolx: reduceboolx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o reduceboolx reduceboolx.o  ${LIBS}
alltoallw2x: alltoallw2x.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o alltoallw2x alltoallw2x.o  ${LIBS}
icscattervx: icscattervx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icscattervx icscattervx.o  ${LIBS}
arcomplex: arcomplex.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o arcomplex arcomplex.o  ${LIBS}
icbcastx: icbcastx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icbcastx icbcastx.o  ${LIBS}
uallreduce: uallreduce.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o uallreduce uallreduce.o  ${LIBS}
uscan: uscan.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o uscan uscan.o  ${LIBS}
icgathervx: icgathervx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icgathervx icgathervx.o  ${LIBS}
icgatherx: icgatherx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icgatherx icgatherx.o  ${LIBS}
uexscan: uexscan.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o uexscan uexscan.o  ${LIBS}
icscatterx: icscatterx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icscatterx icscatterx.o  ${LIBS}
ureduce: ureduce.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o ureduce ureduce.o  ${LIBS}
icbarrierx: icbarrierx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icbarrierx icbarrierx.o  ${LIBS}
icalltoallx: icalltoallx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icalltoallx icalltoallx.o  ${LIBS}
icallgatherx: icallgatherx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icallgatherx icallgatherx.o  ${LIBS}
icreducex: icreducex.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icreducex icreducex.o  ${LIBS}
icallgathervx: icallgathervx.o    ../util/mtest.o
	$(CXX_LINK) $(CXXFLAGS) $(LDFLAGS) -o icallgathervx icallgathervx.o  ${LIBS}
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
SOURCES         =  arcomplex.cxx uallredx.cxx uallreduce.cxx ureduce.cxx uscan.cxx uexscan.cxx \
	icbcastx.cxx icbarrierx.cxx icgatherx.cxx icgathervx.cxx icscatterx.cxx icscattervx.cxx \
	icallgatherx.cxx icallgathervx.cxx icallreducex.cxx icreducex.cxx icalltoallx.cxx \
	alltoallw2x.cxx reduceboolx.cxx 
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------


../util/mtest.o:
	(cd ../util && $(MAKE) mtest.o)

.POSIX:

FORCE_TARGET:

