# This Makefile.in created by simplemake.  Do not edit

# examples/Makefile.  Generated from Makefile.in by configure.

# --------------------------------------------------------------------------
all: all-redirect

SHELL           = /bin/sh
abs_builddir    = /scratch/ptc/ppm/work/mpich2-1.4.1p1/examples
srcdir          = .
abs_srcdir      = /scratch/ptc/ppm/work/mpich2-1.4.1p1/examples
top_srcdir      = ..
MAKE            = make
MPILIBNAME      = mpich
PMPILIBNAME     = pmpich

C_COMPILE       = $(CC) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)
C_LINK         = $(CC)
CC_SHL          = true
C_COMPILE_SHL   = $(CC_SHL) $(DEFS) $(INCLUDES) $(CFLAGS) $(CPPFLAGS)

VPATH = .
MPICHDIR=..
CC       = $(MPICHDIR)/bin/mpicc
CFLAGS   =
CXXFLAGS =
FFLAGS   =
FCFLAGS  =
LDFLAGS  =
CPPFLAGS = 
DEFS     = 
INCLUDES = -I../src/include -I${top_srcdir}/src/include
LIBS   =

EXTRA_PROGRAMS = pmandel pmandel_spawn pmandel_service pmandel_spaserv 	         pmandel_fence hellow icpi parent child srtest 		 spawn_merge_parent spawn_merge_child1 spawn_merge_child2

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
	-@rm -f *.o ${srcdir}/*.o cpi pmandel_service pmandel icpi srtest hellow child pmandel_fence pmandel_spaserv spawn_merge_child2 spawn_merge_child1 parent spawn_merge_parent pmandel_spawn
	-@rm -f cpi.exe pmandel_service.exe pmandel.exe icpi.exe srtest.exe hellow.exe child.exe pmandel_fence.exe pmandel_spaserv.exe spawn_merge_child2.exe spawn_merge_child1.exe parent.exe spawn_merge_parent.exe pmandel_spawn.exe
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
all-redirect: cpi

cpi: cpi.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o cpi cpi.o  -lm ${LIBS}
pmandel_service: pmandel_service.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o pmandel_service pmandel_service.o  -lm ${LIBS}
pmandel: pmandel.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o pmandel pmandel.o  -lm ${LIBS}
icpi: icpi.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o icpi icpi.o  -lm ${LIBS}
srtest: srtest.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o srtest srtest.o  ${LIBS}
hellow: hellow.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o hellow hellow.o  ${LIBS}
child: child.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o child child.o  ${LIBS}
pmandel_fence: pmandel_fence.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o pmandel_fence pmandel_fence.o  -lm ${LIBS}
pmandel_spaserv: pmandel_spaserv.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o pmandel_spaserv pmandel_spaserv.o  -lm ${LIBS}
spawn_merge_child2: spawn_merge_child2.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o spawn_merge_child2 spawn_merge_child2.o  ${LIBS}
spawn_merge_child1: spawn_merge_child1.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o spawn_merge_child1 spawn_merge_child1.o  ${LIBS}
parent: parent.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o parent parent.o  ${LIBS}
spawn_merge_parent: spawn_merge_parent.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o spawn_merge_parent spawn_merge_parent.o  ${LIBS}
pmandel_spawn: pmandel_spawn.o   ../lib/lib${MPILIBNAME}.a
	$(C_LINK) $(CFLAGS) $(LDFLAGS) -o pmandel_spawn pmandel_spawn.o  -lm ${LIBS}
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Dummy target
dependencies:

# --------------------------------------------------------------------------
SOURCES         =  hellow.c srtest.c child.c parent.c cpi.c icpi.c pmandel.c pmandel_spawn.c \
	pmandel_service.c pmandel_spaserv.c pmandel_fence.c spawn_merge_parent.c \
	spawn_merge_child1.c spawn_merge_child2.c 
HEADERFILES     = $(HEADERS)
SOURCEFILES     = $(SOURCES)
# --------------------------------------------------------------------------
#OTHER_DIRS = f77 f90 cxx
# We want to replace the compilation and linking rule used by simplement


# We previously put "../lib/libopa.a" here, but that only covers one of the
# various OPA build/install cases.  Omit it for now.

# pmandel requires a separate set of socket calls (its a long story)
# and may not build on most platforms

# LIBS includes -lmpich and other libraries (e.g., -lpmpich if
# weak symbols are not supported)

# Use make all-programs to test building all example programs
all-programs: cpi $(EXTRA_PROGRAMS) 

.POSIX:

FORCE_TARGET:

