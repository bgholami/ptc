########## MAKEFILE FOR PPM##########
#####################################
# Compiler
CC			=	icc
CFLAGS			=	-O2

CXX			=	icpc
CXXFLAGS		=	-O2

FC			=	ifort
FCFLAGS			=	-O2  

#####################################
######################################
# FOLDERS
DISTDIR		=	$(PWD)/distfiles
WRKDIR		=	$(PWD)/work


#########################################
##           MACHINE
#########################################
KNOWN_SYSTEM=no

# SUPERMUC #############################
ifeq ($(MYMACHINE), SUPERMUC)
#ZLIB_VERSION       = 1.2.5
#HDF5_VERSION       = 1.8.10-snap0
H5PART_VERSION     = 1.6.6
METIS_VERSION      = 4.0
MAKEDEPF90_VERSION = 2.8.8
PPM_VERSION        = 1.2.2

PREFIX          = /home/hpc/pr32ma/lu96mah2/prefix-supermuc

#ZLIB_BASE       = ${PREFIX}/zlib/${ZLIB_VERSION}
#HDF5_BASE       = ${PREFIX}/hdf5/${HDF5_VERSION}
H5PART_BASE     = ${PREFIX}/h5part/${H5PART_VERSION}
METIS_BASE      = ${PREFIX}/metis/${METIS_VERSION}
PPM_BASE        = ${PREFIX}/ppm/${PPM_VERSION}

MPIFC  = mpif90 
MPICC  = mpicc 
MPICXX = mpiCC

BUILD_MAKE_FLAGS     =    -j 8
KNOWN_SYSTEM=yes
endif

# MIYUKO ##############################
ifeq ($(MYMACHINE), miyuko)
ZLIB_VERSION       = 1.2.5
HDF5_VERSION       = 1.8.10-snap0
H5PART_VERSION     = 1.6.6
METIS_VERSION      = 4.0
MAKEDEPF90_VERSION = 2.8.8
PPM_VERSION        = 1.2.2

PREFIX=/home/adami/prefix-O2

ZLIB_BASE       = ${PREFIX}/zlib/${ZLIB_VERSION}
HDF5_BASE       = ${PREFIX}/hdf5/${HDF5_VERSION}
H5PART_BASE     = ${PREFIX}/h5part/${H5PART_VERSION}
METIS_BASE      = ${PREFIX}/metis/${METIS_VERSION}
PPM_BASE        = ${PREFIX}/ppm/${PPM_VERSION}

MPIFC  = mpif90 
MPICC  = mpicc 
MPICXX = mpicxx

BUILD_MAKE_FLAGS     =    -j4
KNOWN_SYSTEM=yes
endif

# MITSUKO ###########################################
ifeq ($(MYMACHINE), mitsuko)
#ZLIB_VERSION       = 1.2.5
#HDF5_VERSION       = 1.8.10-snap0
#H5PART_VERSION     = 1.6.6
METIS_VERSION      = 4.0
MAKEDEPF90_VERSION = 2.8.8
PPM_VERSION        = 1.2.2
MPICH2_VERSION     = 1.4.1p1

PREFIX=/scratch/ptc/prefix-O2

#ZLIB_BASE       = ${PREFIX}/zlib/${ZLIB_VERSION}
#HDF5_BASE       = ${PREFIX}/hdf5/${HDF5_VERSION}
#H5PART_BASE     = ${PREFIX}/h5part/${H5PART_VERSION}
METIS_BASE      = ${PREFIX}/metis/${METIS_VERSION}
PPM_BASE        = ${PREFIX}/ppm/${PPM_VERSION}
MPI_BASE        = ${PREFIX}/mpich2/${MPICH2_VERSION}

MPIFC  = mpif90
MPICC  = mpicc
MPICXX = mpicxx

BUILD_MAKE_FLAGS     =   
KNOWN_SYSTEM=yes
endif

##########################################
# ABORT
ifeq ($(KNOWN_SYSTEM), no)
abort:
	@echo MYMACHINE not defined or unknown!
endif


# LINUX CLUSTER
#-include Makefile.linuxcluster
#-include Makefile.miyuko-dguerster
#-include Makefile.miyuko
#-include Makefile.laptop

# SUPER MUC
#-include Makefile.supermuc
#-include Makefile.mpp
#-include Makefile.ice

#####################################
#####################################
# EXTRACT COMAND
EXTRACT_COMAND	=	scripts/extract


#######################################
#####################################
# TARGETS
all:
	@echo "Use one of the following:"
	@echo "make install_mpich2"
	@echo "make install_zlib"
	@echo "make install_hdf5"
	@echo "make install_h5part"
	@echo "make install_metis"
	@echo "make install_makedepf90"
	@echo "make install_ppm"
	@echo "make install_all"
	@echo "make install_allbutmpi"
	@echo "make install_fftw"


install_mpich2		:	mpich2_extract mpich2_configure mpich2_build mpich2_install
install_zlib		:	zlib_extract zlib_configure zlib_build zlib_install	
install_hdf5_serial     :       hdf5_extractserial hdf5_configureserial hdf5_buildserial hdf5_installserial
install_h5part_serial   :       h5part_extractserial h5part_configureserial h5part_buildserial h5part_installserial
install_hdf5		:	hdf5_extract hdf5_configure hdf5_build hdf5_install
install_h5part		:	h5part_extract h5part_configure h5part_build h5part_install
install_metis		:	metis_extract metis_configure metis_build metis_install	
install_ppm		:	ppm_extract ppm_configure ppm_build ppm_install
install_ppmnumerics	:	ppmnumerics_extract ppmnumerics_configure ppmnumerics_build ppmnumerics_install
install_makedepf90	:	makedepf90_extract makedepf90_configure makedepf90_build makedepf90_install
install_all		:	install_mpich2 install_zlib install_hdf5 install_h5part install_metis install_ppm install_makedepf90
install_allbutmpi	:	install_zlib install_hdf5 install_h5part install_metis install_ppm install_makedepf90
install_fftw		:	fftw_extract fftw_configure fftw_build fftw_install fftw_configure2 fftw_build fftw_install


######################################
#######################################
# INCLUDE FILES
include Makefile.mpich2
include Makefile.zlib
include Makefile.hdf5serial
include Makefile.h5partserial
include Makefile.hdf5
include Makefile.h5part
include Makefile.metis
include Makefile.ppm
include Makefile.ppmnumerics
include Makefile.makedepf90
include Makefile.fftw
