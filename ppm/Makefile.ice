METIS_VERSION      = 4.0
MAKEDEPF90_VERSION = 2.8.8
PPM_VERSION        = 1.2.2
ZLIB_VERSION       = 1.2.5
HDF5_VERSION       = 1.8.10-snap0
H5PART_VERSION     = 1.6.6
#MPICH2_VERSION     = 1.4.1p1

PREFIX          = /home/cluster/t7841/lu96mah4/prefix-mpt

METIS_BASE      = ${PREFIX}/metis/${METIS_VERSION}
ZLIB_BASE       = ${PREFIX}/zlib/${ZLIB_VERSION}
HDF5_BASE       = ${PREFIX}/hdf5/${HDF5_VERSION}
H5PART_BASE     = ${PREFIX}/h5part/${H5PART_VERSION}
PPM_BASE        = ${PREFIX}/ppm/${PPM_VERSION}

#MPI_BASE        = ${PREFIX}/mpich2/${MPICH2_VERSION}
#MPI_LIB         = "-L${MPI_BASE}/lib"
#MPI_INC         = "-I${MPI_BASE}/include"
#MPI_BASE = system one used here
#MPI_LIB = system one "-L${MPI_BASE}/lib"
#MPI_INC = system one

MPIFC           = mpif90
MPICC           = mpicc
MPICXX          = mpicxx 

BUILD_MAKE_FLAGS     =    -j 8
