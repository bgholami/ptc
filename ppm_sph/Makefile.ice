#---------------------------------------------------------------------
# CHOICES AT COMPILE-TIME
#--------------------------------------------------------------------
# Dimensions of problem (2D or 3D)
DIM     ?= 3

#Kernel 
# 1 = Quintic spline,  2 = B-Spline, 3 = Cubic-Spline 4=Wendland kernel
KERNEL  ?= 1

ZLIB_VERSION            = 1.2.3
H5PART_VERSION          = 1.6.2
MAKEDEPF90_VERSION      = 2.8.8
METIS_VERSION           = 4.0
PPM_VERSION             = 1.2.1

PREFIX          =/home/cluster/t7841/lu96mah4/prefix-mpt

PPM_BASE        = ${PREFIX}/ppm/$(PPM_VERSION)
MAKEDEPF90_BASE = ${PREFIX}/makedepf90/$(MAKEDEPF90_VERSION)
H5PART_BASE     = $(PREFIX)/h5part/$(H5PART_VERSION)
METIS_BASE      = ${PREFIX}/metis/${METIS_VERSION}
ZLIB_BASE       = $(PREFIX)/zlib/$(ZLIB_VERSION)

FC              = mpif90
LIBS2           = ${HDF5_F90_SHLIB} ${MPI_F90_LIB}

#COMPILEROPTION         = -O0 -ffree-form -cpp -g
#COMPILEROPTION         = -O2 -ffast-math -ffree-form -cpp
#COMPILEROPTION         = -O0 -ffree-form -cpp -g
#COMPILEROPTION         = -O0 -ftz -g -fpp -check all
#SWITCHWARNINGS         =  -Wall -Waliasing -pedantic-errors -Wunderflow
#SWITCHWARNINGS         = -warn all
#COMPILEROPTION         = -O0 -g -ftz -fpp -check all -traceback
#COMPILEROPTION         = -O1 -ftz -fpp -g -check all -traceback
#COMPILEROPTION         = -O1 -fpe0 -fpp -check all -traceback
#COMPILEROPTION         = -O0 -fpp -ftz -check all -traceback -g
#COMPILEROPTION         = -O0 -fpp -fpe0 -check bounds -traceback
#COMPILEROPTION         = -O2 -fpp -g -check all -traceback -ftz
#COMPILEROPTION         = -O0 -fpp -ftz -g -traceback -check all
#COMPILEROPTION         = -O1 -fpp -fpe0 -g -traceback
#COMPILEROPTION         = -O0 -fpp -ftz -fpe0 -g -traceback -check all -heap-array
COMPILEROPTION          = -O2 -fpp -ftz 

