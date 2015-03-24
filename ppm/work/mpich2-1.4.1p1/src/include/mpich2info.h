/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2007 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

/* This file creates strings for the most important configuration options.
   These are then used in the file src/mpi/init/initthread.c to initialize
   global variables that will then be included in both the library and 
   executables, providing a way to determine what version and features of
   MPICH2 were used with a particular library or executable. 
*/
#ifndef MPICH2INFO_H_INCLUDED
#define MPICH2INFO_H_INCLUDED

#define MPICH2_CONFIGURE_ARGS_CLEAN "--prefix=/scratch/ptc/prefix-O2/mpich2/1.4.1p1 FC=ifort FCFLAGS=-O2 CC=icc CFLAGS=-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2 CXX=icpc CXXFLAGS=-O2 --enable-fc --enable-fast=defopt --enable-timer-type=gettimeofday --with-mpe RSHCOMMAND=/usr/bin/ssh --disable-dependencies"
#define MPICH2_VERSION_DATE "Thu Sep  1 13:53:02 CDT 2011"
#define MPICH2_DEVICE "ch3:nemesis"
#define MPICH2_COMPILER_CC "icc -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O2   -O2"
#define MPICH2_COMPILER_CXX "icpc -O2  -O2"
#define MPICH2_COMPILER_F77 "ifort   -O2"
#define MPICH2_COMPILER_FC "ifort -O2    -O2"

#endif
