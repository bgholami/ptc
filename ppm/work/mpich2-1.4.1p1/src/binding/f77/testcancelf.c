/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface 
 * DO NOT EDIT
 */
#include "mpi_fortimpl.h"


/* Begin MPI profiling block */
#if defined(USE_WEAK_SYMBOLS) && !defined(USE_ONLY_MPI_NAMES) 
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_TEST_CANCELLED( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled_( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak MPI_TEST_CANCELLED = PMPI_TEST_CANCELLED
#pragma weak mpi_test_cancelled__ = PMPI_TEST_CANCELLED
#pragma weak mpi_test_cancelled_ = PMPI_TEST_CANCELLED
#pragma weak mpi_test_cancelled = PMPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_TEST_CANCELLED = pmpi_test_cancelled__
#pragma weak mpi_test_cancelled__ = pmpi_test_cancelled__
#pragma weak mpi_test_cancelled_ = pmpi_test_cancelled__
#pragma weak mpi_test_cancelled = pmpi_test_cancelled__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_TEST_CANCELLED = pmpi_test_cancelled_
#pragma weak mpi_test_cancelled__ = pmpi_test_cancelled_
#pragma weak mpi_test_cancelled_ = pmpi_test_cancelled_
#pragma weak mpi_test_cancelled = pmpi_test_cancelled_
#else
#pragma weak MPI_TEST_CANCELLED = pmpi_test_cancelled
#pragma weak mpi_test_cancelled__ = pmpi_test_cancelled
#pragma weak mpi_test_cancelled_ = pmpi_test_cancelled
#pragma weak mpi_test_cancelled = pmpi_test_cancelled
#endif



#elif defined(HAVE_PRAGMA_WEAK)

#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_TEST_CANCELLED( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak MPI_TEST_CANCELLED = PMPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled__( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_test_cancelled__ = pmpi_test_cancelled__
#elif !defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_test_cancelled = pmpi_test_cancelled
#else
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled_( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_test_cancelled_ = pmpi_test_cancelled_
#endif

#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#if defined(F77_NAME_UPPER)
#pragma _HP_SECONDARY_DEF PMPI_TEST_CANCELLED  MPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _HP_SECONDARY_DEF pmpi_test_cancelled__  mpi_test_cancelled__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _HP_SECONDARY_DEF pmpi_test_cancelled  mpi_test_cancelled
#else
#pragma _HP_SECONDARY_DEF pmpi_test_cancelled_  mpi_test_cancelled_
#endif

#elif defined(HAVE_PRAGMA_CRI_DUP)
#if defined(F77_NAME_UPPER)
#pragma _CRI duplicate MPI_TEST_CANCELLED as PMPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _CRI duplicate mpi_test_cancelled__ as pmpi_test_cancelled__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _CRI duplicate mpi_test_cancelled as pmpi_test_cancelled
#else
#pragma _CRI duplicate mpi_test_cancelled_ as pmpi_test_cancelled_
#endif
#endif /* HAVE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */
/* End MPI profiling block */


/* These definitions are used only for generating the Fortran wrappers */
#if defined(USE_WEAK_SYMBOLS) && defined(HAVE_MULTIPLE_PRAGMA_WEAK) && \
    defined(USE_ONLY_MPI_NAMES)
extern FORT_DLL_SPEC void FORT_CALL MPI_TEST_CANCELLED( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled_( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak mpi_test_cancelled__ = MPI_TEST_CANCELLED
#pragma weak mpi_test_cancelled_ = MPI_TEST_CANCELLED
#pragma weak mpi_test_cancelled = MPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_TEST_CANCELLED = mpi_test_cancelled__
#pragma weak mpi_test_cancelled_ = mpi_test_cancelled__
#pragma weak mpi_test_cancelled = mpi_test_cancelled__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_TEST_CANCELLED = mpi_test_cancelled_
#pragma weak mpi_test_cancelled__ = mpi_test_cancelled_
#pragma weak mpi_test_cancelled = mpi_test_cancelled_
#else
#pragma weak MPI_TEST_CANCELLED = mpi_test_cancelled
#pragma weak mpi_test_cancelled__ = mpi_test_cancelled
#pragma weak mpi_test_cancelled_ = mpi_test_cancelled
#endif

#endif

/* Map the name to the correct form */
#ifndef MPICH_MPI_FROM_PMPI
#if defined(USE_WEAK_SYMBOLS) && defined(HAVE_MULTIPLE_PRAGMA_WEAK)
/* Define the weak versions of the PMPI routine*/
#ifndef F77_NAME_UPPER
extern FORT_DLL_SPEC void FORT_CALL PMPI_TEST_CANCELLED( MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_2USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_test_cancelled__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_test_cancelled_( MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER
extern FORT_DLL_SPEC void FORT_CALL pmpi_test_cancelled( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif

#if defined(F77_NAME_UPPER)
#pragma weak pmpi_test_cancelled__ = PMPI_TEST_CANCELLED
#pragma weak pmpi_test_cancelled_ = PMPI_TEST_CANCELLED
#pragma weak pmpi_test_cancelled = PMPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak PMPI_TEST_CANCELLED = pmpi_test_cancelled__
#pragma weak pmpi_test_cancelled_ = pmpi_test_cancelled__
#pragma weak pmpi_test_cancelled = pmpi_test_cancelled__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak PMPI_TEST_CANCELLED = pmpi_test_cancelled_
#pragma weak pmpi_test_cancelled__ = pmpi_test_cancelled_
#pragma weak pmpi_test_cancelled = pmpi_test_cancelled_
#else
#pragma weak PMPI_TEST_CANCELLED = pmpi_test_cancelled
#pragma weak pmpi_test_cancelled__ = pmpi_test_cancelled
#pragma weak pmpi_test_cancelled_ = pmpi_test_cancelled
#endif /* Test on name mapping */
#endif /* Use multiple pragma weak */

#ifdef F77_NAME_UPPER
#define mpi_test_cancelled_ PMPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_test_cancelled_ pmpi_test_cancelled__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_test_cancelled_ pmpi_test_cancelled
#else
#define mpi_test_cancelled_ pmpi_test_cancelled_
#endif /* Test on name mapping */

/* This defines the routine that we call, which must be the PMPI version
   since we're renaming the Fortran entry as the pmpi version.  The MPI name
   must be undefined first to prevent any conflicts with previous renamings,
   such as those put in place by the globus device when it is building on
   top of a vendor MPI. */
#undef MPI_Test_cancelled
#define MPI_Test_cancelled PMPI_Test_cancelled 

#else

#ifdef F77_NAME_UPPER
#define mpi_test_cancelled_ MPI_TEST_CANCELLED
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_test_cancelled_ mpi_test_cancelled__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_test_cancelled_ mpi_test_cancelled
/* Else leave name alone */
#endif


#endif /* MPICH_MPI_FROM_PMPI */

/* Prototypes for the Fortran interfaces */
#include "fproto.h"
FORT_DLL_SPEC void FORT_CALL mpi_test_cancelled_ ( MPI_Fint *v1, MPI_Fint *v2, MPI_Fint *ierr ){
    int l2;
    *ierr = MPI_Test_cancelled( (MPI_Status *)(v1), &l2 );
    *v2 = MPIR_TO_FLOG(l2);
}
