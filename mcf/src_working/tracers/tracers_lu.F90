!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
!MODULE tracers_LU

!CONTAINS

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
Subroutine LUDCMP(A,N,INDX,D,CODE) 

  IMPLICIT NONE

  INTEGER, PARAMETER       :: NMAX = 100
  REAL(8), PARAMETER       :: TINY = 1.5D-16
  REAL(8)                  :: AMAX,DUM, SUM
  REAL(8), DIMENSION(:, :) :: A(N,N)
  REAL(8), DIMENSION(:)    :: VV(NMAX)
  INTEGER CODE, D
  INTEGER, DIMENSION(:)    :: INDX(N)

  ! I add
  INTEGER I, J, K, N, IMAX 

  D=1; CODE=0

  DO I=1,N
     AMAX=0.d0
     DO J=1,N
        IF (DABS(A(I,J)).GT.AMAX)  AMAX=DABS(A(I,J))
     END DO ! j loop
     IF(AMAX.LT.TINY) THEN
        CODE = 1
        RETURN
     END IF
     VV(I) = 1.d0 / AMAX
  END DO ! i loop

  DO J=1,N
     DO I=1,J-1
        SUM = A(I,J)
        DO K=1,I-1
           SUM = SUM - A(I,K)*A(K,J) 
        END DO ! k loop
        A(I,J) = SUM
     END DO ! i loop
     AMAX = 0.d0
     DO I=J,N
        SUM = A(I,J)
        DO K=1,J-1
           SUM = SUM - A(I,K)*A(K,J) 
        END DO ! k loop
        A(I,J) = SUM
        DUM = VV(I)*DABS(SUM)
        IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
        END IF
     END DO ! i loop  

     IF(J.NE.IMAX) THEN
        DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
        END DO ! k loop
        D = -D
        VV(IMAX) = VV(J)
     END IF

     INDX(J) = IMAX
     IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

     IF(J.NE.N) THEN
        DUM = 1.d0 / A(J,J)
        DO I=J+1,N
           A(I,J) = A(I,J)*DUM
        END DO ! i loop
     END IF
  END DO ! j loop

  RETURN
END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
Subroutine LUBKSB(A,N,INDX,B) 

  IMPLICIT NONE

  REAL(8), DIMENSION(:, :) :: A(N,N) 
  REAL(8), DIMENSION(:)    :: B(N)
  REAL(8)                  :: SUM
  INTEGER, DIMENSION(:)    :: INDX(N)
  INTEGER                  :: N

  ! I add
  INTEGER I, J, LL, II

  II = 0

  DO I=1,N
     LL = INDX(I)
     SUM = B(LL)
     B(LL) = B(I)
     IF(II.NE.0) THEN
        DO J=II,I-1
           SUM = SUM - A(I,J)*B(J)
        END DO ! j loop
     ELSE IF(SUM.NE.0.d0) THEN
        II = I
     END IF
     B(I) = SUM
  END DO ! i loop

  DO I=N,1,-1
     SUM = B(I)
     IF(I < N) THEN
        DO J=I+1,N
           SUM = SUM - A(I,J)*B(J)
        END DO ! j loop
     END IF
     B(I) = SUM / A(I,I)
  END DO ! i loop

  RETURN
END subroutine LUBKSB

!END MODULE tracers_LU

! end of file lu.f90
