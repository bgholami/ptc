! a wendland quartic kernel calculation module 2D
! Input: dij
!normedDist = dij / smoothinglength
normedDist = dij /smoothinglength !cutoff_inv*ratio_cutoff_smoothing

ss3 =  normedDist*normedDist*normedDist
ss2 =  normedDist*normedDist

IF (normedDist .LT. 2.0_MK) THEN

  w = factorW*(1._MK + 2._MK*normedDist)*((2._MK-normedDist)**4._MK)

  gradW = factorGradW*( 2._MK*(2._MK - normedDist)**4._MK - &
  4._MK*(1._MK+2._MK*normedDist)*(2._MK-normedDist)**3._MK)

ELSE
  w = 0.0_MK
  gradW = 0.0_MK
ENDIF
