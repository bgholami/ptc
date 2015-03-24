! a cubic spline kernel calculation module 2D
! see Gray, Monaghan and Swegle (2001)
! Input: dij
normedDist = dij * cutoff_inv * ratio_cutoff_smoothing

IF (normedDist .LE. 1.0) THEN
  ss3 = normedDist**3._mk
  ss2 = normedDist**2._mk
  w = factorW*(1._mk - 1.5_mk*ss2 + 0.75_mk*ss3)
  gradW = factorGradW*(-3._mk*normedDist + 9._mk*0.25_mk*ss2)

ELSE IF (normedDist .le. 2.0_MK) THEN
  w = factorW * 0.25_mk*(2._mk - normedDist)**3._mk
  gradW = -factorGradW*0.25_mk*(2._mk - normedDist)**2._mk

ELSE
   w = 0.0_MK
   gradW = 0.0_MK
ENDIF
