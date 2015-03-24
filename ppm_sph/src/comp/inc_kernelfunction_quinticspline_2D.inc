! a quintic spline kernel calculation module 2D
! see Morris, Fox and Zhu (1997) doi:10.1016/j.jcp.2005.02.001
! Input: dij
normedDist = dij * cutoff_inv * 3._mk

IF (normedDist .LE. 1.0_mk) THEN
  ss3 = (3._mk - normedDist)
  ss2 = (2._mk - normedDist)
  ss1 = (1._mk - normedDist)

   w = factorW*(ss3*ss3*ss3*ss3*ss3 - 6._mk*ss2*ss2*ss2*ss2*ss2 + 15._mk*ss1*ss1*ss1*ss1*ss1)
   gradW = factorGradW * (ss3*ss3*ss3*ss3 - 6._mk*ss2*ss2*ss2*ss2 + 15._mk*ss1*ss1*ss1*ss1)

ELSE IF (normedDist .le. 2.0_MK) THEN
  ss3 = (3._mk - normedDist)
  ss2 = (2._mk - normedDist)

   w = factorW * (ss3*ss3*ss3*ss3*ss3 - 6.0_MK*ss2*ss2*ss2*ss2*ss2)
   gradW = factorGradW * (ss3*ss3*ss3*ss3 - 6.0_MK*ss2*ss2*ss2*ss2)

ELSE IF (normedDist .le. 3.0_MK) THEN
  ss3 = (3._mk - normedDist)

   w = factorW * ss3*ss3*ss3*ss3*ss3
   gradW = factorGradW * ss3*ss3*ss3*ss3

ELSE
   w = 0.0_mk
   gradW = 0.0_mk

ENDIF
