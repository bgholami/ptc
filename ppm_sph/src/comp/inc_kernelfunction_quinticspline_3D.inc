! a quintic spline kernel calculation module 2D
! see Morris, Fox and Zhu (1997) doi:10.1016/j.jcp.2005.02.001
! Input: dij
normedDist = dij * cutoff_inv * 3._mk

IF (normedDist .LE. 1.0) THEN
  ss3 = (3._mk - normedDist)
  ss2 = (2._mk - normedDist)
  ss1 = (1._mk - normedDist)

  !   w = factorW*(ss3*ss3*ss3*ss3*ss3 - 6._mk*ss2*ss2*ss2*ss2*ss2 + 15._mk*ss1*ss1*ss1*ss1*ss1)
  !   gradW = factorGradW * (ss3*ss3*ss3*ss3 - 6._mk*ss2*ss2*ss2*ss2 + 15._mk*ss1*ss1*ss1*ss1)
   w = factorW*(ss3**5 - 6._mk*ss2**5 + 15._mk*ss1**5)
   gradW = factorGradW * (ss3**4 - 6._mk*ss2**4 + 15._mk*ss1**4)

ELSE IF (normedDist .le. 2.0_MK) THEN
  ss3 = (3._mk - normedDist)
  ss2 = (2._mk - normedDist)

  !   w = factorW * (ss3*ss3*ss3*ss3*ss3 - 6.0_MK*ss2*ss2*ss2*ss2*ss2)
  !   gradW = factorGradW * (ss3*ss3*ss3*ss3 - 6.0_MK*ss2*ss2*ss2*ss2)
   w = factorW * (ss3**5 - 6.0_MK*ss2**5)
   gradW = factorGradW * (ss3**4 - 6.0_MK*ss2**4)

ELSE IF (normedDist .le. 3.0_MK) THEN
  ss3 = (3._mk - normedDist)

  !   w = factorW * ss3*ss3*ss3*ss3*ss3
  !   gradW = factorGradW * ss3*ss3*ss3*ss3
   w = factorW * ss3**5
   gradW = factorGradW * ss3**4

ELSE
   w = 0.0_MK
   gradW = 0.0_MK
ENDIF
