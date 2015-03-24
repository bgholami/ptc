do ip = 1, Mpart
  !< Store ip-values to avoid nneigh-times array access
  xip(1) = xp(1,ip)
  xip(2) = xp(2,ip)
#if __DIM==3
  xip(3) = xp(3,ip)
#endif

  do nneigh = 1, verletlist(1)%nvlist(ip)
    jp = verletlist(1)%vlist(nneigh,ip)

    !< Calculate dist-vector and vel-diff
    dx(1) = xp(1,jp) - xip(1)
    dx(2) = xp(2,jp) - xip(2)
#if __DIM==3
    dx(3) = xp(3,jp) - xip(3)
#endif

#if __DIM==2
    dij = sqrt(dx(1)*dx(1)+dx(2)*dx(2))
#elif __DIM==3
    dij = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
#endif
    rdij = 1._mk / (dij + myeps)

#include "inc_kernelfunction.inc"

    !< number density
    drhodt(1,ip) = drhodt(1,ip) + w
#ifdef __SYMMETRY
    drhodt(1,jp) = drhodt(1,jp) + w
#endif
  enddo

enddo
