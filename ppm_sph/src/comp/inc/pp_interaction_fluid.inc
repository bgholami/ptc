do ip = 1, Mpart
  !< Avoid array access nneigh-times for ip
  xip(1) = xp(1,ip)
  xip(2) = xp(2,ip)
  vip(1) = dxpdt(1,ip)
  vip(2) = dxpdt(2,ip)
#if __DIM==3  
  xip(3) = xp(3,ip)
  vip(3) = dxpdt(3,ip)
#endif

  rhoip = pdata(rhoidx,ip)
  volip = pdata(volidx,ip)
  pressip = pdata(pressidx,ip)
  etaip   = phase(ap(phaseidx,ip))%visc_dyn

  !< All neighbours of ip
  do nneigh = 1, verletlist(1)%nvlist(ip)
    jp = verletlist(1)%vlist(nneigh,ip)

    !< Calculate dist-vector and vel-diff
    dx(1) = xp(1,jp) - xip(1)
    dx(2) = xp(2,jp) - xip(2)
    dvel(1) = vip(1) - dxpdt(1,jp)
    dvel(2) = vip(2) - dxpdt(2,jp)
#if __DIM==3    
    dx(3) = xp(3,jp) - xip(3)
    dvel(3) = vip(3) - dxpdt(3,jp)
#endif
    
#if __DIM==2
    dij = sqrt(dx(1)*dx(1)+dx(2)*dx(2))
#elif __DIM==3
    dij = sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
#endif
    rdij = 1._mk / (dij)

#include "inc_kernelfunction.inc"

    !< Pressure term
    rhojp = pdata(rhoidx,jp)
    voljp = pdata(volidx,jp)
    pressjp = pdata(pressidx,jp)
    !< Volume term
    Vi2Vj2 = volip*volip + voljp*voljp

    !< Pressure force
    part_dens_ip = rhoip / (rhoip + rhojp)
    part_dens_jp = 1._MK - part_dens_ip
    p_avg = part_dens_jp*pressip + part_dens_ip*pressjp
    pforce(1) = Vi2Vj2*p_avg*gradW*dx(1)*rdij
    pforce(2) = Vi2Vj2*p_avg*gradW*dx(2)*rdij
#if __DIM==3
    pforce(3) = Vi2Vj2*p_avg*gradW*dx(3)*rdij
#endif

    !< Viscous term
    eta_avg = etaip
    pforce(1) = pforce(1) + Vi2Vj2*eta_avg*rdij*gradW*dvel(1)
    pforce(2) = pforce(2) + Vi2Vj2*eta_avg*rdij*gradW*dvel(2)
#if __DIM==3
    pforce(3) = pforce(3) + Vi2Vj2*eta_avg*rdij*gradW*dvel(3)
#endif

    !< Other forces....



    !< Store interaction-force to pForce_ip -> add to pForce_total after all pp
    pForce_total(1,ip) = pForce_total(1,ip) + pforce(1)
    pForce_total(2,ip) = pForce_total(2,ip) + pForce(2)
#if __DIM==3
    pForce_total(3,ip) = pForce_total(3,ip) + pForce(3)
#endif
#ifdef __SYMMETRY
    pForce_total(1,jp) = pForce_total(1,jp) - pforce(1)
    pForce_total(2,jp) = pForce_total(2,jp) - pforce(2)
#if __DIM==3
    pForce_total(3,jp) = pForce_total(3,jp) - pforce(3)
#endif
#endif


  enddo
enddo
