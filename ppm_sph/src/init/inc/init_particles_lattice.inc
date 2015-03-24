#if __DIM==2
DO j = 1,Np(2)
  y = min_compbox(2) + (j-0.5_mk)*dpx
  DO i=1,Np(1)
    x = min_compbox(1) + (i-0.5_mk)*dpx
    !-------------------------------------------------
    !  Initialize particles in cube
    !-------------------------------------------------
    Npart = Npart + 1

    IF (MOD(Npart,Nproc) .EQ. rank) THEN
      Npart_local = Npart_local + 1
      ! Position
      xp(1,Npart_local) = x
      xp(2,Npart_local) = y
      ! Phase ID
      ap(phaseidx,Npart_local) = 1
      ! Particle ID
      ap(partIDidx,Npart_local) = Npart
      ! Global number of particles
      Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
    ENDIF
  ENDDO
ENDDO
#elif __DIM==3
! The particle volumes are constant
DO k=1,Np(3)
  z = min_compbox(3) + (k-0.5_MK)*dpx
  DO j=1,Np(2)
    y = min_compbox(2) + (j-0.5_MK)*dpx
    DO i=1,Np(1)
      x = min_compbox(1) + (i-0.5_MK)*dpx
      !-------------------------------------------------
      !  Initialize particles in cube
      !-------------------------------------------------
      Npart = Npart + 1

      IF (MOD(Npart,Nproc) .EQ. rank) THEN
        Npart_local = Npart_local + 1
        ! Position
        xp(1,Npart_local) = x
        xp(2,Npart_local) = y
        xp(3,Npart_local) = z
        ! Phase ID
        ap(phaseidx,Npart_local) = 1
        ! Particle ID
        ap(partIDidx,Npart_local) = Npart
        ! Global number of particles
        Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
      ENDIF
    ENDDO
  ENDDO
ENDDO
#endif
