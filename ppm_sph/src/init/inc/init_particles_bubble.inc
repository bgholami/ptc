#if __DIM==2
  DO j = 1,Np(2)
    y = min_compbox(2) + ((REAL(j,MK)-0.5_MK)*dpx(2))
    DO i=1,Np(1)
       x = min_compbox(1) + ((REAL(i,MK)-0.5_MK)*dpx(1))
       !-------------------------------------------------
       !  Initialize particles in cube
       !-------------------------------------------------
       rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK)

       !< HACKED
       !< only ring around bubble interface
       !if (abs(rad-bubble1_rad) .gt. 2.*cutoff) cycle

       Npart = Npart + 1
       IF (MOD(Npart,Nproc) .EQ. rank) THEN
          Npart_local = Npart_local + 1
          ! Position
          xp(1,Npart_local) = x
          xp(2,Npart_local) = y
          ! Phase ID
          ap(phaseidx,Npart_local) = 1
          IF (rad .LE. bubble1_rad) then
             if (bubble1_phase .lt. 1) then
                ap(phaseidx,Npart_local) = N_PHASES
             else
                ap(phaseidx,Npart_local) = bubble1_phase
             endif
          ENDIF
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
  z = min_compbox(3) + ((REAL(k,MK)-0.5_MK)*dpx(3))
  DO j=1,Np(2)
    y = min_compbox(2) + ((REAL(j,MK)-0.5_MK)*dpx(2))
    DO i=1,Np(1)
      x = min_compbox(1) + ((REAL(i,MK)-0.5_MK)*dpx(1))

      !-------------------------------------------------
      !  Initialize particles in cube
      !-------------------------------------------------
      !< HACKED
      !< only ring around bubble interface
      rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK + (z-bubble1_pos(3))**2._MK)
      if (abs(rad-bubble1_rad) .gt. 2.*cutoff) cycle

      Npart = Npart + 1
      IF (MOD(Npart,Nproc) .EQ. rank) THEN
        Npart_local = Npart_local + 1
        ! Position
        xp(1,Npart_local) = x
        xp(2,Npart_local) = y
        xp(3,Npart_local) = z
        ! Phase ID
        ap(phaseidx,Npart_local) = 1
        rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK+ (z-bubble1_pos(3))**2._MK)
        IF ( rad .LE. bubble1_rad) then
          if (bubble1_phase .lt. 1) then
            ap(phaseidx,Npart_local) = N_PHASES
          else
            ap(phaseidx,Npart_local) = bubble1_phase
          endif
        ENDIF
        ! Particle ID
        ap(partIDidx,Npart_local) = Npart
        ! Global number of particles
        Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
      ENDIF
    ENDDO
  ENDDO
ENDDO
#endif

