#if __DIM==2
! Beam
nx_beam = INT( beam_size(1)/dpx(1) + 0.5_MK)
ny_beam = INT( 0.5_MK*(beam_size(2)/dpx(2) + 0.5_MK)) * 2 ! Ensure even ny_beam
DO i = 1, nx_beam
  x = beam_pos(1) + (i-1)*dpx(1)
  DO j = 1, ny_beam
    y = beam_pos(2) - 0.5_MK*(ny_beam-1)*dpx(2) + (j-1)*dpx(2)

    !-------------------------------------------------
    !  Initialize particles
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

! Clamp beam part
nx_beam = INT( beam_support/dpx(1) + 0.5_MK)
DO i = 1, nx_beam
  x = beam_pos(1) - i*dpx(1)
  DO j = 1, ny_beam
    y = beam_pos(2) - 0.5_MK*(ny_beam-1)*dpx(2) + (j-1)*dpx(2)
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
      ! Fix clamp
      if (i .gt. nx_beam - 3) pdata(fixidx,Npart_local) = param_special_fix

    ENDIF
  ENDDO
ENDDO

! Clamp wall part
ny_clamp = INT(clamp_thickness/dpx(1) + 0.5_MK)
DO i = 1, nx_beam
  x = beam_pos(1) - i*dpx(1)
  DO j = 1, ny_clamp
    y = beam_pos(2) - 0.5_MK*(ny_beam-1)*dpx(2) - j*dpx(2)
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
      ! Fix clamp
      pdata(fixidx,Npart_local) = param_special_fix
    ENDIF
  ENDDO
  DO j = 1, ny_clamp
    y = beam_pos(2) + 0.5_MK*(ny_beam-1)*dpx(2) + (j)*dpx(2)
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
      ! Fix clamp
      pdata(fixidx,Npart_local) = param_special_fix
    ENDIF
  ENDDO
ENDDO


#else
call abortmessage(myname,'beam not yet implemented for __DIM=3')
#endif
