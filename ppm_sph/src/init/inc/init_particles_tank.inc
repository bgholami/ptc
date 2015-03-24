#if __DIM==2
DO j = 1,Np(2)
  y = min_compbox(2) + ((REAL(j,mk)-0.5_mk)*dpx(2))
  DO i=1,Np(1)
    x = min_compbox(1) + ((REAL(i,mk)-0.5_mk)*dpx(1))

    ! left wall of tank
    if (x .le. tankbox_min(1) .and. x .ge. tankbox_min(1)-wallthickness(1)) then
      if (y .le. tankbox_max(2) .and. y .ge. tankbox_min(2)-wallthickness(2)) then
        call place_particle(xpos=x,ypos=y,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
        cycle
      endif
    endif

    ! right wall of tank
    if (x .ge. tankbox_max(1) .and. x .le. tankbox_max(1)+wallthickness(1)) then
      if (y .le. tankbox_max(2) .and. y .ge. tankbox_min(2)-wallthickness(2)) then
        call place_particle(xpos=x,ypos=y,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
        cycle
      endif
    endif

    ! lower wall of tank
    if (x .gt. tankbox_min(1) .and. x .lt. tankbox_max(1)) then
      if (y .le. tankbox_min(2) .and. y .ge. tankbox_min(2)-wallthickness(2)) then
        call place_particle(xpos=x,ypos=y,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
        cycle
      endif
    endif

    ! gate
    if (x .gt. gate_min(1) .and. x .lt. gate_max(1)) then
      if (y .gt. gate_min(2) .and. y .lt. gate_max(2)) then
        call place_particle(xpos=x,ypos=y,pID=gate_phase,Np_loc=Npart_local,created=prtl_created)
        cycle
      endif
    endif

    ! water
    if (x .gt. watercolumn_min(1) .and. x .lt. watercolumn_max(1)) then
      if (y .gt. watercolumn_min(2) .and. y .lt. watercolumn_max(2)) then
        rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK)
        if ( rad .le. bubble1_rad) then
          call place_particle(xpos=x,ypos=y,pID=bubble1_phase,Np_loc=Npart_local,created=prtl_created)
        else
          call place_particle(xpos=x,ypos=y,pID=watercolumn_phase,Np_loc=Npart_local,created=prtl_created)
        endif
        cycle
      endif
    endif


  ENDDO
ENDDO

!<----------------- 3D--------------------------
#else
DO k=1,Np(3)
  z = min_compbox(3) + ((REAL(k,MK)-0.5_MK)*dpx(3))
  DO j = 1,Np(2)
    y = min_compbox(2) + ((REAL(j,mk)-0.5_mk)*dpx(2))
    DO i=1,Np(1)
      x = min_compbox(1) + ((REAL(i,mk)-0.5_mk)*dpx(1))

      ! left wall of tank
      if (x .le. tankbox_min(1) .and. x .ge. tankbox_min(1)-wallthickness(1)) then
        if (y .le. tankbox_max(2)+wallthickness(2) .and. y .ge. tankbox_min(2)-wallthickness(2)) then
          if (z .le. tankbox_max(3) .and. z .ge. tankbox_min(3) - wallthickness(3)) then
            call place_particle(xpos=x,ypos=y,zpos=z,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
            cycle
          endif
        endif
      endif

      ! right wall of tank
      if (x .ge. tankbox_max(1) .and. x .le. tankbox_max(1)+wallthickness(1)) then
        if (y .le. tankbox_max(2)+wallthickness(2) .and. y .ge. tankbox_min(2)-wallthickness(2)) then
          if (z .le. tankbox_max(3) .and. z .ge. tankbox_min(3) - wallthickness(3)) then
            call place_particle(xpos=x,ypos=y,zpos=z,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
            cycle
          endif
        endif
      endif

      ! front wall of tank
      if (x .gt. tankbox_min(1) .and. x .lt. tankbox_max(1)) then
        if (y .le. tankbox_min(2) .and. y .ge. tankbox_min(2)-wallthickness(2)) then
          if (z .le. tankbox_max(3) .and. z .ge. tankbox_min(3) - wallthickness(3)) then
            call place_particle(xpos=x,ypos=y,zpos=z,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
            cycle
          endif
        endif
      endif

      ! rear wall of tank
      if (x .gt. tankbox_min(1) .and. x .lt. tankbox_max(1)) then
        if (y .ge. tankbox_max(2) .and. y .le. tankbox_max(2)+wallthickness(2)) then
          if (z .le. tankbox_max(3) .and. z .ge. tankbox_min(3) - wallthickness(3)) then
            call place_particle(xpos=x,ypos=y,zpos=z,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
            cycle
          endif
        endif
      endif

      ! lower wall of tank
      if (x .gt. tankbox_min(1) .and. x .lt. tankbox_max(1)) then
        if (y .gt. tankbox_min(2) .and. y .lt. tankbox_max(2)) then
          if (z .le. tankbox_min(3) .and. z .gt. tankbox_min(3)-wallthickness(3)) then
            call place_particle(xpos=x,ypos=y,zpos=z,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
            cycle
          endif
        endif
      endif

      ! gate
      if (x .ge. gate_min(1) .and. x .le. gate_max(1)) then
        if (y .ge. gate_min(2) .and. y .le. gate_max(2)) then
          if (z .ge. gate_min(3) .and. z .le. gate_max(3)) then
            call place_particle(xpos=x,ypos=y,zpos=z,pID=gate_phase,Np_loc=Npart_local,created=prtl_created)
            cycle
          endif
        endif
      endif

      ! water
      if (x .gt. watercolumn_min(1) .and. x .lt. watercolumn_max(1)) then
        if (y .gt. watercolumn_min(2) .and. y .lt. watercolumn_max(2)) then
          if (z .gt. watercolumn_min(3) .and. z .lt. watercolumn_max(3)) then
            rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK + (z-bubble1_pos(3))**2._mk)
            if ( rad .le. bubble1_rad) then
              call place_particle(xpos=x,ypos=y,zpos=z,pID=bubble1_phase,Np_loc=Npart_local,created=prtl_created)
            else
              call place_particle(xpos=x,ypos=y,zpos=z,pID=watercolumn_phase,Np_loc=Npart_local,created=prtl_created)
            endif
            cycle
          endif
        endif
      endif

    enddo
  enddo
enddo

#endif
