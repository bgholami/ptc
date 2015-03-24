#if __DIM==2
DO j = 1,Np(2)
  y = min_compbox(2) + ((REAL(j,mk)-0.5_mk)*dpx(2))
  DO i=1,Np(1)
    x = min_compbox(1) + ((REAL(i,mk)-0.5_mk)*dpx(1))

    !< tank
    if ( (x .ge. tankbox_min(1)-wallthickness(1) .and. x .le. tankbox_max(1)+wallthickness(1)) .and. &
      (y .ge. tankbox_min(2)-wallthickness(2) .and. y .le. tankbox_max(2)+wallthickness(2))) then
      if ( (x .gt. tankbox_min(1) .and. x .lt. tankbox_max(1)) .and. &
        (y .gt. tankbox_min(2) .and. y .lt. tankbox_max(2))) then
        ! nothing
      else
        !< wall
        call place_particle(xpos=x,ypos=y,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif

    !< bubble 2
    rad = sqrt((x-bubble2_pos(1))**2._MK + (y-bubble2_pos(2))**2._MK)
    if ( rad .le. bubble2_rad) then
      call place_particle(xpos=x,ypos=y,pID=bubble2_phase,Np_loc=Npart_local,created=prtl_created)
      if (prtl_created) cycle
    endif

    !< body 1
    xm = 0.5_mk*(body1_min(1)+body1_max(1))
    ym = 0.5_mk*(body1_min(2)+body1_max(2))
    xsize = body1_max(1)-body1_min(1)
    ysize = body1_max(2)-body1_min(2)
    xbody = (x-xm)*cos(body1_angle/180._mk*pi) - (y-ym)*sin(body1_angle/180._mk*pi)
    ybody = (x-xm)*sin(body1_angle/180._mk*pi) + (y-ym)*cos(body1_angle/180._mk*pi)
    if (xbody .gt. -0.5_mk*xsize .and. xbody .lt. 0.5_mk*xsize) then
      if (ybody .gt. -0.5_mk*ysize .and. ybody .lt. 0.5_mk*ysize) then
        call place_particle(xpos=x,ypos=y,pID=body1_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif
    !< body 2
    xm = 0.5_mk*(body2_min(1)+body2_max(1))
    ym = 0.5_mk*(body2_min(2)+body2_max(2))
    xsize = body2_max(1)-body2_min(1)
    ysize = body2_max(2)-body2_min(2)
    xbody = (x-xm)*cos(body2_angle/180._mk*pi) - (y-ym)*sin(body2_angle/180._mk*pi)
    ybody = (x-xm)*sin(body2_angle/180._mk*pi) + (y-ym)*cos(body2_angle/180._mk*pi)
    if (xbody .gt. -0.5_mk*xsize .and. xbody .lt. 0.5_mk*xsize) then
      if (ybody .gt. -0.5_mk*ysize .and. ybody .lt. 0.5_mk*ysize) then
        call place_particle(xpos=x,ypos=y,pID=body2_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif
    !< body 3
    xm = 0.5_mk*(body3_min(1)+body3_max(1))
    ym = 0.5_mk*(body3_min(2)+body3_max(2))
    xsize = body3_max(1)-body3_min(1)
    ysize = body3_max(2)-body3_min(2)
    xbody = (x-xm)*cos(body3_angle/180._mk*pi) - (y-ym)*sin(body3_angle/180._mk*pi)
    ybody = (x-xm)*sin(body3_angle/180._mk*pi) + (y-ym)*cos(body3_angle/180._mk*pi)
    if (xbody .gt. -0.5_mk*xsize .and. xbody .lt. 0.5_mk*xsize) then
      if (ybody .gt. -0.5_mk*ysize .and. ybody .lt. 0.5_mk*ysize) then
        call place_particle(xpos=x,ypos=y,pID=body3_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif
    !< body 4
    xm = 0.5_mk*(body4_min(1)+body4_max(1))
    ym = 0.5_mk*(body4_min(2)+body4_max(2))
    xsize = body4_max(1)-body4_min(1)
    ysize = body4_max(2)-body4_min(2)
    xbody = (x-xm)*cos(body4_angle/180._mk*pi) - (y-ym)*sin(body4_angle/180._mk*pi)
    ybody = (x-xm)*sin(body4_angle/180._mk*pi) + (y-ym)*cos(body4_angle/180._mk*pi)
    if (xbody .gt. -0.5_mk*xsize .and. xbody .lt. 0.5_mk*xsize) then
      if (ybody .gt. -0.5_mk*ysize .and. ybody .lt. 0.5_mk*ysize) then
        call place_particle(xpos=x,ypos=y,pID=body4_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif
    !< body 5
    xm = 0.5_mk*(body5_min(1)+body5_max(1))
    ym = 0.5_mk*(body5_min(2)+body5_max(2))
    xsize = body5_max(1)-body5_min(1)
    ysize = body5_max(2)-body5_min(2)
    xbody = (x-xm)*cos(body5_angle/180._mk*pi) - (y-ym)*sin(body5_angle/180._mk*pi)
    ybody = (x-xm)*sin(body5_angle/180._mk*pi) + (y-ym)*cos(body5_angle/180._mk*pi)
    if (xbody .gt. -0.5_mk*xsize .and. xbody .lt. 0.5_mk*xsize) then
      if (ybody .gt. -0.5_mk*ysize .and. ybody .lt. 0.5_mk*ysize) then
        call place_particle(xpos=x,ypos=y,pID=body5_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif

    ! water
    if (x .ge. watercolumn_min(1) .and. x .le. watercolumn_max(1)) then
      if (y .ge. watercolumn_min(2) .and. y .le. watercolumn_max(2)) then
        rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK)
        if ( rad .le. bubble1_rad) then
          call place_particle(xpos=x,ypos=y,pID=bubble1_phase,Np_loc=Npart_local,created=prtl_created)
        else
          call place_particle(xpos=x,ypos=y,pID=watercolumn_phase,Np_loc=Npart_local,created=prtl_created)
        endif
        if (prtl_created) cycle
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

    !< tank
    if ( (x .ge. tankbox_min(1)-wallthickness(1) .and. x .le. tankbox_max(1)+wallthickness(1)) .and. &
      (y .ge. tankbox_min(2)-wallthickness(2) .and. y .le. tankbox_max(2)+wallthickness(2)) .and. &
      (z .ge. tankbox_min(3)-wallthickness(3) .and. z .le. tankbox_max(3)+wallthickness(3))) then
      if ( (x .gt. tankbox_min(1) .and. x .lt. tankbox_max(1)) .and. &
        (y .gt. tankbox_min(2) .and. y .lt. tankbox_max(2)) .and. &
        (z .gt. tankbox_min(3) .and. z .lt. tankbox_max(3))) then
        ! nothing
      else
        !< wall
        call place_particle(xpos=x,ypos=y,zpos=z,pID=tankbox_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
    endif


    !< body 1
    if (x .gt. body1_min(1) .and. x .lt. body1_max(1)) then
      if (y .gt. body1_min(2) .and. y .lt. body1_max(2)) then
        if (z .gt. body1_min(3) .and. z .lt. body1_max(3)) then
        call place_particle(xpos=x,ypos=y,zpos=z,pID=body1_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
        endif
      endif
    endif
    !< body 2
    if (x .gt. body2_min(1) .and. x .lt. body2_max(1)) then
      if (y .gt. body2_min(2) .and. y .lt. body2_max(2)) then
      if (z .gt. body2_min(3) .and. z .lt. body2_max(3)) then
        call place_particle(xpos=x,ypos=y,zpos=z,pID=body2_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
        endif
      endif
    endif
    !< body 3
    if (x .gt. body3_min(1) .and. x .lt. body3_max(1)) then
      if (y .gt. body3_min(2) .and. y .lt. body3_max(2)) then
      if (z .gt. body3_min(3) .and. z .lt. body3_max(3)) then
        call place_particle(xpos=x,ypos=y,zpos=z,pID=body3_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
      endif
    endif
    !< body 4
    if (x .gt. body4_min(1) .and. x .lt. body4_max(1)) then
      if (y .gt. body4_min(2) .and. y .lt. body4_max(2)) then
      if (z .gt. body4_min(3) .and. z .lt. body4_max(3)) then
        call place_particle(xpos=x,ypos=y,zpos=z,pID=body4_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
      endif
    endif
    !< body 5
    if (x .gt. body5_min(1) .and. x .lt. body5_max(1)) then
      if (y .gt. body5_min(2) .and. y .lt. body5_max(2)) then
      if (z .gt. body5_min(3) .and. z .lt. body5_max(3)) then
        call place_particle(xpos=x,ypos=y,zpos=z,pID=body5_phase,Np_loc=Npart_local,created=prtl_created)
        if (prtl_created) cycle
      endif
      endif
    endif

    ! water
    if (x .ge. watercolumn_min(1) .and. x .le. watercolumn_max(1)) then
      if (y .ge. watercolumn_min(2) .and. y .le. watercolumn_max(2)) then
      if (z .ge. watercolumn_min(3) .and. z .le. watercolumn_max(3)) then
        rad = sqrt((x-bubble1_pos(1))**2._MK + (y-bubble1_pos(2))**2._MK + (z-bubble1_pos(2))**2._mk)
        if ( rad .le. bubble1_rad) then
          call place_particle(xpos=x,ypos=y,zpos=z,pID=bubble1_phase,Np_loc=Npart_local,created=prtl_created)
        else
          call place_particle(xpos=x,ypos=y,zpos=z,pID=watercolumn_phase,Np_loc=Npart_local,created=prtl_created)
        endif
        if (prtl_created) cycle
      endif
      endif
    endif


    enddo
  enddo
enddo

#endif
