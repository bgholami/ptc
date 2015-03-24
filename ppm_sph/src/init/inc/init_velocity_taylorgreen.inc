#if __DIM==2
  dxpdt(1,j) = -vel_ref*cos(2._MK*pi/length_ref * xp(1,j))*sin(2._MK*pi/length_ref * xp(2,j))
  dxpdt(2,j) =  vel_ref*sin(2._MK*pi/length_ref * xp(1,j))*cos(2._MK*pi/length_ref * xp(2,j))
#elif __DIM==3

#endif
