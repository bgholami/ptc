#if __DIM==2
akl =1.875
akx = akl/(beam_size(1))
cosk = cos(akl)
sink = sin(akl)
coshk = cosh(akl)
sinhk = sinh(akl)

! Gray, Monaghan, Swegle
fk1 = cosk + coshk
fk2 = sink - sinhk
denom = fk1*(cosh(akl)-cos(akl)) + fk2*(sinh(akl) - sin(akl))
    ddx = xp(1,j) - beam_pos(1)
    IF (ddx .LT. 0._MK) CYCLE
    ak = akx*ddx
    aax = fk1*(cosh(ak)-cos(ak))+fk2*(sinh(ak)-sin(ak))
    dxpdt(2,j) = vel_ref*aax/denom*pdata(fixidx,j)
        
#elif __DIM==3
   dxpdt(1,j) = vel_ref
   dxpdt(2,j) = vel_ref
   dxpdt(3,j) = vel_ref
#endif

   
