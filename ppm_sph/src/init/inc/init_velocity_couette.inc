#if __DIM==2
     !dxpdt(1,j) = vel_ref/length_ref*(xp(2,j) - 0.5_MK*(max_phys(2)+min_phys(2)))
     dxpdt(1,j) = vel_ref/length_ref*(xp(2,j) - 0.5_MK*(max_compbox(2)+min_compbox(2)))
#elif __DIM==3
     !dxpdt(1,j) = vel_ref/length_ref*(xp(3,j) - 0.5_MK*(max_compbox(3)+min_compbox(3)))
     dxpdt(1,j) = vel_ref/length_ref*(xp(3,j) - 0.5_MK*(max_compbox(3)+min_compbox(3)))
#endif
