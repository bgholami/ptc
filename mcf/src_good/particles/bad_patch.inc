  !----------------------------------------------------
  ! interpolation loop
  !----------------------------------------------------       
  
  DO it = 1, dim

     f0 = 0.0d0
     f1 = 0.0d0
     r0 = 0.0d0
     r1 = 0.0d0

     r2 = 0.0d0     

     x_t(1:dim) = this%x(1:dim, i)

     DO jp = 1,num

        !! only fluid SPH particles
        if (this%id(this%sid_idx, jp) == 0) then

           x_p(1:dim) = this%x(1:dim, jp)
           r(1:dim) = x_t(1:dim) - x_p(1:dim)
           r_norm = DOT_PRODUCT(r(1:dim), r(1:dim))

           IF (r_norm <= cut_off2) THEN

              r_norm = SQRT(r_norm)

              CALL kernel_kernel(this%kern, r_norm, w, stat_info_sub)
              rhs_density_type = 2
              IF (rhs_density_type == 1) THEN
                 mass_factor = this%m(jp)
              ELSE
                 mass_factor = 1.0_MK
              END IF

              ! function to be interpolated
              f_value =  this%v(it, jp)
              ! density
              d_value = this%rho(jp) / mass_factor


              !! setup system
              f0 = f0 + W / d_value * f_value
              f1 = f1 + W / d_value * f_value * r(1:dim)

              r0 = r0 + W / d_value
              r1 = r1 + W / d_value * r(1:dim)
              do ii = 1, 2
                 r2(ii,:) = r2(ii,:) + W / d_value * r(ii) * r(1:dim)
              end do
!!!!!!!!!!!!!!!!!!!


           END IF

        end if

     end DO

     !! solve system
     a(1,1) = r0 
     a(1,2) = -r1(1)
     a(1,3) = -r1(2) 
     a(2,1) = r1(1) 
     a(2,2) = -r2(1,1)
     a(2,3) = -r2(2,1)
     a(3,1) = r1(2) 
     a(3,2) = -r2(1,2)
     a(3,3) = -r2(2,2)
     f(1) = f0
     f(2) = f1(1)
     f(3) = f1(2) 


     ! solves to find the first unknown of a system of linear equations of size 3.

     ! [aij] * [xj] = [fi]

     !       a12 a23 f3 - a12 a33 f2 - a13 a22 f3 + a13 a32 f2 + a22 a33 f1 - a23 a32 f1
     ! x1 = ---------------------------------------------------------------------------------
     !      a11 a22 a33 - a11 a23 a32 - a12 a21 a33 + a12 a23 a31 + a13 a21 a32 - a13 a22 a31x


     res1 = a(1,2) * a(2,3) * f(3) - &
          a(1,2) * a(3,3) * f(2) - &
          a(1,3) * a(2,2) * f(3) + &
          a(1,3) * a(3,2) * f(2) + &
          a(2,2) * a(3,3) * f(1) - &
          a(2,3) * a(3,2) * f(1)


     res1 = res1 / ( &
          a(1,1) * a(2,2) * a(3,3) - &
          a(1,1) * a(2,3) * a(3,2) - &
          a(1,2) * a(2,1) * a(3,3) + &
          a(1,2) * a(2,3) * a(3,1) + &
          a(1,3) * a(2,1) * a(3,2) - &
          a(1,3) * a(2,2) * a(3,1))
!!!!!!!!!!!!!!!!!!!!

     this%v(it, i) = res1

  END DO

