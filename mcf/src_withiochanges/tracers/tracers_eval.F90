!--------------------------------------------------
! Subroutine  : tracers_eval_*
!--------------------------------------------------
!
! Purpose     : evaluates auxiliary functions 
!               needed to calculate forces acting on tracers 
!
! Reference   : Longest et al.
!               2004 Computers and Fluids
!
!               Cherukat and McLaughlin
!               1994 J Fluid Mechanics
!
! Remark      : This does not really mimic object-oriented 
!               ideas as the rest of the code does. Therefore, 
!               it will be modified in future to maintain
!               consistency.
!
! Revisions   : V0.1 16.05.2011, original version.
!
!--------------------------------------------------
! Author      : Babak Gholami
! Contact     : babak.gholami@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------

      SUBROUTINE tracers_eval_sgn(vnt, vnt_i, flag, sgni)
        
        REAL(MK)              , INTENT(IN)  :: vnt
        REAL(MK), DIMENSION(:), INTENT(IN)  :: vnt_i 
        INTEGER               , INTENT(IN)  :: flag 
        REAL(MK), DIMENSION(:), INTENT(OUT) :: sgni
        
        INTEGER                             :: i
        
        IF (flag == 1) THEN

           DO i = 1, SIZE(vnt_i, 1)

              IF ( vnt < 0.0_MK ) THEN
                 sgni(i) = 0.0_MK
              ELSE
                 IF ( vnt_i(i) < 0.0_MK ) THEN
                    sgni(i) = 1.0_MK
                 ELSE
                    sgni(i) = -1.0_MK
                 END IF
              END IF

           END DO

        ELSE

           DO i = 1, SIZE(vnt_i, 1)
              IF ( vnt_i(i) < 0.0_MK ) THEN
                 sgni(i) = 1.0_MK
              ELSE
                 sgni(i) = -1.0_MK
              END IF
           END DO
           
        END IF
        
        
        RETURN
        
      END SUBROUTINE tracers_eval_sgn

     
      REAL(MK) FUNCTION tracers_eval_lift_function(kappa, lambda)
       
        REAL(MK)                    :: kappa
        REAL(MK)                    :: lambda
        REAL(MK)                    :: temp

        
        temp = 0.0_MK

        temp = 1.76310_MK + 0.35610_MK * kappa - 1.18370_MK * kappa**2.0_MK + &
               0.8451630_MK * kappa**3.0_MK

        temp = temp - &
               (3.241390_MK / kappa + 2.6760_MK + 0.82480_MK * kappa - &
                0.46160_MK * kappa**2.0_MK) * lambda

        temp = temp + &
               (1.80810_MK + 0.87960_MK * kappa - 1.90090_MK * kappa**2.0_MK + &
                0.981490_MK * kappa**3.0_MK) * lambda**2.0_MK 


        tracers_eval_lift_function = temp       
        
        
        RETURN
        
      END FUNCTION tracers_eval_lift_function



      SUBROUTINE tracers_interface_sort_r(data,n1,n2,d)
        ! sorts data w.r.t its dimension d (insertion sort)
        
        REAL(MK), DIMENSION(n1,n2)  :: data
        INTEGER                     :: n1, n2, d
        REAL(MK), DIMENSION(n1)     :: temp
        INTEGER                     :: i, j
 
        DO i = 2, n2
           j = i - 1
           temp = data(:,i)
           DO WHILE (j>=1)
              IF (data(d,j)>temp(d)) THEN
                 data(:,j+1) = data(:,j)
                 j = j - 1    
              ELSE
                 GOTO 10
              END IF
           END DO

10         CONTINUE

           data(:,j+1) = temp
        END DO

        !RETURN

      END SUBROUTINE tracers_interface_sort_r



      SUBROUTINE tracers_interface_sort_i(data,d_data,n1,n2)
        ! sorts data and reorders d_data accordingly

        INTEGER,  DIMENSION(n2)     :: data 
        REAL(MK), DIMENSION(n1,n2)  :: d_data  
        INTEGER                     :: n1, n2
        INTEGER                     :: temp
        REAL(MK), DIMENSION(n1)     :: temp_d
        INTEGER                     :: i, j
 
        DO i = 2, n2
           j = i - 1
           temp = data(i)
           temp_d = d_data(:,i)
           DO WHILE (j>=1)
              IF (data(j)>temp) THEN
                 data(j+1) = data(j)
                 d_data(:,j+1) = d_data(:,j)
                 j = j - 1
              ELSE
                 GOTO 10
              END IF
           END DO

10         CONTINUE

           data(j+1) = temp
           d_data(:,j+1) = temp_d
        END DO

        !RETURN

      END SUBROUTINE tracers_interface_sort_i
      
