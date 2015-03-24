      SUBROUTINE colloid_integrate_angle(this,dt,accuracy_order,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_integrate_angle
        !----------------------------------------------------
        !
        ! Purpose     : Integrate the angles of rotation.
        !
        ! Remark      : Colloid are modeled as rigid body.
        !
        !
        !  Revision   : V0.1  14.10.2009, original.
        !
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)      :: this
        REAL(MK), INTENT(IN)            :: dt
        INTEGER, INTENT(IN)             :: accuracy_order
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim,i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim           = this%num_dim

#if __POSITION_FIXED
#else        
        IF ( this%rotate ) THEN

           !-------------------------------------------------
           ! Select different accuracy oder:
           ! 1 velocity contribution;
           ! 2 velocity + acceleration contribution.
           !-------------------------------------------------
           
           SELECT CASE (accuracy_order)
              
           CASE(1)
           
              DO i = 1, this%num_colloid
                 
                 this%phi(:,i) = this%omega(:,i) * dt
                 
                 !-------------------------------------------
                 ! Currently only usefull for 2D,
                 ! since 3D roation angle can not be simply
                 ! added up for accumulation.
                 !-------------------------------------------
                 
                 this%theta(:,i) = this%theta(:,i) + &
                      this%phi(:,i)
                 
              END DO
              
           CASE(2)
              
              DO i = 1, this%num_colloid
                 
                 this%phi(:,i) = this%omega(:,i) * dt + &
                      0.5_MK * this%alpha(:,i) * dt**2
                 
                 !-------------------------------------------
                 ! Currently only usefull for 2D,
                 ! since 3D roation angle can not be simply
                 ! added up for accumulation.
                 !-------------------------------------------
             
                 this%theta(:,i) = this%theta(:,i) + &
                      this%phi(:,i)
                 
              END DO

              
           END SELECT ! accuracy_order

           
           !-------------------------------------------------
           ! Compute rotation matrix for each colloid.
           !-------------------------------------------------
           
           
           CALL colloid_compute_rotation_matrix(this,stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              
              PRINT *, "colloid_integrate_angle : ", &
                   "Computing rotaiton matrix failed "
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! rotate
        
#endif
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_integrate_angle
      
