      SUBROUTINE colloid_compute_acceleration(this,time,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_acceleration
        !----------------------------------------------------
        !
        ! Purpose     : Compute colloid accelerations,
        !               i.e., translation and rotation.
        !               
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1 19.11.2010, original version.
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
    
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)              :: this
        REAL(MK), INTENT(IN)                    :: time
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: dim,num,i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim   = this%num_dim
        num   = this%num_colloid
        
        
#ifdef __COLLOID_NOACCE
        
        !----------------------------------------------------
        ! No acceleration
        !        
        ! Solution to keep colloid moving/rotating
        ! in a constant velocity, i.e., no translational or
        ! rotational acceleration.
        !----------------------------------------------------
        
        this%f(1:dim,1:num)   = 0.0_MK
        this%alpha(1:3,1:num) = 0.0_MK
        
#else
        
        !----------------------------------------------------
        ! Calculate translating accelerations.
        !----------------------------------------------------
        
        IF( this%translate ) THEN
           
           DO i = 1, dim
              
              this%f(i,1:num) = &
                   this%drag(i,1:num) / this%m(1:num)
              
           END DO
           
        ELSE
           
           this%f(1:dim,1:num) = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! No movement in y-direction.
        !----------------------------------------------------
        
        !IF ( time < 30 ) THEN
        !   this%f(2,1:num) = 0.0_MK
        !END IF
     
        !----------------------------------------------------
        ! Calculate rotating accelerations.
        ! Note that torque and its accleration are 3D.
        !----------------------------------------------------
        
        IF( this%rotate ) THEN
           
           DO i = 1, 3
              
              this%alpha(i,1:num) = &
                   this%torque(i,1:num) / this%mmi(1:num)
              
           END DO
           
        ELSE
           
           this%alpha(1:3,1:num) = 0.0_MK
           
        END IF
        
        
        !----------------------------------------------------
        ! In context of symmetry boundaries, we have to 
        ! eliminate the force on the direction normal to 
        ! symmetry boundaries and torque.
        !----------------------------------------------------
        
        Do i = 1, dim
           
           IF ( this%bcdef(2*i-1) == mcf_bcdef_symmetry .AND.&
                this%bcdef(2*i) == mcf_bcdef_symmetry ) THEN
              
              this%f(i,1:num) = 0.0_MK
              
              this%alpha(1:3,1:num) = 0.0_MK              
              
           END IF
           
        END DO ! i = 1, dim
        
#endif        
        
        
9999    CONTINUE
        
        
        RETURN          
        
      END SUBROUTINE colloid_compute_acceleration
      
      
      
