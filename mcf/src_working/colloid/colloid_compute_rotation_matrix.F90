      SUBROUTINE colloid_compute_rotation_matrix(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_rotation_matrix
        !----------------------------------------------------
        !
        ! Purpose     : Compute rotation matrix for each
        !               colloid.
        !
        ! Referecen   : Chen et al. 
        !               Physics of Fluids, 18, 103605, 2006.
        !
        ! Revision    : V.01  23.08.2010
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
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim,i
        REAL(MK), DIMENSION(3)          :: phi
        REAL(MK)                        :: len
        REAL(MK), DIMENSION(3,3)        :: rot_matrix
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0        
        dim           = this%num_dim
        
        
        IF ( this%rotate ) THEN
           
           DO i = 1, this%num_colloid
           
              phi(1:3) = this%phi(1:3,i)
              
              !----------------------------------------------
              ! Calculate rotated angle from time t to t+dt.
              !----------------------------------------------
              
              len = SQRT(DOT_PRODUCT(phi(1:3),phi(1:3)))
           
              !----------------------------------------------
              ! Make the rotation angle unit vecotr.
              !----------------------------------------------
              
              IF ( len < mcf_machine_zero ) THEN
                 len      = 0.0_MK
                 phi(1)   = 1.0_MK
                 phi(2:3) = 0.0_MK
              ELSE
                 phi(1:3) = phi(1:3) / len
              END IF
           
              rot_matrix(1:dim,1:dim) = 0.0_MK
              
              CALL tool_rotation_matrix(this%tool,&
                   dim, phi(1:3),len,&
                   rot_matrix(1:dim,1:dim),stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "colloid_compute_rotation_matrix ; ", &
                      "Using tool_rotation_matrix failed ! "
                 stat_info = -1
                 GOTO 9999                 
              END IF
              
              this%rot_matrix(1:dim,1:dim,i) = &
                   rot_matrix(1:dim,1:dim)
              
           END DO ! i = 1, num_colloid
           
        END IF ! rotate
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_rotation_matrix
      
