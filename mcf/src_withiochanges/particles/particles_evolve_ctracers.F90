      SUBROUTINE particles_evolve_ctracers(this, &
           num,dt,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_evolve_ctracers
        !----------------------------------------------------
        !
        ! Purpose     : evolve c_tracers on SPH particles
        !
        ! Reference   :
        !
        ! Remark      : 
        !
        !
        ! Revision    : V0.1  24.05.2013, original version.
        !               (from particles_integrate_position)
        !
        !----------------------------------------------------
        ! Author      : Babak Gholami
        ! Contact     : babak.gholami@aer.mw.tum.de
        ! 
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.        
        !----------------------------------------------------
        
        !----------------------------------------------------
        !  Arguments
        !
        !  this           : an object of Particles Class.
        !  num            : number of particles needed to be updated,
        !                   i.e. first num particles in this%x 
        !                   are operated.
        !  dt             : time step.
        !  stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        REAL(MK),INTENT(IN)                     :: dt
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim           = &
             physics_get_num_dim(this%phys,stat_info_sub)

        
        DO i = 1, num
              
           IF ( this%id(this%sid_idx,i) == mcf_particle_type_fluid) THEN
              
              this%c_tracers(i) = &
                   this%c_tracers(i) + &
                   this%dn(i) * dt
           END IF
           
        END DO
           
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_evolve_ctracers
