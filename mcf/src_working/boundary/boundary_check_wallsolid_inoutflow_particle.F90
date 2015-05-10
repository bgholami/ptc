      SUBROUTINE boundary_check_wallsolid_inoutflow_particle(this,&
           p_x,l_w,p_sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : boundary_check_wallsolid_inoutflow_particle
        !----------------------------------------------------
        !
        ! Purpose     : Check if p_x is a solid wall/inflow/outflow 
        !               boundary particle. 
        !       
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 14.12 2009, original version.
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
        
        !----------------------------------------------------
        ! Arguments :
        !
        ! dx        : initial distance between particles.
        ! p_x       : position.
        ! l_w       : indicate inside wall or not.
        ! p_sid     : species ID.
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(Boundary), INTENT(INOUT)           :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: p_x
        LOGICAL, INTENT(OUT)                    :: l_w
        INTEGER, INTENT(OUT)                    :: p_sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: num_dim
        INTEGER                                 :: i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        l_w       = .FALSE.
        num_dim   = this%num_dim
                
        !----------------------------------------------------
        ! Loop over each boundary direction.
        !----------------------------------------------------
        
        DO i = 1, num_dim
           
           IF ( p_x(i) <= this%min_phys(i) .AND. &
                p_x(i) >= this%min_phys_t(i) ) THEN
              
              IF ( this%bcdef(2*i-1) == &
                   mcf_bcdef_wall_solid ) THEN
                 
                 l_w = .TRUE.
                 p_sid = 1-2*i
                 
                 EXIT

              ELSE IF ( (this%bcdef(2*i-1) == mcf_bcdef_inflow) .OR. &
                   (this%bcdef(2*i-1) == mcf_bcdef_outflow) .OR. &
                   ( this%bcdef(2*i-1) == mcf_bcdef_dirichlet) )THEN
                 
                 l_w  = .FALSE.
                 p_sid = mcf_particle_type_fluid
                 
                 EXIT 

              ELSE
                 
                 PRINT *, p_x(1:num_dim), &
                      "boundary_check_wallsolid_inoutflow_particle : ", &
                      "particle is not recognized !"
                 stat_info = -1
                 GOTO 9999
                 
              END IF ! bcdef
              
           ELSE IF ( p_x(i) >= this%max_phys(i) .AND. &
                p_x(i) < this%max_phys_t(i) ) THEN
              
              IF ( this%bcdef(2*i) == &
                   mcf_bcdef_wall_solid ) THEN
                 
                 l_w = .TRUE.
                 p_sid = -2*i
                 
                 EXIT   


              ELSE IF ( (this%bcdef(2*i) == mcf_bcdef_inflow) .OR. &
                   (this%bcdef(2*i) == mcf_bcdef_outflow) .OR. &
                   ( this%bcdef(2*i) == mcf_bcdef_dirichlet) )THEN

                 ! for now fluid
                 l_w  = .FALSE.
                 p_sid = mcf_particle_type_fluid
                                  
                 EXIT  
                 
              ELSE
                 
                 PRINT *, p_x(1:num_dim), &
                      "boundary_check_wallsolid_inoutflow_particle : ", &
                      "particle is not recognized !"
                 stat_info = -1
                 GOTO 9999
                 
              END IF ! bcdef
              
           END IF ! p_x
           
        END DO ! i
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE boundary_check_wallsolid_inoutflow_particle
      
