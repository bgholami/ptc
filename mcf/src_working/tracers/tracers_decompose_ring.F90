      SUBROUTINE tracers_decompose_ring(this,&
           l_map_x,    l_map_v,     &
           l_map_id,   l_map_f,     &
           stat_info)
        !----------------------------------------------------
        ! Program     :  tracers_decompose_ring
        !----------------------------------------------------
        !
        ! Purpose     :  Communicate variables between
        !                neibouring processes.
        !                What to comminicate depends on
        !                the arguments.
        !
        ! Reference   :
        !
        ! Remark      :  Positions are always needed for 
        !                decomposing at the moment.
        !                
        !
        ! Revisions   :  V0.2 02.10.2013, with ppm 1.2.2
        !
        !                V0.1 15.07.2009, original version.
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
        ! Modules used.
        !----------------------------------------------------

        USE ppm_module_map_part_util
        USE ppm_module_map_part
        
      	!----------------------------------------------------
      	!  Arguments
      	!----------------------------------------------------
        
        TYPE(Particles),INTENT(INOUT)   :: this
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_x
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_v
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_id
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_f
        INTEGER, INTENT(OUT)            :: stat_info
        
      	!----------------------------------------------------
      	! Local variables
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        
        
	!----------------------------------------------------
      	! Local/MAPPING
      	!----------------------------------------------------
	
        LOGICAL                         :: map_x
        LOGICAL                         :: map_v
        LOGICAL                         :: map_id
        LOGICAL                         :: map_f
        INTEGER                         :: rank


      	!----------------------------------------------------
      	! Initialize variables
      	!----------------------------------------------------
        
        stat_info       = 0
        stat_info_sub   = 0
        
        map_x     = .FALSE.
        map_v     = .FALSE.
        map_id    = .FALSE.
        map_f     = .FALSE.
        

      	!----------------------------------------------------
        ! Change flag of each variable according to input
        ! arguments. The flags indicate whether or not
        ! the variables needs to be communicated.
      	!----------------------------------------------------
        
        IF( PRESENT(l_map_x) ) THEN	   
           map_x = l_map_x
        END IF
        IF( PRESENT(l_map_v) ) THEN	   
           map_v = l_map_v
        END IF
        IF( PRESENT(l_map_id) ) THEN
           map_id = l_map_id
        END IF
        IF( PRESENT(l_map_f) ) THEN	   
           map_f = l_map_f
        END IF
        
        num_dim = physics_get_num_dim(this%phys,stat_info_sub)
        rank    = technique_get_rank(this%tech,stat_info_sub)
        
	!----------------------------------------------------
      	! Decomposition starts
      	!----------------------------------------------------


	!----------------------------------------------------
      	! Choose mapping
      	!----------------------------------------------------
        
        IF ( map_x ) THEN
           
           CALL ppm_map_part_eqdistrib(this%x, &
                this%num_part_real, stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Starting decomposion x failed!'	   
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

	!----------------------------------------------------
     	! Push the data into buffer
     	!----------------------------------------------------
        
        IF ( map_v ) THEN  

           CALL ppm_map_part_push(this%v, num_dim, &
                this%num_part_real, stat_info_sub)

           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Pushing v failed!'	   
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_id) THEN  

           CALL ppm_map_part_push(this%id, this%num_id, &
                this%num_part_real, stat_info_sub)

           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Pushing id failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( map_f) THEN
           
           CALL ppm_map_part_push(this%f, num_dim, &
                this%num_part_real, stat_info_sub)

           IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Pushing f failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF


	!----------------------------------------------------
     	! Send the buffer
     	!----------------------------------------------------
                
        CALL ppm_map_part_send(this%num_part_real, &
             this%num_part_all, stat_info_sub) 
        
        IF (stat_info_sub /=0 ) THEN
           PRINT *, 'tracers_decompose_ring : ', &
                'Sending data failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
	!----------------------------------------------------
     	! Pop data from the buffer
     	!----------------------------------------------------
    
        IF( map_f) THEN 

           CALL ppm_map_part_pop(this%f, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Poping f failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_id) THEN

           CALL ppm_map_part_pop(this%id, this%num_id, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Poping id failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
   
        IF ( map_v) THEN
            
           CALL ppm_map_part_pop(this%v, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'tracers_decompose_ring : ',&
                   'Poping v failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        

        CALL ppm_map_part_pop(this%x, num_dim, &
             this%num_part_real, &
             this%num_part_all, stat_info_sub)
        
        IF (stat_info_sub /=0 ) THEN
           PRINT *, 'tracers_decompose_ring : ',&
                'Poping x failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Update particles number after mapping
        !----------------------------------------------------
        
        this%num_part_real     = this%num_part_all 
        this%num_part_ghost    = 0
        this%num_part_fluid    = 0
        this%num_part_sym      = 0
        this%num_part_wall_sym = 0
        this%num_part_wall_solid_real  = 0
        this%num_part_wall_solid_ghost = 0
        this%num_part_colloid = 0
        
        
9999    CONTINUE	


        RETURN 

      END SUBROUTINE tracers_decompose_ring

