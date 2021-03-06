      SUBROUTINE particles_map_ghost_put(this,&
           l_map_x,    l_map_v,    &
           l_map_rho,  l_map_m,    &
           l_map_id,   l_map_f,    &
           l_map_vgt,  l_map_evgt, &
           l_map_eval, l_map_aeval,&
           l_map_evec, l_map_aevec,&
           l_map_ct,   l_map_act,  &
           l_map_u,    l_map_au,   &
#if __TRACER 
           l_map_c_tracers,        & 
           l_map_dn,               &
#endif
           stat_info)        
        !----------------------------------------------------
        !  Program     :   particles_map_ghost_put
        !----------------------------------------------------
        !
        !  Purpose     :   Swap the send and receive buffers,
        !                  in order to send the contribution 
        !                  back to their host processes.
        !                 
        !
        !  Routines    :
        !
        !  References  :  Sbalzarini et al.
        !                 2006 Journal of Computational Physics.
        !
        !  Remarks     : Positions are always used for ghost
        !                layers.
        !
        !  Revisions   :  V0.2 30.07.2009, 
        !                 add vgt  evgt, eval, aeval, 
        !                 evec, aevec,ct, act,
        !                 which are usefull for  viscoelastic 
        !                 Oldroyd-B model.
        ! 
        !                 V0.1 16.06.2009, original.
        !
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
        ! Modules
      	!----------------------------------------------------
        
        ! USE ppm_module_map

        USE ppm_module_map_part_ghost
        USE ppm_module_map_part 
        USE ppm_module_map_part_util
        
      	!----------------------------------------------------
        ! Arguments
      	!----------------------------------------------------
        
        TYPE(Particles),INTENT(INOUT)   :: this
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_x
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_v
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_rho
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_m
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_id
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_f
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_vgt
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_evgt
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_eval
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_aeval
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_evec
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_aevec
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_ct
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_act
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_u
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_au
#if __TRACER 
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_c_tracers 
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_dn
#endif
        INTEGER,INTENT(OUT)             :: stat_info
        
      	!----------------------------------------------------
      	! Local Variables
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: symmetry
        INTEGER                         :: isymm

        INTEGER                         :: num_dim
        REAL(MK)                        :: ghost_size
        
      	!----------------------------------------------------
        ! Logical variables, 
        ! indicate which is sent or received.
      	!----------------------------------------------------
        
        INTEGER                         ::  map_type
        LOGICAL                         ::  map_x
        LOGICAL                         ::  map_v
        LOGICAL                         ::  map_rho
        LOGICAL                         ::  map_m
        LOGICAL                         ::  map_id
        LOGICAL                         ::  map_f
        LOGICAL                         ::  map_vgt
        LOGICAL                         ::  map_evgt
        LOGICAL                         ::  map_eval
        LOGICAL                         ::  map_aeval
        LOGICAL                         ::  map_evec
        LOGICAL                         ::  map_aevec  
        LOGICAL                         ::  map_ct
        LOGICAL                         ::  map_act
        LOGICAL                         ::  map_u
        LOGICAL                         ::  map_au
#if __TRACER 
        LOGICAL                         :: map_c_tracers 
        LOGICAL                         :: map_dn
#endif
        INTEGER                         :: topo_id

#ifdef __DEBUG
        INTEGER                        :: debug_flag
        INTEGER                        :: rank
        REAL(MK)                       :: time_routine_start
#endif 
        
      	!--------------------------------
      	! Initialize variables
      	!--------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
#ifdef __DEBUG
        
        rank = technique_get_rank(this%tech,stat_info_sub)
        debug_flag = debug_get_flag(global_debug,stat_info_sub)
        IF (debug_flag == 2 ) THEN
           CALL debug_substart(global_debug,rank,&
                "particles_map_ghost_put",&
                time_routine_start,stat_info)
        END IF
#endif
        
        map_x     = .FALSE.
        map_v     = .FALSE.
        map_rho   = .FALSE.
        map_m     = .FALSE.
        map_id    = .FALSE.
        map_f     = .FALSE.
        map_vgt   = .FALSE.
        map_evgt  = .FALSE.
        map_eval  = .FALSE.
        map_aeval = .FALSE.
        map_evec  = .FALSE.
        map_aevec = .FALSE.
        map_ct    = .FALSE.
        map_act   = .FALSE.
        map_u     = .FALSE.
        map_au    = .FALSE.
#if __TRACER 
        map_c_tracers = .FALSE.
        map_dn        = .FALSE.
#endif
        
        IF( PRESENT(l_map_x) ) THEN
           map_x = l_map_x
        END IF
        IF( PRESENT(l_map_v) ) THEN
           map_v = l_map_v
        END IF
        IF( PRESENT(l_map_rho) ) THEN
           map_rho = l_map_rho
        END IF
        IF( PRESENT(l_map_m) ) THEN
           map_m = l_map_m
        END IF
        IF( PRESENT(l_map_id) ) THEN
           map_id = l_map_id
        END IF
        IF( PRESENT(l_map_f) ) THEN
           map_f = l_map_f
        END IF
        IF( PRESENT(l_map_vgt) ) THEN
           map_vgt = l_map_vgt
        END IF
        IF( PRESENT(l_map_evgt) ) THEN
           map_evgt = l_map_evgt           
        END IF
        IF( PRESENT(l_map_eval) ) THEN
           map_eval = l_map_eval           
        END IF
        IF( PRESENT(l_map_aeval) ) THEN
           map_aeval = l_map_aeval           
        END IF
        IF( PRESENT(l_map_evec) ) THEN
           map_evec = l_map_evec           
        END IF
        IF( PRESENT(l_map_aevec) ) THEN
           map_aevec = l_map_aevec           
        END IF
        IF( PRESENT(l_map_ct) ) THEN
           map_ct = l_map_ct
        END IF
        IF( PRESENT(l_map_act) ) THEN
           map_act = l_map_act
        END IF
        IF( PRESENT(l_map_u) ) THEN
           map_u = l_map_u
        END IF
        IF( PRESENT(l_map_au) ) THEN
           map_au = l_map_au
        END IF
#if __TRACER
        IF( PRESENT(l_map_c_tracers) ) THEN
           map_c_tracers = l_map_c_tracers
        END IF  
        IF( PRESENT(l_map_dn) ) THEN
           map_dn = l_map_dn
        END IF
#endif  

        !----------------------------------------------------
        ! get topology ID
        !----------------------------------------------------

        topo_id = &
             technique_get_topo_id(this%tech,stat_info_sub) 

        !----------------------------------------------------
        ! Check if use symmetry inter-process communication.
        !----------------------------------------------------

        symmetry = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        
        IF (symmetry) THEN 
           isymm = 1
        ELSE
           isymm = 0
        END IF
        
        num_dim    = physics_get_num_dim(this%phys,stat_info_sub)
        
        ghost_size = technique_get_ghost_size(this%tech,stat_info_sub)
        
        !--------------------------------
      	! Mapping ghosts
      	!--------------------------------
        
        IF ( ghost_size <= 0.0_MK) THEN
           PRINT * ,'particles_map_ghost_put : ',&
                'Ghost size should be positive !'
           stat_info = -1
           GOTO 9999
        END IF
        
        !map_type = ppm_param_map_ghost_put   
        
        IF ( map_x ) THEN
           
           !CALL ppm_map_part_ghost(this%x,num_dim,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)  

           CALL ppm_map_part_ghost_put(topo_id, stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *,'particles_map_ghost_put',&
                'Failed to put ghosts!'           
              stat_info = -1
              GOTO 9999           
           END IF
        
        END IF
        
        !map_type = ppm_param_map_push

IF ( map_v ) THEN
           
           !CALL ppm_map_part_ghost(this%v,num_dim,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)           

           CALL ppm_map_part_push(this%v, num_dim, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        
        IF ( map_rho ) THEN
           
           !CALL  ppm_map_part_ghost(this%rho,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)   

           CALL ppm_map_part_push(this%rho, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_m ) THEN
           
           !CALL  ppm_map_part_ghost(this%m,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)       

           CALL ppm_map_part_push(this%m,&
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_id ) THEN
           
           !CALL  ppm_map_part_ghost(this%id,this%num_id,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)        

           CALL ppm_map_part_push(this%id, this%num_id, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_f ) THEN
           
           !CALL  ppm_map_part_ghost(this%f,num_dim,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)  

           CALL ppm_map_part_push(this%f, num_dim, &
                this%num_part_real, stat_info_sub)  

           ! also background pressure force
           CALL ppm_map_part_push(this%f_bp, num_dim, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_vgt ) THEN
           
           !CALL  ppm_map_part_ghost(this%vgt,num_dim**2,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)   

           CALL ppm_map_part_push(this%vgt, num_dim**2, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_evgt ) THEN
           
           !CALL  ppm_map_part_ghost(this%evgt,num_dim**2,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)   

           CALL ppm_map_part_push(this%evgt, num_dim**2, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_eval ) THEN
           
           !CALL  ppm_map_part_ghost(this%eval,num_dim,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)   

           CALL ppm_map_part_push(this%eval, num_dim, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_aeval ) THEN
           
           !CALL  ppm_map_part_ghost(this%aeval,num_dim,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)        

           CALL ppm_map_part_push(this%aeval, num_dim, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_evec ) THEN
           
           !CALL  ppm_map_part_ghost(this%evec,num_dim**2,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)      

           CALL ppm_map_part_push(this%evec, num_dim**2, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_aevec ) THEN
           
           !CALL  ppm_map_part_ghost(this%aevec,num_dim**2,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)

           CALL ppm_map_part_push(this%aevec, num_dim**2, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_ct ) THEN
           
           !CALL  ppm_map_part_ghost(this%ct,num_dim**2,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)     

           CALL ppm_map_part_push(this%ct, num_dim**2, &
                this%num_part_real, stat_info_sub)
           
        END IF

        IF ( map_act ) THEN
           
           !CALL  ppm_map_part_ghost(this%act,num_dim**2,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)     

           CALL ppm_map_part_push(this%act, num_dim**2, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_u ) THEN
           
           !CALL  ppm_map_part_ghost(this%u,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub)    

           CALL ppm_map_part_push(this%u, &
                this%num_part_real, stat_info_sub)
           
        END IF
        
        IF ( map_au ) THEN
           
           !CALL  ppm_map_part_ghost(this%au,&
           !     this%num_part_real,this%num_part_all,&
           !     isymm,ghost_size,map_type,stat_info_sub) 

           CALL ppm_map_part_push(this%au, &
                this%num_part_real, stat_info_sub)
           
        END IF

#if __TRACER

        IF ( map_c_tracers) THEN
           
           !CALL ppm_map_part_ghost(this%c_tracers,&
           !     this%num_part_real,this%num_part_all,& 
           !     isymm,ghost_size,map_type,stat_info_sub)      

           CALL ppm_map_part_push(this%c_tracers, &
                this%num_part_real, stat_info_sub)
                
        END IF 

        IF ( map_dn) THEN  

           !CALL ppm_map_part_ghost(this%dn,&
           !     this%num_part_real,this%num_part_all,& 
           !     isymm,ghost_size,map_type,stat_info_sub)    

           CALL ppm_map_part_push(this%dn, &
                this%num_part_real, stat_info_sub)
           
        END IF

#endif
        
        IF (stat_info_sub /= 0) THEN
           PRINT *,'particles_map_ghost_put : ',&
                'Failed to push buffers !'
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !map_type = ppm_param_map_send
        
        !CALL ppm_map_part_ghost(this%x,num_dim,&
        !     this%num_part_real,this%num_part_ghost,&
        !     isymm,ghost_size,map_type,stat_info_sub) 

        CALL ppm_map_part_send(this%num_part_real, &
             this%num_part_all, stat_info_sub)
          
        
        IF (stat_info_sub /= 0) THEN
           PRINT *, 'particles_map_ghost_put : ',&
                'Failed to send ghosts !'
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !map_type = ppm_param_map_pop


#if __TRACER  

        IF ( map_dn) THEN

           CALL ppm_map_part_ghost_pop(this%dn, 1, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF

        IF ( map_c_tracers ) THEN  


           CALL ppm_map_part_ghost_pop(this%c_tracers, 1, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF

#endif
        
        
        IF ( map_au ) THEN
           
           CALL ppm_map_part_ghost_pop(this%au, 1, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_u ) THEN
           
           CALL ppm_map_part_ghost_pop(this%u, 1, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_act ) THEN
           
           CALL ppm_map_part_ghost_pop(this%act, num_dim**2, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_ct ) THEN
           
           CALL ppm_map_part_ghost_pop(this%ct, num_dim**2, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_aevec ) THEN
           
           CALL  ppm_map_part_ghost_pop(this%aevec, num_dim**2, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_evec ) THEN
           
           CALL  ppm_map_part_ghost_pop(this%evec, num_dim**2, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_aeval ) THEN
   
           CALL  ppm_map_part_ghost_pop(this%aeval, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_eval ) THEN
           
           CALL  ppm_map_part_ghost_pop(this%eval, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_evgt ) THEN
           
           CALL  ppm_map_part_ghost_pop(this%evgt, num_dim**2, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_vgt ) THEN
           
           CALL ppm_map_part_ghost_pop(this%vgt, num_dim**2, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_f ) THEN 

           ! background pressure force
           CALL ppm_map_part_ghost_pop(this%f_bp, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
           CALL ppm_map_part_ghost_pop(this%f, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_id ) THEN
           
           CALL ppm_map_part_ghost_pop(this%id, this%num_id, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_m ) THEN
           
           CALL ppm_map_part_ghost_pop(this%m, 1, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        IF ( map_rho ) THEN
           
           CALL ppm_map_part_ghost_pop(this%rho, 1, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)

        END IF
        
        IF ( map_v ) THEN
           
           CALL ppm_map_part_ghost_pop(this%v, num_dim, &
                this%num_part_real, &
                this%num_part_all, stat_info_sub)
           
        END IF
        
        !CALL ppm_map_part_ghost_pop(this%x, num_dim, &
        !     this%num_part_real, &
        !     this%num_part_all, stat_info_sub)
        !print*, stat_info_sub, this%num_part_real, this%num_part_all
        
        IF (stat_info_sub /= 0) THEN          
           PRINT *,'particles_map_ghost_put : ',&
                'Failed to pop buffers !'
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! load the stored particle mapping
        !----------------------------------------------------
        
        !CALL ppm_map_part_load(stat_info_sub)

        
        
9999    CONTINUE
        
        
#ifdef __DEBUG
        IF ( debug_flag == 2 ) THEN
           CALL debug_substop(global_debug,rank,&
                "particles_map_ghost_put",&
                time_routine_start,stat_info_sub)
        END IF
        
#endif
        
        RETURN                  
        
      END SUBROUTINE particles_map_ghost_put
      
      
