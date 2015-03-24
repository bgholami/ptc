 SUBROUTINE particles_create_delete_particles(this, to_be_created, to_be_deleted, stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_create_delete_particles
        !----------------------------------------------------
        !
        ! Purpose     : creates/deletes particles for inflow/outflow BC
        !                  
        !
        ! Routines    :  
        !
        ! References  :  
        !
        !
        ! Remarks     :  
        !
        ! Revisions   : 
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
      	! Modules :
      	!----------------------------------------------------
        
        
        
	!----------------------------------------------------
      	! Arguments :
      	!----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this        
        INTEGER                         :: to_be_created
        INTEGER                         :: to_be_deleted
        INTEGER, INTENT(OUT)		:: stat_info	
        
      	!----------------------------------------------------
      	! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: ip

        !----------------------------------------------------
        ! Physics parameters :
        !
        ! num_dim        : number of dimension.
        !----------------------------------------------------
        
        INTEGER                         :: num_dim

        !----------------------------------------------------
        ! Number of all tracers
     	!----------------------------------------------------
        
        INTEGER                         :: num_part_real

        !----------------------------------------------------
        ! local variables:
        !  
        !----------------------------------------------------

        REAL(MK), DIMENSION(:,:),POINTER:: x_t, v_t, id_t
        REAL(MK), DIMENSION(:)  ,POINTER:: m_t 
#if __TRACER 
        REAL(MK), DIMENSION(:)  ,POINTER:: c_tracers_t  
        INTEGER                         :: tracers_c_factor
#endif
        INTEGER                         :: net_change, num 
        REAL(MK), DIMENSION(:), POINTER :: length, length_t 
        INTEGER , DIMENSION(:), POINTER :: bcdef  
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK)                        :: io_length

  	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0    
        
        NULLIFY(x_t)    
        NULLIFY(v_t)        
        NULLIFY(m_t)    
        NULLIFY(id_t)
#if __TRACER 
        NULLIFY(c_tracers_t)
#endif  
        NULLIFY(min_phys)
        NULLIFY(length) 
        NULLIFY(length_t)
        NULLIFY(bcdef)

        
        !----------------------------------------------------
        ! Physics parameters :
        !
        ! from a object of Physics class.
        !
        !----------------------------------------------------
        
        num_dim     = &
             physics_get_num_dim(this%phys,stat_info_sub)   

        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_length(this%phys,length,stat_info_sub)
        CALL physics_get_length_t(this%phys,length_t,stat_info_sub)
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)

        io_length = (length_t(1) - length(1)) / 2.0_MK ! in x direction

#if __TRACER
        tracers_c_factor = &
             physics_get_tracers_c_factor(this%phys,stat_info_sub)
#endif


        !----------------------------------------------------
        ! Number of real particles
        !----------------------------------------------------
        num_part_real = this%num_part_real

        net_change = to_be_created - to_be_deleted


        IF ((to_be_created > 0) .OR. (to_be_deleted > 0)) THEN
           
           ALLOCATE(x_t(num_dim,num_part_real+net_change), &
                STAT=stat_info_sub)         
           ALLOCATE(v_t(num_dim,num_part_real+net_change), &
                STAT=stat_info_sub)   
           ALLOCATE(m_t(num_part_real+net_change), &
                STAT=stat_info_sub) 
           ALLOCATE(id_t(this%num_id,num_part_real+net_change), &
                STAT=stat_info_sub)  
#if __TRACER  
           ALLOCATE(c_tracers_t(num_part_real+net_change), &
                STAT=stat_info_sub)
#endif
           
           
           num = 0
           DO ip = 1, num_part_real
              
              IF (this%id(this%pid_idx, ip) < 0) THEN ! delete

                 CONTINUE

              ELSE

                 num = num + 1
                
                 ! rest of particles: copy   
                 x_t(1:num_dim,num) = this%x(1:num_dim,ip)
                 v_t(1:num_dim,num) = this%v(1:num_dim,ip)
                 m_t(num)           = this%m(ip)
                 id_t(1:this%num_id,num) = this%id(1:this%num_id,ip)
#if __TRACER 
                 c_tracers_t(num)   = this%c_tracers(ip) 
#endif

                 IF ((bcdef(1) == mcf_bcdef_inflow) .AND. &
                      (this%id(this%bid_idx, ip) == 1) .AND. &
                      (this%x(1, ip) > min_phys(1))) THEN
                    
                    ! an inflow particle that just entered interior: duplicate (inflow)

                    ! first modify id and c_tracers of the particles itself
                    id_t(this%bid_idx, num) = 0
#if __TRACER 
                    c_tracers_t(num)   = DBLE(tracers_c_factor)
#endif

                    ! now create a new one
                    num = num + 1

                    x_t(1:num_dim,num) = this%x(1:num_dim,ip)
                    x_t(1,num)         = x_t(1,num) - io_length ! assumes flow in x direction
                    v_t(1:num_dim,num) = this%v(1:num_dim,ip) ! velocity will be set to inlet velocity
                    m_t(num)           = this%m(ip)
                    id_t(this%pid_idx,num) = 1  ! could be made unique
                    id_t(this%sid_idx,num) = 0  ! fluid particle
                    id_t(this%bid_idx,num) = 1  ! new inflow particle
#if __TRACER 
                    c_tracers_t(num)   = DBLE(tracers_c_factor) ! just to not mess with the diffusion
#endif

                 END IF

              END IF
                 
              
           END DO
           
           num_part_real = num_part_real + net_change
           this%num_part_real = num_part_real 
           this%num_part_all = this%num_part_all + net_change
           
           IF (ASSOCIATED(this%x)) THEN
              DEALLOCATE(this%x,STAT=stat_info_sub)
           END IF
           IF (ASSOCIATED(this%v)) THEN
              DEALLOCATE(this%v,STAT=stat_info_sub)
           END IF
           IF (ASSOCIATED(this%m)) THEN
              DEALLOCATE(this%m,STAT=stat_info_sub)
           END IF
           IF (ASSOCIATED(this%id)) THEN
              DEALLOCATE(this%id,STAT=stat_info_sub)
           END IF 
#if __TRACER 
           IF (ASSOCIATED(this%c_tracers)) THEN
              DEALLOCATE(this%c_tracers,STAT=stat_info_sub)
           END IF
#endif
           
           ALLOCATE(this%x(num_dim,num_part_real), &
                STAT=stat_info_sub)         
           ALLOCATE(this%v(num_dim,num_part_real), &
                STAT=stat_info_sub)   
           ALLOCATE(this%m(num_part_real), &
                STAT=stat_info_sub) 
           ALLOCATE(this%id(this%num_id,num_part_real), &
                STAT=stat_info_sub) 
#if __TRACER 
           ALLOCATE(this%c_tracers(num_part_real), &
                STAT=stat_info_sub)
#endif
           
           
           
           this%x = x_t        
           this%v = v_t   
           this%m = m_t    
           this%id = id_t
#if __TRACER 
           this%c_tracers = c_tracers_t
#endif

        END IF
        
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE  
        
        
        
        IF (ASSOCIATED(x_t)) THEN
           DEALLOCATE(x_t,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(v_t)) THEN
           DEALLOCATE(v_t,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(m_t)) THEN
           DEALLOCATE(m_t,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(id_t)) THEN
           DEALLOCATE(id_t,STAT=stat_info_sub)
        END IF  
#if __TRACER 
        IF (ASSOCIATED(c_tracers_t)) THEN
           DEALLOCATE(c_tracers_t,STAT=stat_info_sub)
        END IF
#endif 

        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys) 
        END IF
                
        IF(ASSOCIATED(length)) THEN
           DEALLOCATE(length) 
        END IF
        
        IF(ASSOCIATED(length_t)) THEN
           DEALLOCATE(length_t) 
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef) 
        END IF
        
        
        
        RETURN
        
      END SUBROUTINE particles_create_delete_particles
      
