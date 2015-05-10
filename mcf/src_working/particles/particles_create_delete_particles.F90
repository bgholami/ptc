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

  TYPE(Particles), INTENT(INOUT)    :: this        
  INTEGER                           :: to_be_created
  INTEGER                           :: to_be_deleted
  INTEGER, INTENT(OUT)              :: stat_info	

  !----------------------------------------------------
  ! Local variables start here :
  !----------------------------------------------------

  INTEGER                           :: stat_info_sub
  INTEGER                           :: ip  
  TYPE(Boundary), POINTER           :: d_boundary  
  REAL(MK), DIMENSION(:,:), POINTER :: iopatch_n
  REAL(MK)                          :: bwidth, distance
  INTEGER                           :: patch_id

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !----------------------------------------------------

  INTEGER                           :: num_dim

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

  NULLIFY(d_boundary)  
  NULLIFY(iopatch_n)  


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub)   

  CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)  
  CALL boundary_get_iopatch_n(d_boundary, iopatch_n, stat_info_sub) 
  bwidth = boundary_get_bwidth(d_boundary, stat_info_sub)

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

           IF (this%id(this%bid_idx, ip) > 0) THEN

              ! look for inflow particles that just entered interior
              distance = 0.0_MK
              patch_id = 0
              CALL boundary_check_particle_stat(d_boundary, &
                   this%x(1:num_dim, ip), 1, distance, patch_id, stat_info_sub)

              IF (distance == 0.0_MK) THEN ! means it is located at the interior
                 ! It just has entered the interior, so should be duplicated. First
                 ! save its corresponding patch
                 patch_id = this%id(this%bid_idx, ip)

                 ! Now modify id and c_tracers of the particles itself
                 id_t(this%bid_idx, num) = 0
#if __TRACER 
                 c_tracers_t(num)   = DBLE(tracers_c_factor)
#endif

                 ! And finally create a new one
                 num = num + 1

                 x_t(1:num_dim,num) = this%x(1:num_dim,ip) - &
                      bwidth * iopatch_n(patch_id, 1:num_dim)
                 v_t(1:num_dim,num) = this%v(1:num_dim,ip) ! velocity will be set to inlet velocity
                 m_t(num)           = this%m(ip)
                 id_t(this%pid_idx,num) = 1  ! could be made unique
                 id_t(this%sid_idx,num) = 0  ! fluid particle
                 id_t(this%bid_idx,num) = patch_id  ! new inflow particle
#if __TRACER 
                 c_tracers_t(num)   = DBLE(tracers_c_factor) ! just to not mess with the diffusion
#endif

              END IF

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

9999 CONTINUE  



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

  IF(ASSOCIATED(iopatch_n)) THEN
     DEALLOCATE(iopatch_n)
  END IF



  RETURN

END SUBROUTINE particles_create_delete_particles

