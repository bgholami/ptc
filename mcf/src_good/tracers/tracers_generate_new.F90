SUBROUTINE tracers_generate_new(this, t_num_gen, box_dim, p_velocity, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_coupling_generate
  !----------------------------------------------------
  !
  ! Purpose     : generates new tracers in the overlapping 
  !               zone in the near-wall region. Initializes 
  !               position, velocity, and ID of each new tracer.
  !
  !                  
  !
  ! Routines    :  
  !
  ! References  :  
  !
  !
  ! Remarks     :  
  !
  ! Revisions   : V0.2 05.02 2013, 3D and arbitrary geometry
  ! 
  !               V0.1 03.08 2011, original
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

  TYPE(Particles)          , INTENT(INOUT) :: this        
  INTEGER                  , INTENT(IN)    :: t_num_gen
  REAL(MK), DIMENSION(:,:) , INTENT(IN)    :: box_dim
  REAL(MK), DIMENSION(:)   , INTENT(IN)    :: p_velocity
  INTEGER                  , INTENT(OUT)   :: stat_info 


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !----------------------------------------------------

  INTEGER                         :: num_dim	


  !----------------------------------------------------
  ! Number of tracers
  !----------------------------------------------------

  INTEGER                         :: num_part_real  

  !----------------------------------------------------
  ! Colloid parameters :
  !
  !----------------------------------------------------

  TYPE(Colloid),POINTER               :: colloids
  INTEGER                             :: num_colloid 

  !----------------------------------------------------
  ! local variables:
  !  
  !----------------------------------------------------

  TYPE(Technique)                     :: d_tech
  REAL(MK), DIMENSION(:,:),POINTER    :: x_t, v_t, id_t  
  REAL(MK)                            :: uni_rnd
  INTEGER                             :: id_factor, local_id
  INTEGER                             :: large_constant_id
  INTEGER                             :: ip, i, it 
  REAL(MK)                            :: distance
  INTEGER             	              :: col_index, fi 
  REAL(MK), DIMENSION(:), ALLOCATABLE :: v_normal
  INTEGER                             :: stat_info_sub

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0    

  NULLIFY(x_t)    
  NULLIFY(v_t)        
  NULLIFY(id_t)   

  NULLIFY(colloids)


  !----------------------------------------------------
  ! object of class Technique
  !----------------------------------------------------
  d_tech = this%tech  
  
  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------
  
  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub) 
  

  !----------------------------------------------------
  ! Number of real tracers
  !----------------------------------------------------
  num_part_real = &
       tracers_get_num_part_real(this,stat_info_sub)    

  !----------------------------------------------------
  ! Get colloid
  !----------------------------------------------------
  
  num_colloid = &
       physics_get_num_colloid(this%phys,stat_info_sub)
  
  IF ( num_colloid > 0 ) THEN
     
     CALL physics_get_colloid(this%phys, &
          colloids,stat_info_sub)
     
  END IF


  !----------------------------------------------------
  ! Allocate memory for the temporary arrays:
  ! x_t, v_t, id_t, v_normal
  !----------------------------------------------------
  
  ALLOCATE(x_t(num_dim,num_part_real+t_num_gen), &
       STAT=stat_info_sub)         
  ALLOCATE(v_t(num_dim,num_part_real+t_num_gen), &
       STAT=stat_info_sub)   
  ALLOCATE(id_t(this%num_id,num_part_real+t_num_gen), &
       STAT=stat_info_sub)

  ALLOCATE(v_normal(num_dim))

  !----------------------------------------------------
  ! variables for producing unique tracer ID
  !----------------------------------------------------

  ! order of number of processes
  id_factor = FLOOR(LOG10(REAL(d_tech%num_proc))) + 1

  ! initialize the local_id to the last local 
  ! id (largest) this process has given
  local_id = this%last_local_id 

  !----------------------------------------------------
  ! for large constant id
  !---------------------------------------------------- 
  large_constant_id = SIZE(this%num_tracer_deposited, 2) + 1

  !----------------------------------------------------
  ! Create new tracers:
  ! initialize position, velocity, and ID
  !----------------------------------------------------
  DO ip = num_part_real+1, num_part_real+t_num_gen

     ! initialize position according to a uniform 
     ! distribution within the box.
     uni_rnd = random_random_uniform(this%random, stat_info_sub)

     x_t(1:num_dim, ip) = tool_arbitrary_point_picking(this%tool, &
          uni_rnd, box_dim(1:num_dim, 1:2*num_dim), stat_info_sub)
     !x_t(1:num_dim, ip) = box_dim(1:num_dim, 1) ! exactly at particle's position

     ! initialize velocity:
     !v_t(1:num_dim, ip) = p_velocity(1:num_dim) 
     col_index = 0 ! this means all colloids will be searched
     CALL colloid_arbitrary_distance(colloids, &
          .TRUE., & ! use sort data
          col_index, & 
          x_t(1:num_dim, ip), &  ! input position
          distance, &  ! output distance
          v_normal, &  ! output normal vector
          fi,       &  ! output facet
          stat_info_sub)

     v_t(1:num_dim, ip) = p_velocity(1:num_dim) * distance / this%near_wall

     ! initialize ID
     ! Every process "signs" the tracers it generates 
     ! by writing its rank to the right most digits 
     ! of the ID. This guarantees the uniquness of tracer 
     ! IDs generated by different processes at the same time. 
     local_id = local_id + 1
     !id_t(1,ip) = d_tech%rank + local_id * 10**id_factor ! unique ids throughout the simulation
     id_t(1,ip) = large_constant_id  ! same for all, but larger than number of facets
     id_t(2:this%num_id,ip) = 0

  END DO

  !----------------------------------------------------
  ! store the value of the last local ID
  !----------------------------------------------------

  this%last_local_id = local_id


  !----------------------------------------------------
  ! Copy the data from the already existing tracers
  !----------------------------------------------------
  IF (num_part_real > 0) THEN
     x_t(1:num_dim,1:num_part_real) = &
          this%x(1:num_dim, 1:num_part_real)  
     v_t(1:num_dim,1:num_part_real) = &
          this%v(1:num_dim, 1:num_part_real)  
     id_t(1:this%num_id,1:num_part_real) = &
          this%id(1:this%num_id, 1:num_part_real) 
  END IF

  !----------------------------------------------------
  ! update tracer count
  !----------------------------------------------------
  num_part_real = num_part_real + t_num_gen
  this%num_part_real = num_part_real 
  this%num_part_all = num_part_real 


  !----------------------------------------------------
  ! Allocate memory for the permanent arrays
  !----------------------------------------------------
  IF (ASSOCIATED(this%x)) THEN
     DEALLOCATE(this%x,STAT=stat_info_sub)
  END IF
  IF (ASSOCIATED(this%v)) THEN
     DEALLOCATE(this%v,STAT=stat_info_sub)
  END IF
  IF (ASSOCIATED(this%id)) THEN
     DEALLOCATE(this%id,STAT=stat_info_sub)
  END IF

  ALLOCATE(this%x(num_dim,num_part_real), &
       STAT=stat_info_sub)         
  ALLOCATE(this%v(num_dim,num_part_real), &
       STAT=stat_info_sub)   
  ALLOCATE(this%id(this%num_id,num_part_real), &
       STAT=stat_info_sub)


  !----------------------------------------------------
  ! Copy everything to permanent arrays
  !---------------------------------------------------  
  this%x = x_t        
  this%v = v_t   
  this%id = id_t


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

  IF (ASSOCIATED(id_t)) THEN
     DEALLOCATE(id_t,STAT=stat_info_sub)
  END IF 


  RETURN


END SUBROUTINE tracers_generate_new

