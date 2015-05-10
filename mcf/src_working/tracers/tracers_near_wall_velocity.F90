SUBROUTINE tracers_near_wall_velocity(this, VI, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_near_wall_velocity
  !----------------------------------------------------
  !
  ! Purpose     : interpolates interface's velocity to 
  !               tracers positions
  !                  
  !
  ! Routines    : 
  !
  !
  ! Remarks     :  For now works only with non-symmetric
  !                inter-processor communication (i.e. it 
  !                is nearly as slow as it can get!).
  !
  ! Revisions   : V0.1 21.07 2011, original
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
  ! Arguments :
  !----------------------------------------------------

  TYPE(Particles)         , INTENT(IN)  :: this 
  REAL(MK), DIMENSION(:,:), INTENT(OUT) :: VI
  INTEGER                 , INTENT(OUT) :: stat_info


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !----------------------------------------------------

  INTEGER                                 :: num_dim  
  REAL(MK), DIMENSION(:), POINTER         :: dx  

  !----------------------------------------------------
  ! Colloid parameters :
  !----------------------------------------------------

  TYPE(Colloid),POINTER                   :: colloids
  INTEGER                                 :: num_colloid  
  INTEGER , DIMENSION(:,:,:), POINTER     :: coll_f_vlist 
!!$  REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_v
  
  !----------------------------------------------------
  ! local variables
  !----------------------------------------------------

  REAL(MK), DIMENSION(:)  , ALLOCATABLE   :: temp_VI, vecn, vect
  REAL(MK), DIMENSION(:)  , ALLOCATABLE   :: p_point
  REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: x_point, v_point 
  REAL(MK), DIMENSION(:)  , ALLOCATABLE   :: lambda
  INTEGER , DIMENSION(:)  , ALLOCATABLE   :: i_index
  
  REAL(MK)                        :: a_p, d_near_wall
  REAL(MK)                        :: vn
  REAL(MK)                        :: m, dk, a, b, c, d      
  INTEGER                         :: num_part_real, vperf, ip, i
  INTEGER             	          :: stat_info_sub

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0 

  NULLIFY(dx)
  NULLIFY(colloids)
  NULLIFY(coll_f_vlist) 

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub)   

  CALL physics_get_dx(this%phys,dx,stat_info_sub)

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
  ! initialize colloid variables
  !----------------------------------------------------  
  
  CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub)  
!!$  CALL colloid_get_coll_v(colloids, coll_v, stat_info_sub) 

  !----------------------------------------------------
  ! Tracers' parameters
  !----------------------------------------------------

  ! number of tracers
  num_part_real = &
       tracers_get_num_part_real(this,stat_info_sub)   

  ! tracers' radius
  a_p = this%radius

  ! width of near-wall region
  d_near_wall = this%near_wall
  !----------------------------------------------------
  ! initialize interpolation array and other variables
  !----------------------------------------------------
  VI(:,:) = 0.0_MK

  vperf = num_dim
  ALLOCATE(i_index(vperf)) 
  ALLOCATE(x_point(num_dim, vperf), v_point(num_dim, vperf))
  ALLOCATE(p_point(num_dim))
  ALLOCATE(lambda(vperf)) 

  ALLOCATE(temp_VI(num_dim))  
  ALLOCATE(vecn(num_dim))
  ALLOCATE(vect(num_dim))

  !----------------------------------------------------
  ! interpolation loop
  !----------------------------------------------------       

  DO ip = 1, num_part_real

     ! since tracers with a non-zero id(2, ip) do not 
     ! have a valid wall distance, the near-wall 
     ! velocity is calculated only for those with 
     ! id(2, ip) equal to zero. (same as in 
     ! tracers_compute_wall_distance.F90)
     IF (this%id(2, ip) == 0) THEN

        ! find and copy corresponding interface points
        i_index(1:vperf) = coll_f_vlist(this%col_i(ip), 1:vperf, this%facet(ip))

        DO i = 1, vperf
           x_point(1:num_dim, i) = this%x_interface(this%col_i(ip), 1:num_dim, i_index(i))
           v_point(1:num_dim, i) = this%vc_interface(this%col_i(ip), 1:num_dim, i_index(i))
        END DO


        ! Method1: accurate but sensative to concave surfaces (and not curves)
!!$        ! find projection point (in normal direction) on interface
!!$        IF (num_dim == 2) THEN
!!$           ! in 2D it is intersection of two lines
!!$
!!$           ! distance from tracer to its projection
!!$           dk = ((this%x(1, ip) - x_point(1, 1))*(x_point(2, 2) - x_point(2, 1)) &
!!$                - (this%x(2, ip) - x_point(2, 1))*(x_point(1, 2) - x_point(1, 1))) / &
!!$                (this%n_vec(2, ip)*(x_point(1, 2) - x_point(1, 1)) - &
!!$                this%n_vec(1, ip)*(x_point(2, 2) - x_point(2, 1)))
!!$
!!$           
!!$        ELSE 
!!$
!!$           ! in 3D it is intersection of a line with a plane
!!$
!!$           ! find the plane (ax + by + cz + d = 0)
!!$           a =  (x_point(2, 2) - x_point(2, 1)) * (x_point(3, 3) - x_point(3, 1)) - &
!!$                (x_point(2, 3) - x_point(2, 1)) * (x_point(3, 2) - x_point(3, 1)) 
!!$           
!!$           b =  (x_point(3, 2) - x_point(3, 1)) * (x_point(1, 3) - x_point(1, 1)) - &
!!$                (x_point(3, 3) - x_point(3, 1)) * (x_point(1, 2) - x_point(1, 1))  
!!$           
!!$           c =  (x_point(1, 2) - x_point(1, 1)) * (x_point(2, 3) - x_point(2, 1)) - &
!!$                (x_point(1, 3) - x_point(1, 1)) * (x_point(2, 2) - x_point(2, 1))
!!$           
!!$           d = a * x_point(1, 1) + b * x_point(2, 1) + c * x_point(3, 1)
!!$           
!!$           ! distance from tracer to its projection
!!$           dk= -(a * this%x(1, ip) + b * this%x(2, ip) + c * this%x(3, ip) + d) / &
!!$                (a * this%n_vec(1, ip) + b * this%n_vec(2, ip) + c * this%n_vec(3, ip)) 
!!$           
!!$        END IF
!!$           
!!$        ! projection point
!!$        p_point = this%x(1:num_dim, ip) + dk * this%n_vec(1:num_dim, ip)
!!$
!!$
!!$        ! find barycentric coordinates to interpolate
!!$        CALL tool_barycentric_coordinate(this%tool, &
!!$             x_point, p_point, lambda, stat_info_sub)
!!$        
!!$        ! interpolate to projection point at interface
!!$        DO i = 1, num_dim
!!$           temp_VI(i) = DOT_PRODUCT(lambda(1:vperf), v_point(i, 1:vperf))
!!$        END DO
!!$
!!$        ! interpolate to determine the tracer's velocity ( in normal direction)
!!$        VI(1:num_dim, ip) = temp_VI(1:num_dim) * this%h_p(ip) / (this%h_p(ip) + dk)



        ! Method2: not that accurate but safe
        ! Velocity is simply averaged at interface and interpolated 
        ! in normal direction assuming all interface points are d_near_wall
        ! away from wall.
        DO i = 1, num_dim
           temp_VI(i) = SUM(v_point(i, 1:vperf))
        END DO
        temp_VI = temp_VI / DBLE(vperf)

        !VI(1:num_dim, ip) = temp_VI(1:num_dim) * (this%h_p(ip) - a_p) / (d_near_wall - a_p)
        VI(1:num_dim, ip) = temp_VI(1:num_dim) * this%h_p(ip) / d_near_wall


     END IF

  END DO

           
  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE

  IF(ASSOCIATED(dx)) THEN
     DEALLOCATE(dx)
  END IF
  
!!$  IF(ASSOCIATED(coll_v))THEN
!!$     DEALLOCATE(coll_v)
!!$  END IF

  IF(ASSOCIATED(coll_f_vlist))THEN
     DEALLOCATE(coll_f_vlist)
  END IF


  RETURN

END SUBROUTINE tracers_near_wall_velocity



