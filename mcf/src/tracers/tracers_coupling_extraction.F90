SUBROUTINE tracers_coupling_extraction(this,d_particles,stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_coupling_extraction
  !----------------------------------------------------
  !
  ! Purpose     : Coupling of tracers and particles
  !               information are exchanged at 
  !               the interface (extraction)
  !                  
  !
  ! Routines    :  
  !
  ! References  :  
  !
  !
  ! Remarks     :  
  !
  !        
  ! Revisions   : V0.3 11.07.2012, 3D and arbitrary geometry 
  ! 
  !               V0.2 11.07.2012, revise for reuse 
  ! 
  !               V0.1 02.08.2011, original
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
  TYPE(Particles), INTENT(INOUT)  :: d_particles  
  INTEGER, INTENT(OUT)	          :: stat_info	

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !----------------------------------------------------

  INTEGER                         :: num_dim
  REAL(MK), DIMENSION(:), POINTER :: dx 
  REAL(MK)                        :: dt
  !----------------------------------------------------
  ! Colloid parameters :
  !
  !----------------------------------------------------
  
  TYPE(Colloid),POINTER               :: colloids
  INTEGER                             :: num_colloid 
  REAL(MK), DIMENSION(:,:,:), POINTER :: coll_v 
  INTEGER , DIMENSION(:,:,:), POINTER :: coll_f_vlist  
  REAL(MK), DIMENSION(:,:,:), POINTER :: coll_n

  !----------------------------------------------------
  ! Technique parameters :
  !---------------------------------------------------
  
  REAL(MK), DIMENSION(:), POINTER         :: min_sub
  REAL(MK), DIMENSION(:), POINTER         :: max_sub
  
  !----------------------------------------------------
  ! Number of particles and tracers
  !----------------------------------------------------

  INTEGER                         :: num_part_real_particles
  INTEGER                         :: num_part_real_tracers

  !----------------------------------------------------
  ! local variables:
  !  
  !----------------------------------------------------   

  REAL(MK), DIMENSION(:,:), ALLOCATABLE :: box_dim 
  REAL(MK)                              :: a_p, dx_diff
  REAL(MK)                              :: d_near_wall
  INTEGER                               :: t_num_gen, vperf
  INTEGER                               :: ip, i, j, vi
  INTEGER                               :: stat_info_sub 

  REAL(MK), DIMENSION(:), ALLOCATABLE   :: dumv
  INTEGER                         :: rank
  INTEGER                         :: COMM

!CALL ppm_time(t0,stat_info_sub)
  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0   
 
  NULLIFY(dx)
  NULLIFY(colloids) 
  NULLIFY(coll_v)
  NULLIFY(coll_f_vlist)  
  NULLIFY(coll_n) 
  NULLIFY(min_sub)
  NULLIFY(max_sub)

  !----------------------------------------------------
  ! Tracers' parameters :
  !----------------------------------------------------

  ! tracers' radius
  a_p = this%radius

  ! width of the near_wall region
  d_near_wall = this%near_wall


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub)  
  CALL physics_get_dx(this%phys,dx,stat_info_sub)
  dt = physics_get_dt(this%phys, stat_info_sub)

  dx_diff = (2.0_MK * dt * this%dp) ** 0.50_MK

  !----------------------------------------------------
  ! Get colloid
  !----------------------------------------------------
  
  num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
  
  IF ( num_colloid > 0 ) THEN
     
     CALL physics_get_colloid(this%phys, &
          colloids,stat_info_sub)
     
  END IF

  CALL colloid_get_coll_v(colloids, coll_v, stat_info_sub)  
  CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub) 
  CALL colloid_get_coll_n(colloids, coll_n, stat_info_sub)

  !----------------------------------------------------
  ! Get the boundary of this sub-domain. 
  !----------------------------------------------------
  
  CALL technique_get_min_sub(this%tech,min_sub,stat_info_sub)
  CALL technique_get_max_sub(this%tech,max_sub,stat_info_sub)


  !----------------------------------------------------
  ! Number of real particles and tracers
  !----------------------------------------------------
  num_part_real_particles = &
       particles_get_num_part_real(d_particles,stat_info_sub)    

  num_part_real_tracers   = &
       tracers_get_num_part_real(this,stat_info_sub)


  !-------------------------------------------------
  ! Allocate memory and initialize
  !-------------------------------------------------

  vperf = num_dim  

  ALLOCATE(box_dim(num_dim, 2*vperf))

  ALLOCATE(dumv(num_dim))
!CALL ppm_time(t1,stat_info_sub)
  !-------------------------------------------------
  ! Calculate wall parameters
  ! for example: h_p (wall distance)
  !              n_vec (wall normal unit vector)
  !              facet (index of facet a particle belongs to)
  !              col_i (colloid index particle belongs to)
  !-------------------------------------------------

  CALL tracers_compute_wall_distance(d_particles, &
       0, stat_info_sub)
  
!CALL ppm_time(t2,stat_info_sub)
  !---------------------------------------------------
  ! Loop over all particles
  ! if it is close enough, determine how many 
  !    tracers must be generated and...
  ! else, set num_gen_tracers to zero
  !---------------------------------------------------

  DO ip = 1, num_part_real_particles   
     
     IF ((d_particles%id(2, ip) == 0) .AND. (d_particles%h_p(ip) <= d_near_wall)) THEN


        ! number of tracers to be generated
        t_num_gen = INT(d_particles%c_tracers(ip))


        ! We continue only if a positive number of 
        ! tracers are to be generated.
        IF (t_num_gen > 0) THEN

           !---------------------------------------------------
           ! Determine the dimensions of the box in 
           ! which new tracers should be created.
           !---------------------------------------------------
           
           ! box is a quadrilateral in 2D  
           ! The order is going from 1 to vperf index in coll_v (base
           ! located on colloid wall) and then 1 to vperf index in
           ! interface (the other base).
           ! (the order and tracer_coupling_generate must be consistent.)
           dumv = 0.0_MK
           DO i = 1, vperf
              vi = coll_f_vlist(d_particles%col_i(ip), i, d_particles%facet(ip))
              box_dim(1:num_dim, i) = coll_v(d_particles%col_i(ip), 1:num_dim, vi) + &
                   2.0_MK * a_p * coll_n(d_particles%col_i(ip), 1:num_dim, d_particles%facet(ip))
              box_dim(1:num_dim, i+vperf) = &
                   this%x_interface(d_particles%col_i(ip), 1:num_dim, vi) 


              dumv = dumv + this%vn_interface(d_particles%col_i(ip), 1:num_dim, vi)
!!$              DO j = 1, num_dim
!!$                 box_dim(j, i) = MAX(box_dim(j, i), min_sub(j) + a_p)
!!$                 box_dim(j, i) = MIN(box_dim(j, i), max_sub(j) - a_p)
!!$
!!$                 box_dim(j, i+vperf) = MAX(box_dim(j, i+vperf), min_sub(j) + a_p)
!!$                 box_dim(j, i+vperf) = MIN(box_dim(j, i+vperf), max_sub(j) - a_p)
!!$              END DO

           END DO
           dumv      = dumv      / DBLE(vperf) 
           !dumv      = dumv      / 2.0_MK ! middle of box_dim 
           !dumv      = dumv * h_p(ip) / d_near_wall  ! refer to tracers_near_wall_vel.

           !---------------------------------------------------
           ! Generate new tracers.
           ! position, velocity and IDs will be calculated. 
           !---------------------------------------------------
           !box_dim(1:num_dim, 1) = d_particles%x(1:num_dim, ip) ! to extract exactly at particle's position
           !box_dim(1:num_dim, 1) = d_particles%x(1:num_dim, ip) + &
           !     n_vec(1:num_dim, ip) * (d_near_wall - dx_diff/2.0_MK - h_p(ip))
           CALL tracers_generate_new(this, &
                t_num_gen, &
                box_dim(1:num_dim, 1:2*vperf), &
                !d_particles%v(1:num_dim, ip), &
                dumv(1:num_dim), &
                stat_info_sub)  

           IF (stat_info_sub /= 0 ) THEN
              PRINT *, "tracers_particle_coupling : ",&
                   "Generating new tracers failed ! "
              stat_info = -1
              GOTO 9999
           END IF

           ! Update the number of tracers carried by this 
           ! SPH particle. 
           d_particles%c_tracers(ip) = &
                d_particles%c_tracers(ip) - t_num_gen 
           
        END IF

     END IF
  
  END DO
!CALL ppm_time(t3,stat_info_sub)
  !----------------------------------------------------
  ! update the tracer count (new tracers are added)
  !----------------------------------------------------
  !num_part_real_tracers   = &
  !     tracers_get_num_part_real(this,stat_info_sub)  
 
! rank     = technique_get_rank(this%tech,stat_info_sub)
! COMM     = technique_get_comm(this%tech,stat_info_sub)
!call MPI_Barrier(COMM, stat_info_sub)
!if (rank == 0) then
!print*, 'crank:', rank, t1-t0, t2-t0, t3-t0
!end if

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE  


  IF(ASSOCIATED(dx)) THEN
     DEALLOCATE(dx)
  END IF
  
   IF(ASSOCIATED(coll_v))THEN
      DEALLOCATE(coll_v)
   END IF

  IF(ASSOCIATED(coll_f_vlist))THEN
     DEALLOCATE(coll_f_vlist)
  END IF 

  IF(ASSOCIATED(coll_n))THEN
     DEALLOCATE(coll_n)
  END IF

  IF(ASSOCIATED(min_sub)) THEN
     DEALLOCATE(min_sub)
  END IF
  
  IF(ASSOCIATED(max_sub)) THEN
     DEALLOCATE(max_sub)
  END IF
  

  RETURN

END SUBROUTINE tracers_coupling_extraction
