SUBROUTINE tracers_integrate(this,d_particles,time, dt,stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_integrate
  !----------------------------------------------------
  !
  ! Purpose     : Time integration of tracers
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
  ! Revisions   : V0.2 21.11.2011, two-way coupling
  !
  !               V0.1 21.07 2011, original
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
  REAL(MK)       , INTENT(IN)     :: time
  REAL(MK)       , INTENT(IN)     :: dt
  INTEGER        , INTENT(OUT)	  :: stat_info	


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !----------------------------------------------------

  INTEGER                         :: num_dim  


  !----------------------------------------------------
  ! MPI parameters.
  !----------------------------------------------------
  
  INTEGER                                 :: rank
  INTEGER                                 :: COMM
  INTEGER                                 :: MPI_PREC 


  !----------------------------------------------------
  ! Tracers' parameters
  !----------------------------------------------------

  INTEGER                             :: num_part_real    
  INTEGER, DIMENSION(:) , ALLOCATABLE :: num_interface  


  !----------------------------------------------------
  ! local variables:
  !
  ! num_steps      : number of intermediate time steps
  ! dt_inter       : dt corresponding to tracers
  !----------------------------------------------------

  INTEGER                         :: num_steps
  REAL(MK)                        :: dt_inter  
  INTEGER                         :: it, stat_info_sub
  INTEGER                         :: col_num
  REAL(MK), DIMENSION(:,:,:) , ALLOCATABLE  :: V_INT 



  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0    


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
  ! MPI parameters.
  !----------------------------------------------------
  
  rank     = technique_get_rank(this%tech,stat_info_sub)
  MPI_PREC = technique_get_MPI_PREC(this%tech,stat_info_sub)
  COMM     = technique_get_comm(this%tech,stat_info_sub)


  !----------------------------------------------------
  ! determine the number of intermediate steps and dt_inter
  !----------------------------------------------------
  num_steps = FLOOR(dt / this%tau_p) + 1
  dt_inter  = dt / DBLE(num_steps)


  !----------------------------------------------------
  ! Evaluate new velocity at the interface
  !----------------------------------------------------    
  col_num = SIZE(this%x_interface, 1) 

  ALLOCATE(num_interface(col_num))
  num_interface = this%num_interface  

  ALLOCATE(V_INT(col_num, num_dim, MAXVAL(num_interface)))
  V_INT = 0.0_MK
  this%vn_interface = 0.0_MK

  CALL tracers_interpolate_velocity(this, &
       d_particles, &
       V_INT, &
       stat_info_sub)

  ! Allreduce interface velocity to have it global
  CALL MPI_ALLREDUCE(V_INT(1:col_num, 1:num_dim, 1:MAXVAL(num_interface)), &
       this%vn_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface)), &
       col_num*num_dim*MAXVAL(num_interface), &
       MPI_PREC, &
       MPI_SUM, &
       COMM, &
       stat_info_sub)


  !----------------------------------------------------
  ! Coupling of tracers and particles at the interface
  !----------------------------------------------------

  CALL tracers_coupling_extraction(this, &
       d_particles, &
       stat_info_sub)

  !----------------------------------------------------
  ! Apply the collision term: displacements are superimposed.
  !----------------------------------------------------  


  num_part_real =  this%num_part_real  

  IF (ASSOCIATED(this%h_p)) THEN
     DEALLOCATE(this%h_p)
  END IF

  IF (ASSOCIATED(this%n_vec)) THEN
     DEALLOCATE(this%n_vec)
  END IF

  IF (ASSOCIATED(this%facet)) THEN
     DEALLOCATE(this%facet)
  END IF

  IF (ASSOCIATED(this%col_i)) THEN
     DEALLOCATE(this%col_i)
  END IF

  ALLOCATE(this%h_p(num_part_real))
  ALLOCATE(this%n_vec(num_dim, num_part_real))
  ALLOCATE(this%facet(num_part_real))  
  ALLOCATE(this%col_i(num_part_real)) 


  CALL tracers_compute_wall_distance(this, &
       1, stat_info_sub)

  CALL tracers_compute_force_collision(this, &
       dt, stat_info_sub)  

  ! tracers have moved, so...
!  CALL tracers_adjust_tracers(this,&
!       num_part_real,stat_info_sub)  
  CALL tracers_adjust_tracers_special(this,&
       num_part_real,stat_info_sub)


  if (MOD(INT(time/dt), 10) == 1) then

     CALL tracers_decompose_ring( this,&
          l_map_x    = .TRUE., l_map_v  = .TRUE., &
          l_map_id   = .TRUE., &
          stat_info  = stat_info_sub )  

  end if

  num_part_real =  this%num_part_real 

  ! wall distance again (so we do not calculate it anymore 
  ! in the integration loop, but only do a partial update,
  ! i.e. an approximation based on old and new position and
  ! under the assumption that the tracer belongs to the same
  ! facet for the entire integration loop.

  IF (ASSOCIATED(this%h_p)) THEN
     DEALLOCATE(this%h_p)
  END IF

  IF (ASSOCIATED(this%n_vec)) THEN
     DEALLOCATE(this%n_vec)
  END IF

  IF (ASSOCIATED(this%facet)) THEN
     DEALLOCATE(this%facet)
  END IF

  IF (ASSOCIATED(this%col_i)) THEN
     DEALLOCATE(this%col_i)
  END IF

  ALLOCATE(this%h_p(num_part_real))
  ALLOCATE(this%n_vec(num_dim, num_part_real))
  ALLOCATE(this%facet(num_part_real))  
  ALLOCATE(this%col_i(num_part_real))   

  
  CALL tracers_compute_wall_distance(this, &
       1, stat_info_sub)


  !----------------------------------------------------
  !----------------------------------------------------
  ! Integration loop
  !----------------------------------------------------
  !----------------------------------------------------

  DO it = 1, num_steps

     !-------------------------------------------------
     ! calculate interface velocity at the current step 
     !-------------------------------------------------
     this%vc_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface)) = &
          (this%vn_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface)) - &
          this%vo_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface))) * &
          DBLE(it) / DBLE(num_steps) + &
          this%vo_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface)) 


     !-------------------------------------------------
     ! integrate the current step 
     !-------------------------------------------------
     CALL tracers_integrate_step(this, &
          dt_inter, &
          stat_info_sub)  

     IF ( stat_info_sub /= 0 ) THEN
        PRINT *,"tracers_integrate : ", &
             "tracers_integrate_step failed !"
        stat_info = -1
        GOTO 9999
     END IF


  END DO

  !----------------------------------------------------
  ! keep the new interface velocity for the next timestep
  !----------------------------------------------------
  this%vo_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface)) = &
       this%vn_interface(1:col_num, 1:num_dim, 1:MAXVAL(num_interface))

  
  !----------------------------------------------------
  ! Coupling of tracers and particles at the interface
  !----------------------------------------------------
  CALL tracers_coupling_insertion(this, d_particles, stat_info_sub)


  !----------------------------------------------------
  ! Return.
  !---------------------------------------------------- 

9999 CONTINUE  


  RETURN

END SUBROUTINE tracers_integrate

