SUBROUTINE particles_set_inflow_outflow_velocity(this, num, stat_info)
  !----------------------------------------------------
  ! Subroutine  : particles_set_inflow_outflow_velocity
  !----------------------------------------------------
  !
  ! Purpose     : set particle velocities at inflow and 
  !               outflow boundaries
  !                  
  !                  
  ! Reference   :
  !
  ! Remark      : 
  !
  ! Revisions   : 
  !
  !----------------------------------------------------
  ! Author      : 
  ! Contact     : 
  ! 
  ! Dr. Marco Ellero's Emmy Noether Group,
  ! Prof. Dr. N. Adams' Chair of Aerodynamics,
  ! Faculty of Mechanical Engineering,
  ! Technische Universitaet Muenchen, Germany.
  !----------------------------------------------------  


  !----------------------------------------------------
  ! Argumetns :
  !----------------------------------------------------

  TYPE(Particles),INTENT(INOUT)           :: this  
  INTEGER, INTENT(IN)                     :: num
  INTEGER,INTENT(OUT)                     :: stat_info  


  !----------------------------------------------------
  ! Physics, boundary parameters :
  !----------------------------------------------------

  INTEGER                                :: num_dim  
  REAL(MK)                               :: time
  TYPE(Boundary), POINTER                :: d_boundary  
  INTEGER                                :: num_time
  INTEGER , DIMENSION(:)  , POINTER      :: total_num_grid
  REAL(MK), DIMENSION(:)  , POINTER      :: BC_time
  INTEGER , DIMENSION(:,:)  , POINTER    :: num_grid 
  REAL(MK), DIMENSION(:,:,:), POINTER    :: tvec   
  REAL(MK), DIMENSION(:,:)  , POINTER    :: dtvec      
  REAL(MK), DIMENSION(:,:)  , POINTER    :: grido 
  REAL(MK), DIMENSION(:,:,:,:), POINTER  :: BC_vel
  REAL(MK)                               :: tstart

  !----------------------------------------------------
  !  Local variables
  !----------------------------------------------------

  INTEGER                                 :: stat_info_sub
  INTEGER                                 :: ip, j, i
  INTEGER                                 :: indx, indt, indtl, indtr
  REAL(MK)                                :: dt_pulse, t_add
  INTEGER , DIMENSION(2)                  :: ind
  REAL(MK), DIMENSION(3)                  :: r, vt, vl, vr



  !----------------------------------------------------
  !  Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0    


  NULLIFY(d_boundary)    
  NULLIFY(total_num_grid)
  NULLIFY(BC_time)
  NULLIFY(num_grid)
  NULLIFY(tvec)
  NULLIFY(dtvec)
  NULLIFY(grido)
  NULLIFY(BC_vel)


  !----------------------------------------------------
  ! Physics_parameters :
  !----------------------------------------------------

  num_dim = this%num_dim   
  time = physics_get_time_current(this%phys, stat_info_sub)  

  !----------------------------------------------------
  ! Boundary parameters :
  !----------------------------------------------------

  CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)
  num_time = boundary_get_num_time(d_boundary, stat_info_sub)
  CALL boundary_get_total_num_grid(d_boundary, total_num_grid, stat_info_sub)
  CALL boundary_get_BC_time(d_boundary, BC_time, stat_info_sub)  
  CALL boundary_get_num_grid(d_boundary, num_grid, stat_info_sub) 
  CALL boundary_get_tvec(d_boundary, tvec, stat_info_sub)
  CALL boundary_get_dtvec(d_boundary, dtvec, stat_info_sub) 
  CALL boundary_get_grido(d_boundary, grido, stat_info_sub)
  CALL boundary_get_BC_vel(d_boundary, BC_vel, stat_info_sub)
  tstart = boundary_get_tstart(d_boundary, stat_info_sub)

  dt_pulse = BC_time(2) - BC_time(1)

  !----------------------------------------------------
  ! Update velcoity of fluid particles.
  !----------------------------------------------------

  DO ip = 1, num

     IF ( (this%id(this%sid_idx, ip) == mcf_particle_type_fluid ) .AND. &
          (this%id(this%bid_idx, ip) /= 0)) THEN ! i.e. inflow/outflow boundary particle

        this%v(1:num_dim, ip) = 0.0_MK
        ! and just to be safe
        this%f(1:num_dim, ip) = 0.0_MK
        ! and this (since no periodic BC)
        this%f_bp(1:num_dim, ip) = 0.0_MK

        !----------------------------------------------------
        ! evaluate velocity (interpolated in time and space) 
        !----------------------------------------------------

        ! patch index
        j = ABS(this%id(this%bid_idx, ip))

        ! calculate space and time indices
        r(1:num_dim) = this%x(1:num_dim, ip) - grido(j, 1:num_dim)

        ! space indices
        DO i = 1, 2
           ind(i) = FLOOR(DOT_PRODUCT(r(1:num_dim), tvec(j, i, 1:num_dim)) / dtvec(j, i)) + 1
        END DO
        ! translate to index in BC_vel
        indx = (ind(1) - 1) * num_grid(j, 2) + ind(2)

        ! time index
        indt  = FLOOR((time - tstart) / dt_pulse) + 1     
        t_add = (time-tstart) - (indt - 1) * dt_pulse
        indt  = MOD(indt-1, num_time) + 1 ! take care of periodicity of time


        ! simple averaging over space
        indtl = indt
        vl(1:num_dim) = (BC_vel(j, 1:num_dim, indtl, indx) + &
             BC_vel(j, 1:num_dim, indtl, indx+1) + &
             BC_vel(j, 1:num_dim, indtl, indx+num_grid(j, 2)) + &      
             BC_vel(j, 1:num_dim, indtl, indx+num_grid(j, 2)+1)) / 4.0_MK        

        indtr = indt + 1
        indtr = MOD(indtr-1, num_time) + 1
        vr(1:num_dim) = (BC_vel(j, 1:num_dim, indtr, indx) + &
             BC_vel(j, 1:num_dim, indtr, indx+1) + &
             BC_vel(j, 1:num_dim, indtr, indx+num_grid(j, 2)) + &      
             BC_vel(j, 1:num_dim, indtr, indx+num_grid(j, 2)+1)) / 4.0_MK

        ! linear interpolation over time
        vt(1:num_dim) = vl(1:num_dim) * (1.0_MK - t_add/dt_pulse) + &
             vr(1:num_dim) * (t_add/dt_pulse)


        IF (time >= tstart) THEN

#if __TRACER
           CALL physics_set_tracers_c_factor(this%phys, 5, stat_info_sub)
#endif  
           ! velocity profile
           this%v(1:num_dim, ip) = vt(1:num_dim)

        ELSE

           this%v(1:num_dim, ip) = 0.0_MK * vt(1:num_dim) * time / tstart


        END IF

     END IF
     
  END DO ! ip =1, num

9999 CONTINUE     


  IF (ASSOCIATED(total_num_grid)) THEN   
     DEALLOCATE(total_num_grid)  
  END IF

  IF (ASSOCIATED(BC_time)) THEN
     DEALLOCATE(BC_time)  
  END IF

  IF (ASSOCIATED(num_grid)) THEN
     DEALLOCATE(num_grid)  
  END IF

  IF (ASSOCIATED(tvec)) THEN
     DEALLOCATE(tvec)  
  END IF

  IF (ASSOCIATED(dtvec)) THEN
     DEALLOCATE(dtvec)  
  END IF

  IF (ASSOCIATED(grido)) THEN
     DEALLOCATE(grido)  
  END IF

  IF (ASSOCIATED(BC_vel)) THEN
     DEALLOCATE(BC_vel)  
  END IF


  RETURN

END SUBROUTINE particles_set_inflow_outflow_velocity

