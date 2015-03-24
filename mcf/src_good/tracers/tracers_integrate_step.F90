      SUBROUTINE tracers_integrate_step(this,dt,stat_info)
        !----------------------------------------------------
        ! Subroutine  : tracers_integrate_step
        !----------------------------------------------------
        !
        ! Purpose     : Performs one step of time integration 
        !               for tracers with the corresponding dt.
        !                  
        !
        ! Routines    :  
        !
        ! References  :  
        !
        !
        ! Remarks     : Previously, different parts of this 
        !               routine were implemented in 
        !               marching_integrate_VV. 
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
      	! Modules :
      	!----------------------------------------------------
              
        
	!----------------------------------------------------
      	! Arguments :
      	!----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this 
        REAL(MK), INTENT(IN)            :: dt
        INTEGER, INTENT(OUT)		:: stat_info	
        
              
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
        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE :: x_old
        REAL(MK)                              :: dhp, vnmax, vn, vcap
        REAL(MK), DIMENSION(:)  , ALLOCATABLE :: vt
        INTEGER                               :: stat_info_sub
        INTEGER                               :: i, ip, it 

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

        if (num_part_real > 0) then

        
        !----------------------------------------------------
        ! save a copy of the current position to use for 
        ! updating wall distance 
        !----------------------------------------------------
        ALLOCATE(x_old(num_dim, num_part_real))  
        ALLOCATE(vt(num_dim))
        x_old = this%x


        !----------------------------------------------------
        ! Compute forces on tracers at this time.
        !----------------------------------------------------


        CALL tracers_compute_force(this, &  
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"tracers_integrae_step : ", &
                "Computing tracers' forces failed !"
           stat_info = -1
           GOTO 9999
        END IF 
!!$        IF (ASSOCIATED(this%f)) THEN
!!$           DEALLOCATE(this%f,STAT=stat_info_sub)
!!$        END IF
!!$        ALLOCATE(this%f(num_dim,num_part_real), STAT=stat_info_sub)
!!$        this%f = 0.0_MK


        !----------------------------------------------------
	! Position integration for tracers
	!----------------------------------------------------  

        CALL tracers_integrate_position(this,&
             num_part_real,dt,1,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *, "tracers_integrate_step : ",&
                "Integrating tracers' position failed ! "
           stat_info = -1
           GOTO 9999
        END IF


        !----------------------------------------------------
	! Calculate new velocity for tracers
        ! using lamda=1.0 as coefficient in front 
        ! of acceleration.
       	!----------------------------------------------------  

        
        CALL tracers_integrate_velocity(this,&
             num_part_real,dt,1.0_MK,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "tracers_integrae_step : ",&
                "Updating tracers' velocity failed  !"
           stat_info = -1     
           GOTO 9999
        END IF   



        !----------------------------------------------------
        ! update wall distance,
        ! i.e. hp_new = hp_old + (X_new - X_old) . normal_vec
        !----------------------------------------------------

        DO ip = 1, num_part_real

           ! evaluate normal displacement
           dhp =  DOT_PRODUCT(this%x(1:num_dim, ip) - x_old(1:num_dim, ip), &
                this%n_vec(1:num_dim, ip))

           IF ((this%h_p(ip) + dhp) >= 2.0_MK * this%radius) THEN
              
              this%h_p(ip) = this%h_p(ip) + dhp

           ELSE ! in this case, do not change h_p

              ! correct position (move parallel to wall from old position)
              this%x(1:num_dim, ip) = this%x(1:num_dim, ip) + &
                   ABS(dhp) * this%n_vec(1:num_dim, ip)

              ! correct velocity (set normal component to zero)
              !this%v(1:num_dim, ip) = this%v(1:num_dim, ip) - &
              !     DOT_PRODUCT(this%v(1:num_dim, ip), this%n_vec(1:num_dim, ip)) * this%n_vec(1:num_dim, ip)

           END IF

!!$           ! velocity cap (based on wall distance)
!!$           vnmax = (this%h_p(ip) - this%radius) / dt
!!$
!!$           vn = DOT_PRODUCT(this%v(1:num_dim, ip), this%n_vec(1:num_dim, ip))
!!$           vt = this%v(1:num_dim, ip) - vn * this%n_vec(1:num_dim, ip)
!!$
!!$           vcap = SIGN(MIN(vnmax, ABS(vn)), vn)
!!$
!!$           this%v(1:num_dim, ip) = vt + vcap * this%n_vec(1:num_dim, ip)

!!$           if ((vnmax < abs(vn)) .and. &
!!$                (((this%v(1,ip)**2.0_MK + this%v(2,ip)**2.0_MK + this%v(3,ip)**2.0_MK)**0.50_MK)/this%h_p(ip) > 1000.0_MK))then
!!$              print*, this%x(:, ip), this%v(:, ip), vt+vn*this%n_vec(1:num_dim, ip), this%h_p(ip), vn, vnmax, vcap, this%n_vec(1:num_dim, ip), vt
!!$           end if
!!$
!!$           if (this%v(1, ip) > 1000.0_MK) then
!!$              print*, this%v(:, ip), this%x(:, ip), this%h_p(ip)
!!$           end if

        END DO


        !----------------------------------------------------
        ! Adjust real tracers' r/v after motion
        ! according to boundary conditions,
	! in case they go out of the physical domain.
	!----------------------------------------------------

        CALL tracers_adjust_tracers(this,&
             num_part_real,stat_info_sub)

        IF (stat_info_sub /= 0 ) THEN
           PRINT *, "tracers_integrae_step : ",&
                "Adjusting tracers' r or v failed ! "
           stat_info = -1
           GOTO 9999
        END IF 
        
        

!!$        !----------------------------------------------------
!!$        ! Tracers:
!!$        ! Decompose partially, since positions have changed.
!!$        !
!!$        ! Only velocity and ID are commmunicated for tracers.
!!$        !--------------------------------------------------------
!!$
!!$        CALL tracers_decompose_partial( this,&
!!$             l_map_x    = .TRUE., l_map_v  = .TRUE., &
!!$             l_map_id   = .TRUE., &
!!$             stat_info  = stat_info_sub )
!!$
!!$
!!$        IF ( stat_info_sub /= 0 ) THEN
!!$           PRINT *,"tracers_integrae_step : ", &
!!$                "Decomposing tracers partially failed !"
!!$           stat_info = -1
!!$           GOTO 9999
!!$        END IF
!!$
!!$        !----------------------------------------------------
!!$        ! After decompostion, number of real tracers
!!$        ! on each process may have changed.
!!$        !----------------------------------------------------
!!$        
!!$        num_part_real = &
!!$             tracers_get_num_part_real(this,stat_info_sub)



!!$        !----------------------------------------------------
!!$        ! Compute forces on tracers at this time.
!!$        !----------------------------------------------------
!!$        
!!$        CALL tracers_compute_force(this, &  
!!$             stat_info_sub)
!!$        
!!$        IF ( stat_info_sub /= 0 ) THEN
!!$           PRINT *,"tracers_integrae_step : ", &
!!$                "Computing tracers' forces failed !"
!!$           stat_info = -1
!!$           GOTO 9999
!!$        END IF

     end if
        
        !----------------------------------------------------
        ! Return.
     	!----------------------------------------------------
        
9999	CONTINUE        
        
        
        RETURN
        
      END SUBROUTINE tracers_integrate_step
      
