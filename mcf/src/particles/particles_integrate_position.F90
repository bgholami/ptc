      SUBROUTINE particles_integrate_position(this,&
           num,dt,accuracy_order,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_position
        !----------------------------------------------------
        !
        ! Purpose     : Integrate position of particles
        !               with required accuracy.
        !
        ! Reference   :
        !
        ! Remark      : colloidal boundary particle may rotate,
        !               needs to be done seperately.
        !               boundary condition boundary particles 
        !               can be done here.
        !
        !
        ! Revision    : V0.3 24.08.2010, integrate position of
        !               particles, except colloidal boundary 
        !               particles.
        !
        !               V0.2  08.07.2009,
        !               check again the work flow and supply
        !               with more comments.
        !
        !               V0.1  01.04.2009, original version.
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
        !  Arguments
        !
        !  this           : an object of Particles Class.
        !  num            : number of particles needed to be updated,
        !                   i.e. first num particles in this%x 
        !                   are operated.
        !  dt             : time step.
        !  accuracy_order : accuracy required.
        !  stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        REAL(MK),INTENT(IN)                     :: dt
        INTEGER, INTENT(IN)                     :: accuracy_order
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        INTEGER                                 :: stat_info_sub 
        TYPE(Boundary), POINTER                 :: d_boundary
        REAL(MK), DIMENSION(:), POINTER         :: min_phys
        REAL(MK), DIMENSION(:), POINTER         :: max_phys
        INTEGER                                 :: dim, i, j, patch_id
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: new_x
        REAL(MK)                                :: distance
        LOGICAL                                 :: eliminate_particle

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------

        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(d_boundary)
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        
        dim           = &
             physics_get_num_dim(this%phys,stat_info_sub)
        CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)
        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        ALLOCATE(new_x(dim))

#ifdef __POSITION_FIXED

#else      
        !----------------------------------------------------
        ! Select different accuracy oder
        ! 1 only velocity contribution.
        ! 2 velocity + acceleration contribution. (transport velocity)
        !----------------------------------------------------

        Do i = 1, num

           eliminate_particle = .FALSE.
           IF ( this%id(this%sid_idx,i) == mcf_particle_type_fluid) THEN

              SELECT CASE (accuracy_order)

              CASE (1)

                 new_x(1:dim) = &
                      this%x(1:dim,i) + &
                      this%v(1:dim,i) * dt

              CASE (2)

                 new_x(1:dim) = &
                      this%x(1:dim,i) + &
                      this%v(1:dim,i) * dt  + &
                      0.5_MK * this%f_bp(1:dim,i) * dt**2

                 !-------------------------------------------------
                 ! Other integration not available.
                 !-------------------------------------------------

              CASE DEFAULT
                 PRINT *, "particles_integrate_position : ", &
                      "Order of accuracy not available !"
                 stat_info = -1
                 GOTO 9999

              END SELECT ! accuracy_order


              ! getting rid of those buffer particles that have
              ! unphysical force/velocity/position 
              IF (mcf_eliminate_non_physical_buffer_particles) THEN
               
                 distance = 0.0_MK
                 patch_id = 0
                 CALL boundary_check_particle_stat(d_boundary, &
                      this%x(1:dim, i), 0, distance, patch_id, stat_info_sub)
                 ! only buffer particles
                 IF (distance /= 0.0_MK) THEN
                    DO j = 1, num
                       IF ((new_x(j) .LT. min_phys(j)) .OR. &
                            (new_x(j) .GT. max_phys(j))) THEN
                          ! eliminate
                          eliminate_particle = .TRUE.
                          this%id(this%pid_idx, i) = -ABS(this%id(this%pid_idx, i))
                       END IF
                       EXIT
                    END DO
                 END IF

              END IF


              IF (.NOT. eliminate_particle) THEN
                 this%x(1:dim,i) = new_x(1:dim)
              END IF

           END IF

        END DO

#endif
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_integrate_position
