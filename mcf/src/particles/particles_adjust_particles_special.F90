SUBROUTINE particles_adjust_particles_special(this,num,stat_info)
  !----------------------------------------------------
  ! Subroutine  : particles_adjust_particles_special
  !----------------------------------------------------
  !
  ! Purpose     : Adjust particles positions and 
  !               velocities according to boundary
  !               conditions.
  !                
  !	 	      	 
  ! Reference   :
  !              
  ! Remark      : 
  !
  ! Revisions   : 
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
  ! Arguments
  !
  ! this       : an object of Particles Class.
  ! num        : number of particles needed to be updated,
  !              i.e. first 'num' particles in this%x are 
  !              operated.
  ! stat_info  : return flag of status.
  !----------------------------------------------------

  TYPE(Particles), INTENT(INOUT)  :: this
  INTEGER, INTENT(IN)             :: num
  INTEGER, INTENT(OUT)	          :: stat_info

  !----------------------------------------------------
  ! Local variables
  !----------------------------------------------------

  INTEGER                         :: num_dim
  TYPE(Boundary), POINTER         :: d_boundary
  REAL(MK)                        :: bwidth, distance
  INTEGER                         :: j, patch_id
  INTEGER                         :: num_create, num_delete
  INTEGER                         :: stat_info_sub
#if __TRACER
  INTEGER                         :: tracers_c_factor
#endif

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0


  NULLIFY(d_boundary)

  num_create = 0
  num_delete = 0


  !----------------------------------------------------
  ! Check if num is in the range of particles' number,
  ! we are supposed to check only the real particles.
  !----------------------------------------------------

  IF( num > this%num_part_real ) THEN
     PRINT *, "particles_adjust_particles_adjust :", &
          "Num > number of real particles !"
     stat_info = -1
     GOTO 9999
  END IF

  !----------------------------------------------------
  ! Get physics quantities including boundary.
  !----------------------------------------------------

  num_dim = physics_get_num_dim(this%phys,stat_info_sub)
  CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)
  bwidth = boundary_get_bwidth(d_boundary, stat_info_sub)
#if __TRACER
  tracers_c_factor = &
       physics_get_tracers_c_factor(this%phys,stat_info_sub)
#endif
  


  !----------------------------------------------------
  ! Loop over the first num particles.
  !----------------------------------------------------

  DO j = 1, num    

     IF (this%id(this%sid_idx, j) == mcf_particle_type_fluid) THEN

        distance = 0.0_MK
        patch_id = 0
        CALL boundary_check_particle_stat(d_boundary, &
             this%x(1:num_dim, j), 1, distance, patch_id, stat_info_sub)


        IF (distance == 0.0_MK) THEN

           ! interior
           IF (this%id(this%bid_idx, j) > 0) THEN 
              ! former inlet that just entered interior
              num_create = num_create + 1
           END IF

        ELSE IF (distance > 0.0_MK) THEN

           IF (ABS(distance) <= bwidth) THEN   

              ! inlet
              ! optional? to cover interiors that come back to inlet buffer...
              !this%id(this%bid_idx, j) = patch_id
              CONTINUE

           ELSE

              ! delete  
              ! mark with change in id value (delete all in the end)
              this%id(this%pid_idx, j) = -ABS(this%id(this%pid_idx, j))
              num_delete = num_delete + 1

           END IF

        ELSE IF (distance < 0.0_MK) THEN

           IF (ABS(distance) <= bwidth) THEN   

              ! outlet
              this%id(this%bid_idx, j) = patch_id

           ELSE

              ! delete  
              ! mark with change in id value (delete all in the end)
              this%id(this%pid_idx, j) = -ABS(this%id(this%pid_idx, j))
              num_delete = num_delete + 1

           END IF

        END IF

     END IF

  END DO ! j = 1, num of particles


  ! take care of create/delete (for inflow/outflow BC)
  CALL particles_create_delete_particles(this, num_create, num_delete, stat_info_sub)


9999 CONTINUE	


  RETURN

END SUBROUTINE particles_adjust_particles_special
