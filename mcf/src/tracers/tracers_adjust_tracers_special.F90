SUBROUTINE tracers_adjust_tracers_special(this,num,stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_adjust_tracers_special
  !----------------------------------------------------
  !
  ! Purpose     : Adjust tracers positions and 
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
  INTEGER, INTENT(OUT)	        :: stat_info

  !----------------------------------------------------
  ! Local variables
  !----------------------------------------------------

  INTEGER                         :: num_dim
  TYPE(Boundary), POINTER         :: d_boundary  
  REAL(MK)                        :: distance
  INTEGER                         :: j, patch_id
  INTEGER                         :: stat_info_sub

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0

  NULLIFY(d_boundary)


  !----------------------------------------------------
  ! Check if num is in the range of particles' number,
  ! we are supposed to check only the real particles.
  !----------------------------------------------------

  IF( num > this%num_part_real ) THEN   
     PRINT *, "tracers_adjust_tracers_special :", &
          "Num > number of real particles !"
     stat_info = -1
     GOTO 9999
  END IF


  !----------------------------------------------------
  ! Get physics quantities including boundary.
  !----------------------------------------------------

  num_dim = physics_get_num_dim(this%phys,stat_info_sub)
  CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)


  !----------------------------------------------------
  ! Loop over the first num particles.
  !----------------------------------------------------

  DO j = 1, num    

     IF (this%id(2, j) == 0) THEN     

        distance = 0.0_MK
        patch_id = 0
        CALL boundary_check_particle_stat(d_boundary, &
             this%x(1:num_dim, j), 0, distance, patch_id, stat_info_sub)


        IF (distance == 0.0_MK) THEN

           ! interior
           CONTINUE

        ELSE 

           ! outside the interior
           this%id(1, j) = 0  ! delete it
           this%id(2, j) = -1 ! and delete without insertion into SPH particle

        END IF

     END IF

  END DO ! j = 1, num of particles


9999 CONTINUE	


  RETURN

END SUBROUTINE tracers_adjust_tracers_special
