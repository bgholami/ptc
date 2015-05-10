SUBROUTINE boundary_check_particle_stat(this, &
     p_x, mode, distance, patch, stat_info)
  !----------------------------------------------------
  ! Subroutine  : boundary_check_particle_stat
  !----------------------------------------------------
  !
  ! Purpose     : Check if p_x is a interior/inflow/outflow 
  !               boundary particle.
  !               returns distance = 0 for interior (and 0 for patch ID)
  !                       distance > 0 for inlet (and its patch ID) 
  !                       distance < 0 for outlet (and its patch ID)
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
  ! Arguments :
  !
  ! p_x       : position.
  ! distance  : distance
  ! stat_info : status of this routine.
  !----------------------------------------------------

  TYPE(Boundary), INTENT(INOUT)           :: this
  REAL(MK), DIMENSION(:), INTENT(IN)      :: p_x
  INTEGER, INTENT(IN)                     :: mode
  REAL(MK), INTENT(OUT)                   :: distance  
  INTEGER , INTENT(OUT)                   :: patch
  INTEGER, INTENT(OUT)                    :: stat_info


  !----------------------------------------------------
  ! Local variables
  !----------------------------------------------------

  INTEGER                                 :: num_dim, min_patch
  REAL(MK)                                :: dist, max_dist, min_dist
  REAL(MK)                                :: proj_dist, normal_dist, wall_width
  REAL(MK), DIMENSION(3)                  :: r, distList
  REAL(MK), ALLOCATABLE, DIMENSION(:)     :: proj_vec
  INTEGER                                 :: j

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info = 0
  distance = 0.0_MK
  patch = 0
  num_dim   = this%num_dim
  wall_width = this%bwidth / mcf_IOghost_layer_coeff  * mcf_wall_layer_coeff

  IF (this%num_inout > 0) THEN

     ALLOCATE(proj_vec(this%num_inout))
     IF (mode == -1) THEN
        ! find out if projection in patch
        min_dist = 10000.0_MK * this%bwidth
        min_patch = 0
        DO j = 1, this%num_inout
           r = p_x(1:num_dim) - this%ref_point(j, 1:num_dim)
           proj_vec(1:num_dim) = r(1:num_dim) - &
                this%iopatch_n(j, 1:num_dim) * DOT_PRODUCT(r(1:num_dim), this%iopatch_n(j, 1:num_dim))
           proj_dist =  SQRT(DOT_PRODUCT(proj_vec(1:num_dim), proj_vec(1:num_dim)))
           
           IF (proj_dist .LE. (this%clength(j) + wall_width)) THEN
              normal_dist = DOT_PRODUCT(r(1:num_dim), this%iopatch_n(j, 1:num_dim))
              ! we are looking for the smallest negative distance
              IF ((normal_dist .LT. 0.0_MK) .AND. (ABS(normal_dist) .LT. ABS(min_dist))) THEN
                 min_dist = normal_dist
                 min_patch = j
              END IF
           END IF

        END DO

        IF (min_patch == 0) THEN
           distance = 0.0_MK
           patch = 0
        ELSE
           ! this make the output distance positive for inlet, negative for outlet
           distance = SIGN(min_dist, DBLE(this%patch_id(min_patch)))
           patch    = this%patch_id(min_patch)

           IF (ABS(distance) <= this%ewidth) THEN
              distance = 0.0_MK
              patch    = 0
           ELSE
              distance = SIGN(ABS(distance) - this%ewidth, & 
                   DBLE(this%patch_id(min_patch)))
           END IF

        END IF

     ELSE

        !----------------------------------------------------
        ! find which patch it blongs to
        !----------------------------------------------------

        min_patch = 0
        DO j = 1, this%num_inout

           r = p_x(1:num_dim) - this%ref_point(j, 1:num_dim)
           distList(j) = SQRT(DOT_PRODUCT(r(1:num_dim), r(1:num_dim)))

        END DO

        j = MINLOC(distList, 1)
        !max_dist = SQRT((2.0_MK * this%clength(j) + 2.0_MK * this%bwidth)**2.0_MK + (this%bwidth + this%ewidth)**2.0_MK)
        max_dist = SQRT((this%clength(j) + wall_width)**2.0_MK + (this%bwidth + this%ewidth)**2.0_MK)
        
        IF (distList(j)  <= max_dist) THEN  
           min_patch = j
        END IF

        IF (min_patch == 0) THEN

           distance = 0.0_MK  
           patch    = 0

        ELSE

           r = p_x(1:num_dim) - this%ref_point(min_patch, 1:num_dim)
           dist = DOT_PRODUCT(r(1:num_dim), this%iopatch_n(min_patch, 1:num_dim))

           IF (dist >= 0.0_MK) THEN       

              distance = 0.0_MK
              patch    = 0

           ELSE
              ! this make the output distance positive for inlet, negative for outlet
              distance = SIGN(dist, DBLE(this%patch_id(min_patch)))
              patch    = this%patch_id(min_patch)

              IF (mode > 0) THEN
                 IF (ABS(distance) <= this%ewidth) THEN
                    distance = 0.0_MK
                    patch    = 0
                 ELSE
                    distance = SIGN(ABS(distance) - this%ewidth, & 
                         DBLE(this%patch_id(min_patch)))
                 END IF
              END IF

           END IF

        END IF

     END IF

  END IF


9999 CONTINUE


  RETURN

END SUBROUTINE boundary_check_particle_stat

