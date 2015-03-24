FUNCTION tool_triangle_point_picking(this, uni_rnd, tri_dim, stat_info)
  !----------------------------------------------------
  ! Function    : tool_triangle_point_picking
  !----------------------------------------------------
  !
  ! Purpose     : return a point in a triangle
  !
  ! Routines    :
  !
  ! Remarks     : 
  !
  ! References  :
  !
  ! Revisions   : V0.1 05.02 2013, original version.
  !
  !----------------------------------------------------
  ! Author       : Babak Gholami
  ! Contact      : babak.gholami@aer.mw.tum.de
  !
  ! Dr. Marco Ellero's Emmy Noether Group,
  ! Prof. Dr. N. Adams' Chair of Aerodynamics,
  ! Faculty of Mechanical Engineering,
  ! Technische Universitaet Muenchen, Germany.
  !----------------------------------------------------

  !----------------------------------------------------
  ! Arguments     
  !----------------------------------------------------
  
  TYPE(Tool), INTENT(IN)                  :: this
  REAL(MK)    ,INTENT(IN)                 :: uni_rnd
  REAL(MK), DIMENSION(:,:),INTENT(IN)     :: tri_dim
  INTEGER, INTENT(OUT)                    :: stat_info

  REAL(MK), DIMENSION(SIZE(tri_dim, 1))   :: tool_triangle_point_picking

  INTEGER                                 :: num_dim, i 
  REAL(MK), ALLOCATABLE, DIMENSION(:)     :: lambda, x
  REAL(MK)                                :: a0, a1
  LOGICAL                                 :: in_test
  INTEGER                                 :: stat_info_sub
  !----------------------------------------------------
  ! Initialization
  !
  ! This is supposed to be used, otherwise,
  ! compiler complains that it is not used.
  !----------------------------------------------------

  stat_info = this%flag
  stat_info = 0


  num_dim = SIZE(tri_dim, 1)

  ALLOCATE(lambda(3)) ! 3 for triangle and 2 for line
  ALLOCATE(x(num_dim))

  ! generate two uniform random numbers
  CALL RANDOM_NUMBER(a0)
  CALL RANDOM_NUMBER(a1)

  ! compute point x
  x = tri_dim(1:num_dim, 1) + &
       a0 * (tri_dim(1:num_dim, 2) - tri_dim(1:num_dim, 1)) + &
       a1 * (tri_dim(1:num_dim, 3) - tri_dim(1:num_dim, 1))

  ! calculate barycentric coordinates
  CALL tool_barycentric_coordinate(this, tri_dim, x, lambda, stat_info_sub)

  ! determine whether x is inside the triangle
  in_test = .FALSE.
  DO i = 1, num_dim
     IF ( (lambda(i) > 0.0_MK) .AND. (lambda(i) < 1.0_MK) ) THEN
        in_test = .TRUE.
     ELSE
        in_test = .FALSE.
        EXIT
     END IF
  END DO


  IF (in_test) THEN

     tool_triangle_point_picking(1:num_dim) = x(1:num_dim)

  ELSE

     tool_triangle_point_picking(1:num_dim) = &
          tri_dim(1:num_dim, 2) + tri_dim(1:num_dim, 3) - x

  END IF


  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE

  RETURN

END FUNCTION tool_triangle_point_picking
