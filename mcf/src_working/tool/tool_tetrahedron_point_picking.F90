FUNCTION tool_tetrahedron_point_picking(this, uni_rnd, tet_dim, stat_info)
  !----------------------------------------------------
  ! Function    : tool_tetrahedron_point_picking
  !----------------------------------------------------
  !
  ! Purpose     : return a point in a tetrahedron
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
  REAL(MK), DIMENSION(:,:),INTENT(IN)     :: tet_dim
  INTEGER, INTENT(OUT)                    :: stat_info

  REAL(MK), DIMENSION(SIZE(tet_dim, 1))   :: tool_tetrahedron_point_picking

  INTEGER                                 :: num_dim 
  REAL(MK), ALLOCATABLE, DIMENSION(:)     :: v0, v1, v2, v3
  REAL(MK)                                :: s, t, u, a, tmp
  INTEGER                                 :: stat_info_sub
  !----------------------------------------------------
  ! Initialization
  !
  ! This is supposed to be used, otherwise,
  ! compiler complains that it is not used.
  !----------------------------------------------------

  stat_info = this%flag
  stat_info = 0


  num_dim = SIZE(tet_dim, 1)

  ALLOCATE(v0(num_dim), v1(num_dim), v2(num_dim), v3(num_dim))


  v0 = tet_dim(1:num_dim, 1)
  v1 = tet_dim(1:num_dim, 2)
  v2 = tet_dim(1:num_dim, 3)
  v3 = tet_dim(1:num_dim, 4)

  ! generate three uniform random numbers
  CALL RANDOM_NUMBER(s)
  CALL RANDOM_NUMBER(t)
  CALL RANDOM_NUMBER(u)


  IF (s+t > 1.0_MK) THEN ! cut and fold the cube into a prism

     s = 1.0_MK - s
     t = 1.0_MK - t

  END IF

  IF (t+u > 1.0_MK) THEN ! cut and fold the prism into a tetrahedron

     tmp = u
     u = 1.0_MK - s - t
     t = 1.0_MK - tmp

  ELSE IF (s+t+u>1.0_MK) THEN

     tmp = u
     u = s + t + u - 1.0_MK
     s = 1.0_MK - t - tmp

  END IF

  a = 1.0_MK-s-t-u  ! a,s,t,u are the barycentric coordinates of the random point.

  tool_tetrahedron_point_picking(1:num_dim) =  v0(1:num_dim) * a + &
       v1(1:num_dim) * s + v2(1:num_dim) * t + v3(1:num_dim) * u


  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE

  RETURN

END FUNCTION tool_tetrahedron_point_picking
