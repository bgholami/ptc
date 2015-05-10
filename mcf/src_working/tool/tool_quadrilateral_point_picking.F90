FUNCTION tool_quadrilateral_point_picking(this, uni_rnd, box_dim, stat_info)
  !----------------------------------------------------
  ! Function    : tool_quadrilateral_point_picking
  !----------------------------------------------------
  !
  ! Purpose     : return a point in a quadrilateral
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
  REAL(MK), DIMENSION(:,:),INTENT(IN)     :: box_dim
  INTEGER, INTENT(OUT)                    :: stat_info

  REAL(MK), DIMENSION(SIZE(box_dim, 1))   :: tool_quadrilateral_point_picking

  INTEGER                                 :: num_dim 
  REAL(MK), ALLOCATABLE, DIMENSION(:,:,:) :: tc
  REAL(MK), ALLOCATABLE, DIMENSION(:)     :: t_area
  REAL(MK)                                :: rand
  INTEGER                                 :: choice_index
  INTEGER                                 :: stat_info_sub
  !----------------------------------------------------
  ! Initialization
  !
  ! This is supposed to be used, otherwise,
  ! compiler complains that it is not used.
  !----------------------------------------------------

  stat_info = this%flag
  stat_info = 0


  num_dim = SIZE(box_dim, 1)

  ALLOCATE(tc(2, num_dim, 3))
  ALLOCATE(t_area(2))

  ! create two triangles
  tc(1, 1:num_dim, 1) = box_dim(1:num_dim, 1)
  tc(1, 1:num_dim, 2) = box_dim(1:num_dim, 2)
  tc(1, 1:num_dim, 3) = box_dim(1:num_dim, 4)
  
  tc(2, 1:num_dim, 1) = box_dim(1:num_dim, 1)
  tc(2, 1:num_dim, 2) = box_dim(1:num_dim, 3)
  tc(2, 1:num_dim, 3) = box_dim(1:num_dim, 4)  

  ! compute surface area of each triangle
  t_area(1:2) = ABS( tc(1:2, 1, 1) * (tc(1:2, 2, 2)-tc(1:2, 2, 3)) + &
       tc(1:2, 1, 2) * (tc(1:2, 2, 3)-tc(1:2, 2, 1)) + &
       tc(1:2, 1, 3) * (tc(1:2, 2, 1)-tc(1:2, 2, 2)) ) / 2.0_MK
  ! and normalize area to later use for choosing
  t_area = t_area / SUM(t_area)

  ! choose one of the triangles
  CALL RANDOM_NUMBER(rand)
  
  IF (rand <= t_area(1)) THEN
     choice_index = 1
  ELSE
     choice_index = 2
  END IF

  ! create random point inside this triangle
  tool_quadrilateral_point_picking = tool_triangle_point_picking(this, &
       uni_rnd, tc(choice_index, 1:num_dim, 1:3), stat_info_sub)
  

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE

  RETURN

END FUNCTION tool_quadrilateral_point_picking
