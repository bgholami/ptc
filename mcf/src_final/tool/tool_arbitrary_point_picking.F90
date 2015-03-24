FUNCTION tool_arbitrary_point_picking(this, uni_rnd, box_dim, stat_info)
  !----------------------------------------------------
  ! Function    : tool_arbitrary_point_picking
  !----------------------------------------------------
  !
  ! Purpose     : return a point in a quadrilateral/wedge
  !               according to a uniform random number
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

  REAL(MK), DIMENSION(SIZE(box_dim, 1))   :: tool_arbitrary_point_picking

  INTEGER                                 :: num_dim
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

  IF (num_dim == 2) THEN

     tool_arbitrary_point_picking = tool_quadrilateral_point_picking(this, &
          uni_rnd, box_dim(1:num_dim, 1:2*num_dim), stat_info_sub)

  ELSE

     tool_arbitrary_point_picking = tool_wedge_point_picking(this, &
          uni_rnd, box_dim(1:num_dim, 1:2*num_dim), stat_info_sub)

  END IF

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE

  RETURN

END FUNCTION tool_arbitrary_point_picking


