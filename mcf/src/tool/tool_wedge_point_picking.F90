FUNCTION tool_wedge_point_picking(this, uni_rnd, box_dim, stat_info)
  !----------------------------------------------------
  ! Function    : tool_wedge_point_picking
  !----------------------------------------------------
  !
  ! Purpose     : return a point in a wedge
  !
  ! Routines    :
  !
  ! Remarks     : wedge is a triangular prism with non-parallel bases
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

  REAL(MK), DIMENSION(SIZE(box_dim, 1))   :: tool_wedge_point_picking

  INTEGER                                 :: num_dim 
  REAL(MK), ALLOCATABLE, DIMENSION(:,:,:) :: tc
  REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: dum_mat
  REAL(MK), ALLOCATABLE, DIMENSION(:)     :: t_vol
  REAL(MK)                                :: rand, det
  INTEGER                                 :: choice_index, i
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

  ALLOCATE(tc(3, num_dim, 4))
  ALLOCATE(t_vol(3))
  ALLOCATE(dum_mat(num_dim, 3))

  ! create three tetrahedra (assumption: vertices {1,2,3} and {4,5,6} form
  ! two bases. vertex 1 is connected to 4, v2 to v5, and v3 to v6.
  tc(1, 1:num_dim, 1) =  box_dim(1:num_dim, 1)
  tc(1, 1:num_dim, 2) =  box_dim(1:num_dim, 2)  
  tc(1, 1:num_dim, 3) =  box_dim(1:num_dim, 3)  
  tc(1, 1:num_dim, 4) =  box_dim(1:num_dim, 6)

  tc(2, 1:num_dim, 1) =  box_dim(1:num_dim, 1)
  tc(2, 1:num_dim, 2) =  box_dim(1:num_dim, 4)  
  tc(2, 1:num_dim, 3) =  box_dim(1:num_dim, 5)  
  tc(2, 1:num_dim, 4) =  box_dim(1:num_dim, 6) 

  tc(3, 1:num_dim, 1) =  box_dim(1:num_dim, 1)
  tc(3, 1:num_dim, 2) =  box_dim(1:num_dim, 2)  
  tc(3, 1:num_dim, 3) =  box_dim(1:num_dim, 5)  
  tc(3, 1:num_dim, 4) =  box_dim(1:num_dim, 6)


  !compute volume of each tetrahedron
  DO i = 1, 3  
     dum_mat(1:num_dim, 1) = tc(i, 1:num_dim, 1) - tc(i, 1:num_dim, 4)
     dum_mat(1:num_dim, 2) = tc(i, 1:num_dim, 2) - tc(i, 1:num_dim, 4)
     dum_mat(1:num_dim, 3) = tc(i, 1:num_dim, 3) - tc(i, 1:num_dim, 4)

     ! calculate the determinant
     det = dum_mat(1, 1) * (dum_mat(2, 2)*dum_mat(3, 3) - dum_mat(2, 3)*dum_mat(3, 2)) - &
          dum_mat(1, 2) * (dum_mat(2, 1)*dum_mat(3, 3) - dum_mat(2, 3)*dum_mat(3, 1)) + &
          dum_mat(1, 3) * (dum_mat(2, 1)*dum_mat(3, 2) - dum_mat(2, 2)*dum_mat(3, 1))

     t_vol(i) = ABS(det) / 6.0_MK
  END DO

  ! and normalize volume to later use for choosing
  t_vol = t_vol / SUM(t_vol)

  ! choose one
  CALL RANDOM_NUMBER(rand)
  IF (rand <= t_vol(1)) THEN
     choice_index = 1
  ELSE IF (rand <= t_vol(1) + t_vol(2)) THEN
     choice_index = 2
  ELSE
     choice_index = 3
  END IF

  ! create random point inside this tetrahedron
  tool_wedge_point_picking = tool_tetrahedron_point_picking(this, &
       uni_rnd, tc(choice_index, 1:num_dim, 1:4), stat_info_sub)


  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE

  RETURN

END FUNCTION tool_wedge_point_picking
