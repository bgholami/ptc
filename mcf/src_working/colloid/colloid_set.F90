!--------------------------------------------------
! Subroutine  : colloid_set_*
!--------------------------------------------------
!
! Purpose     : Set routines of Class colloid.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
! Author      : Xin Bian
! Contact     : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------

SUBROUTINE colloid_set_tech(this,d_tech,stat_info)
  TYPE(Colloid), INTENT(INOUT)            :: this
  TYPE(Technique),INTENT(IN),TARGET       :: d_tech
  INTEGER, INTENT(OUT)                    :: stat_info

  stat_info = 0
  this%tech => d_tech

  RETURN
END SUBROUTINE colloid_set_tech


SUBROUTINE colloid_set_num_dim(this,d_num_dim, stat_info)
  !-----------------------------------------
  ! Set the number of dimension of colloids.
  !-----------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_num_dim
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  !---------------------------------------
  ! Only 2D, 3D are supported
  !---------------------------------------

  IF(d_num_dim < 2 .OR. d_num_dim > 3 ) THEN
     PRINT *, "colloid_set_num_dim : ", &
          "Dimension is not supported !"
     stat_info = -1
     GOTO 9999
  END IF

  !---------------------------------------
  ! If reset the dimension, 
  ! the memory has to be reallocated
  !---------------------------------------
  IF( d_num_dim /= this%num_dim ) THEN

     this%num_dim = d_num_dim

     IF(ASSOCIATED(this%x))THEN
        DEALLOCATE(this%x)
     END IF
     ALLOCATE(this%x(d_num_dim,this%num_colloid))

     IF(ASSOCIATED(this%v))THEN
        DEALLOCATE(this%v)
     END IF
     ALLOCATE(this%v(d_num_dim,this%num_colloid))

#if __DRAG_PART
     IF(ASSOCIATED(this%drag_lub))THEN
        DEALLOCATE(this%drag_lub)
     END IF
     ALLOCATE(this%drag_lub(d_num_dim,this%num_colloid))

     IF(ASSOCIATED(this%drag_repul))THEN
        DEALLOCATE(this%drag_repul)
     END IF
     ALLOCATE(this%drag_repul(d_num_dim,this%num_colloid))
#endif

     IF(ASSOCIATED(this%drag))THEN
        DEALLOCATE(this%drag)
     END IF
     ALLOCATE(this%drag(d_num_dim,this%num_colloid))

     IF(ASSOCIATED(this%body_force))THEN
        DEALLOCATE(this%body_force)
     END IF
     ALLOCATE(this%body_force(d_num_dim))

     IF(ASSOCIATED(this%f))THEN
        DEALLOCATE(this%f)
     END IF
     ALLOCATE(this%f(d_num_dim,this%num_colloid))   

     IF(ASSOCIATED(this%mom))THEN
        DEALLOCATE(this%mom)
     END IF
     ALLOCATE(this%mom(d_num_dim,this%num_colloid))  

     IF(ASSOCIATED(this%min_phys)) THEN
        DEALLOCATE(this%min_phys)
     END IF
     ALLOCATE(this%min_phys(d_num_dim))

     IF(ASSOCIATED(this%max_phys)) THEN
        DEALLOCATE(this%max_phys)
     END IF
     ALLOCATE(this%max_phys(d_num_dim))

     IF(ASSOCIATED(this%bcdef)) THEN
        DEALLOCATE(this%bcdef)
     END IF
     ALLOCATE(this%bcdef(2*d_num_dim))

  END IF

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_num_dim


SUBROUTINE colloid_set_num_colloid(this,d_num_colloid, stat_info)
  !-----------------------------------------
  ! Set the number of colloids.
  !-----------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_num_colloid
  INTEGER, INTENT(OUT)            :: stat_info

  INTEGER                         :: dim


  stat_info = 0

  dim = this%num_dim

  !---------------------------------------
  ! If reset the number of colloids, 
  ! the memory has to be reallocated
  !---------------------------------------

  IF( d_num_colloid /= this%num_colloid ) THEN

     this%num_colloid = d_num_colloid

     IF(ASSOCIATED(this%shape))THEN
        DEALLOCATE(this%shape)
     END IF
     ALLOCATE(this%shape(d_num_colloid))

     IF(ASSOCIATED(this%radius))THEN
        DEALLOCATE(this%radius)
     END IF
     ALLOCATE(this%radius(dim,d_num_colloid))

     IF(ASSOCIATED(this%freq))THEN
        DEALLOCATE(this%freq)
     END IF
     ALLOCATE(this%freq(d_num_colloid))

     IF(ASSOCIATED(this%m))THEN
        DEALLOCATE(this%m)
     END IF
     ALLOCATE(this%m(d_num_colloid))

     IF(ASSOCIATED(this%mmi))THEN
        DEALLOCATE(this%mmi)
     END IF
     ALLOCATE(this%mmi(d_num_colloid))

     IF(ASSOCIATED(this%x))THEN
        DEALLOCATE(this%x)
     END IF
     ALLOCATE(this%x(dim,d_num_colloid))

     IF(ASSOCIATED(this%v))THEN
        DEALLOCATE(this%v)
     END IF
     ALLOCATE(this%v(dim,d_num_colloid))

#if __DRAG_PART
     IF(ASSOCIATED(this%drag_lub))THEN
        DEALLOCATE(this%drag_lub)
     END IF
     ALLOCATE(this%drag_lub(dim,d_num_colloid))

     IF(ASSOCIATED(this%drag_repul))THEN
        DEALLOCATE(this%drag_repul)
     END IF
     ALLOCATE(this%drag_repul(dim,d_num_colloid))
#endif           

     IF(ASSOCIATED(this%drag))THEN
        DEALLOCATE(this%drag)
     END IF
     ALLOCATE(this%drag(dim,d_num_colloid))

     IF(ASSOCIATED(this%phi))THEN
        DEALLOCATE(this%phi)
     END IF
     ALLOCATE(this%phi(3,d_num_colloid))

     IF(ASSOCIATED(this%theta))THEN
        DEALLOCATE(this%theta)
     END IF
     ALLOCATE(this%theta(3,d_num_colloid))

     IF(ASSOCIATED(this%omega))THEN
        DEALLOCATE(this%omega)
     END IF
     ALLOCATE(this%omega(3,d_num_colloid))

     IF(ASSOCIATED(this%torque))THEN
        DEALLOCATE(this%torque)
     END IF
     ALLOCATE(this%torque(3,d_num_colloid))

     IF(ASSOCIATED(this%num_physical_part))THEN
        DEALLOCATE(this%num_physical_part)
     END IF
     ALLOCATE(this%num_physical_part(d_num_colloid))

     IF(ASSOCIATED(this%num_numerical_part))THEN
        DEALLOCATE(this%num_numerical_part)
     END IF
     ALLOCATE(this%num_numerical_part(d_num_colloid))

     IF(ASSOCIATED(this%f))THEN
        DEALLOCATE(this%f)
     END IF
     ALLOCATE(this%f(dim,d_num_colloid))

     IF(ASSOCIATED(this%alpha))THEN
        DEALLOCATE(this%alpha)
     END IF
     ALLOCATE(this%alpha(3,d_num_colloid))

     IF(ASSOCIATED(this%k_energy))THEN
        DEALLOCATE(this%k_energy)
     END IF
     ALLOCATE(this%k_energy(d_num_colloid))

     IF(ASSOCIATED(this%mom))THEN
        DEALLOCATE(this%mom)
     END IF
     ALLOCATE(this%mom(dim,d_num_colloid))

  END IF

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_num_colloid


SUBROUTINE colloid_set_rho(this,d_rho,stat_info)
  !-----------------------------------------
  ! Set the colloid density.
  !-----------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_rho
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_rho < 0.0_MK  ) THEN
     PRINT *, "colloid_set_rho : ", &
          "Wrong rho !"
     stat_info = -1
  END IF

  this%rho = d_rho

  RETURN

END SUBROUTINE colloid_set_rho


SUBROUTINE colloid_set_rho_type(this,d_rho_type, stat_info)
  !-----------------------------------------
  ! Set the colloid density type.
  !-----------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_rho_type
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_rho_type < 0 .OR. d_rho_type > 1 ) THEN
     PRINT *, "colloid_set_rho_type : ", &
          "Wrong rho type !"
     stat_info = -1
  END IF

  this%rho_type = d_rho_type

  RETURN

END SUBROUTINE colloid_set_rho_type


SUBROUTINE colloid_set_translate(this,d_translate,stat_info)
  !---------------------------------------
  ! Set if colloids can translate or not.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  LOGICAL, INTENT(IN)             :: d_translate
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%translate = d_translate

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_translate


SUBROUTINE colloid_set_rotate(this,d_rotate,stat_info)
  !---------------------------------------
  ! Set if colloids can rotate or not.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  LOGICAL, INTENT(IN)             :: d_rotate
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%rotate = d_rotate

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_rotate


SUBROUTINE colloid_set_place(this,d_place,stat_info)
  !---------------------------------------
  ! Set the type of noslip condition on
  ! the surface of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_place
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%place = d_place

9999 CONTINUE        
  RETURN       

END SUBROUTINE colloid_set_place


SUBROUTINE colloid_set_noslip_type(this,d_noslip_type,stat_info)
  !---------------------------------------
  ! Set the type of noslip condition on
  ! the surface of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_noslip_type
  INTEGER, INTENT(OUT)            :: stat_info


  stat_info = 0

  this%noslip_type = d_noslip_type

9999 CONTINUE        
  RETURN       

END SUBROUTINE colloid_set_noslip_type


SUBROUTINE colloid_set_body_force_type(this,d_type,stat_info)
  !---------------------------------------
  ! Set the type of body force on
  ! colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_type
  INTEGER, INTENT(OUT)            :: stat_info


  stat_info = 0

  this%body_force_type = d_type

9999 CONTINUE        
  RETURN       

END SUBROUTINE colloid_set_body_force_type


SUBROUTINE colloid_set_body_force(this,d_body_force,stat_info)
  !---------------------------------------
  ! Set the body force on all colloids
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:), INTENT(IN)      :: d_body_force
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim

  stat_info = 0

  !---------------------------------------
  ! Check if the input dimension matches.
  !---------------------------------------

  dim = SIZE(d_body_force,1)

  IF( dim /= this%num_dim) THEN
     PRINT *, "colloid_set_body_force : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  this%body_force(1:dim) = d_body_force(1:dim)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_body_force



SUBROUTINE colloid_set_cc_lub_type(this,d_type,stat_info)
  !----------------------------------------------------
  ! Set lubrication interaction type between 
  ! colloid and colloid.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_type
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_type < 0 .OR. d_type > 2 ) THEN
     PRINT *, "colloid_set_cc_lub_type : ", &
          "Wrong lub type !"
     stat_info = -1
     GOTO 9999
  END IF

  this%cc_lub_type = d_type

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_lub_type


SUBROUTINE colloid_set_cc_repul_type(this,d_type,stat_info)
  !----------------------------------------------------
  ! Set repulsive force interaction type between 
  ! colloid and colloid.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_type
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_type < 0 .OR. d_type > 2 ) THEN
     PRINT *, "colloid_set_cc_repul_type : ", &
          "Wrong repul type !"
     stat_info = -1
     GOTO 9999
  END IF

  this%cc_repul_type = d_type

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_repul_type


SUBROUTINE colloid_set_cc_lub_cut_off(this,d_cut_off,stat_info)
  !----------------------------------------------------
  ! Set cut off for lubrication interaction 
  ! between colloid and colloid.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_off
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cc_lub_cut_off = d_cut_off

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_lub_cut_off


SUBROUTINE colloid_set_cc_lub_cut_on(this,d_cut_on,stat_info)
  !----------------------------------------------------
  ! Set cut on for lubrication interaction 
  ! between colloid and colloid.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_on
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cc_lub_cut_on = d_cut_on

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_lub_cut_on


SUBROUTINE colloid_set_cc_repul_cut_off(this,d_cut_off,stat_info)
  !----------------------------------------------------
  ! Set cut off of repulsive interaction 
  ! between colloid and colloid.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_off
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cc_repul_cut_off = d_cut_off

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_repul_cut_off


SUBROUTINE colloid_set_cc_repul_cut_on(this,d_cut_on,stat_info)
  !----------------------------------------------------
  ! Set cut on of repulsive interaction 
  ! between colloid and colloid.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_on
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cc_repul_cut_on = d_cut_on

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_repul_cut_on


SUBROUTINE colloid_set_cc_repul_F0(this,d_F0,stat_info)
  !----------------------------------------------------
  ! Set maximum repulsive force between colloids.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_F0
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cc_repul_F0 = d_F0

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cc_repul_F0


SUBROUTINE colloid_set_cw_lub_type(this,d_type, stat_info)
  !----------------------------------------------------
  ! Set lubrication interaction type between 
  ! colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_type
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_type < 0 .OR. d_type > 2 ) THEN
     PRINT *, "colloid_set_cw_lub_type : ", &
          "Wrong lub type !"
     stat_info = -1
     GOTO 9999
  END IF

  this%cw_lub_type = d_type

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_lub_type


SUBROUTINE colloid_set_cw_repul_type(this,d_type,stat_info)
  !----------------------------------------------------
  ! Set repulsive force interaction type between 
  ! colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_type
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_type < 0 .OR. d_type > 2 ) THEN
     PRINT *, "colloid_set_cw_repul_type : ", &
          "Wrong repul type !"
     stat_info = -1
     GOTO 9999
  END IF

  this%cw_repul_type = d_type

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_repul_type


SUBROUTINE colloid_set_cw_lub_cut_off(this,d_cut_off, stat_info)
  !----------------------------------------------------
  ! Set cut off for lubrication interaction 
  ! between colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_off
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cw_lub_cut_off = d_cut_off

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_lub_cut_off


SUBROUTINE colloid_set_cw_lub_cut_on(this,d_cut_on,stat_info)
  !----------------------------------------------------
  ! Set cut on for lubrication interaction 
  ! between colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_on
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cw_lub_cut_on = d_cut_on

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_lub_cut_on


SUBROUTINE colloid_set_cw_repul_cut_off(this,d_cut_off,stat_info)
  !----------------------------------------------------
  ! Set cut off of repulsive interaction 
  ! between colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_off
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cw_repul_cut_off = d_cut_off

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_repul_cut_off


SUBROUTINE colloid_set_cw_repul_cut_on(this,d_cut_on,stat_info)
  !----------------------------------------------------
  ! Set cut on of repulsive interaction 
  ! between colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_on
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cw_repul_cut_on = d_cut_on

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_repul_cut_on


SUBROUTINE colloid_set_cw_repul_F0(this,d_F0,stat_info)
  !----------------------------------------------------
  ! Set maximum repulsive force between colloid and wall.
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_F0
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  this%cw_repul_F0 = d_F0

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cw_repul_F0


SUBROUTINE colloid_set_shape(this, d_shape,stat_info)
  !---------------------------------------
  ! Set the shape of a colloid object
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  INTEGER, DIMENSION(:), INTENT(IN)       :: d_shape
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num

  stat_info = 0

  num = SIZE(d_shape,1)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_shape : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%shape(:) = d_shape(1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_shape


SUBROUTINE colloid_set_radius(this,d_ra,stat_info)
  !---------------------------------------
  ! Set the radius of a colloid object 
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_ra
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim,num

  stat_info = 0

  dim = SIZE(d_ra,1)
  num = SIZE(d_ra,2)

  IF (dim /= this%num_dim) THEN
     PRINT *, "colloid_set_radius : ", &
          "Number of dimension doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_ra : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%radius(:,:) = d_ra(:,:)

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_radius


SUBROUTINE colloid_set_freq(this,d_freq,stat_info)
  !---------------------------------------
  ! Set frequency of surface roughness.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  INTEGER, DIMENSION(:),INTENT(IN)        :: d_freq
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num

  stat_info = 0

  num = SIZE(d_freq,1)

  IF ( num /= this%num_colloid ) THEN
     PRINT *, "colloid_set_freq : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%freq(:) = d_freq(:)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_freq


SUBROUTINE colloid_set_m(this, d_m, stat_info)
  !---------------------------------------
  ! Set the total mass of a colloid object
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:), INTENT(IN)      :: d_m
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num


  stat_info = 0

  num = SIZE(d_m,1)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_m : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%m(:) = d_m(1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_m


SUBROUTINE colloid_set_mmi(this, d_mmi, stat_info)
  !---------------------------------------
  ! Set the mass momentum inertia of colloids
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:), INTENT(IN)      :: d_mmi
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num


  stat_info = 0

  num = SIZE(d_mmi,1)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_mmi : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%mmi(:) = d_mmi(1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_mmi


SUBROUTINE colloid_set_x(this,d_x,stat_info)
  !---------------------------------------
  ! Set the positions of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_x
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num

  stat_info = 0

  !---------------------------------------
  ! Check if the input position's dimension
  ! and num of colloids match.
  !---------------------------------------
  dim = SIZE(d_x,1)
  num = SIZE(d_x,2)

  IF( dim /= this%num_dim) THEN
     PRINT *, "colloid_set_x : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_x : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%x(:,:) = d_x(1:dim,1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_x


SUBROUTINE colloid_set_v(this,d_v,stat_info)
  !---------------------------------------
  ! Set the velocity of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_v
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num

  stat_info = 0

  !---------------------------------------
  ! Check if the input velocity's dimension
  ! and num of colloids match.
  !---------------------------------------
  dim = SIZE(d_v,1)
  num = SIZE(d_v,2)

  IF( dim /= this%num_dim) THEN
     PRINT *, "colloid_set_v : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_v : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%v(:,:) = d_v(1:dim,1:num)

9999 CONTINUE
  RETURN       

END SUBROUTINE colloid_set_v


#if __DRAG_PART
SUBROUTINE colloid_set_drag_lub(this,d_drag,stat_info)
  !---------------------------------------
  ! Set the drag/force lubrication part
  ! exerted on a colloid object
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_drag
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num


  stat_info = 0

  !---------------------------------------
  ! Check if the input force's dimension and
  ! num of colloids match.
  !---------------------------------------

  dim = SIZE(d_drag,1)
  num = SIZE(d_drag,2)

  IF( dim /= this%num_dim) THEN
     PRINT *, "colloid_set_drag_lub : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_drag_lub : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%drag_lub(:,:) = d_drag(1:dim,1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_drag_lub


SUBROUTINE colloid_set_drag_repul(this,d_drag,stat_info)
  !---------------------------------------
  ! Set the drag/force repulsive part
  ! exerted on a colloid object
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_drag
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num


  stat_info = 0

  !---------------------------------------
  ! Check if the input force's dimension and
  ! num of colloids match.
  !---------------------------------------

  dim = SIZE(d_drag,1)
  num = SIZE(d_drag,2)

  IF( dim /= this%num_dim) THEN
     PRINT *, "colloid_set_drag : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_drag : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%drag_repul(:,:) = d_drag(1:dim,1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_drag_repul
#endif


SUBROUTINE colloid_set_drag(this,d_drag,stat_info)
  !---------------------------------------
  ! Set the drag/force exerted on a colloid object
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_drag
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num


  stat_info = 0

  !---------------------------------------
  ! Check if the input force's dimension and
  ! num of colloids match.
  !---------------------------------------

  dim = SIZE(d_drag,1)
  num = SIZE(d_drag,2)

  IF( dim /= this%num_dim) THEN
     PRINT *, "colloid_set_drag : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_drag : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%drag(:,:) = d_drag(1:dim,1:num)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_drag


SUBROUTINE colloid_add_drag(this, d_drag,stat_info)
  !---------------------------------------
  ! Accumulate the drag/force on colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_drag
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num

  stat_info = 0

  dim = SIZE(d_drag,1)
  num = SIZE(d_drag,2)

  !---------------------------------------
  ! check if the input force's dimension matches
  !---------------------------------------

  IF (dim /= this%num_dim) THEN
     PRINT *, "colloid_add_drag : ", &
          "Dimension doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF
  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_add_drag : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%drag(1:dim,1:num) =  &
       this%drag(1:dim,1:num) + d_drag(1:dim,1:num)

9999 CONTINUE
  RETURN       

END SUBROUTINE colloid_add_drag


SUBROUTINE colloid_set_rotation_matrix(this,d_matrix,stat_info)
  !---------------------------------------
  ! Set rotation matrix.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:,:),INTENT(IN)   :: d_matrix
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim1,dim2, num

  stat_info = 0

  dim1 = SIZE(d_matrix,1)
  dim2 = SIZE(d_matrix,2)
  num  = SIZE(d_matrix,3)

  IF ( dim1 /= this%num_dim .OR. &
       dim2 /= this%num_dim ) THEN
     PRINT *, "colloid_set_rotation_matrix : ", &
          "Number of dimension doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  IF ( num /= this%num_colloid ) THEN
     PRINT *, "colloid_set_rotation_matrix : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%rot_matrix(:,:,:) = d_matrix(:,:,:)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_rotation_matrix


SUBROUTINE colloid_set_phi(this,d_phi,stat_info)
  !---------------------------------------
  ! Set rotated angle.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_phi
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num

  stat_info = 0

  num = SIZE(d_phi,2)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_phi : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%phi(:,:) = d_phi(:,:)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_phi


SUBROUTINE colloid_set_theta(this,d_theta,stat_info)
  !---------------------------------------
  ! Set rotated angle.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:),INTENT(IN)       :: d_theta
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num

  stat_info = 0

  num = SIZE(d_theta,2)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_theta : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%theta(:,:) = d_theta(:,:)

9999 CONTINUE

  RETURN       

END SUBROUTINE colloid_set_theta


SUBROUTINE colloid_set_omega(this,d_omega,stat_info)
  !---------------------------------------
  ! Set the rotating velocity of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_omega
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num

  stat_info = 0

  !---------------------------------------
  ! Check if the input rotating velocity's
  ! dimension is 3.
  !---------------------------------------
  dim = SIZE(d_omega,1)
  num = SIZE(d_omega,2)

  IF( dim /= 3) THEN
     PRINT *, "colloid_set_oemga : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_omega : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%omega(1:3,1:num) = d_omega(1:3,1:num)

9999 CONTINUE
  RETURN       

END SUBROUTINE colloid_set_omega


SUBROUTINE colloid_set_torque(this,d_torque,stat_info)
  !---------------------------------------
  ! Set the torque of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_torque
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: dim
  INTEGER                                 :: num

  stat_info = 0

  !---------------------------------------
  ! Check if the input torque's
  ! dimension is 3.
  !---------------------------------------
  dim = SIZE(d_torque,1)
  num = SIZE(d_torque,2)

  IF( dim /= 3) THEN
     PRINT *, "colloid_set_torque : ", "Wrong Dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  IF( num /= this%num_colloid) THEN
     PRINT *, "colloid_set_torque : ", "Wrong number of colloids !"
     stat_info = -1
     GOTO 9999
  END IF

  this%torque(1:3,1:num) = d_torque(1:3,1:num)

9999 CONTINUE
  RETURN       

END SUBROUTINE colloid_set_torque


SUBROUTINE colloid_set_num_physical_part(this,&
     d_num_physical_part,stat_info)
  !---------------------------------------
  ! Set the number of physical particles,
  ! which consitute colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  INTEGER, DIMENSION(:), INTENT(IN)       :: d_num_physical_part
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num, i


  stat_info = 0

  num = SIZE(d_num_physical_part,1)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_num_physical_part : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF

  this%num_physical_part(:) = &
       d_num_physical_part(1:num)

  this%num_physical_part_tot = 0

  DO i =1, this%num_colloid

     this%num_physical_part_tot =  &
          this%num_physical_part_tot + &
          this%num_physical_part(i)                   
  END DO

9999 CONTINUE  

  RETURN       

END SUBROUTINE colloid_set_num_physical_part


SUBROUTINE colloid_set_num_numerical_part(this,&
     d_num_numerical_part,stat_info)
  !---------------------------------------
  ! Set the number of numerical particles,
  ! which consitute colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)            :: this
  INTEGER, DIMENSION(:), INTENT(IN)       :: d_num_numerical_part
  INTEGER, INTENT(OUT)                    :: stat_info

  INTEGER                                 :: num,i

  stat_info = 0

  num = SIZE(d_num_numerical_part,1)

  IF (num /= this%num_colloid) THEN
     PRINT *, "colloid_set_num_numerical_part : ", &
          "Number of colloids doesn't match !"
     stat_info = -1
     GOTO 9999
  END IF
  this%num_numerical_part(:) = &
       d_num_numerical_part(:)

  this%num_numerical_part_tot = 0

  DO i =1, this%num_colloid

     this%num_numerical_part_tot =  &
          this%num_numerical_part_tot + &
          this%num_numerical_part(i)                   
  END DO

9999 CONTINUE        
  RETURN       

END SUBROUTINE colloid_set_num_numerical_part


SUBROUTINE colloid_set_min_phys(this,d_min_phys,stat_info)

  TYPE(Colloid), INTENT(INOUT)             :: this
  REAL(MK), DIMENSION(:)                   :: d_min_phys
  INTEGER, INTENT(OUT)                     :: stat_info

  INTEGER                                  :: dim

  stat_info = 0
  dim = SIZE(d_min_phys)
  IF (dim /= this%num_dim ) THEN
     PRINT *, "colloid_set_min_phys : ", "Wrong dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  this%min_phys(1:dim) = d_min_phys(1:dim)

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_min_phys


SUBROUTINE colloid_set_max_phys(this,d_max_phys,stat_info)

  TYPE(Colloid), INTENT(INOUT)             :: this
  REAL(MK), DIMENSION(:)                   :: d_max_phys
  INTEGER, INTENT(OUT)                     :: stat_info

  INTEGER                                  :: dim

  stat_info = 0
  dim = SIZE(d_max_phys)
  IF (dim /= this%num_dim ) THEN
     PRINT *, "colloid_set_max_phys : ", "Wrong dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  this%max_phys(1:dim) = d_max_phys(1:dim)

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_max_phys


SUBROUTINE colloid_set_min_phys_t(this,d_min_phys_t,stat_info)

  TYPE(Colloid), INTENT(INOUT)             :: this
  REAL(MK), DIMENSION(:)                   :: d_min_phys_t
  INTEGER, INTENT(OUT)                     :: stat_info

  INTEGER                                  :: dim

  stat_info = 0

  dim = SIZE(d_min_phys_t)

  IF (dim /= this%num_dim ) THEN
     PRINT *, "colloid_set_min_phys_t : ", "Wrong dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  this%min_phys_t(1:dim) = d_min_phys_t(1:dim)

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_min_phys_t


SUBROUTINE colloid_set_max_phys_t(this,d_max_phys_t,stat_info)

  TYPE(Colloid), INTENT(INOUT)             :: this
  REAL(MK), DIMENSION(:)                   :: d_max_phys_t
  INTEGER, INTENT(OUT)                     :: stat_info

  INTEGER                                  :: dim

  stat_info = 0

  dim = SIZE(d_max_phys_t)

  IF (dim /= this%num_dim ) THEN
     PRINT *, "colloid_set_max_phys_t : ", "Wrong dimension !"
     stat_info = -1
     GOTO 9999
  END IF

  this%max_phys_t(1:dim) = d_max_phys_t(1:dim)

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_max_phys_t


SUBROUTINE colloid_set_cut_off(this,d_cut_off,stat_info)

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_cut_off
  INTEGER, INTENT(OUT)            :: stat_info

  REAL(MK)                        :: din
  INTEGER                         :: stat_info_sub

  stat_info     = 0
  stat_info_sub = 0

  IF (d_cut_off <= 0.0_MK) THEN
     PRINT *, "colloid_set_cut_off : ", &
          "cut off should be non-negative !"
     stat_info = -1
     GOTO 9999
  END IF

  this%cut_off = d_cut_off

  din = mcf_colloid_inner_layer_coeff*d_cut_off

  IF ( this%rho_type == mcf_colloid_rho_type_dynamic ) THEN

     din = din * 2.0_MK

  END IF

  CALL colloid_set_din(this,din, stat_info_sub)

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_cut_off


SUBROUTINE colloid_set_dout(this,d_dout,stat_info)
  !---------------------------------------
  ! Set the minimal distance of a fluid
  ! particle from the surface of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_dout
  INTEGER, INTENT(OUT)            :: stat_info


  stat_info = 0

  IF (d_dout < 0.0_MK) THEN
     PRINT *, "colloid_set_dout : ", &
          "dout should be positive !"
     stat_info = -1
     GOTO 9999
  END IF

  this%dout = d_dout

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_dout


SUBROUTINE colloid_set_din(this,d_din, stat_info)
  !---------------------------------------
  ! Set the maximum distance of a boundary
  ! particle from the surface of colloids.
  !---------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_din
  INTEGER, INTENT(OUT)            :: stat_info


  stat_info = 0

  IF ( d_din < 0.0_MK ) THEN
     PRINT *, "colloid_set_din : ", &
          "din should be positive !"
     stat_info = -1
     GOTO 9999
  END IF

  this%din = d_din

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_din


SUBROUTINE colloid_set_eta(this,d_eta,stat_info)
  !-----------------------------------------
  ! Set the fluid dynamics viscosity.
  !-----------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  REAL(MK), INTENT(IN)            :: d_eta
  INTEGER, INTENT(OUT)            :: stat_info

  stat_info = 0

  IF ( d_eta < 0.0_MK ) THEN
     PRINT *, "colloid_set_eta : ", &
          "Wrong eta !"
     stat_info = -1
  END IF

  this%eta = d_eta

  RETURN

END SUBROUTINE colloid_set_eta


SUBROUTINE colloid_set_bcdef(this,d_bcdef,stat_info)

  TYPE(Colloid), INTENT(INOUT)          :: this
  INTEGER, DIMENSION(:)                 :: d_bcdef
  INTEGER, INTENT(OUT)                  :: stat_info


  INTEGER                               :: bcdef_dim

  stat_info = 0

  bcdef_dim = SIZE(d_bcdef)

  IF( bcdef_dim /= 2*this%num_dim ) THEN
     PRINT *, "colloid_set_bcdef : ", &
          "Wrong dimension !"
     stat_info = -1
     GOTO 9999

  END IF

  this%bcdef(1:bcdef_dim) = d_bcdef(1:bcdef_dim)


9999 CONTINUE

  RETURN

END SUBROUTINE  colloid_set_bcdef


SUBROUTINE colloid_set_boundary(this,d_boundary,stat_info)

  TYPE(Colloid), INTENT(INOUT)    :: this
  TYPE(Boundary), TARGET          :: d_boundary
  INTEGER, INTENT(OUT)            :: stat_info

  INTEGER                         :: stat_info_sub
  INTEGER                         :: num_dim

  stat_info     = 0
  stat_info_sub = 0

  num_dim = &
       boundary_get_num_dim(d_boundary,stat_info_sub)

  IF(num_dim /= this%num_dim) THEN
     PRINT *, "colloid_set_boundary : ",&
          "Boundarys dimension doesn't match with colloid !"
     stat_info = -1
     GOTO 9999
  END IF

  this%boundary => d_boundary

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_boundary


#if 0
SUBROUTINE colloid_set_image(this,stat_info)
  !----------------------------------------------------
  ! Remark      :
  !               In case of periodic or Lees-Edwards
  !               boundaries, the images of colloid's
  !               center has to be taken into account
  !               to decide if a potential particle 
  !               is inside the geometry of a colloid.
  !
  !               For 2D, maximum 3**2=9=8 images;
  !               For 3D, maximum 3**3=27=26 images
  !               (including itself).
  !----------------------------------------------------

  TYPE(Colloid), INTENT(INOUT)          :: this
  INTEGER, INTENT(OUT)                  :: stat_info


  INTEGER                               :: num_dim
  INTEGER                               :: num
  INTEGER, DIMENSION(3)                 :: istart
  INTEGER, DIMENSION(3)                 :: iend
  INTEGER                               :: i,j,k       
  REAL(MK), DIMENSION(3,27)             :: image



  stat_info = 0


  num_dim = this%num_dim

  istart(:) = 0
  iend(:)   = 0

  DO i = 1, num_dim

     !-------------------------------------------------
     ! If one side is periodic, 
     ! the opposite must be periodic also.
     ! If one side is Lees-Edwards, 
     ! the opposite must be Lees-Edwards also.
     !-------------------------------------------------

     IF ( this%bcdef(2*i-1) == mcf_bcdef_periodic .OR. &
          this%bcdef(2*i-1) == mcf_bcdef_LE ) THEN

        istart(i) = -1
        iend(i)   = 1

     END IF

  END DO

  num = 1
  image(1:num_dim,1) = 0.0_MK

  DO k = istart(3), iend(3)

     DO j = istart(2), iend(2)

        DO i = istart(1), iend(1)

           !-------------------------------------------
           ! Box(cell) itself doesn't count,
           ! we exclude.
           !-------------------------------------------

           IF ( i == 0 .AND. &
                j == 0 .AND. &
                k == 0 ) THEN

              CYCLE

           END IF

           num = num + 1
           image(1, num) = i * &
                ( this%max_phys(1) -this%min_phys(1) )
           image(2, num) = j * &
                ( this%max_phys(2) -this%min_phys(2) )

           IF ( num_dim == 3 ) THEN
              image(3, num) = k * &
                   ( this%max_phys(3) -this%min_phys(3) )
           END IF

        END DO ! i

     END DO ! j

  END DO !  k

  this%num_image = num

  IF (ASSOCIATED(this%image)) THEN
     DEALLOCATE(this%image)
  END IF

  ALLOCATE(this%image(num_dim,num))

  this%image(1:num_dim,1:num) = &
       image(1:num_dim,1:num)


9999 CONTINUE

  RETURN

END SUBROUTINE  colloid_set_image
#endif


SUBROUTINE colloid_set_arbitrary(this, d_arbitrary_num, in_file, stat_info)
  !-----------------------------------------
  ! Set up the arbitrary shape colloids.
  !-----------------------------------------

  TYPE(Colloid), INTENT(INOUT)    :: this
  INTEGER, INTENT(IN)             :: d_arbitrary_num
  CHARACTER(LEN=MAX_CHAR)         :: in_file
  INTEGER, INTENT(OUT)            :: stat_info  

  INTEGER                         :: ilenread
  LOGICAL                         :: lExist
  INTEGER                         :: num_dim, num_vperf, max_degree
  INTEGER                         :: i, j


  stat_info = 0

  this%coll_arbitrary_num  = d_arbitrary_num

  IF (d_arbitrary_num > 0) THEN


     ! check the file
     ilenread = LEN_TRIM(in_file)

     !----------------------------------------------------
     ! Check if a file name is given.
     !----------------------------------------------------   

     IF ( ilenread < 1 ) THEN
        PRINT *,'coll_arbitrary_file : ',&
             'No file name is given !'
        stat_info = -1
        GOTO 9999

     ELSE

        !----------------------------------------------------
        ! Check if the file exists. 
        !----------------------------------------------------

        INQUIRE(FILE=in_file,EXIST=lExist)

        IF (.NOT.lExist) THEN 
           PRINT *,'NO such file : ',&
                in_file
           stat_info = -1
           GOTO 9999
        END IF
     END IF


     NULLIFY(this%coll_v) 
     NULLIFY(this%coll_v_flist) 
     NULLIFY(this%coll_f_vlist)  
     NULLIFY(this%coll_n) 
     NULLIFY(this%coll_v_num) 
     NULLIFY(this%coll_f_num) 
     NULLIFY(this%coll_inout)
     NULLIFY(this%coll_sid)


     num_dim = this%num_dim 
     num_vperf = num_dim

     ALLOCATE(this%coll_v_num(d_arbitrary_num))  
     ALLOCATE(this%coll_f_num(d_arbitrary_num)) 
     ALLOCATE(this%coll_inout(d_arbitrary_num))

     ! read input file
     OPEN(10,  FILE=in_file)

     DO j = 1, d_arbitrary_num

        READ(10, *)  this%coll_v_num(j), max_degree
        IF ( j == 1) THEN
           ALLOCATE(this%coll_v(d_arbitrary_num, num_dim, MAXVAL(this%coll_v_num)))
           ALLOCATE(this%coll_v_flist(d_arbitrary_num, max_degree, &
                MAXVAL(this%coll_v_num)))
        END IF

        DO i = 1, this%coll_v_num(j)
           READ(10, *)  this%coll_v(j, 1:num_dim, i)

!!$           ! §§§
!!$           ! quick fix: avoid those at outlet to be moved to inlet 
!!$           ! due to periodic BC (I assume periodic in x direction).
!!$           IF (this%coll_v(j, 1, i) == this%max_phys(1)) THEN
!!$              this%coll_v(j, 1, i) = this%coll_v(j, 1, i) - &
!!$                   this%cut_off / 1000.0_MK
!!$           END IF
!!$           ! and vice versa
!!$           IF (this%coll_v(j, 1, i) == this%min_phys(1)) THEN
!!$              this%coll_v(j, 1, i) = this%coll_v(j, 1, i) + &
!!$                   this%cut_off / 1000.0_MK
!!$           END IF

        END DO

        DO i = 1, this%coll_v_num(j)
           READ(10, *)  this%coll_v_flist(j, 1:max_degree, i)
        END DO

        READ(10, *)  this%coll_f_num(j)
        IF ( j == 1) THEN 
           ALLOCATE(this%coll_f_vlist(d_arbitrary_num, num_vperf, &
                MAXVAL(this%coll_f_num)))
           ALLOCATE(this%coll_n(d_arbitrary_num, num_dim, MAXVAL(this%coll_f_num)))
        END IF

        DO i = 1, this%coll_f_num(j)
           READ(10, *)  this%coll_f_vlist(j, 1:num_vperf, i)
        END DO

        ! surface normal for each facet
        DO i = 1, this%coll_f_num(j)
           READ(10, *)  this%coll_n(j, 1:num_dim, i)
        END DO

        ! sorting data
        READ(10, *)  this%coll_sdx
        READ(10, *)  this%coll_snum, this%coll_sco
        ALLOCATE(this%coll_sid(this%coll_snum))
        DO i = 1, this%coll_snum
           READ(10, *)  this%coll_sid(i)
        END DO


     END DO

     CLOSE(10)

  END IF

9999 CONTINUE

  RETURN

END SUBROUTINE colloid_set_arbitrary



SUBROUTINE colloid_arbitrary_distance(this, sort, col_index, p_x, out_distance, out_normal, fi, stat_info)
  !-----------------------------------------
  ! Set arbitrary colloid distance
  !-----------------------------------------

  TYPE(Colloid)         , INTENT(IN)    :: this
  LOGICAL               , INTENT(IN)    :: sort
  INTEGER               , INTENT(INOUT) :: col_index
  REAL(MK), DIMENSION(:), INTENT(IN)    :: p_x
  REAL(MK)              , INTENT(OUT)   :: out_distance
  REAL(MK), DIMENSION(:), INTENT(OUT)   :: out_normal
  INTEGER               , INTENT(OUT)   :: fi
  INTEGER               , INTENT(OUT)   :: stat_info 

  INTEGER                               :: num_dim, num_vperf, max_degree
  INTEGER                               :: px_dim, N_interface, min_dist
  INTEGER                               :: i, j, closest_i, ii, jj, vi
  INTEGER                               :: istart, iend, sloc
  REAL(MK)                              :: closest_d, d, distance 
  REAL(MK)                              :: bary_u, bary_v
  REAL(MK), ALLOCATABLE, DIMENSION(:,:) :: x_point  
  REAL(MK), ALLOCATABLE, DIMENSION(:)   :: projection_point, lambda
  REAL(MK), ALLOCATABLE, DIMENSION(:)   :: mid_point, d_temp
  REAL(MK), allocatable, dimension(:)   :: e, normal_distance
  REAL(MK), allocatable, dimension(:)   :: proj_test
  INTEGER                               :: stat_info_sub


  stat_info = 0
  stat_info_sub = 0

  num_dim = this%num_dim  
  num_vperf = num_dim
  px_dim = SIZE(p_x)

  IF( px_dim /= num_dim ) THEN
     PRINT *, "colloid_arbitrary_distance : ", &
          "Wrong dimension !"
     stat_info = -1
     GOTO 9999

  END IF


  max_degree =  SIZE(this%coll_v_flist, 2)

  allocate(mid_point(num_dim), d_temp(num_dim))  
  allocate(x_point(num_dim, num_vperf)) 
  allocate(projection_point(num_dim))
  allocate(lambda(num_vperf))
  allocate(e(max_degree), normal_distance(max_degree))
  allocate(proj_test(max_degree))


  IF (sort) THEN
     ! Used sort data to limit search based on location
     ! (this assumes only one colloid. in fact there is no
     !  need to define than one. For now I assume only one 
     !  arbitary colloid exists. coll_arbitrary_num will be
     !  removed in future.)
     sloc = FLOOR((p_x(this%coll_sco) - this%min_phys(this%coll_sco)) / this%coll_sdx) + 1
     sloc = MAX(sloc, 1)
     istart = this%coll_sid(sloc) + 1
     iend   = this%coll_sid(sloc+1)
  ELSE
     istart = 1
     iend   = MAXVAL(this%coll_v_num)
  END IF
  
  ! if we know which colloid is closer (or of interest)..
  IF (col_index == 0) THEN 

     closest_d = 100000.0_MK * this%cut_off    
     DO ii = 1, this%coll_arbitrary_num

        DO i = istart, iend

           d_temp =  this%coll_v(ii, 1:num_dim, i) - p_x(1:num_dim)
           d = DOT_PRODUCT(d_temp, d_temp)
           d = d ** 0.50_MK
           IF (d < closest_d) THEN
              closest_d = d
              closest_i = i
              col_index = ii
           END IF
        END DO

     END DO

  END IF


  ! choose the best facet
  e = 1000.0_MK 
  normal_distance = 1000.0_MK
  out_distance = 1000.0_MK
  proj_test = 1000.0_MK
  
  DO ii = 1, max_degree
     fi = this%coll_v_flist(col_index, ii, closest_i)

     IF (fi > 0 ) THEN
        ! for each facet...
        d_temp = 0.0d0
        DO jj = 1, num_vperf
           vi = this%coll_f_vlist(col_index, jj, fi)
           d_temp = d_temp + this%coll_v(col_index, 1:num_dim, vi)
        END DO
        mid_point = d_temp / DBLE(num_vperf)

        ! vector from the mid_point to grid point
        d_temp = p_x(1:num_dim) - mid_point(1:num_dim)

        ! distance to mid_point
        e(ii) = DOT_PRODUCT(d_temp, d_temp) ** 0.50_MK
        ! dot product gives the normal distance
        normal_distance(ii) = DOT_PRODUCT(d_temp, this%coll_n(col_index, 1:num_dim, fi))

        ! We still need to make sure that the normal 
        ! distance is the closest distance (in case
        ! of convex curves/surfaces it always is, but
        ! for concave curves/surfaces it might not be.)
        ! The idea is to check if the projection point
        ! belongs to the curve/surface.

        ! calculate projection
        ! (the reason for negative sign is the direction 
        ! of normal vector (which is towards the interior), 
        ! notice that distance is positive for interior and 
        ! negative for points outside fluid domain.
        projection_point = p_x(1:num_dim) - normal_distance(ii) * this%coll_n(col_index, 1:num_dim, fi)


        DO i = 1, num_vperf
           vi = this%coll_f_vlist(col_index, i, fi)
           x_point(1:num_dim, i) = this%coll_v(col_index, 1:num_dim, vi)
        END DO

        ! find barycentric coordinates 
        CALL tool_barycentric_coordinate(this%tool, &
             x_point, projection_point, lambda, stat_info_sub) 


        ! determine whether projection point in inside or not
        ! if inside it is 1.0, if outside larger (depending on how much).
        proj_test(ii) = SUM(ABS(lambda(1:num_vperf)))

     END IF

  END DO


  IF (ABS(MINVAL(proj_test) - 1.0_MK) <= mcf_machine_zero) THEN

     ! one or more of projection points are inside
     DO ii = 1, max_degree
        IF ((proj_test(ii) - 1.0_MK) <= mcf_machine_zero) THEN
           IF (ABS(normal_distance(ii)) < ABS(out_distance)) THEN

              fi = this%coll_v_flist(col_index, ii, closest_i)
              out_distance = normal_distance(ii)
              out_normal   = this%coll_n(col_index, 1:num_dim, fi)
              
              
           END IF
        END IF
     END DO

  ELSE  

     ! none of the projection points are inside: we choose
     ! the facet with smallest proj_test, i.e. the facet with 
     ! smallest sum of absolute value of barycentric coordinates.

     min_dist = MINLOC(proj_test, 1)
     fi = this%coll_v_flist(col_index, min_dist, closest_i)
     out_distance = normal_distance(min_dist) ! normal distance
     out_normal   = this%coll_n(col_index, 1:num_dim, fi)

  END IF


9999 CONTINUE

     RETURN

   END SUBROUTINE colloid_arbitrary_distance
