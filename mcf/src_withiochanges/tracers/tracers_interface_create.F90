SUBROUTINE tracers_interface_create(this, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_interface_create
  !----------------------------------------------------
  !
  ! Purpose     : creates interface using wall/colloid 
  !               particles
  !
  ! 
  ! Routines    : 
  !               
  !               
  !
  !
  ! Remarks     : 
  !
  !
  ! References  :
  !
  ! Revisions   : V0.2  04.02 2013, for arbitrary geometry
  !               V0.1  08.11 2011, original version.
  !
  !----------------------------------------------------
  ! Author      : Babak Gholami
  ! Contact     : babak.gholami@aer.mw.tum.de
  !
  ! Dr. Marco Ellero's Emmy Noether Group,
  ! Prof. Dr. N. Adams' Chair of Aerodynamics,
  ! Faculty of Mechanical Engineering,
  ! Technische Universitaet Muenchen, Germany.
  !----------------------------------------------------

  !----------------------------------------------------
  ! Modules :
  !----------------------------------------------------


  !----------------------------------------------------
  ! Arguments :
  !----------------------------------------------------

  TYPE(Particles), INTENT(INOUT)          :: this
  INTEGER, INTENT(OUT)	                  :: stat_info  


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim    : number of dimension. 
  !----------------------------------------------------

  INTEGER                                 :: num_dim        
  REAL(MK), DIMENSION(:), POINTER         :: dx    

  !----------------------------------------------------
  ! Technique parameters :
  !
  ! min_sub    : dimension of the subdomain (minimum) 
  ! max_sub    : dimension of the subdomain (maximum)
  !----------------------------------------------------

  REAL(MK), DIMENSION(:), POINTER         :: min_sub, max_sub

  !----------------------------------------------------
  ! Colloid parameters :
  !
  ! coll_x     : center of the colloid
  ! coll_rad   : dimensions of the colloid
  !
  !----------------------------------------------------

  TYPE(Colloid),POINTER                   :: colloids
  INTEGER                                 :: num_colloid  
  INTEGER                                 :: coll_arbitrary_num 
  INTEGER , DIMENSION(:)    , POINTER     :: coll_v_num
  REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_v 
  INTEGER , DIMENSION(:,:,:), POINTER     :: coll_v_flist
  REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_n


  !----------------------------------------------------
  ! Local variables start here :
  !----------------------------------------------------

  INTEGER                                 :: max_degree, f_counter, fi, local_check
  INTEGER                                 :: i, j, k, stat_info_sub
  REAL(MK)                                :: d_near_wall, a_p, dum
  REAL(MK), DIMENSION(:), ALLOCATABLE     :: mean_normal
  LOGICAL                                 :: count_flag

  CHARACTER(len=100)                      :: in_file  = 'interface_velocity.data'
  REAL(MK)                                :: dum1, dum2

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0

  NULLIFY(colloids)
  NULLIFY(coll_v_num) 
  NULLIFY(coll_v) 
  NULLIFY(coll_v_flist) 
  NULLIFY(coll_n) 

  NULLIFY(dx)
  NULLIFY(min_sub) 
  NULLIFY(max_sub)


  !----------------------------------------------------
  ! Get physics variables,
  ! species, dimension and which lattice.
  !----------------------------------------------------

  num_dim  = &
       physics_get_num_dim(this%phys,stat_info_sub)  

  !----------------------------------------------------
  ! get dimensions of this subdomain
  !---------------------------------------------------- 

  CALL technique_get_min_sub(this%tech, min_sub, stat_info_sub) 
  CALL technique_get_max_sub(this%tech, max_sub, stat_info_sub)

  !----------------------------------------------------
  ! Get colloid
  !----------------------------------------------------

  num_colloid = &
       physics_get_num_colloid(this%phys,stat_info_sub)

  coll_arbitrary_num = 0 

  IF ( num_colloid > 0 ) THEN

     CALL physics_get_colloid(this%phys, &
          colloids,stat_info_sub) 

     coll_arbitrary_num = colloid_get_arbitrary_num(colloids, stat_info_sub)

  END IF


  IF (coll_arbitrary_num > 0) THEN
     !----------------------------------------------------
     ! initialize colloid variables
     !----------------------------------------------------  

     CALL colloid_get_v_num(colloids, coll_v_num, stat_info_sub)
     CALL colloid_get_coll_v(colloids, coll_v, stat_info_sub)  
     CALL colloid_get_coll_v_flist(colloids, coll_v_flist, stat_info_sub)  
     CALL colloid_get_coll_n(colloids, coll_n, stat_info_sub)


     !----------------------------------------------------
     ! initialize tracer variables
     !----------------------------------------------------  

     ! tracers' radius
     a_p = this%radius

     ! width of near-wall region     
     d_near_wall = this%near_wall

     CALL physics_get_dx(this%phys,dx,stat_info_sub)

     !----------------------------------------------------
     ! allocate memory for interface
     !----------------------------------------------------  

     ALLOCATE(this%x_interface(coll_arbitrary_num, num_dim, MAXVAL(coll_v_num)))
     ALLOCATE(this%vo_interface(coll_arbitrary_num, num_dim, MAXVAL(coll_v_num)))
     ALLOCATE(this%vc_interface(coll_arbitrary_num, num_dim, MAXVAL(coll_v_num)))
     ALLOCATE(this%vn_interface(coll_arbitrary_num, num_dim, MAXVAL(coll_v_num)))
     ALLOCATE(this%local_interface(coll_arbitrary_num, MAXVAL(coll_v_num)))
     ALLOCATE(this%insertion_list(coll_arbitrary_num, MAXVAL(coll_v_num)))
     ALLOCATE(this%num_interface(coll_arbitrary_num))

     IF( stat_info_sub /= 0 ) THEN
        PRINT *, &
             "tracers_create_interface : ", &
             "Allocating memory for variables has problem !"
        stat_info = -1
        GOTO 9999
     END IF

     this%x_interface = 0.0_MK
     this%vo_interface= 0.0_MK
     this%vc_interface= 0.0_MK 
     this%vn_interface= 0.0_MK
     this%local_interface = .FALSE.
     this%insertion_list = 0
     this%num_interface = 0 

     ! read input file
!     OPEN(10,  FILE=in_file)

     !----------------------------------------------------
     ! calculate position of interface points
     !--------------------------------------------------- 

     ALLOCATE(mean_normal(num_dim))

     DO j = 1, coll_arbitrary_num

        max_degree =  SIZE(coll_v_flist, 2)

        DO i = 1, coll_v_num(j)

           mean_normal(1:num_dim) = 0.0_MK
           f_counter = 0
           count_flag = .TRUE.

           DO WHILE ((count_flag) .AND. (f_counter < max_degree))

              fi = coll_v_flist(j, f_counter + 1, i)  
              IF (fi > 0 ) THEN

                 mean_normal(1:num_dim) = mean_normal(1:num_dim) + &
                      coll_n(j, 1:num_dim, fi) 

                 f_counter = f_counter + 1

              ELSE

                 count_flag = .FALSE.

              END IF

           END DO

           mean_normal(1:num_dim) = mean_normal(1:num_dim) / dble(f_counter)

           ! normalize the final normal vector

           dum = (DOT_PRODUCT(mean_normal, mean_normal))**0.50_MK

           ! normalized to unity (interface point exactly d_near_wall away from
           ! the corresponding surface point.
           mean_normal = mean_normal / dum 

           ! this enforces the interface to look exactly like surface (at the 
           ! expense of some interface points not exactly d_near_wall from 
           ! the surface. So, in case of sharpe edges the interface still looks
           ! like an offset of surface. Since the actual position of interface
           ! points is used only for velocity interpolation, as long as the
           ! interpolation scheme handles data out of the defined range 
           ! (extrapolation) these two approaches won't make a big difference. 
           ! Otherwise, i.e. no extrapolation, this second approach must be used.
           mean_normal = mean_normal / dum

           this%x_interface(j, 1:num_dim, i) = coll_v(j, 1:num_dim, i) + &
                mean_normal(1:num_dim) * d_near_wall

           ! read interface velocity
!           READ(10, *) dum1, dum2, this%vn_interface(j, 1:num_dim, i)

           ! check if it is local
           local_check = 0
           DO k = 1, num_dim
              IF ((this%x_interface(j, k, i) >= min_sub(k)) .AND. &
                   (this%x_interface(j, k, i) <= max_sub(k))) THEN

                 local_check = local_check + 1

              END IF
           END DO
           IF ( local_check == num_dim ) THEN
              this%local_interface(j, i) = .TRUE.
           END IF


        END DO

     END DO



     ! set the number of interface points
     this%num_interface = coll_v_num  


     ! set number of deposited tracers  
     IF (ASSOCIATED(this%num_tracer_deposited)) THEN
        DEALLOCATE(this%num_tracer_deposited,STAT=stat_info_sub)
     END IF
     
     ALLOCATE(this%num_tracer_deposited(coll_arbitrary_num, SIZE(coll_n, 3)), &
          STAT=stat_info_sub)      
     
     this%num_tracer_deposited = 0
 

     this%vc_interface = this%vn_interface
     this%vo_interface = this%vn_interface

!     CLOSE(10)

  END IF


9999 CONTINUE 


  IF(ASSOCIATED(dx)) THEN
     DEALLOCATE(dx)
  END IF

  IF(ASSOCIATED(min_sub)) THEN
     DEALLOCATE(min_sub)
  END IF

  IF(ASSOCIATED(max_sub)) THEN
     DEALLOCATE(max_sub)
  END IF

  IF(ASSOCIATED(coll_v_num)) THEN
     DEALLOCATE(coll_v_num)
  END IF

  IF(ASSOCIATED(coll_v)) THEN
     DEALLOCATE(coll_v)
  END IF

  IF(ASSOCIATED(coll_v_flist)) THEN
     DEALLOCATE(coll_v_flist)
  END IF

  IF(ASSOCIATED(coll_n)) THEN
     DEALLOCATE(coll_n)
  END IF


  RETURN

END SUBROUTINE tracers_interface_create




