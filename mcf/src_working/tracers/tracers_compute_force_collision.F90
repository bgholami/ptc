SUBROUTINE tracers_compute_force_collision(this, dt, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_compute_force_collision
  !----------------------------------------------------
  !
  ! Purpose     : Applying collision force,
  !               which means displacements will be
  !               superimposed on the current position.
  !                  
  !
  ! Routines    :  
  !
  ! References  :  Longest et al.
  !                2004 Computers and Fluids
  !                Fogelson
  !                1992 J. of Computational Physics
  !
  !
  ! Remarks     :  
  !
  ! Revisions   : V0.1 07.07 2011, original
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

  TYPE(Particles), INTENT(INOUT)      :: this        
  REAL(MK)       , INTENT(IN)         :: dt
  INTEGER        , INTENT(OUT)	      :: stat_info  

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !----------------------------------------------------

  INTEGER                             :: num_dim        

  !----------------------------------------------------
  ! Number of all tracers
  !----------------------------------------------------

  INTEGER                             :: num_part_real  

  !----------------------------------------------------
  ! Local variables start here :
  ! d_p       : diffusion coefficient
  ! gauss_rnd : gaussian random variable 
  !----------------------------------------------------
  
  REAL(MK), DIMENSION(:), ALLOCATABLE :: nv, t1, t2, delta_x, dxt 
  INTEGER                             :: i, ip, it
  REAL(MK)                            :: d_p, a_p
  REAL(MK)                            :: gauss_rnd
  REAL(MK)                            :: fn, ft, z, h
  INTEGER                             :: stat_info_sub

  !----------------------------------------------------


  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub) 

  !----------------------------------------------------
  ! Number of real tracers
  !----------------------------------------------------
  num_part_real = &
       tracers_get_num_part_real(this,stat_info_sub) 

  !----------------------------------------------------
  ! Allocate memory 
  !----------------------------------------------------
  ALLOCATE(nv(num_dim), t1(num_dim))
  IF (num_dim == 3) THEN
     ALLOCATE(t2(num_dim))
  END IF
  ALLOCATE(delta_x(num_dim), dxt(num_dim))


  !----------------------------------------------------
  ! tracers' parameters  
  !----------------------------------------------------
  a_p = this%radius  
  d_p = this%dp


  ! loop over tracers
  DO ip = 1, num_part_real

     ! only if a normal tracer
     IF (this%id(1,ip) > 0) THEN


         ! normal vec (nv)
        nv = this%n_vec(1:num_dim, ip)

        ! one/two vectors in tangent plane (t1, t2)
        ! (we assume nv is not 0)
        it = 1
        DO WHILE ((it <= num_dim) .AND. (nv(it) == 0.0_MK))
           it = it + 1
        END DO
        t1 = 1.0_MK
        t1(it) = (-SUM(nv) + nv(it)) / nv(it)
        t1 = t1 / (DOT_PRODUCT(t1, t1) ** 0.50_MK)
        
        IF (num_dim == 3) THEN
           ! use cross product of nv and t1
           t2(1) = nv(2) * t1(3) - nv(3) * t1(2) 
           t2(2) = nv(3) * t1(1) - nv(1) * t1(3) 
           t2(3) = nv(1) * t1(2) - nv(2) * t1(1)
        END IF

        
        ! evaluate diffusion is hindrance factors
        fn = 0.0_MK
        ft = 0.0_MK
        IF (this%h_p(ip) > 0.0_MK) THEN

           ! to avoid confusion (following notation of kazoe2011measurements)
           z = this%h_p(ip)
           h = z - a_p
           
           fn = (6.0_MK * h**2.0_MK + 2.0_MK * a_p * h) / &
                (6.0_MK * h**2.0_MK + 9.0_MK * a_p * h + 2.0_MK * a_p**2.0_MK)
           
           ft = 1.0_MK - 9.0_MK/16.0_MK * (a_p/z) + 1.0_MK/8.0_MK * (a_p/z)**3.0_MK - &
                45.0_MK/256.0_MK * (a_p/z)**4.0_MK - 1.0_MK/16.0_MK * (a_p/z)**5.0_MK

        END IF


        ! Once for each dimension (independent)
        DO i = 1, num_dim
           
           gauss_rnd    = random_random_Gaussian2(this%random, stat_info_sub) 
           delta_x(i)   = gauss_rnd * (2.0d0*dt*d_p)**0.50d0

        END DO
        
        ! hinder diffusion (first index corresponds to normal direction)
        delta_x(1) = delta_x(1) * fn
        delta_x(2:num_dim) = delta_x(2:num_dim) * ft


        ! present displacement vector in xyz coordinates
        dxt(1:num_dim) = &
             delta_x(1) * nv(1:num_dim) + &
             delta_x(2) * t1(1:num_dim)
        IF (num_dim == 3) THEN  
           dxt(1:num_dim) = dxt(1:num_dim) + &
                delta_x(3) * t2(1:num_dim)
        END IF

        ! superimpose the displacement on current position
        this%x(1:num_dim,ip) = this%x(1:num_dim,ip) + dxt(1:num_dim)


     END IF

  END DO

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE        


  RETURN

END SUBROUTINE tracers_compute_force_collision

