SUBROUTINE tracers_compute_wall_distance(this, func, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_compute_wall_distance
  !----------------------------------------------------
  !
  ! Purpose     : Computes wall distance for all of 
  !               the tracers
  !
  !
  ! Routines    : 
  !
  !
  ! Remarks     :  This is completely scenario-specific.
  !                The routine is implemented to calculate 
  !                tracers' wall distances and works only for 
  !                the backward-facing step geometry.
  !                       
  !
  ! Revisions   : V0.1 12.05 2011, original
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
  ! Arguments :
  !----------------------------------------------------

  TYPE(Particles)         , INTENT(INOUT) :: this 
  ! func:  0 particle, 
  !        1 tracer with id modification for those that leave, 
  !        2 tracer with id modification for those that leave and those that deposit  
  INTEGER ,                 INTENT(IN)    :: func
  INTEGER ,                 INTENT(OUT)	  :: stat_info

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  !
  !----------------------------------------------------

  INTEGER                         :: num_dim        
  REAL(MK), DIMENSION(:), POINTER :: dmin_phys
  REAL(MK), DIMENSION(:), POINTER :: dmax_phys 
  REAL(MK), DIMENSION(:), POINTER :: dx          


  !----------------------------------------------------
  ! Colloid parameters :
  !
  !----------------------------------------------------

  TYPE(Colloid),POINTER               :: colloids
  INTEGER                             :: num_colloid, coll_arbitrary_num
  INTEGER , DIMENSION(:,:,:), POINTER :: coll_f_vlist
  !----------------------------------------------------
  ! local variables
  !----------------------------------------------------        

  REAL(MK)                            :: a_p, d_near_wall
  REAL(MK)                            :: d_far, distance, d_wall, r
  INTEGER             	              :: col_index, ip, it, fi, i
  INTEGER                             :: in, im, nbond
  REAL(MK)                            :: K0, baa, bab, bac, bad
  REAL(MK)                            :: Pa, kB, T, f, Bk, t1, t2
  REAL(MK)                            :: vt, vn, S, FS, eta
  INTEGER                             :: vperf
  REAL(MK)                            :: uni_rnd  
  INTEGER                             :: stat_info_sub
  REAL(MK), DIMENSION(:), ALLOCATABLE :: v_point, v_normal, ft
  INTEGER, DIMENSION(:), ALLOCATABLE  :: i_index

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0        

  NULLIFY(dmin_phys)
  NULLIFY(dmax_phys) 
  NULLIFY(dx)
  NULLIFY(colloids)  
  NULLIFY(coll_f_vlist) 


  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from an object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub)
  CALL physics_get_dx(this%phys,dx,stat_info_sub)
  CALL physics_get_min_phys(this%phys,dmin_phys,stat_info_sub)
  CALL physics_get_max_phys(this%phys,dmax_phys,stat_info_sub) 


  !----------------------------------------------------
  ! Get colloid
  !----------------------------------------------------

  num_colloid = &
       physics_get_num_colloid(this%phys,stat_info_sub)

  IF ( num_colloid > 0 ) THEN

     CALL physics_get_colloid(this%phys, &
          colloids,stat_info_sub)

     coll_arbitrary_num = colloid_get_arbitrary_num(colloids, stat_info_sub)

  END IF

  IF (coll_arbitrary_num > 0) THEN

     CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub)  

  END IF

  !----------------------------------------------------
  ! Tracers' parameters
  !----------------------------------------------------

  ! tracers' radius
  a_p = this%radius

  ! width of near-wall region
  d_near_wall = this%near_wall  

  !----------------------------------------------------
  ! adhesion parameters
  !----------------------------------------------------
  nbond = this%nbond
  K0    = this%K0
  baa   = this%baa  
  bab   = this%bab  
  bac   = this%bac  
  bad   = this%bad

  kB    = 1.38064880e-23 * 1.0e7 ! to convert to (cm2 * g / s2 / K)
  T     = 298 ! K

  !----------------------------------------------------
  ! prepare the wall distance array
  !----------------------------------------------------

  ! A measure of 'far'
  d_far  = 100_MK * MAXVAL(ABS(dmax_phys-dmin_phys))

  ! if we check the tracer/particle to be fluid, 
  ! we need to initialize h_p like this:
  this%h_p = d_far
  
  ! allocate memory
  ALLOCATE(v_normal(num_dim))
  ALLOCATE(ft(num_dim)) 
  vperf = num_dim 
  ALLOCATE(i_index(vperf))   
  ALLOCATE(v_point(num_dim))  

  
  ! tracer aspect ratio coefficient
  FS = 1.66810_MK ! i.e sphere
  
  ! viscosity
  eta = physics_get_eta(this%phys,stat_info_sub) 

  !----------------------------------------------------
  ! Calculation of wall distance
  !----------------------------------------------------    

  DO ip = 1, this%num_part_real 


     ! if fluid tracer/particle 
     !(§§§ see if it is possible to get rid of this.
     ! in any case, change tracers_near_wall_velocity 
     ! accordingly.)
     ! and if interior tracer/particles (id(3) = 0)
     IF ((this%id(2,ip) == 0) .AND. (this%id(3,ip) == 0)) THEN


        d_wall = 0.0_MK
        r      = 1.0_MK


        ! This is completely general (no assumptions on geometry)
        ! Using geometry knowledge (above) limits the number of SPH particles for hp calculation even more, but 1. it is not general 2. with sorted arbitrary colloid the difference is small
        IF ( (func >= 1) .OR. & ! i.e. tracer
             ((func == 0) .AND. (this%facet(ip) == 0)) & ! i.e. if particle then only those marked before (in tracers_interpolate_velocity)
             ) THEN

           col_index = 0 ! this means all colloids will be searched
           CALL colloid_arbitrary_distance(colloids, &
                .TRUE., & ! use sort data
                col_index, & 
                this%x(1:num_dim, ip), &  ! input position
                distance, &  ! output distance
                v_normal, &  ! output normal vector
                fi,       &  ! output facet
                stat_info_sub)

           this%h_p(ip) = distance
           this%n_vec(1:num_dim, ip) = v_normal
           this%facet(ip) = fi
           this%col_i(ip) = col_index


           IF (func >= 1) THEN


              ! avoid wall penetration
              IF (this%h_p(ip) < 0.0_MK) THEN    

                 this%id(1, ip) = 0  ! delete it
                 this%id(2, ip) = -1 ! and delete without insertion into SPH particle

              ELSE IF (this%h_p(ip) < 2.0_MK * a_p) THEN

                 this%x(1:num_dim, ip) = this%x(1:num_dim, ip) + &
                      2.0_MK * (4.0_MK * a_p - this%h_p(ip)) * this%n_vec(1:num_dim, ip)

                 !this%v(1:num_dim, ip) = this%v(1:num_dim, ip) - &
                 !     DOT_PRODUCT(this%v(1:num_dim, ip), this%n_vec(1:num_dim, ip)) * this%n_vec(1:num_dim, ip)

                 this%h_p(ip) = 8.0_MK * a_p - this%h_p(ip)

              END IF


              ! See if any change in ID is required:            
              ! id: positive >> normal
              !     zero     >> delete
              !     negative >> deposited
              IF (this%h_p(ip) > d_near_wall) THEN

                 this%id(1,ip) = 0


              ELSE IF ((func == 2) .AND. & 
                   (this%h_p(ip) >  a_p         ) .AND. &
                   (this%h_p(ip) <= a_p * 4.0_MK) .AND. &
                   (this%id(1,ip) /= this%facet(ip))) THEN ! deposition


                 ! make sure this tracer won't be tried for deposition on this facet after this
                 this%id(1, ip) = this%facet(ip)

                 ! evaluate the probability of adhesion
                 Pa = 0.0_MK
                 t2 = 0.0_MK 

!!$                 ! calculate component of total force in tangent plane (parallel to wall)
!!$                 ft(1:num_dim) = this%f(1:num_dim, ip) - &
!!$                      DOT_PRODUCT(this%f(1:num_dim, ip), this%n_vec(1:num_dim, ip)) * &
!!$                      this%n_vec(1:num_dim, ip)
!!$
!!$                 f = DOT_PRODUCT(ft(1:num_dim), ft(1:num_dim)) ** 0.50_MK 


                 ! calculate component of Stokes drag force in tangent plane (parallel to wall)
                 ! find and copy corresponding interface points
                 i_index(1:vperf) = coll_f_vlist(this%col_i(ip), 1:vperf, this%facet(ip))

                 v_point = 0.0_MK
                 DO i = 1, vperf
                    v_point(1:num_dim) = v_point(1:num_dim) + &
                         this%vc_interface(this%col_i(ip), 1:num_dim, i_index(i))
                 END DO
                 v_point = v_point / dble(vperf)

                 vn = DOT_PRODUCT(v_point, this%n_vec(1:num_dim, ip))
                 vt = DOT_PRODUCT(v_point - vn * this%n_vec(1:num_dim, ip), &
                      v_point - vn * this%n_vec(1:num_dim, ip)) ** 0.5_MK

                 ! shear rate
                 S = vt / d_near_wall  
                 f = 6.0_MK * mcf_pi * a_p * this%h_p(ip) * eta * S * FS


                 DO in = 1, nbond
                    t1 = 1.0_MK
                    DO im = 1, in
                       Bk = K0 * &
                            exp(-(baa * f / dble(im) / kB / T)**bab) / &
                            (1.0_MK + bac* (baa * f / dble(im) / kB / T)**bad)

                       t1 = t1 * Bk / dble(im)
                    END DO
                    t2 = t2 + t1
                 END DO
                 Pa = 1.0_MK - (1.0_MK + t2)**(-1.0_MK)


                 ! determine if adhesion occurs 
                 uni_rnd = random_random_uniform(this%random, stat_info_sub)
                 IF (uni_rnd <= Pa) THEN

                    this%id(1, ip) = -abs(this%id(1, ip))

                 END IF

              END IF

           END IF

        END IF

     END IF

  END DO

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE         


  IF(ASSOCIATED(dmin_phys))THEN
     DEALLOCATE(dmin_phys)
  END IF

  IF(ASSOCIATED(dmax_phys))THEN
     DEALLOCATE(dmax_phys)
  END IF

  IF(ASSOCIATED(dx)) THEN
     DEALLOCATE(dx)
  END IF 

  IF(ASSOCIATED(coll_f_vlist))THEN
     DEALLOCATE(coll_f_vlist)
  END IF
  

  RETURN

END SUBROUTINE tracers_compute_wall_distance
