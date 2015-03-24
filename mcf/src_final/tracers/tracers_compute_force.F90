SUBROUTINE tracers_compute_force(this, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_compute_force
  !----------------------------------------------------
  !
  ! Purpose     : Computing forces for  
  !               non-symmetric inter-processor 
  !               communiction
  !                  
  !
  ! Routines    :  
  !
  ! References  :  Longest et al.
  !                2004 Computers and Fluids
  !
  !
  ! Remarks     :  This should be extended to symmetric
  !                inter-processor communication
  !
  ! Revisions   :  V0.3 04.02.2013
  !                3D and arbitrary geometry
  !
  !                V0.2 21.07.2011
  !                adapted to use interface velocity
  !
  !                V0.1 11.05 2011, original
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

  TYPE(Particles), INTENT(INOUT)  :: this  
  INTEGER, INTENT(OUT)	          ::  stat_info	


  !----------------------------------------------------
  ! Control parameters :
  !
  ! symmetry      : indicate if we use symmetric
  !                 inter-process communication or not.
  !----------------------------------------------------

  LOGICAL                         :: symmetry

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  ! rho            : density
  ! eta            : viscosity
  !----------------------------------------------------

  INTEGER                         :: num_dim
  REAL(MK)                        :: rho
  REAL(MK)                        :: eta       

  !----------------------------------------------------
  ! Number of all tracers
  !----------------------------------------------------

  INTEGER                         :: num_part_real

  !----------------------------------------------------
  ! local variables:
  !  
  ! VI             : interpolated velocity of particles
  ! h_p            : tracers' wall distance
  !----------------------------------------------------

  REAL(MK), DIMENSION(:,:), ALLOCATABLE :: VI 
  REAL(MK), DIMENSION(:)  , ALLOCATABLE :: f_lub, f_lift, t_vec 
  REAL(MK), DIMENSION(:)  , ALLOCATABLE :: u_ni, v_ni, u_ti, v_ti 
  REAL(MK), DIMENSION(:)  , ALLOCATABLE :: sgn_ni, sgn_ti
  REAL(MK)                              :: a_p, tau_p, kappa
  REAL(MK)                              :: tw_sep, m_p, rho_p
  REAL(MK)                              :: v_n, u_n, v_t, u_t
  REAL(MK)                              :: u_s, shear_rate, lift_function
  REAL(MK)                              :: mag
  INTEGER                               :: i, ip, it 
  INTEGER                               :: stat_info_sub

  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0    


  !----------------------------------------------------
  ! Tracers' parameters :
  !----------------------------------------------------

  ! tracers' radius
  a_p = this%radius

  ! tracers' near-wall factor
  ! i.e. for (h_p <= 2.0 * a_p * tw_sep) tracers feel the wall
  tw_sep = 10.0_MK

  ! tracers' density ratio
  rho_p = this%rho_ratio


  !----------------------------------------------------
  ! Control parameters :
  !
  ! Get control variables.
  !----------------------------------------------------

  symmetry  = control_get_symmetry(this%ctrl,stat_info_sub)

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = physics_get_num_dim(this%phys,stat_info_sub) 

  eta         = physics_get_eta(this%phys,stat_info_sub)         

  rho         = physics_get_rho(this%phys,stat_info_sub)  


  !-------------------------------------------------
  ! Calculate tracers' density
  !-------------------------------------------------

  ! first determine tracers' density
  rho_p = rho / rho_p   

  ! momentum response time
  tau_p = 2.0_MK * rho_p * a_p**2.0_MK / 9.0_MK / eta * 1.0_MK
  this%tau_p = tau_p

  ! mass of tracers
  m_p = 4.0d0 / 3.0d0 * 3.14150d0 * a_p**3.0d0 * rho_p
  this%m_p = m_p     


  !----------------------------------------------------
  ! Number of real tracers
  !----------------------------------------------------

  num_part_real = tracers_get_num_part_real(this,stat_info_sub) 

  IF (num_part_real > 0) then


  !----------------------------------------------------
  ! Allocate memory for force,
  !
  ! for non-symmetry : allocate num_part_real
  !
  !----------------------------------------------------

  IF (ASSOCIATED(this%f)) THEN
     DEALLOCATE(this%f,STAT=stat_info_sub)
  END IF

  ALLOCATE(this%f(num_dim,num_part_real), STAT=stat_info_sub)

  this%f = 0.0_MK        


  !-------------------------------------------------
  ! Allocate memory
  !-------------------------------------------------

  ALLOCATE(f_lub(num_dim), f_lift(num_dim))
  ALLOCATE(t_vec(num_dim))  
  ALLOCATE(u_ni(num_dim), v_ni(num_dim)) 
  ALLOCATE(u_ti(num_dim), v_ti(num_dim)) 
  ALLOCATE(sgn_ni(num_dim), sgn_ti(num_dim))

!!$  !-------------------------------------------------
!!$  ! Calculate wall parameters
!!$  !------------------------------------------------- 
!!$
!!$  ! calculate this%h_p and this%n_vec
!!$  CALL tracers_compute_wall_distance(this, &
!!$       1, stat_info_sub)
    

  !-------------------------------------------------
  ! Interpolate interface velocities to the current 
  ! position of tracers
  !-------------------------------------------------        

  ALLOCATE(VI(num_dim, num_part_real))

  CALL tracers_near_wall_velocity(this, VI, stat_info_sub)  

  IF(stat_info_sub /= 0 ) THEN
     PRINT *, "tracers_compute_force:", "near_wall velocity failed. " 
          stat_info = -1
     GOTO 9999
  END IF


  !-------------------------------------------------
  ! Calculate the drag force on tracers
  !-------------------------------------------------
#include "tracers_compute_force_drag.inc"               


  !-------------------------------------------------
  ! Loop too calculate lubrication and lift force
  !-------------------------------------------------  

  DO ip = 1, num_part_real

     ! if deposition, according to: 
     ! id: positive >> normal
     !     zero     >> delete
     !     negative >> deposited
     if ( this%id(1, ip) > 0 ) then

        kappa = a_p / this%h_p(ip)

        ! surface_tangent unit vector in direction of tracer velocity
        t_vec(1:num_dim) = this%v(1:num_dim, ip) - &
             DOT_PRODUCT(this%v(1:num_dim, ip), this%n_vec(1:num_dim, ip)) * &
             this%n_vec(1:num_dim, ip)
        mag = DOT_PRODUCT(t_vec(1:num_dim), t_vec(1:num_dim))
        mag = mag ** 0.50_MK   
        if (mag /= 0.0_MK) then
           t_vec(1:num_dim) = t_vec(1:num_dim) / mag 
        else
           t_vec(1:num_dim) = 0.0_MK
        end if


        ! wall-normal and wall-tangent components of velocities 
        ! (i.e. SPH velocity (u) and tracer velocity (v))

        u_n = DOT_PRODUCT(VI(1:num_dim,ip)    , this%n_vec(1:num_dim, ip))
        v_n = DOT_PRODUCT(this%v(1:num_dim,ip), this%n_vec(1:num_dim, ip))     
        
        u_t = DOT_PRODUCT(VI(1:num_dim,ip)    , t_vec(1:num_dim))
        v_t = DOT_PRODUCT(this%v(1:num_dim,ip), t_vec(1:num_dim))
        
        u_ni(1:num_dim) = u_n * this%n_vec(1:num_dim, ip)
        v_ni(1:num_dim) = v_n * this%n_vec(1:num_dim, ip)
        
        u_ti(1:num_dim) = u_t * t_vec(1:num_dim)
        v_ti(1:num_dim) = v_t * t_vec(1:num_dim)

        ! evaluate sgn_n and sgn_t
        ! The third argument determines if 
        !           sgn_n should be calculated (1), or if
        !           sgn_t should be calculated (2)
        CALL tracers_eval_sgn(v_n, v_ni(1:num_dim), 1, sgn_ni(1:num_dim)) 
        CALL tracers_eval_sgn(v_t, v_ti(1:num_dim), 2, sgn_ti(1:num_dim))


        !-------------------------------------------------
        ! Calculate drag modifications
        !-------------------------------------------------

        ! force cap (if too close to walls)
        ! minimum allowed distance between EDGE of tracer and wall is 5% of radius
        !kappa = MAX(kappa, 1.0_MK/1.05_MK)


        IF ((kappa < 1.0_MK) .AND. &
             (kappa >= 1.0_MK/(2.0_MK*tw_sep))) THEN ! also, if too far from wall no drag modification  
           !-------------------------------------------------
           ! Calculate the lubrication force on tracers
           !-------------------------------------------------                  
#include "tracers_compute_force_lubrication.inc"        
           
           
           !-------------------------------------------------
           ! Calculate the lift force on tracers
           !-------------------------------------------------
#include "tracers_compute_force_lift.inc"           
        END IF
        
                
     else
        this%v(:, ip) = 0.0_MK
        this%f(:, ip) = 0.0_MK

     end if
     
  END DO

END IF



  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------
  
9999 CONTINUE        


  RETURN

END SUBROUTINE tracers_compute_force
