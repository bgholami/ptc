SUBROUTINE particles_integrate_velocity(this,&
     num,dt,lambda,stat_info)
  !----------------------------------------------------
  ! Subroutine  : particles_integrate_velocity
  !----------------------------------------------------
  !
  ! Purpose     : Integrate velocity of particles
  !               using accleration f with lamda, 
  !               which is a prediction or
  !               correction coefficient.
  !                  
  !                  
  ! Reference   :
  !
  ! Remark      : colloidal boundary particle may rotate,
  !               needs to be done seperately.
  !               wall boundary particles may have
  !               different behavior, such as,
  !               oscillating shear, needs to be done
  !               seperately.
  !
  ! Revisions   : V0.4 22.11.2010, integrate velocity of
  !               only fluid particles.
  ! 
  !               V0.3 24.08.2010, integrate velocity of
  !               particles, except colloidal boundary
  !               particles.
  !
  !               V0.2 09.07.2009, 
  !               check again the work flow is correct and
  !               supply with more comments for code.
  !
  !               V0.1 01.04.2009, original version.
  !
  !----------------------------------------------------
  ! Author      : Xin Bian
  ! Contact     : xin.bian@aer.mw.tum.de
  ! 
  ! Dr. Marco Ellero's Emmy Noether Group,
  ! Prof. Dr. N. Adams' Chair of Aerodynamics,
  ! Faculty of Mechanical Engineering,
  ! Technische Universitaet Muenchen, Germany.
  !----------------------------------------------------


  !----------------------------------------------------
  !  Arguments
  !
  !  this       : an object of Particles Class.
  !  num        : first num of particles needed to be updated.
  !  dt         : time step.
  !  lamda      : coefficient of acceleration.
  !  stat_info  : return flag of status.
  !----------------------------------------------------

  TYPE(Particles), INTENT(INOUT)          :: this
  INTEGER, INTENT(IN)                     :: num
  REAL(MK), INTENT(IN)                    :: dt
  REAL(MK), INTENT(IN)                    :: lambda
  INTEGER, INTENT(OUT)                    :: stat_info

  !----------------------------------------------------
  !  Local variables
  !----------------------------------------------------

  INTEGER                                 :: stat_info_sub
  INTEGER                                 :: dim,i
  REAL(MK)                                :: time
  REAL(MK), DIMENSION(2,60)               :: vmx
  INTEGER                                 :: index
  REAL(MK)                                :: vmax, vr, vl, dt_pulse, t_add, tstart, rad, radius

  !----------------------------------------------------
  ! other local variables:
  !
  !----------------------------------------------------
  integer                               :: it, ii, jp
  REAL(MK), DIMENSION(3)          :: x_t, x_p, v_p
  INTEGER                         :: d, rc, rhs_density_type 
  REAL(MK)                        :: f_value, d_value
  REAL(MK), DIMENSION(2)          :: f1, r1
  REAL(MK), DIMENSION(3)          :: f
  REAL(MK), DIMENSION(3,3)        :: a
  REAL(MK), DIMENSION(2,2)        :: r2
  REAL(MK)                        :: f0, r0, res1, mass_factor  
  REAL(MK), DIMENSION(3)          :: r
  REAL(MK)                        :: r_norm
  REAL(MK)                        :: w 
  REAL(MK)                        :: cut_off2, cut_off 

  !----------------------------------------------------
  !  Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0


  dim  = physics_get_num_dim(this%phys,stat_info_sub)  

#include "hinds_vmax.inc"
  dt_pulse = vmx(1, 2) - vmx(1, 1)
  rad = 0.6350_MK / 2.0_MK
  tstart = dt * 30000.0_MK


  cut_off     = &
       physics_get_cut_off(this%phys,stat_info_sub)
  cut_off2 = cut_off * cut_off 

  !----------------------------------------------------
  ! Update velcoity of fluid particles.
  !----------------------------------------------------

  DO i = 1, num

     IF ( this%id(this%sid_idx,i) == mcf_particle_type_fluid ) THEN

        IF (this%id(this%bid_idx, i) == 0) THEN

           this%v(1:dim,i) = &
                this%v(1:dim,i) + &
                lambda * this%f(1:dim,i) * dt

        ELSE ! i.e. inflow/outflow boundary particle

           this%v(1:dim, i) = 0.0_MK
           ! and just to be safe
           this%f(1:dim, i) = 0.0_MK
           ! and this (since no periodic BC)
           this%f_bp(1:dim, i) = 0.0_MK


!!$           ! Poiseuille (morris)
!!$           if ((this%x(2, i) >= 0.0_MK) .AND. & 
!!$                (this%x(2, i) <= 0.0010_MK)) then
!!$
!!$              if (this%id(this%bid_idx, i) == 1) then
!!$                 ! inflow
!!$                 this%v(1, i) = -50.0_MK * this%x(2, i)**2.0 + &
!!$                      0.050_MK * this%x(2, i)
!!$              end if
!!$
!!$           end if

!!$           ! Taper 2D
!!$           radius = ABS(this%x(2, i))
           ! Taper 3D
           radius = (this%x(2, i)**2.0_MK + this%x(3, i)**2.0_MK)**0.50_MK !3D   
           if (radius < 0.31750_MK) then 
              

              !if (this%id(this%bid_idx, i) == 1) then
              ! inflow

!!$              ! steady...
!!$              this%v(1, i) = -25.86920_MK * this%x(2, i)**2.0 - &
!!$                   0.000916610_MK * this%x(2, i) + 2.60810_MK
!!$
!!$              time = physics_get_time_current(this%phys, stat_info_sub)
!!$              
!!$              this%v(1, i) = this%v(1, i) * &
!!$                   (1.0_MK - EXP(-1.0_MK * time / 0.5720_MK))


              ! pulsatile (Hinds)
                    time = physics_get_time_current(this%phys, stat_info_sub)
                    IF (time >= tstart) THEN

#if __TRACER
                       CALL physics_set_tracers_c_factor(this%phys, 5, stat_info_sub)
#endif  
                       ! find corresponding vmax
                       index = FLOOR((time-tstart) / dt_pulse) + 1
                       t_add = (time-tstart) - (index - 1) * dt_pulse
                       index = MOD(index, 60) ! pulsatility

                       if (index /= 0) then
                          vl = vmx(2, index)
                          vr = vmx(2, index + 1)
                       else
                          vl = vmx(2, 60)
                          vr = vmx(2, 1)
                       end if

                       ! linear interpolation between two data points
                       vmax = vr * (t_add/dt_pulse) + vl * (1.0_MK - t_add/dt_pulse)

                       ! velocity profile
                       this%v(1, i) = -vmax * ((radius / rad )**2.0_MK - 1.0_MK)

                    ELSE

                       vmax = vmx(2, 1) * (1.0_MK - EXP(-1.0_MK * time / 1.1440_MK))                    

                       this%v(1, i) = -vmax * ((radius / rad )**2.0_MK - 1.0_MK)                       

                       
                    END IF

                    ! else ! i.e. outflow

                    !#include "bad_patch.inc"

                    !                    end if

           else ! inflow/outflow but not moving at all
              this%f_bp(1:dim, i) = 0.0_MK
           end if

        END IF
        
     END IF

  END DO ! i =1, num

9999 CONTINUE

  RETURN

END SUBROUTINE particles_integrate_velocity

