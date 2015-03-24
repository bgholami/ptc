SUBROUTINE tracers_interpolate_velocity(this, d_particles, V_INT, stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_interpolate_velocity
  !----------------------------------------------------
  !
  ! Purpose     : interpolates SPH velocity to interface
  !                  
  !
  ! Routines    : 
  !
  !
  ! Remarks     :  For now works only with non-symmetric
  !                inter-processor communication (i.e. it 
  !                is nearly as slow as it can get!).
  !
  ! Revisions   :  V0.3 05.02.2013
  !                3D and arbitrary geometry
  ! 
  !                V0.2 18.10.2011
  !                concentration interpolation added
  !
  !                V0.1 11.05.2011, original
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

  USE ppm_module_typedef, ONlY : ppm_t_clist

  !----------------------------------------------------
  ! Arguments :
  !----------------------------------------------------

  TYPE(Particles)             , INTENT(IN)     :: this        
  TYPE(Particles)             , INTENT(INOUT)  :: d_particles        
  REAL(MK), DIMENSION(:, :, :), INTENT(OUT)    :: V_INT
  INTEGER                     , INTENT(OUT)    :: stat_info


  !----------------------------------------------------
  ! Control parameters :
  !----------------------------------------------------

  INTEGER                                      :: rhs_density_type

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! num_dim        : number of dimension.
  ! cut_off        : compact support domain.
  ! cut_off2       : cut_off * cut_off.
  !
  !----------------------------------------------------

  INTEGER                                      :: num_dim
  REAL(MK)                                     :: cut_off
  REAL(MK)                                     :: cut_off2 


  !----------------------------------------------------
  ! Technique parameters :
  !
  ! cell list        : cell list
  ! min_sub, max_sub : min and max extent of subdomain
  ! num_sub          : number of sumdomains on this process
  ! (the rest we don't need here)
  !---------------------------------------------------- 

  TYPE(ppm_t_clist),DIMENSION(:),POINTER       :: cell_list
  REAL(MK), DIMENSION(:), POINTER              :: min_sub
  REAL(MK), DIMENSION(:), POINTER              :: max_sub 
  INTEGER                                      :: num_sub
  INTEGER, DIMENSION(:,:), POINTER             :: sub_bcdef
  INTEGER, DIMENSION(:,:), POINTER             :: inp
  INTEGER, DIMENSION(:,:), POINTER             :: jnp
  INTEGER                                      :: nnp 


  !----------------------------------------------------
  ! cell counters and indices
  !---------------------------------------------------- 

  INTEGER                                      :: n1
  INTEGER                                      :: n2  
  INTEGER                                      :: icstart, icend, ci
  INTEGER                                      :: jcstart, jcend, cj
  INTEGER                                      :: kcstart, kcend, ck
  INTEGER                                      :: ccell   
  INTEGER                                      :: istart, iend, ipart
  INTEGER                                      :: jstart, jend, jpart
  INTEGER                                      :: il, ih, jl, jh, kl, kh


  !----------------------------------------------------
  ! number of particles
  !----------------------------------------------------

  INTEGER                                      :: num_part_real_particles


  !----------------------------------------------------
  ! Indices of particles and tracers
  !----------------------------------------------------

  INTEGER                                      :: ip, jp, it, ic
  INTEGER                                      :: i, j, k


  !----------------------------------------------------
  ! tracer-particle distance.
  !----------------------------------------------------

  REAL(MK), DIMENSION(:), ALLOCATABLE          :: r
  REAL(MK)                                     :: r_norm


  !----------------------------------------------------
  ! kernel parameters.
  !----------------------------------------------------

  REAL(MK)                                     :: w 


  !----------------------------------------------------
  ! local variables
  !----------------------------------------------------

  REAL(MK), DIMENSION(:), ALLOCATABLE          :: x_t, x_p
  REAL(MK), DIMENSION(:), ALLOCATABLE          :: f_value, res1
  REAL(MK)                                     :: d_value
  REAL(MK), DIMENSION(:), ALLOCATABLE          :: f0, r1, f
  REAL(MK), DIMENSION(:,:), ALLOCATABLE        :: f1, r2
  REAL(MK), DIMENSION(:,:), ALLOCATABLE        :: a
  REAL(MK)                                     :: r0, mass_factor, denom  
  INTEGER                                      :: num_col, val_dim, nmh
  INTEGER , DIMENSION(:), ALLOCATABLE          :: cell_index, ngl
  REAL(MK), DIMENSION(:), ALLOCATABLE          :: cell_dx
  INTEGER                                      :: stat_info_sub


  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0

  NULLIFY(min_sub)
  NULLIFY(max_sub)  
  NULLIFY(cell_list)
  NULLIFY(sub_bcdef)        
  NULLIFY(inp)
  NULLIFY(jnp)


  !----------------------------------------------------
  ! Control parameters :
  !
  ! Get control variables.
  !----------------------------------------------------

  rhs_density_type = &
       control_get_rhs_density_type(this%ctrl,stat_info_sub)

  !----------------------------------------------------
  ! Physics parameters :
  !
  ! from a object of Physics class.
  !
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub)
  cut_off     = &
       physics_get_cut_off(this%phys,stat_info_sub)
  cut_off2 = cut_off * cut_off 


  !----------------------------------------------------
  ! Get the boundary of this sub-domain. 
  !----------------------------------------------------

  CALL technique_get_min_sub(this%tech,min_sub,stat_info_sub)
  CALL technique_get_max_sub(this%tech,max_sub,stat_info_sub)


  !----------------------------------------------------
  ! Get particles' cell list and prepare data to rank
  !----------------------------------------------------        

  CALL technique_get_cell_list(d_particles%tech, &
       num_sub, cell_list, &
       inp, jnp, nnp, sub_bcdef, stat_info_sub)


  IF (num_sub .GT. 1) THEN
     print*, 'Tracers_interpolate_velocity : ', &
          'technique_get_min_sub returns one set of dimension only.', &
          'More than one subdomain on each process is not supported!'
     stat_info = -1
     GOTO 9999
  END IF

  ! allocate arrays
  allocate(cell_index(num_dim), ngl(num_dim))
  allocate(cell_dx(num_dim))

  DO i = 1, num_dim 
     nmh = INT((max_sub(i) - min_sub(i)) / d_particles%tech%ghost_size)

     ngl(i)     = cell_list(num_sub)%nm(i) - nmh - 1
     cell_dx(i) = (max_sub(i) - min_sub(i)) / REAL(nmh, MK)
  END DO


  ! set cell index limits
  ! symmetry, so cell indices starts from 0
  icstart = 0              
  jcstart = 0
  kcstart = 0

  icend = cell_list(num_sub)%nm(1)-2
  jcend = cell_list(num_sub)%nm(2)-2 

  n1 = cell_list(num_sub)%nm(1)

  IF ( num_dim == 2 ) THEN

     kcend = kcstart
     n2 = 0

  ELSE

     kcend = cell_list(num_sub)%nm(3)-2      
     n2 = cell_list(num_sub)%nm(1) * &
          cell_list(num_sub)%nm(2)

  END IF


  !----------------------------------------------------
  ! prepare the interpolation array
  !----------------------------------------------------

  num_col = SIZE(this%x_interface, 1)

  V_INT = 0.0_MK 

  ! dimension of quantity to by interpolated
  val_dim = num_dim


  !----------------------------------------------------
  ! allocate interpolation arrays
  !----------------------------------------------------

  allocate(a(num_dim+1, num_dim+1))
  allocate(f(num_dim+1))
  allocate(f0(val_dim), f_value(val_dim))
  allocate(x_t(num_dim), x_p(num_dim), res1(num_dim))
  allocate(r(num_dim), r1(num_dim))
  allocate(f1(val_dim, num_dim), r2(num_dim, num_dim))


  !----------------------------------------------------
  ! allocate wall distance data for particles
  ! (those close to walls are marked here and 
  !  used later.)
  !----------------------------------------------------

  IF (ASSOCIATED(d_particles%h_p)) THEN
     DEALLOCATE(d_particles%h_p)
  END IF

  IF (ASSOCIATED(d_particles%n_vec)) THEN
     DEALLOCATE(d_particles%n_vec)
  END IF

  IF (ASSOCIATED(d_particles%facet)) THEN
     DEALLOCATE(d_particles%facet)
  END IF

  IF (ASSOCIATED(d_particles%col_i)) THEN
     DEALLOCATE(d_particles%col_i)
  END IF

  num_part_real_particles = d_particles%num_part_real

  ALLOCATE(d_particles%h_p(num_part_real_particles))
  ALLOCATE(d_particles%n_vec(num_dim, num_part_real_particles))
  ALLOCATE(d_particles%facet(num_part_real_particles))  
  ALLOCATE(d_particles%col_i(num_part_real_particles))  

  ! allocate facet to -1 since those close to walls
  ! will be marked with facet(jp) = 0 in this routine.
  d_particles%facet = -1

  !----------------------------------------------------
  ! interpolation loop
  !----------------------------------------------------       

  DO ic = 1, num_col

     DO ip = 1, this%num_interface(ic) 

        !! only if it is local
        !IF ( this%local_interface(ic, ip)) THEN

        f0 = 0.0d0
        f1 = 0.0d0
        r0 = 0.0d0
        r1 = 0.0d0

        r2 = 0.0d0     


        x_t(1:num_dim) = this%x_interface(ic, 1:num_dim, ip)

        ! find the cell it belongs to
        cell_index(1:num_dim) = FLOOR((x_t(1:num_dim) - min_sub(1:num_dim)) / &
             cell_dx(1:num_dim)) + ngl(1:num_dim)

        ! enforce right limits for cell index
        il = MAX(cell_index(1)-1, icstart)  
        ih = MIN(cell_index(1)+1, icend)   
        jl = MAX(cell_index(2)-1, jcstart)  
        jh = MIN(cell_index(2)+1, jcend)
        IF ( num_dim == 2 ) THEN
           kl = kcstart
           kh = kcend
        ELSE
           kl = MAX(cell_index(3)-1, kcstart)  
           kh = MIN(cell_index(3)+1, kcend)
        END IF

        ! loop over neighboring cells
        DO ck = kl, kh
           DO cj = jl, jh
              DO ci = il, ih

                 ! cell number
                 ccell = ci + 1 + n1 * cj + n2 * ck 

                 ! Get pointers of first and last particles in ccell
                 istart = cell_list(num_sub)%lhbx(ccell)
                 iend   = cell_list(num_sub)%lhbx(ccell+1)-1

                 ! loop over particles in ccell
                 DO ipart = istart, iend

                    ! index of the particles (in data array)
                    jp = cell_list(num_sub)%lpdx(ipart)

                    !! only fluid SPH particles
                    !if (d_particles%id(2, jp) == 0) then

                    x_p(1:num_dim) = d_particles%x(1:num_dim, jp)
                    r(1:num_dim) = x_t(1:num_dim) - x_p(1:num_dim)
                    r_norm = DOT_PRODUCT(r(1:num_dim), r(1:num_dim))


                    IF (r_norm <= this%near_wall**2.0_MK) THEN

                       ! mark this particle for wall distance calculation
                       ! (only these particle might be close enough to walls)
                       d_particles%facet(jp) = 0 ! otherwise -1

                    END IF


                    !IF (r_norm <= cut_off2) THEN
                    IF ((this%local_interface(ic, ip)) .AND. (r_norm <= cut_off2)) THEN 

                       r_norm = SQRT(r_norm)

                       CALL kernel_kernel(d_particles%kern, r_norm, w, stat_info_sub)

                       IF (rhs_density_type == 1) THEN
                          mass_factor = d_particles%m(jp)
                       ELSE
                          mass_factor = 1.0_MK
                       END IF

                       ! function to be interpolated
                       f_value(1:num_dim) =  d_particles%v(1:num_dim, jp)
                       ! density
                       d_value            = d_particles%rho(jp) / mass_factor

#include "setup_system.inc" 

                    END IF

                    !end if

                 END DO

              END DO

           END DO

        END DO
        
        IF ( this%local_interface(ic, ip)) THEN   

#include "solve_system.inc"

           V_INT(ic, 1:num_dim, ip) = res1(1:num_dim)

        END IF

        !END IF  ! only local

     END DO

  END DO

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE 


  IF(ASSOCIATED(min_sub)) THEN
     DEALLOCATE(min_sub)
  END IF

  IF(ASSOCIATED(max_sub)) THEN
     DEALLOCATE(max_sub)
  END IF

  IF(ASSOCIATED(cell_list)) THEN
     DEALLOCATE(cell_list)
  END IF

  IF(ASSOCIATED(sub_bcdef)) THEN
     DEALLOCATE(sub_bcdef)
  END IF

  IF(ASSOCIATED(inp)) THEN
     DEALLOCATE(inp)
  END IF

  IF(ASSOCIATED(jnp)) THEN
     DEALLOCATE(jnp)
  END IF


  RETURN

END SUBROUTINE tracers_interpolate_velocity
