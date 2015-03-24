!--------------------------------
!< Calc constants as timestep, ref pressure,...
!--------------------------------
subroutine comp_constants

  use mod_data_global
  use mod_data_sph
  use mod_data_physics
  use mod_data_ctrl
  use mod_data_prtl

  use mod_util

  use mod_io

  use ppm_module_data, only : ppm_pi_d, ppm_pi_s, ppm_kind_double

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)    :: myname = 'comp_constants'
  character(maxchar)    :: cbuf
  real(mk)              :: force
  real(mk)              :: dt_min

  real(mk)              :: dij, w, gradW
  integer               :: i, j

  !-----------------------------
  !< Body of routine
  !-----------------------------


  !-----------------------------
  !< PI
  !-----------------------------
  if (mk .eq. ppm_kind_double) then
    pi = ppm_pi_d
  else
    pi = ppm_pi_s
  endif
  pi2 = pi*pi
  pi_inv = 1._mk/pi

  !-----------------------------
  !< Performance analysis
  !-----------------------------
  allocate(gather_nproc(Nproc))
  gather_nproc = 0.


  !-----------------------------
  !< Define particle spacing and smoothinglength...
  !-----------------------------
  smoothinglength = cutoff / ratio_cutoff_smoothing
  smoothinglength2 = smoothinglength*smoothinglength
  smoothinglength3 = smoothinglength2*smoothinglength
  cutoff_inv = 1._mk / cutoff
  cutoff2 = cutoff*cutoff
  dpx = smoothinglength / ratio_smoothing_dx

  !-------------------------------------
  !< Tolerance level
  !-------------------------------------
  myeps = 1.e-30_mk

  !-------------------------------------
  !< Skin size of verlet list -> check manually if work has to be done
  !-------------------------------------
  verlet_skin = verlet_skinfactor * dpx
  verlet_skin2 = verlet_skin*verlet_skin

  !-----------------------------
  !< Box size and ghost boxes
  !-----------------------------
  sph_bsize = cutoff + 2._mk*verlet_skin

  !-----------------------------
  !< Verlet list
  !-----------------------------
  !< all + N_PHASES
  NVERLETLIST = 1 + N_PHASES
  allocate(verletlist(NVERLETLIST))

  !-----------------------------
  !< KERNEL PARAMETER
  !-----------------------------
#ifndef PARAM_KERNEL_TYPE_QUINTICSPLINE
#define PARAM_KERNEL_TYPE_QUINTICSPLINE 1
#endif
#ifndef PARAM_KERNEL_TYPE_BSPLINE
#define PARAM_KERNEL_TYPE_BSPLINE 2
#endif
#ifndef PARAM_KERNEL_TYPE_CUBIC
#define PARAM_KERNEL_TYPE_CUBIC 3
#endif
#ifndef PARAM_KERNEL_TYPE_WENDLANDQUARTIC
#define PARAM_KERNEL_TYPE_WENDLANDQUARTIC 4
#endif


#if __KERNEL==PARAM_KERNEL_TYPE_QUINTICSPLINE
#if __DIM==2
  norm = 63.0_mk/478.0_mk*pi_inv
  factorW = norm * (cutoff_inv * cutoff_inv)
  factorGradW = -15._mk*norm*(cutoff_inv*cutoff_inv*cutoff_inv)
#elif __DIM==3
  norm = 27.0_MK/120.0_mk*pi_inv
  factorW = norm * (cutoff_inv*cutoff_inv*cutoff_inv)
  factorGradW = -15.0_mk*norm*(cutoff_inv*cutoff_inv*cutoff_inv*cutoff_inv)
#endif

#elif __KERNEL==PARAM_KERNEL_TYPE_CUBIC
#if __DIM==2
  factorW = 10._MK / 7._MK *pi_inv * (cutoff_inv * ratio_cutoff_smoothing)**2._MK
  factorGradW = 10._MK/ 7._MK *pi_inv *(cutoff_inv * ratio_cutoff_smoothing)**3._MK
#elif __DIM==3
  call abortmessage(myname,'factorW for cubic 3d not yet implemented.')
#endif

#elif __KERNEL==PARAM_KERNEL_TYPE_BSPLINE
#if __DIM==2
  norm = 40.0_MK /  7.0_MK * pi_inv
#elif __DIM==3
  norm = 8._MK*pi_inv
#endif

#elif __KERNEL==PARAM_KERNEL_TYPE_WENDLANDQUARTIC
#if __DIM==2
  factorW = 7._MK / 64._MK / pi / (smoothinglength**2._MK)
  factorGradW = 7._MK/ 64._MK / pi /(smoothinglength**3._MK)
#endif

#endif

  !------------------------------------------
  ! Calc kernel W(x=0) for density summation
  !-------------------------------------------
  dij = 0._mk
#include "src/comp/inc_kernelfunction.inc"
  w_zero = w

  dij = dpx
#include "src/comp/inc_kernelfunction.inc"
  w_dx = w

  !----------------------------------------
  !< Calc n_entries for array of pointers for mapping
  !--------------------------------------
  allocate(realmaparray(ldamaparray))
  allocate(intmaparray(ldamaparray))
  call reset_ptr_to_list(realmaparray)
  call reset_ptr_to_list(intmaparray)

  !----------------------------------------
  !< Calculate some quantities
  !----------------------------------------
  do i = 1, N_PHASES
    phase(i)%visc_kin = phase(i)%visc_dyn / phase(i)%density 
  enddo

  !----------------------------------------
  !<  Sound speeds
  !----------------------------------------
  DO i = 1, N_PHASES
    phase(i)%soundspeed(:) = 0.0_MK
    !< Sound speed -> velocity
    phase(i)%soundspeed(idx_sound_vel) = sqrt(abs(vel_ref)**2.0_MK / phase(i)%delta_dens)
    !< Sound speed -> viscosity
    phase(i)%soundspeed(idx_sound_visc) = &
    sqrt(phase(i)%visc_kin*abs(vel_ref)/length_ref/phase(i)%delta_dens)
    !< Sound speed -> body force
    force = bodyforce(1)**2.0_MK + bodyforce(2)**2.0_MK
#if __DIM==3
    force = force + bodyforce(3)**2.0_MK
#endif
    force = sqrt(force)
    phase(i)%soundspeed(idx_sound_force) = sqrt(force*length_ref / phase(i)%delta_dens)
  END DO !< Loop over all phases


  !----------------------------------------
  !<  Maximum Sound speeds
  !----------------------------------------
  DO i = 1, N_PHASES
    phase(i)%cmax = phase(i)%soundspeed(NCSOUNDCRITERIA)
    DO j = 1, NCSOUNDCRITERIA-1
      phase(i)%cmax = MAX(phase(i)%cmax,phase(i)%soundspeed(j))
    END DO
  ENDDO

  !----------------------------------------
  ! Reference pressure
  !----------------------------------------
  DO i = 1, N_PHASES
    SELECT CASE (phase(i)%type)
      CASE (param_type_fluid)
        phase(i)%press_ref = phase(i)%density*phase(i)%cmax*phase(i)%cmax/&
        phase(i)%gamma
      CASE DEFAULT
        phase(i)%press_ref = 0.
    END SELECT
    phase(i)%press_b = phase(i)%press_b*phase(i)%press_ref
  ENDDO

  !----------------------------------------
  !  timesteps
  !----------------------------------------
  DO i = 1, N_PHASES
    phase(i)%timestep(:) = -1._mk
    ! dt -> velocity condition
    phase(i)%timestep(idx_dt_vel) = 0.25_MK*smoothinglength / (phase(i)%cmax + abs(vel_ref))
    ! dt -> viscous condition
    if (phase(i)%visc_kin .gt. myeps) phase(i)%timestep(idx_dt_visc) = 0.125_MK * smoothinglength2 / (phase(i)%visc_kin + myeps)
    ! dt -> body force condition
    if (force .gt. myeps) then
      phase(i)%timestep(idx_dt_force) = 0.25_MK * sqrt(smoothinglength / (force + myeps))
    endif
  ENDDO

  !----------------------------------------
  !< Minimum timestep
  !----------------------------------------
  dt_min = 1.e99_MK
  DO i = 1,N_PHASES
    DO j = 1, NDTCRITERIA
      if (phase(i)%timestep(j) .gt. 0._MK) dt_min = MIN(dt_min, phase(i)%timestep(j))
    END DO
  END DO

  dt_global = dt_min
  rdt_global = 1._MK / dt_global


  !------------------------------------------------------------------
  !< Restart filename -> Must be after read_ctrl...that's why not
  !< in defaults...
  !------------------------------------------------------------------
  select case (choice_restart)
    case (param_restart_ascii)
      write(cbuf,'(A)') '/restart.ascii'

    case (param_restart_binary)
      write(cbuf,'(A)') '/restart.ufmt'

    case (param_restart_none)
      !< do nothing
      write(cbuf,'(A)') '/norestart'

    case default
      call abortmessage(myname,'Wrong choice_restart, not existing')

  end select

  i = len_trim(trim(thispath)) + len_trim(cbuf)
  call check_char_length(i,'comp_constants','restartfile')
  write(restartfile,'(A)') trim(thispath)//trim(cbuf)
  cbuf = trim(cbuf)//'old'
  i = len_trim(trim(thispath)) + len_trim(cbuf)
  call check_char_length(i,'comp_constants','oldrestartfile')
  write(oldrestartfile,'(A)') trim(thispath)//trim(cbuf)

  !-----------------------------
  !< Remove old output
  !-----------------------------
  if (.not. restart .and. rank .eq. 0) then
    write(cbuf,'(A)') 'rm -f '//trim(outputfile)//'*'
    call system(trim(cbuf))
    write(cbuf,'(A)') 'rm -f '//trim(ghostfile)//'*'
    call system(trim(cbuf))
  endif

end subroutine comp_constants
