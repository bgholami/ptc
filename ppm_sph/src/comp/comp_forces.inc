!------------------------------------
!< Subroutine that calculates forces
!< and density change and puts everything
!< to the other procs...
!------------------------------------
subroutine comp_forces

  use mod_data_global
  use mod_data_prtl
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_physics

  use mod_io
  use mod_mapping
  use mod_util


  use ppm_module_util_dbg
  use ppm_module_time
  use ppm_module_map_part_util

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  !------------------------------
  !< Arguments
  !------------------------------

  !-------------------------------
  !< Local Variables
  !--------------------------------
  character(maxchar)        :: myname = 'comp_forces'
  integer   :: info

  integer   :: ip
  real(mk)  :: rmass

#ifdef __MPI
  integer   :: ierr
#endif

  !---------------------------------
  !< PARTICLE INTERACTIONS PER PHASE
  !---------------------------------
  call timer_start(wct_timestep(forceidx),info)
  call comp_pp_interactions

  !--------------------------------
  !< Calc acceleration from forces -> include ghost particles
  !< and put them to the other procs
  !--------------------------------
  do ip = 1, Mpart
    rmass = 1._mk / pdata(massidx,ip)
    d2xpdt2(1,ip) = pForce_total(1,ip)*rmass
    d2xpdt2(2,ip) = pForce_total(2,ip)*rmass
#if __DIM==3
    d2xpdt2(3,ip) = pForce_total(3,ip)*rmass
#endif
  enddo
  call timer_stop(wct_timestep(forceidx),info)

  !------------------------------------
  !< Put ghost contributions to procs
  !------------------------------------
  call timer_start(wct_timestep(putidx),info)
  call map_ghost_put(myname,realvec1=d2xpdt2)
  call ppm_map_part_load(info)
  call timer_stop(wct_timestep(putidx),info)

  !------------------------------------
  !< Body force -> absolut / relative force
  !------------------------------------
  call timer_start(wct_timestep(bodyforceidx),info)
  do ip = 1, Npart
    d2xpdt2(1,ip) = d2xpdt2(1,ip) - bodyforce(1)
    d2xpdt2(2,ip) = d2xpdt2(2,ip) - bodyforce(2)
#if __DIM==3
    d2xpdt2(3,ip) = d2xpdt2(3,ip) - bodyforce(3)
#endif
  enddo
  call timer_stop(wct_timestep(bodyforceidx),info)

  if (sph_debug .gt. 0) call logmessage(myname,'done...')

  return

end subroutine comp_forces


