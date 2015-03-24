!--------------------------------
!< Subroutine that computes the first particle interactions:
!< Density change
!--------------------------------
subroutine comp_particleloop_1

  use mod_data_global
  use mod_data_ctrl
  use mod_data_prtl
  use mod_data_sph

  use mod_mapping
  use mod_io
  use mod_util

  use ppm_module_map_part

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'comp_particleloop_1'
  integer                     :: info

  integer                     :: ip, jp
  integer                     :: nneigh
  real(mk)                    :: dij, rdij
  real(mk)                    :: w, gradW

  real(mk), dimension(__DIM)  :: xip
  real(mk), dimension(__DIM)  :: dx

  !<-----------------------------
  !<-----------------------------
  call timer_start(wct_timestep(ploop1idx),info)

  !<-----------------------------
  !< Density summation
  !<-----------------------------
#define __DENSITYSUMMATION
  if (symmetry) then
#define __SYMMETRY
#include "inc/particleloop_1.inc"
#undef __SYMMETRY
    else
#include "inc/particleloop_1.inc"
  endif
#undef __DENSITYSUMMATION

  !<-----------------------------
  !<-----------------------------
  call timer_stop(wct_timestep(ploop1idx),info)


  !<-----------------------------
  !< PUT
  !< drhodt used also for sum_j Wij
  !<-----------------------------
  call timer_start(wct_timestep(putidx),info)
  call map_ghost_put(myname,realvec1=drhodt)
  call ppm_map_part_load(info) ! After a put load mapping !!!
  call timer_stop(wct_timestep(putidx),info)

end subroutine comp_particleloop_1
