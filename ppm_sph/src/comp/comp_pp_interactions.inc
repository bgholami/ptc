!--------------------------------
!< Subroutine where force calculation is performed
!--------------------------------
subroutine comp_pp_interactions

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_prtl
  use mod_data_physics

  use mod_io
  use mod_util

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'comp_pp_interactions'
  !character(maxchar)          :: cbuf
  integer                     :: ip, nneigh, jp, j
  real(mk)                    :: dij, rdij
  real(mk), dimension(__DIM)  :: xip, vip
  real(mk), dimension(__DIM)  :: dx, dvel

  real(mk)                    :: w, gradW
  real(mk)                    :: p_avg
  real(mk)                    :: rhoip, rhojp
  real(mk)                    :: volip, voljp
  real(mk)                    :: pressip, pressjp
  real(mk)                    :: part_dens_ip, part_dens_jp
  real(mk)                    :: etaip
  real(mk)                    :: eta_avg
  real(mk)                    :: Vi2Vj2

  real(mk), dimension(__DIM)  :: pforce

  !----------------------------
  !< Loop: FOR EACH PURE PHASE
  !----------------------------
  if (symmetry) then
#define __SYMMETRY
#include "inc/pp_interaction_fluid.inc"
#undef __SYMMETRY
  else
#include "inc/pp_interaction_fluid.inc"
  endif

end subroutine comp_pp_interactions
