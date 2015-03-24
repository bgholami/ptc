!--------------------------------
!< Subroutine: initial number density and interface flags calculated
!< Define volume of a particle here
!--------------------------------
subroutine comp_init_state

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_prtl
  use mod_data_physics

  use mod_io

  use mod_mapping

  use ppm_module_map_part_util

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'comp_init_state'
  integer                     :: ip, nneigh, jp
  real(mk), dimension(__DIM)  :: xip
  real(mk), dimension(__DIM)  :: dx
  real(mk)                    :: dij
  real(mk)                    :: w, gradW
  real(mk)                    :: rdij
  integer                     :: info

  !-----------------------------
  !< Body of routine
  !-----------------------------
  !< Debug output information
  if (sph_debug .gt. 0) call debugmessage(trim(myname),'START')

  if (symmetry) then
#define __SYMMETRY
#include "inc/sum_kernel.inc"
#undef __SYMMETRY
  else
#include "inc/sum_kernel.inc"
  endif

  !< ip-ip
  do ip = 1, Npart
    drhodt(1,ip) = drhodt(1,ip) + w_zero
  enddo

  !< Put ghost contribution
  call map_ghost_put(myname,realvec1=drhodt)
  call ppm_map_part_load(info) ! After a put load mapping !!!

  !< Volume of particle
  do ip = 1, Npart
    pdata(volidx,ip) = 1._mk/drhodt(1,ip)
  enddo

  !< Debug output information
  if (sph_debug .gt. 0) call debugmessage(trim(myname),'END')

end subroutine comp_init_state
