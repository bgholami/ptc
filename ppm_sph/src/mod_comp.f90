!--------------------------------
!< Module that contains all computations
!--------------------------------------
module mod_comp

  implicit none

contains

#include "src/comp/press_EOS.inc"
#include "src/comp/comp_forces.inc"
#include "src/comp/comp_pp_interactions.inc"
#include "src/comp/comp_particleloop_1.inc"
#include "src/comp/comp_init_state.inc"
#include "src/comp/comp_phasevlist.inc"
!#include "src/comp/comp_stresstensor.inc"
!#include "src/comp/comp_interfaceloop_1.inc"
!#include "src/comp/comp_interfaceloop_2.inc"
#include "src/comp/comp_boundarycondition.inc"
#include "src/comp/comp_constants.inc"

!----------------------------------------
!< Get density from pressure
!----------------------------------------
!real(mk) function density_from_pressure(press,ip)
!
!  use mod_data_global
!  use mod_data_physics
!  use mod_data_prtl
!
!  real(mk), intent(in) :: press
!  integer,  intent(in) :: ip
!
!  integer :: pip
!
!  pip = ap(wallphaseidx,ip)
!
!  !< Check that result exists
!  if (((press - phase(pip)%press_b)/phase(pip)%press_ref + 1._mk) .lt. 0._mk) then
!    density_from_pressure = phase(pip)%density
!  else
!    density_from_pressure = ( (press - phase(pip)%press_b)/phase(pip)%press_ref + 1._mk)**&
!    (1._mk / phase(pip)%gamma)*phase(pip)%density
!  endif
!
!  !< Boundary particles at fringes might have zero density as only one fluid
!  !< particle at distance r_c interacts due to weighting wij
!  !< anyway no influence, but code was somehow unstable otherwise..
!  if (density_from_pressure .lt. myeps) density_from_pressure = phase(pip)%density
!
!end function density_from_pressure

end module mod_comp
