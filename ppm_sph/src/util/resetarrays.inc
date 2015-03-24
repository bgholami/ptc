!--------------------------------
!< Subroutine that resets all arrays
!--------------------------------
subroutine resetarrays

  use mod_data_global
  use mod_data_prtl
  use mod_data_physics, only : N_PHASES
  use mod_data_ctrl

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'reshapearrays'

  !-----------------------------
  !< Reset arrays
  !-----------------------------
  pForce_total = 0._mk

  d2xpdt2 = 0._mk
  drhodt = 0._mk

end subroutine resetarrays
