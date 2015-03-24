!--------------------------------
!< Subroutine that shifts velocity of particles
!--------------------------------
subroutine shift_velocity(dt_scheme)

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_prtl
  use mod_data_physics

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  real(mk), intent(in)    :: dt_scheme

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'shift_velocity'
  integer                     :: ip

  !-----------------------------
  !< Shift velocity
  !-----------------------------
  do ip = 1, Npart
    dxpdt(1,ip) = dxpdt(1,ip) + dt_scheme*d2xpdt2(1,ip)
    dxpdt(2,ip) = dxpdt(2,ip) + dt_scheme*d2xpdt2(2,ip)
#if __DIM==3
    dxpdt(3,ip) = dxpdt(3,ip) + dt_scheme*d2xpdt2(3,ip)
#endif
  enddo

end subroutine shift_velocity
