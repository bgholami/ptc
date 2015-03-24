!--------------------------------
!< Subroutine to kick particle positions
!--------------------------------
subroutine shift_position(dt_scheme)

  use mod_data_global
  use mod_data_prtl
  use mod_data_sph
  use mod_data_ctrl

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  real(mk), intent(in)  :: dt_scheme

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'shift_position'
  integer                     :: ip

  !-----------------------------
  !< Kick position
  !-----------------------------
  do ip = 1, Npart
      xp(1,ip) = xp(1,ip) + dt_scheme*dxpdt(1,ip)
      xp(2,ip) = xp(2,ip) + dt_scheme*dxpdt(2,ip)
      dx_vlist(1,ip) = dx_vlist(1,ip) + dt_scheme*dxpdt(1,ip)
      dx_vlist(2,ip) = dx_vlist(2,ip) + dt_scheme*dxpdt(2,ip)
#if __DIM==3
      xp(3,ip) = xp(3,ip) + dt_scheme*dxpdt(3,ip)
      dx_vlist(3,ip) = dx_vlist(3,ip) + dt_scheme*dxpdt(3,ip)
#endif
  enddo

end subroutine shift_position
