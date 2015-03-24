!--------------------------------
!< Subroutine that updates density -> summation form
!--------------------------------
subroutine shift_density(dt_scheme)

  use mod_data_global
  use mod_data_prtl
  use mod_data_sph
  use mod_data_ctrl
  use mod_data_physics

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  real(mk), intent(in)     :: dt_scheme

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'shift_density'
  integer                     :: ip

  !-----------------------------
  !< Density summation form
  !-----------------------------
  do ip = 1, Npart
    !< Vol_i = 1. / sum W_ij
    pdata(volidx,ip) = 1._mk / (drhodt(1,ip) + w_zero)
    !< rho_i = m_i * sum W_ij
    pdata(rhoidx,ip) = pdata(massidx,ip) / pdata(volidx,ip)
  enddo

end subroutine shift_density
