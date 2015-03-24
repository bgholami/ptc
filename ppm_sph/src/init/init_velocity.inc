!--------------------------------
!< Subroutine to define initial velocity
!--------------------------------
subroutine init_velocity

  use mod_data_global
  use mod_data_physics
  use mod_data_prtl
  use mod_data_ctrl
  use mod_data_sph


  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  integer                     :: ip


  !-----------------------------
  !< 
  !-----------------------------
  do ip = 1, Npart
    dxpdt(1:NDXPDT,ip) = 0._mk
  enddo

end subroutine init_velocity
