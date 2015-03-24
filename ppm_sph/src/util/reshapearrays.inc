!--------------------------------
!< Subroutine that resets all arrays
!--------------------------------
subroutine reshapearrays

  use mod_data_global
  use mod_data_prtl
  use mod_data_physics, only : N_PHASES
  use mod_data_sph, only : Mpart
  use mod_data_ctrl, only : ndim

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'resetarrays'
  integer                     :: info
  integer                     :: newsize

  !-----------------------------
  !< Reset arrays
  !-----------------------------
  if (size(pForce_total,2) .lt. Mpart) then
    !newsize = 2*size(pforce_total,2)
    newsize = Mpart+10

    deallocate(d2xpdt2,STAT=info)
    deallocate(drhodt)
    deallocate(pForce_total)

    allocate(d2xpdt2(ndim,newsize))
    allocate(drhodt(NDRHODT,newsize))
    allocate(pForce_total(ndim,newsize))

  endif


end subroutine reshapearrays
