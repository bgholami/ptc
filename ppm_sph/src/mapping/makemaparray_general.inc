!--------------------------------
!< Subroutine mapping everything on the topology
!< Input: global or partial
!--------------------------------
subroutine makemaparray_general(realmaparray_loc,intmaparray_loc,info)

  use mod_data_global
  use mod_data_prtl
  use mod_data_ctrl

  use mod_util
  use mod_io

  implicit none

  !-----------------------------
  !< Arguments
  !-----------------------------
  type(t_realarray), dimension(:), pointer  :: realmaparray_loc
  type(t_intarray), dimension(:), pointer   :: intmaparray_loc
  integer, intent(out)                      :: info

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'makemaparray'
  !character(len=maxchar)      :: cbuf

  !< Everything will go well ;-)
  info = 0

  !< Reset array list
  call reset_ptr_to_list(realmaparray_loc)
  call reset_ptr_to_list(intmaparray_loc)

  !< Int arrays
  call add_ptr_to_list(intmaparray_loc,ap)

  !< Real arrays
  call add_ptr_to_list(realmaparray_loc,pdata)
  call add_ptr_to_list(realmaparray_loc,dxpdt)

end subroutine makemaparray_general
