!--------------------------------
!< Writing output using h5part in parallel
!--------------------------------
subroutine write_ppmdbg(intstep)

  use mod_data_global
  use mod_data_prtl
  use mod_data_ctrl
  use mod_data_sph

  use mod_util

  use ppm_module_util_dbg


  implicit none


  !-----------------------------
  !< Arguments
  !-----------------------------
  integer, intent(in)           :: intstep

  !-----------------------------
  !< Local variables
  !-----------------------------
  integer                     :: info
  character(len=maxchar)      :: myname = 'write_ppmdbg'

  !-----------------------------
  ! write output
  !-----------------------------
  call ppm_dbg_print(topo_id,sph_bsize,intstep,1,info,xp,Npart,Mpart)
  if (info .ne. 0) call abortmessage(myname,'ppm_dbg_print failed')

  !-----------------------------
  !< Confirm successful output
  !-----------------------------
  call logmessage(myname,'ppm_dbg_print Output written')


end subroutine write_ppmdbg
