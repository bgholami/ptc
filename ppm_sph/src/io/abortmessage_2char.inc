!-----------------------------------
!< Subroutine to abort safely
!-----------------------------------
subroutine abortmessage_2char(routine,mesg)
  use mod_data_global
  use mod_util, only : finalize
  use ppm_module_write


  implicit none

  !-----------------------------
  !< Arguments variables
  !-----------------------------
  character(len=*), intent(in   ) :: routine, mesg

  !-----------------------------
  !< Local variables
  !-----------------------------
  integer                 :: info

  !-----------------------------
  !< Body
  !-----------------------------
  if (rank .eq. 0) then
    call ppm_write(rank,trim(routine),trim(mesg),info) !< stdout
    call ppm_write(rank,trim(routine),trim(mesg),info,unit_log) !< log output
  endif

  call finalize(0,info)

  stop

end subroutine abortmessage_2char
