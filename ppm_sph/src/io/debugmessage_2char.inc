!-----------------------------------
!< Subroutine to abort safely
!-----------------------------------
subroutine debugmessage_2char(routine,mesg)

  use mod_data_global

  use ppm_module_write

  implicit none

  !-----------------------------
  !< Arguments variables
  !-----------------------------
  character(len=*), intent(in   ) :: routine, mesg

  character(len=maxchar)  :: text

  !-----------------------------
  !< Local variables
  !-----------------------------
  integer                 :: info

  !-----------------------------
  !< Body
  !-----------------------------

  text = trim(mesg)//"  (DEBUGMODE)"

  call ppm_write(rank,trim(routine),trim(mesg),info) !< stdout
  call ppm_write(rank,trim(routine),trim(mesg),unit_proclog) !< proclog-file

end subroutine debugmessage_2char
