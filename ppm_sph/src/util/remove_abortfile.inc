!------------------------------------
!< Subroutine to remove the abortfile
!------------------------------------
subroutine remove_abortfile

  use mod_data_global

  use ppm_module_write

  implicit none

  !------------------------------
  !< Arguments
  !------------------------------

  !-------------------------------
  !< Local Variables
  !--------------------------------
  character(maxchar)      :: myname = 'remove_abortfile'
  logical     :: lexist
  character(len=maxchar) :: cbuf
  integer                :: info
  !---------------------------------
  !< Body of routine
  !---------------------------------

  !----------------------------------------------------------------------
  !<  Check if the ABORT file exists
  !----------------------------------------------------------------------
  IF (rank .EQ. 0) THEN
    INQUIRE(FILE=trim(abortfile),EXIST=lexist)
    IF (lexist) then
      WRITE(cbuf,'(2A)') 'rm -f ',TRIM(abortfile)
      CALL system(TRIM(cbuf))
      WRITE(cbuf,'(2A)') 'Deleting ',trim(abortfile)
      if (rank .eq. 0) then
        call ppm_write(rank,trim(myname),trim(cbuf),info)
        call ppm_write(rank,trim(myname),trim(cbuf),info,unit_log)
      endif
    ENDIF
  ENDIF

end subroutine remove_abortfile


