!--------------------------------
!< Subroutine to check if abort calculation due to ABORTFILE
!--------------------------------
subroutine check_abortfile(lexist)

  use mod_data_global

  use ppm_module_write

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  !-----------------------------
  !< Arguments
  !-----------------------------
  logical, intent(out)  :: lexist

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'check_abortfile'
  integer                     :: info, ios
  character(maxchar)          :: cbuf

  !-----------------------------
  !< Check if file exists
  !-----------------------------
  lexist = .false.
  ios = 0
  if (rank .eq. 0) inquire(FILE=abortfile,EXIST=lexist,IOSTAT=ios)

#ifdef __MPI
  !-------------------------------------------------------------------------
  !  Broad cast the result
  !-------------------------------------------------------------------------
  CALL MPI_BCast(ios,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(lexist,1,MPI_LOGICAL,0,comm,info)
#endif

  if (ios .ne. 0) then
    lexist = .true.
    write(cbuf,'(A)') 'Problem with inquire FILE=abortfile (ios .ne. 0) !'
    call ppm_write(rank,trim(myname),trim(cbuf),info)
    call ppm_write(rank,trim(myname),trim(cbuf),info,unit_log)
    return
  endif

  if (lexist) then
    write(cbuf,'(A)') 'ABORTFILE found, terminating calculation....'
    if (rank .eq. 0) call ppm_write(rank,trim(myname),trim(cbuf),info)
    if (rank .eq. 0) call ppm_write(rank,trim(myname),trim(cbuf),info,unit_log)
  endif


  return

end subroutine check_abortfile
