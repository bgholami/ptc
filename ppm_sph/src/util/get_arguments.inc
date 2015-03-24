!------------------------------------------
!< Routine to get the command line arguments
!------------------------------------------
subroutine get_arguments
  use mod_data_global, only: maxchar, commandname, prgname, rank, ctrlfile

  implicit none


  !----------------------------
  !< Local variables
  !----------------------------
  integer            :: ilenprg
  integer            :: info
  integer            :: nargc


  !----------------------------
  !< Get entire command line
  !----------------------------
  call get_command(commandname,info)
  if (info .lt. 1) then
    if (rank .eq. 0) write(*,*) __FILE__,__LINE__,':get_command failed'
    call finalize(1,info)
  endif

  !----------------------------
  ! Get program name
  !----------------------------
  call getarg(0,prgname)
  ilenprg = len_trim(prgname)

  !---------------------------
  !< Read arguments
  !---------------------------
  nargc = command_argument_count()
  if (nargc .lt. 1) then
    if (rank .eq. 0) write(*,*) "No inputfile passed to ",trim(prgname),"!"
    call finalize(0,info)
    stop
  endif
  call getarg(1,ctrlfile)

end subroutine get_arguments
