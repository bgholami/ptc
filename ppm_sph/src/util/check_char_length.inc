subroutine check_char_length(final_length,caller,varname)

  use mod_data_global

  implicit none

  !-----------------------------
  !< Arguments
  !-----------------------------
  integer, intent(in)        :: final_length
  character(LEN=*), intent(in)     :: caller
  character(LEN=*), intent(in)     :: varname

  !-----------------------------
  !< Local variables
  !-----------------------------
  integer                     :: info

  !-----------------------------
  !< Body
  !-----------------------------
  if (rank .eq. 0) then
    if (final_length .gt. maxchar) then
      write(*,*) "Length of '",trim(varname),"' in routine '",trim(caller),"' is going to be too long"
      call finalize(0,info)
      stop
    endif
  endif

end subroutine check_char_length
