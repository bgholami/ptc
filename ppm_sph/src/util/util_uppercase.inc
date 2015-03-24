!--------------------------------
!< Subroutine to convert char-variable to uppercase
!--------------------------------
subroutine util_uppercase(string,ilen,info)

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  CHARACTER(LEN=*), INTENT(INOUT) :: string
  INTEGER         , INTENT(IN   ) :: ilen
  INTEGER         , INTENT(  OUT) :: info

  !-------------------------------------------------------------------------
  !<  Local variables
  !-------------------------------------------------------------------------
  INTEGER          :: i,j
  INTEGER          :: i1,i2,i3,iadd

  !-----------------------------
  !< Body of routine
  !-----------------------------
   info = 0

  !-------------------------------------------------------------------------
  !  Convert to upper case
  !-------------------------------------------------------------------------
  i1   = IACHAR('a') - 1
  i2   = IACHAR('z') + 1
  i3   = IACHAR('A')
  iadd = i3 - i1 - 1
  DO i=1,ilen
     j = IACHAR(string(i:i))
     IF (j.GT.i1.AND.j.LT.i2) THEN
        string(i:i) = CHAR(j+iadd)
     ENDIF
  ENDDO

  return

end subroutine util_uppercase
