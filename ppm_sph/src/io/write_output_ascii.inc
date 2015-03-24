!--------------------------------
!< Writing output using ascii
!--------------------------------
subroutine write_output_ascii(fname,prtl_1,prtl_N)

  use mod_data_global
  use mod_data_prtl
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_physics

  use mod_util
  use ppm_module_io
  use ppm_module_data

  implicit none

  !-----------------------------
  !< Arguments
  !-----------------------------
  character(len=*), intent(in)  :: fname
  integer, intent(in)           :: prtl_1
  integer, intent(in)           :: prtl_N

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'write_output_ascii'
  integer                     :: info
  character(len=maxchar)      :: cbuf
  character(len=maxchar)      :: filename, iofmt
  real(mk)                    :: rfac
  real(mk), dimension(:), allocatable :: outvec

  integer                                       :: i, ii

  !-----------------------------
  !< Construct filename
  !-----------------------------
  call check_char_length(len_trim(trim(fname)//'-000000.out'),'write_output_ascii','filename')
  write(filename,'(A,I6.6,A)') trim(fname)//'-',current_step,'.out'

  !-------------------------------------------------------------------------
  !<  OPEN/APPEND FILE
  !-------------------------------------------------------------------------
  call ppm_io_unused_unit(unit_writeascii,info)
  call ppm_io_open(IUNIT=unit_writeascii,FILENAME=filename,&
  ACTN=ppm_param_io_write,POSN=ppm_param_io_replace,&
  PFMT=ppm_param_io_ascii,PMODE=ppm_param_io_centralized,INFO=info)
  if (info .ne. 0) call abortmessage(myname,'failed to open '//trim(filename))

  !---------------------------------
  !< Write first line
  !---------------------------------
  if (choice_output_nondim) then
    write(cbuf,'(A,F17.7)') '#Snapshot t = ',current_time/length_ref*vel_ref
  else
    write(cbuf,'(A,F17.7)') '#Snapshot t = ',current_time
  endif

  cbuf = TRIM(cbuf)
  CALL ppm_io(unit_writeascii,cbuf,ppm_param_io_write,&
  ppm_param_io_root,STAT=info)
  IF (info .NE. 0) CALL abortmessage(myname,'Writing first line failed')

  ldaoutvec_ascii = ndim
  write(cbuf,'(A)') '# 01:x 02:y'
  write(iofmt,'(A)') '(F12.6,1x,F12.6'
  if (ndim .eq. 3) then
    write(cbuf,'(A)') trim(cbuf)//' 03:z'
    write(iofmt,'(A)') trim(iofmt)//',1x,F12.6'
  endif
  if (output_velocity) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':Ux'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+2,':Uy'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    if (ndim .eq. 3) then
      write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+3,':Uz'
      write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    endif
    ldaoutvec_ascii = ldaoutvec_ascii + ndim
  endif
  if (output_density) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':rho'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_pressure) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':p'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_volume) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':V'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_mass) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':m'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_h) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':h'
    write(iofmt,'(A)') trim(iofmt)//',1x,E13.6'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_phaseidx) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':Phase'
    write(iofmt,'(A)') trim(iofmt)//',1x,F2.0'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_id) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':ID'
    write(iofmt,'(A)') trim(iofmt)//',1x,F2.0'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  if (output_rank) then
    write(cbuf,'(A,I2.2,A)') trim(cbuf)//' ',ldaoutvec_ascii+1,':rank'
    write(iofmt,'(A)') trim(iofmt)//',1x,F5.0'
    ldaoutvec_ascii = ldaoutvec_ascii + 1
  endif
  write(iofmt,'(A)') trim(iofmt)//')'

  call ppm_io(IUNIT=unit_writeascii,ADATA=cbuf,ACTN=ppm_param_io_write,&
  DIST=ppm_param_io_root,IOFMT='(A)',STAT=info)

  !-----------------------------
  !< Allocate outvec
  !-----------------------------
  allocate(outvec(ldaoutvec_ascii))
  outvec= 0._mk

  !-----------------------------------
  !< Store everything in outvec
  !-----------------------------------
  ldaoutvec_ascii = 0 !< here used as counter for outvec(ldaoutvec_ascii,:)
  if (choice_output_nondim .eqv. .true.) then
#define __NONDIMOUTPUT
#include "inc/output_ascii.inc"
#undef __NONDIMOUTPUT
  else
#include "inc/output_ascii.inc"
  endif

  !-----------------------------
  !< Deallocate outvec
  !-----------------------------
  deallocate(outvec)

  !---------------------------------
  !< Close file
  !---------------------------------
  call ppm_io_close(unit_writeascii,info)
  if (info .ne. 0) call abortmessage(myname,'closing '//trim(filename)//' failed...')

  !-----------------------------
  !< Confirm successful output
  !-----------------------------
  write(cbuf,fmt_char_char) 'Output written to file: ',trim(filename)
  call logmessage('write_output_ascii',trim(cbuf))


end subroutine write_output_ascii
