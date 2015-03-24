!--------------------------------------------
!< Subroutine writing some vlist results
!--------------------------------------------
subroutine write_vlistresult

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_prtl

  use mod_util

  use ppm_module_io
  use ppm_module_data

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'write_vlistresult'
  integer                     :: info
  character(maxchar)          :: cbuf
  integer                     :: lenchar
  real(mk), dimension(:,:), allocatable :: write_vec
  integer                               :: lda_write_vec
  integer                     :: i, ii
  character(maxchar)          :: fmt_write_vec

  !-----------------------------
  !< Write some output of vlist
  !-----------------------------
  !< Debug output information
  if (sph_debug .gt. 0) call debugmessage(trim(myname),'START')

  !< Filename
  write(cbuf,'(A,I9.9,A)') '/vlistdetails_',current_step,'.debug'
  lenchar = len_trim(thispath)+len_trim(cbuf)
  call check_char_length(lenchar,trim(myname),'filename')
  write(cbuf,'(A)') trim(thispath)//trim(cbuf)

  !< Open file
  call ppm_io_unused_unit(unit_debug,info)
  call ppm_io_open(unit_debug,trim(cbuf),ppm_param_io_write,ppm_param_io_replace,&
  ppm_param_io_ascii,ppm_param_io_centralized,info)
  if (info .ne. 0) call abortmessage(trim(myname),'fileopen failed!')


  !< Write first line
  if (ndim .eq. 2) then
    write(cbuf,'(A)') '# 1:x 2:y 3:rank 4:nvlist'
    lda_write_vec = 4
  else
    write(cbuf,'(A)') '# 1:x 2:y 3:z 4:rank 5:nvlist'
    lda_write_vec = 5
  endif
  call ppm_io(IUNIT=unit_debug,ADATA=cbuf,ACTN=ppm_param_io_write,&
  DIST=ppm_param_io_root,STAT=info)
  if (info .ne. 0) call abortmessage(trim(myname),'write first line failed')

  !< Create 2d array with information
  allocate(write_vec(lda_write_vec,Npart))
  do i = 1, Npart
    do ii = 1, ndim
      write_vec(ii,i) = xp(ii,i)
    enddo
    write_vec(ndim + 1,i) = real(rank,mk)
    write_vec(ndim + 2,i) = real(verletlist(1)%nvlist(i),mk)
  enddo

  !< Format specifier
  write(fmt_write_vec,'(A)') '('
  do i = 1, ndim
    write(fmt_write_vec,'(A)') trim(fmt_write_vec)//'F10.5,1x,'
  enddo
  write(fmt_write_vec,'(A)') trim(fmt_write_vec)//'F6.0,1x,F5.0)'


  !< Write result
  call ppm_io(IUNIT=unit_debug,ADATA=write_vec,&
  ACTN=ppm_param_io_write,DIST=ppm_param_io_concat,&
  IOFMT=trim(fmt_write_vec),STAT=info)
  if (info .ne. 0) call abortmessage(trim(myname),'write vlistresult failed')

  !< Deallocate write_vec
  deallocate(write_vec)

  !< Close file
  call ppm_io_close(unit_debug,info)
  if (info .ne. 0) call abortmessage(trim(myname),'could not close file')

  !< Debug output information
  if (sph_debug .gt. 0) call debugmessage(trim(myname),'END')

  return


end subroutine write_vlistresult
