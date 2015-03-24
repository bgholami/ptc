!--------------------------------
!< Subroutine to define default settings / setup folder...
!--------------------------------
subroutine defaults

  use mod_data_global
  use mod_data_ctrl

  use mod_util
  use mod_io

  use ppm_module_write
  use ppm_module_io

  implicit none

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'defaults'
  integer                     :: len_ctrlfile
  integer                     :: i
  integer                     :: info
  character(len=maxchar)      :: cbuf

  !-----------------------------
  !< Get current path thispath
  !-----------------------------
  !< Relative directories
  thispath=''
  !< Absolut path
  !call getcwd(thispath)
  !thispath=trim(thispath)//'/'

  i = 1
  do while (i .le. len_trim(ctrlfile))
    if (ctrlfile(i:i) .eq. ".") then
      len_ctrlfile = i-1
      goto 9
    endif
    i = i + 1
  enddo
  call ppm_write(rank,myname,'ctrlfile without fileextension!',info)
  call finalize(0,info)
  stop

9 continue

  i = len_trim(trim(thispath)) + len_ctrlfile
  call check_char_length(i,'defaults','thispath')

  write(thispath,'(A)') trim(thispath)//ctrlfile(1:len_ctrlfile)

  !-------------------------------------------------------------------------
  !<  Default log filename
  !-------------------------------------------------------------------------
  if (rank .eq. 0) then
    i = len_trim(trim(thispath)) + len_trim('/ppm_sph.log')
    call check_char_length(i,'defaults','logfile')
    write(logfile,'(2A)') TRIM(thispath),'/ppm_sph.log'
  endif

  !-------------------------------------------------------------------------
  !< Conservation file
  !-------------------------------------------------------------------------
  i = len_trim(trim(thispath)) + len_trim('/conservation.txt')
  call check_char_length(i,'defaults','conservationfile')
  WRITE(conservationfile,'(A)') TRIM(thispath)//'/conservation'

  !-------------------------------------------------------------------------
  !< logfile for every proc
  !-------------------------------------------------------------------------
  if (Nproc .lt. 10) then
    WRITE(proclogfile,'(A,A,I1.1,A)') TRIM(thispath),'/proc_',rank,'.log'
  elseif (Nproc .lt. 100) then
    WRITE(proclogfile,'(A,A,I2.2,A)') TRIM(thispath),'/proc_',rank,'.log'
  elseif (Nproc .lt. 1000) then
    WRITE(proclogfile,'(A,A,I3.3,A)') TRIM(thispath),'/proc_',rank,'.log'
  elseif (Nproc .lt. 10000) then
    WRITE(proclogfile,'(A,A,I4.4,A)') TRIM(thispath),'/proc_',rank,'.log'
  elseif (Nproc .lt. 100000) then
    WRITE(proclogfile,'(A,A,I5.5,A)') TRIM(thispath),'/proc_',rank,'.log'
  else
    if (rank .eq. 0) then
      write(*,*) 'Really use more than 100000 procs???'
    endif
    call finalize(0,info)
  endif

  !-------------------------------------------------------------------------
  !<  Abort file
  !-------------------------------------------------------------------------
  WRITE(abortfile,'(2A)') TRIM(thispath),'/ABORTFILE'

  !-------------------------------------------------------------------------
  !< Create folder
  !-------------------------------------------------------------------------
  WRITE(cbuf,'(2A)') 'mkdir -p ',TRIM(thispath)
  CALL system(TRIM(cbuf))
  IF (rank .EQ. 0) THEN
    WRITE(cbuf,'(5A)') 'cp ',TRIM(ctrlfile),' ',TRIM(thispath),'/Used.inputfile'
    CALL system(TRIM(cbuf))

    WRITE(cbuf,'(2A,1x,A)') 'mv -f ',TRIM(logfile),TRIM(logfile)//"old"
    CALL system(TRIM(cbuf))
  ENDIF

  if (sph_debug .gt. 0) then
    WRITE(cbuf,'(2A,1x,A)') 'mv -f ',TRIM(proclogfile),TRIM(proclogfile)//"old"
    CALL system(TRIM(cbuf))
  endif

  !-----------------------------
  !< Outputfilenames
  !-----------------------------
  write(outputfile,'(2A)') trim(thispath),'/particles'
  write(ghostfile,'(2A)') trim(thispath),'/ghosts'

  !-----------------------------
  !< Print names for timings
  !-----------------------------
  wct_name(stepidx)       = 'timestep'
  wct_name(velidx)        = 'shiftvel'
  wct_name(posidx)        = 'shiftpos'
  wct_name(densidx)       = 'density'
  wct_name(propidx)       = 'properties'
  wct_name(pushidx)       = 'push'
  wct_name(forceidx)      = 'force'
  wct_name(vlistidx)      = 'vlist'
  wct_name(bcidx)         = 'bc'
  wct_name(resetidx)      = 'reset'
  wct_name(putidx)        = 'put'
  wct_name(ploop1idx)     = 'ploop1'
  wct_name(checkidx)      = 'checkvlist'
  wct_name(ghostidx)      = 'ghostxp'
  wct_name(bodyforceidx)  = 'bodyforce'

end subroutine defaults
