subroutine finalize(errstat,info)

	use mod_data_global
	use mod_data_ctrl
	use mod_data_physics
	use mod_data_prtl
	use mod_data_sph

	use ppm_module_finalize

	implicit none

#ifdef __MPI
	include 'mpif.h'
#endif

	!---------------------------------
	!< Arguments
	!---------------------------------
	integer, intent(in)     :: errstat
	integer, intent(  out)  :: info

  !---------------------------------
  !< Local variables
  !---------------------------------
  logical                 :: opened
  integer                 :: i

  !--------------------------
  !< Deallocate arrays
  !--------------------------
  if (allocated(pForce_total)) deallocate(pForce_total)

  do i = 1, NVERLETLIST
    if (associated(verletlist(i)%vlist)) deallocate(verletlist(i)%vlist)
    if (associated(verletlist(i)%nvlist)) deallocate(verletlist(i)%nvlist)
    if (allocated(verletlist(i)%pidx)) deallocate(verletlist(i)%pidx)
  enddo
  if (allocated(verletlist)) deallocate(verletlist)

  if (allocated(gather_nproc)) deallocate(gather_nproc)

  !--------------------------
  !< Free pointer
  !--------------------------
  if (associated(Ntotal)) deallocate(Ntotal)

  if (associated(xp)) deallocate(xp)
  if (associated(dxpdt)) deallocate(dxpdt)
  if (associated(pdata)) deallocate(pdata)
  if (associated(ap)) deallocate(ap)
  if (associated(d2xpdt2)) deallocate(d2xpdt2)
  if (associated(drhodt)) deallocate(drhodt)
  if (associated(dx_vlist)) deallocate(dx_vlist)

  if (associated(cost)) deallocate(cost)
  if (associated(min_sub)) deallocate(min_sub)
  if (associated(max_sub)) deallocate(max_sub)
  if (associated(sub2proc)) deallocate(sub2proc)
  if (associated(pcost)) deallocate(pcost)
  if (associated(ident)) deallocate(ident)

  if (associated(realmaparray)) deallocate(realmaparray)
  if (associated(intmaparray)) deallocate(intmaparray)

 	!---------------------------------
	!< Finalize ppm
	!---------------------------------
	call ppm_finalize(info)

  !---------------------------------
  !< Close files
  !---------------------------------
  inquire(unit=unit_log,opened=opened)
  if (opened) close(unit_log,iostat=info)
  inquire(unit=unit_proclog,opened=opened)
  if (opened) close(unit_proclog,iostat=info)
  inquire(unit=unit_writeascii,opened=opened)
  if (opened) close(unit_writeascii,iostat=info)
  inquire(unit=unit_debug,opened=opened)
  if (opened) close(unit_debug,iostat=info)
  inquire(unit=unit_ppmlog,opened=opened)
  if (opened) close(unit_ppmlog,iostat=info)

#ifdef __MPI
	!---------------------------------
	! Finalize mpi
	!---------------------------------
	if (errstat .EQ. 0) THEN
	  !---------------------------------------------------------------------
	  !  Finalize MPI
	  !---------------------------------------------------------------------
	  call MPI_Finalize(info)
	  IF (info .NE. 0) THEN
	    WRITE(*,*) 'FAILED TO FINALIZE MPI. BAILING OUT!'
	    STOP
	  ENDIF
	else
	  !---------------------------------------------------------------------
	  !  Abort MPI
	  !---------------------------------------------------------------------
	  CALL MPI_Abort(comm,errstat,info)
	  IF (info .NE. 0) THEN
	    WRITE(*,*) 'FAILED TO ABORT MPI!'
	  ENDIF
    STOP
	endif
#endif

end subroutine finalize

