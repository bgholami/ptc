!--------------------------------
!< Subroutine that gives min, max, mean and var of numbers on different procs
!--------------------------------
subroutine minmaxmeanvar_r8(realin,outvec,info)

  use mod_data_global

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  !-----------------------------
  !< Arguments
  !-----------------------------
  real(mk), intent(in)  :: realin
  real(mk), dimension(n_procbal), intent(out) :: outvec
  integer                             :: info

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'minmaxmeanvar'
  integer                     :: ierr
  real(mk)                    :: sumvar
  integer                     :: j

  !-----------------------------
  !< Get min, max, mean and var
  !-----------------------------
#ifdef __MPI
  gather_nproc = 0.
  call MPI_GATHER(realin,1,MPPREC,gather_nproc,1,MPPREC,0,comm,ierr)

  outvec(minidx) = minval(gather_nproc)
  outvec(maxidx) = maxval(gather_nproc)
  outvec(meanidx) = sum(gather_nproc)/Nproc

  sumvar = 0._mk
  do j = 1, Nproc
    sumvar = sumvar + (outvec(meanidx) - gather_nproc(j))**2
  enddo
  outvec(varidx) = sumvar/Nproc
#else
  outvec(minidx) = realin
  outvec(maxidx) = realin
  outvec(meanidx) = realin
  outvec(varidx)  = 0._mk
#endif

end subroutine minmaxmeanvar_r8
