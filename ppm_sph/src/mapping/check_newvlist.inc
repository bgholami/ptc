!---------------------------------------------------------------
!< Subroutine that checks if a new vlist is required
!---------------------------------------------------------------
subroutine check_newvlist(lnew)

  use mod_data_global
  use mod_data_sph
  use mod_data_prtl
  use mod_data_ctrl
  use mod_data_physics

  use ppm_module_data

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif


  !-----------------------------
  !< Arguments
  !-----------------------------
  logical, intent(out)    :: lnew

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'check_newvlist'
  !integer                     :: info
  real(mk)                    :: dist
  integer                     :: i
  integer                     :: counter
  integer                     :: rcvbuf, ierr

  !----------------------------
  !< Initialize
  !----------------------------
  lnew = .false.

  !-----------------------------
  !< Check if a particle moved further than verlet_skin
  !-----------------------------
  counter = 0

  do i = 1, Npart
    dist = dx_vlist(1,i)**2._mk + dx_vlist(2,i)**2._mk
#if __DIM==3
    dist = dist + dx_vlist(3,i)**2._mk
#endif
    !< Particle moved to far -> new vlist required (0.5*skin)**2
    if (dist .ge. verlet_skin2 - myeps) then
      counter = 1
      exit
    endif
  enddo

  !----------------------------------------
  !< Check if symmetry-bc -> every step new vlist
  !< BUG in ppm: TODO: ghosts of symm-bc not updated with push, only by
  !< ghost_get....
  !-------------------------------------------
!  do i = 1, 2*ndim
!    if (bcdef_ppm(i) .eq. ppm_param_bcdef_symmetry) counter = 1
!  enddo

!  do i = 1, Npart
!    dist = (xp(1,i)-xp_vlist(1,i))**2._mk + (xp(2,i)-xp_vlist(2,i))**2._mk
!#if __DIM==3
!    dist = dist + (xp(3,i)-xp_vlist(3,i))**2._mk
!#endif
!    !< Particle moved to far -> new vlist required (0.5*skin)**2
!    if (dist .ge. verlet_skin2 - myeps) then
!      counter = 1
!      exit
!    endif
!  enddo

#ifdef __MPI
    rcvbuf = 0
    call MPI_ALLREDUCE(counter,rcvbuf,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    counter = rcvbuf
#endif

  if (counter .gt. 0) lnew = .true.

  !if (mod(current_step,1) .eq. 0) lnew = .true.

  !< Always new list
  !lnew = .true.

  return

end subroutine check_newvlist
