!--------------------------------
!< Subroutine that creates one single particle if on proc
!--------------------------------
#if __DIM==2
subroutine place_particle(xpos,ypos,pID,Np_loc,created)
#else
subroutine place_particle(xpos,ypos,zpos,pID,Np_loc,created)
#endif

  use mod_data_global
  use mod_data_sph
  use mod_data_prtl
  use mod_data_physics
  use mod_data_ctrl
  use mod_io, only: abortmessage

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  real(mk), intent(in   ) :: xpos
  real(mk), intent(in   ) :: ypos
#if __DIM==3
  real(mk), intent(in   ) :: zpos
#endif
  integer , intent(in   ) :: pID
  integer , intent(inout) :: Np_loc
  logical , intent(  out) :: created

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'place_particle'

  !-----------------------------
  !< Check vars
  !-----------------------------
  if (pID .gt. N_PHASES) call abortmessage(myname,'pID > N_PHASES!')
  if (pID .lt. 1) then
    created = .false.
    return
    !call abortmessage(myname,'pID < 1')
  endif

  !-----------------------------
  !< define particle
  !-----------------------------
  Npart = Npart + 1
  created = .true.

  if (mod(Npart,Nproc) .eq. rank) then
    Np_loc = Np_loc + 1
    ! Position
    xp(1,Np_loc) = xpos
    xp(2,Np_loc) = ypos
#if __DIM==3
    xp(3,Np_loc) = zpos
#endif

    ! Phase ID
    ap(phaseidx,Np_loc) = pID

    ! Particle ID
    ap(partIDidx,Np_loc) = Npart
  endif

end subroutine place_particle
