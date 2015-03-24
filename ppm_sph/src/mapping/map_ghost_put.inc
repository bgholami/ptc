!--------------------------------
!< Subroutine to put ghost contributions to procs
!--------------------------------
subroutine map_ghost_put(caller,intvec1,realvec1,realvec2,realvec3,realvec4)

  use mod_data_global
  use mod_data_ctrl
  use mod_data_prtl
  use mod_data_sph
  use mod_data_physics

  use mod_io
  use mod_util

  use ppm_module_map_part_ghost
  use ppm_module_map_part
  use ppm_module_write
  use ppm_module_data

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  character(maxchar)            :: caller
  integer , optional, dimension(:,:), pointer :: intvec1
  real(mk), optional, dimension(:,:), pointer :: realvec1
  real(mk), optional, dimension(:,:), pointer :: realvec2
  real(mk), optional, dimension(:,:), pointer :: realvec3
  real(mk), optional, dimension(:,:), pointer :: realvec4

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'map_ghost_put'
  integer                     :: info
  character(maxchar)          :: cbuf

  integer                     :: ldaintvec1  = 0
  integer                     :: ldarealvec1 = 0
  integer                     :: ldarealvec2 = 0
  integer                     :: ldarealvec3 = 0
  integer                     :: ldarealvec4 = 0

  integer                     :: ip
  integer                     :: j

  !--------------------------------
  !< Debug output
  !---------------------------------
  if (sph_debug .gt. 0) then
    call logmessage(trim(myname),'Put ghost particles starting...')
    write(cbuf,'(A,I2,A,I6,A,I6,A,I2,1x,I6)') 'VORHER: rank: ',rank,&
    ' Npart: ',Npart,' Mpart: ',Mpart,' size(xp): ',size(xp,1),size(xp,2)
    call ppm_write(rank,trim(myname),trim(cbuf),info)
  endif

  !--------------------------------
  !< Check if symmetric-interactions -> only then use put
  !--------------------------------
  if (.not. symmetry) return
  !call abortmessage(myname,'are you sure you want to run in symmetry=.true. mode ?')

  !--------------------------------
  !< Defaults
  !--------------------------------
  if (present(intvec1))  ldaintvec1  = size(intvec1,1)
  if (present(realvec1)) ldarealvec1 = size(realvec1,1)
  if (present(realvec2)) ldarealvec2 = size(realvec2,1)
  if (present(realvec3)) ldarealvec3 = size(realvec3,1)
  if (present(realvec4)) ldarealvec4 = size(realvec4,1)

  !----------------------------------------------------
  !< Modify ghost values for symmetric ghosts
  !< Maybe obsolete with new ppm library 1.2.2 -> ask omar...
  !----------------------------------------------------
  do ip = Npart + 1, Mpart
    do j = 1, ndim
      !<-------------------------------
      !< Upper boundary
      !<-------------------------------
      if (xp(j,ip) .gt. max_compbox(j)) then
        if (bcdef_ppm(2*j) .eq. ppm_param_bcdef_symmetry .or. &
        bcdef_ppm(2*j) .eq. ppm_param_bcdef_neumann .or. &
        bcdef_ppm(2*j) .eq. ppm_param_bcdef_dirichlet) then
          if (present(intvec1)) intvec1(:,ip) = 0._mk
          if (present(realvec1)) realvec1(:,ip) = 0._mk
          if (present(realvec2)) realvec2(:,ip) = 0._mk
          if (present(realvec3)) realvec3(:,ip) = 0._mk
          if (present(realvec4)) realvec4(:,ip) = 0._mk
        endif
      !<-------------------------------
      !< Lower boundary
      !<-------------------------------
      elseif (xp(j,ip) .lt. min_compbox(j)) then
        if (bcdef_ppm(2*j-1) .eq. ppm_param_bcdef_symmetry .or. &
        bcdef_ppm(2*j-1) .eq. ppm_param_bcdef_neumann .or. &
        bcdef_ppm(2*j-1) .eq. ppm_param_bcdef_dirichlet) then
          if (present(intvec1)) intvec1(:,ip) = 0._mk
          if (present(realvec1)) realvec1(:,ip) = 0._mk
          if (present(realvec2)) realvec2(:,ip) = 0._mk
          if (present(realvec3)) realvec3(:,ip) = 0._mk
          if (present(realvec4)) realvec4(:,ip) = 0._mk
        endif
      endif
    enddo
  enddo

  !--------------------------------
  !< PUT
  !--------------------------------
  call ppm_map_part_ghost_put(topo_id,info)
  if (info .ne. 0) call abortmessage(trim(myname),'part_ghost_put failed')

  !--------------------------------
  !< PUSH
  !--------------------------------
  if (present(intvec1)) call ppm_map_part_push(intvec1,ldaintvec1,Npart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': push intvec1 failed')
  if (present(realvec1)) call ppm_map_part_push(realvec1,ldarealvec1,Npart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': push realvec1 failed')
  if (present(realvec2)) call ppm_map_part_push(realvec2,ldarealvec2,Npart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': push realvec2 failed')
  if (present(realvec3)) call ppm_map_part_push(realvec3,ldarealvec3,Npart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': push realvec3 failed')
  if (present(realvec4)) call ppm_map_part_push(realvec4,ldarealvec4,Npart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': push realvec4 failed')

  !--------------------------------
  !< SEND
  !--------------------------------
  call ppm_map_part_send(Npart,Mpart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'part_send failed')

  !--------------------------------
  !< POP
  !--------------------------------
  if (present(realvec4)) call ppm_map_part_ghost_pop(realvec4,ldarealvec4,Npart,Mpart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': pop realvec4 failed')
  if (present(realvec3)) call ppm_map_part_ghost_pop(realvec3,ldarealvec3,Npart,Mpart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': pop realvec3 failed')
  if (present(realvec2)) call ppm_map_part_ghost_pop(realvec2,ldarealvec2,Npart,Mpart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': pop realvec2 failed')
  if (present(realvec1)) call ppm_map_part_ghost_pop(realvec1,ldarealvec1,Npart,Mpart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': pop realvec1 failed')
  if (present(intvec1)) call ppm_map_part_ghost_pop(intvec1,ldaintvec1,Npart,Mpart,info)
  if (info .ne. 0) call abortmessage(trim(myname),'called from '//trim(caller)//': pop intvec1 failed')

  !--------------------------------
  ! Debug output
  !--------------------------------
  if (sph_debug .gt. 0) then
    call logmessage(trim(myname),'Put ghost particles with data done')
    write(cbuf,'(A,I2,A,I6,A,I6,A,I2,1x,I6)') 'NACHHER: rank: ',rank,&
    ' Npart: ',Npart,' Mpart: ',Mpart,' size(xp): ',size(xp,1),size(xp,2)
    call ppm_write(rank,trim(myname),trim(cbuf),info)
  endif

end subroutine map_ghost_put
