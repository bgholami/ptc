!--------------------------------
!< Subroutine generic mapping particles/ghosts on the topology
!< Input: global or partial or ghost or ghostpushpp
!--------------------------------
subroutine mapping(caller,maptype,realmaparray_loc,intmaparray_loc)
  use mod_data_global, only: maxchar, t_realarray, t_intarray, mk, &
  sph_debug, param_maptype_global, topo_id, param_maptype_partial, param_maptype_getghost,&
  isymm, param_maptype_ghostupdate, param_maptype_ghostupdatenoxp, lastrealmaparray, &
  lastintmaparray, fmt_char_int10, rank, unit_proclog
  use mod_data_prtl, only: xp
  use mod_data_sph, only: npart, cutoff, verlet_skin, mpart
  use mod_data_ctrl, only: ndim

  use mod_util
  use mod_io

  use ppm_module_map
  use ppm_module_map_part_global
  use ppm_module_map_part_partial
  use ppm_module_write
  use ppm_module_data
  use ppm_module_substart
  use ppm_module_substop

  implicit none


  !-----------------------------
  !< Arguments
  !-----------------------------
  character(maxchar)          :: caller
  integer, intent(in)         :: maptype
  type(t_realarray), dimension(:), pointer :: realmaparray_loc
  type(t_intarray), dimension(:), pointer :: intmaparray_loc

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'mapping'
  integer                     :: info
  character(len=maxchar)      :: cbuf
  integer                     :: j
  real(mk)                    :: ctime

  !-----------------------------
  !< Debug
  !-----------------------------
  if (sph_debug .gt. 0) then
    call substart(myname,ctime,info)
  endif

  !-----------------------------
  !< Map everything on topology
  !-----------------------------
  select case (maptype)

  !< GLOBAL MAPPING OF PARTICLES
  case (param_maptype_global)
    CALL ppm_map_part_global(topo_id,xp,Npart,info) ! positions
    if (info .ne. 0) call abortmessage(trim(myname),'map_part_global xp failed')

  !< PARTIAL MAPPING OF PARTICLES
  case (param_maptype_partial)
    CALL ppm_map_part_partial(topo_id,xp,Npart,info) ! positions
    if (info .ne. 0) call abortmessage(trim(myname),'map_part_global xp failed')

  case (param_maptype_getghost)
    call ppm_map_part_ghost_get(topo_id,xp,ndim,Npart,isymm,cutoff+2._mk*verlet_skin,info)
    if (info .ne. 0) call abortmessage(trim(myname),'part_ghost_get failed')
    if (sph_debug .gt. 0) call logmessage(myname,'ppm_map_part_ghost_get done')

  case (param_maptype_ghostupdate)
    call ppm_map_part_push(xp,__DIM,Npart,info,pushpp=.true.)
    if (info .ne. 0) call abortmessage(myname,'called from '//trim(caller)//': push xp failed')

  case (param_maptype_ghostupdatenoxp)
    !< do nothing -> only properties are pushed...

  case default
    call abortmessage(trim(myname),'Wrong maptype when calling mapping')
  end select

  !--------------------------
  !< PUSH
  !--------------------------
  do j = 1, lastrealmaparray
    call ppm_map_part_push(realmaparray_loc(j)%ptr_data,realmaparray_loc(j)%lda,Npart,info)
    if (info .ne. 0) call abortmessage(myname,'called from '//trim(caller)//': push realmaparray_loc failed')
  enddo
  do j = 1, lastintmaparray
    call ppm_map_part_push(intmaparray_loc(j)%ptr_data,intmaparray_loc(j)%lda,Npart,info)
    if (info .ne. 0) call abortmessage(myname,'called from '//trim(caller)//': push intmaparray_loc failed')
  enddo

  !--------------------------
  !< SEND
  !--------------------------
  CALL ppm_map_part_send(Npart,Mpart,info) ! send
  if (info .ne. 0) call abortmessage(trim(myname),'send failed')

  !--------------------------
  !< POP
  !--------------------------
  do j = lastintmaparray, 1, -1
    call ppm_map_part_pop(intmaparray_loc(j)%ptr_data,intmaparray_loc(j)%lda,Npart,Mpart,info)
    if (info .ne. 0) call abortmessage(myname,'called from '//trim(caller)//': pop intmaparray_loc failed')
  enddo
  do j = lastrealmaparray, 1, -1
    call ppm_map_part_pop(realmaparray_loc(j)%ptr_data,realmaparray_loc(j)%lda,Npart,Mpart,info)
    if (info .ne. 0) call abortmessage(myname,'called from '//trim(caller)//': pop realmaparray_loc failed')
  enddo

  if (maptype .ne. param_maptype_ghostupdatenoxp) then
    CALL ppm_map_part_pop(xp,ndim,Npart,Mpart,info) ! positions
    if (info .ne. 0) call abortmessage(trim(myname),'pop xp failed')
  endif

  !< Update Npart only for global or partial mapping...
  if (maptype .eq. param_maptype_global .or. maptype .eq. param_maptype_partial) then
    Npart = Mpart
  endif

  if (sph_debug .gt. 0) then
      !< Some details of map
    call logmessage(trim(myname),'done mapping particles on topo')
    write(cbuf,fmt_char_int10) 'Local number of particles is now: ',Npart
    call logmessage(trim(myname),trim(cbuf))
    if (rank .ne. 0) call ppm_write(rank,trim(myname),trim(cbuf),info)

    if (sph_debug .gt. 1) then
      call ppm_write(rank,trim(myname),'done mapping particles on topo',info,unit_proclog)
      write(cbuf,fmt_char_int10) 'Local number of particles is now: ',Npart
      call ppm_write(rank,trim(myname),trim(cbuf),info,unit_proclog)
    endif

    !< Write successful end of mapping
    if (maptype .eq. ppm_param_map_global) then
      call logmessage(trim(myname),'global mapping done')
    elseif (maptype .eq. ppm_param_map_partial) then
      call logmessage(trim(myname),'partial mapping done')
    endif
  endif

  !-----------------------------
  !< Debug
  !-----------------------------
  if (sph_debug .gt. 0) then
    call substop(myname,ctime,info)
  endif


end subroutine mapping
