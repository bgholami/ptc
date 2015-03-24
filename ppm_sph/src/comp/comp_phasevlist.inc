!------------------------------------
!< Subroutine that creates interface lists
!------------------------------------
subroutine comp_phasevlist

  use mod_data_global
  use mod_data_sph
  use mod_data_prtl
  use mod_data_global
  use mod_data_ctrl
  use mod_data_physics

  use mod_io
  use mod_mapping
  use mod_util
	
  use ppm_module_neighlist
  use ppm_module_map_part_util

  implicit none

  !------------------------------
  !< Arguments
  !------------------------------

  !-------------------------------
  !< Local Variables
  !--------------------------------
  character(maxchar)	:: myname = 'comp_phasevlist'
  integer             :: info
  integer             :: j

  !---------------------------------
  !< Body of routine
  !---------------------------------
  if (sph_debug .gt. 0) call logmessage(myname,'entering...')

  !< Now get number of interface particles
  do j = 2,NVERLETLIST
    verletlist(j)%ldapidx = 0
  enddo

  !< Now allocate pidx
  do j = 2, NVERLETLIST
    deallocate(verletlist(j)%pidx)
    allocate(verletlist(j)%pidx(verletlist(j)%ldapidx))
  enddo

  !< null ldapidx and use as counter
  do j = 2, NVERLETLIST
    verletlist(j)%ldapidx = 0
  enddo

  !< Make lists for every phase + interface
  do j = 2, NVERLETLIST
    if (verletlist(j)%ldapidx .eq. Mpart) then
      call ppm_neighlist_vlist(topoid=topo_id,xp=xp,np=Mpart,cutoff=cutoff,&
      skin=2._mk*verlet_skin,lsymm=symmetry,vlist=verletlist(j)%vlist,&
      nvlist=verletlist(j)%nvlist,info=info)

    elseif (verletlist(j)%ldapidx .gt. 0) then
      call ppm_neighlist_vlist(topoid=topo_id,xp=xp,np=Mpart,cutoff=cutoff,&
      skin=2._mk*verlet_skin,lsymm=symmetry,vlist=verletlist(j)%vlist,&
      nvlist=verletlist(j)%nvlist,info=info,pidx=verletlist(j)%pidx,&
      npidx=verletlist(j)%ldapidx)
    endif

    if (info .ne. 0) call abortmessage(trim(myname),'multiple ppm_neighlist_vlist failed')
    if (sph_debug .gt. 0) call logmessage(trim(myname),'create multiple neighlist_vlist done')
  enddo

  return

end subroutine comp_phasevlist


