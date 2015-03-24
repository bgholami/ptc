!--------------------------------
!< Subroutine: define restart state of particles
!--------------------------------
subroutine init_state

  use mod_data_global
  use mod_data_ctrl
  use mod_data_physics
  use mod_data_sph
  use mod_data_prtl

  use mod_io
  use mod_comp
  use mod_mapping
  use mod_util

  use ppm_module_map_part
  use ppm_module_map_part_util
  use ppm_module_map_part_ghost
  use ppm_module_neighlist

  use ppm_module_substart
  use ppm_module_substop


  implicit none

  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(len=maxchar)    :: myname = 'init_state'
  integer                   :: i
  integer                   :: info

  !-----------------------------
  !< Set TIME
  !-----------------------------
  current_time = start_time
  current_step = 0
  restart_step_start = 0
  step_lastrebalance = 0

  !---------------------------------------
  !< Get ghost particles before building the verlet list
  !---------------------------------------
  call reset_ptr_to_list(realmaparray)
  call reset_ptr_to_list(intmaparray)
  call add_ptr_to_list(intmaparray,ap)
  call mapping(myname,param_maptype_getghost,realmaparray,intmaparray)
  call get_ptr_to_list(intmaparray,ap)
  call logmessage(trim(myname),'mapping getghost successfull')

  !----------------------------------------
  !< Store the mapping buffer
  !----------------------------------------
  call ppm_map_part_store(info)
  if (info .ne. 0) call abortmessage(myname,'ppm_map_part_store failed')

  !------------------------------------------
  !< Create list
  !------------------------------------------
  !< Before first time
  !call ppm_clist_destroy(clist,info)
  call ppm_neighlist_vlist(topo_id,xp,Mpart,cutoff,2._mk*verlet_skin,symmetry,&
  verletlist(1)%vlist,verletlist(1)%nvlist,info)!,clist=clist)
  if (info .ne. 0) call abortmessage(trim(myname),'ppm_neighlist_vlist failed')
  call logmessage(trim(myname),'create neighlist_vlist done')

  !------------------------------------------
  !< Store positions when vlist was created -> to detect when new vlist required
  !------------------------------------------
  allocate(dx_vlist(__DIM,Npart))
  dx_vlist = 0._mk

  !----------------------------------------------
  !< Show vlistresult after rebalancing if sph_debug > 0
  !----------------------------------------------
  if (sph_debug .gt. 0) call write_vlistresult

  !-------------------------------------------------------------
  !< Calc particle number density and set interface flags
  !< and calc/define particle volume
  !--------------------------------------------------------------
  allocate(drhodt(NDRHODT,Mpart))
  drhodt = 0._MK
  call comp_init_state

  !-----------------------------
  !< Set / calculate density, mass and h
  !-----------------------------
  do i = 1, Npart
    pdata(massidx,i) = pdata(volidx,i) * phase(ap(phaseidx,i))%density
    pdata(rhoidx,i) = phase(ap(phaseidx,i))%density
    pdata(hidx,i) = smoothinglength
  enddo

  !-----------------------------
  !< Set initial pressure
  !-----------------------------
  call press_EOS

  !--------------------------------
  !< Set initial velocity
  !--------------------------------
  allocate(dxpdt(NDXPDT,Mpart))
  dxpdt(1:NDXPDT,1:Mpart) = 0._mk

  call init_velocity
  call logmessage('init_state','Initial velocity defined')

  !---------------------------------------
  !< now get everything for the ghosts
  !---------------------------------------
  call reset_ptr_to_list(realmaparray)
  call reset_ptr_to_list(intmaparray)
  call add_ptr_to_list(realmaparray,dxpdt)
  call add_ptr_to_list(realmaparray,pdata)
  call mapping(myname,param_maptype_ghostupdatenoxp,realmaparray,intmaparray)
  call get_ptr_to_list(realmaparray,pdata)
  call get_ptr_to_list(realmaparray,dxpdt)
  call logmessage(trim(myname),'map_ghost_push successfull')

  !-----------------------------------------
  !< Apply velocity boundary condition
  !-----------------------------------------
  call comp_boundarycondition('dxpdt')

  !-------------------------------------------
  !< Initialize forces
  !-------------------------------------------
  allocate(pForce_total(ndim,Mpart))
  pForce_total = 0._mk
  allocate(d2xpdt2(ndim,Mpart))
  d2xpdt2 = 0._mk
  drhodt = 0._mk

  !-----------------------------------------
  !< Calculate forces and acceleration at t=0
  !-----------------------------------------
  call comp_forces
  call logmessage(myname,'Forces computed')

end subroutine init_state

