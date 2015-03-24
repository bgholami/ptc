!--------------------------------
!< Subroutine to do one step of scheme: kickdriftkick
!--------------------------------
subroutine run_kickdriftkick(dt_kickdriftkick)

  use mod_data_global
  use mod_data_ctrl
  use mod_data_prtl
  use mod_data_physics
  use mod_data_sph

  use mod_io
  use mod_util

  use mod_comp
  use mod_mapping

  use ppm_module_impose_part_bc
  use ppm_module_data
  use ppm_module_time
  use ppm_module_neighlist
  use ppm_module_topo
  use ppm_module_util_dbg
  use ppm_module_map_part_util

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  !-----------------------------
  !< Arguments
  !-----------------------------
  real(mk), intent(in)  :: dt_kickdriftkick

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'run_kickdriftkick'
  integer                     :: info
  logical                     :: topo_ok

  integer                     :: ierr
   
  !-----------------------------
  !< Velocity at half step u(n+1/2) - "KICK"
  !-----------------------------
  call timer_start(wct_timestep(velidx),info)
  call shift_velocity(0.5_mk*dt_kickdriftkick)
  call timer_stop(wct_timestep(velidx),info)

  !------------------------------------
  !< Position at next step x(n+1) - "DRIFT"
  !------------------------------------
  call timer_start(wct_timestep(posidx),info)
  call shift_position(dt_kickdriftkick)
  call timer_stop(wct_timestep(posidx),info)


  !------------------------------------
  !< Check if new vlist required -> check if particle moved out of skin
  !------------------------------------
  call timer_start(wct_timestep(checkidx),info)
  !< check if new verlet list required
  call check_newvlist(lnew=make_new_vlist)

  !< If new vlist: a) impose_part_bc (if not yet in rebalancing)
  !<               b) map partial (if not yet globally in rebalancing)
  if (make_new_vlist) then
    !< "Beam" particles if perodicity is used
    call ppm_impose_part_bc(topo_id,xp,Npart,info)
    call ppm_impose_part_bc(topo_id,xp,Npart,info) !< Omar: twice needed (bug)
    if (info .ne. 0) call abortmessage(trim(myname),'ppm_impose_part_bc failed.')
    !< Create list of arrays to map
    call makemaparray(realmaparray,intmaparray,info)
    !< Map everything
    call mapping(trim(myname),param_maptype_partial,realmaparray,intmaparray)
    !< Update local pointers
    call getmaparray(realmaparray,intmaparray,info)

    if (sph_debug .gt. 1) then
      !---------------------------------
      !< Check topology
      !---------------------------------
      call ppm_topo_check(topo_id,xp,Npart,topo_ok,info)
      if (info .ne. 0) call abortmessage(myname,'1 call ppm_topo_check failed')
      if (.not. topo_ok) call abortmessage(myname,'1 ppm_topo_check: topology check failed, something wrong...')
      call logmessage(myname,'topo_check OK')
    endif

  endif
  call timer_stop(wct_timestep(checkidx),info)


  !------------------------------------------------
  !< Get new particle data at ghosts
  !-------------------------------------------------
  call timer_start(wct_timestep(ghostidx),info)
  if (make_new_vlist) then
    !< Create list of arrays to map
    call makemaparray(realmaparray,intmaparray,info)
    !< New ghosts with properties
    call mapping(myname,param_maptype_getghost,realmaparray,intmaparray)
    !< Update local pointers from arraylist
    call getmaparray(realmaparray,intmaparray,info)
    call ppm_map_part_store(info)
    if (info .ne. 0) call abortmessage(myname,'ppm_map_part_store failed')

  else
    call ppm_map_part_load(info)
    if (info .ne. 0) call abortmessage(myname,'ppm_map_part_load failed')
    !< Create list of arrays to map
    call makemaparray(realmaparray,intmaparray,info)
    !< update ghosts
    call mapping(myname,param_maptype_ghostupdate,realmaparray,intmaparray)
    !< update local pointer
    call getmaparray(realmaparray,intmaparray,info)
  endif
  call timer_stop(wct_timestep(ghostidx),info)

  !-------------------------------------------------------------
  !< Apply boundary condition: velocity
  !-------------------------------------------------------------
  call timer_start(wct_timestep(bcidx),info)
  call comp_boundarycondition('dxpdt')
  call timer_stop(wct_timestep(bcidx),info)

  !-----------------------------------
  !< Build new vlist
  !-----------------------------------
  call timer_start(wct_timestep(vlistidx),info)
  if (make_new_vlist) then
    call ppm_neighlist_vlist(topo_id,xp,Mpart,cutoff,2._mk*verlet_skin,symmetry,&
    verletlist(1)%vlist,verletlist(1)%nvlist,info)!,clist=clist)
    if (info .ne. 0) call abortmessage(myname,'ppm_neighlist_vlist failed')
    !< counter
    count_newvlist = count_newvlist + 1
    ! next timestep do not make vlist
    make_new_vlist = .false.

    !----------------------------------------------
    !< Show vlistresult if sph_debug > 0
    !----------------------------------------------
    if (sph_debug .gt. 1) call write_vlistresult

    !------------------------------------------------
    !< Store positions when vlist is created
    !< to check movement < skin
    !----------------------------------------------
    if (size(dx_vlist,2) .lt. Npart) then
      deallocate(dx_vlist)
      allocate(dx_vlist(__DIM,2*Npart))
    endif
    dx_vlist = 0._mk

    !-------------------------
    !< Reshape pointers / arrays
    !-------------------------
    call reshapearrays

  endif
  call timer_stop(wct_timestep(vlistidx),info)

  !-------------------------------
  !< Null forces and drhodt
  !-------------------------------
  call timer_start(wct_timestep(resetidx),info)
  call resetarrays
  call timer_stop(wct_timestep(resetidx),info)

  !----------------------------------------------
  !< Calculate PP-INTERACTIONS "1"
  !< Density summation
  !-----------------------------------------------
  !< Timing inside...
  call comp_particleloop_1

  !---------------------------------------------
  !< Shift/Update density
  !---------------------------------------------
  call timer_start(wct_timestep(densidx),info)  !call ppm_time(wct_shiftdensity_s,info)
  call shift_density(dt_kickdriftkick)
  call timer_stop(wct_timestep(densidx),info)  !call ppm_time(wct_shiftdensity_s,info)

  !-----------------------
  !< Properties
  !----------------------
  call timer_start(wct_timestep(propidx),info)
  call press_EOS
  call timer_stop(wct_timestep(propidx),info)

  !------------------------------------
  !< Push density, volume and pressure now
  !------------------------------------
  call timer_start(wct_timestep(pushidx),info)
  call reset_ptr_to_list(realmaparray)
  call reset_ptr_to_list(intmaparray)
  call add_ptr_to_list(realmaparray,pdata)
  call mapping(myname,param_maptype_ghostupdatenoxp,realmaparray,intmaparray)
  call get_ptr_to_list(realmaparray,pdata)
  call timer_stop(wct_timestep(pushidx),info)

  !-------------------------------------------------------------
  !< PARTICLE INTERACTIONS
  !< Force calculation
  !-------------------------------------------------------------
  !< Timing inside
  call comp_forces

  !-------------------------------------------------------------
  !< PUSH
  !-------------------------------------------------------------
  call timer_start(wct_timestep(pushidx),info)
  call reset_ptr_to_list(realmaparray)
  call reset_ptr_to_list(intmaparray)
  call add_ptr_to_list(realmaparray,d2xpdt2)
  call mapping(myname,param_maptype_ghostupdatenoxp,realmaparray,intmaparray)
  call get_ptr_to_list(realmaparray,d2xpdt2)
  call timer_stop(wct_timestep(pushidx),info)

  !-----------------------------------
  !< Velocity at full step u(n+1) - "KICK"
  !----------------------------------
  call timer_start(wct_timestep(velidx),info)
  call shift_velocity(0.5_mk*dt_kickdriftkick)
  call timer_stop(wct_timestep(velidx),info)

  return

end subroutine run_kickdriftkick
