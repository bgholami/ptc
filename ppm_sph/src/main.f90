! ============================================================================
!< Name        : Main program
!< Author      : Stefan Adami and Omar Awile
!< Version     :
!< Copyright   :
!< Description :
! ============================================================================
program main

  use mod_data_global
  use mod_data_sph
  use mod_data_physics
  use mod_data_prtl
  use mod_data_ctrl

  use mod_util
  use mod_init
  use mod_mapping
  use mod_io
  use mod_run

  use ppm_module_data
  use ppm_module_topo
  use ppm_module_write
  use ppm_module_find_duplicates
  use ppm_module_time
  use ppm_module_loadbal
  use ppm_module_util_dbg
  use ppm_module_map

  use ppm_module_substart
  use ppm_module_substop

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  !------------------------------------------
  !< Local variables
  !------------------------------------------
  character(maxchar)      :: myname = 'main'
  integer                 :: info                ! < info variable
  logical                 :: topo_ok

  !-------------------------------------------------------------------------
  !< Measure time until start of steps
  !-------------------------------------------------------------------------
  call timer_reset(wct_prestep,info)
  call timer_start(wct_prestep,info)

  !-------------------------------------------------------------------------
  !<  Initialise
  !--------------------------------------------------------------------------
  info = 0
  call initialize(info)


  !-------------------------------------------------------------------------
  !<  Create particles
  !-------------------------------------------------------------------------
  if (.not. restart) then
    call create_particles
  else
    !call restart_particles
  endif

  !-------------------------------------------------------------------------
  !<  MAKE TOPOLOGY
  !-------------------------------------------------------------------------
  topo_id = 0
  if (restart) then
    !< First create purely geometry based topo
    call ppm_mktopo(topo_id,choice_decomp,choice_assign,min_compbox,&
    max_compbox,bcdef_ppm,sph_bsize,cost,info)
    !----------------------------------------------------------------------
    !< Particle positions in restart file might
    !< exceed max_compbox/min_compbox due to verlet lists...
    !< Impose part_bc to ensure that they are inside domain
    !----------------------------------------------------------------------
    call ppm_impose_part_bc(topo_id,xp,Npart,info)
    call ppm_impose_part_bc(topo_id,xp,Npart,info) !< Omar: twice needed (bug)
  endif

  call ppm_mktopo(topo_id,xp,Npart,choice_decomp,choice_assign,min_compbox,&
  max_compbox,bcdef_ppm,sph_bsize,cost,info)
  if (info .ne. 0) call abortmessage('main','ppm_mktopo failed !!!')
  call logmessage('main','PPM MAKE TOPOLOGY SUCCESSFUL!')


  !-------------------------------------------------------------------------
  !< Map everything on the topology globally
  !-------------------------------------------------------------------------
  if (restart) then
  !    call map_part_everything(myname,maptype=ppm_param_map_global,&
  !    intvec1=ap,realvec1=pdata,realvec2=dxpdt,realvec3=d2xpdt2,&
  !    realvec4=drhodt,realvec5=surfspec,realvec6=dmGdtdiff,&
  !    realvec7=tensor,realvec8=dtensordt,realvec9=dmCdtdiff,&
  !    realvec10=bulkspec)
  else
    call reset_ptr_to_list(realmaparray)
    call reset_ptr_to_list(intmaparray)
    call add_ptr_to_list(realmaparray,pdata)
    call add_ptr_to_list(intmaparray,ap)
    call mapping(myname,param_maptype_global,realmaparray,intmaparray)
    call get_ptr_to_list(intmaparray,ap)
    call get_ptr_to_list(realmaparray,pdata)
  endif
  call logmessage('main','Mapped on topology!')


  !---------------------------------
  !< Check topology
  !---------------------------------
  call ppm_topo_check(topo_id,xp,Npart,topo_ok,info)
  if (info .ne. 0) call abortmessage('main','call ppm_topo_check failed')
  if (.not. topo_ok) call abortmessage('main','ppm_topo_check: topology check failed, something wrong...')
  call logmessage('main','Topology check ok!')

  !-----------------------------------
  !< Check uniqueness of particles
  !-----------------------------------
  call ppm_find_duplicates(xp,NXP,Npart,nident,ident,info)
  if (info .ne. 0) call abortmessage('main','call ppm_find_duplicates failed')
  if (associated(ident)) deallocate(ident)
  if (nident .ne. 0) call abortmessage('main','ppm_find_duplicates: collocating particles found...')
  call logmessage('main','No duplicates found')

  !-----------------------------------
  !< Define initial state
  !-----------------------------------
  if (.not. restart) then
    call init_state
    call logmessage(myname,'Initial state defined')
  else
    !call restart_state
    !call logmessage(myname,'Restart state defined')
    call abortmessage(myname,'Restart not yet available...')
  endif

  !--------------------------
  !< Write first output
  !--------------------------
  if (.not. restart) then
    call write_output
  else
    !< If restart only write output noe if reset steps used
    if (restart_reset_steps) call write_output
  endif

  !-------------------------------
  !< Measure time until starting steps...
  !-------------------------------
  call timer_stop(wct_prestep,info)

  !-------------------------------
  !< Do the work
  !-------------------------------
  call run_steps

  !-------------------------------------------------------------------------
  !< Finalize
  !-------------------------------------------------------------------------
  call abortmessage('main','PROGRAM FINISHED SUCCESSFULLY')
  stop

end program main
