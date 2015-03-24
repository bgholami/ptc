!------------------------------------
!< Subroutine check ctrlfile
!------------------------------------
subroutine check_ctrlfile

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_physics
  use mod_io

  use ppm_module_data

  implicit none

  !------------------------------
  !< Arguments
  !------------------------------

  !-------------------------------
  !< Local Variables
  !--------------------------------
  character(maxchar)        :: myname = 'check_ctrlfile'
  real(mk)                  :: boxvol, minvol
  integer                   :: i

  !---------------------------------
  !< Body of routine
  !---------------------------------
  
  !< Check ndim
  if (ndim .ne. __DIM) call abortmessage(trim(myname),'Check input: ndim .ne. __DIM!!!')
  
  !< Check minimum size for one proc
  minvol = 1._mk
  boxvol = 1._mk
  do i = 1, ndim
    boxvol = boxvol*(max_domain(i)-min_domain(i))
    minvol = minvol*cutoff*2
  enddo
  if (int(boxvol/minvol) .lt. Nproc) then
    call abortmessage(trim(myname),'Cutoff too large for that number of procs')
  endif

  !---------------------------------------------------------------
  !< Change bcdef according to ppm_param.h
  !---------------------------------------------------------------
  do i = 1, 2*ndim
    !<1 = periodic
    if (bcdef(i) .eq. param_bcdef_periodic) then
      bcdef_ppm(i) = ppm_param_bcdef_periodic
      cycle
    endif
    !<2 = open boundary / free space
    if (bcdef(i) .eq. param_bcdef_openboundary) then
     bcdef_ppm(i) = ppm_param_bcdef_freespace
      cycle
    endif
    !<3 = symmetry
    if (bcdef(i) .eq. param_bcdef_symmetry) then
     bcdef_ppm(i) = ppm_param_bcdef_symmetry
      cycle
    endif
    !<4 = neumann
    if (bcdef(i) .eq. param_bcdef_neumann) then
     bcdef_ppm(i) = ppm_param_bcdef_neumann
      cycle
    endif
    !<5 = dirichlet
    if (bcdef(i) .eq. param_bcdef_dirichlet) then
     bcdef_ppm(i) = ppm_param_bcdef_dirichlet
      cycle
    endif
    !<6 = robin
    if (bcdef(i) .eq. param_bcdef_robin) then
     bcdef_ppm(i) = ppm_param_bcdef_robin
      cycle
    endif
    !<7 = mirror wall
    if (bcdef(i) .eq. param_bcdef_mirrorwall) then
     bcdef_ppm(i) = ppm_param_bcdef_symmetry
      cycle
    endif
    if (bcdef(i) .eq. param_bcdef_injection) then
     call abortmessage(myname,'bcdef 8 not defined')
      cycle
    endif

    if (bcdef(i) .eq. param_bcdef_solidwall) then
     call abortmessage(myname,'bcdef 9 not defined')
      cycle
    endif
    !<10 = antisymmetry
    if (bcdef(i) .eq. param_bcdef_antisymmetry) then
     bcdef_ppm(i) = ppm_param_bcdef_antisymmetry
      cycle
    endif
  enddo

  !---------------------------------------------------------------
  !< Change decomp according to ppm_param.h
  !---------------------------------------------------------------
  if (choice_decomp .eq. param_decomp_tree) choice_decomp = ppm_param_decomp_tree
  if (choice_decomp .eq. param_decomp_pruned_cell) choice_decomp = ppm_param_decomp_pruned_cell
  if (choice_decomp .eq. param_decomp_bisection) choice_decomp = ppm_param_decomp_bisection
  if (choice_decomp .eq. param_decomp_xpencil) choice_decomp = ppm_param_decomp_xpencil
  if (choice_decomp .eq. param_decomp_ypencil) choice_decomp = ppm_param_decomp_ypencil
  if (choice_decomp .eq. param_decomp_zpencil) choice_decomp = ppm_param_decomp_zpencil
  if (choice_decomp .eq. param_decomp_cuboid) choice_decomp = ppm_param_decomp_cuboid
  if (choice_decomp .eq. param_decomp_user_defined) choice_decomp = ppm_param_decomp_user_defined
  if (choice_decomp .eq. param_decomp_xy_slab) choice_decomp = ppm_param_decomp_xy_slab
  if (choice_decomp .eq. param_decomp_xz_slab) choice_decomp = ppm_param_decomp_xz_slab
  if (choice_decomp .eq. param_decomp_yz_slab) choice_decomp = ppm_param_decomp_yz_slab
  if (choice_decomp .eq. param_decomp_cartesian) choice_decomp = ppm_param_decomp_cartesian

  !---------------------------------------------------------------
  !< Change assign according ppm_param.h
  !---------------------------------------------------------------
  if (choice_assign .eq. param_assign_internal) choice_assign = ppm_param_assign_internal
  if (choice_assign .eq. param_assign_nodal_cut) choice_assign = ppm_param_assign_nodal_cut
  if (choice_assign .eq. param_assign_nodal_comm) choice_assign = ppm_param_assign_nodal_comm
  if (choice_assign .eq. param_assign_dual_cut) choice_assign = ppm_param_assign_dual_cut
  if (choice_assign .eq. param_assign_dual_comm) choice_assign = ppm_param_assign_dual_comm
  if (choice_assign .eq. param_assign_user_defined) choice_assign = ppm_param_assign_user_defined

end subroutine check_ctrlfile


