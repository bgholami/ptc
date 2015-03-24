!--------------------------------
!< Subroutine that defines boundary condition on particles Npart+1 - Mpart
!--------------------------------
subroutine comp_boundarycondition(char_in)

  use mod_data_global
  use mod_data_physics
  use mod_data_prtl
  use mod_data_sph
  use mod_data_ctrl


  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------
  character, intent(in), optional  :: char_in*(*)

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(maxchar)          :: myname = 'comp_boundarycondition'
  !integer                     :: info
  integer                     :: j
  integer                     :: ip

  !-----------------------------
  !< Do work
  !-----------------------------
  !< Velocity boundary condigion

  select case (trim(char_in))

    !-----------------------------
    !< VELOCITY BOUNDARY CONDITION
    !-----------------------------
    case ('dxpdt','DXPDT')
      do ip = Npart+1, Mpart
        do j = 1, ndim
          !<-------------------------------
          !< Upper boundary
          if (xp(j,ip) .ge. max_compbox(j)) then
            !< Symmetry
            if (bcdef(2*j) .eq. param_bcdef_symmetry) then
              dxpdt(j,ip) = -dxpdt(j,ip)
            !< Mirror wall
            elseif (bcdef(2*j) .eq. param_bcdef_mirrorwall) then
              dxpdt(1:ndim,ip) = 2._mk*bcvel(2*j,1:ndim) - dxpdt(1:ndim,ip)
            !< Antisymmetry
            elseif (bcdef(2*j) .eq. param_bcdef_antisymmetry) then
              !< do nothing
            endif

          !<-------------------------------
          !< Lower boundary
          elseif (xp(j,ip) .le. min_compbox(j)) then
            !< Symmetry
            if (bcdef(2*j-1) .eq. param_bcdef_symmetry) then
              dxpdt(j,ip) = -dxpdt(j,ip)
            !< Mirror wall
            elseif (bcdef(2*j-1) .eq. param_bcdef_mirrorwall) then
              dxpdt(1:ndim,ip) = 2._mk*bcvel(2*j-1,1:ndim) - dxpdt(1:ndim,ip)
            !< Anti-symmetry
            elseif (bcdef(2*j-1) .eq. param_bcdef_antisymmetry) then
              !< do nothing
            endif
          endif
        enddo
      enddo

  end select

end subroutine comp_boundarycondition
