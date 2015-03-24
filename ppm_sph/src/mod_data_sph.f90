!---------------------------
!< Module that contains variables for sph stuff
!--------------------------------------------
module mod_data_sph

  use mod_data_global, only : mk
!  use ppm_module_data, only : ppm_t_clist

  implicit none

  !< cutoff
  REAL(MK) :: cutoff, cutoff_inv, cutoff2

  !< Ratio cutoff / dx (DEFAULT = 3)
  REAL(MK) :: ratio_cutoff_smoothing = 3.
  REAL(MK) :: ratio_smoothing_dx = 1.

  !----------------------------
  !< Smoothinglength
  !-----------------------------
  REAL(MK) :: smoothinglength = 0.
  REAL(MK) :: smoothinglength2
  REAL(MK) :: smoothinglength3

  !--------------------------
  !< Verlet list
  !--------------------------
  type list
    !< Verlet lists
    integer, dimension(:,:), pointer    :: vlist => null()
    integer, dimension(:  ), pointer    :: nvlist => null()
    !< Pidx for special neighlist with interface particles only
    integer, dimension(:), allocatable  :: pidx
    integer                             :: ldapidx = 0
  end type list

  !< List variables
  type(list), dimension(:), allocatable :: verletlist
  integer                               :: NVERLETLIST

  real(mk)                          :: verlet_skinfactor
  real(mk)                          :: verlet_skin !< enlarge pp in verlet
  real(mk)                          :: verlet_skin2

  !-----------------------------
  !< Number of particles and discretization
  !-------------------------------
  !< Particle distancing (initial)
  integer, dimension(__DIM)   :: Np
  real(mk)                    :: dpx
  real(mk)                    :: sph_bsize

  !< Number of particles of each phase
  integer, dimension(:), pointer  :: Ntotal => null()

  !< Number of particles + ghosts
  integer                         :: Npart, Mpart

  !< Total number of particles + ghosts
  integer                         :: sumNtotal, sumMtotal

  !-------------------------------------------
  !< Kernel variables
  !---------------------------------------------
  REAL(MK) :: norm
  REAL(MK) :: factorW, factorGradW, normedDist
  REAL(MK) :: ss3, ss2, ss1
  !< Zero contribution to kernel
  real(mk)  :: w_zero = 0._mk
  !< Kernel at dx -> artificial stress
  real(mk)  :: w_dx = 0._mk

end module mod_data_sph
