!---------------------------------
!< Module containing all physics variables
!---------------------------------
module mod_data_prtl
  use mod_data_global, only : mk, t_realarray, t_intarray

  implicit none

  !----------------------
  !< Particle data
  !----------------------
  !< Array of pointers -> mapping...
  type(t_realarray), dimension(:), pointer :: realmaparray => null()
  type(t_intarray), dimension(:), pointer :: intmaparray => null()


  !< Position
  real(mk), dimension(:,:), pointer :: xp => null()

  !< Velocity
  real(mk), dimension(:,:), pointer :: dxpdt => null()

  !< Pdata
  real(mk), dimension(:,:), pointer :: pdata => null()

  !< ID's "attributes of particles"
  integer, dimension(:,:), pointer  :: ap => null()

  !< Acceleration
  real(mk), dimension(:,:), pointer :: d2xpdt2 => null()

  !< Density change rate
  real(mk), dimension(:,:),   pointer :: drhodt => null()

  !< Position when vlist is created
  real(mk), dimension(:,:), pointer :: dx_vlist => null()

  !--------------------------
  !< Forces + change rates
  !-----------------------
  real(mk), dimension(:,:), allocatable :: pForce_total

  !------------------------
  !< Size of vectors
  !------------------------
  integer, parameter  :: NXP = __DIM
  integer, parameter  :: NDXPDT = __DIM
  integer, parameter  :: ND2XPDT2 = __DIM

  integer, parameter  :: NPDATA = 5
  integer, parameter  :: rhoidx = 1
  integer             :: pressidx = 2
  integer             :: volidx = 3
  integer             :: massidx = 4
  integer             :: hidx = 5

  integer, parameter  :: NAP = 2
  integer             :: phaseidx = 1
  integer             :: partIDidx = 2

  integer, parameter  :: NDRHODT = 1

end module mod_data_prtl
