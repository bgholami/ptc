!---------------------------------
!< Module containing all prtl parameter
!---------------------------------
module mod_data_physics
  use mod_data_global
  implicit none

  !< compbox size
  real(mk), dimension(__DIM)    :: min_compbox, max_compbox
  real(mk), dimension(__DIM)    :: min_domain, max_domain


  !< boundaries
  integer, dimension(2*__DIM)         :: bcdef, bcdef_ppm
  real(mk), dimension(2*__DIM,__DIM)  :: bcvel

  !< Body force
  real(mk), dimension(__DIM)    :: bodyforce

  !< Ref values
  real(mk)                      :: length_ref, vel_ref, density_ref

  !< Phases
  integer, parameter            :: param_max_n_phases = 1
  integer                       :: N_PHASES

  !---------------------------------
  !< Type definition of phase
  !---------------------------------
  TYPE material
    CHARACTER(LEN=MAXCHAR) :: name = 'NONAME'
    REAL(MK) :: density = 1._MK
    INTEGER  :: type = 2
    REAL(MK) :: visc_dyn = 1._MK
    REAL(MK) :: visc_kin = 0._MK
    REAL(MK) :: gamma = 7._MK
    REAL(MK) :: delta_dens = 0.01_mk
    REAL(MK), DIMENSION(NDTCRITERIA) :: timestep = 1.e99_MK
    REAL(MK), DIMENSION(NCSOUNDCRITERIA) :: soundspeed = 0._MK
    REAL(MK) :: cmax = 0._MK
    REAL(MK) :: press_ref = 0._MK
    REAL(MK) :: press_b = 1._MK
    INTEGER :: output = 0
  END TYPE material

  !------------------------------------------------------
  !< Phases
  !------------------------------------------------------
  type(material), dimension(param_max_n_phases) :: phase

end module mod_data_physics
