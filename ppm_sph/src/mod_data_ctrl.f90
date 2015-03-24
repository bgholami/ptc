!---------------------------------
!< Module containing all ctrl file parameter
!---------------------------------
module mod_data_ctrl

  use mod_data_global, only : mk, maxchar

  implicit none

  !< Number of dimensions
  integer         :: ndim

  !< simcase
  integer                       :: simcase

  !< Flags (implicit - explicit)
  integer                       :: flag_density_equation = 0

  !< Restart
  logical                       :: restart = .false.
  logical                       :: restart_reset_steps = .false.
  logical                       :: restart_reset_velocity = .false.
  integer                       :: choice_restart

  !< Timestepping and output
  integer                       :: timesteppingscheme
  real(mk)                      :: start_time, end_time
  real(mk)                      :: timeoutput
  integer                       :: choice_output
  logical                       :: choice_output_nondim

  !< Get which output
  logical                       :: output_velocity
  logical                       :: output_density
  logical                       :: output_pressure
  logical                       :: output_volume
  logical                       :: output_mass
  logical                       :: output_h
  logical                       :: output_phaseidx
  logical                       :: output_id
  logical                       :: output_rank

  !< restart
  integer                       :: freqrestart, freqlogoutput

  !< max iterations
  integer                       :: maxitsteps

  !< model switches
  logical                       :: switch_showtimings

  !< PPM - Decomposition and assignment to procs
  integer                       :: choice_assign
  integer                       :: choice_decomp

end module mod_data_ctrl
