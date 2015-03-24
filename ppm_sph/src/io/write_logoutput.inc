!--------------------------------
!< Subroutine that does all the logoutput
!--------------------------------
subroutine write_logoutput

  use mod_data_global
  use mod_data_sph, only: sumNtotal
  use mod_data_ctrl
  use mod_util, only : allreduce_mean


  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(maxchar)          :: myname = 'write_logoutput'
  character(maxchar)          :: cbuf

  real(8)                     :: avg_wctstep, avgwctstep
  real(8)                     :: wct_togo_hours, wct_togo_min


  !-----------------------------
  !< 
  !-----------------------------
  !< Get avgwct time now
  avg_wctstep = wct_step%delta / (current_step - step_wctstep_s)

  !< Get mean performance
  avgwctstep = 0.
  call allreduce_mean(avg_wctstep,avgwctstep)

  write(cbuf,'(A)') '_________________________________________________________________'
  call logmessage(myname,cbuf)

  WRITE(cbuf,'(A,I12,A,I12,A,E17.7)') 'Completed steps: ',current_step,&
  ' (of ',nsteps,') with dt: ',dt_global
  CALL logmessage(trim(myname),TRIM(cbuf))

  WRITE(cbuf,fmt_char_real) 'Avg. time for a time step (s): ',avgwctstep
  CALL logmessage(trim(myname),TRIM(cbuf))

  WRITE(cbuf,fmt_char_exp) 'Performance (Time/Particle/Iteration) : ',&
  avgwctstep/real(sumNtotal,KIND=mk)*real(Nproc,KIND=mk)
  CALL logmessage(trim(myname),TRIM(cbuf))

  wct_togo_hours = INT(REAL((nsteps-current_step),MK)*&
  avgwctstep/60._mk/60._mk)
  wct_togo_min = REAL((nsteps-current_step),MK)*&
  avgwctstep/60._MK - wct_togo_hours*60._MK
  WRITE(cbuf,'(3A,F10.2,A,F10.2,A)') 'WCT to finish (approx.) ',&
  TRIM(ctrlfile),' : ',wct_togo_hours,' h and ',wct_togo_min,' min.'
  CALL logmessage(trim(myname),TRIM(cbuf))

  WRITE(cbuf,'(A,E17.7,A,E17.7,A)') 'Calculated time: ',current_time,' of ',&
  end_time,' sec.'
  CALL logmessage(trim(myname),TRIM(cbuf))

end subroutine write_logoutput
