!--------------------------------
!< Subroutine: here the actual work is done !!
!--------------------------------
subroutine run_steps
  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_prtl

  use mod_io
  use mod_util
  use mod_mapping
  use mod_util

  use ppm_module_time
  use ppm_module_map_part_util
  use ppm_module_util_dbg

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
  character(maxchar)          :: myname = 'run_steps'
  integer                     :: info
  character(maxchar)          :: cbuf
  real(mk)              :: t_next_output
  logical               :: abortexist = .false.

  real(mk)              :: tsum

  integer               :: sum_Ntotal, sum_Mtotal

#ifdef __MPI
  integer               :: ierr
#endif

  !-----------------------------------------------------
  !< Total number of iterations to be done
  !-----------------------------------------------------
  nsteps = INT((end_time - start_time) / dt_global + 0.5_MK)

  WRITE(cbuf,*) '----------------------------------'
  CALL logmessage('run_steps',cbuf)
  WRITE(cbuf,*) '|  Starting main time loop...'
  CALL logmessage('run_steps',cbuf)
  WRITE(cbuf,*) '----------------------------------'
  CALL logmessage('run_steps',cbuf)

  WRITE(cbuf,fmt_char_int10) 'Number of timesteps: ',nsteps
  CALL logmessage('run_steps',TRIM(cbuf))
  WRITE(cbuf,fmt_char_real) 'Timestep: ',dt_global
  CALL logmessage('run_steps',TRIM(cbuf))
  WRITE(cbuf,fmt_char_real) 'Start time of simulation: ',current_time
  CALL logmessage('run_steps',TRIM(cbuf))
  WRITE(cbuf,fmt_char_real) 'End time: ',end_time
  CALL logmessage('run_steps',TRIM(cbuf))

  !-----------------------------
  !< Null timer
  !-----------------------------
  call timer_reset(wct_step,info)

  !-----------------------------
  !< Gather timing information
  !-----------------------------
  !< Measure total time from start to end of run_steps
  call MPI_Barrier(comm,ierr)
  call timer_start(wct_step,info)
  step_wctstep_s = current_step

  !-----------------------------
  !< MAIN LOOP
  !-----------------------------
  do while (current_step .lt. nsteps .and. &
  (current_step .lt. maxitsteps))

    !-------------------------------------
    !< Update new stepnumber and actual time
    !-------------------------------------
    current_step = current_step + 1
    current_time = current_time + dt_global

    !----------------------------------
    !< Performance timer
    !----------------------------------
    call timer_reset(wct_timestep,info)
    call timer_start(wct_timestep(stepidx),info)

    !----------------------------------
    !< DO ONE TIMESTEP
    !----------------------------------
    select case (timesteppingscheme)
      case (param_timestepping_kickdriftkick)
        call run_kickdriftkick(dt_global)

      case default
        write(cbuf,fmt_char_int2) 'Non-existing timesteppingscheme: ',&
        timesteppingscheme
        call abortmessage('run_steps',trim(cbuf))

    end select

    !< Accumulate runtime of timestep
    call timer_stop(wct_timestep(stepidx),info)

    !---------------------------------
    !< Write output
    !---------------------------------
    t_next_output = abs(NINT(current_time/(timeoutput))*timeoutput-current_time)
    if (t_next_output .lt. 0.5_mk*dt_global) then
      call timer_start(wct_output,info)
      call write_output
      if (ppm_debug .eq. 1) call ppm_dbg_print(topo_id,sph_bsize,current_step,1,info)
      if (ppm_debug .eq. 2) call ppm_dbg_print(topo_id,sph_bsize,current_step,1,info,xp,Npart,Mpart)
      call timer_start(wct_output,info)
    endif

    !----------------------------
    !< Check abortfile
    !----------------------------
    call check_abortfile(abortexist)
    if (abortexist) exit

    !-----------------------
    !< Write restart
    !-----------------------
    !if (choice_restart .ne. 0 .and. mod(current_step,freqrestart) .eq. 0) then
    !  call write_restart
    !endif

    !---------------------------------
    !< Write conservation
    !---------------------------------
    !    if (mod(current_step,freqlogoutput) .eq. 0) then
    !      call write_conservation
    !    endif

    !---------------------------------
    !< Write log from time to time
    !---------------------------------
    if (current_step .eq. 1 .or. mod(current_step,freqlogoutput) .eq. 0) then
      call timer_stop(wct_step,info)
      call write_logoutput
      !------------------------------------
      !< Show detailed performance analyse
      !< of every code section
      !-------------------------------------
      if (switch_showtimings) then
        call show_timings
        !call write_timings
      endif
    endif

  !< End of main step loop
  enddo

  !<---------------------------------------------
  !< Stop timing
  !<---------------------------------------------
  call MPI_Barrier(comm,ierr)
  call timer_stop(wct_step,info)

  !-----------------------------------------------
  !< Write output if ABORTED or MAXITSTEPS met
  !< and not an outputtime anyway
  !< Write restart file at the end if not already written this step
  !-----------------------------------------------
  t_next_output = abs(NINT(current_time/timeoutput)*timeoutput-current_time)
  if (t_next_output .ge. 0.5_mk*dt_global) then
    call write_output
    if (ppm_debug .eq. 1) call ppm_dbg_print(topo_id,0._mk,current_step,1,info)
    if (ppm_debug .eq. 2) call ppm_dbg_print(topo_id,sph_bsize,current_step,1,info,xp,Npart,Mpart)
    !call write_restart
  elseif (mod(current_step,freqrestart) .ne. 0) then
    !call write_restart
  endif

  !----------------------------------
  !< Timing results and calc details
  !----------------------------------
  sumNtotal = Npart
#ifdef  __MPI
  call MPI_ALLREDUCE(sumNtotal,sum_Ntotal,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  sumNtotal = sum_Ntotal
#endif

  WRITE(cbuf,'(A)') '_________________________________________________________'
  CALL logmessage(trim(myname),TRIM(cbuf))
  WRITE(cbuf,fmt_char_int10) 'Global number of particles: ', sumNtotal
  CALL logmessage(trim(myname),TRIM(cbuf))

  !< Ghosts
  sumMtotal = Mpart-Npart
#ifdef  __MPI
  call MPI_ALLREDUCE(sumMtotal,sum_Mtotal,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  sumMtotal = sum_Mtotal
#endif
  WRITE(cbuf,fmt_char_int10) 'Global number of ghost particles: ', sumMtotal
  CALL logmessage(trim(myname),TRIM(cbuf))

  WRITE(cbuf,fmt_char_int10) 'Number of completed time steps: ',&
  current_step
  CALL logmessage(trim(myname),TRIM(cbuf))

  write(cbuf,fmt_char_real) 'Current time: ',current_time
  call logmessage(myname,cbuf)

  IF (current_step .GE. maxitsteps) THEN
    WRITE(cbuf,fmt_char) 'STOPPED BECAUSE MAXITSTEPS MET !'
    CALL logmessage(trim(myname),TRIM(cbuf))
  ENDIF
  IF (abortexist) THEN
    WRITE(cbuf,'(A)') 'STOPPED BECAUSE ABORTFILE FOUND !'
    CALL logmessage(trim(myname),TRIM(cbuf))
  ENDIF

  !<--------------------------------------------------------
  !< Total time since start of simulation
  !<--------------------------------------------------------
  tsum = 0._mk
#ifdef __MPI
  call MPI_ALLREDUCE(wct_step%delta,tsum,1,MPPREC,MPI_SUM,comm,ierr)
  tsum = tsum / Nproc
#else
  tsum = wct_step%delta
#endif

  if (current_step - restart_step_start .gt. 0) then
    WRITE(cbuf,fmt_char_real) 'Avg. seconds per time step: ',&
    tsum/(REAL(current_step - restart_step_start,mk)+myeps)
    CALL logmessage(trim(myname),cbuf)

    WRITE(cbuf,fmt_char_exp) 'Performance (Time/Particle/Iteration): ',&
    tsum/REAL(current_step - restart_step_start,mk)/&
    real(sumNtotal,KIND=mk) * real(Nproc,KIND=mk)
    CALL logmessage(trim(myname),TRIM(cbuf))
  endif

  write(cbuf,fmt_char_int10) 'Number of rebalancings: ',count_rebalancing
  call logmessage(myname,cbuf)

  if (count_newvlist .lt. 1) then
    write(cbuf,fmt_char) 'Verlet list only once created!'
    call logmessage(myname,cbuf)
  else
    write(cbuf,fmt_char_int10) 'Avg. steps between make_new_vlist: ',&
    int((current_step - restart_step_start) / count_newvlist)
    call logmessage(myname,cbuf)
  endif

  !< Wct-time for total initialization
  if (wct_prestep%delta .lt. 60._mk) then
    WRITE(cbuf,fmt_char_real_char) 'Time to setup timesteps: ',&
    wct_prestep%delta,' sec. '
  else
    WRITE(cbuf,fmt_char_real_char) 'Time to setup timesteps: ',&
    wct_prestep%delta,' sec. '
  endif
  CALL logmessage(trim(myname),cbuf)

  !< Calculation time for steps
  if (tsum .lt. 60._mk) then
    WRITE(cbuf,fmt_char_real_char) 'All timesteps: ',&
    tsum,' sec. '
  else
    WRITE(cbuf,fmt_char_real_char) 'All timesteps: ',&
    tsum/60._mk,' min. '
  endif
  CALL logmessage(trim(myname),cbuf)

end subroutine run_steps
