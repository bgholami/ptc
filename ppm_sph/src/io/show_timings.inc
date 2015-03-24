!--------------------------------
!< Subroutine that does all the logoutput
!--------------------------------
subroutine show_timings

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_util

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
  character(maxchar)          :: myname = 'show_timings'
  character(maxchar)          :: cbuf
  character(maxchar)          :: fmt_local
  integer                     :: info
  real(mk)                    :: rsteps

  real(mk), dimension(n_procbal)      :: procperf !minmaxmeanvar
  real(mk)                            :: timestepmean
  real(mk)                            :: ttotal
  integer                             :: j

  integer                             :: sumpp

  !-----------------------------
  !< Check total time of one timestep
  !-----------------------------
  write(cbuf,'(A)') '****************************************************************************'
  call logmessage(myname,cbuf)
  write(cbuf,'(A)') '************************* SHOW DETAILED TIMINGS ****************************'
  call logmessage(myname,cbuf)
  write(cbuf,'(A)') 'Timings are:       min      |     max      |     mean     |     var      |   %    |  calls'
  call logmessage(myname,cbuf)

  call minmaxmeanvar(wct_timestep(stepidx)%t_sum,procperf,info)
  write(fmt_local,'(A)') '(A15,E12.5," | ",E12.5," | ",E12.5," | ",E12.5)'
  write(cbuf,trim(fmt_local)) 'timestep (s):', procperf(minidx),procperf(maxidx),&
  procperf(meanidx),procperf(varidx)
  call logmessage(myname,cbuf)
  write(cbuf,'(A)') '_____________________________________________________________________________'
  call logmessage(myname,cbuf)
  timestepmean = procperf(meanidx)

  !< Show all subtimings:
  ttotal = 0._mk
  do j = 2, ntimings
    call minmaxmeanvar(wct_timestep(j)%t_sum,procperf,info)
    ttotal = ttotal + procperf(meanidx)
    write(fmt_local,'(A)') '(A10," (s):",E12.5," | ",E12.5," | ",E12.5," | ",E12.5," | ",F6.2," | ",I2)'
    write(cbuf,trim(fmt_local)) trim(wct_name(j)),procperf(minidx),procperf(maxidx),procperf(meanidx),&
    procperf(varidx),procperf(meanidx)/timestepmean*100._mk,wct_timestep(j)%calls
    call logmessage(myname,cbuf)
  enddo
  ttotal = timestepmean - ttotal
 
  !< Sum interactions
  sumpp = 0
  do j = 1, Npart
    sumpp = sumpp + verletlist(1)%nvlist(j)
  enddo
  call minmaxmeanvar(real(sumpp,mk),procperf,info)
  write(fmt_local,'(A)') '(A10,"    :",E12.5," | ",E12.5," | ",E12.5," | ",E12.5)'
  write(cbuf,trim(fmt_local)) 'proc-pp',procperf(minidx),procperf(maxidx),procperf(meanidx),&
  procperf(varidx)
  call logmessage(myname,cbuf)

  write(cbuf,'(A)') '_____________________________________________________________________________'
  call logmessage(myname,cbuf)
  !< Write remaining time to see if something is missing...
  write(fmt_local,'(A)') '(A10," (s):",E12.5," | ",E12.5," | ",E12.5," | ",E12.5," | ",F6.2," | ",I2,")")'
  write(cbuf,trim(fmt_local)) "(Rest: ",ttotal,ttotal,ttotal,0._mk,ttotal/timestepmean*100._mk,0
  call logmessage(myname,cbuf)

  write(cbuf,'(A)') '****************************************************************************'
  call logmessage(myname,cbuf)

end subroutine show_timings
