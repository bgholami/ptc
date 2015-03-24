!--------------------------------
!< Give an overview about all parameter
!--------------------------------
subroutine dump_setup

  use mod_data_global
  use mod_data_ctrl
  use mod_data_physics
  use mod_data_sph

  use mod_io

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  character(len=maxchar)      :: cbuf
  character(len=maxchar)      :: fmt_local
  integer                     :: i

  !-----------------------------
  !< Body of routine
  !-----------------------------

  !---------------------------
  !< Symmetric interactions?
  !---------------------------
  write(cbuf,trim(fmt_char_boolean)) 'Interactions symmetric: ',symmetry
  CALL logmessage('dump_setup',TRIM(cbuf))
  write(cbuf,trim(fmt_char_int2)) 'isymm: ',isymm
  CALL logmessage('dump_setup',TRIM(cbuf))

  !-----------------------------
  !< Write case
  !-----------------------------
  WRITE(cbuf,trim(fmt_char_int2)) 'SIMCASE: ',simcase
  CALL logmessage('dump_setup',TRIM(cbuf))

  !-----------------------------------------------------------
  !< WRITE INFORMATION
  !-----------------------------------------------------------
  WRITE(cbuf,trim(fmt_char_real)) 'CUTOFF: ',cutoff
  CALL logmessage('dump_setup',TRIM(cbuf))
  WRITE(cbuf,trim(fmt_char_real)) 'SMOOTHINGLENGTH: ',smoothinglength
  CALL logmessage('dump_setup',TRIM(cbuf))
  WRITE(cbuf,trim(fmt_char_real)) 'PARTICLE SPACING : ',dpx
  CALL logmessage('dump_setup',TRIM(cbuf))
  WRITE(cbuf,trim(fmt_char_real)) 'VERLET_SKIN / DX : ',verlet_skin/dpx
  CALL logmessage('dump_setup',TRIM(cbuf))

  !----------------------------------------
  !< Format for output
  !----------------------------------------
  WRITE(fmt_local,'(A)') '(A,A,E15.5)'

  !----------------------------------------
  !< SOUNDSPEEDS
  !----------------------------------------
  WRITE(cbuf,fmt_char) 'Sound speeds: '
  DO i = 1,N_PHASES
    WRITE(cbuf,'(3A)') TRIM(cbuf),'  ',TRIM(phase(i)%name)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'Cs (vel): '
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%soundspeed(idx_sound_vel)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'Cs (visc): '
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%soundspeed(idx_sound_visc)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'Cs (force): '
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%soundspeed(idx_sound_force)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'Max. Cs (used):'
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%cmax
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'Mach number:'
  DO i = 1,N_PHASES

    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',abs(vel_ref)/phase(i)%cmax
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'Reynolds number:'
  DO i = 1,N_PHASES
    if (phase(i)%visc_dyn .gt. myeps) then
      WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%density*abs(vel_ref)*length_ref/(phase(i)%visc_dyn+myeps)
    else
      WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',0._mk
    endif
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  !---------------------------------------------------------
  ! Output pressure
  !---------------------------------------------------------
  WRITE(cbuf,fmt_char) 'Equation of state: '
  DO i = 1,N_PHASES
    WRITE(cbuf,'(3A)') TRIM(cbuf),'  ',TRIM(phase(i)%name)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'p0: '
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%press_ref
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'b (constant): '
  DO i = 1,N_PHASES
    write(cbuf,trim(fmt_local)) trim(cbuf),'  ',phase(i)%press_b
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  !---------------------------------------------------------
  ! OUTPUT TIMESTEPS
  !---------------------------------------------------------
  WRITE(cbuf,fmt_char) 'Timesteps: '
  DO i = 1,N_PHASES
    WRITE(cbuf,'(A,A,A)') TRIM(cbuf),'  ',TRIM(phase(i)%name)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'dt (vel): '
  DO i = 1,N_PHASES
    WRITE(cbuf,trim(fmt_local)) TRIM(cbuf),'  ',phase(i)%timestep(idx_dt_vel)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'dt (visc): '
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%timestep(idx_dt_visc)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,fmt_char) 'dt (force): '
  DO i = 1,N_PHASES
    WRITE(cbuf,fmt_local) TRIM(cbuf),'  ',phase(i)%timestep(idx_dt_force)
  ENDDO
  CALL logmessage('dump_setup',TRIM(cbuf))

  WRITE(cbuf,trim(fmt_char_real)) 'Min. timestep (used): ', dt_global
  CALL logmessage('dump_setup',TRIM(cbuf))

end subroutine dump_setup
