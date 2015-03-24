!--------------------------------
!< Routine that defines all input variables
!--------------------------------
subroutine define_inputvariables

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_physics

  use ppm_module_ctrl
  use mod_io, only : logmessage

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------

  !-----------------------------
  !< Geometry
  !-----------------------------
  call arg_group('Geometry')

  call arg(ndim,'DIMENSIONS',long_flag='--dim',&
	ctrl_name = 'ndim', default = __DIM, min = 2, max = __DIM)

	call arg(simcase,'SIMCASE',long_flag='--simcase',&
	ctrl_name='simcase',default=1,min=1,max=param_max_simcase)

#if __DIM==2
	call arg(min_domain,'MIN_DOMAIN',long_flag='--min_domain',&
	ctrl_name='min_domain',default=(/-2._mk, -2._mk/))
	call arg(max_domain,'MAX_DOMAIN',long_flag='--max_domain',&
	ctrl_name='max_domain',default=(/2._mk, 2._mk/))
#elif __DIM==3
	call arg(min_domain,'MIN_DOMAIN',long_flag='--min_domain',&
	ctrl_name='min_domain',default=(/-2._mk,-2._mk,-2._mk/))
	call arg(max_domain,'MAX_DOMAIN',long_flag='--max_domain',&
	ctrl_name='max_domain',default=(/2._mk,2._mk,2._mk/))
#endif

#if __DIM==2
  call arg(bcdef(1:2),'BC-X: min/max',long_flag='--boundarycondition-x',&
  ctrl_name='bc-x',default=(/1,1/),min=1,max=10)
  call arg(bcdef(3:4),'BC-Y: min/max',long_flag='--boundarycondition-y',&
  ctrl_name='bc-y',default=(/7,7/),min=1,max=10)
#elif __DIM==3
  call arg(bcdef(1:2),'BC-X: min/max',long_flag='--boundarycondition-x',&
  ctrl_name='bc-x',default=(/1,1/),min=1,max=10)
  call arg(bcdef(3:4),'BC-Y: min/max',long_flag='--boundarycondition-y',&
  ctrl_name='bc-y',default=(/1,1/),min=1,max=10)
  call arg(bcdef(5:6),'BC-Z: min/max',long_flag='--boundarycondition-z',&
  ctrl_name='bc-z',default=(/7,7/),min=1,max=10)
#endif

#if __DIM==2
  call arg(bcvel(1,1:__DIM),'WALLVEL_X_MIN',long_flag='--wallvel-x-min',&
  ctrl_name='wallvel-x-min',default=(/0._mk,0._mk/))
  call arg(bcvel(2,1:__DIM),'WALLVEL_X_MAX',long_flag='--wallvel-x-max',&
  ctrl_name='wallvel-x-max',default=(/0._mk,0._mk/))
  call arg(bcvel(3,1:__DIM),'WALLVEL_Y_MIN',long_flag='--wallvel-y-min',&
  ctrl_name='wallvel-y-min',default=(/-2._mk,0._mk/))
  call arg(bcvel(4,1:__DIM),'WALLVEL_Y_MAX',long_flag='--wallvel-y-max',&
  ctrl_name='wallvel-y-max',default=(/2._mk,0._mk/))
#elif __DIM==3
  call arg(bcvel(1,1:__DIM),'WALLVEL_X_MIN',long_flag='--wallvel-x-min',&
  ctrl_name='wallvel-x-min',default=(/0._mk,0._mk,0._mk/))
  call arg(bcvel(2,1:__DIM),'WALLVEL_X_MAX',long_flag='--wallvel-x-max',&
  ctrl_name='wallvel-x-max',default=(/0._mk,0._mk,0._mk/))
  call arg(bcvel(3,1:__DIM),'WALLVEL_Y_MIN',long_flag='--wallvel-y-min',&
  ctrl_name='wallvel-y-min',default=(/0._mk,0._mk,0._mk/))
  call arg(bcvel(4,1:__DIM),'WALLVEL_Y_MAX',long_flag='--wallvel-y-max',&
  ctrl_name='wallvel-y-max',default=(/0._mk,0._mk,0._mk/))
  call arg(bcvel(5,1:__DIM),'WALLVEL_Z_MIN',long_flag='--wallvel-z-min',&
  ctrl_name='wallvel-z-min',default=(/-2._mk,0._mk,0._mk/))
  call arg(bcvel(6,1:__DIM),'WALLVEL_Z_MAX',long_flag='--wallvel-z-max',&
  ctrl_name='wallvel-z-max',default=(/2._mk,0._mk,0._mk/))
#endif

  !-----------------------------
  !< Physics
  !-----------------------------
  call arg_group('Physics')

#if __DIM==2
  call arg(bodyforce,'BODYFORCE',flag='-g',long_flag='--bodyforce',&
  ctrl_name='bodyforce',default=(/0._mk,0._mk/))
#else
  call arg(bodyforce,'BODYFORCE',flag='-g',long_flag='--bodyforce',&
  ctrl_name='bodyforce',default=(/0._mk,0._mk,0._mk/))
#endif
  call arg(length_ref,'LENGTH_REF',long_flag='--length-reference',&
  ctrl_name='length_ref',default=1._mk,min=0._mk)
  call arg(density_ref,'DENSITY_REF',long_flag='--density-reference',&
  ctrl_name='density_ref',default=1._mk,min=0._mk)
  call arg(vel_ref,'VELOCITY_REF',long_flag='--velocity-reference',&
  ctrl_name='vel_ref',default=1._mk,min=0._mk)
  call arg(N_PHASES,'N_PHASES',long_flag='--n-phases',&
  ctrl_name='N_PHASES',default=1,min=1,max=param_max_n_phases)

  !-----------------------------
  !< Phases
  !< hard coded to use parse_args...
  !-----------------------------
  call arg_group('PHASE 1')

  call arg(phase(1)%name,'PHASE_NAME (1)',long_flag='--phase1-name',&
  ctrl_name='NAME',default='FLUID 1')
  call arg(phase(1)%type,'PHASE_TYPE (1)',long_flag='--phase1-type',&
  ctrl_name='TYPE',default=param_type_fluid)
  call arg(phase(1)%density,'PHASE_DENSITY (1)',long_flag='--phase1-density',&
  ctrl_name='DENSITY',default=1._mk)
  call arg(phase(1)%visc_dyn,'PHASE_VISCOSITY (dynamic) (1)',long_flag='--phase1-dynvisc',&
  ctrl_name='VISC_DYN',default=1._mk)
  call arg(phase(1)%gamma,'PHASE_GAMMA (1)',long_flag='--phase1-gamma',&
  ctrl_name='GAMMA',default=7._mk,min=0._mk)
  call arg(phase(1)%delta_dens,'PHASE_DELTADENS (1)',long_flag='--phase1-deltadens',&
  ctrl_name='DELTA_DENS',default=0.01_mk,min=0._mk)
  call arg(phase(1)%press_b,'BACKGROUNDPRESSURE (*p0) (1)',long_flag='--pressre-background',&
  ctrl_name='PRESS_BACKGROUND',default=1._mk,min=0._mk)
  call arg(phase(1)%output,'PHASE_OUTPUT(1)',long_flag='--phase1-output',&
  ctrl_name='OUTPUT',default=1,min=0,max=1)


  !-----------------------------
  !< SPH parameter
  !-----------------------------
  CALL arg_group('SPH parameter')

  CALL arg(cutoff, 'cutoff',long_flag = '--cutoff',&
  ctrl_name = 'cutoff',default = 0.12_mk, min=0._mk)

  call arg(ratio_cutoff_smoothing,'RATIO-CUTOFF-SMOOTHINGLENGTH',long_flag='--ratio-cutoff_smoothing',&
  ctrl_name='ratio_cutoff_smoothing',default=3._mk)

  call arg(ratio_smoothing_dx,'RATIO-SMOOTHING-DX',long_flag='--ratio_smoothing_dx',&
  ctrl_name='ratio_smoothing_dx',default=1._mk)

  call arg(verlet_skinfactor,'SKINFACTOR',long_flag='--verlet-skinfactor',&
  ctrl_name='verlet_skinfactor',default=0.1_mk,min=0._mk,max=3._mk)


  !-----------------------------
  !< Numerics
  !-----------------------------
  call arg_group('NUMERICS')

  call arg(choice_assign,'SUB-TO-PROC ASSIGN',long_flag='--choice_assign',&
  ctrl_name='choice_assign',default=param_assign_internal,max=param_max_assign)
  call arg(choice_decomp,'DOMAIN DECOMP.',long_flag='--choice_decomp',&
  ctrl_name='choice_decomp',default=param_decomp_bisection,max=param_max_decomp)

  !-----------------------------
  !< SWITCHES
  !-----------------------------
  call arg_group('SWITCHES')

  call arg(switch_showtimings,'SHOW TIMINGS',long_flag='--switch_showtimings',&
  ctrl_name='switch_showtimings',default=.false.)

  !-----------------------------
  !< Timestepping
  !-----------------------------
  call arg_group('TIMESTEPPING')
  call arg(timesteppingscheme,'TIME-SCHEME',long_flag='--timesteppingscheme',&
  ctrl_name='timesteppingscheme',default=param_timestepping_kickdriftkick,min=1,max=param_max_timestepping)
  call arg(start_time,'START-TIME',long_flag='--starttime',&
  ctrl_name='Start-Time',default=0._mk)
  call arg(end_time,'END-TIME',long_flag='--endtime',&
  ctrl_name='End-Time',default=3._mk)

  !-----------------------------
  !< Output switch
  !-----------------------------
  call arg_group('Output')

  call arg(choice_output,'Switch output type',long_flag='--choice_output',&
  ctrl_name = 'choice_output', default = 0, min = 0,max = param_max_choice_output)

  call arg(timeoutput,'Output times',long_flag='--time-output',&
  ctrl_name='timeoutput',default=1._mk)

  call arg(choice_output_nondim,'Output non-dim.',long_flag='--output_nondim',&
  ctrl_name = 'output_nondim',default=.false.)

  call arg(output_velocity,'output_velocity',long_flag='--output_velocity',&
  ctrl_name='output_velocity',default=.true.)

  call arg(output_h,'output_h',long_flag='--output-h',&
  ctrl_name='output_h',default=.false.)

  call arg(output_density,'output_density',long_flag='--output_density',&
  ctrl_name='output_density',default=.false.)

  call arg(output_pressure,'output_pressure',long_flag='--output_pressure',&
  ctrl_name='output_pressure',default=.false.)

  call arg(output_mass,'output_mass',long_flag='--output_mass',&
  ctrl_name='output_mass',default=.false.)

  call arg(output_volume,'output_volume',long_flag='--output_volume',&
  ctrl_name='output_volume',default=.false.)

  call arg(output_phaseidx,'output_phase',long_flag='--output_phase',&
  ctrl_name='output_phase',default=.false.)

  call arg(output_id,'output_id',long_flag='--output_id',&
  ctrl_name='output_id',default=.false.)

  call arg(output_rank,'output_rank',long_flag='--output_rank',&
  ctrl_name='output_rank',default=.true.)

  !-----------------------------
  !< Log
  !-----------------------------
  call arg_group('LOGGING')

  call arg(freqlogoutput,'FREQLOGOUTPUT',long_flag='--freqlogoutput',&
  ctrl_name='freqlogoutput',default=5,min=1)
  call arg(maxitsteps,'MAXITSTEPS',long_flag='--maxitsteps',&
  ctrl_name='MAXITSTEPS',default=10,min=1)
  call arg(sph_debug,'SPH_DEBUG',long_flag='--sph_debug',&
  ctrl_name='sph_debug',default=0,min=0,max=param_max_sph_debug)
  call arg(ppm_debug,'PPM_DEBUG',long_flag='--ppm_debug',&
  ctrl_name='ppm_debug',default=0,min=0,max=param_max_ppm_debug)

  !-----------------------------
  !< Restart
  !-----------------------------
  call arg_group('RESTART')

  CALL arg(choice_restart,'choice_restart',long_flag='--choice_restart',&
  ctrl_name='choice_restart',default= param_restart_none,min= 0,max=param_max_restart)

  call arg(restart_reset_steps,'reset_steps',long_flag='--restart_reset_steps',&
  ctrl_name='restart_reset_steps',default=.false.)
  call arg(restart_reset_velocity,'reset_velocity',long_flag='--restart_reset_velocity',&
  ctrl_name='restart_reset_velocity',default=.false.)

  CALL arg(freqrestart,'freqrestart',long_flag='--freqrestart',&
  ctrl_name = 'freqrestart', default= 10000000)

end subroutine define_inputvariables
