!--------------------------------------
!<  Subroutine to initialize simulation
!--------------------------------------
subroutine initialize(info)

  use mod_data_global
  use mod_data_sph
  use mod_data_ctrl
  use mod_data_physics
  use mod_io
  use mod_util
  use mod_comp, only : comp_constants

  use ppm_module_init
  use ppm_module_io
  use ppm_module_data, only : ppm_kind_single, ppm_kind_double, ppm_comm
  use ppm_module_ctrl

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif


  !-------------------------------------------
  !< Arguments
  !-------------------------------------------
  integer, intent(inout) :: info

  !-------------------------------------------
  !< Local variables
  !-------------------------------------------
  character(len=maxchar)   :: myname = 'initialize'
  character(len=maxchar)   :: mesg
  character(len=maxchar)   :: cbuf
  integer                 :: tolexp

#ifdef __MPI
  !-------------------------------------------------------------------------
  !  Initialise MPI
  !-------------------------------------------------------------------------
  info = 0
  comm = MPI_COMM_WORLD
  call MPI_Init(info)
  if (info .ne. 0) then
    write(*,*) 'FAILED TO INITIALIZE MPI. ABORTING!'
    stop
  ENDIF

  !-------------------------------------------------------------------------
  !  Get rank of processor
  !-------------------------------------------------------------------------
  CALL MPI_Comm_Rank(comm,rank,info)
  !-------------------------------------------------------------------------
  !  Get total number of processors
  !-------------------------------------------------------------------------
  CALL MPI_Comm_Size(comm,Nproc,info)
  !-------------------------------------------------------------------------
  !  Get name of processor
  !-------------------------------------------------------------------------
  CALL MPI_Get_Processor_Name(procname,len_procname,info)
#else
  !-------------------------------------------------------------------------
  !<  Initialise number of processors
  !-------------------------------------------------------------------------
  Nproc = 1
  rank  = 0
  info = HOSTNM(procname)
#endif

  !-------------------------------------------------------------------------
  !< Get input
  !-------------------------------------------------------------------------
  ppm_comm = comm
  call define_inputvariables
  !  !call disable_ctrl()
  call parse_args(info)
  if (info .ne. 0) then
    call finalize(1,info)
    stop
  endif

!  simcase = 1
!  N_PHASES=1
!  ndim=2
!  min_domain=0._mk
!  max_domain=21.e-3_mk
!  bcdef=1
!  bcvel=0._mk
!  bcvel(4,1) = 2.e-5_mk
!  cutoff=3.e-3
!  ratio_cutoff_smoothing = 3._mk
!  ratio_smoothing_dx = 1._mk
!  verlet_skin=0.1_mk
!  phase(1)%name='FLUID 1'
!  phase(1)%density=1000.
!  phase(1)%type = 2
!  phase(1)%visc_dyn = 1.e-3
!  phase(1)%gamma = 7._mk
!  phase(1)%delta_dens = 0.01_mk
!  phase(1)%press_b = 1._mk
!  phase(1)%output = 1
!  bodyforce = 0._mk
!  length_ref = 1.e-3_mk
!  vel_ref = 2.e-5_mk
!  density_ref = 1000._mk
!  choice_assign = 1
!  choice_decomp = 2
!  !  switch_showtimings = .true.
!  timesteppingscheme = 1
!  start_time = 0._mk
!  end_time = 2._mk
!  choice_output=1
!  !output_velocity = .true.
!  output_rank = .true.
!  freqlogoutput = 2
!  maxitsteps = 10000
!  sph_debug = 0
!  ppm_debug = 0
!  choice_restart = .false.

  !  print*, maxitsteps
  !  print*, phase(1)%press_b
  !  print*, phase(1)%delta_dens
  !  print*, phase(1)%gamma
  !  print*, phase(1)%visc_dyn
  !  print*, phase(1)%type
  !  print*, phase(1)%density
  !  print*, trim(phase(1)%name)
  !  print*, verlet_skin
  !  print*, ratio_smoothing_dx
  !  print*, ratio_cutoff_smoothing
  !  print*, bcvel
  !  print*, simcase
  !  print*, N_PHASES
  !  print*, ndim
  !  print*, min_domain
  !  print*, max_domain
  !  print*, bcdef
  !  print*, cutoff

  !-------------------------------------------------------------------------
  !< Get arguments of command line
  !-------------------------------------------------------------------------
  call get_arguments

  !-------------------------------------------------------------------------
  !< set defaults
  !-------------------------------------------------------------------------
  call defaults

  !-------------------------------------------------------------------------
  !  Open log file
  !-------------------------------------------------------------------------
  IF (rank.EQ.0) THEN
    call ppm_io_unused_unit(unit_log,info)
    OPEN(unit_log,FILE=TRIM(logfile),ACTION='WRITE',IOSTAT=info)
  ENDIF

  !-------------------------------------------------------------------------
  !  Print a first sign of life
  !-------------------------------------------------------------------------
  WRITE(mesg,'(2A)') 'PROGRAM STARTED ON: ',trim(procname)
  CALL logmessage(trim(myname),mesg)
  WRITE(mesg,'(2A)') 'using the command: ',trim(commandname)
  CALL logmessage(trim(myname),mesg)
  WRITE(mesg,'(A,I5,A)') 'RUNNING WITH ',Nproc,' Processors'
  CALL logmessage(trim(myname),mesg)

  !-------------------------------------------------------------------------
  !<  Initialise the ppm library
  !-------------------------------------------------------------------------
  call ppm_io_unused_unit(unit_ppmlog,info)
  tolexp = LOG10(EPSILON(1._mk))
  CALL ppm_init(dim=ndim,prec=mk,tolexp=tolexp,comm=comm,&
  debug=ppm_debug,info=info,logfile=unit_ppmlog,stderr=0,stdout=6)
  IF (info .NE. 0) call abortmessage(trim(myname),'ppm_init failed')

  !-------------------------------------------------------------------------
  !<  Check ctrl file data
  !-------------------------------------------------------------------------
  call check_ctrlfile
  call logmessage(trim(myname),'ctrlfile ok')

  !-------------------------------------------------------------------------
  !< Debug log file / proc log file
  !-------------------------------------------------------------------------
  if (sph_debug .gt. 0) then
    call ppm_io_unused_unit(unit_proclog,info)

    if (unit_proclog .gt. 255) then
      call abortmessage(trim(myname),'unit_proclog > 255!')
    endif
    open(unit_proclog,file=trim(proclogfile),action='write',iostat=info)
  endif

#ifdef __MPI
  !-------------------------------------------------------------------------
  !  Define MPI precision
  !-------------------------------------------------------------------------
  IF (mk .EQ. ppm_kind_single) THEN
    cbuf = 'Using single precision'
    MPPREC = MPI_REAL
  ELSEIF (mk .EQ. ppm_kind_double) THEN
    cbuf = 'Using double precision'
    MPPREC = MPI_DOUBLE_PRECISION
  ELSE
    WRITE(cbuf,'(A,I4)')'Unknown precision: ',MK
    call logmessage('sph_init',trim(cbuf))
    call finalize(0,info)
    stop
  ENDIF
  CALL logmessage('sph_init',TRIM(cbuf))
#endif

  !-------------------------------------------------------------------------
  !< Remove abortfile
  !---------------------------------------------------
  call remove_abortfile

  !-------------------------------------------------------------------------
  !<  Calc constants (timesteps, sound speed, reference pressure ,...
  !-------------------------------------------------------------------------
  call comp_constants

  !---------------------------------------------------
  !< Dump parameter -> logfile / stdout
  !---------------------------------------------------
  call dump_setup

  return

end subroutine initialize
