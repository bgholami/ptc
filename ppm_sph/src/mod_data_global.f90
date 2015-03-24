!============================================
!< Module containing all global data
!============================================
module mod_data_global

  use ppm_module_data, only : ppm_char, ppm_kind_double

  use iso_c_binding, only : c_int64_t, c_int8_t, c_double

  implicit none

  !< Precision
  integer, parameter :: mk = ppm_kind_double
  !< Max. Character length
  integer, parameter :: maxchar = ppm_char


  !------------------------------
  !< MPI variables
  !------------------------------
  integer           :: comm             !< communicator
  integer           :: rank             !< ID of actual proc
  integer           :: Nproc            !< Number of procs
  character(maxchar):: procname         !< Name of proc
  integer           :: len_procname     !< len of procname
  integer           :: MPPREC


  !------------------------------
  ! < program
  !------------------------------
  character(maxchar)    :: prgname
  character(maxchar)    :: commandname

  !------------------------------
  ! < files and folder
  !------------------------------
  character(maxchar)    :: ctrlfile       !< control file
  character(maxchar)    :: thispath       !< current path
  character(maxchar)    :: logfile        !< log file
  character(maxchar)    :: abortfile      !< abortfile
  character(maxchar)    :: outputposfile  !< outputposfile for special output
  character(maxchar)    :: proclogfile    !< logfile for every proc.
  character(maxchar)    :: outputfile, ghostfile !< filenames of output
  character(maxchar)    :: restartfile, oldrestartfile !< filename of restart
  character(maxchar)    :: conservationfile !< output conservation props

  !-----------------------------------
  !< file units
  !-----------------------------------
  integer               :: unit_log                   = 10 !< set by ppm internally
  integer               :: unit_writeascii            = 10 !< set by ppm internally
  integer               :: unit_debug                 = 10 !< set by ppm internally
  integer               :: unit_proclog               = 10 !< set by ppm internally
  integer               :: unit_ppmlog                = 10 !< set by ppm internally

  !-----------------------------------
  !< Debug level
  !-----------------------------------
  integer, parameter    :: param_max_sph_debug = 2
  integer, parameter    :: param_max_ppm_debug = 2
  integer               :: sph_debug
  integer               :: ppm_debug

  !-----------------------------------
  !< Number of limiting criteria
  !< Viscousity / Velocity / Body force / Surface tension / Diffusion
  !-----------------------------------
  INTEGER, PARAMETER :: NDTCRITERIA  = 3
  INTEGER, PARAMETER :: idx_dt_vel   = 1
  INTEGER, PARAMETER :: idx_dt_visc  = 2
  INTEGER, PARAMETER :: idx_dt_force = 3

  !-----------------------------------
  !< Viscousity / Velocity / Body force / Surface tension
  !-----------------------------------
  INTEGER, PARAMETER :: NCSOUNDCRITERIA   = 3
  INTEGER, PARAMETER :: idx_sound_vel   = 1
  INTEGER, PARAMETER :: idx_sound_visc  = 2
  INTEGER, PARAMETER :: idx_sound_force = 3

  !-----------------------------------
  !< Format specifies
  !-----------------------------------
  character(len=8)              :: fmt_char_int2      = '(A40,I2)'
  character(len=9)              :: fmt_char_int10     = '(A40,I10)'
  character(len=12)             :: fmt_char_real      = '(A40,F20.10)'
  character(len=11)             :: fmt_char_exp       = '(A40,E17.7)'
  character(len=7)              :: fmt_char_char      = '(A40,A)'
  character(len=7)              :: fmt_char_boolean   = '(A40,L)'
  character(len=5)              :: fmt_char           = '(A40)'
  character(len=14)             :: fmt_char_real_char = '(A40,F20.10,A)'

  !---------------------------------
  !< PI
  !--------------------------------
  real(mk)                      :: pi, pi2, pi_inv

  !--------------------------------
  !< Tolerance level
  !--------------------------------
  real(mk)                      :: myeps

  !---------------------------------
  !< Global time marching
  !---------------------------------
  real(mk)                      :: current_time
  real(mk)                      :: dt_global !< Global timestep
  real(mk)                      :: rdt_global !< 1/dt_global
  integer                       :: nsteps !< steps to be done
  integer                       :: current_step
  integer                       :: restart_step_start !< step at restart

  !------------------------------------------
  !< Output details
  !-------------------------------------------
  !< For conservation output and specialoutput
  integer                             :: ldaoutvec_conservation
  integer                             :: ldaoutvec_ascii

  integer                             :: ldaspecialoutvec
  real(mk), dimension(:), allocatable :: specialoutvec

  !------------------------------------
  !< PPM details
  !------------------------------------
  integer                           :: topo_id
  real(mk), dimension(:), pointer   :: cost => null()
  real(mk), dimension(:,:), pointer :: min_sub => null(), max_sub => null()
  integer                           :: nsubs
  integer, dimension(:), pointer    :: sub2proc => null()
  real(mk), dimension(:), pointer   :: pcost => null()
  integer, dimension(:,:), pointer  :: ident => null()
  integer                           :: nident

  !< Symmetric interactions ?
  integer                           :: isymm = 1 !< 1= symmetry, 0 = no symmetry
  logical                           :: symmetry = .true.

  !< Rebalance
  logical                           :: rebalance = .false.
  logical                           :: lflush = .false.
  integer                           :: count_rebalancing = 0

  !< New verlet list
  logical                           :: make_new_vlist = .false.
  integer                           :: count_newvlist = 0

  !< Map array for variable mapping
  integer, parameter                        :: ldamaparray = 20 !< -> not more than 20 vecs to communicate...
  integer                                   :: lastintmaparray = 0
  integer                                   :: lastrealmaparray = 0

  !-----------------------------
  !< Timing results
  !-----------------------------
  integer                           :: step_wctstep_s = 0
  integer                           :: step_lastrebalance

  !-----------------------------
  !< Type definition for timer
  !< 1:tic, 2:toc, 3:toc-tic, 4:sum(toc-tic), 5:calls
  !-----------------------------
  type t_timer
    real(8)                             :: tstart
    real(8)                             :: tend
    real(8)                             :: delta !> = tend - tstart
    real(8)                             :: t_sum !> = t_sum + tend - tstart
    integer                             :: calls !> number of calls for this timer
  end type

  type(t_timer)                        :: wct_prestep !> everything before real calc starts
  type(t_timer)                        :: wct_step    !> time of whole step routine
  type(t_timer)                        :: wct_output  !> time of output

  integer, parameter                   :: ntimings = 15   !> subtimings

  integer, parameter                    :: stepidx      = 1
  integer, parameter                    :: forceidx     = 2
  integer, parameter                    :: ploop1idx    = 3
  integer, parameter                    :: velidx       = 4
  integer, parameter                    :: posidx       = 5
  integer, parameter                    :: densidx      = 6
  integer, parameter                    :: propidx      = 7
  integer, parameter                    :: pushidx      = 8
  integer, parameter                    :: vlistidx     = 9
  integer, parameter                    :: bcidx        = 10
  integer, parameter                    :: resetidx     = 11
  integer, parameter                    :: putidx       = 12
  integer, parameter                    :: checkidx     = 13
  integer, parameter                    :: ghostidx     = 14
  integer, parameter                    :: bodyforceidx = 15


  type(t_timer), dimension(ntimings)          :: wct_timestep !> detailed timing of one timestep (without output)
  character(len=maxchar), dimension(ntimings) :: wct_name

  !------------------------------------
  !< Vars used for minmaxmeanvar -> utils
  !------------------------------------
  real(8), dimension(:), allocatable    :: gather_nproc
  integer, parameter                    :: n_procbal = 4
  integer, parameter                    :: minidx = 1
  integer, parameter                    :: maxidx = 2
  integer, parameter                    :: meanidx= 3
  integer, parameter                    :: varidx = 4

  !------------------------------------
  !< Setup parameter: init_particles
  !------------------------------------
  integer, parameter :: param_max_simcase                   = 2
  integer, parameter :: param_simcase_lattice               = 1

  !--------------------------
  !< Setup parameter: density equation choice
  !--------------------------
  INTEGER, PARAMETER :: param_density_summation = 1

  !--------------------------
  !< Setup parameter: output choices
  !-----------------------------
  integer, parameter :: param_max_choice_output         = 3
  integer, parameter :: param_output_none               = 0
  integer, parameter :: param_output_ascii              = 1
  integer, parameter :: param_output_unformatted        = 2
  integer, parameter :: param_output_ppmdbg             = 3

  !--------------------------
  !< Setup parameter: restart choices
  !-----------------------------
  integer, parameter  :: param_max_restart    = 2
  integer, parameter  :: param_restart_none   = 0
  integer, parameter  :: param_restart_ascii  = 1
  integer, parameter  :: param_restart_binary = 2

  !-----------------------------------------
  !< Particle type definitions
  !------------------------------------------
  INTEGER, PARAMETER :: param_n_types = 2
  INTEGER, PARAMETER :: param_type_solidboundary = 1
  INTEGER, PARAMETER :: param_type_fluid         = 2

  !----------------------------------------
  !< Setup parameter: Choice time-stepping scheme
  !-------------------------------------------------
  integer, parameter  :: param_max_timestepping           = 1
  INTEGER, PARAMETER :: param_timestepping_kickdriftkick  = 1

  !-----------------------------------------------
  !< Boundary definition
  !-----------------------------------------------
  integer, parameter :: param_bcdef_periodic      = 1
  integer, parameter :: param_bcdef_openboundary  = 2
  integer, parameter :: param_bcdef_symmetry      = 3
  integer, parameter :: param_bcdef_neumann       = 4
  integer, parameter :: param_bcdef_dirichlet     = 5
  integer, parameter :: param_bcdef_robin         = 6
  integer, parameter :: param_bcdef_mirrorwall    = 7
  integer, parameter :: param_bcdef_injection     = 8
  integer, parameter :: param_bcdef_solidwall     = 9
  integer, parameter :: param_bcdef_antisymmetry  = 10

  !---------------------------------------------------
  !< Prescribe velocity
  !---------------------------------------------------
  integer, parameter :: param_prescribemotion_sin      = 1
  integer, parameter :: param_prescribemotion_acc      = 2
  integer, parameter :: param_prescribemotion_jump     = 3
  integer, parameter :: param_prescribemotion_velocity = 4

  !---------------------------------------------------
  !< Map types
  !---------------------------------------------------
  integer, parameter :: param_maptype_partial         = 1
  integer, parameter :: param_maptype_global          = 2
  integer, parameter :: param_maptype_getghost        = 3
  integer, parameter :: param_maptype_ghostupdate     = 4
  integer, parameter :: param_maptype_ghostupdatenoxp = 5

  !---------------------------------------------------
  !< Subdomain-to-processor assignment schemes
  !---------------------------------------------------
  integer, parameter :: param_max_assign              = 6
  integer, parameter :: param_assign_internal         = 1
  integer, parameter :: param_assign_nodal_cut        = 2
  integer, parameter :: param_assign_nodal_comm       = 3
  integer, parameter :: param_assign_dual_cut         = 4
  integer, parameter :: param_assign_dual_comm        = 5
  integer, parameter :: param_assign_user_defined     = 6

  !---------------------------------------------------
  !< Domain decomposition techniques
  !---------------------------------------------------
  integer, parameter :: param_max_decomp          = 13
  INTEGER, PARAMETER :: param_decomp_tree         = 1
  INTEGER, PARAMETER :: param_decomp_pruned_cell  = 2
  INTEGER, PARAMETER :: param_decomp_bisection    = 3
  INTEGER, PARAMETER :: param_decomp_xpencil      = 4
  INTEGER, PARAMETER :: param_decomp_ypencil      = 5
  INTEGER, PARAMETER :: param_decomp_zpencil      = 6
  INTEGER, PARAMETER :: param_decomp_cuboid       = 7
  INTEGER, PARAMETER :: param_decomp_user_defined = 8
  INTEGER, PARAMETER :: param_decomp_xy_slab      = 10
  INTEGER, PARAMETER :: param_decomp_xz_slab      = 11
  INTEGER, PARAMETER :: param_decomp_yz_slab      = 12
  INTEGER, PARAMETER :: param_decomp_cartesian    = 13

  !---------------------------------------------------
  !< Derived type for array of pointers... -> used for mapping
  !---------------------------------------------------
  type t_realarray
    real(mk), dimension(:,:), pointer :: ptr_data => null()
    integer                           :: lda
  end type
  type t_intarray
    integer, dimension(:,:), pointer  :: ptr_data => null()
    integer                           :: lda
  end type


end module mod_data_global

