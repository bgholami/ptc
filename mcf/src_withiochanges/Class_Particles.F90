      MODULE Class_Particles
        !----------------------------------------------------
        ! Class	      :	Particles
        !----------------------------------------------------
        !
        !  Purpose    :
        !> \brief       Variables and corresponding operations
        !>              for Particles infomation.
        !>	   	
        !>              The variable memebers are public 
        !  Remarks    :
        !
        !  References :
        !
        !  Revisions  : V0.1 03.03.2009 
        !
        !----------------------------------------------------
        !  Author     : Xin Bian
        !  Contact    : xin.bian@aer.mw.tum.de
        !  Dr. Marco Ellero's Emmy Noether Group,
        !  Prof. Dr. N. Adams' Chair of Aerodynamics,
        !  Faculty of Mechanical Engineering,
        !  Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        USE mcf_header
        USE Class_Debug
        USE Class_Tool
        
        USE Class_Control
        USE Class_Physics
        USE Class_Colloid
        USE Class_Rhs
        USE Class_StateEquation
        USE Class_Technique
        USE Class_Kernel
	USE Class_Random
        
        
        IMPLICIT NONE
        SAVE
        
        TYPE Particles
           
           TYPE(Control), POINTER               :: ctrl
           TYPE(Physics), POINTER               :: phys
           TYPE(Rhs), POINTER                   :: rhs
           TYPE(StateEquation), POINTER         :: stateEquation
           TYPE(Kernel), POINTER                :: kern
           TYPE(Technique), POINTER             :: tech
           LOGICAL                              :: pp_interact_cc
           LOGICAL                              :: pp_interact_cw
           
           !-------------------------------------------------
           ! x     : position
           ! v     : velocity
           ! rho   : mass / number density
           ! m     : mass
           ! p     : pressure
           ! id    : (1,:) particle ID
           !         (2,:) species  ID, i.e.,
           !          0 for fluid, > 0 for colloids
           !          < 0 for walls
           ! f     : force per unit mass / acceleration 
           ! f_bp  : background pressure force per unit 
           !         mass / acceleration (Transport velocity formulation)
           ! fa_min: minimum force acceleration.
           ! fa_max: maximum force acceleration.
           !           
           ! u     : thermal energy.
           ! au    : acceleration of thermal energy.
           !
           ! Remark :
           !         1) Using PPM library, if there are
           !         ghosts particles, they are always
           !         in the end of the data struture
           !          continously.
           !-------------------------------------------------
           
           INTEGER                              :: num_dim
           REAL(MK)                             :: h
           REAL(MK), DIMENSION(:,:), POINTER    :: x
           REAL(MK), DIMENSION(:,:), POINTER    :: v
           REAL(MK), DIMENSION(:), POINTER      :: rho
           REAL(MK), DIMENSION(:), POINTER      :: rho_norm
           REAL(MK)                             :: rho_min
           REAL(MK)                             :: rho_max
           REAL(MK), DIMENSION(:), POINTER      :: m
           REAL(MK), DIMENSION(:), POINTER      :: p
           INTEGER, DIMENSION(:,:), POINTER     :: id
           REAL(MK), DIMENSION(:,:), POINTER    :: f 
           REAL(MK), DIMENSION(:,:), POINTER    :: f_bp
           REAL(MK)                             :: fa_min
           REAL(MK)                             :: fa_max
           REAL(MK)                             :: dt_f
           
           REAL(MK), DIMENSION(:), POINTER      :: u
           REAL(MK), DIMENSION(:), POINTER      :: au
           
           !-------------------------------------------------
           ! Quantites more related to non-Newtonian fluid.
           !
           ! vgt   : velocity gradient tensor
           !
           ! evgt  : egenvector dynamics vgt
           ! eval  : egenvalues
           ! aeval : acceleration of egenvalues
           ! evec  : egenvectors
           ! aevec : acceleration of egenvectors           
           ! 
           ! ct    : conformation tensor
           ! act   : acceleation of conformation tensor

           ! Remark :
           !         2) Since PPM doesn't support,
           !         tensor/matrix for particles,
           !         i.e., it doesn't support 
           !         rank=3 dimensional data,
           !         we have to use dim**2 array
           !         for each particle, if matrix
           !         is needed. 
           !         In context of Fortran language,
           !         column dominating convention
           !         is used, i.e., index change
           !         first by rows then columns.

           !
           !
           ! vgt,evgt, ct or act's
           ! matrix notation, 2D(3D):
           !  Cxx   Cxy  (Cxz)
           !  Cyx   Cyy  (Cyz)
           ! (Czx) (Czy) (Czz)
           !
           ! array notation of 2D in order:
           ! Cxx, Cyx, Cxy, Cyy;
           ! array notation of 3D in order:
           ! Cxx, Cyx, Czx, Cxy, Cyy, Czy,
           ! Cxz, Cyz, Czz.
           !
           ! egenvalue 2D(3D):
           ! e1,e2(,e3).
           !
           ! egenvector matrix notation, 2D(3D):
           !
           !  ev1_x   ev2_x  (ev3_x)
           !  ev1_y   ev2_y  (ev3_y)
           ! (ev1_z) (ev2_z) (ev3_z)
           !
           ! 2D array notation in order :
           ! ev1_x, ev1_y, ev2_x, ev2_y
           ! 3D array notation in order :
           ! ev1_x, ev1_y, ev1_z, ev2_x, ev2_y, ev2_z,
           ! ev3_x, ev3_y, ev3_z
           !
           ! However, pressure tensor(pt) is not needed
           ! to communicate between processors,
           ! therefore is is matrxi notation.
           !-------------------------------------------------
           
           REAL(MK), DIMENSION(:,:), POINTER    :: vgt
           REAL(MK), DIMENSION(:,:), POINTER    :: evgt
           REAL(MK), DIMENSION(:,:), POINTER    :: eval
           REAL(MK), DIMENSION(:,:), POINTER    :: aeval
           REAL(MK), DIMENSION(:,:), POINTER    :: evec
           REAL(MK), DIMENSION(:,:), POINTER    :: aevec
           REAL(MK), DIMENSION(:,:), POINTER    :: ct
           REAL(MK), DIMENSION(:,:), POINTER    :: act
           REAL(MK),DIMENSION(:,:,:),POINTER    :: pt
           
           ! Num of ids for each particle
           INTEGER                              :: num_id
           ! Particle id index
           INTEGER                              :: pid_idx
           ! Spieces id index
           INTEGER                              :: sid_idx
           ! Boundary id index (in case of inflow/outflow boundary)
           INTEGER                              :: bid_idx
            
           
           !----------------------------------------------------
           ! real :  real particles on local process;
           ! all  :  all particles including ghost;
           ! ghost:  ghost particles.
           !
           ! all = real + ghost
           !
           ! fluid   : number of real fluid particles;
           ! sym     : number of symmetry boundary particles.
           !           ( are all ghosts particles)
           ! wall_sym: number of wall boundary particles created
           !           by PPM using symmetry/mirror particles.
           !           ( are all ghosts particles)
           ! wall_solid: inital number of wall boundary real particles created
           !           by MCF using solid particles.
           ! wall_solid_real : number of wall boundary ghost particles
           !           created by MCF using solid particles.
           !           ( are real particles)
           ! wall_solid_ghost : number of wall boundary ghost particles
           !           created by MCF using solid particles.
           !           ( are ghosts particles)
           ! le      : number of Lees-Edwards boundary particles.
           !           ( are ghost particles)
           !
           ! shear   : number of wall boundary particles,
           !           created by ppm using symmetry/mirror,
           !           or by MCF as solid wall, or number of
           !           Lees-Edwards boundary particles.
           !
           !
           ! colloid : number of colloid boundary particles.
           !
           !
           ! real = fluid + wall_solid + colloid
           !-----------------------------------------------------
           
           INTEGER                              :: num_part_real
           INTEGER                              :: num_part_all
           INTEGER                              :: num_part_ghost
           
           INTEGER                              :: num_part_fluid
           INTEGER                              :: num_part_sym
           INTEGER                              :: num_part_wall_sym
           INTEGER                              :: num_part_wall_solid
           INTEGER                              :: num_part_wall_solid_real
           INTEGER                              :: num_part_wall_solid_ghost
           INTEGER                              :: num_part_le
           INTEGER                              :: num_part_colloid  

           
           !----------------------------------------------------
           ! Record the index and species ID of each :
           !
           ! symmetry boundary particle;
           ! ( are all ghost particles)
           ! wall_using symmetry boundary particles;
           ! ( are all ghost particles)
           ! wall using solid boundary particles;
           ! ( part of are real, part of are ghost particles)
           ! Lees-Edwards boundary particles.
           ! ( are all ghost particles)
           ! colloid boundary particle;
           ! ( part of are real, part of are ghost particles)
           !----------------------------------------------------
           
           INTEGER, DIMENSION(:,:), POINTER     :: part_sym_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_walL_sym_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_wall_solid_real_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_wall_solid_ghost_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_le_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_colloid_list
           
           TYPE(Tool)                           :: tool
           TYPE(Random), POINTER                :: random           

#if __TRACER 

           !----------------------------------------------------
           ! Additional parameters concerning tracers 
           ! (used by SPH fluid particles):
           ! 
           ! c_tracers     : number of tracers on a SPH particle
           ! dn            : evolution of c_tracers
           !----------------------------------------------------
           REAL(MK), DIMENSION(:)    , POINTER  :: c_tracers  
           REAL(MK), DIMENSION(:)    , POINTER  :: dn
           
           !----------------------------------------------------
           ! Additional parameters that ONLY Tracers use:
           !
           ! radius        : radius of the tracers
           ! rho_ratio     : density ratio (rho_particles/rho_tracers)
           ! c_factor      : integer factor to convert concentration 
           !                 to tracer number
           ! dp            : tracers diffusion coefficient
           ! tau_p         : momentum response time
           ! m_p           : mass of tracers
           ! near_wall     : near-wall region distance
           ! last_local_id : the last local tracer ID generated on this process
           ! x_interface   : position of interface
           ! vo_interface  : old velocity of interface
           ! vc_interface  : current velocity of interface
           ! vn_interface  : new velocity of interface
           ! local_interface 
           !               : indicates if an interface point is local or not
           ! insertion_list 
           !               : number of tracers to be inserted into SPH particles
           !                 at this interface point
           ! num_interface : number of interface points
           ! num_tracer_deposited 
           !               : current number of deposited tracers in the near-wall 
           !                 control volume
           ! h_p           : particle/tracer wall distance
           ! n_vec         : particle/tracer wall normal
           ! facet         : particle/tracer corresponding wall facer
           ! col_i         : particle/tracer corresponding wall colloid
           ! adhesion parameters according to Piper1998
           ! nbond         : adhesion; number of bonds
           ! K0            : adhesion; (= mr * ml * Ac * K0a)
           ! baa           : adhesion; binding affinity a  
           ! bab           : adhesion; binding affinity b
           ! bac           : adhesion; binding affinity c
           ! bad           : adhesion; binding affinity d
           !----------------------------------------------------
                      
           REAL(MK)                             :: radius
           REAL(MK)                             :: rho_ratio           
           INTEGER                              :: c_factor
           REAL(MK)                             :: dp
           REAL(MK)                             :: tau_p
           REAL(MK)                             :: m_p
           REAL(MK)                             :: near_wall
           INTEGER                              :: last_local_id
           REAL(MK), DIMENSION(:,:,:), POINTER  :: x_interface
           REAL(MK), DIMENSION(:,:,:), POINTER  :: vo_interface  
           REAL(MK), DIMENSION(:,:,:), POINTER  :: vc_interface  
           REAL(MK), DIMENSION(:,:,:), POINTER  :: vn_interface
           LOGICAL , DIMENSION(:,:)  , POINTER  :: local_interface
           INTEGER , DIMENSION(:,:)  , POINTER  :: insertion_list
           INTEGER , DIMENSION(:)    , POINTER  :: num_interface
           INTEGER , DIMENSION(:,:)  , POINTER  :: num_tracer_deposited 

           REAL(MK), DIMENSION(:)    , POINTER  :: h_p 
           REAL(MK), DIMENSION(:,:)  , POINTER  :: n_vec
           INTEGER , DIMENSION(:)    , POINTER  :: facet
           INTEGER , DIMENSION(:)    , POINTER  :: col_i  

           INTEGER                              :: nbond
           REAL(MK)                             :: K0  
           REAL(MK)                             :: baa
           REAL(MK)                             :: bab
           REAL(MK)                             :: bac
           REAL(MK)                             :: bad
           
#endif
           
        END TYPE Particles
        
        
        INTERFACE particles_new
           MODULE PROCEDURE particles_init
        END INTERFACE

        
      CONTAINS
        
#include "src/particles/particles_new.F90"
#include "src/particles/particles_finalize.F90"
#include "src/particles/particles_get.F90"
#include "src/particles/particles_set.F90"
#include "src/particles/particles_reset.F90"
#include "src/particles/particles_init_global_inter.F90"
#include "src/particles/particles_init_global_exter.F90"
#include "src/particles/particles_init_partial_inter.F90"
#include "src/particles/particles_init_partial_inter_adjust.F90"
#include "src/particles/particles_init_partial_exter.F90"
#include "src/particles/particles_set_colloid_on_lattice.F90"
#include "src/particles/particles_decompose_global.F90"
#include "src/particles/particles_decompose_partial.F90"
#include "src/particles/particles_map_ghost_get.F90"
#include "src/particles/particles_map_ghost_put.F90"
#include "src/particles/particles_compute_mass.F90"
#include "src/particles/particles_compute_density.F90"
#include "src/particles/particles_normalize_density.F90"
#include "src/particles/particles_find_density_extreme.F90"
#include "src/particles/particles_compute_pressure.F90"
#include "src/particles/particles_compute_interaction.F90"
#include "src/particles/particles_find_force_extreme.F90"
#include "src/particles/particles_compute_dt_f.F90"
#include "src/particles/particles_apply_body_force.F90"
#include "src/particles/particles_compute_act.F90"
#include "src/particles/particles_compute_evgt.F90"
#include "src/particles/particles_compute_aeval.F90"
#include "src/particles/particles_compute_aevec.F90"
#include "src/particles/particles_compute_ct.F90"
#include "src/particles/particles_compute_pressure_tensor.F90"
#include "src/particles/particles_collect_colloid_interaction.F90"
#include "src/particles/particles_collect_boundary_interaction.F90"
#include "src/particles/particles_integrate_position.F90"
#include "src/particles/particles_integrate_velocity.F90"
#include "src/particles/particles_integrate_colloid_position.F90"
#include "src/particles/particles_integrate_boundary_position.F90"
#include "src/particles/particles_adjust_particles.F90"
#include "src/particles/particles_create_delete_particles.F90"
#include "src/particles/particles_integrate_eval.F90"
#include "src/particles/particles_integrate_evec.F90"
#include "src/particles/particles_integrate_ct.F90"
#include "src/particles/particles_integrate_potential_energy.F90"
#include "src/particles/particles_set_boundary_ghost_id.F90"
#include "src/particles/particles_set_boundary_velocity.F90"
#include "src/particles/particles_set_boundary_ghost_velocity.F90"
#include "src/particles/particles_reset_boundary_velocity.F90"
#include "src/particles/particles_reset_boundary_interaction.F90"
#include "src/particles/particles_reset_boundary_ghost_interaction.F90"
#include "src/particles/particles_set_colloid_velocity.F90"
#include "src/particles/particles_reset_colloid_velocity.F90"
#include "src/particles/particles_reset_colloid_interaction.F90"
#include "src/particles/particles_set_flow_developed.F90"
#if __TRACER 
#include "src/particles/particles_evolve_ctracers.F90"
#endif

      END MODULE Class_Particles
      
      
