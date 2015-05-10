      MODULE Class_Tracers
        !----------------------------------------------------
        ! Class	      :	Tracers
        !----------------------------------------------------
        !
        !  Purpose    :
        !> \brief       Variables and corresponding operations
        !>              for Tracers information.
        !>	   	
        !>              The variable members are public 
        !  Remarks    :
        !
        !  References :
        !
        !  Revisions  : V0.1 28.04.2011 
        !
        !----------------------------------------------------
        !  Author     : Babak Gholami
        !  Contact    : babak.gholami@aer.mw.tum.de
        !  Dr. Marco Ellero's Emmy Noether Group,
        !  Prof. Dr. N. Adams' Chair of Aerodynamics,
        !  Faculty of Mechanical Engineering,
        !  Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
#if __TRACER 

        USE Class_Particles        
        USE Class_Kernel

        
        IMPLICIT NONE
        SAVE
        
        INTERFACE tracers_new
           MODULE PROCEDURE particles_init
        END INTERFACE

        INTERFACE tracers_display_parameters
           MODULE PROCEDURE particles_display_parameters
        END INTERFACE

        INTERFACE tracers_decompose_global
           MODULE PROCEDURE particles_decompose_global
        END INTERFACE

        INTERFACE tracers_get_num_part_real
           MODULE PROCEDURE particles_get_num_part_real
        END INTERFACE        

        INTERFACE tracers_get_x
           MODULE PROCEDURE particles_get_x
        END INTERFACE        

        INTERFACE tracers_get_v
           MODULE PROCEDURE particles_get_v
        END INTERFACE        

        INTERFACE tracers_get_id
           MODULE PROCEDURE particles_get_id
        END INTERFACE 

        INTERFACE tracers_get_num_dim
           MODULE PROCEDURE particles_get_num_dim
        END INTERFACE   

        INTERFACE tracers_get_num_id
           MODULE PROCEDURE particles_get_num_id
        END INTERFACE      
       
        INTERFACE tracers_integrate_position
           MODULE PROCEDURE particles_integrate_position
        END INTERFACE        

        INTERFACE tracers_integrate_velocity
           MODULE PROCEDURE particles_integrate_velocity
        END INTERFACE        

        INTERFACE tracers_decompose_partial
           MODULE PROCEDURE particles_decompose_partial
        END INTERFACE

        INTERFACE tracers_finalize
           MODULE PROCEDURE particles_finalize
        END INTERFACE     

      CONTAINS

#include "src/tracers/tracers_interpolate_velocity.F90" 
#include "src/tracers/tracers_compute_force.F90"
#include "src/tracers/tracers_compute_wall_distance.F90"
#include "src/tracers/tracers_compute_force_collision.F90"
#include "src/tracers/tracers_delete_tracers.F90"
#include "src/tracers/tracers_eval.F90"
#include "src/tracers/tracers_integrate.F90"
#include "src/tracers/tracers_integrate_step.F90"
#include "src/tracers/tracers_adjust_tracers.F90"
#include "src/tracers/tracers_near_wall_velocity.F90"
#include "src/tracers/tracers_coupling_extraction.F90"
#include "src/tracers/tracers_coupling_insertion.F90"
#include "src/tracers/tracers_generate_new.F90"
#include "src/tracers/tracers_interface_create.F90"
#include "src/tracers/tracers_lu.F90"
#include "src/tracers/tracers_decompose_ring.F90"
#include "src/tracers/tracers_adjust_tracers_special.F90"

#endif

      END MODULE Class_Tracers
      
      
