!--------------------------------
!< Module that contains all routines required for initialization
!--------------------------------
module mod_init


  implicit none
contains

#include "src/init/create_particles.inc"
#include "src/init/place_particle.inc"
#include "src/init/initialize.inc"
#include "src/init/init_state.inc"
#include "src/init/init_velocity.inc"
#include "src/init/init_defaults.inc"
#include "src/init/dump_setup.inc"
#include "src/init/check_ctrlfile.inc"
#include "src/init/define_inputvariables.inc"

end module mod_init
