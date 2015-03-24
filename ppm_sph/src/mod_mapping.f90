!-------------
!< Module contains all mapping routines
!-----------------------------
module mod_mapping
    implicit none

    interface makemaparray
      module procedure makemaparray_general
    end interface

    interface getmaparray
      module procedure getmaparray_general
    end interface

contains

#include "src/mapping/map_ghost_put.inc"

#include "src/mapping/check_newvlist.inc"

#include "src/mapping/makemaparray_general.inc"
#include "src/mapping/getmaparray_general.inc"

#include "src/mapping/mapping.inc"

end module mod_mapping
