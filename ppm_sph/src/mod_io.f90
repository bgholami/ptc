!-----------------------------------------------
!< Module that contains all output routines
!-----------------------------------------------
module mod_io
  implicit none

  interface logmessage
    module procedure logmessage_2char
  end interface

  interface abortmessage
    module procedure abortmessage_2char
  end interface

  interface debugmessage
    module procedure debugmessage_2char
  end interface

contains

#include "src/io/write_output.inc"
#include "src/io/write_output_ascii.inc"
#include "src/io/write_ppmdbg.inc"

#include "src/io/logmessage_2char.inc"
#include "src/io/abortmessage_2char.inc"
#include "src/io/debugmessage_2char.inc"

#include "src/io/write_logoutput.inc"
#include "src/io/show_timings.inc"

#include "src/io/write_vlistresult.inc"

end module mod_io
