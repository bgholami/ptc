!--------------------------------
!< Module that contains all routines for utilities
!--------------------------------
module mod_util

  use mod_data_global, only : mk, t_realarray, t_intarray

  implicit none

  INTERFACE trace
    MODULE PROCEDURE trace_matrix
  END INTERFACE

  interface abs_vec
    module procedure absolut_vector
  end interface

  interface reset_ptr_to_list
    module procedure reset_ptr_to_list_r
    module procedure reset_ptr_to_list_i
  end interface
  interface add_ptr_to_list
    module procedure add_ptr_to_list_r
    module procedure add_ptr_to_list_i
  end interface
  interface get_ptr_to_list
    module procedure get_ptr_to_list_r
    module procedure get_ptr_to_list_i
  end interface

  interface allreduce_mean
    module procedure allreduce_mean_r8
  end interface

  interface minmaxmeanvar
    module procedure minmaxmeanvar_r8
!    module procedure minmaxmeanvar_i8
  end interface

  interface timer_reset
    module procedure timer_reset_array
    module procedure timer_reset_single
  end interface

contains

  !< Reset pointer list -> REAL
  subroutine reset_ptr_to_list_r(orig_array)
    use mod_data_global

    implicit none

    integer :: i

    type(t_realarray), dimension(:), pointer :: orig_array
    forall (i=1:ldamaparray) orig_array(i)%ptr_data => null()
    lastrealmaparray = 0
  end subroutine reset_ptr_to_list_r

  !< Reset pointer list -> INTEGER
  subroutine reset_ptr_to_list_i(orig_array)
    use mod_data_global

    implicit none

    integer :: i

    type(t_intarray), dimension(:), pointer :: orig_array
    forall (i=1:ldamaparray) orig_array(i)%ptr_data => null()
    lastintmaparray = 0
  end subroutine reset_ptr_to_list_i

  !< ADD something to pointer list -> REAL
  subroutine add_ptr_to_list_r(orig_array,add_element)
    use mod_data_global

    implicit none

    type(t_realarray), dimension(:), pointer :: orig_array
    real(mk), dimension(:,:), pointer :: add_element

    if (lastrealmaparray .lt. ldamaparray) then
      lastrealmaparray = lastrealmaparray + 1
      orig_array(lastrealmaparray)%ptr_data => add_element
      orig_array(lastrealmaparray)%lda = size(add_element,1)

    else
      print*,__FILE__,__LINE__
      stop
    endif
  end subroutine add_ptr_to_list_r

  !< ADD something to pointer list -> INTEGER
  subroutine add_ptr_to_list_i(orig_array,add_element)
    use mod_data_global

    implicit none

    type(t_intarray), dimension(:), pointer :: orig_array
    integer, dimension(:,:), pointer :: add_element

    if (lastintmaparray .lt. ldamaparray) then
      lastintmaparray = lastintmaparray + 1
      orig_array(lastintmaparray)%ptr_data => add_element
      orig_array(lastintmaparray)%lda = size(add_element,1)
    else
      print*,__FILE__,__LINE__
      stop
    endif
  end subroutine add_ptr_to_list_i


  !< GET new local pointer from list -> REAL
  subroutine get_ptr_to_list_r(orig_array,get_element)
    use mod_data_global

    implicit none

    type(t_realarray), dimension(:), pointer :: orig_array
    real(mk), dimension(:,:), pointer :: get_element

    if (lastrealmaparray .gt. 0) then
      get_element => orig_array(lastrealmaparray)%ptr_data
      orig_array(lastrealmaparray)%lda = 0
      orig_array(lastrealmaparray)%ptr_data => null()
      lastrealmaparray = lastrealmaparray - 1
    else
      print*,__FILE__,__LINE__
      stop
    endif
  end subroutine get_ptr_to_list_r

  !< GET new local pointer from list -> INTEGER
  subroutine get_ptr_to_list_i(orig_array,get_element)
    use mod_data_global

    implicit none

    type(t_intarray), dimension(:), pointer :: orig_array
    integer, dimension(:,:), pointer :: get_element

    if (lastintmaparray .gt. 0) then
      get_element => orig_array(lastintmaparray)%ptr_data
      orig_array(lastintmaparray)%lda = 0
      orig_array(lastintmaparray)%ptr_data => null()
      lastintmaparray = lastintmaparray - 1
    else
      print*,__FILE__,__LINE__
      stop
    endif
  end subroutine get_ptr_to_list_i


  subroutine timer_start(timer,info)
    use mod_data_global
    use ppm_module_time

    implicit none

    !< Arguments
    type(t_timer), intent(inout) :: timer
    integer, intent(out)                 :: info

    call ppm_time(timer%tstart,info)
    if (info .ne. 0) then
      print*, __FILE__,__LINE__
      stop
    endif
    timer%calls = timer%calls + 1

  end subroutine timer_start

  subroutine timer_stop(timer,info)
    use mod_data_global
    use ppm_module_time

    implicit none

    !< Arguments
    type(t_timer), intent(inout) :: timer
    integer, intent(out)                 :: info

    call ppm_time(timer%tend,info)
    if (info .ne. 0) then
      print*, __FILE__,__LINE__
      stop
    endif
    timer%delta = timer%tend - timer%tstart
    timer%t_sum = timer%t_sum + timer%delta

  end subroutine timer_stop

  subroutine timer_reset_array(timer,info)
    use mod_data_global
    use ppm_module_time

    implicit none

    !< Arguments
    type(t_timer), dimension(:), intent(inout) :: timer
    integer, intent(out)                 :: info
    !< Local vars
    integer :: ntimer
    integer :: j

    ntimer = size(timer,1)
    do j = 1, ntimer
      timer(j)%tstart = 0._mk
      timer(j)%tend   = 0._mk
      timer(j)%delta  = 0._mk
      timer(j)%t_sum   = 0._mk
      timer(j)%calls  = 0
    enddo
    info = 0
  end subroutine timer_reset_array

  subroutine timer_reset_single(timer,info)
    use mod_data_global
    use ppm_module_time

    implicit none

    !< Arguments
    type(t_timer), intent(inout) :: timer
    integer, intent(out)                 :: info

    timer%tstart = 0._mk
    timer%tend   = 0._mk
    timer%delta  = 0._mk
    timer%t_sum   = 0._mk
    timer%calls  = 0
    info = 0

  end subroutine timer_reset_single


  subroutine allreduce_mean_r8(realin,realout)
    use mod_data_global

    implicit none
#ifdef __MPI
    include 'mpif.h'
#endif
    real(8), intent(in)  :: realin
    real(8), intent(out) :: realout

    integer :: ierr

    realout = 0.

#ifdef __MPI
    call MPI_ALLREDUCE(realin,realout,1,MPPREC,MPI_SUM,comm,ierr)
    realout = realout / Nproc
#else
    realout = realin
#endif

  end subroutine allreduce_mean_r8


#include "src/util/check_char_length.inc" 
#include "src/util/util_uppercase.inc"
#include "src/util/remove_abortfile.inc"
#include "src/util/check_abortfile.inc"
#include "src/util/finalize.inc"
#include "src/util/get_arguments.inc"
#include "src/util/resetarrays.inc"
#include "src/util/reshapearrays.inc"
#include "src/util/minmaxmeanvar_r8.inc"

  !----------------------------------------------------------------------
  !< TRACE OF A MATRIX
  !----------------------------------------------------------------------
  REAL(mk) FUNCTION trace_matrix(inputmatrix)

    REAL(mk), DIMENSION(__DIM,__DIM), INTENT(IN) :: inputmatrix
    INTEGER :: i

    trace_matrix = 0._MK
    DO i = 1, __DIM
    trace_matrix = trace_matrix + inputmatrix(i,i)
  ENDDO

END FUNCTION trace_matrix

!----------------------------------------------------------------------
! LENGTH OF A VECTOR
!----------------------------------------------------------------------
REAL(MK) FUNCTION absolut_vector(inputvec)
  USE mod_data_global

  REAL(MK), DIMENSION(__DIM), INTENT(IN) :: inputvec
  INTEGER :: i

  absolut_vector = 0._MK
  DO i = 1, __DIM
  absolut_vector = absolut_vector + inputvec(i)**2._MK
ENDDO

absolut_vector = sqrt(absolut_vector)

END FUNCTION absolut_vector

end module mod_util
