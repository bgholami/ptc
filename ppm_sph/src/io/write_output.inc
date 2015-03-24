!--------------------------------
!< Subroutine for output
!--------------------------------
subroutine write_output

  use mod_data_global
  use mod_data_ctrl
  use mod_data_sph
  use mod_data_prtl
  use mod_data_physics

  use mod_util

  use ppm_module_map_part_util
  use ppm_module_io
  use ppm_module_data

  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  !character(len=maxchar)      :: myname = 'write_output'
  !integer                     :: info
  character(len=maxchar)      :: cbuf

  !character(len=maxchar)      :: filename, iofmt
  !real(mk)    :: znorm1, znorm2, znorm3, norm_dev
  !real(mk)    :: znorm1_all, znorm2_all, znorm3_all, norm_dev_all
  !integer     :: ierr
  !integer     :: k, j, i, ip
  !real(mk)    :: g_analyt, z
  !real(mk)    :: sumarea, sum_area,r0

  !-----------------------------
  !< Choose output
  !-----------------------------
  select case (choice_output)

    case (param_output_none)
      ! nothing to do

    case (param_output_ascii)
      call write_output_ascii(outputfile,1,Npart)
      if (sph_debug .gt. 0) call write_output_ascii(ghostfile,Npart+1,Mpart)

    case (param_output_ppmdbg)
      call write_ppmdbg(current_step)

    case default
      write(cbuf,fmt_char_int2) 'Wrong choice_output: ',choice_output
      call abortmessage('write_output',trim(cbuf))
  end select

end subroutine write_output
