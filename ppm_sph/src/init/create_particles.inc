!--------------------------------
!< Subroutine to create particles
!--------------------------------
subroutine create_particles

  use mod_data_global
  use mod_data_sph
  use mod_data_ctrl
  use mod_data_physics
  use mod_data_prtl

  use mod_io
  use ppm_module_write

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
  character(maxchar)          :: myname = 'create_particles'
  integer                     :: Ntot
  integer                     :: i
  character(len=maxchar)      :: cbuf

  integer                     :: sum_Ntotal
  integer                     :: Npart_local
  integer                     :: array_size
  integer                     :: ierr, info
  integer                     :: k,j
  real(mk)                    :: x,y
#if __DIM==3
  real(mk)                    :: z
#endif


  !-----------------------------
  !< Estimate number of particles
  !-----------------------------
  do i = 1, ndim
    Np(i)  = int((max_domain(i) - min_domain(i))/dpx+0.5_mk)
    min_compbox(i) = min_domain(i) - myeps
    max_compbox(i) = min_domain(i) + dpx*Np(i) -myeps
  enddo

  Ntot = 1
  do i = 1, ndim
    Ntot = Np(i)*Ntot
  enddo

  write(cbuf,fmt_char_real) 'Initial dx: ',dpx
  call logmessage('create_particles',trim(cbuf))

  write(cbuf,fmt_char_int10) 'Initial points in x: ',Np(1)
  call logmessage('create_particles',trim(cbuf))
  write(cbuf,fmt_char_int10) 'Initial points in y: ',Np(2)
  call logmessage('create_particles',trim(cbuf))
#if __DIM==3
  write(cbuf,fmt_char_int10) 'Initial points in z: ',Np(3)
  call logmessage('create_particles',trim(cbuf))
#endif

  write(cbuf,fmt_char_real) 'Comp. box size in x: ',max_compbox(1) - min_compbox(1)
  call logmessage('create_particles',trim(cbuf))
  write(cbuf,fmt_char_real) 'Comp. box size in y: ',max_compbox(2) - min_compbox(2)
  call logmessage('create_particles',trim(cbuf))
#if __DIM==3
  write(cbuf,fmt_char_real) 'Comp. box size in z: ',max_compbox(3) - min_compbox(3)
  call logmessage('create_particles',trim(cbuf))
#endif

  write(cbuf,fmt_char_int10) 'Total number of lattice points: ',Ntot
  call logmessage('create_particles',trim(cbuf))



  !----------------------------------------
  !< Nullify and allocate arrays
  !----------------------------------------
  allocate(Ntotal(N_PHASES))
  Ntotal = 0

  array_size = int(Ntot / Nproc + 2.5)

  !< Position
  allocate(xp(NXP,array_size))
  xp = 0._mk

  !< pdata
  allocate(pdata(NPDATA,array_size))
  pdata = 0._mk

  !< integer attributes
  allocate(ap(NAP,array_size))
  ap = 0

  !-------------------------------
  !< Create particles
  !-------------------------------
  Npart = 0
  Npart_local = 0

  !-----------------------------------------
  !<  INITIAL PARTICLE POSITIONS
  !-----------------------------------------
  SELECT CASE(simcase)

    CASE (param_simcase_lattice)
      ! Initialize particles on regular lattice0
#include "inc/init_particles_lattice.inc"

    CASE DEFAULT
      ! Wrong simcase
      CALL abortmessage('create_particles','Wrong input of simcase!')
  END SELECT

  !---------------------------
  !< Total number of particles for each phase and globally
  !---------------------------
  Ntotal = 0
  do k = 1, Npart_local
    Ntotal(ap(phaseidx,k)) = Ntotal(ap(phaseidx,k)) + 1
  enddo

  do k = 1, N_PHASES
    sumNtotal = Ntotal(k)
#ifdef  __MPI
    call MPI_ALLREDUCE(sumNtotal,sum_Ntotal,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    sumNtotal = sum_Ntotal
#endif
    Ntotal(k) = sumNtotal
  enddo
  sumNtotal = sum(Ntotal)

  WRITE(cbuf,'(A)') '__________________________________________'
  CALL logmessage('create_particles',TRIM(cbuf))
  WRITE(cbuf,fmt_char_int10) 'Global number of particles: ', sumNtotal
  CALL logmessage('create_particles',TRIM(cbuf))
  WRITE(cbuf,'(A)') '__________________________________________'
  CALL logmessage('create_particles',TRIM(cbuf))

  WRITE(cbuf,'(A,I5,A,I15)') 'Initialized particles locally on Proc. ',rank,': ',Npart_local
  call logmessage('create_particles',trim(cbuf))
  if (rank .gt. 0) call ppm_write(rank,'create_particles',trim(cbuf),info)

  !<Every proc to its logfile
  if (sph_debug .gt. 0) call ppm_write(rank,'create_particles',trim(cbuf),info,unit_proclog)

  !------------------------
  !< Abort if no particles
  !------------------------
  if (sumNtotal .lt. 1) call abortmessage(myname,'No particles initialized...')

  !------------------------
  !< Now use Npart on each proc
  !--------------------------
  Npart = Npart_local

end subroutine create_particles
