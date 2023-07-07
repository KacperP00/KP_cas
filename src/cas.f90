!--------------------------------------------------------------------------------!
! CAS: Cross-sectionally Averaged Spray                                          !
!                                                                                ! 
! This program solves cross-sectionally averaged one-D multiphase equations      !
! describing spray formation                                                     !
!                                                                                !
! Author: Abhishek Y. Deshmukh                                                   !
! Major contributors: Avijit Saha, Carsten Giefer                                !
!--------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com or a.saha at itv.rwth-aachen.de           !
!--------------------------------------------------------------------------------!

program cas

  use parser
  use spray_defs
  use spray_func
#ifdef MDUC_MPI
  use MPI
#endif
  implicit none

  !include 'mpif.h'

  !--------------------------------------
  type(spray_t), pointer :: spray
  type(pc_t), pointer :: pc
  real(WP) :: tstart, tend
  !--------------------------------------

  integer :: ierr, ntasks, myid

#ifdef MDUC_MPI
  ! I n i t i a l i z a t i o n of the MPI environment
  call MPI_INIT ( ierr )
  call MPI_COMM_RANK ( MPI_COMM_WORLD , myid , ierr )
  call MPI_COMM_SIZE ( MPI_COMM_WORLD , ntasks , ierr )
#endif

  ! Allocate spray object
  call allocate_spray(spray)

#ifdef MDUC_MPI
  spray%nranks = ntasks
  spray%irank = myid
  spray%comm = MPI_COMM_WORLD
#endif

  ! Get input filename from command line if present
  call getarg(1,spray%inp_fname)

  ! If input filename not provided, choose default filename
  if (spray%inp_fname == '') then
     write(*,*) 'Input filename not provided. Reading default cas.in ...'
     spray%inp_fname = "cas.in"
  end if

  ! Read parameters from input file
  call read_inputs(spray)

  ! Initialize spray
  call init_spray(spray)

  call cpu_time(tstart)

  ! Solve the equations
  call run_spray(spray) 

  call cpu_time(tend)

  write(*,*) 'Runtime: ',tend-tstart,' seconds'

#ifdef MDUC_MPI
  if(myid .eq. 0) then
#endif

  ! Write output
  call write_output(spray,spray%step,spray%ndtime)

#ifdef MDUC_MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$  if(spray%twoD) then
!!$     call spray_combust_poststep_mrif_2D(spray)
!!$  else
  call spray_combust_poststep_mrif(spray)
!!$  end if

  end if
#endif

  ! Deallocate spray object
  call deallocate_spray(spray)

#ifdef MDUC_MPI
  call MPI_FINALIZE ( ierr )
#endif

end program cas
