!--------------------------------------------------!
! CATS: Cross-sectionally Averaged Turbulent Spray !
! This program solves cross-sectionally averaged   !
! one-D multiphase equations describing spray      !
! formation along with turbulence modeling         !
!                                                  !
! Author: Abhishek Y. Deshmukh                     !
! Contact: a.deshmukh@itv.rwth-aachen.de           !
!--------------------------------------------------!

program cats

  use parser
  use spray_defs
  use spray_func
!!$  use output
  implicit none

  !--------------------------------------
  type(spray_t), pointer :: spray
  type(pc_t), pointer :: pc
  real(WP) :: tstart, tend
  !--------------------------------------

  ! Allocate spray object
  call allocate_spray(spray)

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

  ! Write output
  call write_output(spray,spray%step,spray%ndtime)

  ! Deallocate spray object
  call deallocate_spray(spray)

end program cats
