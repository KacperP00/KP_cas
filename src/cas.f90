!----------------------------------------!
! CAS: Cross-sectionally Averaged Spray  !
! This program solves one-D multiphase   !
! equations describing spray formation   !
!----------------------------------------!

program cas

  use parser
  use spray_defs
  use spray_func
!!$  use output
  
  implicit none

  !--------------------------------------
  type(spray_t), pointer :: spray
  type(pc_t), pointer :: pc
  character(len=128) :: fname
  !--------------------------------------
  
  ! Get input filename from command line if present
  call getarg(1,fname)

  ! If input filename not provided, choose default filename
  if (fname == '') then
     write(*,*) 'Input filename not provided. Reading default cas.in ...'
     fname = "cas.in"
  end if

  ! Allocate spray object
  call allocate_spray(spray)

  ! Read parameters from input file
  call read_inputs(spray,fname)

  ! Initialize spray
  call init_spray(spray)

  ! Solve the equations
  call run_spray(spray) 

  ! Write output
  call write_output(spray,spray%step,spray%ndtime)

  ! Deallocate spray object
  call deallocate_spray(spray)
end program cas
