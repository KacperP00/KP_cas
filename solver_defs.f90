module solver_defs
  use precision
  use math
  use rk_defs

  implicit none

  type solver_t

     ! Scheme
     character(len=32) :: scheme = 'LF'

     ! Order
     integer :: order = 1

     ! RK Solver
     type(rk_solver_t), pointer :: rk
     
     ! Convective scheme
     real(WP), dimension(:,:), pointer :: divc, weno5p, weno5m

     ! Optimal weights for weno5
     real(WP) :: a0, a1, a2

     ! Stencils for weno5
     real(WP), dimension(:), pointer :: S0_p, S1_p, S2_p, S0_m, S1_m, S2_m, f_p, f_m, f_t
     
     ! State and flux vectors
     real(WP), dimension(:,:), pointer :: W, Wold, F, S, Res

     ! Flux, Residual, Wave speeds(liquid and gas phase)
     real(WP), dimension(:), pointer :: Flux, alpha_l, alpha_g    

  end type solver_t

contains

  subroutine allocate_solver(solver,nzo)
    implicit none
    
    ! ---------------------------------
    type(solver_t), pointer, intent(inout) :: solver
    integer, intent(in) :: nzo
    ! ---------------------------------

!!$    nullify(solver%rk, &
!!$            solver%divc, &
!!$            solver%weno5p, solver%weno5m, solver%S0_p, solver%S1_p, solver%S2_p, &
!!$            solver%S0_m, solver%S1_m, solver%S2_m, solver%f_p, solver%f_m, solver%f_t, &
!!$            solver%W, solver%Wold, solver%F, solver%S, solver%Res, &
!!$            solver%Flux, solver%alpha_l, solver%alpha_g)

    allocate(solver%W(8,nzo)); solver%W = 0.0_WP
    allocate(solver%Wold(8,nzo)); solver%Wold = 0.0_WP
    allocate(solver%F(8,nzo)); solver%F = 0.0_WP
    allocate(solver%S(8,nzo)); solver%S = 0.0_WP
    allocate(solver%Res(8,nzo)); solver%Res = 0.0_WP

    allocate(solver%Flux(nzo)); solver%Flux = 0.0_WP
    allocate(solver%alpha_l(nzo)); solver%alpha_l = 0.0_WP
    allocate(solver%alpha_g(nzo)); solver%alpha_g = 0.0_WP

    allocate(solver%rk%RK(8,nzo)); solver%rk%RK = 0.0_WP
    allocate(solver%rk%dRK(8,nzo)); solver%rk%dRK = 0.0_WP

  end subroutine allocate_solver

  subroutine deallocate_solver(solver)
    implicit none
    
    ! ---------------------------------
    type(solver_t), pointer, intent(inout) :: solver

    ! ---------------------------------

    call deallocate_rk(solver%rk)

    deallocate(solver%divc)

    if(associated(solver%weno5p)) then
       deallocate(solver%weno5p, solver%weno5m, solver%S0_p, solver%S1_p, solver%S2_p, solver%S0_m, solver%S1_m, solver%S2_m, solver%f_p, solver%f_m, solver%f_t)
    end if

    deallocate(solver%W, solver%Wold, solver%F, solver%S)
    
    deallocate(solver%Flux, solver%Res, solver%alpha_l, solver%alpha_g)
    
    deallocate(solver%rk%RK, solver%rk%dRK)

    deallocate(solver)

  end subroutine deallocate_solver

end module solver_defs
