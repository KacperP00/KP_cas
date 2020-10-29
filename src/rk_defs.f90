module rk_defs
  use precision
  use math

  implicit none

  type rktvd_solver_t
     character(len=32) :: scheme = 'noname'
     integer :: stage = 0, order = 0, niter, nzo, neq, nsc
     real(WP), dimension(:,:), pointer :: coeff(:,:) => null()
     real(WP), dimension(:,:,:), pointer :: RK => null(), Res => null()
   end type rktvd_solver_t

contains

  subroutine allocate_rktvd(solver)
    implicit none

    type(rktvd_solver_t), pointer, intent(inout) :: solver
    ! ---------------------------------

    allocate(solver%coeff(solver%stage,2*solver%order))
    allocate(solver%RK(solver%stage,solver%neq+solver%nsc,solver%nzo))
    allocate(solver%Res(solver%stage,solver%neq+solver%nsc,solver%nzo))

  end subroutine allocate_rktvd

  subroutine deallocate_rktvd(solver)
    implicit none

    type(rktvd_solver_t), pointer, intent(inout) :: solver
    ! ---------------------------------

    deallocate(solver%coeff)
    deallocate(solver%RK)
    deallocate(solver%Res)

  end subroutine deallocate_rktvd  

end module rk_defs
