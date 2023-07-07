module solver_func
  use precision
  use math
  use solver_defs
  use spray_defs
  use rk_defs
  use rk_func

  implicit none

contains

  subroutine solver_init(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    call rk_init(spray)

    call computeDivc1(spray)

    select case(spray%solver%scheme)
    case ('WENO3')

       ! Compute Residual
       call init_weno3(spray)

    case ('WENO5')

       ! Compute Residual
       call init_weno5(spray)

    case default

    end select

  end subroutine solver_init

  subroutine solver_run(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    if (spray%solver%strang) then
       ! Solver with Strang splitting
       call solverStrangSplit(spray)       
    else       
       call solverNormal(spray)
    end if
    
  end subroutine solver_run

  subroutine solverStrangSplit(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    ! Symmetric Strang splitting
    
    ! Step 1: Solve 1/2 time step for convective part
    select case(spray%solver%scheme)
    case ('UPWIND1')

       ! Compute Residual
       call computeResidual_UPWIND1(spray)       

    case ('LF')

       ! Compute Residual
       call computeResidual_LF(spray)

    case ('WENO3')

       ! Compute Residual
       call computeResidual_WENO3(spray)

    case ('WENO5')

       ! Compute Residual
       call computeResidual_WENO5(spray)

    case default

    end select

    ! Update solution
    call updateSolution(spray)

    !call updatePrimVars(spray)
    call updateFlowVariables(spray)

    !call buildSourceVector(spray)


    ! Step 2: Add source terms 1 time step
    call addSource(spray)

    ! Update solution
    call updateSolution(spray)

    !call updatePrimVars(spray)
    call updateFlowVariables(spray)

    call buildFluxVector(spray)

    call buildStateVector(spray)

    ! Step 3: Solve 1/2 time step for convective part
    select case(spray%solver%scheme)
    case ('UPWIND1')

       ! Compute Residual
       call computeResidual_UPWIND1(spray)       

    case ('LF')

       ! Compute Residual
       call computeResidual_LF(spray)

    case ('WENO3')

       ! Compute Residual
       call computeResidual_WENO3(spray)

    case ('WENO5')

       ! Compute Residual
       call computeResidual_WENO5(spray)

    case default

    end select

    ! Update solution
    call updateSolution(spray)
    
  end subroutine solverStrangSplit

  subroutine solverNormal(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    select case(spray%solver%scheme)
    case ('UPWIND1')

       ! Compute Residual
       call computeResidual_UPWIND1(spray)       

    case ('LF')

       ! Compute Residual
       call computeResidual_LF(spray)

    case ('WENO3')

       ! Compute Residual
       call computeResidual_WENO3(spray)

    case ('WENO5')

       ! Compute Residual
       call computeResidual_WENO5(spray)

    case default

    end select

    ! Update solution
    call updateSolution(spray)

  end subroutine solverNormal

  ! Compute first order divergence operator
  subroutine computeDivc1(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    real(WP), dimension(:,:), pointer :: divc => null()

    allocate(divc(0:1,spray%kmino:spray%kmaxo))

    spray%solver%divc => divc

    divc(0,:) = -1.0_WP/spray%dz
    divc(1,:) = 1.0_WP/spray%dz

  end subroutine computeDivc1

  subroutine init_weno3(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(solver_t), pointer :: solver=>null()

    solver => spray%solver

    ! Allocate stencils
    allocate(solver%S0_p(-2:+0))
    allocate(solver%S1_p(-2:+0))
    allocate(solver%S0_m(-1:+1))
    allocate(solver%S1_m(-1:+1))

    ! Allocate WENO5 operators
    allocate(solver%weno3p(-2:0,spray%kmino:spray%kmaxo))
    allocate(solver%weno3m(-1:1,spray%kmino:spray%kmaxo))

    ! Allocate flux vectors
    allocate(solver%f_p(spray%kmino:spray%kmaxo))
    allocate(solver%f_m(spray%kmino:spray%kmaxo))
    allocate(solver%f_t(spray%kmino:spray%kmaxo))

    ! Interpolation at the faces
    solver%S0_p = 0.0_WP; solver%S0_m = 0.0_WP
    solver%S1_p = 0.0_WP; solver%S1_m = 0.0_WP

    solver%S0_p(-2) = real(-1.0/2.0,WP); solver%S0_p(-1) = real( 3.0/2.0,WP)
    solver%S0_m( 0) = real( 3.0/2.0,WP); solver%S0_m(+1) = real(-1.0/2.0,WP)

    solver%S1_p(-1) = real( 1.0/2.0,WP); solver%S1_p( 0) = real(1.0/2.0,WP)
    solver%S1_m(-1) = real( 1.0/2.0,WP); solver%S1_m( 0) = real(1.0/2.0,WP)

    ! Optimal weights
    solver%a0 = real(1.0/3.0,WP); solver%a1 = real(2.0/3.0,WP)

  end subroutine init_weno3

  subroutine computeWeno3Coeff(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(solver_t), pointer :: solver => null()
    integer :: i
    real(WP) :: alpha,r,c, a0, a1, a2
    real(WP), parameter :: epsilon = 1.0e-6_WP
    real(WP) :: ISa, ISb, IS0_p, IS1_p, IS0_m, IS1_m, sum_inv
    real(WP), dimension(:), pointer :: F_p=>null(), F_m=>null(), &
         S0_p=>null(), S1_p=>null(), &
         S0_m=>null(), S1_m=>null()
    ! ---------------------------------

     solver => spray%solver    

     F_p => solver%f_p; F_m => solver%f_m
     S0_p => solver%S0_p; S1_p => solver%S1_p
     S0_m => solver%S0_m; S1_m => solver%S1_m
  
     ! compute WENO5 operators...
     do i=spray%kmin,spray%kmax

        ! Direction x - U>0
        ! compute smoothness indicators at regular faces...
        IS0_p = (F_p(i-1)-F_p(  i))**2
        IS1_p = (F_p(i-2)-F_p(i-1))**2

        ! Compute new local weights
        a0 = solver%a0/(epsilon+IS0_p)**2
        a1 = solver%a1/(epsilon+IS1_p)**2
        
        sum_inv = 1.0_WP/(a0+a1+a2)
        solver%weno3p(:,i) = (a0*S0_p(:)+a1*S1_p(:))*sum_inv

        ! Direction x - U<0
        ! compute smoothness indicators at regular faces...
        IS0_m = (F_m(i  )-F_m(i-1))**2
        IS1_m = (F_m(i+1)-F_m(i  ))**2

        ! Compute new local weights
        a0 = solver%a0/(epsilon+IS0_m)**2
        a1 = solver%a1/(epsilon+IS1_m)**2

        sum_inv = 1.0_WP/(a0+a1+a2)
        solver%weno5m(:,i) = (a0*S0_m(:)+a1*S1_m(:))*sum_inv

     end do
   
     ! apply boundary treatment to coefficients
     ! order of scheme is reduced near boundary
     ! to avoid large stencil size
     call apply_weno3_bc()

   contains
     
     subroutine apply_weno3_bc()
       implicit none

       ! Reduced to first order scheme at boundaries
       ! Left boundary
       solver%weno3p(: ,spray%kmino:spray%kmin-1) = 0.0_WP
       solver%weno3p(-1,spray%kmino:spray%kmin-1) = 1.0_WP

       ! Right boundary
       solver%weno3p(: ,spray%kmax+1:spray%kmaxo) = 0.0_WP
       solver%weno3p(-1,spray%kmax+1:spray%kmaxo) = 1.0_WP
       
       ! Left boundary
       solver%weno3m(:,spray%kmino:spray%kmin-1) = 0.0_WP
       solver%weno3m(0,spray%kmino:spray%kmin-1) = 1.0_WP

       ! Right boundary
       solver%weno3m(:,spray%kmax+1:spray%kmaxo) = 0.0_WP
       solver%weno3m(0,spray%kmax+1:spray%kmaxo) = 1.0_WP

     end subroutine apply_weno3_bc
       
  end subroutine computeWeno3Coeff

  subroutine computeWeno3Flux(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(solver_t), pointer :: solver => null()
    integer :: k

    solver => spray%solver

    solver%f_t(spray%kmino:spray%kmin) = solver%f_p(spray%kmino:spray%kmin) &
         + solver%f_m(spray%kmino:spray%kmin)

    solver%f_t(spray%kmax+1:spray%kmaxo) = solver%f_p(spray%kmax+1:spray%kmaxo)&
         + solver%f_m(spray%kmax+1:spray%kmaxo)

    do k=spray%kmin+1,spray%kmax
       solver%f_t(k) = sum(solver%weno3p(:,k)*solver%f_p(k-2:k  )) &
                     + sum(solver%weno3m(:,k)*solver%f_m(k-1:k+1))
    end do

  end subroutine computeWeno3Flux

  subroutine init_weno5(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(solver_t), pointer :: solver=>null()

    solver => spray%solver

    ! Allocate stencils
    allocate(solver%S0_p(-3:+1))
    allocate(solver%S1_p(-3:+1))
    allocate(solver%S2_p(-3:+1))
    allocate(solver%S0_m(-2:+2))
    allocate(solver%S1_m(-2:+2))
    allocate(solver%S2_m(-2:+2))

    ! Allocate WENO5 operators
    allocate(solver%weno5p(-3:1,spray%kmino:spray%kmaxo))
    allocate(solver%weno5m(-2:2,spray%kmino:spray%kmaxo))

    ! Allocate flux vectors
    allocate(solver%f_p(spray%kmino:spray%kmaxo))
    allocate(solver%f_m(spray%kmino:spray%kmaxo))
    allocate(solver%f_t(spray%kmino:spray%kmaxo))

    ! Interpolation at the faces
    solver%S0_p = 0.0_WP; solver%S0_m = 0.0_WP
    solver%S1_p = 0.0_WP; solver%S1_m = 0.0_WP
    solver%S2_p = 0.0_WP; solver%S2_m = 0.0_WP

    solver%S0_p(-1) = real( 1.0/3.0,WP); solver%S0_p( 0) = real(5.0/6.0,WP);  solver%S0_p(+1) = real(-1.0/6.0,WP)
    solver%S0_m(-2) = real(-1.0/6.0,WP); solver%S0_m(-1) = real(5.0/6.0,WP);  solver%S0_m( 0) = real( 1.0/3.0,WP)

    solver%S1_p(-2) = real(-1.0/6.0,WP); solver%S1_p(-1) = real(5.0/6.0,WP);  solver%S1_p( 0) = real( 1.0/3.0,WP)
    solver%S1_m(-1) = real( 1.0/3.0,WP); solver%S1_m( 0) = real(5.0/6.0,WP);  solver%S1_m(+1) = real(-1.0/6.0,WP)

    solver%S2_p(-3) = real( 1.0/3.0,WP); solver%S2_p(-2) = real(-7.0/6.0,WP); solver%S2_p(-1) = real(11.0/6.0,WP)
    solver%S2_m( 0) = real(11.0/6.0,WP); solver%S2_m(+1) = real(-7.0/6.0,WP); solver%S2_m(+2) = real( 1.0/3.0,WP)

    ! Optimal weights
    solver%a0 = 0.3_WP; solver%a1 = 0.6_WP; solver%a2 = 0.1_WP

  end subroutine init_weno5

  subroutine computeWeno5Coeff(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(solver_t), pointer :: solver => null()
    integer :: i
    real(WP) :: alpha,r,c, a0, a1, a2
    real(WP), parameter :: epsilon = 1.0e-6_WP
    real(WP) :: ISa, ISb, IS0_p, IS1_p, IS2_p,IS0_m, IS1_m, IS2_m, sum_inv
    real(WP), dimension(:), pointer :: F_p=>null(), F_m=>null(), &
         S0_p=>null(), S1_p=>null(), S2_p=>null(), &
         S0_m=>null(), S1_m=>null(), S2_m=>null()
    ! ---------------------------------

     solver => spray%solver    

     F_p => solver%f_p; F_m => solver%f_m
     S0_p => solver%S0_p; S1_p => solver%S1_p; S2_p => solver%S2_p
     S0_m => solver%S0_m; S1_m => solver%S1_m; S2_m => solver%S2_m

     ISa = real(13.0/12.0,WP); ISb = real(1.0/4.0,WP)
  
     ! compute WENO5 operators...
     do i=spray%kmin,spray%kmax

        ! Direction x - U>0
        ! compute smoothness indicators at regular faces...
        IS0_p = ISa*(F_p(i+1)-2.0_WP*F_p(i)+F_p(i-1))**2 + ISb*(F_p(i+1)-4.0_WP*F_p(i)+3.0_WP*F_p(i-1))**2
        IS1_p = ISa*(F_p(i-2)-2.0_WP*F_p(i-1)+F_p(i))**2 + ISb*(F_p(i-2)-F_p(i))**2
        IS2_p = ISa*(F_p(i-3)-2.0_WP*F_p(i-2)+F_p(i-1))**2 + ISb*(F_p(i-3)-4.0_WP*F_p(i-2)+3.0_WP*F_p(i-1))**2

        ! Compute new local weights
        a0 = solver%a0/(epsilon+IS0_p)**2
        a1 = solver%a1/(epsilon+IS1_p)**2
        a2 = solver%a2/(epsilon+IS2_p)**2
        
        sum_inv = 1.0_WP/(a0+a1+a2)
        solver%weno5p(:,i) = (a0*S0_p(:)+a1*S1_p(:)+a2*S2_p(:))*sum_inv

        ! Direction x - U<0
        ! compute smoothness indicators at regular faces...
        IS0_m = ISa*(F_m(i-2)-2.0_WP*F_m(i-1)+F_m(i))**2 + ISb*(F_m(i-2)-4.0_WP*F_m(i-1)+3.0_WP*F_m(i))**2
        IS1_m = ISa*(F_m(i+1)-2.0_WP*F_m(i)+F_m(i-1))**2 + ISb*(F_m(i+1)-F_m(i-1))**2
        IS2_m = ISa*(F_m(i+2)-2.0_WP*F_m(i+1)+F_m(i))**2 + ISb*(F_m(i+2)-4.0_WP*F_m(i+1)+3.0_WP*F_m(i))**2

        ! Compute new local weights
        a0 = solver%a0/(epsilon+IS0_m)**2
        a1 = solver%a1/(epsilon+IS1_m)**2
        a2 = solver%a2/(epsilon+IS2_m)**2
        sum_inv = 1.0_WP/(a0+a1+a2)
        solver%weno5m(:,i) = (a0*S0_m(:)+a1*S1_m(:)+a2*S2_m(:))*sum_inv

     end do
   
     ! apply boundary treatment to coefficients
     ! order of scheme is reduced near boundary
     ! to avoid large stencil size
     call apply_weno5_bc()

   contains
     
     subroutine apply_weno5_bc()
       implicit none

       ! Reduced to first order scheme at boundaries
       ! Left boundary
       solver%weno5p(: ,spray%kmino:spray%kmin-1) = 0.0_WP
       solver%weno5p(-1,spray%kmino:spray%kmin-1) = 1.0_WP

       ! Right boundary
       solver%weno5p(: ,spray%kmax+1:spray%kmaxo) = 0.0_WP
       solver%weno5p(-1,spray%kmax+1:spray%kmaxo) = 1.0_WP
       
       ! Left boundary
       solver%weno5m(:,spray%kmino:spray%kmin-1) = 0.0_WP
       solver%weno5m(0,spray%kmino:spray%kmin-1) = 1.0_WP

       ! Right boundary
       solver%weno5m(:,spray%kmax+1:spray%kmaxo) = 0.0_WP
       solver%weno5m(0,spray%kmax+1:spray%kmaxo) = 1.0_WP

     end subroutine apply_weno5_bc
       
  end subroutine computeWeno5Coeff

  subroutine computeWeno5Flux(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(solver_t), pointer :: solver => null()
    integer :: k

    solver => spray%solver

    solver%f_t(spray%kmino:spray%kmin) = solver%f_p(spray%kmino:spray%kmin) &
         + solver%f_m(spray%kmino:spray%kmin)

    solver%f_t(spray%kmax+1:spray%kmaxo) = solver%f_p(spray%kmax+1:spray%kmaxo)&
         + solver%f_m(spray%kmax+1:spray%kmaxo)

    do k=spray%kmin+1,spray%kmax
       solver%f_t(k) = sum(solver%weno5p(:,k)*solver%f_p(k-3:k+1)) &
                     + sum(solver%weno5m(:,k)*solver%f_m(k-2:k+2))
    end do

  end subroutine computeWeno5Flux

  subroutine buildStateVector(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:,:), pointer :: W=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_a=>null(), Y_g=>null(), &
                                       u_l=>null(), u_g=>null(), d2=>null(), d3=>null(), dm=>null(), dvar=>null(), &
                                       Td=>null(), Tg=>null(), b=>null(), k_g=>null(), eps_g=>null(), zvar_g=>null()
    real(WP), dimension(spray%nzo) :: b2

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_a => spray%Y_a; Y_g => spray%Y_g
    u_l => spray%u_l; u_g => spray%u_g; d3 => spray%d3; d2 => spray%d2; dm => spray%dm; dvar => spray%dvar
    Td => spray%Td; Tg => spray%Tg; b => spray%b
    k_g => spray%k_g; eps_g => spray%eps_g; zvar_g=>spray%zvar_g
    W => spray%solver%W

    W = 0.0_WP

    b2 = b**2

    W(1,:) = rho*Y_a*b2
    W(2,:) = rho*Y_v*b2
    W(3,:) = rho*Y_g*u_g*b2
    W(4,:) = rho*Y_g*k_g*b2
    W(5,:) = rho*Y_g*eps_g*b2
    W(6,:) = rho*Y_g*zvar_g*b2

    W(7,:) = rho*Y_l*b2
    W(8,:) = rho*Y_l*u_l*b2
    W(9,:) = rho*Y_l*Td*b2
    W(10,:) = rho*Y_l*dm*b2
    W(11,:) = rho*Y_l*dvar*b2
    W(12,:) = rho*Y_l*d2*b2
    W(13,:) = rho*Y_l*d3*b2
    W(14,:) = rho*Y_g*Tg*b2

  end subroutine buildStateVector

  subroutine buildFluxVector(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:,:), pointer :: F=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_a=>null(), Y_g=>null(), &
                                       u_l=>null(), u_g=>null(), d3=>null(), d2=>null(), dm=>null(), dvar=>null(), &
                                       Td=>null(), Tg=>null(), b=>null(), k_g=>null(), eps_g=>null(), zvar_g=>null()
    real(WP), dimension(spray%nzo) :: b2

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_a => spray%Y_a; Y_g => spray%Y_g
    u_l => spray%u_l; u_g => spray%u_g; d3 => spray%d3; d2 => spray%d2; dm => spray%dm; dvar => spray%dvar;
    Td => spray%Td; Tg => spray%Tg; b => spray%b
    k_g => spray%k_g; eps_g => spray%eps_g; zvar_g => spray%zvar_g

    F => spray%solver%F

    Y_g = 1.0_WP - Y_l

    F = 0.0_WP

    b2 = b**2

    F(1,:) = rho*Y_a*u_g*b2
    F(2,:) = rho*Y_v*u_g*b2
    F(3,:) = rho*Y_g*u_g*u_g*b2
    F(4,:) = rho*Y_g*k_g*u_g*b2
    F(5,:) = rho*Y_g*eps_g*u_g*b2
    F(6,:) = rho*Y_g*zvar_g*u_g*b2

    F(7,:) = rho*Y_l*u_l*b2
    F(8,:) = rho*Y_l*u_l*u_l*b2
    F(9,:) = rho*Y_l*Td*u_l*b2
    F(10,:) = rho*Y_l*dm*u_l*b2
    F(11,:) = rho*Y_l*dvar*u_l*b2
    F(12,:) = rho*Y_l*d2*u_l*b2
    F(13,:) = rho*Y_l*d3*u_l*b2
    F(14,:) = rho*Y_g*Tg*u_g*b2

  end subroutine buildFluxVector

  subroutine buildSourceVector(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:,:), pointer :: S=>null()
    real(WP), dimension(:), pointer :: u_l=>null(), d2=>null(), d3=>null(), dm=>null(), Td=>null(), &
                                       Tg=>null(), Tv=>null(), b=>null(), &
                                       omega_ent=>null(), omega_vap=>null(), omega_vapdm=>null(), &
                                       omega_vapd2=>null(), omega_vapd3=>null(), f_drag=>null(), &
                                       omega_bre1=>null(), omega_bre2=>null(), omega_bre3=>null(), &
                                       omega_exp=>null(), omega_expdm=>null(), omega_expd2=>null(), omega_expd3=>null(), &
                                       omega_T=>null(), omega_k_g_p=>null(), omega_k_g_d=>null(), &
                                       omega_eps_g_p=>null(), omega_eps_g_d=>null(), omega_zvar_g_p, &
                                       omega_zvar_g_d
    real(WP), dimension(spray%nzo) :: b2
    real(WP) :: Ta, Tfuel

    u_l => spray%u_l; d2 => spray%d2; dm => spray%dm; Td => spray%Td; Tv => spray%Tv; Tg => spray%Tg; b => spray%b

    omega_ent => spray%omega_ent; omega_vap => spray%omega_vap; omega_vapdm => spray%omega_vapdm; 
    omega_vapd2 => spray%omega_vapd2; omega_vapd3 => spray%omega_vapd3; f_drag => spray%f_drag; 
    omega_exp => spray%omega_exp; omega_expdm => spray%omega_expdm; 
    omega_expd2 => spray%omega_expd2; omega_expd3 => spray%omega_expd3; 
    omega_bre1 => spray%omega_bre1; omega_bre2 => spray%omega_bre2; omega_bre3 => spray%omega_bre3;
    omega_T => spray%omega_T; omega_k_g_p => spray%omega_k_g_p; omega_k_g_d => spray%omega_k_g_d
    omega_eps_g_p => spray%omega_eps_g_p; omega_eps_g_d => spray%omega_eps_g_d
    omega_zvar_g_p => spray%omega_zvar_g_p; omega_zvar_g_d => spray%omega_zvar_g_d

    !omega_exp = 0.0_WP; omega_expdm = 0.0_WP; omega_expd2 = 0.0_WP; omega_expd3 = 0.0_WP

    Ta = spray%T_a; Tfuel = spray%T_fuel

    S => spray%solver%S

    S = 0.0_WP

    b2 = b**2

    S(1,:) = omega_ent
    S(2,:) = (omega_vap+omega_exp)*b2
    S(3,:) = (-f_drag + (omega_vap+omega_exp)*u_l)*b2
    S(4,:) = omega_k_g_p - omega_k_g_d
    S(5,:) = omega_eps_g_p - omega_eps_g_d
    S(6,:) = omega_zvar_g_p - omega_zvar_g_d

    S(7,:) = -(omega_vap+omega_exp)*b2
    S(8,:) = (f_drag - (omega_vap+omega_exp)*u_l)*b2
    S(9,:) = (omega_T - (omega_vap+omega_exp)*Td)*b2
    S(10,:) = (-omega_bre1 - 4.0_WP/3.0_WP*omega_vapdm - 2.0_WP/3.0_WP*omega_expdm)*b2
    S(11,:) = (-omega_bre2 - 5.0_WP/3.0_WP*omega_vapd2)*b2 &
               - 2.0_WP*dm*S(10,:) + S(7,:)*dm**2
    S(12,:) = (-omega_bre2 - 5.0_WP/3.0_WP*omega_vapd2 - 1.0_WP/3.0_WP*omega_expd2)*b2
    S(13,:) = (-omega_bre3 - 2.0_WP       *omega_vapd3)*b2

    S(14,:) = -(spray%C_l/spray%Cp_g)*(omega_T)*b2+(spray%Cp_v/spray%Cp_g)*(omega_vap+omega_exp)*Tv*b2 &
                        +omega_ent*(spray%Cp_a/spray%Cp_g)*(Ta/Tfuel)

  end subroutine buildSourceVector

  subroutine buildScalars(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:,:), pointer :: W=>null(), F=>null(), S=>null(), Sc=>null(), Wt=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_a=>null(), Y_g=>null(), &
                                       u_l=>null(), u_g=>null(), d2=>null(), d3=>null(), dm=>null(), dvar=>null(), &
                                       Td=>null(), b=>null(), k_g=>null(), eps_g=>null(), zvar_g=>null(), omega_vap=>null()
    real(WP), dimension(spray%nzo) :: b2
    real(WP) :: curr_vap_mass, factor
    integer :: i

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_a => spray%Y_a; Y_g => spray%Y_g
    u_l => spray%u_l; u_g => spray%u_g; d3 => spray%d3; d2 => spray%d2; dm => spray%dm; dvar => spray%dvar; 
    Td => spray%Td; b => spray%b
    k_g => spray%k_g; eps_g => spray%eps_g; zvar_g=>spray%zvar_g
    W => spray%solver%W
    F => spray%solver%F
    S => spray%solver%S

    Sc => spray%Sc
    Wt => spray%Wt
    omega_vap => spray%omega_vap

    W(spray%neq+1:,:) = 0.0_WP
    F(spray%neq+1:,:) = 0.0_WP
    S(spray%neq+1:,:) = 0.0_WP

    b2 = b**2

    curr_vap_mass = sum(rho*spray%rho_l*Y_v*spray%dz*Pi*b2*spray%D_eff**3)
    factor = int(curr_vap_mass/(spray%inj_mass/spray%nsc))
    do i = 1,spray%nsc
       W(spray%neq+i,:) = rho*Sc(i,:)*b2
       F(spray%neq+i,:) = rho*Sc(i,:)*u_g*b2
       if(factor==i-1) then
          S(spray%neq+i,:) = omega_vap*b2
       end if
    end do

  end subroutine buildScalars

  subroutine computeResidual_UPWIND1(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: i, k
    integer, pointer :: kmin, kmax, kmino, kmaxo
    real(WP), dimension(:,:), pointer :: divc=>null(), W=>null(), Wold=>null(), F=>null(), S=>null(), Res=>null()
    real(WP), dimension(:), pointer :: u_l=>null(), u_g=>null(), alpha_l=>null(), alpha_g=>null(), &
         Flux=>null()
    real(WP) :: a_l, a_r, FL, FR
    type(solver_t), pointer :: solver=>null()

    solver => spray%solver

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    divc => spray%solver%divc; W => spray%solver%W; Wold => spray%solver%Wold; 
    F => spray%solver%F; S => spray%solver%S

    u_l => spray%u_l; u_g => spray%u_g; alpha_l => spray%solver%alpha_l; alpha_g => spray%solver%alpha_g
    Flux => spray%solver%Flux; Res => spray%solver%Res

    do i = 1,6

       if(i >= abs(spray%skip_turb)) cycle

       do k = kmin,kmax
          
          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*F(i,k-1:k))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
          end if

       end do

    end do

    do i = 7,14

       if(i == spray%skip_d2) cycle
       
       if(i == spray%skip_d3) cycle

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*F(i,k-1:k))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
          end if

       end do

    end do

    if(spray%nsc > 0) then
       do i = spray%neq+1,spray%neq+spray%nsc

          do k = kmin,kmax

             if(spray%solver%strang) then
                Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*F(i,k-1:k))) ! - S(i,k))
             else
                Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
             end if

          end do

       end do
    end if

    Wold = W

  end subroutine computeResidual_UPWIND1

  subroutine computeResidual_LF(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: i, k
    integer, pointer :: kmin, kmax, kmino, kmaxo
    real(WP), dimension(:,:), pointer :: divc=>null(), W=>null(), Wold=>null(), F=>null(), S=>null(), Res=>null()
    real(WP), dimension(:), pointer :: u_l=>null(), u_g=>null(), alpha_l=>null(), alpha_g=>null(), &
                                       Flux=>null(), alpha_m=>null()
    real(WP) :: a_l, a_r, FL, FR

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    divc => spray%solver%divc; W => spray%solver%W; Wold => spray%solver%Wold; 
    F => spray%solver%F; S => spray%solver%S

    u_l => spray%u_l; u_g => spray%u_g; alpha_l => spray%solver%alpha_l; alpha_g => spray%solver%alpha_g;
    alpha_m => spray%solver%alpha_m

    Flux => spray%solver%Flux; Res => spray%solver%Res

    do k=kmino+1,kmaxo
       alpha_l(k) = max(abs(u_l(k-1)),abs(u_l(k))); !hx/dt;%max(abs(u_l(i-1)),abs(u_l(i)));
       alpha_g(k) = max(abs(u_g(k-1)),abs(u_g(k))); !hx/dt;%max(abs(u_l(i)),abs(u_l(i+1)));
       alpha_m(k) = max(alpha_l(k),alpha_g(k))
    end do

    do i = 1,6

       if(i >= abs(spray%skip_turb)) cycle

       Flux = 0.5_WP*(F(i,1:kmaxo-1)+F(i,2:kmaxo)) + 0.5_WP*alpha_m*(W(i,1:kmaxo-1)-W(i,2:kmaxo))

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*Flux(k-1:k))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*Flux(k-1:k)) - S(i,k))
          end if

       end do

    end do

    do i = 7,14

       if(i == spray%skip_d2) cycle

       if(i == spray%skip_d3) cycle

       Flux = 0.5_WP*(F(i,1:kmaxo-1)+F(i,2:kmaxo)) + 0.5_WP*alpha_m*(W(i,1:kmaxo-1)-W(i,2:kmaxo))

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*Flux(k-1:k))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*Flux(k-1:k)) - S(i,k))
          end if

       end do

    end do

    if(spray%nsc > 0) then
       do i = spray%neq+1,spray%neq+spray%nsc

          Flux = 0.5_WP*(F(i,1:kmaxo-1)+F(i,2:kmaxo)) + 0.5_WP*alpha_m*(W(i,1:kmaxo-1)-W(i,2:kmaxo))

          do k = kmin,kmax
             if(spray%solver%strang) then
                Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*Flux(k-1:k))) ! - S(i,k))
             else
                Res(i,k) = -spray%dt*(sum(divc(:,k)*Flux(k-1:k)) - S(i,k))
             end if

          end do

       end do
    end if

    Wold = W

    spray%solver%rktvd%RK(spray%solver%rktvd%niter,:,:) = W
    spray%solver%rktvd%Res(spray%solver%rktvd%niter,:,:) = Res

  end subroutine computeResidual_LF

  subroutine computeResidual_WENO3(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: i, k
    integer, pointer :: kmin, kmax, kmino, kmaxo
    real(WP), dimension(:,:), pointer :: divc=>null(), W=>null(), Wold=>null(), F=>null(), S=>null(), Res=>null()
    real(WP), dimension(:), pointer :: u_l=>null(), u_g=>null(), alpha_l=>null(), alpha_g=>null(), alpha_m=>null(),&
         Flux=>null()
    real(WP) :: a_l, a_r, FL, FR
    real(WP), dimension(spray%nzo) :: c_s
    type(solver_t), pointer :: solver=>null()

    solver => spray%solver

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    divc => spray%solver%divc; W => spray%solver%W; Wold => spray%solver%Wold; 
    F => spray%solver%F; S => spray%solver%S

    u_l => spray%u_l; u_g => spray%u_g; alpha_l => spray%solver%alpha_l; alpha_g => spray%solver%alpha_g;
    alpha_m => spray%solver%alpha_m
    Flux => spray%solver%Flux; Res => spray%solver%Res

    ! Global LF splitting
    alpha_l = maxval(u_l); alpha_g = maxval(u_g)
    alpha_m = max(alpha_l,alpha_g)

    do i = 1,6

       if(i >= abs(spray%skip_turb)) cycle

       solver%f_p = 0.5_WP*(F(i,:) + alpha_g*W(i,:))
       solver%f_m = 0.5_WP*(F(i,:) - alpha_g*W(i,:))

       call computeWeno3Coeff(spray)

       call computeWeno3Flux(spray)

       !where(solver%f_t < 0.0_WP) solver%f_t = 0.0_WP

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*solver%f_t(k:k+1))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
          end if

       end do

    end do

    do i = 7,14

       if(i == spray%skip_d2) cycle

       if(i == spray%skip_d3) cycle

       solver%f_p = 0.5_WP*(F(i,:) + alpha_l*W(i,:))
       solver%f_m = 0.5_WP*(F(i,:) - alpha_l*W(i,:))

       call computeWeno3Coeff(spray)

       call computeWeno3Flux(spray)

       !where(solver%f_t < 0.0_WP) solver%f_t = 0.0_WP

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*solver%f_t(k:k+1))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
          end if

       end do

    end do

    if(spray%nsc > 0) then
       do i = spray%neq+1,spray%neq+spray%nsc

          solver%f_p = 0.5_WP*(F(i,:) + alpha_g*W(i,:))
          solver%f_m = 0.5_WP*(F(i,:) - alpha_g*W(i,:))

          call computeWeno3Coeff(spray)

          call computeWeno3Flux(spray)

          !where(solver%f_t < 0.0_WP) solver%f_t = 0.0_WP

          do k = kmin,kmax

             if(spray%solver%strang) then
                Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*solver%f_t(k:k+1))) ! - S(i,k))
             else
                Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
             end if

          end do

       end do
    end if

    Wold = W

    spray%solver%rktvd%RK(spray%solver%rktvd%niter,:,:) = W
    spray%solver%rktvd%Res(spray%solver%rktvd%niter,:,:) = Res

  end subroutine computeResidual_WENO3

  subroutine computeResidual_WENO5(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: i, k
    integer, pointer :: kmin, kmax, kmino, kmaxo
    real(WP), dimension(:,:), pointer :: divc=>null(), W=>null(), Wold=>null(), F=>null(), S=>null(), Res=>null()
    real(WP), dimension(:), pointer :: u_l=>null(), u_g=>null(), alpha_l=>null(), alpha_g=>null(), &
         Flux=>null()
    real(WP) :: a_l, a_r, FL, FR
    type(solver_t), pointer :: solver=>null()

    solver => spray%solver

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    divc => spray%solver%divc; W => spray%solver%W; Wold => spray%solver%Wold; 
    F => spray%solver%F; S => spray%solver%S

    u_l => spray%u_l; u_g => spray%u_g; alpha_l => spray%solver%alpha_l; alpha_g => spray%solver%alpha_g
    Flux => spray%solver%Flux; Res => spray%solver%Res

    ! Global LF splitting
    alpha_l = maxval(u_l); alpha_g = maxval(u_g)

    do i = 1,6

       if(i >= abs(spray%skip_turb)) cycle

       solver%f_p = 0.5_WP*(F(i,:) + alpha_g*W(i,:))
       solver%f_m = 0.5_WP*(F(i,:) - alpha_g*W(i,:))

       call computeWeno5Coeff(spray)

       call computeWeno5Flux(spray)

       !where(solver%f_t < 0.0_WP) solver%f_t = 0.0_WP

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*solver%f_t(k:k+1))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
          end if

       end do

    end do

    do i = 7,14

       if(i == spray%skip_d2) cycle

       if(i == spray%skip_d3) cycle

       solver%f_p = 0.5_WP*(F(i,:) + alpha_l*W(i,:))
       solver%f_m = 0.5_WP*(F(i,:) - alpha_l*W(i,:))

       call computeWeno5Coeff(spray)

       call computeWeno5Flux(spray)

       !where(solver%f_t < 0.0_WP) solver%f_t = 0.0_WP

       do k = kmin,kmax

          if(spray%solver%strang) then
             Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*solver%f_t(k:k+1))) ! - S(i,k))
          else
             Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
          end if

       end do

    end do

    if(spray%nsc > 0) then
       do i = spray%neq+1,spray%neq+spray%nsc

          solver%f_p = 0.5_WP*(F(i,:) + alpha_g*W(i,:))
          solver%f_m = 0.5_WP*(F(i,:) - alpha_g*W(i,:))

          call computeWeno5Coeff(spray)

          call computeWeno5Flux(spray)

          !where(solver%f_t < 0.0_WP) solver%f_t = 0.0_WP

          do k = kmin,kmax

             if(spray%solver%strang) then
                Res(i,k) = -0.5_WP*spray%dt*(sum(divc(:,k)*solver%f_t(k:k+1))) ! - S(i,k))
             else
                Res(i,k) = -spray%dt*(sum(divc(:,k)*F(i,k-1:k)) - S(i,k))
             end if

          end do

       end do
    end if

    Wold = W

    spray%solver%rktvd%RK(spray%solver%rktvd%niter,:,:) = W
    spray%solver%rktvd%Res(spray%solver%rktvd%niter,:,:) = Res

  end subroutine computeResidual_WENO5

  subroutine updateSolution(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(rktvd_solver_t), pointer :: rktvd
    real(WP), dimension(:,:,:), pointer :: RK=>null(), Res=>null()
    real(WP), dimension(:,:), pointer :: W=>null()
    real(WP), dimension(:,:), pointer :: coeff
    integer :: k, istage
    integer, pointer :: kmino, kmaxo, kmin, kmax

    kmin => spray%kmin; kmax => spray%kmax; kmino => spray%kmino; kmaxo => spray%kmaxo

    W => spray%solver%W
    rktvd => spray%solver%rktvd

    RK => rktvd%RK; Res => rktvd%Res; coeff => rktvd%coeff

    istage = rktvd%niter

    select case (rktvd%order)

    case(1)

       do k = kmin,kmax

          W(:,k) = coeff(istage,1)*RK(istage,:,k) + coeff(istage,2)*Res(istage,:,k)

       end do

    case(2)

       if(istage == 1) then
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage,:,k) + coeff(istage,2)*Res(istage,:,k)

          end do
       else
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage-1,:,k) + coeff(istage,2)*RK(istage,:,k) &
                  + coeff(istage,3)*Res(istage,:,k)

          end do
       end if

    case(3)

       if(istage == 1) then
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage,:,k) + coeff(istage,2)*Res(istage,:,k)

          end do
       else if(istage == 2) then
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage-1,:,k) + coeff(istage,2)*RK(istage,:,k) &
                  + coeff(istage,3)*Res(istage,:,k)

          end do
       else
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage-2,:,k) + coeff(istage,2)*RK(istage,:,k) &
                  + coeff(istage,3)*Res(istage,:,k)

          end do
       end if

    case(4)
       write(*,*) 'Warning! Scheme not working properly...'
       if(istage == 1) then
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage,:,k) + coeff(istage,2)*Res(istage,:,k)

          end do
       else if(istage == 2) then
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage-1,:,k) - coeff(istage,2)*Res(istage-1,:,k) &
                  + coeff(istage,3)*RK(istage,:,k)   + coeff(istage,4)*Res(istage,:,k)

          end do
       else if(istage == 3) then
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage-2,:,k) - coeff(istage,2)*Res(istage-2,:,k) &
                  + coeff(istage,3)*RK(istage-1,:,k) - coeff(istage,4)*Res(istage-1,:,k) &
                  + coeff(istage,5)*RK(istage,:,k)   + coeff(istage,6)*Res(istage,:,k)

          end do
       else 
          do k = kmin,kmax

             W(:,k) = coeff(istage,1)*RK(istage-3,:,k) + coeff(istage,2)*Res(istage-3,:,k) &
                  + coeff(istage,3)*RK(istage-2,:,k) + coeff(istage,4)*Res(istage-2,:,k) &
                  + coeff(istage,5)*RK(istage-1,:,k)                                     &
                  + coeff(istage,6)*RK(istage,:,k)   + coeff(istage,7)*Res(istage,:,k)

          end do
       end if

    case default
       write(*,*) 'Scheme not implemented. Order of the scheme: ',rktvd%order
    end select

  end subroutine updateSolution

  subroutine addSource(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: i, k
    integer, pointer :: kmin, kmax, kmino, kmaxo
    real(WP), dimension(:,:), pointer :: divc=>null(), W=>null(), Wold=>null(), F=>null(), S=>null(), Res=>null()
    real(WP), dimension(:), pointer :: u_l=>null(), u_g=>null(), alpha_l=>null(), alpha_g=>null(), &
                                       Flux=>null(), alpha_m=>null()

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    W => spray%solver%W; Wold => spray%solver%Wold; 
    S => spray%solver%S
    Res => spray%solver%Res
    do i = 1,6

       if(i >= abs(spray%skip_turb)) cycle

       do k = kmin,kmax

          Res(i,k) = spray%dt*S(i,k)

       end do

    end do

    do i = 7,14

       if(i == spray%skip_d2) cycle

       if(i == spray%skip_d3) cycle

       do k = kmin,kmax

          Res(i,k) = spray%dt*S(i,k)

       end do

    end do

    if(spray%nsc > 0) then
       do i = spray%neq+1,spray%neq+spray%nsc

          do k = kmin,kmax

             Res(i,k) = spray%dt*S(i,k)

          end do

       end do
    end if

    Wold = W

    spray%solver%rktvd%RK(spray%solver%rktvd%niter,:,:) = W
    spray%solver%rktvd%Res(spray%solver%rktvd%niter,:,:) = Res

  end subroutine addSource

  subroutine updateFlowVariables(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l
    real(WP), dimension(:,:), pointer :: W=>null(), Sc=>null(), Tdi=>null(), dsd=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_a=>null(), Y_g=>null(), &
                                       u_l=>null(), u_g=>null(), d3=>null(), d2=>null(), dm=>null(), dvar=>null(), &
                                       Td=>null(), Tg=>null(), Tv=>null(), b=>null(), k_g=>null(), eps_g=>null(), zvar_g=>null(), &
                                       mu_t_g=>null(), chi_g=>null(), DRg=>null()

    real(WP), pointer :: DRv=>null(), DRa=>null()
    real(WP), dimension(spray%nzo) :: T_low, T_high, Tg_low, Tg_high
    real(WP) :: nu, eps = 1.0E-16_WP, C_scale
    integer :: k, i

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_a => spray%Y_a; Y_g => spray%Y_g
    u_l => spray%u_l; u_g => spray%u_g; d3 => spray%d3; d2 => spray%d2; dm => spray%dm; dvar => spray%dvar; 
    Td => spray%Td; Tv => spray%Tv; Tg => spray%Tg; Tdi => spray%Tdi; dsd => spray%dsd; b => spray%b
    k_g => spray%k_g; eps_g => spray%eps_g; zvar_g => spray%zvar_g; mu_t_g => spray%mu_t_g; chi_g => spray%chi_g
    DRv => spray%DRv; DRa => spray%DRa; DRg => spray%DRg
    W => spray%solver%W
    Sc => spray%Sc

    T_low = spray%MP/spray%T_fuel
    Tg_low = 298.15_WP/spray%T_fuel
    Tg_high = 1.0_WP  !900.0_WP/spray%T_fuel

    !rho = 0.0_WP; u_l = 0.0_WP; 
 
    !Y_l = 0.0_WP; Y_a = 0.0_WP; Y_v = 0.0_WP;
    !dm = 0.0_WP; d2 = 0.0_WP; Td = 0.0_WP; dvar = 0.0_WP

    dm = 0.0_WP; dvar = 0.0_WP; d2 = 0.0_WP

    do k = spray%kmino,spray%kmaxo

       if ( (W(1,k)+W(2,k)) .gt. 0.0_WP ) then
          u_g(k) = min(1.0_WP,max(0.0_WP,W(3,k)/(W(1,k)+W(2,k))))
          Tg(k) = max(Tg_low(k),min(Tg_high(k),W(14,k)/(W(1,k)+W(2,k))))

          ! Check if the gas temperature is out of the range (250.0_WP,600.0_WP) PS: The range is selected by intuition
          ! One should change this accoridng to the operating conditions
          !if(W(14,k)/(W(1,k)+W(2,k)) .gt. Tg_high(k) .or. W(14,k)/(W(1,k)+W(2,k)) .lt. Tg_low(k))then
          !  write(*,*)'Gas temperature is out of the selected range:', spray%T_fuel*(W(14,k)/(W(1,k)+W(2,k))), spray%Td(k)*spray%T_fuel
          !  stop
          !end if
       end if

       Y_l(k) = W(7,k)/(W(1,k)+W(2,k)+W(7,k))
       if ( Y_l(k) <= 1.0E-06_WP ) Y_l(k) = 0.0_WP
       Y_v(k) = max(0.0_WP,W(2,k)/(W(1,k)+W(2,k)+W(7,k)))
       Y_a(k) = max(0.0_WP,1.0_WP - Y_l(k) - Y_v(k))
       Y_g(k) = 1.0_WP - Y_l(k)
       rho(k) = 1.0_WP/(Y_l(k) + DRv*Y_v(k) + DRa*Y_a(k))
       spray%b_old(k) = b(k)
       b(k) = max(0.5_WP,sqrt((W(1,k)+W(2,k)+W(7,k))/rho(k)))

       if(Y_l(k) > 1.0E-06_WP) then

          u_l(k) = max(0.0_WP,W(8,k)/W(7,k))

          if(u_l(k) < 0.0_WP .or. u_l(k) > 1.0_WP) then
             write(*,*) 'ul, ug:', u_l(k), u_g(k)
             stop
          end if
          if(u_g(k) < 0.0_WP .or. u_g(k) > 1.0_WP) then
             write(*,*) 'ul, ug:', u_l(k), u_g(k)
             stop
          end if

          dm(k) = max(0.0_WP,W(10,k)/W(7,k))
          !dm(k) = W(10,k)/rho(k)/Y_l(k)/b(k)**2

          dvar(k) = max(0.0_WP,W(11,k)/W(7,k))
          !dvar(k) = W(11,k)/rho(k)/Y_l(k)/b(k)**2

          !d2(k) = min(1.0_WP,max(0.0_WP,W(12,k)/W(7,k)))

          !d3(k) = min(1.0_WP,max(0.0_WP,W(13,k)/rho(k)/Y_l(k)/b(k)**2))

          d2(k) = dm(k)**2 + dvar(k)
          !dvar(k) = d2(k) - dm(k)**2

          if(spray%Td(k)*spray%T_fuel.ge.spray%T_boil)then
            Tv(k) = (spray%C_l(k)*Td(k)*spray%T_fuel/spray%Cp_v-spray%L_fb/spray%Cp_v)/spray%T_fuel
          else
            Tv(k) = (spray%C_l(k)*Td(k)*spray%T_fuel/spray%Cp_v-spray%L_f(k)/spray%Cp_v)/spray%T_fuel
          end if

          if (spray%flash_boiling) then
            T_high(k) =1.0_WP
          else
            T_high(k) = spray%T_sat(k)/spray%T_fuel
          end if
          Td(k) = max(T_low(k),min(T_high(k),W(9,k)/W(7,k)))         

          ! Correct temperature of diameter classes with mean Temperature
          C_scale = 1.0_WP
          if(Td(k)>0.0_WP) then
             C_scale = sum(Tdi(:,k)*dsd(:,k))/Td(k)
          end if
          if(C_scale > 0.0_WP) then
             Tdi(:,k) = Tdi(:,k)/C_scale
          else
             Tdi(:,k) = Td(k)
          end if
          where(Tdi(:,k) < T_low(k)) Tdi(:,k) = T_low(k)
          where(Tdi(:,k) > T_high(k)) Tdi(:,k) = T_high(k)
          where(Tdi(:,k) /= Tdi(:,k)) Tdi(:,k) = Td(k)
       end if

       if(Y_l(k) <= 1e-06) then
          Y_l(k) = 0.0_WP
          u_l(k) = 0.0_WP
          dm(k) = 0.0_WP
          d2(k) = 0.0_WP
          dvar(k) = 0.0_WP
          d3(k) = 0.0_WP
          Td(k) = 0.0_WP
       end if

       !if(abs(u_g(k)) > 0.0_WP .and. Y_g(k) > 0.0_WP) then
       ! TODO: Is this condition necessary???
       !if(Y_g(k) > 0.01_WP) then
       if(Y_g(k) > eps) then
          if(W(4,k) > 0.0_WP) then
             k_g(k) = max(0.0_WP,W(4,k)/rho(k)/Y_g(k)/b(k)**2)
          end if
          if(W(5,k) > 0.0_WP) then
             eps_g(k) = max(0.0_WP,W(5,k)/rho(k)/Y_g(k)/b(k)**2)
          end if
          zvar_g(k) = max(0.0_WP,W(6,k)/rho(k)/Y_g(k)/b(k)**2)
          if(k_g(k) /= k_g(k)) then
             write(*,*) 'k_g is NaN'
          end if
          if(k_g(k) > 1.0e+100_WP) then
             write(*,*) 'k_g is NaN'
          end if          
          if(eps_g(k) /= eps_g(k)) then
             write(*,*) 'eps_g is NaN'
          end if
          if(zvar_g(k) /= zvar_g(k)) then
             write(*,*) 'zvar_g is NaN'
          end if

          if(eps_g(k) > 0.0_WP) then
             !mu_t_g(k) = spray%c_mu*rho(k)*sqrt(Y_g(k)*Y_a(k))*k_g(k)**2/eps_g(k)
             if(DRg(k) /= DRg(k)) then
                mu_t_g(k) = 0.0_WP
             else
                !mu_t_g(k) = spray%c_mu*sqrt(1.0_WP/DRa/DRg(k))*k_g(k)**2/eps_g(k)
                !mu_t_g(k) = spray%c_mu/DRg(k)*k_g(k)**2/eps_g(k)
                mu_t_g(k) = spray%c_mu*sqrt(spray%rho_ra(k)*spray%rho_l/DRg(k))*k_g(k)**2/eps_g(k)*(0.5_WP*spray%U_inj**2*spray%tau)
             end if
             if(mu_t_g(k) /= mu_t_g(k)) then
                write(*,*) 'mu_t_g is NaN'
             end if
          end if
          if(k_g(k) > 0.0_WP) then
             chi_g(k) = spray%c_zvar*eps_g(k)/k_g(k)*zvar_g(k)
          else
             chi_g(k) = min(1000.0_WP*spray%tau,2.0_WP*b(k)/max(eps,u_g(k)))
          end if
       end if

       if(spray%nsc > 0) then
          do i=1,spray%nsc-1
             Sc(i,k) = max(0.0_WP,W(i+spray%neq,k)/(W(1,k)+W(2,k)+W(7,k)))
          end do
          Sc(spray%nsc,k) = spray%Y_v(k) - sum(Sc(1:spray%nsc-1,k))
       end if

    end do

    !call compute_b(spray)
    ! Constant local gas temperature
    !spray%Tg(spray%kmino:spray%kmin-1) = spray%T_a/spray%T_fuel
    ! Remark:
    ! Choice of local gas temperature calculation 
    ! does not affect liquid penetration significantly
    ! ~ 1 mm difference for OME1
    ! Local gas temperature adjusted by vapor temperature
    ! Warning: This approach results in jumps in gas and reference temperature, 
    ! which cause jumps in gaseous properties.
    !where(spray%Td > 0.0_WP) spray%Tg = (spray%Y_v*(spray%Td - spray%De*spray%CR*spray%LR) &
    !                                  +  spray%Y_a*spray%T_a/spray%T_fuel)/spray%Y_g

    !!$$where(spray%Td > 0.0_WP) spray%Tg = (spray%Y_v*spray%Tv/spray%T_fuel &
    !!$$                                  +  spray%Y_a*spray%T_a/spray%T_fuel)/spray%Y_g
    !do k = spray%kmino,spray%kmaxo
    !  if (spray%Tg(k)*spray%T_fuel.lt.260.0_WP) then
    !     write(*,*)'Gas is colling down, clipping is required: Td, Tg, Ta:', Td(k)*spray%T_fuel, spray%Tg(k)*spray%T_fuel, spray%T_a
    !  end if
    !end do
    !where(spray%Td > 0.0_WP) spray%Tg = &
    !     (spray%Y_v*spray%Td/spray%T_fuel + spray%Y_a*spray%T_a/spray%T_fuel)/spray%Y_g

    !spray%mu_t_g = spray%c_mu*spray%k_g**2/spray%eps_g
    !spray%c_mu*sqrt(1.0_WP/spray%DRa/spray%DRg)*spray%k_g**2/spray%eps_g

    !spray%Tg = spray%T_a/spray%T_fuel
    !First cell is a BC
    spray%mu_t_g(spray%kmin) = 0.0_WP
    ! n-dodecane
    !nu = 3.475402137_WP
    ! n-heptane
    !nu = 3.51283381901_WP 

    spray%Zmix_g = (spray%stoic_coeff*spray%Y_v/spray%Y_g-spray%Y_O2*spray%Y_a/spray%Y_g+spray%Y_O2)/(spray%stoic_coeff+spray%Y_O2)

    if(spray%nsc > 0) then
       do i=1,spray%nsc
          spray%Wt(i,:) = spray%Sc(i,:)/spray%Y_v !(spray%stoic_coeff*spray%Sc(i,:)/spray%Y_g-spray%Y_O2*spray%Y_a/spray%Y_g+spray%Y_O2)/(spray%stoic_coeff+spray%Y_O2)/spray%Zmix_g
          !spray%Wt(i,:) = (spray%stoic_coeff*spray%Sc(i,:)/spray%Y_g-spray%Y_O2*spray%Y_a/spray%Y_g/spray%nsc+spray%Y_O2/spray%nsc)/(spray%stoic_coeff+spray%Y_O2)/spray%Zmix_g
       end do
       where(spray%Wt /= spray%Wt .or. spray%Wt < 0.0_WP) spray%Wt = 0.0_WP
    end if

    !spray%zmix_g = 12*12.0107_WP*spray%Y_v/spray%MW_f + 26*1.00794_WP*spray%Y_v/spray%MW_f & ! Fuel
    !             + 12.0107_WP*spray%Y_a*0.093629/44.01_WP & ! CO2
    !             + 2*1.016*spray%Y_a*0.022292/18.01528_WP ! H2O

  end subroutine updateFlowVariables

end module solver_func
