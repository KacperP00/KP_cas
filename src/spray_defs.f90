module spray_defs
  use precision
  use pc_defs
  use solver_defs

  implicit none

  ! Definition of spray object
  type spray_t

     ! Number of grid points
     integer, pointer :: Nz, nzo, step
     integer, pointer :: kmin, kmax, kmino, kmaxo

     ! Numerics
     real(WP), pointer :: dz, Lz
     real(WP), pointer :: dt, ndtime, final_time, ndftime, tau
     real(WP), pointer :: MaxCFL, CFL

     ! Non-dimensional numbers
     real(WP), pointer :: Re, We, WR, De, DRa, DRv, VRa, VRv
     real(WP), dimension(:), pointer :: DRg, VRg, LR, CR
     
     ! Geometric variables
     real(WP), pointer :: noz_D, noz_LD, noz_rD, noz_Dsac
     
     ! Flow variables
     real(WP), dimension(:), pointer :: z, rho, Y_l, Y_v, Y_a, Y_g, u_l, u_g, d2, dm, d3, Td, b, Tg

     ! Source terms
     real(WP), dimension(:), pointer :: omega_ent, omega_vap, omega_vapdm, omega_vapd2, omega_vapd3, &
                                        f_drag, omega_bre1, omega_bre2, omega_bre3, omega_T

     ! Fuel name
     character(len=128), pointer :: Fuel

     ! Pure chemical object for liquid phase
     type(pc_t), dimension(:), pointer :: pc_l
     
     ! Pure chemical object for vapor phase
     type(pc_t), dimension(:), pointer :: pc_v
 
     ! Fuel properties
     real(WP), pointer :: T_fuel, sigma, rho_l, visc_l, MW_f, MP, NBP
     real(WP), dimension(:), pointer :: L_f, C_l, p_vap

     ! Fuel properties from table
     character(len=128), pointer :: LFPTname, VFPTname
     real(WP), dimension(:,:), pointer :: LFPT, VFPT
     
     ! Fuel vapor properties
     real(WP), pointer :: rho_v, visc_v, lambda_v, Cp_v

     ! Ambient properties
     real(WP), pointer :: MW_a, Z_a, P_a, T_a, rho_a, visc_a, lambda_a, Cp_a

     ! Gas phase mixture properties
     real(WP), dimension(:), pointer :: rho_g, visc_g, lambda_g, Cp_g, Sc_g, Pr_g

     ! Reference variables for evaporation model
     real(WP), dimension(:), pointer :: T_ref, Y_ref, rho_rv, visc_rv, lambda_rv, Cp_rv, G_rv, rho_ra, visc_ra, lambda_ra, Cp_ra

     ! Injector/Injection parameters
     character(len=128), pointer :: roi_file
     real(WP), dimension(:,:), pointer :: roi
     real(WP), pointer :: P_inj, C_d, U_inj

     ! Spray angle parameters
     character(len=128), pointer :: spray_angle_model
     real(WP), pointer :: theta, beta, Cnoz, C_theta

     ! Nozzle flow constants/variables
     real(WP), pointer :: K_in, K_exp, Cc0, Cc, D_eff, const_inj_vel

     ! Droplet size distribution
     character(len=128), pointer :: init_dsd_name
     integer, pointer :: nd
     real(WP), pointer :: h, init_dm, init_d2, init_d3, d10, d20, d30
     real(WP), dimension(:,:), pointer :: dsd, CD, Red, Shd, Nud, dsdlam
     real(WP), dimension(:), pointer :: di

     ! Universal constants
     real(WP), pointer :: R_gas

     ! Solver
     type(solver_t), pointer :: solver

     ! Output
     integer, pointer :: outfreq
     character(len=128) :: datafilename, plfilename

     ! Postprocessing
     real(WP), dimension(:), pointer :: time, LPL, VPL

     ! Fixed Non-dimensional numbers
     real(WP), pointer :: fixed_Re, fixed_We, fixed_DRa, fixed_DRv, fixed_VRa, fixed_VRv, fixed_De
  end type spray_t


  ! Definition of SolverVec object
  type SolverVec_t

     ! Flow variables
     real(WP), dimension(:), pointer :: z, rho, Y_l, Y_v, Y_a, Y_g, u_l, u_g, d2, dm, d3, Td, b, Tg

     ! Source terms
     real(WP), dimension(:), pointer :: omega_ent, omega_vap, omega_vapdm, omega_vapd2, omega_vapd3, &
                                        f_drag, omega_bre1, omega_bre2, omega_bre3, omega_T

     ! Solver
     type(solver_t), pointer :: solver

     ! Number of grid points
     integer, pointer :: Nz, nzo, step
     integer, pointer :: kmin, kmax, kmino, kmaxo

     ! Numerics
     real(WP), pointer :: dz, Lz

  end type SolverVec_t   

contains

  ! Allocate spray with default values
  subroutine allocate_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    !call nullify_spray(spray)

    allocate(spray)
    allocate(spray%Nz); spray%Nz = -9999
    allocate(spray%nzo); spray%nzo = -9999
    allocate(spray%step); spray%step = -9999
    allocate(spray%kmin); spray%kmin = -9999
    allocate(spray%kmax); spray%kmax = -9999
    allocate(spray%kmino); spray%kmino = -9999
    allocate(spray%kmaxo); spray%kmaxo = -9999

    allocate(spray%Lz); spray%Lz = -9999.0_WP
    allocate(spray%dz); spray%dz = -9999.0_WP
    allocate(spray%dt); spray%dt = -9999.0_WP
    allocate(spray%ndtime); spray%ndtime = -9999.0_WP
    allocate(spray%final_time); spray%final_time = -9999.0_WP
    allocate(spray%ndftime); spray%ndftime = -9999.0_WP
    allocate(spray%tau); spray%tau = -9999.0_WP
    allocate(spray%MaxCFL); spray%MaxCFL = -9999.0_WP
    allocate(spray%CFL); spray%CFL = -9999.0_WP

    allocate(spray%Re); spray%Re = -9999.0_WP
    allocate(spray%We); spray%We = -9999.0_WP
    allocate(spray%DRa); spray%DRa = -9999.0_WP
    allocate(spray%DRv); spray%DRv = -9999.0_WP
    allocate(spray%VRa); spray%VRa = -9999.0_WP
    allocate(spray%VRv); spray%VRv = -9999.0_WP
    allocate(spray%WR); spray%WR = -9999.0_WP
    allocate(spray%De); spray%De = -9999.0_WP

    allocate(spray%noz_D); spray%noz_D = -9999.0_WP
    allocate(spray%noz_LD); spray%noz_LD = -9999.0_WP
    allocate(spray%noz_rD); spray%noz_rD = -9999.0_WP
    allocate(spray%noz_Dsac); spray%noz_Dsac = -9999.0_WP
    allocate(spray%Cnoz); spray%Cnoz = -9999.0_WP
    allocate(spray%C_theta); spray%C_theta = -9999.0_WP
    allocate(spray%theta); spray%theta = -9999.0_WP
    allocate(spray%beta); spray%beta = -9999.0_WP
    allocate(spray%spray_angle_model); spray%spray_angle_model = 'noname'

    allocate(spray%Fuel); spray%Fuel = 'noname'
    allocate(spray%LFPTname); spray%LFPTname = 'noname'
    allocate(spray%VFPTname); spray%VFPTname = 'noname'

    allocate(spray%T_fuel); spray%T_fuel = -9999.0_WP
    allocate(spray%sigma); spray%sigma = -9999.0_WP
    allocate(spray%rho_l); spray%rho_l = -9999.0_WP
    allocate(spray%visc_l); spray%visc_l = -9999.0_WP
    allocate(spray%MW_f); spray%MW_f = -9999.0_WP
    allocate(spray%MP); spray%MP = -9999.0_WP
    allocate(spray%NBP); spray%NBP = -9999.0_WP

    allocate(spray%MW_a); spray%MW_a = -9999.0_WP
    allocate(spray%Z_a); spray%Z_a = -9999.0_WP
    allocate(spray%P_a); spray%P_a = -9999.0_WP
    allocate(spray%T_a); spray%T_a = -9999.0_WP

    allocate(spray%rho_a); spray%rho_a = -9999.0_WP
    allocate(spray%visc_a); spray%visc_a = -9999.0_WP
    allocate(spray%lambda_a); spray%lambda_a = -9999.0_WP
    allocate(spray%Cp_a); spray%Cp_a = -9999.0_WP

    allocate(spray%rho_v); spray%rho_v = -9999.0_WP
    allocate(spray%visc_v); spray%visc_v = -9999.0_WP
    allocate(spray%lambda_v); spray%lambda_v = -9999.0_WP
    allocate(spray%Cp_v); spray%Cp_v = -9999.0_WP

    allocate(spray%roi_file); spray%roi_file = 'noname'
    allocate(spray%P_inj); spray%P_inj = -9999.0_WP
    allocate(spray%C_d); spray%C_d = -9999.0_WP
    allocate(spray%U_inj); spray%U_inj = -9999.0_WP

    allocate(spray%K_in); spray%K_in = -9999.0_WP
    allocate(spray%K_exp); spray%K_exp = -9999.0_WP
    allocate(spray%Cc0); spray%Cc0 = -9999.0_WP
    allocate(spray%Cc); spray%Cc = -9999.0_WP
    allocate(spray%D_eff); spray%D_eff = -9999.0_WP
    allocate(spray%const_inj_vel); spray%const_inj_vel = -9999.0_WP

    allocate(spray%init_dsd_name); spray%init_dsd_name = 'noname'
    allocate(spray%nd); spray%nd = -9999
    allocate(spray%h); spray%h = -9999.0_WP
    allocate(spray%init_dm); spray%init_dm = -9999.0_WP
    allocate(spray%init_d2); spray%init_d2 = -9999.0_WP
    allocate(spray%init_d3); spray%init_d3 = -9999.0_WP

    allocate(spray%R_gas); spray%R_gas = -9999.0_WP

    allocate(spray%solver)

    allocate(spray%solver%rk)

    allocate(spray%solver%nr)
    ! Default initialization
    spray%solver%nr%tol = 1.0E-02_WP
    spray%solver%nr%relax_coeff = 2.0E-02_WP
    spray%solver%nr%alpha = 1.0E-04_WP
    spray%solver%nr%max_count = 2000

    allocate(spray%outfreq); spray%outfreq = -9999

    allocate(spray%fixed_Re); spray%fixed_Re = -9999_WP
    allocate(spray%fixed_We); spray%fixed_We = -9999_WP
    allocate(spray%fixed_DRa); spray%fixed_DRa = -9999_WP
    allocate(spray%fixed_DRv); spray%fixed_DRv = -9999_WP
    allocate(spray%fixed_VRa); spray%fixed_VRa = -9999_WP
    allocate(spray%fixed_VRv); spray%fixed_VRv = -9999_WP
    allocate(spray%fixed_De); spray%fixed_De = -9999_WP

  end subroutine allocate_spray

  ! Allocate spray with default values
  subroutine allocate_SolverVec(spray,SolverVecA,SolverVecB,SolverVecC,SolverVecD)
    implicit none

    ! ---------------------------------
    type(SolverVec_t), pointer, intent(inout) :: SolverVecA,SolverVecB,SolverVecC,SolverVecD
    type(spray_t), pointer, intent(inout) :: spray
    ! ---------------------------------


    allocate(SolverVecA)
    allocate(SolverVecA%dz); SolverVecA%dz = -9999
    allocate(SolverVecA%nzo); SolverVecA%nzo = -9999
    
    SolverVecA%dz = spray%Lz/(spray%Nz-1)
    SolverVecA%nzo = spray%Nz+6
    
    
    allocate(SolverVecA%z(SolverVecA%nzo)); SolverVecA%z = -9999.0_WP
    allocate(SolverVecA%rho(SolverVecA%nzo)); SolverVecA%rho = -9999.0_WP
    allocate(SolverVecA%Y_l(SolverVecA%nzo)); SolverVecA%Y_l = -9999.0_WP
    allocate(SolverVecA%Y_v(SolverVecA%nzo)); SolverVecA%Y_v = -9999.0_WP
    allocate(SolverVecA%Y_a(SolverVecA%nzo)); SolverVecA%Y_a = -9999.0_WP
    allocate(SolverVecA%Y_g(SolverVecA%nzo)); SolverVecA%Y_g = -9999.0_WP
    allocate(SolverVecA%u_l(SolverVecA%nzo)); SolverVecA%u_l = -9999.0_WP
    allocate(SolverVecA%u_g(SolverVecA%nzo)); SolverVecA%u_g = -9999.0_WP
    allocate(SolverVecA%d3(SolverVecA%nzo)); SolverVecA%d3 = -9999.0_WP
    allocate(SolverVecA%d2(SolverVecA%nzo)); SolverVecA%d2 = -9999.0_WP
    allocate(SolverVecA%dm(SolverVecA%nzo)); SolverVecA%dm = -9999.0_WP
    allocate(SolverVecA%Td(SolverVecA%nzo)); SolverVecA%Td = -9999.0_WP
    allocate(SolverVecA%Tg(SolverVecA%nzo)); SolverVecA%Tg = -9999.0_WP
    allocate(SolverVecA%b(SolverVecA%nzo)); SolverVecA%b = -9999.0_WP

    allocate(SolverVecA%omega_ent(SolverVecA%nzo)); SolverVecA%omega_ent = -9999.0_WP
    allocate(SolverVecA%omega_vap(SolverVecA%nzo)); SolverVecA%omega_vap = -9999.0_WP
    allocate(SolverVecA%omega_vapdm(SolverVecA%nzo)); SolverVecA%omega_vapdm = -9999.0_WP
    allocate(SolverVecA%omega_vapd2(SolverVecA%nzo)); SolverVecA%omega_vapd2 = -9999.0_WP
    allocate(SolverVecA%omega_vapd3(SolverVecA%nzo)); SolverVecA%omega_vapd3 = -9999.0_WP
    allocate(SolverVecA%f_drag(SolverVecA%nzo)); SolverVecA%f_drag = -9999.0_WP
    allocate(SolverVecA%omega_bre1(SolverVecA%nzo)); SolverVecA%omega_bre1 = -9999.0_WP
    allocate(SolverVecA%omega_bre2(SolverVecA%nzo)); SolverVecA%omega_bre2 = -9999.0_WP
    allocate(SolverVecA%omega_bre3(SolverVecA%nzo)); SolverVecA%omega_bre3 = -9999.0_WP
    allocate(SolverVecA%omega_T(SolverVecA%nzo)); SolverVecA%omega_T = -9999.0_WP

    allocate(SolverVecA%solver)

    !call allocate_solver(SolverVecA%solver,spray%nzo)

    allocate(SolverVecA%solver%rk)

    allocate(SolverVecA%solver%nr)
    ! Default initialization
    SolverVecA%solver%nr%tol = 1.0E-02_WP
    SolverVecA%solver%nr%relax_coeff = 2.0E-02_WP
    SolverVecA%solver%nr%alpha = 1.0E-04_WP
    SolverVecA%solver%nr%max_count = 2000


    allocate(SolverVecB)
    allocate(SolverVecB%dz); SolverVecB%dz = -9999
    allocate(SolverVecB%nzo); SolverVecB%nzo = -9999
 
    SolverVecB%dz = spray%Lz/(spray%Nz-1)
    SolverVecB%nzo = spray%Nz+6
    
    allocate(SolverVecB%z(SolverVecB%nzo)); SolverVecB%z = -9999.0_WP
    allocate(SolverVecB%rho(SolverVecB%nzo)); SolverVecB%rho = -9999.0_WP
    allocate(SolverVecB%Y_l(SolverVecB%nzo)); SolverVecB%Y_l = -9999.0_WP
    allocate(SolverVecB%Y_v(SolverVecB%nzo)); SolverVecB%Y_v = -9999.0_WP
    allocate(SolverVecB%Y_a(SolverVecB%nzo)); SolverVecB%Y_a = -9999.0_WP
    allocate(SolverVecB%Y_g(SolverVecB%nzo)); SolverVecB%Y_g = -9999.0_WP
    allocate(SolverVecB%u_l(SolverVecB%nzo)); SolverVecB%u_l = -9999.0_WP
    allocate(SolverVecB%u_g(SolverVecB%nzo)); SolverVecB%u_g = -9999.0_WP
    allocate(SolverVecB%d3(SolverVecB%nzo)); SolverVecB%d3 = -9999.0_WP
    allocate(SolverVecB%d2(SolverVecB%nzo)); SolverVecB%d2 = -9999.0_WP
    allocate(SolverVecB%dm(SolverVecB%nzo)); SolverVecB%dm = -9999.0_WP
    allocate(SolverVecB%Td(SolverVecB%nzo)); SolverVecB%Td = -9999.0_WP
    allocate(SolverVecB%Tg(SolverVecB%nzo)); SolverVecB%Tg = -9999.0_WP
    allocate(SolverVecB%b(SolverVecB%nzo)); SolverVecB%b = -9999.0_WP

    allocate(SolverVecB%omega_ent(SolverVecB%nzo)); SolverVecB%omega_ent = -9999.0_WP
    allocate(SolverVecB%omega_vap(SolverVecB%nzo)); SolverVecB%omega_vap = -9999.0_WP
    allocate(SolverVecB%omega_vapdm(SolverVecB%nzo)); SolverVecB%omega_vapdm = -9999.0_WP
    allocate(SolverVecB%omega_vapd2(SolverVecB%nzo)); SolverVecB%omega_vapd2 = -9999.0_WP
    allocate(SolverVecB%omega_vapd3(SolverVecB%nzo)); SolverVecB%omega_vapd3 = -9999.0_WP
    allocate(SolverVecB%f_drag(SolverVecB%nzo)); SolverVecB%f_drag = -9999.0_WP
    allocate(SolverVecB%omega_bre1(SolverVecB%nzo)); SolverVecB%omega_bre1 = -9999.0_WP
    allocate(SolverVecB%omega_bre2(SolverVecB%nzo)); SolverVecB%omega_bre2 = -9999.0_WP
    allocate(SolverVecB%omega_bre3(SolverVecB%nzo)); SolverVecB%omega_bre3 = -9999.0_WP
    allocate(SolverVecB%omega_T(SolverVecB%nzo)); SolverVecB%omega_T = -9999.0_WP

    allocate(SolverVecB%solver)

    !call allocate_solver(SolverVecB%solver,spray%nzo)

    allocate(SolverVecB%solver%rk)

    allocate(SolverVecB%solver%nr)
    ! Default initialization
    SolverVecB%solver%nr%tol = 1.0E-02_WP
    SolverVecB%solver%nr%relax_coeff = 2.0E-02_WP
    SolverVecB%solver%nr%alpha = 1.0E-04_WP
    SolverVecB%solver%nr%max_count = 2000

    allocate(SolverVecC)
    allocate(SolverVecC%dz); SolverVecC%dz = -9999
    allocate(SolverVecC%nzo); SolverVecC%nzo = -9999

    SolverVecC%dz = spray%Lz/(spray%Nz-1)
    SolverVecC%nzo = spray%Nz+6

    allocate(SolverVecC%z(SolverVecC%nzo)); SolverVecC%z = -9999.0_WP
    allocate(SolverVecC%rho(SolverVecC%nzo)); SolverVecC%rho = -9999.0_WP
    allocate(SolverVecC%Y_l(SolverVecC%nzo)); SolverVecC%Y_l = -9999.0_WP
    allocate(SolverVecC%Y_v(SolverVecC%nzo)); SolverVecC%Y_v = -9999.0_WP
    allocate(SolverVecC%Y_a(SolverVecC%nzo)); SolverVecC%Y_a = -9999.0_WP
    allocate(SolverVecC%Y_g(SolverVecC%nzo)); SolverVecC%Y_g = -9999.0_WP
    allocate(SolverVecC%u_l(SolverVecC%nzo)); SolverVecC%u_l = -9999.0_WP
    allocate(SolverVecC%u_g(SolverVecC%nzo)); SolverVecC%u_g = -9999.0_WP
    allocate(SolverVecC%d3(SolverVecC%nzo)); SolverVecC%d3 = -9999.0_WP
    allocate(SolverVecC%d2(SolverVecC%nzo)); SolverVecC%d2 = -9999.0_WP
    allocate(SolverVecC%dm(SolverVecC%nzo)); SolverVecC%dm = -9999.0_WP
    allocate(SolverVecC%Td(SolverVecC%nzo)); SolverVecC%Td = -9999.0_WP
    allocate(SolverVecC%Tg(SolverVecC%nzo)); SolverVecC%Tg = -9999.0_WP
    allocate(SolverVecC%b(SolverVecC%nzo)); SolverVecC%b = -9999.0_WP

    allocate(SolverVecC%omega_ent(SolverVecC%nzo)); SolverVecC%omega_ent = -9999.0_WP
    allocate(SolverVecC%omega_vap(SolverVecC%nzo)); SolverVecC%omega_vap = -9999.0_WP
    allocate(SolverVecC%omega_vapdm(SolverVecC%nzo)); SolverVecC%omega_vapdm = -9999.0_WP
    allocate(SolverVecC%omega_vapd2(SolverVecC%nzo)); SolverVecC%omega_vapd2 = -9999.0_WP
    allocate(SolverVecC%omega_vapd3(SolverVecC%nzo)); SolverVecC%omega_vapd3 = -9999.0_WP
    allocate(SolverVecC%f_drag(SolverVecC%nzo)); SolverVecC%f_drag = -9999.0_WP
    allocate(SolverVecC%omega_bre1(SolverVecC%nzo)); SolverVecC%omega_bre1 = -9999.0_WP
    allocate(SolverVecC%omega_bre2(SolverVecC%nzo)); SolverVecC%omega_bre2 = -9999.0_WP
    allocate(SolverVecC%omega_bre3(SolverVecC%nzo)); SolverVecC%omega_bre3 = -9999.0_WP
    allocate(SolverVecC%omega_T(SolverVecC%nzo)); SolverVecC%omega_T = -9999.0_WP

    allocate(SolverVecC%solver)

    !call allocate_solver(SolverVecC%solver,spray%nzo)

    allocate(SolverVecC%solver%rk)

    allocate(SolverVecC%solver%nr)
    ! Default initialization
    SolverVecC%solver%nr%tol = 1.0E-02_WP
    SolverVecC%solver%nr%relax_coeff = 2.0E-02_WP
    SolverVecC%solver%nr%alpha = 1.0E-04_WP
    SolverVecC%solver%nr%max_count = 2000

    allocate(SolverVecD)
    allocate(SolverVecD%dz); SolverVecD%dz = -9999
    allocate(SolverVecD%nzo); SolverVecD%nzo = -9999

    SolverVecD%dz = spray%Lz/(spray%Nz-1)
    SolverVecD%nzo = spray%Nz+6

    allocate(SolverVecD%z(SolverVecD%nzo)); SolverVecD%z = -9999.0_WP
    allocate(SolverVecD%rho(SolverVecD%nzo)); SolverVecD%rho = -9999.0_WP
    allocate(SolverVecD%Y_l(SolverVecD%nzo)); SolverVecD%Y_l = -9999.0_WP
    allocate(SolverVecD%Y_v(SolverVecD%nzo)); SolverVecD%Y_v = -9999.0_WP
    allocate(SolverVecD%Y_a(SolverVecD%nzo)); SolverVecD%Y_a = -9999.0_WP
    allocate(SolverVecD%Y_g(SolverVecD%nzo)); SolverVecD%Y_g = -9999.0_WP
    allocate(SolverVecD%u_l(SolverVecD%nzo)); SolverVecD%u_l = -9999.0_WP
    allocate(SolverVecD%u_g(SolverVecD%nzo)); SolverVecD%u_g = -9999.0_WP
    allocate(SolverVecD%d3(SolverVecD%nzo)); SolverVecD%d3 = -9999.0_WP
    allocate(SolverVecD%d2(SolverVecD%nzo)); SolverVecD%d2 = -9999.0_WP
    allocate(SolverVecD%dm(SolverVecD%nzo)); SolverVecD%dm = -9999.0_WP
    allocate(SolverVecD%Td(SolverVecD%nzo)); SolverVecD%Td = -9999.0_WP
    allocate(SolverVecD%Tg(SolverVecD%nzo)); SolverVecD%Tg = -9999.0_WP
    allocate(SolverVecD%b(SolverVecD%nzo)); SolverVecD%b = -9999.0_WP

    allocate(SolverVecD%omega_ent(SolverVecD%nzo)); SolverVecD%omega_ent = -9999.0_WP
    allocate(SolverVecD%omega_vap(SolverVecD%nzo)); SolverVecD%omega_vap = -9999.0_WP
    allocate(SolverVecD%omega_vapdm(SolverVecD%nzo)); SolverVecD%omega_vapdm = -9999.0_WP
    allocate(SolverVecD%omega_vapd2(SolverVecD%nzo)); SolverVecD%omega_vapd2 = -9999.0_WP
    allocate(SolverVecD%omega_vapd3(SolverVecD%nzo)); SolverVecD%omega_vapd3 = -9999.0_WP
    allocate(SolverVecD%f_drag(SolverVecD%nzo)); SolverVecD%f_drag = -9999.0_WP
    allocate(SolverVecD%omega_bre1(SolverVecD%nzo)); SolverVecD%omega_bre1 = -9999.0_WP
    allocate(SolverVecD%omega_bre2(SolverVecD%nzo)); SolverVecD%omega_bre2 = -9999.0_WP
    allocate(SolverVecD%omega_bre3(SolverVecD%nzo)); SolverVecD%omega_bre3 = -9999.0_WP
    allocate(SolverVecD%omega_T(SolverVecD%nzo)); SolverVecD%omega_T = -9999.0_WP

    allocate(SolverVecD%solver)

    !call allocate_solver(SolverVecD%solver,spray%nzo)

    allocate(SolverVecD%solver%rk)

    allocate(SolverVecD%solver%nr)
    ! Default initialization
    SolverVecD%solver%nr%tol = 1.0E-02_WP
    SolverVecD%solver%nr%relax_coeff = 2.0E-02_WP
    SolverVecD%solver%nr%alpha = 1.0E-04_WP
    SolverVecD%solver%nr%max_count = 2000

    call allocate_solver(SolverVecA%solver,spray%nzo)
    call allocate_solver(SolverVecB%solver,spray%nzo)
    call allocate_solver(SolverVecC%solver,spray%nzo)
    call allocate_solver(SolverVecD%solver,spray%nzo)

   end subroutine allocate_SolverVec

   subroutine assign_SolverVec(SolverVecA,SolverVecB,SolverVecC,SolverVecD,spray)
     implicit none

     ! ---------------------------------
     type(spray_t), pointer, intent(inout) :: spray
     type(SolverVec_t), pointer, intent(inout) :: SolverVecA,SolverVecB,SolverVecC,SolverVecD

     IF (spray%solver%rk%niter == 1) THEN
        SolverVecA%z = spray%z
        SolverVecA%rho = spray%rho
        SolverVecA%Y_l = spray%Y_l
        SolverVecA%Y_v = spray%Y_v
        SolverVecA%Y_a = spray%Y_a
        SolverVecA%Y_g = spray%Y_g
        SolverVecA%u_l = spray%u_l
        SolverVecA%u_g = spray%u_g
        SolverVecA%d3 = spray%d3
        SolverVecA%d2 = spray%d2
        SolverVecA%dm = spray%dm
        SolverVecA%Td = spray%Td
        SolverVecA%Tg = spray%Tg
        SolverVecA%b = spray%b

        SolverVecA%omega_ent = spray%omega_ent
        SolverVecA%omega_vap = spray%omega_vap
        SolverVecA%omega_vapdm = spray%omega_vapdm
        SolverVecA%omega_vapd2 = spray%omega_vapd2
        SolverVecA%omega_vapd3 = spray%omega_vapd3
        SolverVecA%f_drag = spray%f_drag
        SolverVecA%omega_bre1 = spray%omega_bre1
        SolverVecA%omega_bre2 = spray%omega_bre2
        SolverVecA%omega_bre3 = spray%omega_bre3
        SolverVecA%omega_T = spray%omega_T
        !SolverVecA%solver = spray%solver

     ELSE IF (spray%solver%rk%niter == 2) THEN

        SolverVecB%z = spray%z
        SolverVecB%rho = spray%rho
        SolverVecB%Y_l = spray%Y_l
        SolverVecB%Y_v = spray%Y_v
        SolverVecB%Y_a = spray%Y_a
        SolverVecB%Y_g = spray%Y_g
        SolverVecB%u_l = spray%u_l
        SolverVecB%u_g = spray%u_g
        SolverVecB%d3 = spray%d3
        SolverVecB%d2 = spray%d2
        SolverVecB%dm = spray%dm
        SolverVecB%Td = spray%Td
        SolverVecB%Tg = spray%Tg
        SolverVecB%b = spray%b

        SolverVecB%omega_ent = spray%omega_ent
        SolverVecB%omega_vap = spray%omega_vap
        SolverVecB%omega_vapdm = spray%omega_vapdm
        SolverVecB%omega_vapd2 = spray%omega_vapd2
        SolverVecB%omega_vapd3 = spray%omega_vapd3
        SolverVecB%f_drag = spray%f_drag
        SolverVecB%omega_bre1 = spray%omega_bre1
        SolverVecB%omega_bre2 = spray%omega_bre2
        SolverVecB%omega_bre3 = spray%omega_bre3
        SolverVecB%omega_T = spray%omega_T
        !SolverVecB%solver = spray%solver

     ELSE IF (spray%solver%rk%niter == 3) THEN

        SolverVecC%z = spray%z
        SolverVecC%rho = spray%rho
        SolverVecC%Y_l = spray%Y_l
        SolverVecC%Y_v = spray%Y_v
        SolverVecC%Y_a = spray%Y_a
        SolverVecC%Y_g = spray%Y_g
        SolverVecC%u_l = spray%u_l
        SolverVecC%u_g = spray%u_g
        SolverVecC%d3 = spray%d3
        SolverVecC%d2 = spray%d2
        SolverVecC%dm = spray%dm
        SolverVecC%Td = spray%Td
        SolverVecC%Tg = spray%Tg
        SolverVecC%b = spray%b

        SolverVecC%omega_ent = spray%omega_ent
        SolverVecC%omega_vap = spray%omega_vap
        SolverVecC%omega_vapdm = spray%omega_vapdm
        SolverVecC%omega_vapd2 = spray%omega_vapd2
        SolverVecC%omega_vapd3 = spray%omega_vapd3
        SolverVecC%f_drag = spray%f_drag
        SolverVecC%omega_bre1 = spray%omega_bre1
        SolverVecC%omega_bre2 = spray%omega_bre2
        SolverVecC%omega_bre3 = spray%omega_bre3
        SolverVecC%omega_T = spray%omega_T
        !SolverVecC%solver = spray%solver

     ELSE IF (spray%solver%rk%niter == 4) THEN

        SolverVecD%z = spray%z
        SolverVecD%rho = spray%rho
        SolverVecD%Y_l = spray%Y_l
        SolverVecD%Y_v = spray%Y_v
        SolverVecD%Y_a = spray%Y_a
        SolverVecD%Y_g = spray%Y_g
        SolverVecD%u_l = spray%u_l
        SolverVecD%u_g = spray%u_g
        SolverVecD%d3 = spray%d3
        SolverVecD%d2 = spray%d2
        SolverVecD%dm = spray%dm
        SolverVecD%Td = spray%Td
        SolverVecD%Tg = spray%Tg
        SolverVecD%b = spray%b

        SolverVecD%omega_ent = spray%omega_ent
        SolverVecD%omega_vap = spray%omega_vap
        SolverVecD%omega_vapdm = spray%omega_vapdm
        SolverVecD%omega_vapd2 = spray%omega_vapd2
        SolverVecD%omega_vapd3 = spray%omega_vapd3
        SolverVecD%f_drag = spray%f_drag
        SolverVecD%omega_bre1 = spray%omega_bre1
        SolverVecD%omega_bre2 = spray%omega_bre2
        SolverVecD%omega_bre3 = spray%omega_bre3
        SolverVecD%omega_T = spray%omega_T
        !SolverVecD%solver = spray%solver


     END IF


   end subroutine assign_SolverVec

  subroutine allocate_spray_grid_vars(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l, pc_v
    integer :: k

    spray%dz = spray%Lz/(spray%Nz-1)
    spray%nzo = spray%Nz+6
    spray%kmino = 1; spray%kmaxo = spray%nzo
    spray%kmin  = 4; spray%kmax  = spray%Nz

    allocate(spray%pc_l(spray%nzo))
    allocate(spray%pc_v(spray%nzo))

    do k = spray%kmino,spray%kmaxo
       pc_l => spray%pc_l(k)
       pc_v => spray%pc_v(k)

       call allocate_pc(pc_l)
       call allocate_pc(pc_v)

       spray%pc_l(k) = pc_l
       spray%pc_v(k) = pc_v

       nullify(pc_l,pc_v)

    end do

    allocate(spray%rho_ra(spray%nzo)); spray%rho_ra = -9999.0_WP
    allocate(spray%visc_ra(spray%nzo)); spray%visc_ra = -9999.0_WP
    allocate(spray%lambda_ra(spray%nzo)); spray%lambda_ra = -9999.0_WP
    allocate(spray%Cp_ra(spray%nzo)); spray%Cp_ra = -9999.0_WP

    allocate(spray%rho_rv(spray%nzo)); spray%rho_rv = -9999.0_WP
    allocate(spray%visc_rv(spray%nzo)); spray%visc_rv = -9999.0_WP
    allocate(spray%lambda_rv(spray%nzo)); spray%lambda_rv = -9999.0_WP
    allocate(spray%Cp_rv(spray%nzo)); spray%Cp_rv = -9999.0_WP
    allocate(spray%G_rv(spray%nzo)); spray%G_rv = -9999.0_WP

    allocate(spray%rho_g(spray%nzo)); spray%rho_g = -9999.0_WP
    allocate(spray%visc_g(spray%nzo)); spray%visc_g = -9999.0_WP
    allocate(spray%lambda_g(spray%nzo)); spray%lambda_g = -9999.0_WP
    allocate(spray%Cp_g(spray%nzo)); spray%Cp_g = -9999.0_WP
    allocate(spray%Sc_g(spray%nzo)); spray%Sc_g = -9999.0_WP
    allocate(spray%Pr_g(spray%nzo)); spray%Pr_g = -9999.0_WP

    allocate(spray%C_l(spray%nzo)); spray%C_l = -9999.0_WP
    allocate(spray%p_vap(spray%nzo)); spray%p_vap = -9999.0_WP
    allocate(spray%L_f(spray%nzo)); spray%L_f = -9999.0_WP

    allocate(spray%T_ref(spray%nzo)); spray%T_ref = -9999.0_WP
    allocate(spray%Y_ref(spray%nzo)); spray%Y_ref = -9999.0_WP

    allocate(spray%DRg(spray%nzo)); spray%DRg = -9999.0_WP
    allocate(spray%VRg(spray%nzo)); spray%VRg = -9999.0_WP
    allocate(spray%LR(spray%nzo)); spray%LR = -9999.0_WP
    allocate(spray%CR(spray%nzo)); spray%CR = -9999.0_WP

    allocate(spray%z(spray%nzo)); spray%z = -9999.0_WP
    allocate(spray%rho(spray%nzo)); spray%rho = -9999.0_WP
    allocate(spray%Y_l(spray%nzo)); spray%Y_l = -9999.0_WP
    allocate(spray%Y_v(spray%nzo)); spray%Y_v = -9999.0_WP
    allocate(spray%Y_a(spray%nzo)); spray%Y_a = -9999.0_WP
    allocate(spray%Y_g(spray%nzo)); spray%Y_g = -9999.0_WP
    allocate(spray%u_l(spray%nzo)); spray%u_l = -9999.0_WP
    allocate(spray%u_g(spray%nzo)); spray%u_g = -9999.0_WP
    allocate(spray%d3(spray%nzo)); spray%d3 = -9999.0_WP
    allocate(spray%d2(spray%nzo)); spray%d2 = -9999.0_WP
    allocate(spray%dm(spray%nzo)); spray%dm = -9999.0_WP
    allocate(spray%Td(spray%nzo)); spray%Td = -9999.0_WP
    allocate(spray%Tg(spray%nzo)); spray%Tg = -9999.0_WP
    allocate(spray%b(spray%nzo)); spray%b = -9999.0_WP

    allocate(spray%dsd(spray%nd,spray%nzo)); spray%dsd = -9999.0_WP
    allocate(spray%CD(spray%nd,spray%nzo)); spray%CD = -9999.0_WP
    allocate(spray%Red(spray%nd,spray%nzo)); spray%Red = -9999.0_WP
    allocate(spray%Shd(spray%nd,spray%nzo)); spray%Shd = -9999.0_WP
    allocate(spray%Nud(spray%nd,spray%nzo)); spray%Nud = -9999.0_WP
    allocate(spray%di(spray%nd)); spray%di = -9999.0_WP
    allocate(spray%dsdlam(4,spray%nzo)); spray%dsdlam = -9999.0_WP

    call allocate_solver(spray%solver,spray%nzo)

    allocate(spray%omega_ent(spray%nzo)); spray%omega_ent = -9999.0_WP
    allocate(spray%omega_vap(spray%nzo)); spray%omega_vap = -9999.0_WP
    allocate(spray%omega_vapdm(spray%nzo)); spray%omega_vapdm = -9999.0_WP
    allocate(spray%omega_vapd2(spray%nzo)); spray%omega_vapd2 = -9999.0_WP
    allocate(spray%omega_vapd3(spray%nzo)); spray%omega_vapd3 = -9999.0_WP
    allocate(spray%f_drag(spray%nzo)); spray%f_drag = -9999.0_WP
    allocate(spray%omega_bre1(spray%nzo)); spray%omega_bre1 = -9999.0_WP
    allocate(spray%omega_bre2(spray%nzo)); spray%omega_bre2 = -9999.0_WP
    allocate(spray%omega_bre3(spray%nzo)); spray%omega_bre3 = -9999.0_WP
    allocate(spray%omega_T(spray%nzo)); spray%omega_T = -9999.0_WP

  end subroutine allocate_spray_grid_vars

  subroutine deallocate_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    deallocate(spray%Nz,spray%nzo,spray%step,spray%kmin,spray%kmax,spray%kmino,spray%kmaxo)

    deallocate(spray%z,spray%Lz,spray%dz)

    deallocate(spray%dt,spray%ndtime,spray%final_time,spray%ndftime,spray%tau)

    deallocate(spray%MaxCFL,spray%CFL)

    deallocate(spray%Re,spray%We,spray%DRa,spray%DRv,spray%DRg,spray%VRa,spray%VRv,spray%VRg)
    
    deallocate(spray%WR,spray%De,spray%LR,spray%CR)

    deallocate(spray%noz_D,spray%noz_LD,spray%noz_rD,spray%noz_Dsac)
    
    deallocate(spray%Cnoz,spray%C_theta,spray%theta,spray%beta,spray%spray_angle_model)

    deallocate(spray%rho,spray%Y_l,spray%Y_v,spray%Y_a, spray%Y_g, &
               spray%u_l,spray%u_g, &
               spray%d3,spray%d2,spray%dm, &
               spray%Td,spray%b,spray%Tg)

    deallocate(spray%Fuel,spray%T_fuel,spray%sigma,spray%rho_l,spray%visc_l,spray%C_l, &
               spray%p_vap,spray%MW_f,spray%L_f,spray%MP,spray%NBP)

    deallocate(spray%pc_l)

    deallocate(spray%pc_v)

    if (associated(spray%LFPT)) then
       deallocate(spray%LFPT)
    end if

    if (associated(spray%VFPT)) then
       deallocate(spray%VFPT)
    end if
    
    deallocate(spray%rho_v, spray%visc_v, spray%lambda_v, spray%Cp_v)

    deallocate(spray%rho_rv, spray%visc_rv, spray%lambda_rv, spray%Cp_rv, spray%G_rv)
    
    deallocate(spray%MW_a, spray%Z_a, spray%P_a, spray%T_a, spray%rho_a, spray%visc_a, spray%lambda_a, spray%Cp_a)

    deallocate(spray%rho_ra, spray%visc_ra, spray%lambda_ra, spray%Cp_ra)

    deallocate(spray%rho_g, spray%visc_g, spray%lambda_g, spray%Cp_g, spray%Sc_g, spray%Pr_g)

    deallocate(spray%T_ref, spray%Y_ref)

    deallocate(spray%roi_file, spray%roi)

    deallocate(spray%P_inj, spray%C_d, spray%U_inj)

    deallocate(spray%K_in, spray%K_exp, spray%Cc0, spray%Cc, spray%D_eff)

    deallocate(spray%init_dsd_name)

    deallocate(spray%nd)

    deallocate(spray%h, spray%init_dm, spray%init_d2, spray%init_d3)

    deallocate(spray%dsd, spray%CD, spray%Red, spray%Shd, spray%Nud, spray%dsdlam)

    deallocate(spray%di)

    deallocate(spray%R_gas)

    call deallocate_solver(spray%solver)

    deallocate(spray%omega_ent, spray%omega_vap, spray%omega_vapdm, spray%omega_vapd2, spray%omega_vapd3, &
               spray%f_drag, spray%omega_bre1, spray%omega_bre2, spray%omega_bre3, spray%omega_T)

    deallocate(spray%outfreq)

    deallocate(spray%time, spray%LPL, spray%VPL)

    deallocate(spray%fixed_Re, spray%fixed_We,  spray%fixed_DRa, spray%fixed_DRv, spray%fixed_VRa, spray%fixed_VRv, spray%fixed_De)

    deallocate(spray)

  end subroutine deallocate_spray

  subroutine deallocate_SolverVec(SolverVecA,SolverVecB,SolverVecC,SolverVecD)
    implicit none

    ! ---------------------------------
    type(SolverVec_t), pointer, intent(inout) :: SolverVecA,SolverVecB,SolverVecC,SolverVecD

    ! ---------------------------------

    deallocate(SolverVecA%z,SolverVecB%z)

    deallocate(SolverVecA%rho,SolverVecA%Y_l,SolverVecA%Y_v,SolverVecA%Y_a, SolverVecA%Y_g, &
               SolverVecA%u_l,SolverVecA%u_g, &
               SolverVecA%d3,SolverVecA%d2,SolverVecA%dm, &
               SolverVecA%Td,SolverVecA%b,SolverVecA%Tg)
    
    deallocate(SolverVecA%omega_ent, SolverVecA%omega_vap, SolverVecA%omega_vapdm, SolverVecA%omega_vapd2, SolverVecA%omega_vapd3, &
               SolverVecA%f_drag, SolverVecA%omega_bre1, SolverVecA%omega_bre2, SolverVecA%omega_bre3, SolverVecA%omega_T)

    deallocate(SolverVecB%rho,SolverVecB%Y_l,SolverVecB%Y_v,SolverVecB%Y_a, SolverVecB%Y_g, &
               SolverVecB%u_l,SolverVecB%u_g, &
               SolverVecB%d3,SolverVecB%d2,SolverVecB%dm, &
               SolverVecB%Td,SolverVecB%b,SolverVecB%Tg)
    
    deallocate(SolverVecB%omega_ent, SolverVecB%omega_vap, SolverVecB%omega_vapdm, SolverVecB%omega_vapd2, SolverVecB%omega_vapd3, &
               SolverVecB%f_drag, SolverVecB%omega_bre1, SolverVecB%omega_bre2, SolverVecB%omega_bre3, SolverVecB%omega_T)

    deallocate(SolverVecC%rho,SolverVecC%Y_l,SolverVecC%Y_v,SolverVecC%Y_a, SolverVecC%Y_g, &
               SolverVecC%u_l,SolverVecC%u_g, &
               SolverVecC%d3,SolverVecC%d2,SolverVecC%dm, &
               SolverVecC%Td,SolverVecC%b,SolverVecC%Tg)
    
    deallocate(SolverVecC%omega_ent, SolverVecC%omega_vap, SolverVecC%omega_vapdm, SolverVecC%omega_vapd2, SolverVecC%omega_vapd3, &
               SolverVecC%f_drag, SolverVecC%omega_bre1, SolverVecC%omega_bre2, SolverVecC%omega_bre3, SolverVecC%omega_T)

    deallocate(SolverVecD%rho,SolverVecD%Y_l,SolverVecD%Y_v,SolverVecD%Y_a, SolverVecD%Y_g, &
               SolverVecD%u_l,SolverVecD%u_g, &
               SolverVecD%d3,SolverVecD%d2,SolverVecD%dm, &
               SolverVecD%Td,SolverVecD%b,SolverVecD%Tg)
    
    deallocate(SolverVecD%omega_ent, SolverVecD%omega_vap, SolverVecD%omega_vapdm, SolverVecD%omega_vapd2, SolverVecD%omega_vapd3, &
               SolverVecD%f_drag, SolverVecD%omega_bre1, SolverVecD%omega_bre2, SolverVecD%omega_bre3, SolverVecD%omega_T)

    call deallocate_solver(SolverVecA%solver)
    call deallocate_solver(SolverVecB%solver)
    call deallocate_solver(SolverVecC%solver)
    call deallocate_solver(SolverVecD%solver)

  end subroutine deallocate_SolverVec
end module spray_defs
