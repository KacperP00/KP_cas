module spray_defs
  use precision
  use pc_defs
  use solver_defs

  implicit none

  ! Definition of spray object
  type spray_t
     
     ! Input
     character(len=128), pointer :: inp_fname

     ! Number of grid points
     integer, pointer :: Nz, nzo, step
     integer, pointer :: kmin, kmax, kmino, kmaxo

     ! Numerics
     real(WP), pointer :: dz, Lz
     real(WP), pointer :: dt, ndtime, final_time, ndftime, tau
     real(WP), pointer :: MaxCFL, CFL
     real(WP), dimension(:), pointer :: CFL_conv, CFL_bre, CFL_evap

     ! Non-dimensional numbers
     real(WP), pointer :: Re, We, WR, De, DRa, DRv, VRa, VRv, Pr_l
     real(WP), dimension(:), pointer :: DRg, VRg, VRtg, TCRg, LR, CR
     
     ! Geometric variables
     real(WP), pointer :: noz_D, noz_LD, noz_rD, noz_Dsac
     integer, pointer :: num_noz     

     ! Flow variables
     real(WP), dimension(:), pointer :: z, rho, Y_l, Y_v, Y_a, Y_g, u_l, u_g, d2, dm, d3, dvar, Td, b, Tg, &
                                        k_g, eps_g, mu_t_g, zmix_g, zvar_g, chi_g, chi_g_stl

     ! Source terms
     real(WP), dimension(:), pointer :: omega_ent, omega_vap, omega_vapdm, omega_vapd2, omega_vapd3, &
                                        f_drag, omega_bre1, omega_bre2, omega_bre3, omega_T, omega_k_g_p, &
                                        omega_k_g_d, omega_eps_g_p, omega_eps_g_d, &
                                        omega_zvar_g_p, omega_zvar_g_d

     ! Fuel name
     character(len=128), pointer :: Fuel

     ! Pure chemical object for liquid phase
     type(pc_t), dimension(:), pointer :: pc_l
     
     ! Pure chemical object for vapor phase
     type(pc_t), dimension(:), pointer :: pc_v
 
     ! Fuel properties
     real(WP), pointer :: T_fuel, sigma, rho_l, visc_l, lambda_l, MW_f, MP, NBP, stoic_coeff
     real(WP), dimension(:), pointer :: L_f, C_l, p_vap, T_sat, sigma_loc, rho_l_loc, visc_l_loc, lambda_l_loc

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
     real(WP), pointer :: P_inj, C_d, Anoz, U_inj, inj_mass
     logical :: noz_flow_model = .false.

     ! Spray angle parameters
     character(len=128), pointer :: spray_angle_model
     real(WP), pointer :: theta, beta, Cnoz, C_theta

     ! Nozzle flow constants/variables
     real(WP), pointer :: K_in, K_exp, Cc0, Cc, D_eff, const_inj_vel

     ! Droplet size distribution
     character(len=128), pointer :: init_dsd_name
     integer, pointer :: nd
     real(WP), pointer :: init_dm, init_d2, init_d3, init_dvar, d10, d20, d30, dgf
     real(WP), dimension(:), pointer :: h
     real(WP), dimension(:,:), pointer :: di, dsd, CD, Red, Shd, Nud, dsdlam
     integer, dimension(:), pointer :: dsd_type
     integer, pointer :: skip_d2, skip_d3

     ! Breakup parameters with default values
     real(WP) :: B0 = 0.61_WP, B1 = 10.0_WP, C3 = 2.5_WP, Crel = 1.0_WP

     ! Evaporation parameters
     real(WP) :: Cevap = 1.0_WP

     ! Universal constants
     real(WP), pointer :: R_gas

     ! Turbulence parameters
     logical :: turb_model = .false.
     integer, pointer :: skip_turb
     real(WP), pointer :: c_k, c_mu, c_eps1, c_eps2, c_zvar

     ! Combustion model
     character(len=128), pointer :: combustion_model
     real(WP) :: C_chi = 1.0_WP, Y_O2 = 0.164215_WP
     real(WP), dimension(:), pointer :: zz=>null(), bpdf=>null()

     ! Solver
     type(solver_t), pointer :: solver

     ! Output
     integer, pointer :: outfreq, datafreq
     character(len=128) :: datafilename, plfilename
     logical :: saveDataFile=.false.

     ! Postprocessing
     real(WP), dimension(:), pointer :: time, LPL, VPL, chi_st

     ! For sensitivity analysis
     real(WP) :: f_rho_l=0.0_WP, f_mu_l=0.0_WP, f_lambda_l=0.0_WP, f_Lv=0.0_WP, f_C_l=0.0_WP, f_Pv=0.0_WP, f_sigma=0.0_WP
     real(WP) :: f_rho_v=0.0_WP, f_mu_v=0.0_WP, f_lambda_v=0.0_WP, f_Cp_v=0.0_WP, f_G_v=0.0_WP

     ! Fixed Non-dimensional numbers
     real(WP), pointer :: fixed_Re, fixed_We, fixed_DRa, fixed_DRv, fixed_VRa, fixed_VRv, fixed_De

     ! MPI related
     integer :: nranks=1, irank=0, comm

     ! End of simulation
     logical :: end=.false.

  end type spray_t  

  ! Spray droplet distribution function
  integer, parameter :: &
       type_delta                     = 1, &
       type_rosin_rammler             = 2, &
       type_log_normal                = 3, &
       type_gamma                     = 4, &
       type_inverse_gamma             = 5, &
       type_general_gamma             = 6, &
       type_maximum_entropy_formalism = 7

contains

  ! Allocate spray with default values
  subroutine allocate_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    !call nullify_spray(spray)

    allocate(spray)

    allocate(spray%inp_fname); spray%inp_fname = 'noname'

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
    allocate(spray%num_noz); spray%num_noz = -9999
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
    allocate(spray%lambda_l); spray%lambda_l = -9999.0_WP
    allocate(spray%MW_f); spray%MW_f = -9999.0_WP
    allocate(spray%MP); spray%MP = -9999.0_WP
    allocate(spray%NBP); spray%NBP = -9999.0_WP
    allocate(spray%stoic_coeff); spray%stoic_coeff = -9999.0_WP

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
    allocate(spray%Anoz); spray%Anoz = -9999.0_WP
    allocate(spray%U_inj); spray%U_inj = -9999.0_WP
    allocate(spray%inj_mass); spray%inj_mass = -9999.0_WP

    allocate(spray%K_in); spray%K_in = -9999.0_WP
    allocate(spray%K_exp); spray%K_exp = -9999.0_WP
    allocate(spray%Cc0); spray%Cc0 = -9999.0_WP
    allocate(spray%Cc); spray%Cc = -9999.0_WP
    allocate(spray%D_eff); spray%D_eff = -9999.0_WP
    allocate(spray%const_inj_vel); spray%const_inj_vel = -9999.0_WP

    allocate(spray%init_dsd_name); spray%init_dsd_name = 'noname'
    allocate(spray%nd); spray%nd = -9999
    allocate(spray%dgf); spray%dgf = -9999.0_WP
    allocate(spray%init_dm); spray%init_dm = -9999.0_WP
    allocate(spray%init_d2); spray%init_d2 = -9999.0_WP
    allocate(spray%init_d3); spray%init_d3 = -9999.0_WP
    allocate(spray%init_dvar); spray%init_dvar = -9999.0_WP

    allocate(spray%skip_turb); spray%skip_turb = -9999
    allocate(spray%skip_d2); spray%skip_d2 = -9999
    allocate(spray%skip_d3); spray%skip_d3 = -9999    

    allocate(spray%R_gas); spray%R_gas = -9999.0_WP

    allocate(spray%c_k); spray%c_k = -9999.0_WP
    allocate(spray%c_mu); spray%c_mu = -9999.0_WP
    allocate(spray%c_eps1); spray%c_eps1 = -9999.0_WP
    allocate(spray%c_eps2); spray%c_eps2 = -9999.0_WP
    allocate(spray%c_zvar); spray%c_zvar = -9999.0_WP

    allocate(spray%combustion_model); spray%combustion_model = 'noname'

    allocate(spray%solver)

    allocate(spray%solver%rktvd)

    allocate(spray%solver%nr)
    ! Default initialization
    spray%solver%nr%tol = 1.0E-02_WP
    spray%solver%nr%relax_coeff = 2.0E-02_WP
    spray%solver%nr%alpha = 1.0E-04_WP
    spray%solver%nr%max_count = 2000

    allocate(spray%outfreq); spray%outfreq = -9999
    allocate(spray%datafreq); spray%datafreq = -9999

    allocate(spray%fixed_Re); spray%fixed_Re = -9999_WP
    allocate(spray%fixed_We); spray%fixed_We = -9999_WP
    allocate(spray%fixed_DRa); spray%fixed_DRa = -9999_WP
    allocate(spray%fixed_DRv); spray%fixed_DRv = -9999_WP
    allocate(spray%fixed_VRa); spray%fixed_VRa = -9999_WP
    allocate(spray%fixed_VRv); spray%fixed_VRv = -9999_WP
    allocate(spray%fixed_De); spray%fixed_De = -9999_WP

  end subroutine allocate_spray

  subroutine allocate_spray_grid_vars(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l=>null(), pc_v=>null()
    integer :: k

    spray%dz = spray%Lz/(spray%Nz-1)
    spray%nzo = spray%Nz-1+6
    spray%kmino = 1; spray%kmaxo = spray%nzo
    spray%kmin  = 4; spray%kmax  = spray%Nz+3

    allocate(spray%pc_l(spray%nzo))
    allocate(spray%pc_v(spray%nzo))

    do k = spray%kmino,spray%kmaxo
       pc_l => spray%pc_l(k)
       pc_v => spray%pc_v(k)

       call allocate_pc(pc_l)
       call allocate_pc(pc_v)

       spray%pc_l(k) = pc_l
       spray%pc_v(k) = pc_v

       deallocate(pc_l,pc_v)

    end do
    
    allocate(spray%CFL_conv(spray%nzo)); spray%CFL_conv = -9999.0_WP
    allocate(spray%CFL_bre(spray%nzo)); spray%CFL_bre = -9999.0_WP
    allocate(spray%CFL_evap(spray%nzo)); spray%CFL_evap = -9999.0_WP

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
    allocate(spray%T_sat(spray%nzo)); spray%T_sat = -9999.0_WP
    allocate(spray%L_f(spray%nzo)); spray%L_f = -9999.0_WP
    allocate(spray%rho_l_loc(spray%nzo)); spray%rho_l_loc = -9999.0_WP
    allocate(spray%sigma_loc(spray%nzo)); spray%sigma_loc = -9999.0_WP
    allocate(spray%visc_l_loc(spray%nzo)); spray%visc_l_loc = -9999.0_WP
    allocate(spray%lambda_l_loc(spray%nzo)); spray%lambda_l_loc = -9999.0_WP

    allocate(spray%T_ref(spray%nzo)); spray%T_ref = -9999.0_WP
    allocate(spray%Y_ref(spray%nzo)); spray%Y_ref = -9999.0_WP

    allocate(spray%DRg(spray%nzo)); spray%DRg = -9999.0_WP
    allocate(spray%VRg(spray%nzo)); spray%VRg = -9999.0_WP
    allocate(spray%VRtg(spray%nzo)); spray%VRtg = -9999.0_WP
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
    allocate(spray%dvar(spray%nzo)); spray%dvar = -9999.0_WP
    allocate(spray%dsd_type(spray%nzo)); spray%dsd_type = -9999
    allocate(spray%Td(spray%nzo)); spray%Td = -9999.0_WP
    allocate(spray%Tg(spray%nzo)); spray%Tg = -9999.0_WP
    allocate(spray%k_g(spray%nzo)); spray%k_g = -9999.0_WP
    allocate(spray%eps_g(spray%nzo)); spray%eps_g = -9999.0_WP
    allocate(spray%mu_t_g(spray%nzo)); spray%mu_t_g = -9999.0_WP
    allocate(spray%zmix_g(spray%nzo)); spray%zmix_g = -9999.0_WP
    allocate(spray%zvar_g(spray%nzo)); spray%zvar_g = -9999.0_WP
    allocate(spray%chi_g(spray%nzo)); spray%chi_g = -9999.0_WP
    allocate(spray%chi_g_stl(spray%nzo)); spray%chi_g_stl = -9999.0_WP
    allocate(spray%b(spray%nzo)); spray%b = -9999.0_WP

    allocate(spray%di(spray%nd,spray%nzo)); spray%di = -9999.0_WP
    allocate(spray%dsd(spray%nd,spray%nzo)); spray%dsd = -9999.0_WP
    allocate(spray%h(spray%nzo)); spray%h = -9999.0_WP
    allocate(spray%CD(spray%nd,spray%nzo)); spray%CD = -9999.0_WP
    allocate(spray%Red(spray%nd,spray%nzo)); spray%Red = -9999.0_WP
    allocate(spray%Shd(spray%nd,spray%nzo)); spray%Shd = -9999.0_WP
    allocate(spray%Nud(spray%nd,spray%nzo)); spray%Nud = -9999.0_WP
 
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
    allocate(spray%omega_k_g_p(spray%nzo)); spray%omega_k_g_p = -9999.0_WP
    allocate(spray%omega_k_g_d(spray%nzo)); spray%omega_k_g_d = -9999.0_WP
    allocate(spray%omega_eps_g_p(spray%nzo)); spray%omega_eps_g_p = -9999.0_WP
    allocate(spray%omega_eps_g_d(spray%nzo)); spray%omega_eps_g_d = -9999.0_WP
    allocate(spray%omega_zvar_g_p(spray%nzo)); spray%omega_zvar_g_p = -9999.0_WP
    allocate(spray%omega_zvar_g_d(spray%nzo)); spray%omega_zvar_g_d = -9999.0_WP

  end subroutine allocate_spray_grid_vars

  subroutine deallocate_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    deallocate(spray%inp_fname)

    deallocate(spray%Nz,spray%nzo,spray%step,spray%kmin,spray%kmax,spray%kmino,spray%kmaxo)

    deallocate(spray%z,spray%Lz,spray%dz)

    deallocate(spray%dt,spray%ndtime,spray%final_time,spray%ndftime,spray%tau)

    deallocate(spray%MaxCFL,spray%CFL,spray%CFL_conv,spray%CFL_bre,spray%CFL_evap)

    deallocate(spray%Re,spray%We,spray%DRa,spray%DRv,spray%DRg,spray%VRa,spray%VRv,spray%VRg)
    
    deallocate(spray%VRtg)
    
    deallocate(spray%WR,spray%De,spray%LR,spray%CR)

    deallocate(spray%noz_D,spray%noz_LD,spray%noz_rD,spray%noz_Dsac, spray%num_noz)
    
    deallocate(spray%Cnoz,spray%C_theta,spray%theta,spray%beta,spray%spray_angle_model)

    deallocate(spray%rho,spray%Y_l,spray%Y_v,spray%Y_a, spray%Y_g, &
               spray%u_l,spray%u_g, &
               spray%d3,spray%d2,spray%dm,spray%dvar, &
               spray%Td,spray%b,spray%Tg, &
               spray%k_g,spray%eps_g,spray%mu_t_g, &
               spray%zmix_g, spray%zvar_g, spray%chi_g, spray%chi_g_stl)

    deallocate(spray%Fuel,spray%T_fuel,spray%sigma,spray%rho_l,spray%visc_l, spray%lambda_l,spray%C_l, &
               spray%rho_l_loc, spray%sigma_loc, spray%visc_l_loc, spray%lambda_l_loc, &
               spray%p_vap,spray%T_sat, spray%MW_f,spray%L_f,spray%MP,spray%NBP,spray%stoic_coeff)

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

    deallocate(spray%P_inj, spray%Anoz, spray%C_d, spray%U_inj, spray%inj_mass)

    deallocate(spray%K_in, spray%K_exp, spray%Cc0, spray%Cc, spray%D_eff)

    deallocate(spray%init_dsd_name)

    deallocate(spray%nd)

    deallocate(spray%init_dm, spray%init_d2, spray%init_d3, spray%init_dvar, spray%dgf)

    deallocate(spray%skip_turb, spray%skip_d2, spray%skip_d3)

    deallocate(spray%dsd, spray%dsd_type, spray%CD, spray%Red, spray%Shd, spray%Nud, spray%dsdlam)

    deallocate(spray%di,spray%h)

    deallocate(spray%R_gas)

    deallocate(spray%c_k, spray%c_mu, spray%c_eps1, spray%c_eps2, spray%c_zvar)

#ifdef MDUC_MPI
    deallocate(spray%combustion_model, spray%zz, spray%bpdf)
#endif

    call deallocate_solver(spray%solver)

    deallocate(spray%omega_ent, spray%omega_vap, spray%omega_vapdm, spray%omega_vapd2, spray%omega_vapd3, &
               spray%f_drag, spray%omega_bre1, spray%omega_bre2, spray%omega_bre3, spray%omega_T, &
               spray%omega_k_g_p, spray%omega_k_g_d, spray%omega_eps_g_p, spray%omega_eps_g_d, &
               spray%omega_zvar_g_p, spray%omega_zvar_g_d)

    deallocate(spray%outfreq, spray%datafreq)

    deallocate(spray%time, spray%LPL, spray%VPL, spray%chi_st)

    deallocate(spray%fixed_Re, spray%fixed_We,  spray%fixed_DRa, spray%fixed_DRv, spray%fixed_VRa, spray%fixed_VRv, spray%fixed_De)

    deallocate(spray)

  end subroutine deallocate_spray
end module spray_defs
