!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com or a.saha at itv.rwth-aachen.de           !
!--------------------------------------------------------------------------------!

module spray_defs
  use precision
  use string
  use pc_defs
  use fpt_defs
  use solver_defs

  implicit none

  ! Definition of spray object
  type spray_t
     
     ! Input
     character(len=128), pointer :: inp_fnamea

     ! Number of grid points
     integer, pointer :: Nz, nzo, step, nro
     integer, pointer :: kmin, kmax, kmino, kmaxo

     ! Numerics
     real(WP), pointer :: dz, Lz
     real(WP), pointer :: dt, ndtime, final_time, ndftime, tau
     real(WP), pointer :: MaxCFL, CFL
     real(WP), dimension(:), pointer :: CFL_conv, CFL_bre, CFL_evap

     ! Number of equations and scalars
     integer :: neq, nsc

     ! Non-dimensional numbers
     real(WP), pointer :: Re, We, WR, De, DRa, DRv, VRa, VRv, Pr_l
     real(WP), dimension(:), pointer :: DRl, DRg, VRl, VRg, VRtg, TCRg, LR, CR, SR
     
     ! Geometric variables
     real(WP), pointer :: noz_D, noz_LD, noz_rD, noz_Dsac
     ! Variables to avoid extremely small time-step sizes

     integer, pointer :: num_noz     

     ! Flow variables
     real(WP), dimension(:), pointer :: z, rho, Y_l, Y_v, Y_a, Y_g, u_l, u_g, d2, dm, d3, dvar, Td, b, Tg, Tv, &
                                        k_g, eps_g, mu_t_g, zmix_g, zvar_g, chi_g, chi_g_stl, b_old

     ! Source terms
     real(WP), dimension(:), pointer :: omega_ent, omega_vap, omega_vapdm, omega_vapd2, omega_vapd3, omega_vapres, &
                                        omega_exp, omega_expdm, omega_expd2, omega_expd3, &
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
      real(WP), pointer :: P_bub, L_fb, L_fbub, C_lb, T_fuel, T_boil, T_crit, sigma,  sigmab, sigma_bub, rho_l, rho_lb, &
                                        rho_lbub, visc_l, visc_lb, lambda_l, lambda_lb, MW_f, MVol, MP, NBP, stoic_coeff
      real(WP), dimension(:), pointer :: L_f, C_l, p_vap, T_sat, sigma_loc, rho_l_loc, visc_l_loc, lambda_l_loc

     ! Analytical heat transfer correction factor
     real(WP), dimension(:,:), pointer :: f2

     ! Fuel properties from table
     character(len=128), pointer :: FPTname
     type(fpt_t), pointer :: fpt
     real(WP), dimension(:,:), pointer :: LFPT, VFPT
     
     ! Fuel vapor properties
     real(WP), pointer :: rho_v, visc_v, lambda_v, Cp_v

     ! Ambient properties
     real(WP), pointer :: MW_a, Z_a, P_a, T_a, rho_a, visc_a, lambda_a, Cp_a

     ! Gas phase mixture properties
     real(WP), dimension(:), pointer :: rho_g, visc_g, lambda_g, Cp_g, Sc_g, Pr_g

     ! Reference variables for evaporation model
     real(WP), dimension(:), pointer :: T_ref, Y_ref, rho_rv, visc_rv, lambda_rv, Cp_rv, G_rv, rho_ra, visc_ra, lambda_ra, Cp_ra

     ! Bubble parameter 
     real(WP), dimension(:), pointer :: Rdot, Nbub, db

     ! Injector/Injection parameters
     character(len=128), pointer :: roi_file
     real(WP), dimension(:,:), pointer :: roi
     real(WP), pointer :: P_inj, C_d, Anoz, U_inj, inj_mass
     logical :: noz_flow_model = .false.
     logical :: flash_boiling = .true.

     ! Spray angle parameters
     character(len=128), pointer :: spray_angle_model
     real(WP), pointer ::  Cnoz, C_theta, theta, beta
     !real(WP), dimension(:), pointer :: theta, beta

     ! Saturation P-T file
     character(len=128), pointer :: satPT_file

     ! Entrainment model (default is true)
     logical :: ent_model = .true.

     ! Drag model (default is true)
     logical :: drag_model = .true.

     ! Nozzle flow constants/variables
     real(WP), pointer :: K_in, K_exp, Cc0, Cc, D_eff, const_inj_vel

     ! Droplet size distribution
     character(len=128), pointer :: init_dsd_name
     integer, pointer :: nd
     real(WP), pointer :: init_dm, init_d2, init_d3, init_dvar, d10, d20, d30, dgf
     real(WP), dimension(:), pointer :: h
     real(WP), dimension(:,:), pointer :: di, dsd, CD, Red, Shd, Nud, dsdlam, Tdi
     integer, dimension(:), pointer :: dsd_type
     integer, pointer :: skip_d2, skip_d3

     ! Breakup model (default is true)
     logical :: bre_model = .true.
     ! Breakup parameters with default values
     real(WP) :: B0 = 0.61_WP, B1 = 10.0_WP, C3 = 2.5_WP, Crel = 1.0_WP, D0 = 1.0_WP, D1 = 1.0_WP

     ! Evaporation model (default is true)
     logical :: evap_model = .true.
     ! Evaporation parameters
     real(WP) :: Cevap = 1.0_WP
     real(WP) :: A0 != 1.0_WP
     real(WP) :: A1 != 1.0_WP zmiany do sprawdzenia KP 

     ! Universal constants
     real(WP), pointer :: R_gas

     ! Turbulence parameters
     logical :: turb_model = .false.
     integer, pointer :: skip_turb
     real(WP), pointer :: c_mu, c_eps1, c_eps2, c_zvar
     real(WP), dimension(:), pointer :: c_k
     character(len=128), pointer :: ck_type
     real(WP) :: ck_val = -9999.0_WP

     ! Combustion model
     character(len=128), pointer :: combustion_model
     real(WP) :: C_chi = 1.0_WP, Y_O2 = 0.164215_WP, Zmix_st, flm_int=1.0E-04_WP
     real(WP), dimension(:), pointer :: zz=>null(), bpdf=>null(), alpha=>null(), MassPDF=>null()
     real(WP), dimension(:), pointer :: Tf=>null(), SootNumDen=>null(), SootVolume=>null(), SootMass=>null()
     real(WP), dimension(:,:), pointer :: Sc=>null(), Wt=>null(), Yspec=>null()
     type(srif_t), pointer :: first_rif=>null()
     integer :: nflamelet = -9999
     ! 2D profiles of Zmix and Zvar
     real(WP), dimension(:,:), pointer :: Zmix2D=>null(), Zvar2D=>null(), r=>null(), rho2D=>null(), rho_g2D=>null(), gVF2D=>null()
     real(WP), dimension(:,:), pointer :: Tf2D=>null(), SootNumDen2D=>null(), SootVolume2D=>null(), SootMass2D=>null()
     real(WP), dimension(:,:,:), pointer :: Yspec2D=>null()
     real(WP), dimension(:), pointer :: xi=>null()
     logical :: twoD = .false.
     real(WP) :: TotSootMass

     ! Solver
     type(solver_t), pointer :: solver

     ! Output
     integer, pointer :: outfreq, datafreq
     character(len=128) :: datafilename, plfilename
     logical :: saveDataFile=.false.

     ! Postprocessing
     real(WP), dimension(:), pointer :: time, LPL, VPL, chi_st, chi_st1, utip, dtip
     real(WP), dimension(:), pointer :: chi_st_m

     ! For sensitivity analysis
     real(WP) :: f_rho_l=0.0_WP, f_mu_l=0.0_WP, f_lambda_l=0.0_WP, f_Lv=0.0_WP, f_C_l=0.0_WP, f_Pv=0.0_WP, f_sigma=0.0_WP
     real(WP) :: f_rho_v=0.0_WP, f_mu_v=0.0_WP, f_lambda_v=0.0_WP, f_Cp_v=0.0_WP, f_G_v=0.0_WP

     ! Fixed Non-dimensional numbers
     real(WP), pointer :: fixed_Re, fixed_We, fixed_DRa, fixed_DRv, fixed_VRa, fixed_VRv, fixed_De

     ! MPI related
     integer :: nranks=1, irank=0, comm

     ! End of simulation
     logical :: end=.false.

     character(len=256) :: inp_fname

  end type spray_t

  ! Adapted from CIAO routines by Marco Davidovic
  type mduc_t
     ! ==================================================================================== !
     !!! this derived datatype has a corresponding MPI datatype                           !!!
     !!! any changes must be consistent with create_mpi_.. !!!
     ! ==================================================================================== !
     SEQUENCE

     integer(8) :: mem          ! memory location
     !integer :: index          ! mainly to distinguish variables
     integer :: nx1,nx2         ! number of grid points in each direction
     integer :: hdim            ! enthalpy dimension in mduc model
     integer :: transf2d        ! coordinate transformation type in 2D
     ! activate(1)/deactivate(0) source terms : 
     ! Nucleation,Condensation,Coagulation,SurfaceGrowth,Oxidation,Fragmentation
     integer, dimension(:), pointer :: sootActive => null()
     ! ==================================================================================== !
     !!! this derived datatype has a corresponding MPI datatype                           !!!
     !!! any changes must be consistent with create_mpi_.. !!!
     ! ==================================================================================== !
  end type mduc_t

  type srif_t
     ! flamelet number
     integer :: i
     ! is this flamelet currently solved?
     logical :: active
     ! mduc data
     type(mduc_t), pointer :: mduc
     ! flow data
     !type(rif_flow_t), dimension(:),allocatable :: f
     ! major species fields
     real(WP), dimension(:), pointer :: Yi_SC => null() 
     ! flamelet marker (tracking scalar for MRIF)
     real(WP), dimension(:), pointer :: mark => null()
     ! flamelet weight (MRIF)
     real(WP), dimension(:), pointer :: weight => null()
     ! mask for locally disabling rif
     integer, dimension(:), pointer :: mask => null()
     ! progress variable field
     integer, dimension(:), pointer :: prog => null()
     character(len=str_medium), dimension(:), allocatable :: species,prog_species
     integer, dimension(:), allocatable   :: iprog_species
     real(WP), dimension(:), allocatable :: bpdf, Temp, sootM00, sootM10
     real(WP), dimension(:,:), allocatable :: Yspec
     real(WP), dimension(:,:), allocatable :: Temp2D, sootM002D, sootM102D
     real(WP), dimension(:,:,:), allocatable :: Yspec2D
     !real(WP) :: Temp, sootM00, sootM10
     ! flamelet dimension
     integer :: ndim
     ! number of species in chemical mechanism
     integer :: Nspec
     ! number of non-steady-state species 
     integer :: NspecS 
     ! number of major species
     integer :: NspecT
     ! soot method (local / conditional)
     integer :: soot_method 
     ! number of soot moments
     integer :: nSootEq
     ! soot species names
     character(len=str_medium), dimension(:), pointer :: soot_names => null()
     ! sott species index
     integer, dimension(:), pointer :: soot_ispec => null()
     ! phi threshold
     real(WP) :: phimin_thresh(2) 
     ! weight threshold (MRIF)
     real(WP) :: wmin_thresh
     ! pdf type (recently changed to integer values)
     !character(len=str_medium) :: pdf_type(2)
     integer :: pdf_type(2)
     ! injection status
     integer :: inj
     ! start of chemistry calculation
     real(WP), dimension(:), pointer :: tinit => null()
     ! scalar dissipation rate mode (distribution or global-zero)
     integer, dimension(:), pointer :: chi_mode => null()
     ! scalar dissipation rate functional form
     integer :: chi_form  
     !
     integer, dimension(:), pointer :: chi_norm => null()
     ! calculate reference scalar dissipation to scale functional form?
     logical :: calc_chi_ref
     ! write frequency in steps
     integer :: dump_freq
     ! root process which solves 1d flamelet
     integer :: root, nproc, irank
     ! store mduc grid
     real(WP), dimension(:), pointer :: zz=>null()
     real(WP), dimension(:), pointer :: chi=>null()
     ! flamelet solution arrays
     ! temperature 1d
     real(WP), dimension(:), pointer :: fT_1d => null(), frrho_1d => null()
     ! all species 1d
     real(WP), dimension(:,:), pointer :: fY_1d => null()
     ! major species 1d
     real(WP), dimension(:,:), pointer :: fYT_1d => null()
     ! soot source terms 1d
     real(WP), dimension(:,:), pointer :: fS_1d => null()
     ! temperature 2d
     real(WP), dimension(:,:), pointer :: fT_2d => null()
     ! major species 2d
     real(WP), dimension(:,:,:), pointer :: fY_2d => null()
     ! next flamelet (MRIF)
     type(srif_t), pointer :: next => null()
  end type srif_t

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

    !allocate(character(len=64) :: spray%inp_fname); spray%inp_fname = 'noname'

    allocate(spray%Nz); spray%Nz = -9999
    allocate(spray%nzo); spray%nzo = -9999
    allocate(spray%nro); spray%nro = -9999
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

    allocate(spray%FPTname); spray%FPTname = 'noname'
    allocate(spray%T_fuel); spray%T_fuel = -9999.0_WP
    allocate(spray%T_boil); spray%T_boil = -9999.0_WP
    allocate(spray%T_crit); spray%T_crit = -9999.0_WP
    allocate(spray%P_bub); spray%P_bub = -9999.0_WP
    allocate(spray%sigma); spray%sigma = -9999.0_WP
    allocate(spray%sigmab); spray%sigmab = -9999.0_WP
    allocate(spray%sigma_bub); spray%sigma_bub = -9999.0_WP
    allocate(spray%rho_l); spray%rho_l = -9999.0_WP
    allocate(spray%rho_lb); spray%rho_lb = -9999.0_WP
    allocate(spray%rho_lbub); spray%rho_lbub = -9999.0_WP
    allocate(spray%visc_l); spray%visc_l = -9999.0_WP
    allocate(spray%visc_lb); spray%visc_lb = -9999.0_WP
    allocate(spray%lambda_l); spray%lambda_l = -9999.0_WP
    allocate(spray%lambda_lb); spray%lambda_lb = -9999.0_WP
    allocate(spray%C_lb); spray%C_lb = -9999.0_WP
    allocate(spray%L_fb); spray%L_fb = -9999.0_WP
    allocate(spray%L_fbub); spray%L_fbub = -9999.0_WP
    allocate(spray%MW_f); spray%MW_f = -9999.0_WP
    allocate(spray%MVol); spray%MVol = -9999.0_WP
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
    allocate(spray%satPT_file); spray%satPT_file = 'noname'
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
    allocate(spray%ck_type); spray%ck_type = 'Dynamic'

    allocate(spray%R_gas); spray%R_gas = -9999.0_WP

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

    !allocate(spray%theta(spray%nzo)); spray%theta = 0.0_WP
    !allocate(spray%beta(spray%nzo)); spray%beta = 0.0_WP

    allocate(spray%DRl(spray%nzo)); spray%DRl = -9999.0_WP
    allocate(spray%DRg(spray%nzo)); spray%DRg = -9999.0_WP
    allocate(spray%VRl(spray%nzo)); spray%VRl = -9999.0_WP
    allocate(spray%VRg(spray%nzo)); spray%VRg = -9999.0_WP
    allocate(spray%VRtg(spray%nzo)); spray%VRtg = -9999.0_WP
    allocate(spray%LR(spray%nzo)); spray%LR = -9999.0_WP
    allocate(spray%CR(spray%nzo)); spray%CR = -9999.0_WP
    allocate(spray%SR(spray%nzo)); spray%SR = -9999.0_WP

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
    allocate(spray%Tv(spray%nzo)); spray%Tv = 0.0_WP
    allocate(spray%c_k(spray%nzo)); spray%c_k = -9999.0_WP
    allocate(spray%k_g(spray%nzo)); spray%k_g = -9999.0_WP
    allocate(spray%eps_g(spray%nzo)); spray%eps_g = -9999.0_WP
    allocate(spray%mu_t_g(spray%nzo)); spray%mu_t_g = -9999.0_WP
    allocate(spray%zmix_g(spray%nzo)); spray%zmix_g = -9999.0_WP
    allocate(spray%zvar_g(spray%nzo)); spray%zvar_g = -9999.0_WP
    allocate(spray%chi_g(spray%nzo)); spray%chi_g = -9999.0_WP
    allocate(spray%chi_g_stl(spray%nzo)); spray%chi_g_stl = -9999.0_WP
    allocate(spray%b(spray%nzo)); spray%b = -9999.0_WP
    allocate(spray%b_old(spray%nzo)); spray%b_old = -9999.0_WP
    allocate(spray%alpha(spray%nzo)); spray%alpha = -9999.0_WP

    allocate(spray%di(spray%nd,spray%nzo)); spray%di = -9999.0_WP
    allocate(spray%f2(spray%nd,spray%nzo)); spray%f2 = 1.0_WP
    !allocate(spray%Tdi(spray%nd,spray%nzo)); spray%Tdi = -9999.0_WP
    allocate(spray%dsd(spray%nd,spray%nzo)); spray%dsd = -9999.0_WP
    allocate(spray%h(spray%nzo)); spray%h = -9999.0_WP
    allocate(spray%CD(spray%nd,spray%nzo)); spray%CD = -9999.0_WP
    allocate(spray%Red(spray%nd,spray%nzo)); spray%Red = -9999.0_WP
    allocate(spray%Shd(spray%nd,spray%nzo)); spray%Shd = -9999.0_WP
    allocate(spray%Nud(spray%nd,spray%nzo)); spray%Nud = -9999.0_WP
 
    allocate(spray%dsdlam(4,spray%nzo)); spray%dsdlam = -9999.0_WP
    
    allocate(spray%Rdot(spray%nd)); spray%Rdot = 0.0_WP
    allocate(spray%Nbub(spray%nd)); spray%Nbub = 0.0_WP
    allocate(spray%db(spray%nd)); spray%db = 0.0_WP

    if(spray%twoD) then
       allocate(spray%Zmix2D(spray%nzo,spray%nro+1)); spray%Zmix2D = -9999.0_WP
       allocate(spray%Zvar2D(spray%nzo,spray%nro+1)); spray%Zvar2D = -9999.0_WP
       allocate(spray%r(spray%nzo,spray%nro+1)); spray%r = -9999.0_WP
       allocate(spray%rho2D(spray%nzo,spray%nro+1)); spray%rho2D = -9999.0_WP
       allocate(spray%rho_g2D(spray%nzo,spray%nro+1)); spray%rho_g2D = -9999.0_WP
       allocate(spray%gVF2D(spray%nzo,spray%nro+1)); spray%gVF2D = -9999.0_WP
       allocate(spray%xi(spray%nro+1)); spray%xi = -9999.0_WP
    end if

    call allocate_solver(spray%solver,spray%nzo)

    allocate(spray%omega_ent(spray%nzo)); spray%omega_ent = -9999.0_WP
    allocate(spray%omega_vap(spray%nzo)); spray%omega_vap = -9999.0_WP
    allocate(spray%omega_vapdm(spray%nzo)); spray%omega_vapdm = -9999.0_WP
    allocate(spray%omega_vapd2(spray%nzo)); spray%omega_vapd2 = -9999.0_WP
    allocate(spray%omega_vapd3(spray%nzo)); spray%omega_vapd3 = -9999.0_WP
    allocate(spray%omega_vapres(spray%nzo)); spray%omega_vapres = -9999.0_WP

    allocate(spray%omega_exp(spray%nzo)); spray%omega_exp = -9999.0_WP
    allocate(spray%omega_expdm(spray%nzo)); spray%omega_expdm = -9999.0_WP
    allocate(spray%omega_expd2(spray%nzo)); spray%omega_expd2 = -9999.0_WP
    allocate(spray%omega_expd3(spray%nzo)); spray%omega_expd3 = -9999.0_WP

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

  subroutine list_destroy( list )
    type(srif_t), pointer  :: list

    type(srif_t), pointer  :: current
    type(srif_t), pointer  :: next

    current => list
    do while ( associated(current%next) )
       next => current%next
       deallocate( current%mduc)
       deallocate( current%Yspec )
       deallocate( current%Yspec2D )
       deallocate( current%species )
       deallocate( current%bpdf )
       deallocate( current%Temp )
       deallocate( current%Temp2D )
       if(associated(current%fS_1d)) deallocate( current%sootM00 )
       if(associated(current%fS_1d)) deallocate( current%sootM002D )
       if(associated(current%fS_1d)) deallocate( current%sootM10 )
       if(associated(current%fS_1d)) deallocate( current%sootM102D )
       if(associated(current%fS_1d)) deallocate( current%frrho_1d )
       deallocate( current%prog_species )
       deallocate( current%iprog_species )
       deallocate( current%fY_1d )
       deallocate( current%fT_1d )
       deallocate( current%chi )
       if(associated(current%fS_1d)) deallocate( current%fS_1d )
       deallocate( current )
       current => next
    enddo
  end subroutine list_destroy

  subroutine deallocate_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    !deallocate(spray%inp_fname)

    deallocate(spray%Nz,spray%nzo,spray%nro,spray%step,spray%kmin,spray%kmax,spray%kmino,spray%kmaxo)

    deallocate(spray%z,spray%Lz,spray%dz)

    deallocate(spray%dt,spray%ndtime,spray%final_time,spray%ndftime,spray%tau)

    deallocate(spray%MaxCFL,spray%CFL,spray%CFL_conv,spray%CFL_bre,spray%CFL_evap)

    deallocate(spray%Re,spray%We,spray%DRa,spray%DRv,spray%DRg,spray%VRa,spray%VRv,spray%VRg)
    
    deallocate(spray%VRtg)
    
    deallocate(spray%WR,spray%De,spray%LR,spray%CR,spray%SR,spray%DRl,spray%VRl)

    deallocate(spray%noz_D,spray%noz_LD,spray%noz_rD,spray%noz_Dsac, spray%num_noz)
    
    deallocate(spray%Cnoz,spray%C_theta,spray%theta,spray%beta,spray%spray_angle_model)

    deallocate(spray%rho,spray%Y_l,spray%Y_v,spray%Y_a, spray%Y_g, &
               spray%u_l,spray%u_g, &
               spray%d3,spray%d2,spray%dm,spray%dvar, &
               spray%Td,spray%b, spray%b_old, spray%Tg, spray%Tv, &
               spray%k_g,spray%eps_g,spray%mu_t_g, &
               spray%zmix_g, spray%zvar_g, spray%chi_g, spray%chi_g_stl, &
               spray%alpha)

    deallocate(spray%Fuel,spray%T_fuel,spray%sigma,spray%rho_l,spray%visc_l,spray%lambda_l,spray%C_l, &
               spray%T_boil,spray%T_crit,spray%sigmab,spray%sigma_bub,spray%rho_lb,spray%rho_lbub, &
	       spray%visc_lb,spray%lambda_lb,spray%C_lb,spray%L_fb, &
               spray%rho_l_loc, spray%sigma_loc, spray%visc_l_loc, spray%lambda_l_loc, &
               spray%P_bub,spray%p_vap,spray%T_sat,spray%MW_f,spray%MVol,spray%L_f,spray%MP,spray%NBP,spray%stoic_coeff)

    deallocate(spray%Rdot, spray%Nbub, spray%db)

    deallocate(spray%f2)

    deallocate(spray%pc_l)

    deallocate(spray%pc_v)

    deallocate(spray%FPTname)
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

    deallocate(spray%satPT_file, spray%roi_file, spray%roi)

    deallocate(spray%P_inj, spray%Anoz, spray%C_d, spray%U_inj, spray%inj_mass)

    deallocate(spray%K_in, spray%K_exp, spray%Cc0, spray%Cc, spray%D_eff)

    deallocate(spray%init_dsd_name)

    deallocate(spray%nd)

    deallocate(spray%init_dm, spray%init_d2, spray%init_d3, spray%init_dvar, spray%dgf)

    deallocate(spray%skip_turb, spray%skip_d2, spray%skip_d3, spray%ck_type)

    deallocate(spray%dsd, spray%dsd_type, spray%CD, spray%Red, spray%Shd, spray%Nud, spray%dsdlam)

    deallocate(spray%di,spray%h)!,spray%Tdi)

    deallocate(spray%R_gas)

    deallocate(spray%c_k, spray%c_mu, spray%c_eps1, spray%c_eps2, spray%c_zvar)

#ifdef MDUC_MPI
    deallocate(spray%combustion_model, spray%zz, spray%bpdf, spray%MassPDF)
    deallocate(spray%Tf, spray%SootNumDen, spray%SootVolume, spray%SootMass, spray%Yspec)

    if(spray%twoD) then
       deallocate(spray%Tf2D, spray%SootNumDen2D, spray%SootVolume2D, spray%SootMass2D, spray%Yspec2D)
    end if

    call list_destroy(spray%first_rif)
#endif

    if(spray%twoD) then
       deallocate(spray%Zmix2D, spray%Zvar2D, spray%xi, spray%r, spray%rho2D, spray%rho_g2D, spray%gVF2D)
    end if

    if (associated(spray%Sc)) then
       deallocate(spray%Sc)
    end if

    if (associated(spray%Wt)) then
       deallocate(spray%Wt)
    end if

    call deallocate_solver(spray%solver)

    deallocate(spray%omega_ent, spray%omega_vap, spray%omega_vapdm, spray%omega_vapd2, spray%omega_vapd3, spray%omega_vapres, &
            spray%f_drag, spray%omega_bre1, spray%omega_bre2, spray%omega_bre3, &
	    spray%omega_exp, spray%omega_expdm, spray%omega_expd2, spray%omega_expd3, &
            spray%omega_T, spray%omega_k_g_p, spray%omega_k_g_d, spray%omega_eps_g_p, spray%omega_eps_g_d, &
            spray%omega_zvar_g_p, spray%omega_zvar_g_d)


    deallocate(spray%outfreq, spray%datafreq)

    deallocate(spray%time, spray%LPL, spray%VPL, spray%chi_st, spray%chi_st_m, spray%chi_st1, spray%utip, spray%dtip)

    deallocate(spray%fixed_Re, spray%fixed_We,  spray%fixed_DRa, spray%fixed_DRv, spray%fixed_VRa, spray%fixed_VRv, spray%fixed_De)

    deallocate(spray)

  end subroutine deallocate_spray
end module spray_defs
