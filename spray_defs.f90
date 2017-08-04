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
     real(WP), pointer :: noz_D, noz_LD, noz_rD, noz_DoDi, Cnoz, theta, beta
     
     ! Flow variables
     real(WP), dimension(:), pointer :: z, rho, Y_l, Y_v, Y_a, Y_g, u_l, u_g, d2, dm, Td, b, Tg

     ! Source terms
     real(WP), dimension(:), pointer :: omega_ent, omega_vap, omega_vapdm, omega_vapd2,f_drag, omega_bre1, omega_bre2, omega_T

     ! Fuel name
     character(len=128), pointer :: Fuel

     ! Pure chemical object for liquid phase
     type(pc_t), dimension(:), pointer :: pc_l
     
     ! Pure chemical object for vapor phase
     type(pc_t), dimension(:), pointer :: pc_v
 
     ! Fuel properties
     real(WP), pointer :: T_fuel, sigma, rho_l, visc_l, MW_f
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

     ! Nozzle flow constants/variables
     real(WP), pointer :: K_in, K_exp, Cc0, Cc, D_eff

     ! Droplet size distribution
     character(len=128), pointer :: init_dsd_name
     integer, pointer :: nd
     real(WP), pointer :: h, init_dm, init_d2
     real(WP), dimension(:,:), pointer :: dsd, CD, Red, Shd, Nud
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

contains

  subroutine nullify_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    nullify(spray%Nz,spray%nzo,spray%step,spray%kmin,spray%kmax,spray%kmino,spray%kmaxo)

    nullify(spray%z,spray%Lz,spray%dz)

    nullify(spray%dt,spray%ndtime,spray%final_time,spray%ndftime,spray%tau)

    nullify(spray%MaxCFL,spray%CFL)

    nullify(spray%Re,spray%We,spray%DRa,spray%DRv,spray%DRg,spray%VRa,spray%VRv,spray%VRg)
    
    nullify(spray%WR,spray%De,spray%LR,spray%CR)

    nullify(spray%noz_D,spray%noz_LD,spray%noz_rD,spray%noz_DoDi,spray%Cnoz,spray%theta,spray%beta)

    nullify(spray%rho,spray%Y_l,spray%Y_v,spray%Y_a, spray%Y_g,spray%u_l,spray%u_g,spray%d2,spray%dm,spray%Td,spray%b,spray%Tg)

    nullify(spray%Fuel,spray%T_fuel,spray%sigma,spray%rho_l,spray%visc_l,spray%C_l,spray%p_vap,spray%MW_f,spray%L_f)

    nullify(spray%LFPTname,spray%LFPT,spray%VFPTname,spray%VFPT)

    nullify(spray%pc_l)

    nullify(spray%pc_v)
    
    nullify(spray%rho_v, spray%visc_v, spray%lambda_v, spray%Cp_v)

    nullify(spray%rho_rv, spray%visc_rv, spray%lambda_rv, spray%Cp_rv, spray%G_rv)
    
    nullify(spray%MW_a, spray%Z_a, spray%P_a, spray%T_a, spray%rho_a, spray%visc_a, spray%lambda_a, spray%Cp_a)

    nullify(spray%rho_ra, spray%visc_ra, spray%lambda_ra, spray%Cp_ra)

    nullify(spray%rho_g, spray%visc_g, spray%lambda_g, spray%Cp_g, spray%Sc_g, spray%Pr_g)

    nullify(spray%T_ref, spray%Y_ref)

    nullify(spray%roi_file, spray%roi)

    nullify(spray%P_inj, spray%C_d, spray%U_inj)

    nullify(spray%K_in, spray%K_exp, spray%Cc0, spray%Cc, spray%D_eff)

    nullify(spray%init_dsd_name)

    nullify(spray%nd)

    nullify(spray%h, spray%init_dm, spray%init_d2)

    nullify(spray%dsd, spray%CD, spray%Red, spray%Shd, spray%Nud)

    nullify(spray%di)

    nullify(spray%R_gas)

    nullify(spray%solver)

    nullify(spray%omega_ent, spray%omega_vap, spray%omega_vapdm, spray%omega_vapd2, spray%f_drag, spray%omega_bre1, spray%omega_bre2, spray%omega_T)

    nullify(spray%outfreq)

    nullify(spray)

  end subroutine nullify_spray

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
    allocate(spray%noz_DoDi); spray%noz_DoDi = -9999.0_WP
    allocate(spray%Cnoz); spray%Cnoz = -9999.0_WP
    allocate(spray%theta); spray%theta = -9999.0_WP
    allocate(spray%beta); spray%beta = -9999.0_WP

    allocate(spray%Fuel); spray%Fuel = 'noname'
    allocate(spray%LFPTname); spray%LFPTname = 'noname'
    allocate(spray%VFPTname); spray%VFPTname = 'noname'

    allocate(spray%T_fuel); spray%T_fuel = -9999.0_WP
    allocate(spray%sigma); spray%sigma = -9999.0_WP
    allocate(spray%rho_l); spray%rho_l = -9999.0_WP
    allocate(spray%visc_l); spray%visc_l = -9999.0_WP
    allocate(spray%MW_f); spray%MW_f = -9999.0_WP

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

    allocate(spray%init_dsd_name); spray%init_dsd_name = 'noname'
    allocate(spray%nd); spray%nd = -9999
    allocate(spray%h); spray%h = -9999.0_WP
    allocate(spray%init_dm); spray%init_dm = -9999.0_WP
    allocate(spray%init_d2); spray%init_d2 = -9999.0_WP

    allocate(spray%R_gas); spray%R_gas = -9999.0_WP

    allocate(spray%solver)

    allocate(spray%solver%rk)

    allocate(spray%outfreq); spray%outfreq = -9999

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

    call allocate_solver(spray%solver,spray%nzo)

    allocate(spray%omega_ent(spray%nzo)); spray%omega_ent = -9999.0_WP
    allocate(spray%omega_vap(spray%nzo)); spray%omega_vap = -9999.0_WP
    allocate(spray%omega_vapdm(spray%nzo)); spray%omega_vapdm = -9999.0_WP
    allocate(spray%omega_vapd2(spray%nzo)); spray%omega_vapd2 = -9999.0_WP
    allocate(spray%f_drag(spray%nzo)); spray%f_drag = -9999.0_WP
    allocate(spray%omega_bre1(spray%nzo)); spray%omega_bre1 = -9999.0_WP
    allocate(spray%omega_bre2(spray%nzo)); spray%omega_bre2 = -9999.0_WP
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

    deallocate(spray%noz_D,spray%noz_LD,spray%noz_rD,spray%noz_DoDi,spray%Cnoz,spray%theta,spray%beta)

    deallocate(spray%rho,spray%Y_l,spray%Y_v,spray%Y_a, spray%Y_g,spray%u_l,spray%u_g,spray%d2,spray%dm,spray%Td,spray%b,spray%Tg)

    deallocate(spray%Fuel,spray%T_fuel,spray%sigma,spray%rho_l,spray%visc_l,spray%C_l,spray%p_vap,spray%MW_f,spray%L_f)

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

    deallocate(spray%h, spray%init_dm, spray%init_d2)

    deallocate(spray%dsd, spray%CD, spray%Red, spray%Shd, spray%Nud)

    deallocate(spray%di)

    deallocate(spray%R_gas)

    call deallocate_solver(spray%solver)

    deallocate(spray%omega_ent, spray%omega_vap, spray%omega_vapdm, spray%omega_vapd2, spray%f_drag, spray%omega_bre1, spray%omega_bre2, spray%omega_T)

    deallocate(spray%outfreq)

    deallocate(spray%time, spray%LPL, spray%VPL)

    deallocate(spray%fixed_Re, spray%fixed_We,  spray%fixed_DRa, spray%fixed_DRv, spray%fixed_VRa, spray%fixed_VRv, spray%fixed_De)

    deallocate(spray)

  end subroutine deallocate_spray
end module spray_defs
