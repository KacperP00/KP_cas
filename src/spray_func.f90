module spray_func
  use precision
  use math
  use parser
  use spray_defs
  use pc_defs
  use pc_func
  use fpt_defs
  use fpt_func
  use solver_defs
  use solver_func
#ifdef MDUC_MPI
  use spray_combust
  use MPI
#endif

  implicit none

!  include 'mpif.h'

contains

  subroutine init_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(fpt_t), pointer :: fpt
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, ierr
    real(WP) :: scal, Cm, Cvar, var, nu
    real(WP), dimension(:), pointer :: T_d

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    spray%R_gas = 8.3144598_WP   ! J/K/mol

    ! Set the number of flamelets
    if(spray%combustion_model == 'MRIF') then
       if(spray%nflamelet > 0) then
          spray%nsc = spray%nflamelet
          spray%solver%nsc = spray%nsc
          write(*,*) 'User specified number of flamelets: ',spray%nsc
       else
          spray%nsc = int(spray%final_time/spray%flm_int)
          spray%solver%nsc = spray%nsc
          write(*,*) 'Setting ',spray%nsc,' flamelets...'
       end if
    end if

    ! Initialize grid and grid variables
    call allocate_spray_grid_vars(spray)

    ! Load fuel properties table if provided
    if ( spray%FPTname /= 'noname' ) then
       
       call allocate_fpt(spray%fpt)
       spray%fpt%file = spray%FPTname
       call readFuelProperties(spray%fpt)

    end if

    ! Initialize solvers
    call solver_init(spray)
    spray%neq = spray%solver%neq
    spray%nsc = spray%solver%nsc

    do k=kmino,kmaxo
       spray%z(k) = k*spray%dz - 3.5_WP*spray%dz
    end do

    ! Set sensitivity factors
    ! Default 1.0; Otherwise specified by user
    if(spray%f_rho_l==0.0_WP) spray%f_rho_l = 1.0_WP
    if(spray%f_mu_l==0.0_WP) spray%f_mu_l = 1.0_WP
    if(spray%f_lambda_l==0.0_WP) spray%f_lambda_l = 1.0_WP
    if(spray%f_Lv==0.0_WP) spray%f_Lv = 1.0_WP
    if(spray%f_C_l==0.0_WP) spray%f_C_l = 1.0_WP
    if(spray%f_Pv==0.0_WP) spray%f_Pv = 1.0_WP
    if(spray%f_sigma==0.0_WP) spray%f_sigma = 1.0_WP

    if(spray%f_rho_v==0.0_WP) spray%f_rho_v = 1.0_WP
    if(spray%f_mu_v==0.0_WP) spray%f_mu_v = 1.0_WP
    if(spray%f_lambda_v==0.0_WP) spray%f_lambda_v = 1.0_WP
    if(spray%f_Cp_v==0.0_WP) spray%f_Cp_v = 1.0_WP
    if(spray%f_G_v==0.0_WP) spray%f_G_v = 1.0_WP
    
    ! Liquid fuel properties
    if (associated(spray%FPT)) then
       call getLiquidFuelPropertiesFromLFPT(spray)
    else
       call computeLiquidFuelProperties(spray)
       call computeSatLiquidFuelProperties(spray)
    end if

    ! Ambient properties
    call computeAmbientProperties(spray)

    ! Vapor fuel properties
    if (associated(spray%FPT)) then
       call getVaporFuelPropertiesFromVFPT(spray,spray%T_a)
    else
       call computeVaporFuelProperties(spray,spray%T_a)
    end if

    if( spray%const_inj_vel .gt. 0.0_WP ) then
       spray%C_d = 1.0_WP

       spray%U_inj = spray%const_inj_vel

       spray%D_eff = spray%noz_D
    else if(spray%noz_flow_model) then
       ! Compute nozzle flow
       call nozzle_flow_model(spray)
    else if(spray%inj_mass .gt. 0.0_WP) then
       ! Compute injection parameters
       call injection_params(spray)
    else
       write(*,*) '#### Error!!! Specify one of the following to compute nozzle exit velocity....'
       write(*,*) '#### 1. Use nozzle flow model           : <.true./.false.>'
       write(*,*) '#### 2. Constant injection velocity     : <VALUE>'
       write(*,*) '#### 3. Injected mass                   : <VALUE>'
       write(*,*) '####    Area of nozzle                  : <VALUE>'
       write(*,*) '####    Discharge coefficient of nozzle : <VALUE>'
       call abort
    end if

    ! Compute non-dimensional parameters
    call compute_constNonDparams(spray)

    ! Time scale
    spray%tau = spray%D_eff/spray%U_inj

    spray%ndtime = 0.0_WP; spray%ndftime = spray%final_time/spray%tau
    
    spray%step = 0

    ! Read Rate of Injection profile if provided
    call read_ROI_from_file(spray)

    ! Compute spray half cone angle and spreading coefficient
    call compute_beta(spray)

!!$    ! Scaling parameters of initial DSD (from Spray A primary breakup simulations)
!!$    scal = 2.2278E+09_WP; Cm = 0.009363081_WP; Cvar = 2.38132E-05_WP
!!$
!!$    if(spray%init_dm == -9999.0_WP) then
!!$       spray%init_dm = min(sqrt(1.0E-01_WP),Cm*scal*spray%DRa/(spray%Re*spray%We))
!!$    end if
!!$
!!$    var = Cvar*(scal*spray%DRa/(spray%Re*spray%We))**2
!!$
!!$    if(spray%init_d2 == -9999.0_WP) then
!!$       spray%init_d2 = min(1.12E-01_WP,spray%init_dm**2 + var)
!!$    end if
!!$
!!$    spray%init_dvar = spray%init_d2 - spray%init_dm**2

!!$    if(spray%init_dm == -9999.0_WP) then
!!$       call initDropletSize(spray)
!!$    end if

    ! Skip solving turbulence equations
    if(.not. spray%turb_model) spray%skip_turb = 4
    
    ! Skip solving d2 equation
    !spray%skip_d2 = 12
    
    ! Skip solving d3 equation
    spray%skip_d3 = 13

    ! Initialize droplets
    call init_droplets(spray)
    
    write(*,*) '#####', spray%init_dm, spray%init_d2, spray%init_d3, spray%init_dvar

    ! Initialize non-dimensionalized flow variables and source terms
    spray%b = 0.5_WP !+ spray%z*spray%beta
    spray%b_old = spray%b

    !call compute_b(spray)

    spray%Y_l(kmino:kmin-1) = 1.0_WP; spray%Y_l(kmin:kmaxo) = 0.0_WP
    spray%Y_v = 0.0_WP;               
    spray%Y_a(kmino:kmin-1) = 0.0_WP; spray%Y_a(kmin:kmaxo) = 1.0_WP
    spray%Y_g = spray%Y_v + spray%Y_a

    spray%rho = 1.0_WP/(spray%DRa-spray%Y_l*(spray%DRa-1.0_WP))

    spray%u_l(kmino:kmin-1) = 0.0_WP; spray%u_l(kmin:kmaxo) = 0.0_WP

    spray%u_g = 0.0_WP

    spray%d3(kmino:kmin-1) = spray%init_d3; spray%d3(kmin:kmaxo) = 0.0_WP
    spray%d2(kmino:kmin-1) = spray%init_d2; spray%d2(kmin:kmaxo) = 0.0_WP
    spray%dm(kmino:kmin-1) = spray%init_dm; spray%dm(kmin:kmaxo) = 0.0_WP
    spray%dvar(kmino:kmin-1) = spray%init_dvar; spray%dvar(kmin:kmaxo) = 0.0_WP
 
!    spray%d3 = spray%init_d3*0.5_WP*(1.0_WP-TANH(4.0_WP*spray%z-10.0_WP))
!    spray%d2 = spray%init_d2*0.5_WP*(1.0_WP-TANH(4.0_WP*spray%z-10.0_WP))
!    spray%dm = spray%init_dm*0.5_WP*(1.0_WP-TANH(4.0_WP*spray%z-10.0_WP))

    spray%Td = 0.0_WP
    spray%Td(kmino:kmin-1) = 1.0_WP
    spray%Tdi = 0.0_WP
    spray%Tdi(:,kmino:kmin-1) = 1.0_WP
    spray%Tg = spray%T_a/spray%T_fuel
    !spray%Tg = (spray%Y_v*(1.0_WP - spray%De) + spray%Y_a*spray%T_a/spray%T_fuel)/spray%Y_g
    !spray%Tg(kmino:kmin-1) = spray%T_a/spray%T_fuel

    ! Compute stoichiometric coefficient for fuel
    if (associated(spray%FPT)) then
       spray%stoic_coeff = (spray%fpt%ChemicalFormula%C &
                         +  0.25_WP*spray%fpt%ChemicalFormula%H - 0.5_WP*spray%fpt%ChemicalFormula%O)&
                         *  32.0E-03_WP/spray%fpt%MolecularWeight
    else
       spray%stoic_coeff = (spray%pc_v(spray%kmaxo)%ChemicalFormula%C &
                         +  0.25_WP*spray%pc_v(spray%kmaxo)%ChemicalFormula%H - 0.5_WP*spray%pc_v(spray%kmaxo)%ChemicalFormula%O)&
                         *  32.0E-03_WP/spray%MW_f
    end if

    spray%Zmix_g = (spray%stoic_coeff*spray%Y_v/spray%Y_g-spray%Y_O2*spray%Y_a/spray%Y_g+spray%Y_O2)/(spray%stoic_coeff+spray%Y_O2)

    spray%Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

    ! Initialize turbulence
    ! Standard k-eps
    if(spray%ck_type == 'Dynamic') then
       write(*,*) 'Setting variable C_k...'
       ! Set some value to C_k, this will be anyway recalculated in turbulence model
       spray%c_k = 7.0_WP
    elseif(spray%ck_type == 'Const' .or. spray%ck_type == 'Constant') then
       if(spray%ck_val < 0.0_WP) then          
          write(*,*) 'Error!!! in initializing turbulence model....'
          write(*,*) 'Please provide a value for constant C_k...'
          write(*,*) 'C_k : <value>'
          call abort()
       end if
       write(*,*) 'Setting a constant C_k of ',spray%ck_val, '...'
       spray%c_k = spray%ck_val
    else
       spray%c_k = 7.0_WP
       write(*,*) 'Setting a constant C_k of ',maxval(spray%c_k), '...'
    end if
    spray%c_mu = 0.09_WP; spray%c_eps1 = 1.44_WP; spray%c_eps2 = 1.92_WP
    spray%c_zvar = 2.0_WP
    spray%k_g = 0.0_WP
    spray%zvar_g = 0.0_WP

    !spray%eps_g = sqrt(spray%c_k*spray%c_mu*spray%k_g*spray%rho)*spray%k_g/spray%b
    spray%eps_g = sqrt(spray%c_k*spray%c_mu*spray%k_g)*spray%k_g/spray%b

    !spray%k_g = 1.0E-04_WP
    !spray%eps_g = 1.0E-010_WP

    ! TEST: 
    !spray%mu_t_g = spray%c_mu*spray%rho*sqrt(spray%Y_a*spray%Y_g)*spray%k_g**2/spray%eps_g
    !spray%mu_t_g = spray%c_mu*sqrt(1.0/spray%DRa/spray%DRa)*spray%k_g**2/spray%eps_g
    !spray%mu_t_g = spray%c_mu/spray%DRg*spray%k_g**2/spray%eps_g

    !Average density according to Tamanini 1981
    spray%mu_t_g = spray%c_mu*sqrt(spray%rho_ra*spray%rho_l/spray%DRg)*spray%k_g**2/spray%eps_g*(0.5_WP*spray%U_inj**2*spray%tau)
    
    !Average density according to derivation in CnF2021 paper 
    !-> Tests didn't show any difference compared to Tamanini's suggestion for average density
    !spray%mu_t_g = spray%c_mu*spray%rho*spray%rho_l*spray%k_g**2/spray%eps_g*(0.5_WP*spray%U_inj**2*spray%tau)

    where (spray%eps_g == 0.0_WP) spray%mu_t_g = 0.0_WP

    spray%chi_g = spray%c_zvar*spray%eps_g/spray%k_g*spray%zvar_g
    where (spray%k_g == 0.0_WP) spray%chi_g = 0.0_WP

    spray%chi_g_stl = 0.0_WP

    ! Reference temperature and mass fraction for evaporation model
    call computeRefTemperature(spray)

    ! Compute initial droplet size distribution
    if(spray%dgf < 0.0_WP) spray%dgf = 5.0_WP
    call compute_DSD(spray)

    call computeTimeStep(spray)

    ! Post processing
    allocate(spray%time(1)); spray%time   = 0.0_WP
    allocate(spray%LPL(1));  spray%LPL    = 0.0_WP
    allocate(spray%VPL(1));  spray%VPL    = 0.0_WP
    allocate(spray%chi_st(1));  spray%chi_st = 0.0_WP
    allocate(spray%chi_st1(1));  spray%chi_st1 = 0.0_WP
    allocate(spray%utip(1));  spray%utip    = 0.0_WP
    allocate(spray%dtip(1));  spray%dtip    = 0.0_WP

    if(spray%combustion_model == 'MRIF') then
       allocate(spray%chi_st_m(spray%nsc));  spray%chi_st_m = 0.0_WP
    end if

    ! Compute gas mixture properties at reference temperature
    call computeGasMixtureProperties(spray,spray%T_ref)

    ! Compute non-dimensional parameters
    call compute_varNonDparams(spray)

    ! Get fuel vapor temperature at liquid temperature
    allocate(T_d(spray%nzo)); T_d = -9999.0_WP
    T_d = spray%Td*spray%T_fuel
    !call getVaporFuelTemperature(spray,T_d)
    deallocate(T_d)

    ! Initialize source terms to zeros
    spray%omega_ent  = 0.0_WP; spray%omega_vap  = 0.0_WP; spray%omega_vapdm  = 0.0_WP; 
    spray%omega_vapd2  = 0.0_WP; spray%omega_vapd3  = 0.0_WP; spray%f_drag  = 0.0_WP;
    spray%omega_bre1 = 0.0_WP; spray%omega_bre2 = 0.0_WP; spray%omega_T = 0.0_WP;

    if( spray%fixed_Re > 0.0_WP ) then
       spray%Re = spray%fixed_Re
    elseif( spray%fixed_We > 0.0_WP ) then
       spray%We = spray%fixed_We
    elseif( spray%fixed_DRa > 0.0_WP ) then
       spray%DRa = spray%fixed_DRa
       call computeGasMixtureProperties_fixed_DR_VR(spray,spray%T_ref)
    elseif( spray%fixed_DRv > 0.0_WP ) then
       spray%DRv = spray%fixed_DRv
       call computeGasMixtureProperties_fixed_DR_VR(spray,spray%T_ref)
    elseif( spray%fixed_VRa > 0.0_WP ) then
       spray%VRa = spray%fixed_VRa
       call computeGasMixtureProperties_fixed_DR_VR(spray,spray%T_ref)
    elseif( spray%fixed_VRv > 0.0_WP ) then
       spray%VRv = spray%fixed_VRv
       call computeGasMixtureProperties_fixed_DR_VR(spray,spray%T_ref)
    elseif( spray%fixed_De > 0.0_WP ) then
       spray%De = spray%fixed_De
    end if

    !!$$call compute_alpha(spray)
    if(spray%twoD .and. spray%irank .eq. 0) then
       call createSelfSimilarProfiles(spray)
    end if
    ! For RIF combustion model (MDUC)
#ifdef MDUC_MPI
    call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$    if(spray%twoD) then
!!$       call spray_combust_init_2D(spray)
!!$       call spray_combust_poststep_mrif_2D(spray)
!!$    else
    call spray_combust_init(spray)
    call spray_combust_poststep_mrif(spray)
!!$    end if
    
#endif

    ! Output important parameters
    if(spray%irank .eq. 0) call write_params(spray)

    ! Write initial output file 
    if(spray%irank .eq. 0) call write_output(spray,spray%step,spray%ndtime)

  end subroutine init_spray

  subroutine run_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: niter, ierr
    real(WP), dimension(spray%kmino:spray%kmaxo) :: T_d

    do while (spray%ndtime < spray%ndftime .and. .not.spray%end )
       
       write(*,*) 'step: ', spray%step, 'Time: ', spray%ndtime, 'dt: ', spray%dt !, 'End: ', spray%ndftime

       do niter = 1,spray%solver%rktvd%stage

          spray%solver%rktvd%niter = niter

          ! Apply BCs
          call applyBC(spray)

          ! Compute time step
          call computeTimeStep(spray)

          ! Compute Droplet Size Distribution
          call compute_DSD(spray)

          ! Update models: Entrainment, Drag, Breakup, Evaporation and heating
          call entrainmentTerm(spray)

          call dragModel(spray)

          call breakupModel(spray)

          call evaporationModel(spray)
          
          call turbulenceModel(spray)

          ! Update State Vector
          call buildStateVector(spray)

          ! Update Flux Vector
          call buildFluxVector(spray)

          ! Update Source Vector
          call buildSourceVector(spray)

          ! Solve multiple flamelets
          if(spray%nsc > 0) then
             call buildScalars(spray)
          end if

          ! Solve hyperbolic system of equations without source terms
          call solver_run(spray)

          ! Update flow variables
          call updateFlowVariables(spray)

          ! Update reference temperature
          call computeRefTemperature(spray)

          ! Update liquid fuel properties
          T_d = spray%Td*spray%T_fuel
          if (associated(spray%LFPT)) then
             call updateLiquidFuelPropertiesFromLFPT(spray,T_d)
          else
             call updateLiquidFuelProperties(spray,T_d)
          end if

          ! Get fuel vapor temperature at liquid temperature
          !call getVaporFuelTemperature(spray,T_d)

          ! Update reference properties
          if( spray%fixed_DRa > 0.0_WP .or. spray%fixed_DRv > 0.0_WP .or. &
              spray%fixed_VRa > 0.0_WP .or. spray%fixed_VRv > 0.0_WP ) then
             call computeGasMixtureProperties_fixed_DR_VR(spray,spray%T_ref)
          else
             call computeGasMixtureProperties(spray,spray%T_ref)
          end if

          ! Compute non-dimensional parameters
          call compute_constNonDparams(spray)

          ! Update varying Non-dimensional parameters
          call compute_varNonDparams(spray)
       
       end do

       ! Advance time 
       call advanceTime(spray)

       ! Compute Penetration Lengths
       call getPenetration(spray,spray%step)

       ! Conditional mean scalar dissipation rate 
       if(spray%combustion_model == 'MRIF') then
          call getScalarDissipationMRIF(spray,spray%step)
       else
          call getScalarDissipation(spray,spray%step)
       end if

       !!$$call compute_alpha(spray)
       
       if(spray%twoD .and. spray%irank .eq. 0) then
          call createSelfSimilarProfiles(spray)
       end if

       ! For RIF combustion model, call MDUC
#ifdef MDUC_MPI
       call MPI_barrier(MPI_COMM_WORLD,ierr)
       
!!$       if(spray%twoD) then
!!$          call spray_combust_poststep_mrif_2D(spray)
!!$       else
       call spray_combust_poststep_mrif(spray)
!!$       end if
#endif

       if(spray%irank .eq. 0) call write_output(spray,spray%step,spray%ndtime)

       call modify_inputs(spray)

    end do

    spray%end = .true.

  end subroutine run_spray

  ! Initial droplet size
  subroutine init_droplets(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    select case(spray%init_dsd_name)
       
    ! Blob injection: set to Log-Normal distribution with zero variance
    case ('Blob', 'blob')
        write(*,*) 'Initializing blob injection...'
        write(*,*) 'Log-Normal DSD will be used in the limit of zero variance...'
        write(*,*) 'For non-zero variance, computation will be continued with Log-Normal DSD...'
       spray%dsd_type = type_log_normal
       !spray%init_dm = 1.0_WP
       !spray%init_dvar = 0.001_WP
    ! Delta distribution
    case ('Delta','delta')
       write(*,*) 'Initializing delta (mono-dispersed) distribution for droplet sizes...'
       spray%dsd_type = type_delta
    ! Log-Normal distribution
    case ('Log Normal', 'Log-Normal', 'lognormal', 'LogNormal', 'log normal','LN','L-N')
       write(*,*) 'Initializing Log-Normal distribution for droplets...'
       spray%dsd_type = type_log_normal
    ! Rosin-Rammler distributtttttion
    case ('Rosin Rammler', 'Rosin-Rammler', 'RR', 'R-R')
       write(*,*) 'Initializing Rosin-Rammler distribution for droplet sizes...'
       spray%dsd_type = type_rosin_rammler
    ! Gamma distribution
    case ('Gamma','gamma')
       write(*,*) 'Initializing Gamma distribution for droplet sizes...'
       spray%dsd_type = type_gamma
    ! Generalized gamma distribution
    case ('Generalized Gamma', 'generalized gamma', 'Generalized gamma')
       write(*,*) 'Initializing Generalized Gamma distribution for droplet sizes...'
       spray%dsd_type = type_general_gamma
    ! Maximum-Entropy-Formalism
    case ('MEF','MEP','Maximum-Entropy-Formalism','Maximum-Entropy-Principle')
       write(*,*) 'Initializing Maximum entropy formalism for distribution of droplet sizes...'
       spray%dsd_type = type_maximum_entropy_formalism
    ! No distribution set
    case default
       write(*,*) 'Error!!! in setting initial droplet size distribution....'
       write(*,*) 'Please provide initial droplet size distribution...'
       write(*,*) 'Initial DSD : <value>'
       write(*,*) 'Available options for <value>:'
       write(*,*) '1 : Blob'
       write(*,*) '2 : Delta'
       write(*,*) '3 : Log-Normal'
       write(*,*) '4 : Rosin-Rammler'
       write(*,*) '5 : Gamma'
       write(*,*) '6 : Generalized Gamma'
       write(*,*) '7 : Maximum-Entropy-Formalism'
       call abort()

    end select

    if(spray%init_dm == -9999.0_WP) then
       write(*,*) 'Error!!! in initializing droplets....'
       write(*,*) 'Please provide initial Mean droplet diameter...'
       write(*,*) 'Initial Dm : <value>'
       call abort()
    end if

    if(spray%init_d2 == -9999.0 .and. spray%init_dvar == -9999.0_WP) then
       write(*,*) 'Error!!! in initializing droplets....'
       write(*,*) 'Please provide either initial variance of droplet diameter...'
       write(*,*) 'Initial Dvar : <value>'
       write(*,*) 'or initial Second moment of droplet diameter...'
       write(*,*) 'Initial D2 : <value>'
       call abort()
    end if

    ! Consistency between d2 and dvar
    if(spray%init_dvar >= 0.0_WP) then
       spray%init_d2 = spray%init_dm**2 + spray%init_dvar
    end if
    if(spray%init_d2 > 0.0_WP .and. spray%init_dvar < 0.0_WP) then
       spray%init_dvar = spray%init_d2 - spray%init_dm**2
    end if 

    if(spray%nd == -9999) then
       write(*,*) 'Error!!! in initializing droplets....'
       write(*,*) 'Please provide number of diameter classes for distribution...'
       write(*,*) 'Number of diameter classes : <value>'
       call abort()
    end if

  end subroutine init_droplets

  ! Initial droplet size
  subroutine initDropletSize(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: eps = 1E-16_WP
    real(WP), pointer :: Re=>null(), We=>null()
    real(WP) :: u_rel

    real(WP) :: Weg, Wel, Rel, Z, T, OmegaKH, LambdaKH, tauKH, dKH
    real(WP) :: B0, B1, C1, C2, C3, D10, D20, D30, pw

    Re => spray%Re; We => spray%We
             
    ! Kelvin-Helmoltz WAVE model
    B0 = 0.61_WP; B1 = 10.0_WP
    Weg = 0.5_WP*We/spray%DRa
    Wel = 0.5_WP*We
    Rel = 0.5_WP*Re
    Z = sqrt(Wel)/(Rel)
    T = Z*sqrt(Weg)
    OmegaKH = (0.34_WP+0.38_WP*Weg**1.5_WP)/(1.0_WP+Z)/(1.0_WP+1.4_WP*T**0.6_WP)*sqrt(8.0_WP/We) !/spray%noz_D**3)
    LambdaKH = (9.02_WP/2.0_WP)*(1.0_WP+0.45_WP*sqrt(Z))*(1.0_WP+0.4_WP*T**0.7_WP)/(1.0_WP+0.865_WP*Weg**1.67_WP)**0.6_WP
    tauKH = (3.726_WP*B1/2.0_WP)/(OmegaKH*LambdaKH)
    dKH = 2.0_WP*B0*LambdaKH

    if ( dKH > 1.0_WP ) then
       dKH = min((3.0_WP*PI/OmegaKH)**(1.0_WP/3.0_WP),(1.5_WP*LambdaKH)**(1.0_WP/3.0_WP))
    end if

    C1 = 14.0807_WP !2.0824e+04 !0.1085_WP !14.0206_WP !6.4723_WP 
    C2 = 15.8461_WP !2.3434e+04 !0.1221_WP !15.8105_WP !7.2934_WP 
    C3 = 17.4024_WP !2.5736e+04 !0.1341_WP !17.3763_WP !8.0139_WP

    pw = 1.0_WP !2.0_WP

    D10 = min(1.0_WP,C1*dKH**pw);
    D20 = min(1.0_WP,C2*dKH**pw);
    D30 = min(1.0_WP,C3*dKH**pw);

    spray%init_dm = D10
    spray%init_d2 = D20**2
    spray%init_d3 = D30**3

    write(*,*) 'Initial droplet size D10: ', D10, 'D20: ', D20, 'D30:', D30

  end subroutine initDropletSize

subroutine nozzle_flow_model(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    ! ---------------------------------
    real(WP) :: D, Cd, u_mean, u_eff, d_eff
    real(WP) :: P_inj, P_amb, P_vap, CN
    real(WP) :: Tcrit, Tfuel, Tboil, rhol, viscl, FI, delT_star

    D = spray%noz_D; Cd = spray%C_d

    P_inj = spray%P_inj; P_amb = spray%P_a; P_vap = spray%p_vap(1); rhol = spray%rho_l; viscl = spray%visc_l
    Tcrit = spray%T_crit; Tfuel = spray%T_fuel; Tboil = spray%T_boil

    u_mean = sqrt(2.0_WP*(P_inj-P_amb)/rhol)

    u_eff = Cd*sqrt(2.0_WP*(P_inj-P_amb)/rhol)
    
    CN = (P_inj-P_vap)/(P_inj-P_amb)

    delT_star = (Tfuel-Tboil)/(Tcrit-Tboil)

    FI = delT_star-CN

    d_eff = 1.0e-6_WP*(118.40_WP-28.29_WP*FI)

    spray%D_eff = d_eff

    spray%U_inj = u_eff

  end subroutine nozzle_flow_model

  subroutine nozzle_flow_modelAD(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: Cc0 = 0.61_WP, tol = 1e-05_WP
    real(WP) :: D, rd, ld, Kin, Kexp, Cc, u_mean, u_vena, u_eff
    real(WP) :: p_vena, p_venar, p_vap, p1, Aeff, Cd, Cdold, err, d_eff, Re
    real(WP) :: P_inj, P_amb, rho_l, visc_l, f
    
    D = spray%noz_D
    rd = spray%noz_rD; ld = spray%noz_LD
    Kin = spray%K_in; Kexp = spray%K_exp
    
    P_inj = spray%P_inj; P_amb = spray%P_a
    rho_l = spray%rho_l; visc_l = spray%visc_l

    p_vap = spray%p_vap(1);

    if (rd <= 0.15_WP) then
       Cc = ((1.0_WP/Cc0)**2 - 11.4_WP*rd)**(-0.5_WP)
    else
       Cc = 1.0_WP
    end if

    spray%Cc0 = Cc0; spray%Cc = Cc
    
    ! Initial conditions
    Cd = 1.0
    u_mean = sqrt(2.0_WP*(P_inj-P_amb)/rho_l)
    err = 1.0_WP
    d_eff = D

    do while (err > tol)
       Re = rho_l*u_mean*d_eff/visc_l
       f = max(0.316_WP*Re**(-0.25_WP),64.0_WP/Re)
       Cd = 1.0_WP/sqrt(Kin+f*ld+1.0_WP)
       u_mean = Cd*sqrt(2.0_WP*(P_inj-P_amb)/rho_l)
       u_vena = u_mean/Cc
       p_vena = P_inj - 0.5_WP*rho_l*u_vena**2
       p_venar = P_amb + 0.5_WP*rho_l*u_mean**2 *(1.0_WP-1.0_WP/Cc + Kexp + f*ld)
       if(p_vena < p_vap) then
          write(*,*) '########  Cavitation detected....'
          p1 = p_vap + 0.5_WP*rho_l*u_vena**2
          Cd = Cc*sqrt((p1-p_vap)/(p1-P_amb))
          u_eff = u_vena - (P_amb-p_vap)/rho_l/u_mean;
          Aeff = 0.25_WP*Pi*D**2 *u_mean/u_eff
          u_mean = u_eff
          d_eff = sqrt(4.0_WP*Aeff/Pi);
       end if
       err = abs(Cdold -Cd)
       Cdold = Cd
    end do

    spray%C_d = Cd

    spray%U_inj = u_mean
    
    spray%D_eff = d_eff

  end subroutine nozzle_flow_modelAD

  subroutine injection_params(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP) :: A_eff, inj_mass, peak_mfr
    integer :: i

    A_eff = spray%C_d*spray%Anoz

    spray%D_eff = sqrt(4.0_WP*A_eff/Pi)

    spray%tau = 1.0_WP

    call read_ROI_from_file(spray)

    inj_mass = 0.0_WP
    do i=2,size(spray%roi,1)
        inj_mass = inj_mass + 0.5_WP*(spray%roi(i,2)+spray%roi(i-1,2)) &
                                    *(spray%roi(i,1)-spray%roi(i-1,1))
    end do

    peak_mfr = spray%inj_mass/inj_mass/spray%num_noz

    spray%U_inj = peak_mfr/spray%rho_l/A_eff

    spray%tau = spray%D_eff/spray%U_inj

    spray%roi(:,1) = spray%roi(:,1)/spray%tau

  end subroutine

  subroutine compute_constNonDparams(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    ! Fuel Jet exit Reynolds number
    spray%Re = spray%rho_l*spray%U_inj*spray%D_eff/spray%visc_l
    
    ! Fuel Jet Weber number
    spray%We = spray%rho_l*spray%U_inj**2*spray%D_eff/spray%sigma                 
    ! Density Ratio liquid to ambient gas
    spray%DRa = spray%rho_l/spray%rho_a
    
    ! Density Ratio liquid to vapor
    spray%DRv = spray%rho_l/spray%rho_v
 
    ! Viscosity Ratio liquid to ambient gas
    spray%VRa = spray%visc_l/spray%visc_a
    
    ! Viscosity Ratio liquid to vapor
    spray%VRv = spray%visc_l/spray%visc_v
    
    ! Ratio of Molecular weights of ambient to fuel
    spray%WR = spray%MW_a/spray%MW_f
    
    ! New non-dimensional parameter
    spray%De = spray%L_f(1)/(spray%T_fuel*spray%C_l(1))                              
    ! Prandtl number 
    !spray%Pr_l = spray%visc_l*spray%C_l(1)/spray%lambda_l

  end subroutine compute_constNonDparams

  subroutine compute_varNonDparams(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    ! Density Ratio liquid to gas phase mixture
    spray%DRg = spray%rho_l/spray%rho_g

    ! Viscosity Ratio liquid to gas phase mixture
    spray%VRg = spray%visc_l/spray%visc_g

    ! Ratio of Latent Heats at droplet temperature to fuel temperature
    spray%LR = spray%L_f/spray%L_f(1)

    ! Ratio of Heat Capacities at droplet temperature to fuel temperature
    spray%CR = spray%C_l(1)/spray%C_l

    ! Viscosity Ratio liquid to gas phase turbulent viscosity
    spray%VRtg = spray%mu_t_g/spray%visc_l

    ! Ratio of surface tension at droplet temperature to fuel temperature
    spray%SR = spray%sigma/spray%sigma_loc
    where(spray%sigma_loc==0.0_WP) spray%SR = 0.0_WP

    ! Ratio of surface tension at droplet temperature to fuel temperature
    spray%DRl = spray%rho_l_loc/spray%rho_l

    ! Ratio of viscosity at droplet temperature to fuel temperature
    spray%VRl = spray%visc_l_loc/spray%visc_l
   
  end subroutine compute_varNonDparams

  subroutine compute_b(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    real(WP) :: numerator, denominator, factor
    integer :: k

    do k=spray%kmino,spray%kmaxo

       denominator = sqrt(spray%rho(k)*(spray%solver%F(3,k)+spray%solver%F(5,k)))
       numerator = spray%solver%F(1,k)+spray%solver%F(2,k)+spray%solver%F(4,k)

       if((spray%Y_l(k)*spray%u_l(k) + spray%Y_g(k)*spray%u_g(k)) > 0.0_WP ) then
          factor = sqrt(spray%Y_l(k)*spray%u_l(k)**2 + spray%Y_g(k)*spray%u_g(k)**2)/(spray%Y_l(k)*spray%u_l(k) + spray%Y_g(k)*spray%u_g(k))
       end if
       
       if(denominator > 0.0_WP ) then
          spray%b(k) = max(0.5_WP,numerator/denominator*factor)
       else
          spray%b(k) = 0.5_WP
       end if

    end do

  end subroutine compute_b

  subroutine compute_beta(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP) :: dd, ld, DRa ,theta,beta
    real(WP) :: A, gammaa, fgamma
    real(WP) :: Rp, Rpstar, kb, NdSrf, zeta, a0, a1, a2, a3, k

    dd = spray%noz_D/spray%noz_Dsac; ld = spray%noz_LD
    DRa = spray%DRa

    ! Several models for spray cone angle
    ! Ref.: Spray Atomization Models in Engine Applications,
    !       from Correlations to Direct Numerical Simulations
    ! theta is half cone angle
    ! beta is tan(theta)
    select case (trim(spray%spray_angle_model))
    case ('Reitz-Bracco')

       A = 3.0_WP + 0.28_WP*(ld)

       gammaa = (spray%Re/spray%We)**2*DRa

       fgamma = sqrt(3.0_WP)/6.0_WP*(1.0_WP-exp(-10.0_WP*gammaa))

       spray%beta = 4.0_WP*(Pi/A)*sqrt(1.0_WP/DRa)*fgamma

       spray%theta = atan(spray%beta)*180.0_WP/Pi

    case ('Reitz-Bracco-Simplified')

       A = 3.0_WP + 0.28_WP*(ld)

       fgamma = sqrt(3.0_WP)/6.0_WP

       spray%beta = 4.0_WP*(Pi/A)*sqrt(1.0_WP/DRa)*fgamma

       spray%theta = atan(spray%beta)*180.0_WP/Pi

    case ('Ruiz-Chigier')

       A = 3.0_WP + 0.28_WP*(ld)

       gammaa = (spray%Re/spray%We)**2*DRa

       fgamma = sqrt(3.0_WP)/6.0_WP*(1.0_WP-exp(-10.0_WP*gammaa))

       spray%beta = 4.0_WP*(Pi/A)*sqrt(1.0_WP/DRa)*fgamma*(spray%Re/spray%We)**(-0.25_WP)

       spray%theta = atan(spray%beta)*180.0_WP/Pi

    case ('Arai')

       spray%theta = 0.025_WP*(spray%Re**2/DRa*spray%VRa**2/2.0_WP/spray%C_d)**0.25_WP

       spray%beta = tan(spray%theta*Pi/180.0_WP)

    case ('Hiroyasu-Arai','Hiroyasu')

       spray%Cnoz = 41.75_WP*(ld**(-0.22_WP))*(dd)**0.15_WP

       spray%theta = spray%Cnoz*(1.0_WP/DRa)**0.26_WP

       spray%beta = tan(spray%theta*Pi/180.0_WP)

    case ('Arregle')  

       spray%beta = spray%noz_D**0.508_WP*spray%P_inj**0.00943_WP*spray%rho_a**0.335_WP

       spray%theta = atan(spray%beta)*180.0_WP/Pi

    case ('Siebers')

       if (spray%C_theta .lt. 0.0_WP) spray%C_theta = 1.0_WP

       spray%beta = spray%C_theta*((1.0_WP/DRa)**0.19_WP - 0.0043_WP*(DRa)**0.5_WP)

       spray%theta = atan(spray%beta)*180.0_WP/Pi
    case ('Reddemann')

       spray%beta = 0.212_WP*(1.0_WP/DRa)**0.198_WP*spray%Re**0.0526_WP

       spray%theta = atan(spray%beta)*180.0_WP/Pi

  case ('Kitamura')

       Rp = spray%p_vap(1)/spray%P_a
       Rpstar = (spray%P_inj-spray%P_a)/(spray%p_vap(1)-spray%P_a)

       a0 = ((36.0_WP*Pi)*(spray%MVol)**2.0_WP)**(1.0_WP/3.0_WP)
       kb = 1.380649e-23_WP
       NdSrf = a0*spray%sigma/(kb*spray%T_fuel)
       zeta = log10((NdSrf**3.0_WP)*(Rp**2.0_WP)/(spray%MW_f**2.0_WP))
       spray%theta = 0.5_WP*(zeta) 
       spray%beta = tan(spray%theta*Pi/180.0_WP)

  !!===============================!!
  !!  For varying spray cone angle !!
  !!    From Kamoun et al. 2010    !!
  !!===============================!!
  !!do  k=spray%kmin,spray%kmax
  !!  if(spray%z(k).le.2)then
  !!     a1 = -3.234_WP; a2 = 363.07_WP; a3 = -10021.0_WP
  !!     spray%theta(k) = 0.5_WP*(a1*zeta**2.0_WP+a2*zeta+a3)  !  !0.5_WP*(zeta) !
  !!     spray%beta(k) = tan(spray%theta(k)*Pi/180.0_WP)
  !!  else if(spray%z(k).gt.2 .and. spray%z(k).le.5)then
  !!    a1 = -2.368_WP; a2 = 272.084_WP; a3 = -7649.0_WP
  !!    spray%theta(k) = 0.5_WP*(a1*zeta**2.0_WP+a2*zeta+a3)  !  !0.5_WP*(zeta) !
  !!    spray%beta(k) = tan(spray%theta(k)*Pi/180.0_WP)
  !!  else if(spray%z(k).gt.5 .and. spray%z(k).le.10)then
  !!     a1 = -3.208_WP; a2 = 366.61_WP; a3 = -10324.0_WP
  !!     spray%theta(k) = 0.5_WP*(a1*zeta**2.0_WP+a2*zeta+a3)  !  !0.5_WP*(zeta) !
  !!      spray%beta(k) = tan(spray%theta(k)*Pi/180.0_WP)
  !!  else if(spray%z(k).gt.10 .and. spray%z(k).le.15)then
  !!     a1 = -3.043_WP; a2 = 348.364_WP; a3 = -9829.0_WP
  !!     spray%theta(k) = 0.5_WP*(a1*zeta**2.0_WP+a2*zeta+a3)  !  !0.5_WP*(zeta) !
  !!    spray%beta(k) = tan(spray%theta(k)*Pi/180.0_WP)
  !!  else if(spray%z(k).gt.15 .and. spray%z(k).le.20)then
  !!     a1 = -2.851_WP; a2 = 332.92_WP; a3 = -9559.0_WP
  !!    spray%theta(k) = 0.5_WP*(a1*zeta**2.0_WP+a2*zeta+a3)  !  !0.5_WP*(zeta) !
  !!    spray%beta(k) = tan(spray%theta(k)*Pi/180.0_WP)
  !!  else 
  !!     a1 = -4.353_WP; a2 = 497.69_WP; a3 = -14096.0_WP
  !!    spray%theta(k) = 0.5_WP*(a1*zeta**2.0_WP+a2*zeta+a3)  !  !0.5_WP*(zeta) !
  !!     spray%beta(k) = tan(spray%theta(k)*Pi/180.0_WP)
  !! end if 
  !!end do


    case ('Constant','constant')

       if (spray%beta .lt. 0.0_WP .and. spray%theta .gt. 0.0_WP) then
          spray%beta = tan(spray%theta*Pi/180.0_WP)
          write(*,*) '#### Spray half-cone angle provided by user....'
          write(*,*) '#### Computing tangent of spray half-cone angle....'
          write(*,*) '#### theta = ',spray%theta, 'beta = ', spray%beta
       else if (spray%beta .gt. 0.0_WP .and. spray%theta .lt. 0.0_WP) then
          spray%theta = atan(spray%beta)*180.0_WP/Pi
          write(*,*) '#### Tangent of spray half-cone angle provided by user....'
          write(*,*) '#### Computing spray half-cone angle....'
         write(*,*) '#### theta = ',spray%theta, 'beta = ', spray%beta
       else
          write(*,*) 'Error!!! in computing spray half-cone angle....'
          write(*,*) 'Please provide spray half-cone angle or check the routine to compute it...'
          write(*,*) 'Spray half-cone angle : <value>'
          write(*,*) 'OR'
          write(*,*) 'Tangent of spray half-cone angle : <value>'
          call abort
       end if
       
!!$    case ('Varying')
!!$
!!$       

    case default

       write(*,*) 'Error!!! No spray angle model provided....'
       write(*,*) 'Please choose one of the following spray angle models:'
       write(*,*) 'Spray half-cone angle model : <name>'
       write(*,*) 'Reitz-Bracco, Reitz-Bracco-Simplified, Ruiz-Chigier, Arai,...'
       write(*,*) 'Hiroyasu-Arai, Arregle, Siebers, Constant'
       call abort        

    end select

    !write(*,*) '#### Spray half-cone angle model: ',trim(spray%spray_angle_model)
    !write(*,*) '#### theta = ',spray%theta, 'beta = ', spray%beta


  end subroutine compute_beta

  subroutine compute_DSD(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(spray%nzo):: Dbar, q, ratio, Dm, D2, D3, Dvar, mu, sigma2, B, alpha, beta, gamm, var, skew
    real(WP), dimension(spray%nd):: hh, f, weight, DD
    real(WP), dimension(spray%nd+1):: D
    real(WP), dimension(spray%nd,spray%nzo) :: dsd, dsdln
    real(WP) :: L, h, norm, r_min, r_max, peak, diff, eps = 1.0E-16_WP
    integer :: i, k, n, idx

    if(spray%step==0) write(*,*) 'Computing droplet size distribution...'

    Dm = spray%dm
    D3 = spray%d3
    Dvar = spray%dvar
    !where(Dvar < 0.0_WP) Dvar = 0.0_WP
    !D2 = Dm**2+Dvar

    D2 = spray%d2
    Dvar = D2-Dm**2

    do k=spray%kmino,spray%kmaxo

       if(spray%init_dsd_name .eq. 'blob' .or. spray%init_dsd_name .eq. 'Blob') then
!!$          if(D2(k)/Dm(k)**2 > 1.05_WP .and. D2(k)/Dm(k)**2 < 1.9_WP) then
!!$             spray%dsd_type(k) = type_rosin_rammler
!!$          end if
       end if

       if((trim(spray%init_dsd_name) == trim('Rosin-Rammler')) .and. Dvar(k) < 1E-8_WP .and. Dm(k) > 0.0_WP) then
          spray%dsd_type(k) = type_delta
       elseif (trim(spray%init_dsd_name) == trim('Rosin-Rammler')) then
          spray%dsd_type(k) = type_rosin_rammler
       end if

       select case(spray%dsd_type(k))
       case (type_rosin_rammler)

          !if(spray%step==0) write(*,*) 'Computing Rosin-Rammler DSD...'

          Dbar = 0.0_WP

          if(Dvar(k) > 0.0005) then
             write(*,*) 'Dvar'
          end if

          L = min(1.0_WP,spray%dgf*Dm(k))

          !call RosinRammlerPDF(Dm(k),D2(k),L,spray%nd,DD,f)

          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h

          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)
          hh = h

          !call createDropletGrid(Dm(k),D2(k),L,spray%nd,D,spray%di(:,k))
          !hh = D(2:spray%nd+1)-D(1:spray%nd)
          !spray%di(:,k) = (D(:spray%nd)+0.5_WP*hh)

          dsd(:,k) = 0.0_WP

          r_min = 1.05_WP; r_max = 1.9_WP
          ratio(k) = D2(k)/Dm(k)**2

          if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP) then

             if (ratio(k) > r_min ) then
                ratio(k) = min(r_max,ratio(k))
                call bisection(q(k),ratio(k))
             else
                q(k) = 5.1334_WP
             end if

             Dbar(k) = Dm(k)/gamma(1.0_WP+1.0_WP/q(k))

             f = (q(k)/Dbar(k)**q(k))*spray%di(:,k)**(q(k)-1.0_WP)*exp(-(spray%di(:,k)/Dbar(k))**q(k))

             norm = sum(f*hh)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*hh/norm;
             end if

          end if

          spray%dsd(:,k) = dsd(:,k)

       case (type_rosin_rammler+type_log_normal)

          !if(spray%step==0) write(*,*) 'Computing blended Rosin-Rammler Log-Normal DSD...'

          Dbar = 0.0_WP

          L = min(1.0_WP,spray%dgf*Dm(k))

          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h

          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)

          dsd(:,:) = 0.0_WP

          r_min = 1.05_WP; r_max = 1.9_WP
          ratio(k) = D2(k)/Dm(k)**2

          if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP) then

             if (ratio(k) > r_min ) then
                ratio(k) = min(r_max,ratio(k))
                call bisection(q(k),ratio(k))
             else
                q(k) = 5.1334_WP
             end if

             Dbar(k) = Dm(k)/gamma(1.0_WP+1.0_WP/q(k))

             f = (q(k)/Dbar(k)**q(k))*D**(q(k)-1.0_WP)*exp(-(D/Dbar(k))**q(k))

             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*h/norm;
             end if

          end if

          dsdln(:,:) = 0.0_WP

          mu = log(Dm**2/sqrt(D2))
          sigma2 = log(D2/Dm**2)

          B = 1.0_WP/sqrt(twoPi*sigma2)

          if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP) then

             f = B(k)*exp(-((log(D)-mu(k))**2/2.0_WP/sigma2(k)))/D

             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsdln(:,k) = f/norm;
             end if

          end if

          if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP) then
             peak = maxval(dsdln(:,k))
             weight = dsdln(:,k)/peak
             do n=2,spray%nd
                if(weight(n-1) == 1.0_WP) then
                   weight(n) = 1.0_WP
                end if
             end do
          end if

          !if (k==1) write(*,*) weight

          f = (1.0_WP-weight)*dsdln(:,k) + weight*dsd(:,k)
          norm = sum(f*h)
          if (norm > 0.0_WP) then
             dsd(:,k) = f/norm;
          end if

       case (type_log_normal)

          !if(spray%step==0) write(*,*) 'Computing Log-Normal DSD...'

          L = min(1.0_WP,spray%dgf*Dm(k))

          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h

          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)
          hh = h
          !call createDropletGrid(Dm(k),D2(k),L,spray%nd,D,spray%di(:,k))
          !hh = D(2:spray%nd+1)-D(1:spray%nd)
          !spray%di(:,k) = (D(:spray%nd)+0.5_WP*hh)

          dsd(:,:) = 0.0_WP

          mu(k) = log(Dm(k)**2/sqrt(D2(k)))
          sigma2(k) = log(D2(k)/Dm(k)**2)

          B(k) = 1.0_WP/sqrt(twoPi*sigma2(k))

          if (Dm(k) > 0.0_WP .and. D2(k) >= 0.0_WP) then

             if (sigma2(k) <= 0.0_WP) then

                call getIndex(D,Dm(k),idx)

                f(idx) = 1.0_WP

             else

                f = B(k)*exp(-((log(spray%di(:,k))-mu(k))**2/2.0_WP/sigma2(k)))/spray%di(:,k)

             end if

             norm = sum(f*hh)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*hh/norm;
             end if

          end if

          spray%dsd(:,k) = dsd(:,k)

       case (type_gamma)

          !if(spray%step==0) write(*,*) 'Computing Gamma DSD...'

          L = min(1.0_WP,spray%dgf*Dm(k))

          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h

          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)

          dsd(:,:) = 0.0_WP

          alpha(k) = Dm(k)**2/(D2(k)-Dm(k)**2)
          beta(k) = Dm(k)/(D2(k)-Dm(k)**2)
          B(k) = beta(k)**alpha(k)/gamma(alpha(k))

          !if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP .and. abs(D2(k)-Dm(k)**2) > 0.0_WP ) then
          !if((beta(k)**alpha(k) < 1000.0_WP) .and. (beta(k)**alpha(k) > 1.0E-3_WP) .and. (gamma(alpha(k)) < 1000.0_WP) .and. (gamma(alpha(k)) > 1.0E-3_WP)) then
          if(.not. isnan(B(k)) .and. (B(k)>0.0_WP) .and. .not. isinf(B(k)) ) then

                f = B(k)*spray%di(:,k)**(alpha(k)-1.0_WP)*exp(-beta(k)*spray%di(:,k))

                norm = sum(f*h)
                if (norm > 0.0_WP) then
                   dsd(:,k) = f*h/norm;
                end if

          else if ((D2(k)-Dm(k)**2) > 0.0_WP) then

             mu(k) = Dm(k)
             sigma2(k) = D2(k)-Dm(k)**2
             B(k) = 1.0_WP/sqrt(twoPi*sigma2(k))
             f = B(k)*exp(-((spray%di(:,k))-mu(k))**2/2.0_WP/sigma2(k))
             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*h/norm;
             end if

          else

             call getIndex(D,Dm(k),idx)

             f(idx) = 1.0_WP

             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*h/norm
             end if

          end if

          spray%dsd(:,k) = dsd(:,k)

!!$       case (type_inverse_gamma)
!!$
!!$          !if(spray%step==0) write(*,*) 'Computing Inverse Gamma DSD...'
!!$
!!$          L = min(1.0_WP,spray%dgf*Dm(k))
!!$
!!$          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h
!!$
!!$          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
!!$          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
!!$          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)
!!$
!!$          dsd(:,:) = 0.0_WP
!!$
!!$          alpha(k) = Dm(k)**2/(D2(k)-Dm(k)**2)
!!$          !alpha(k) = max(1.0_WP+1.0E-12_WP,min(100.0_WP,alpha(k)))
!!$          beta(k) = Dm(k)/(D2(k)-Dm(k)**2)
!!$          !beta(k) = max(1.0_WP,min(12.0_WP*alpha(k),beta(k)))
!!$
!!$
!!$          !if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP .and. abs(D2(k)-Dm(k)**2) > 0.0_WP ) then
!!$          if((beta(k)**alpha(k) < 1000.0_WP) .and. (beta(k)**alpha(k) > 1.0E-3_WP) .and. (gamma(alpha(k)) < 1000.0_WP) .and. (gamma(alpha(k)) > 1.0E-3_WP)) then
!!$             B(k) = beta(k)**alpha(k)/gamma(alpha(k))
!!$
!!$             f = B(k)*D**(alpha(k)-1.0_WP)*exp(-beta(k)*D)
!!$
!!$             norm = sum(f*h)
!!$             if (norm > 0.0_WP) then
!!$                dsd(:,k) = f*h/norm;
!!$             end if
!!$
!!$          else if ((D2(k)-Dm(k)**2) > 0.0_WP) then
!!$             mu(k) = Dm(k)
!!$             sigma2(k) = D2(k)-Dm(k)**2
!!$             B(k) = 1.0_WP/sqrt(twoPi*sigma2(k))
!!$             f = B(k)*exp(-((spray%di(:,k))-mu(k))**2/2.0_WP/sigma2(k))
!!$             norm = sum(f*h)
!!$             if (norm > 0.0_WP) then
!!$                dsd(:,k) = f*h/norm;
!!$             end if
!!$          else
!!$
!!$             call getIndex(D,Dm(k),idx)
!!$
!!$             f(idx) = 1.0_WP
!!$
!!$             norm = sum(f*h)
!!$             if (norm > 0.0_WP) then
!!$                dsd(:,k) = f*h/norm
!!$             end if
!!$
!!$          end if
!!$
!!$          spray%dsd(:,k) = dsd(:,k)

       case (type_general_gamma)!'GenGamma','Generalized Gamma')

          !if(spray%step==0) write(*,*) 'Computing Generalized Gamma DSD with three parameters...'

          L = min(1.0_WP,spray%dgf*Dm(k))

          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h

          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)

          dsd(:,:) = 0.0_WP

          var(k) = D2(k)-Dm(k)**2
          skew(k) = (D3(k)-3.0_WP*Dm(k)*D2(k)+2.0_WP*Dm(k)**3)/var(k)**1.5
          alpha(k) = 4.0_WP/skew(k)**2
          beta(k) = sqrt(var(k)/alpha(k))
          gamm(k) = Dm(k) - 2.0_WP*sqrt(var(k))/skew(k)

          if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP) then

             B(k) = 1.0_WP/(beta(k)**alpha(k))/gamma(alpha(k))

             f = B(k)*(D-gamm(k))**(alpha(k)-1.0_WP)*exp(-(D-gamm(k))/beta(k))

             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*h/norm;
             end if

          end if

          spray%dsd(:,k) = dsd(:,k)

       case (type_maximum_entropy_formalism)

          !if(spray%step==0) write(*,*) 'Computing DSD using maximum entropy formalism...'

          dsd(:,:) = 0.0_WP

          if (Dm(k) > 0.0_WP .and. D2(k) > 0.0_WP) then

             call maximumEntropyFormalism(spray,Dm(k),D2(k),D3(k),k,D,f)

             D(1) = 1E-16_WP

             h = D(2)-D(1); spray%h(k) = h

             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*h/norm;
             end if

          end if

          dsd(1,:) = 0.0_WP

       case (type_delta)

          !write(*,*) 'Computing Delta DSD...'

          L = Dm(k) !sqrt(D2(k))

          h = 2.0_WP*L/(2.0_WP*spray%nd-1.0_WP); spray%h(k) = h

          D = (/ (h*real(i,WP),i=0,spray%nd,1) /)
          f = (/ (real(0.0,WP),i=1,spray%nd,1) /)
          spray%di(:,k) = (D(:spray%nd)+0.5_WP*h)

          dsd(:,k) = 0.0_WP

          if (Dm(k) > 0.0_WP) then

             call getIndex(D,Dm(k),idx)

             f(idx) = 1.0_WP

             norm = sum(f*h)
             if (norm > 0.0_WP) then
                dsd(:,k) = f*h/norm
             end if

          end if

          spray%dsd(:,k) = dsd(:,k)

       case default

          write(*,*) 'Error!!! Please select initial droplet size distribution...'
          write(*,*) 'Initial DSD : <value>'
          write(*,*) 'Available options for <value> : '
          write(*,*) 'Delta, Rosin-Rammler, Log-Normal, Gamma, Generalized Gamma, Maximum-Entropy-Formalism'
          call abort()

       end select

    end do

  end subroutine compute_DSD

  subroutine createDropletGrid(xm,x2,xmax,nx,x_,x)
    implicit none

    ! ---------------------------------
    real(WP), intent(in) :: xm, x2, xmax
    integer, intent(in) :: nx
    real(WP), dimension(:), intent(out) :: x_,x
    ! ---------------------------------

    real(WP), dimension(nx+1) :: f, dx_, dx, x_adapt, f_adapt, dx_adapt
    real(WP), dimension(nx) :: invdiff
    real(WP) :: dx_min, dx_max, peak_x, peak_xn
    real(WP) :: dxx, mu, sigma2, B, xbar, r_min, r_max, q, ratio, norm
    real(WP) :: Integral, bin, dummy_bin, dx_guess1, dx_guess2, dx_guess3, err_bin, err_bin2, err_case, x_left, x_right, f_left, f_right, eps
    integer :: i,idx,ind,kk

    dx_min = 0.0001_WP*xm
    dx_max = 10.0_WP*xmax/nx

    dxx = 2.0_WP*xmax/(2.0_WP*nx-1.0_WP)

    x_ = (/ (dxx*real(i,WP),i=0,nx,1) /)
    f = (/ (real(0.0,WP),i=1,nx,1) /)
    x = (x_(:nx)+0.5_WP*dxx)

    mu = log(xm**2/sqrt(x2))
    sigma2 = log(x2/xm**2)
    if(sigma2 <= 0.0_WP .or. isnan(sigma2)) then
       return
    end if

    call lognormalPDF(xm,mu,sigma2,x_,f)

    eps = 1.0E-05_WP

    Integral = 0
    do i=1,nx
       Integral = Integral + 0.5_WP*(f(i)+f(i+1))*dxx
    end do
    bin = Integral/nx

    do i=1,nx
       
       x_left = x_adapt(i)
       f_left = f_adapt(i)

       dx_guess1 = dxx;
       x_right = x_adapt(i) + dx_guess1
       if (x_right == 0) then
          f_right = 0
       else
          call lognormalPDFscalar(mu,sigma2,x_right,f_right)
       end if

       dummy_bin = 0.5*(f_left+f_right)*dx_guess1
       err_bin = dummy_bin - bin

       kk = 0
       do while (abs(err_bin) > eps)
          kk = kk+1
          if (kk==1) then
             dx_guess2 = dx_guess1
             if (err_bin > eps) then
                dx_guess1 = dx_guess2/2            
             end if
             if (err_bin <  -eps) then
                dx_guess1 = dx_guess2*2   
             end if
          else
             err_case = err_bin2*err_bin
             dx_guess3 = dx_guess2
             dx_guess2 = dx_guess1
             if (err_case > 0) then
                if (err_bin > eps) then
                   dx_guess1 = dx_guess2/2
                end if
                if (err_bin <  -eps) then
                   dx_guess1 = dx_guess2*2
                end if
             else
                dx_guess1 = (dx_guess3+dx_guess2)/2
             end if
          end if

          x_right = x_adapt(i) + dx_guess1
          if(x_right == 0) then
             f_right = 0
          else
             call lognormalPDFscalar(mu,sigma2,x_right,f_right)
          end if
          dummy_bin = 0.5*(f_left+f_right)*dx_guess1
          err_bin2 = err_bin
          err_bin = dummy_bin - bin

       end do

       x_adapt(i+1) = x_right
       f_adapt(i+1) = f_right
       dx_adapt(i+1) = dx_guess1
    end do

    x_ = x_adapt

!!$    open(unit=1000,file='DSD.grid',form="formatted",status="unknown",action="write")
!!$    do i=1,nx
!!$       write(1000,FMT='(ES15.5E3,ES15.5E3,ES15.5E3,ES15.5E3)') x_(i),f(i),x_adapt(i),f_adapt(i)
!!$    end do
!!$    close(unit=1000)

!!$    B = 1.0_WP/sqrt(twoPi*sigma2)
!!$
!!$    if (xm > 0.0_WP .and. sigma2 .gt. 0.0_WP) then
!!$       
!!$       f = B*exp(-((log(x)-mu)**2/2.0_WP/sigma2))/x
!!$
!!$       norm = sum(f*dxx)
!!$       if (norm > 0.0_WP) then
!!$          f = f/norm;
!!$       end if
!!$
!!$       ind = maxloc(f);
!!$       peak_x = x(ind);
!!$
!!$       ! Adapt 
!!$       invdiff = -f+maxval(f)+1.0_WP; !1.0_WP/(abs(f(:nx-1)-f(2:nx))+1E-8)
!!$
!!$       dx_adapt(1:nx) = invdiff/maxval(invdiff)
!!$       dx_adapt(1) = dx_adapt(2)
!!$
!!$       dx_ = dx_adapt/sum(dx_adapt)*(xmax+0.5_WP*(x_(nx+1)-x_(nx)))
!!$
!!$       x_ = 0.0_WP
!!$       do i=2,nx+1
!!$          x_(i) = x_(i-1) + dx_(i-1); !min(dx_max,max(dx_min,dx_(i-1)))
!!$       end do
!!$
!!$       !x_ = xmax*x_
!!$
!!$       ! Shift to cell center
!!$       dx = x_(2:nx+1)-x_(1:nx)
!!$
!!$       x = x_(:nx)+0.5_WP*dx
!!$
!!$       if(maxval(x) .gt. 1.0_WP) then
!!$          x = x/maxval(x)
!!$       end if
!!$
!!$    end if
!!$
!!$    peak_xn = x(ind)
!!$
!!$    x = peak_x/peak_xn*x
!!$
!!$    B = 1.0_WP/sqrt(twoPi*sigma2)
!!$
!!$    if (xm > 0.0_WP .and. sigma2 .gt. 0.0_WP) then
!!$       
!!$       f = B*exp(-((log(x)-mu)**2/2.0_WP/sigma2))/x
!!$
!!$       norm = sum(f*dxx)
!!$       if (norm > 0.0_WP) then
!!$          f = f/norm;
!!$       end if
!!$    end if

!!$    if (xm > 0.0_WP .and. x2 > 0.0_WP) then
!!$
!!$       if (ratio > r_min ) then
!!$          ratio = min(r_max,ratio)
!!$          call bisection(q,ratio)
!!$       else
!!$          q = 5.1334_WP
!!$       end if
!!$
!!$       xbar = xm/gamma(1.0_WP+1.0_WP/q)
!!$
!!$       f = (q/xbar**q)*x**(q-1.0_WP)*exp(-(x/xbar)**q)
!!$
!!$       norm = sum(f)
!!$       if (norm > 0.0_WP) then
!!$          f = f/norm;
!!$       end if
!!$
!!$       ! Adapt 
!!$       invdiff = 1.0_WP/abs(f(:nx-1)-f(2:nx))
!!$
!!$       dx_adapt(2:nx) = invdiff/maxval(invdiff)
!!$       dx_adapt(1) = dx_adapt(2)
!!$
!!$       dx_ = dx_adapt/sum(dx_adapt)
!!$
!!$       x_ = 0.0_WP
!!$       do i=2,nx+1
!!$          x_(i) = x_(i-1) + dx_(i-1) 
!!$       end do
!!$
!!$       x_ = xmax*x_
!!$
!!$       dx = x_(2:nx+1)-x_(1:nx)
!!$
!!$       x = x_(:nx)+0.5_WP*dx
!!$
!!$       f = (q/xbar**q)*x**(q-1.0_WP)*exp(-(x/xbar)**q)
!!$
!!$       norm = sum(f*dx)
!!$       if (norm > 0.0_WP) then
!!$          f = f*dx/norm;
!!$       end if
!!$
!!$    end if

 
  end subroutine createDropletGrid

      
  subroutine lognormalPDF(xm,mu,sigma2,x,fy)
    implicit none

    ! ---------------------------------
    real(WP), intent(in) :: xm, mu, sigma2
    real(WP), dimension(:), intent(in) :: x
    real(WP), dimension(:), intent(out) :: fy
    ! ---------------------------------
    real(WP) :: B
    integer :: idx,ind

    B = 1.0_WP/sqrt(twoPi*sigma2)

    if (sigma2 .gt. 0.0_WP) then

       fy = 2.0_WP + B*exp(-((log(x)-mu)**2/2.0_WP/sigma2))/x
    else
       call getIndex(x,xm,idx)
       fy(idx) = 1.0_WP
    end if

    where(x == 0.0_WP) fy = 0.0_WP

  end subroutine lognormalPDF

  subroutine lognormalPDFscalar(mu,sigma2,x,fyy)
    implicit none

    ! ---------------------------------
    real(WP), intent(in) :: mu, sigma2, x
    real(WP), intent(out) :: fyy
    ! ---------------------------------

    real(WP) :: B

    if (sigma2 .gt. 0.0_WP) then

       B = 1.0_WP/sqrt(twoPi*sigma2)

       fyy = 2.0_WP + B*exp(-((log(x)-mu)**2/2.0_WP/sigma2))/x
    else
       fyy = 1.0
    end if

    if(x == 0.0_WP) fyy = 0.0_WP

  end subroutine lognormalPDFscalar
  
  subroutine getIndex(A,val,idx)
    implicit none
    real(WP), dimension(:), intent(in) :: A
    real(WP), intent(in) :: val
    integer, intent(out) :: idx
    ! -------------------------------

    real(WP) :: tmp
    integer :: i, n, flag

    n = size(A)

    do i=1,n-1
       if (val .ge. A(i) .and. val .lt. A(i+1)) then
          idx = i
       end if
    end do

  end subroutine getIndex

  subroutine sortGetIndex(A,val,idx)
    implicit none
    real(WP), dimension(:), intent(inout) :: A
    real(WP), intent(in) :: val
    integer, intent(out) :: idx
    ! -------------------------------

    real(WP) :: tmp
    integer :: i, n, flag

    n = size(A)

    do i=1,n-1
       if (A(i)>A(i+1)) then
          tmp = A(i)
          A(i) = A(i+1)
          A(i+1) = tmp
          idx = i+1
       end if
    end do

  end subroutine sortGetIndex

  subroutine bisection(q_final,ratio)
    implicit none

    real(WP), intent(inout) :: q_final
    real(WP), intent(in) :: ratio
    ! -------------------------------
    real(WP), parameter :: ql0 = 1.0, qh0 = 10.0_WP, tol = 1E-06
    real(WP) :: fs, fl, fh, qs, ql, qh, err
    integer :: iter
    ! -------------------------------

    iter = 0
    err = 1.0_WP
    ql = ql0; qh = qh0; qs = 0.5_WP*(ql+qh)
    do while ( err > tol)
       iter = iter + 1
       fs = gamma(1.0_WP+2.0_WP/qs)/(gamma(1.0_WP+1.0_WP/qs))**2 - ratio
       fl = gamma(1.0_WP+2.0_WP/ql)/(gamma(1.0_WP+1.0_WP/ql))**2 - ratio
       fh = gamma(1.0_WP+2.0_WP/qh)/(gamma(1.0_WP+1.0_WP/qh))**2 - ratio

       if (fl*fs < 0.0_WP) then
          qh = qs; qs = 0.5_WP*(ql+qh)
       elseif (fs*fh < 0.0_WP) then
          ql = qs; qs = 0.5_WP*(ql+qh)
       end if

       err = abs(fs)

    end do
    !write(*,*) 'Bisection iters:', iter
    q_final = qs
  end subroutine bisection

  ! NOT TESTED EXTENSIVELY
  subroutine maximumEntropyFormalism(spray,Dm,D2,D3,k,Dp,pdf)
    implicit none

    ! ---------------------------------

    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(spray%nd), intent(out) :: pdf
    real(WP), intent(in) :: Dm,D2,D3
    integer, intent(in) :: k

    ! ---------------------------------
    real(WP), dimension(spray%nd):: Dp
    real(WP), dimension(spray%nd):: D, g1, g2, g3, g4, p
    real(WP), dimension(4):: lambda, lambdaNew, dlambda, func, pold, pnew, gradf, funcnew
    real(WP), dimension(4,4):: Jac, JacInv, Jacnew, JacInvnew
    real(WP), dimension(2,2):: L12
    real(WP), dimension(2):: gL, coeff
    real(WP) :: L, h, tol, err, relax, relax1, relax2, Eg1, Eg2, Eg3, Eg4, alpha, fold, fnew, gf0, gf1, gfdash0, gfrelax1, gfrelax2
    integer :: i, n, count, max_count

    L = 4.0_WP

    n = spray%nd
    h = L/n

    D = (/ (h*real(i,WP),i=0,n,1) /)

    g1 = (/ (1.0_WP,i=0,n,1) /)
    g2 = D
    g3 = D**2
    g4 = D**3

    tol = spray%solver%nr%tol !1.0E-2_WP
    err = 1.0_WP

    lambda = (/ 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP /)

    p = 3.0_WP*D**2.*exp(-lambda(1)*g1-lambda(2)*g2-lambda(3)*g3-lambda(4)*g4)
    p = p/sum(p*h)

    relax = spray%solver%nr%relax_coeff !1.0E-1_WP
    relax1 = relax 
    relax2 = relax

    Eg1 = 1.0_WP
    Eg2 = Dm/D3**(1.0_WP/3.0_WP)
    Eg3 = D2/D3**(2.0_WP/3.0_WP)
    Eg4 = 1.0_WP

    alpha = spray%solver%nr%alpha !1.0E-04_WP

    count = 0
    max_count = spray%solver%nr%max_count !2000

    do while (abs(err)>tol .and. count<max_count)

       func = (/ sum(p*g1*h)-Eg1, sum(p*g2*h)-Eg2, sum(p*g3*h)-Eg3, sum(p*g4*h)-Eg4 /)
      
       Jac = transpose( reshape( &
             (/ sum(p*g1*(-g1)*h), sum(p*g1*(-g2)*h), sum(p*g1*(-g3)*h), sum(p*g1*(-g4)*h),   &
                sum(p*g2*(-g1)*h), sum(p*g2*(-g2)*h), sum(p*g2*(-g3)*h), sum(p*g2*(-g4)*h),   &
                sum(p*g3*(-g1)*h), sum(p*g3*(-g2)*h), sum(p*g3*(-g3)*h), sum(p*g3*(-g4)*h),   &
                sum(p*g4*(-g1)*h), sum(p*g4*(-g2)*h), sum(p*g4*(-g3)*h), sum(p*g4*(-g4)*h) /),&
             shape(Jac)))

       pold = -matmul(pinv(func),Jac)

       fold = 0.5_WP*sum(func*func)

       gradf = matmul(Jac,func)

       lambdaNew = lambda + relax*pold

       dlambda = lambdaNew-lambda

       err = sqrt(sum(func*func))

       lambda = lambdaNew;

       p = 3.0_WP*D**2*exp(-lambda(1)*g1-lambda(2)*g2-lambda(3)*g3-lambda(4)*g4)
       p = p/sum(p*h)

       funcnew = (/ sum(p*g1*h)-Eg1, sum(p*g2*h)-Eg2, sum(p*g3*h)-Eg3, sum(p*g4*h)-Eg4 /)

       Jacnew = transpose( reshape( &
             (/ sum(p*g1*(-g1)*h), sum(p*g1*(-g2)*h), sum(p*g1*(-g3)*h), sum(p*g1*(-g4)*h),   &
                sum(p*g2*(-g1)*h), sum(p*g2*(-g2)*h), sum(p*g2*(-g3)*h), sum(p*g2*(-g4)*h),   &
                sum(p*g3*(-g1)*h), sum(p*g3*(-g2)*h), sum(p*g3*(-g3)*h), sum(p*g3*(-g4)*h),   &
                sum(p*g4*(-g1)*h), sum(p*g4*(-g2)*h), sum(p*g4*(-g3)*h), sum(p*g4*(-g4)*h) /),&
             shape(Jacnew)))

       pnew = -matmul(pinv(funcnew),Jacnew)

       fnew = 0.5_WP*sum(funcnew*funcnew)

       if(fnew < fold + alpha*sum(gradf*dlambda)) then
             lambdaNew = lambda + relax*pnew
             
             dlambda = lambdaNew - lambda
             err = sqrt(sum(funcnew*funcnew))

             lambda = lambdaNew;

             p = 3.0_WP*D**2*exp(-lambda(1)*g1-lambda(2)*g2-lambda(3)*g3-lambda(4)*g4)
             p = p/sum(p*h)
       end if

       if(fnew > fold + alpha*sum(gradf*dlambda)) then
          relax2 = relax1
          relax1 = relax
          gf0 = fold
          gf1 = fnew
          gfdash0 = sum(gradf*pold)
          gfrelax1 = (gf1-gf0-gfdash0)*relax1**2+gfdash0*relax1+gf0
          gfrelax2 = (gf1-gf0-gfdash0)*relax2**2+gfdash0*relax2+gf0

          if(relax == 1.0_WP) then

             relax = max(0.1_WP,-gfdash0/(2.0_WP*(gf1-gf0-gfdash0)))

          else

             L12 = transpose( reshape( &
                   (/ 1.0_WP/relax1**2, -1.0_WP/relax2**2, -relax2/relax1**2, relax1/relax2**2 /), &
                   shape(L12)))

             gL = (/ gfrelax1-gfdash0*relax1-gf0, gfrelax2-gfdash0*relax2-gf0 /)

             if(relax1 == relax2) relax2 = relax1 + 1.0E-2_WP

             coeff = 1.0_WP/(relax1-relax2)*matmul(L12,gL)
             
             if(coeff(1) == 0.0_WP) coeff(1) = 1.0e-2_WP

             if(coeff(2) == 0.0_WP) coeff(2) = 1.0e-1_WP

             relax = (-coeff(2)+sqrt(max(0.0_WP,coeff(2)**2 - 3.0_WP*coeff(1)*gfdash0)))/3.0_WP/coeff(1)
             relax = max(1.0E-6_WP,max(0.1_WP*relax1,min(0.5_WP*relax1,relax)))
          end if

       end if

       count = count + 1

       if(count == max_count .and. err > 2.0_WP*tol) then

          max_count = min(2*max_count,1000000)
          
       end if

       if(count > max_count) write(*,*) 'Warning! Newton iteration not converged!!! Error: ', err
    end do

    p = 3.0_WP*D**2*exp(-lambda(1)*g1-lambda(2)*g2-lambda(3)*g3-lambda(4)*g4)
    p = p/sum(p*h)

    p = p/D3**(1.0_WP/3.0_WP)

    D = D*D3**(1.0_WP/3.0_WP)

    pdf = p(1:spray%nd)

    Dp = D(1:spray%nd)

    spray%dsdlam(1,k) = lambda(1); spray%dsdlam(2,k) = lambda(2)
    spray%dsdlam(3,k) = lambda(3); spray%dsdlam(4,k) = lambda(4)

    write(*,*) 'Newton itrations: ', count, 'Error: ', err, 'Tol: ', tol
    
  contains
    pure function pinv(A) result(X)
      !! Performs a direct calculation of the inverse of a 4×4 matrix.
      real(WP), intent(in) :: A(4)   !! Vector
      real(WP)             :: X(4)   !! Inverse vector using Moore-Penrose pseudoinverse

      real(WP) :: L2

      L2 = (A(1)**2 + A(2)**2 + A(3)**2 + A(4)**2)

      X = (/ A(1)/L2, A(2)/L2, A(3)/L2, A(4)/L2  /)

    end function pinv
  end subroutine maximumEntropyFormalism

  subroutine computeLiquidFuelProperties(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l => null()
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k
 
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo
    
    do k=kmino,kmaxo

       pc_l => spray%pc_l(k)

       if ( pc_l%name == 'pc' ) then
          call choosePC(pc_l,spray%Fuel)
       end if

       ! Set temperature of liquid and vapor phaseson the grid
       pc_l%T = spray%T_fuel
       pc_l%p = spray%P_a

       ! Compute fuel properties: Liquid phase
       call computeLiqDensity(pc_l)
       call computeLiqViscosity(pc_l)
       call computeLiqThermalConductivity(pc_l)
       call computeLiqHeatCapacity(pc_l)
       call computeVapPressure(pc_l)
       call computeHeatOfVap(pc_l)
       call computeSurfaceTension(pc_l)

       ! Set fuel properties to spray
       spray%C_l(k) = spray%f_C_l*pc_l%liqHeatCapacity
       spray%p_vap(k) = spray%f_Pv*pc_l%vapPressure%val
       spray%L_f(k) = spray%f_Lv*pc_l%HeatOfVap
       spray%T_crit = pc_l%Tcrit
       if(spray%L_f(k) == 0.0_WP) then
          spray%T_sat(k) = pc_l%Tcrit
       else
          !spray%T_sat(k) = min(pc_l%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/pc_l%NormalBoilingPoint - spray%R_gas*log(spray%P_a/101325.0_WP)/pc_l%MolecularWeight*1000.0_WP/spray%L_f(k))))
          spray%T_sat(k) = min(pc_l%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/pc_l%Tcrit - spray%R_gas*log(spray%P_a/pc_l%Pcrit)/pc_l%MolecularWeight*1000.0_WP/spray%L_f(k))))
       end if

       nullify(pc_l)
    end do
    
    call computeBoilingTemperature(spray%satPT_file,spray%P_a,spray%T_boil)

    pc_l => spray%pc_l(1)

    ! Set fuel properties to spray (constant throughout spray simulation)
    spray%sigma = spray%f_sigma*pc_l%SurfaceTension%val
    spray%rho_l = spray%f_rho_l*pc_l%liqDensity
    spray%visc_l = spray%f_mu_l*pc_l%liqViscosity%val
    spray%lambda_l = spray%f_lambda_l*pc_l%liqThermalConductivity%val
    spray%MW_f = pc_l%MolecularWeight/1000.0_WP
    spray%MP = pc_l%MeltingPoint
    spray%NBP = pc_l%NormalBoilingPoint

    spray%sigma_loc = spray%f_mu_l*spray%sigma
    spray%rho_l_loc = spray%f_rho_l*spray%rho_l
    spray%visc_l_loc = spray%f_mu_l*spray%visc_l
    spray%lambda_l_loc = spray%f_lambda_l*spray%lambda_l

  end subroutine computeLiquidFuelProperties

subroutine computeSatLiquidFuelProperties(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l => null()

       pc_l => spray%pc_l(1)

       if ( pc_l%name == 'pc' ) then
          call choosePC(pc_l,spray%Fuel)
       end if

       ! Set temperature of liquid and vapor phaseson the grid
       pc_l%T = spray%T_boil
       pc_l%p = spray%P_a

       ! Compute fuel properties: Liquid phase
       call computeLiqDensity(pc_l)
       call computeLiqViscosity(pc_l)
       call computeLiqThermalConductivity(pc_l)
       call computeLiqHeatCapacity(pc_l)
       call computeVapPressure(pc_l)
       call computeHeatOfVap(pc_l)
       call computeSurfaceTension(pc_l)

       ! Set fuel properties to spray
       spray%C_lb = pc_l%liqHeatCapacity
       spray%L_fb = pc_l%HeatOfVap

       nullify(pc_l)
 
    
    pc_l => spray%pc_l(1)

    ! Set fuel properties to spray (constant throughout spray simulation)
    spray%sigmab = pc_l%SurfaceTension%val
    spray%rho_lb = pc_l%liqDensity
    spray%visc_lb = pc_l%liqViscosity%val
    spray%lambda_lb = pc_l%liqThermalConductivity%val

  end subroutine computeSatLiquidFuelProperties

  subroutine computeSatBubbleProperties(spray,Tbub)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l => null()
    real(WP) :: Tbub

       pc_l => spray%pc_l(1)

       if ( pc_l%name == 'pc' ) then
          call choosePC(pc_l,spray%Fuel)
       end if

       ! Set temperature of liquid and vapor phaseson the grid
       pc_l%T = Tbub
       pc_l%p = spray%P_a

       ! Compute fuel properties: Liquid phase
       call computeVapPressure(pc_l)
       call computeHeatOfVap(pc_l)
       call computeLiqDensity(pc_l)
       call computeSurfaceTension(pc_l)

       ! Set fuel properties to spray
       spray%L_fbub = pc_l%HeatOfVap
       spray%P_bub = spray%f_Pv*pc_l%vapPressure%val
       spray%sigma_bub = pc_l%SurfaceTension%val
       spray%rho_lbub = pc_l%liqDensity
       nullify(pc_l)

    pc_l => spray%pc_l(1)

  end subroutine computeSatBubbleProperties

  subroutine getVaporFuelTemperature(spray,Tl)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: Tl

    ! ---------------------------------
    type(pc_t), pointer :: pc => null()
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, iter
    real(WP) :: dT, T, CpL, CpLold, CpL_integral, CpG, CpGold, CpG_integral
    
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    call allocate_pc(pc)
    if ( pc%name == 'pc' ) then
       call choosePC(pc,spray%Fuel)
    end if

    do k=kmino,kmaxo

       if(Tl(k) == 0.0_WP) then
          cycle
       end if

       ! Set temperature of liquid and vapor phases on the grid
       pc%T = Tl(k)
       pc%p = spray%P_a

       call computeHeatOfVap(pc)
       CpL_integral = spray%f_Lv*pc%HeatOfVap

       call computeLiqHeatCapacity(pc)

       ! Subtracting heat of vaporization from the liquid enthalpy
       CpG_integral = 0.0_WP
       dT = 0.1_WP
       T = Tl(k)  
       CpGold = spray%f_C_l*pc%liqHeatCapacity
       T = T + dT
       iter = 0
       do while(CpG_integral<CpL_integral)
          pc%T = T
          if(T<300.0_WP) pc%T = 300.0_WP
          call computeIG_HeatCapacity(pc)
          CpG = spray%f_Cp_v*pc%IG_HeatCapacity
          CpG_integral = CpG_integral + 0.5_WP*(CpG+CpGold)*dT
          T = T + dT
          if(CpG_integral>CpL_integral) then
             dT = 2.0_WP*(CpG_integral-CpL_integral)/(CpG+CpGold)
             T = T + dT
             spray%Tv(k) = max(10.0_WP,T)
             !write(*,*) 'Fuel vapor temp iter:',iter,spray%Tv(k)
             exit
          end if
          CpGold = CpG
          iter = iter + 1
       end do
    end do
    call deallocate_pc(pc)
  end subroutine getVaporFuelTemperature

  subroutine getLiquidFuelPropertiesFromLFPT(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(fpt_t), pointer :: fpt
    real(WP), dimension(:), pointer :: values
    real(WP) :: T, dummy, T_sat
    integer :: col
    
    fpt => spray%fpt

    allocate(spray%LFPT(size(fpt%lTemp,1),8))
    spray%LFPT(:,1) = fpt%lTemp
    spray%LFPT(:,2) = fpt%liqViscosity
    spray%LFPT(:,3) = fpt%SurfaceTension
    spray%LFPT(:,4) = fpt%HeatOfVap
    spray%LFPT(:,5) = fpt%vapPressure
    spray%LFPT(:,6) = fpt%liqThermalConductivity
    spray%LFPT(:,7) = fpt%liqDensity
    spray%LFPT(:,8) = fpt%liqHeatCapacity

    spray%NBP = fpt%NormalBoilingPoint
    spray%MP = fpt%MeltingPoint
    spray%MW_f = fpt%MolecularWeight

    allocate(values(size(spray%LFPT,2))); values = 0.0_WP

    ! Set temperature of liquid and vapor phases on the grid
    T = max(minval(spray%LFPT(:,1)),min(spray%T_fuel,maxval(spray%LFPT(:,1))))

    ! Linearly interpolate fuel properties: Liquid phase
    ! Use temperature column as a x-value for interpolation
    col = 1
    call interpolate1(spray%LFPT,col,T,values)

    ! Set fuel properties to spray
    spray%C_l = spray%f_C_l*values(8)
    spray%L_f = spray%f_Lv*values(4)

    if(values(5) < 0.0_WP) then
       spray%p_vap = spray%f_Pv*(fpt%pcrit)*exp((spray%L_f*spray%MW_f/spray%R_gas)*(1.0_WP/fpt%Tcrit-1.0_WP/T))
       if(T >= fpt%Tcrit) then
          spray%T_sat = fpt%Tcrit
       else
          spray%T_sat = min(fpt%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/fpt%Tcrit - spray%R_gas*log(spray%P_a/fpt%pcrit)/fpt%MolecularWeight/spray%L_f)))
       end if
    else
       spray%p_vap = spray%f_Pv*values(5)
       if(spray%P_a >= fpt%pcrit) then
          spray%T_sat = fpt%Tcrit
       else
          call interpolate1val(fpt%vapPressure,fpt%lTemp,spray%P_a,T_sat,dummy)
          spray%T_sat = T_sat
       end if
    end if

    ! Set fuel properties to spray (constant throughout spray simulation)
    spray%sigma = spray%f_sigma*values(3)
    spray%rho_l = spray%f_rho_l*values(7)
    spray%visc_l = spray%f_mu_l*values(2)
    spray%lambda_l = spray%f_lambda_l*values(6)

    spray%sigma_loc = spray%f_mu_l*spray%sigma
    spray%rho_l_loc = spray%f_rho_l*spray%rho_l
    spray%visc_l_loc = spray%f_mu_l*spray%visc_l
    spray%lambda_l_loc = spray%f_lambda_l*spray%lambda_l

    if (spray%MW_f == -9999.0_WP) then
       write(*,*) 'Error: Molecular weight of the fuel not provided. Please add following line to input file:'
       write(*,*) 'Fuel molecular weight : <value>'
       call abort
    end if

    if (spray%NBP == -9999.0_WP) then
       write(*,*) 'Error: Normal Boiling Point of the fuel not provided. Please add following line to input file:'
       write(*,*) 'Normal boiling point : <value>'
       call abort
    end if

    if (spray%MP == -9999.0_WP) then
       write(*,*) 'Error: Melting point of the fuel not provided. You may add following line to input file:'
       write(*,*) 'Melting point : <value>'

       spray%MP = 263.0_WP

       write(*,*) 'Setting melting point of the fuel to ', spray%MP, '...'

    end if

    deallocate(values)
  
  end subroutine getLiquidFuelPropertiesFromLFPT

  ! Linear Interpolation routine
  subroutine interpolate1(Table,col,X,Y)
    implicit none
    ! ---------------------------------
    real(WP), dimension(:,:), intent(in) :: Table
    integer, intent(in) :: col
    real(WP), intent(in) :: X
    real(WP), dimension(:), intent(out) :: Y
    ! ---------------------------------
    integer, dimension(2) :: Tshape
    integer :: i, l, h

    Tshape = shape(Table)

    do i=1,Tshape(1)-1
       if (X >= Table(i,col) .and. X <= Table(i+1,col)) then
          Y = (Table(i+1,:)-Table(i,:))/(Table(i+1,col)-Table(i,col))*(X-Table(i,col)) + Table(i,:)
       end if
    end do

  end subroutine interpolate1

  ! Linear Interpolation routine single value
  subroutine interpolate1val(X,Y,xval,yval,dxval)
    implicit none
    ! ---------------------------------
    real(WP), dimension(:), intent(in) :: X, Y
    real(WP), intent(inout) :: xval, yval, dxval
    ! ---------------------------------
    integer :: n, i, l, h

    n = size(X,1)

    do i=1,n-1
       if (xval >= X(i) .and. xval <= X(i+1)) then
          yval = (Y(i+1)-Y(i))/(X(i+1)-X(i))*(xval-X(i)) + Y(i)
          dxval = (X(i+1)-X(i))
       end if
    end do

  end subroutine interpolate1val

  subroutine updateLiquidFuelProperties(spray,T)
    implicit none
    real(WP) :: eps = 1.0E-16_WP
    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    type(pc_t), pointer :: pc_l => null()
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k
 
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    do k=kmino,kmaxo
     !if ( spray%d2(k) > eps*eps .and. spray%dm(k) > eps ) then
       pc_l => spray%pc_l(k)

       if ( pc_l%name == 'pc' ) then
          call choosePC(pc_l,spray%Fuel)
       end if

       ! Set temperature of liquid and vapor phaseson the grid
       pc_l%T = T(k)
       pc_l%p = spray%P_a

       ! Compute fuel properties: Liquid phase
       call computeLiqDensity(pc_l)
       call computeLiqViscosity(pc_l)
       call computeLiqThermalConductivity(pc_l)
       call computeSurfaceTension(pc_l)

       call computeLiqHeatCapacity(pc_l)
       call computeVapPressure(pc_l)
       call computeHeatOfVap(pc_l)

       ! Set fuel properties to spray
       spray%C_l(k) = spray%f_C_l*pc_l%liqHeatCapacity
       spray%p_vap(k) = spray%f_Pv*pc_l%vapPressure%val
       spray%L_f(k) = spray%f_Lv*pc_l%HeatOfVap
       if(spray%L_f(k) == 0.0_WP) then
          spray%T_sat(k) = pc_l%Tcrit
       else
          !spray%T_sat(k) = min(pc_l%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/pc_l%NormalBoilingPoint - spray%R_gas*log(spray%P_a/101325.0_WP)/pc_l%MolecularWeight*1000.0_WP/spray%L_f(k))))
          spray%T_sat(k) = min(pc_l%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/pc_l%Tcrit - spray%R_gas*log(spray%P_a/pc_l%Pcrit)/pc_l%MolecularWeight*1000.0_WP/spray%L_f(k))))
       end if

       spray%rho_l_loc(k) = spray%f_rho_l*pc_l%liqDensity
       spray%visc_l_loc(k) = spray%f_mu_l*pc_l%liqViscosity%val
       spray%lambda_l_loc(k) = spray%f_lambda_l*pc_l%liqThermalConductivity%val
       spray%sigma_loc(k) = spray%f_sigma*pc_l%SurfaceTension%val

       nullify(pc_l)
     !end if
    end do

    ! XXXXXXXXXXXXXXX TESTING 
!!$    spray%sigma = pc_l%SurfaceTension%val
!!$    spray%rho_l = pc_l%liqDensity
!!$    spray%visc_l = pc_l%liqViscosity%val
 
  end subroutine updateLiquidFuelProperties

  subroutine updateLiquidFuelPropertiesFromLFPT(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    type(fpt_t), pointer :: fpt
    real(WP), dimension(:), pointer :: values
    real(WP) :: dummy
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: col, k
 
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    fpt => spray%fpt

    allocate(values(size(spray%LFPT,2))); values = 0.0_WP
    
    do k=kmino,kmaxo

       ! Linearly interpolate fuel properties: Liquid phase
       ! Use temperature column as a x-value for interpolation
       col = 1
       call interpolate1(spray%LFPT,col,T(k),values)

       ! Set fuel properties to spray
       spray%C_l(k) = spray%f_C_l*values(8)
       spray%L_f(k) = spray%f_Lv*values(4)

       if(values(5) < 0.0_WP) then
          spray%p_vap(k) = spray%f_Pv*(fpt%pcrit)*exp((spray%L_f(k)*spray%MW_f/spray%R_gas)*(1.0_WP/fpt%Tcrit-1.0_WP/T(k)))
          if(T(k) >= fpt%Tcrit) then
             spray%T_sat(k) = fpt%Tcrit
          else
             spray%T_sat(k) = min(fpt%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/fpt%Tcrit - spray%R_gas*log(spray%P_a/fpt%pcrit)/fpt%MolecularWeight/spray%L_f(k))))
          end if
       else
          spray%p_vap(k) = spray%f_Pv*values(5)
          if(spray%P_a >= fpt%pcrit) then
             spray%T_sat(k) = fpt%Tcrit
          else
             call interpolate1val(fpt%vapPressure,fpt%lTemp,spray%P_a,spray%T_sat(k),dummy)
          end if
       end if
    end do
  
    deallocate(values)
  end subroutine updateLiquidFuelPropertiesFromLFPT

  ! Ambient properties are fixed to air properties
  ! TODO: Different composition of ambient gas
  subroutine computeAmbientProperties(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), pointer :: T_a=>null()

    T_a => spray%T_a

    spray%rho_a = spray%P_a*spray%MW_a/(spray%Z_a*spray%R_gas*T_a)

    spray%visc_a = - 0.000000539820785e-5_WP*T_a**2 &
                   + 0.003926564011368e-5_WP*T_a + 0.761680769346030e-5_WP

    spray%lambda_a = -0.000000006622335_WP*T_a**2 + 0.000065566840252_WP*T_a + 0.006983319214040_WP

    spray%Cp_a = 1.9327E-10_WP*T_a**4 - 7.9999E-07_WP*T_a**3 + 1.1407E-03_WP*T_a**2 &
               - 4.4890E-01_WP*T_a + 1.0575E+03_WP

  end subroutine computeAmbientProperties

  subroutine computeVaporFuelProperties(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), intent(in) :: T

    ! ---------------------------------
    type(pc_t), pointer :: pc_v => null()

    pc_v => spray%pc_v(spray%kmaxo)

    if ( pc_v%name == 'pc' ) then
       call choosePC(pc_v,spray%Fuel)
    end if

    ! Set temperature and pressure of vapor phases on the grid      
    pc_v%T = T
    pc_v%p = spray%P_a

    ! Compute fuel properties: Vapor phase
    call computeIG_vapDensity(pc_v)
    call computeVapViscosity(pc_v)
    call computeVapThermalConductivity(pc_v)
    call computeIG_HeatCapacity(pc_v)
    
    ! Set fuel vapor properties in spray object
    spray%rho_v = spray%f_rho_v*pc_v%IG_vapDensity
    spray%visc_v = spray%f_mu_v*pc_v%vapViscosity%val
    spray%lambda_v = spray%f_lambda_v*pc_v%vapThermalConductivity%val
    spray%Cp_v = spray%f_Cp_v*pc_v%IG_HeatCapacity

  end subroutine computeVaporFuelProperties

  subroutine getVaporFuelPropertiesFromVFPT(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), intent(in):: T

    ! ---------------------------------
    type(pc_t), pointer :: pc_v => null()
    type(fpt_t), pointer :: fpt => null()
    real(WP), dimension(:), pointer :: values
    integer :: col
    real(WP) :: R_gas

!!$    pc_v => spray%pc_v(spray%kmaxo)
!!$    if ( pc_v%name == 'pc' ) then
!!$       call choosePC(pc_v,spray%Fuel)
!!$    end if
!!$
!!$    ! Set temperature and pressure of vapor phases on the grid      
!!$    pc_v%T = T
!!$    pc_v%p = spray%P_a

    fpt => spray%fpt

    if(.not. associated(spray%VFPT)) then
       allocate(spray%VFPT(size(fpt%vTemp,1),6))
       spray%VFPT(:,1) = fpt%vTemp
       spray%VFPT(:,2) = fpt%vapViscosity
       spray%VFPT(:,3) = fpt%vapThermalConductivity
       spray%VFPT(:,4) = fpt%vapDensity
       spray%VFPT(:,5) = fpt%vapHeatCapacity
       spray%VFPT(:,6) = fpt%DiffusionCoefficient
    end if

    allocate(values(size(spray%VFPT,2))); values = 0.0_WP
        
    if (spray%MW_f == -9999.0_WP) then
       write(*,*) 'Error: Molecular weight of the fuel not provide. Please add following line to input file:'
       write(*,*) 'Fuel molecular weight : <value>'
       call abort
    end if

    R_gas = spray%R_gas   ! J/K/mol

    ! Linearly interpolate fuel properties: vapor phase
    ! Use temperature column as a x-value for interpolation
    col = 1
    call interpolate1(spray%VFPT,col,T,values)

    if(values(4) < 0.0_WP) then
       spray%rho_v = spray%f_rho_v*spray%P_a*spray%MW_f/(R_gas*T)
    else
       spray%rho_v = spray%f_rho_v*values(4)
    end if
    spray%visc_v = spray%f_mu_v*values(2)
    spray%lambda_v = spray%f_lambda_v*values(3)
    spray%Cp_v = spray%f_Cp_v*values(5)

    deallocate(values)
 
  end subroutine getVaporFuelPropertiesFromVFPT

  ! Remark:
  ! Choice of calculation of reference temperature 
  ! has little influence on liquid penetration
  subroutine computeRefTemperature(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k

    ! Reference temperature for evaporation model
    spray%T_ref = spray%Tg*spray%T_fuel 
    ! Warning: This approach results in jumps in gas and reference temperature, 
    ! which cause jumps in gaseous properties.
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo
   
   if(spray%flash_boiling)then
   do k=kmino,kmaxo
    if(spray%Td(k)*spray%T_fuel.gt.spray%T_boil)then
      spray%T_ref(k) = (spray%Tg(k) + 2.0_WP*spray%T_boil/spray%T_fuel)/3.0_WP*spray%T_fuel
    else 
      if(spray%Td(k) > 0.0_WP)then
        spray%T_ref(k) = (spray%Tg(k) + 2.0_WP*spray%Td(k))/3.0_WP*spray%T_fuel
      end if
    end if
   end do
   else 
      where(spray%Td > 0.0_WP) spray%T_ref = (spray%Tg + 2.0_WP*spray%Td)/3.0_WP*spray%T_fuel
   end if
    !spray%T_ref = spray%T_fuel
  end subroutine computeRefTemperature

  subroutine computeRefAmbientProperties(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    real(WP), dimension(spray%nzo) :: denom

    denom = spray%Z_a*spray%R_gas*T    

    spray%rho_ra = spray%P_a*spray%MW_a/denom

    spray%visc_ra = - 0.000000539820785e-5_WP*T**2 + 0.003926564011368e-5_WP*T + 0.761680769346030e-5_WP

    spray%lambda_ra = -0.000000006622335_WP*T**2 + 0.000065566840252_WP*T + 0.006983319214040_WP

    spray%Cp_ra = 1.9327E-10_WP*T**4 - 7.9999E-07_WP*T**3 + 1.1407E-03_WP*T**2 &
               - 4.4890E-01_WP*T + 1.0575E+03_WP

  end subroutine computeRefAmbientProperties

  subroutine computeRefVaporFuelProperties(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    type(pc_t), pointer :: pc_v => null()
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k
 
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo
   
    do k=kmino,kmaxo
       pc_v => spray%pc_v(k)

       if ( pc_v%name == 'pc' ) then
          call choosePC(pc_v,spray%Fuel)
       end if

       ! Set temperature and pressure of vapor phases on the grid      
       pc_v%T = T(k)
       pc_v%p = spray%P_a

       ! Compute fuel properties: Vapor phase
       call computeIG_vapDensity(pc_v)
       call computeVapViscosity(pc_v)
       call computeVapThermalConductivity(pc_v)
       call computeIG_HeatCapacity(pc_v)
       !call computeDiffusionCoefficientFuller(pc_v)
       call computeDiffusionCoefficientWilkeLee(pc_v)
       !call computeDiffusionCoefficientHighPres(pc_v,spray%Y_a(k),spray%Y_v(k))
       !call computeDiffusionCoefficientTheory(pc_v)

       ! Set fuel vapor properties in spray object
       spray%rho_rv(k) = spray%f_rho_v*pc_v%IG_vapDensity
       spray%visc_rv(k) = spray%f_mu_v*pc_v%vapViscosity%val
       spray%lambda_rv(k) = spray%f_lambda_v*pc_v%vapThermalConductivity%val
       spray%Cp_rv(k) = spray%f_Cp_v*pc_v%IG_HeatCapacity
       !spray%G_rv(k) = spray%f_G_v*pc_v%DiffusionCoefficientFuller
       spray%G_rv(k) = spray%f_G_v*pc_v%DiffusionCoefficientWilkeLee
       !spray%G_rv(k) = spray%f_G_v*pc_v%DiffusionCoefficientHighPres
       !spray%G_rv(k) = spray%f_G_v*pc_v%DiffusionCoefficientTheory
       !spray%G_rv(k) = spray%lambda_rv(k)/spray%rho_rv(k)/spray%Cp_rv(k)
       !spray%G_rv(k) = 1.293_WP*4.16e-06_WP/spray%rho_rv(k)*(T(k)/273.0_WP)**(1.6_WP-1.0_WP)
       !spray%G_rv(k) = 1.293_WP*2.49e-05_WP/spray%rho_rv(k)*(T(k)/273.0_WP)**(1.8_WP-1.0_WP)

       nullify(pc_v)
    end do
  end subroutine computeRefVaporFuelProperties

  subroutine getRefVaporFuelPropertiesFromVFPT(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in):: T

    ! ---------------------------------
    type(fpt_t), pointer :: fpt => null()
    real(WP), dimension(:), pointer :: values
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: col, k
    real(WP) :: R_gas

    fpt => spray%fpt

    if(.not. associated(spray%VFPT)) then
       allocate(spray%VFPT(size(fpt%vTemp,1),6))
       spray%VFPT(:,1) = fpt%vTemp
       spray%VFPT(:,2) = fpt%vapViscosity
       spray%VFPT(:,3) = fpt%vapThermalConductivity
       spray%VFPT(:,4) = fpt%vapDensity
       spray%VFPT(:,5) = fpt%vapHeatCapacity
       spray%VFPT(:,6) = fpt%DiffusionCoefficient
       fpt%DiffusionCoefficient = 0.0_WP
    end if
    
    allocate(values(size(spray%VFPT,2))); values = 0.0_WP
    
    if (spray%MW_f == -9999.0_WP) then
       write(*,*) 'Error: Molecular weight of the fuel not provide. Please add following line to the input file:'
       write(*,*) 'Fuel molecular weight : <value>'
       call abort
    end if

    R_gas = spray%R_gas   ! J/K/kmol
    
    allocate(values(size(spray%VFPT,2))); values = 0.0_WP

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo
   
    do k=kmino,kmaxo
       ! Linearly interpolate fuel properties: Liquid phase
       ! Use temperature column as a x-value for interpolation
       col = 1
       call interpolate1(spray%VFPT,col,T(k),values)

       if(values(4) < 0.0_WP) then
          spray%rho_rv(k) = spray%f_rho_v*spray%P_a*spray%MW_f/(R_gas*T(k))
       else
          spray%rho_rv(k) = spray%f_rho_v*values(4)
       end if

       spray%visc_rv(k) = spray%f_mu_v*values(2)
       spray%lambda_rv(k) = spray%f_lambda_v*values(3)
       spray%Cp_rv(k) = spray%f_Cp_v*values(5)
    
       if(values(6) < 0.0_WP) then
          !call computeDiffusionCoefficientFullerFPT(fpt,spray%P_a,T(k),spray%G_rv(k))
          call computeDiffusionCoefficientWilkeLeeFPT(fpt,spray%P_a,T(k),spray%G_rv(k))
          !call computeDiffusionCoefficientTheoryFPT(fpt,spray%P_a,T(k),spray%G_rv(k))
          spray%G_rv(k) = spray%f_G_v*spray%G_rv(k)
          !spray%G_rv(k) = 1.293_WP*4.16e-06_WP/spray%rho_rv(k)*(T(k)/273.0_WP)**(1.6_WP-1.0_WP)
       else
          spray%G_rv(k) = spray%f_G_v*values(6)/spray%rho_rv(k)
          !spray%G_rv(k) = 1.293_WP*4.16e-06_WP/spray%rho_rv(k)*(T(k)/273.0_WP)**(1.6_WP-1.0_WP)
       end if
       
    end do
    
    deallocate(values)

  end subroutine getRefVaporFuelPropertiesFromVFPT

  subroutine computeGasMixtureProperties(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    real(WP), dimension(:), pointer :: Y_v=>null(), Y_a=>null()
    real(WP), dimension(spray%nzo) :: phi_va, phi_av

    Y_v => spray%Y_v; Y_a => spray%Y_a

    if (associated(spray%VFPT)) then
       call getRefVaporFuelPropertiesFromVFPT(spray,T)
    else
       call computeRefVaporFuelProperties(spray,T)
    end if

    call computeRefAmbientProperties(spray,T)

    spray%rho_g = (Y_a + Y_v)/(Y_a/spray%rho_ra + Y_v/spray%rho_rv)
    spray%rho_g(spray%kmino:spray%kmin-1) = MIN_REAL_WP !spray%rho_g(spray%kmin)

    ! Wilke's formula for mixture viscosity
    phi_va = ((1.0_WP+sqrt(spray%visc_rv/spray%visc_ra)*(spray%WR)**0.25_WP)**2) &
                     /sqrt(8.0_WP*(1.0_WP+1.0_WP/spray%WR))
    phi_av = ((1.0_WP+sqrt(spray%visc_ra/spray%visc_rv)*(1.0_WP/spray%WR)**0.25_WP)**2) &
                     /sqrt(8.0_WP*(1.0_WP+spray%WR))

    spray%visc_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%visc_rv &
                 + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%visc_ra
    spray%visc_g(spray%kmino:spray%kmin-1) = MIN_REAL_WP !spray%visc_g(spray%kmin)

    ! The same mixture rule for viscosity can be used for thermal conductivity: Poling 2001
    spray%lambda_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%lambda_rv &
                   + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%lambda_ra
    spray%lambda_g(spray%kmino:spray%kmin-1) = MIN_REAL_WP !spray%lambda_g(spray%kmin)

    spray%Cp_g = (Y_v*spray%Cp_rv + Y_a*spray%Cp_ra)/(Y_v+Y_a)
    spray%Cp_g(spray%kmino:spray%kmin-1) = MIN_REAL_WP !spray%Cp_g(spray%kmin)

    !phi_av = ((1.0_WP+(1.0_WP/spray%WR)**0.25_WP)**2) &
    !                 /sqrt(8.0_WP*(1.0_WP+spray%WR))
    
    !phi_va = (1.0_WP-Y_v/spray%Y_g)*spray%G_rv
    !spray%Sc_g = spray%visc_g/(spray%rho_g*phi_va)

    spray%Sc_g = spray%visc_g/(spray%rho_g*spray%G_rv)
    spray%Pr_g = spray%visc_g*spray%Cp_g/spray%lambda_g

  end subroutine computeGasMixtureProperties

  subroutine computeGasMixtureProperties_fixed_DR_VR(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    real(WP), dimension(:), pointer :: Y_v=>null(), Y_a=>null()
    real(WP), dimension(spray%nzo) :: phi_va, phi_av

    Y_v => spray%Y_v; Y_a => spray%Y_a

    spray%rho_g = (Y_a + Y_v)/(Y_a*spray%DRa + Y_v*spray%DRv)*spray%rho_l
    spray%rho_g(spray%kmino:spray%kmin-1) = spray%rho_g(spray%kmin)

    phi_va = ((1.0_WP+sqrt(spray%VRa/spray%VRv)*(spray%WR)**0.25_WP)**2) &
                     /sqrt(8.0_WP*(1.0_WP+1.0_WP/spray%WR))
    phi_av = ((1.0_WP+sqrt(spray%VRv/spray%VRa)*(1.0_WP/spray%WR)**0.25_WP)**2) &
                     /sqrt(8.0_WP*(1.0_WP+spray%WR))

    spray%visc_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%visc_l/spray%VRv + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%visc_l/spray%VRa
    spray%visc_g(spray%kmino:spray%kmin-1) = spray%visc_g(spray%kmin)

    spray%lambda_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%lambda_rv + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%lambda_ra
    spray%lambda_g(spray%kmino:spray%kmin-1) = spray%lambda_g(spray%kmin)

    spray%Cp_g = (Y_v*spray%Cp_rv + Y_a*spray%Cp_ra)/(Y_v+Y_a)
    spray%Cp_g(spray%kmino:spray%kmin-1) = spray%Cp_g(spray%kmin)

    phi_va = (1.0_WP-Y_v/spray%Y_g)*spray%G_rv

    spray%Sc_g = spray%visc_g/(spray%rho_g*phi_va)
    spray%Pr_g = spray%visc_g*spray%Cp_g/spray%lambda_g

  end subroutine computeGasMixtureProperties_fixed_DR_VR

  subroutine entrainmentTerm(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:), pointer :: u_g=>null(), omega_ent=>null()
    real(WP), pointer :: DRa=>null()

    if(.not. spray%ent_model) return
    
    u_g => spray%u_g; omega_ent => spray%omega_ent;  DRa => spray%DRa

    omega_ent(spray%kmin:spray%kmax) = spray%beta*u_g(spray%kmin:spray%kmax)/DRa*spray%b(spray%kmin:spray%kmax)
    !omega_ent(spray%kmin:spray%kmax) = 7.0_WP/spray%Re*(1.0_WP/(spray%VRg(spray%kmin:spray%kmax))+spray%VRtg(spray%kmin:spray%kmax))
    !omega_ent(spray%kmin:spray%kmax) = spray%c_k/spray%Re*(1.0_WP/spray%VRtg(spray%kmin:spray%kmax))
    
    !omega_ent(spray%kmin:spray%kmax) = 4.0_WP*spray%beta/DRa*spray%z(spray%kmin:spray%kmax)/sqrt((1.0_WP/DRa-1.0_WP)**2+4.0_WP*spray%b(spray%kmin:spray%kmax)**2/DRa)

  end subroutine entrainmentTerm

subroutine evaporationModel(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: eps = 1E-16_WP
    real(WP), dimension(:,:), pointer :: di=>null(), dsd=>null(), Red=>null(), Shd=>null(), Nud=>null(), f2=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_ref=>null(), &
                                       h=>null(), b=>null(), &
                                       DRg=>null(), VRg=>null(), Sc_g=>null(), Pr_g=>null(), &
                                       Td=>null(), Tg=>null(), CR=>null(), LR=>null()
    real(WP), pointer :: Re=>null(), De=>null()
    real(WP), pointer :: viscl=>null(), c_l=>null(), lambda_l=>null(), Lv_b=>null(), Tfuel=>null(), Tboil=>null()
    real(WP), dimension(spray%nzo) :: Lv, Lk, Xeq, T_d, Ystar_fe, Bdeq, zeta, norm_var, tmp
    real(WP), dimension(spray%nd) :: zetta, Xneq, Ystar_f, Yref, Bd, coeff, K_vap_ext, K_vap_ht1, K_vap_ht2, K_vap_swl, Nu_int, Qd, Vd, dst_fl
    real(WP), dimension(spray%nd,spray%nzo) :: D
    real(WP) :: coeff_a, coeff_full, TBd, K_T, x, sigma
    real(WP) :: HtCoef, Pr_l, HRg1, HRg2, f3, deltaT
    integer :: k
    
    if(.not. spray%evap_model) return

    spray%omega_vap = 0.0_WP; spray%omega_vapdm = 0.0_WP; spray%omega_vapd2 = 0.0_WP; 
    spray%omega_vapd3 = 0.0_WP
    spray%omega_exp = 0.0_WP; spray%omega_expdm = 0.0_WP; spray%omega_expd2 = 0.0_WP; 
    spray%omega_expd3 = 0.0_WP
    spray%omega_T = 0.0_WP; spray%Shd = 0.0_WP; spray%Nud = 0.0_WP; 

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_ref => spray%Y_ref; Td => spray%Td; Tg => spray%Tg
    viscl => spray%visc_l; c_l => spray%C_l(1); lambda_l => spray%lambda_l; Lv_b => spray%L_fb; Tboil => spray%T_boil; Tfuel => spray%T_fuel 
    di => spray%di; b => spray%b; h => spray%h; dsd => spray%dsd; Red => spray%Red; Shd => spray%Shd; Nud => spray%Nud
    DRg => spray%DRg; VRg => spray%VRg; Sc_g => spray%Sc_g; Pr_g => spray%Pr_g; f2 => spray%f2
    Re => spray%Re; CR => spray%CR; LR => spray%LR; De => spray%De

    Lv = spray%L_f; TBd = spray%NBP; T_d = spray%T_fuel*spray%Td; D = spray%D_eff*di

   do k = spray%kmin,spray%kmax

       if ( spray%d2(k) > eps*eps .and. spray%dm(k) > eps ) then

        if (spray%flash_boiling .and. T_d(k).gt.Tboil) then
          Pr_l = viscl*c_l/lambda_l  
          HRg1 = (c_l*Tfuel)/Lv_b 
          HRg2 = (c_l*Tfuel)/Lv(k)
          deltaT = T_d(k) - Tboil

          ! Calculation of heat transfer coefficient suggested by Adachi et al.
          if ((deltaT .gt. 0.0) .and. (deltaT .lt. 5.0)) then
             HtCoef = 760.0_WP*(deltaT)**0.26_WP
          else if ((deltaT .gt. 5.0) .and. (deltaT .lt. 25.0)) then
             HtCoef = 27.0_WP*(deltaT)**2.33_WP
          else 
             HtCoef = 13800.0_WP*(deltaT)**0.39_WP
          end if

          Nu_int = HtCoef*D(:,k)/lambda_l
          K_vap_ht1 = 4.0_WP*Nu_int*HRg1*(Td(k)-(Tboil/Tfuel))/Pr_l
      
          Nud(:,k) = (2.0_WP + 0.552_WP*(Red(:,k)**0.5_WP)*(Pr_g(k)**(1.0_WP/3.0_WP)))

          if(T_d(k).le.Tboil+5.0_WP)then       
            call computeSatBubbleProperties(spray,Tboil+5.0_WP)           !(spray%P_a/spray%P_bub.gt.0.97_WP)then
            f3 = 1.0_WP/(1.0_WP-spray%P_a/spray%P_bub)
          else 
            call computeSatBubbleProperties(spray,T_d(k)) 
            f3 = 1.0_WP/(1.0_WP-spray%P_a/spray%P_bub)
          end if

          K_vap_ht2 = 4.0_WP*Nud(:,k)*HRg1*f3*abs(Tboil/Tfuel-Tg(k))/Pr_l

          K_vap_ext = K_vap_ht1+K_vap_ht2

          K_vap_swl = (1.0_WP/spray%DRv)*4.0_WP*di(:,k)*spray%Nbub*(spray%Rdot/spray%U_inj)  
        
          zetta = (0.5_WP*Pr_g(k)*VRg(k)*HRg1/Pr_l)*(Nu_int*(Td(k)-Tboil/Tfuel)+f3*Nud(:,k)*abs(Tboil/Tfuel-Tg(k)))
          if(zetta(1).le.25.0_WP)then
            f2(:,k) = zetta(1)/(exp(zetta(1))-1.0_WP+eps)
          else 
            f2(:,k) = 1.0_WP
          end if

          ! Only evaporation. No condensation modeled
          spray%omega_vap(k)   = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext/di(:,k)**2))/Re)
          spray%omega_vapdm(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext/di(:,k)))/Re)
          spray%omega_vapd2(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext))/Re)
          spray%omega_vapd3(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext*di(:,k)))/Re)

          spray%omega_exp(k)   = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl/di(:,k)**2)))
          spray%omega_expdm(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl/di(:,k))))
          spray%omega_expd2(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl)))
          spray%omega_expd3(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_swl*di(:,k))))

          Qd = (1.0_WP/(Pr_g(k)*VRg(k)))*(spray%Cp_g(k)/spray%C_l(k))*f2(:,k)*(Tg(k)-(Tboil/Tfuel))*Nud(:,k)/di(:,k)
          K_T = sum((6.0_WP*Qd/di(:,k)-1.5_WP*CR(k)*(K_vap_ext/HRg1+K_vap_swl*Re/HRg2)/(di(:,k)**2))*dsd(:,k))/Re

          spray%omega_T(k) = K_T*rho(k)*Y_l(k)

          if(isnan(spray%omega_vap(k)).or.isnan(spray%omega_vapdm(k)).or.isnan(spray%omega_vapd2(k)).or.isnan(spray%omega_vapd3(k))) then
             write(*,*) 'Omega_vap is NaN:',spray%omega_vap(k),spray%omega_vapdm(k),spray%omega_vapd2(k),spray%omega_vapd3(k)
             stop
          end if
       
          if(isnan(spray%omega_exp(k)).or.isnan(spray%omega_expdm(k)).or.isnan(spray%omega_expd2(k)).or.isnan(spray%omega_expd3(k))) then
             write(*,*) 'Omega_exp is NaN:',spray%omega_exp(k),spray%omega_expdm(k),spray%omega_expd2(k),spray%omega_expd3(k)
             stop
          end if

          if(isnan(spray%omega_T(k))) then
             write(*,*) 'Omega_T is NaN:',spray%omega_T(k)
             stop
          end if

        else 
    ! Clausius-Clapeyron relation has to be used carefully after comparing with measured vapor
    ! pressure. It is better to use vapor pressure directly if available.
    !Xeq = spray%f_Pv*(101325.0_WP/spray%P_a)*exp((Lv*spray%MW_f/spray%R_gas)*(1.0_WP/TBd-1.0_WP/T_d))

    Xeq = spray%f_Pv*(spray%p_vap/spray%P_a)

    Lk = spray%visc_ra*sqrt(2.0_WP*Pi*T_d*spray%R_gas/spray%MW_f)/(Sc_g*spray%P_a)

    !Ystar_fe = 1.0_WP/(1.0_WP + spray%WR*(spray%P_a/spray%p_vap-1.0_WP))

    Ystar_fe = Xeq/(Xeq + (1.0_WP-Xeq)*spray%WR)
    
    where(Ystar_fe < 0.0_WP) Ystar_fe = 0.0_WP

    !spray%Y_ref = spray%Cevap*spray%Y_v
    !spray%Y_ref = spray%Cevap*spray%Y_v*spray%Y_v/spray%Y_g
    !spray%Y_ref = (spray%Y_v/spray%Y_g + 2.0_WP*Ystar_fe)/3.0_WP

    ! Remark:
    ! This has large influence on liquid penetration of 
    ! different fuels. So far 1/3rd rule couldn't capture 
    ! trends. The following model for reference vapor mass 
    ! fraction captures these trends better. The droplets
    ! encounter hot ambient gas from the perifery due to
    ! mixing. To replicate these conditions, the avg. vapor
    ! mass fraction Y_v is scaled down by a coefficient,
    ! which decreases as the spray becomes wider mimicking
    ! more hot gas entrainment.
    !spray%Y_ref = (-1.5_WP/spray%b)*spray%Y_v
    !spray%Y_ref = (0.5_WP/spray%b)*(spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP
    spray%Y_ref = (0.5_WP/spray%b)*(spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP
    !spray%Y_ref = (spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP
    !spray%Y_ref = spray%Y_v

!!$    do k = spray%kmin,spray%kmax
!!$       spray%Cevap = max(0.0_WP,rand_normal(spray%Zmix_g(k),sqrt(spray%Zvar_g(k))))
!!$       if(spray%Cevap < spray%Zmix_g(k)) then
!!$          spray%Cevap = 0.0_WP
!!$       else
!!$          spray%Cevap = 1.0_WP
!!$       end if
!!$       spray%Y_ref(k) = spray%Cevap*spray%Y_v(k)
!!$    end do
!!$       spray%Y_ref(k) = min(spray%Cevap/spray%Zmix_g(k),1.0_WP)*spray%Y_v(k)
!!$       if(spray%Zmix_g(k) > 0.0_WP) then
!!$          write(*,*) spray%Y_v(k), spray%Y_ref(k)
!!$       end if
!!$       spray%Y_ref(k) = spray%Cevap*spray%Y_v(k)
!!$       if(spray%Zmix_g(k) > 0.0_WP) then
!!$          write(*,*) spray%Cevap, spray%Y_v(k), spray%Y_ref(k)
!!$       end if
!!$    end do

    !spray%Y_ref = spray%Y_v - 3.0_WP*spray%Zvar_g**0.5/spray%Zmix_g
    !where(spray%Y_ref < 0.0_WP) spray%Y_ref = 0.0_WP
    
    !write(*,*) spray%Y_v, spray%Zvar_g**0.5, spray%Y_ref

    Bdeq = (Ystar_fe - Y_ref)/(1.0_WP-Ystar_fe)

    where(Bdeq < 0.0_WP) Bdeq = 0.0_WP

          zeta = 0.5_WP*(Pr_g/Sc_g)*log(1.0_WP+Bdeq)

          Shd(:,k) = (2.0_WP + 0.552_WP*((Red(:,k))**0.5_WP)*((Sc_g(k))**(1.0_WP/3.0_WP)))
         
          zetta = zeta(k)*Shd(:,k)
          
          Xneq = Xeq(k) - 2.0_WP*(Lk(k)/D(:,k))*zetta
          where(Xneq <0.0_WP) Xneq = 0.0_WP

          Ystar_f = Xneq/(Xneq + (1.0_WP-Xneq)*spray%WR)
          where(Ystar_f < 0.0_WP) Ystar_f = 0.0_WP

          !Yref = spray%Y_ref(k) !spray%Y_v(k)
          !Yref = (spray%Y_v(k)/spray%Y_g(k) + 2.0_WP*Ystar_f)/3.0_WP
          Yref = (0.5_WP/spray%b(k))*(spray%Y_v(k) + 2.0_WP*Ystar_f)/3.0_WP
          !Yref = (spray%Y_v(k) + 2.0_WP*Ystar_f)/3.0_WP

          Bd = (Ystar_f - Yref)/(1.0_WP-Ystar_f)
          where(Bd < 0.0_WP) Bd = 0.0_WP

          ! Constant profile for Y_ref
          coeff = log(1.0_WP+Bd)

          zetta =  0.5_WP*(Pr_g(k)/Sc_g(k))*coeff*Shd(:,k)

          !call bisection(zetta)

          K_vap_ext = 4.0_WP*coeff*Shd(:,k)/(VRg(k)*Sc_g(k))

          ! Only evaporation. No condensation modeled
          spray%omega_vap(k)   = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext/di(:,k)**2))/Re)
          spray%omega_vapdm(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext/di(:,k)))/Re)
          spray%omega_vapd2(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext))/Re)
          spray%omega_vapd3(k) = max(0.0_WP,(1.5_WP*rho(k)*Y_l(k))*sum((dsd(:,k)*K_vap_ext*di(:,k)))/Re)

          Nud(:,k) = (2.0_WP + 0.552_WP*((Red(:,k))**0.5_WP)*((Pr_g(k))**(1.0_WP/3.0_WP)))
          !Qd = (1.0_WP/(Pr_g(k)*VRg(k)))*(spray%Cp_g(k)/spray%C_l(k))*(zetta/(exp(zetta)-1.0_WP+eps))*(spray%T_a/spray%T_fuel-Td(k))*Nud(:,k)/di
          Qd = (1.0_WP/(Pr_g(k)*VRg(k)))*(spray%Cp_g(k)/spray%C_l(k))*(zetta/(exp(zetta)-1.0_WP+eps))*(Tg(k)-Td(k))*Nud(:,k)/di(:,k);
          K_T = sum((6.0_WP*Qd/di(:,k)-1.5_WP*De*CR(k)*LR(k)*K_vap_ext/di(:,k)**2)*dsd(:,k))/Re

          spray%omega_T(k) = K_T*rho(k)*Y_l(k)

          if(spray%omega_vap(k) /= spray%omega_vap(k) .or. spray%omega_T(k) /= spray%omega_T(k) ) then
             write(*,*) 'Omega_vap/T NaN'
          end if
       end if
     end if

    end do

  contains

    subroutine getDist(A,B,Bmax)
      implicit none
      real(WP), dimension(:), intent(in) :: A
      real(WP), dimension(:), intent(out) :: B
      real(WP), intent(in) :: Bmax
      ! -------------------------------
      integer :: i, n,Aidx

      n = size(A)

      Aidx = maxloc(A,1)

      B(1:Aidx) = 2.0_WP*Bmax*(/ (real(i,WP),i=0,Aidx,1) /)/real(Aidx,WP)
      B(Aidx+1:n) = 2.0_WP*Bmax*(real(n,WP) -(/ (real(i,WP),i=Aidx+1,n,1)/))/real(n-Aidx,WP)

      !write(*,*) Aidx, B(Aidx), B(Aidx+1)

    end subroutine getDist

  end subroutine evaporationModel

   subroutine createSelfSimilarProfiles(spray)
     implicit none

     ! ---------------------------------
     type(spray_t), intent(inout) :: spray

     ! ---------------------------------
     integer :: k, ierr, i
     real(WP), dimension(:,:), pointer :: r=>null()
     real(WP), dimension(:), pointer :: xi=>null(), alpha=>null()
     real(WP), dimension(spray%nzo) :: Vf_a, Vf_v, Vf_l
     real(WP), dimension(spray%nro+1) :: rho, phi
     real(WP) :: Vfc_a, Vfc_v, Vfc_l, sum_V, bR2, rho_hat, phi_hat, phi_c, rho_c, h, integral
     real(WP) :: bratio, alphaa, betaa, RR
     character(len=512) :: fname,fname1

     r => spray%r; xi => spray%xi; alpha => spray%alpha

     Vf_a = spray%rho*spray%Y_a*spray%DRa
     Vf_v = spray%rho*spray%Y_v*spray%DRv
     Vf_l = spray%rho*spray%Y_l

     h = 1.0_WP/spray%nro
     xi = (/ (h*real(i,WP),i=0,spray%nro,1) /)
     !write(*,*) xi
     do k = spray%kmin,spray%kmax
        
        ! Assuming self-similar profiles for vapor and liquid phases
        Vfc_v = Vf_v(k)*(alpha(k)+1.0_WP)*(alpha(k)+2.0_WP)/alpha(k)**2
        Vfc_l = Vf_l(k)*(alpha(k)+1.0_WP)*(alpha(k)+2.0_WP)/alpha(k)**2

        ! Balance the rest with air
        Vfc_a = 1.0_WP - Vfc_v - Vfc_l
        
        ! Reconfirm summation of all volume fractions to 1.0_WP
        sum_V = Vfc_a + Vfc_v + Vfc_l
        Vfc_a = Vfc_a/sum_V
        Vfc_v = Vfc_v/sum_V
        Vfc_l = Vfc_l/sum_V

        ! Build gas (ambient+vapor) volume fraction profile
        ! This will be later used in computation of soot mass
        spray%gVF2D(k,:) = 1.0_WP - Vfc_l*(1.0_WP-xi**alpha(k))**2

        ! Build gas (ambient+vapor) density profile
        ! This will be later used in computation of soot mass
        spray%rho_g2D(k,:) = (1.0_WP - Vfc_v*(1.0_WP-xi**alpha(k))**2 - Vfc_l*(1.0_WP-xi**alpha(k))**2)*spray%rho_ra &
                           + Vfc_v*(1.0_WP-xi**alpha(k))**2*spray%rho_rv
        spray%rho_g2D(k,:) = spray%rho_g2D(k,:)/spray%gVF2D(k,:)
        spray%rho_g2D(spray%kmino:spray%kmin-1,:) = MIN_REAL_WP

        ! Build mixture density profile
        rho = 1.0_WP/spray%DRa +((1.0_WP/spray%DRv-1.0_WP/spray%DRa)*Vfc_v+(1.0_WP-1.0_WP/spray%DRa)*Vfc_l)*((1.0_WP-xi**alpha(k))**2)
        spray%rho2D(k,:) = rho

        integral = 2.0_WP*sum(rho*xi*h)
        ! R from Musculus et al. (2009) is not the same as could be computed from jet half-width
        RR = spray%b(k)*1.545_WP !spray%z(k)*spray%beta*1.545_WP
        
        ! Testing: NOT FINAL
        !RR = spray%b(k)*2.67*(1.0/spray%DRa)**0.16_WP

        bR2 = (spray%b(k)/RR)**2!(1.0_WP-sqrt(0.5_WP))**(2.0_WP/alpha(k))
        ! Since there will be numerical error in the evaluation of integral, 
        ! it makes sense to use rho_hat computed from the profile here
        rho_hat = (integral-1.0_WP/spray%DRa)/bR2+1.0_WP/spray%DRa
        !write(*,*) h, integral, 1.0_WP/spray%DRa, spray%rho(k), rho_hat, spray%rho(k)/rho_hat
        
        phi = (1.0_WP - xi**alpha(k))**2
        integral = 2.0_WP*sum(rho*phi*xi*h) 
        phi_c = rho_hat*spray%Zmix_g(k)*bR2/integral
        spray%Zmix2D(k,:) = phi_c*phi

        bratio = spray%b(k)/spray%b_old(k)
        alphaa = 2.0_WP*((1+bratio**2)**2/(1.0_WP+2.0_WP*bratio**2)**1.5_WP)**((1.0_WP+2.0_WP*bratio**2)/(2.0_WP*bratio**2));
        betaa = ((1.0_WP+2.0_WP*bratio**2)/(1.0_WP+bratio**2))-1.0_WP;
        !alphaa = 1.0_WP
        betaa = 1.0_WP

        phi = betaa**2*spray%Zmix2D(k,:)*(alphaa*phi_c-spray%Zmix2D(k,:))
        integral = 2.0_WP*sum(rho*phi*xi*h)
        phi_c = rho_hat*spray%Zvar_g(k)*bR2/integral
        spray%Zvar2D(k,:) = phi_c*phi
        spray%r(k,:) = xi*spray%b(k)/sqrt(bR2)
     end do
     
     if(spray%irank .eq. 0) then

        if (mod(spray%step,spray%datafreq) == 0 .or. spray%end) then
           
           write(fname,"(A,A1,I0.6)") trim('Zmix2D.out'), '_', spray%step

           open(unit=100,file=trim(fname),form="formatted",status="replace",action="write")

           write(100,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro,(spray%xi(i), i=1,spray%nro+1)

           write(fname1,"(A,A1,I0.6)") trim('Zvar2D.out'), '_', spray%step

           open(unit=101,file=trim(fname1),form="formatted",status="replace",action="write")

           write(101,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro,(spray%xi(i), i=1,spray%nro+1)

           do k = spray%kmin,spray%kmax
              
              write(100, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Zmix2D(k,i), i=1,spray%nro+1)

              write(101, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Zvar2D(k,i), i=1,spray%nro+1)

           end do

           close(unit=100)
           close(unit=101)
        end if
     end if
    
   end subroutine createSelfSimilarProfiles

  subroutine compute_alpha(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP) :: ratio,alpha
    integer :: k

    do k = spray%kmin,spray%kmax

       ratio = (spray%b(k)/0.5_WP)**2

       call NewtonRoot(alpha,ratio,spray%DRa)

       if(alpha > 1.5_WP) then
          spray%alpha(k) = alpha
       else
          spray%alpha(k) = max(1.5_WP,alpha)
       end if
       
    end do
    
  contains
    subroutine NewtonRoot(xfinal,AA0,DR)
      implicit none

      real(WP), intent(inout) :: xfinal
      real(WP), intent(in) :: DR, AA0
      ! -------------------------------
      real(WP), parameter :: tol = 1E-06
      real(WP) :: x0, x1, fx0, fdashx0, err
      integer :: iter
      ! -------------------------------

      iter = 0
      err = 1.0_WP
      x0 = 100.0_WP
      do while ( err > tol)
         iter = iter + 1
         fx0 = 6.0_WP*(AA0-1.0_WP)*x0**4-(7.0_WP+18.0_WP*DR)*x0**3 &
             - (2.0_WP+33.0_WP*DR)*x0**2 - 20.0_WP*DR*x0 - 4.0_WP*DR
         fdashx0 = 24.0_WP*(AA0-1.0_WP)*x0**3-3.0_WP*(7.0_WP+18.0_WP*DR)*x0**2 &
                - 2.0_WP*(2.0_WP+33.0_WP*DR)*x0 - 20.0_WP*DR
         x1 = x0 - fx0/fdashx0

         err = abs(x1-x0)

         x0 = x1
      end do
      !write(*,*) 'Newton root iters:', iter
      xfinal = x1

    end subroutine NewtonRoot
  end subroutine compute_alpha
  subroutine dragModel(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: eps = 1E-16_WP
    real(WP), dimension(:,:), pointer ::  di=>null(), dsd=>null(), Red=>null(), CD=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), u_l=>null(), u_g=>null(), &
                                       DRg=>null(), VRg=>null(), h=>null(), &
                                       DRl=>null(), VRl=>null()
    real(WP), pointer :: Re=>null()
    real(WP), dimension(spray%nzo) :: u_rel, u_reli

    real(WP), dimension(spray%nd) :: Hvs
    real(WP) :: CDd, Ref
    integer :: k

    if(.not. spray%drag_model) return

    spray%f_drag = 0.0_WP; spray%Red = 0.0_WP; spray%CD = 0.0_WP

    rho => spray%rho; Y_l => spray%Y_l; u_l => spray%u_l; u_g => spray%u_g; 
    di => spray%di; h => spray%h; dsd => spray%dsd; Red => spray%Red
    Re => spray%Re; DRg => spray%DRg; VRg => spray%VRg; CD => spray%CD
    DRl => spray%DRl; VRl => spray%VRl

    u_reli = spray%Crel*u_g-u_l !spray%Crel*abs(u_g-u_l)
    u_rel = abs(u_reli) !spray%Crel*abs(u_g-u_l)

    Hvs = 0.0_WP

    do k = spray%kmin,spray%kmax

       ! Factors to scale Re_l for local droplet condition
       Ref = 1.0_WP !DRl(k)/VRl(k)
       
       Red(:,k) = (Re*VRg(k)/DRg(k))*u_rel(k)*di(:,k) *Ref

       !Red(1,k) = Red(2,k)

       where(Red(:,k) <= 1000.0_WP) Hvs = 1.0_WP
       
       if ( spray%d2(k) > eps .and. spray%dm(k) > eps .and. u_rel(k) > eps) then

          CD(:,k) = (24.0_WP/(Red(:,k)+eps) + 4.0_WP/((Red(:,k)+eps)**(1.0_WP/3.0_WP)))*Hvs + 0.424_WP*(1.0_WP-Hvs)
          CDd = sum(CD(:,k)/di(:,k)*dsd(:,k))

          ! Only the drag force acting against droplets considered
          spray%f_drag(k) = min(0.0_WP,0.75_WP*CDd*rho(k)*Y_l(k)*u_rel(k)*u_reli(k)/DRg(k))
         
       end if
      
    end do
    !spray%Red = Red
    !spray%CD = CD
  end subroutine dragModel

   subroutine breakupModel(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: eps = 1E-16_WP
    real(WP), dimension(:,:), pointer :: di=>null(), dsd=>null(), Red=>null(), CD=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), u_l=>null(), u_g=>null(), &
                                       DRg=>null(), VRg=>null(), h=>null(), Td=>null(), &
                                       SR=>null(), DRl=>null(), VRl=>null(), Rdot=>null(), Nbub=>null(), db=>null()
    real(WP), pointer :: Re=>null(), We=>null()
    real(WP), dimension(spray%nzo) :: T_d, u_rel, u_rel2, K_bre1, K_bre2, K_bre3, K_bref1, K_bref2, K_bref3

    real(WP), dimension(spray%nd) :: Weg, Wel, Rel, Z, T, OmegaKH, LambdaKH, tauKH, dKH, dKH1, dKH2, dKH21, dKH22, &
                                     OmegaRT, LambdaRT, tauRT, dRT, critKH, critKH1, critRT, tb, dst, diff, Rdot_new, X_coef, Y_coef, Z_coef, &
                                     dRddt, tauTHM, eta_thm, omega_thm
    real(WP) :: L_breakup, Wef, Ref
    real(WP), pointer :: rhol=>null(), rhov=>null(), viscv=>null(), viscl=>null(), sigma=>null(), Tfuel=>null(), Tboil=>null()
    real(WP), dimension(:,:), pointer :: Rb=>null()
    real(WP), dimension(:), pointer :: rhog=>null(), viscg=>null()
    real(WP) :: Rcrit, deltaT, Nden, error, Vd, Pdif        

    integer :: k

    if(.not. spray%bre_model) return

    u_rel = 0.0_WP; u_rel2 = 0.0_WP; K_bre1 = 0.0_WP; K_bre2 = 0.0_WP; K_bre3 = 0.0_WP;
    K_bref1 = 0.0_WP; K_bref2 = 0.0_WP; K_bref3 = 0.0_WP;

    Weg = 0.0_WP; Wel = 0.0_WP; Rel = 0.0_WP; Z = 0.0_WP; T = 0.0_WP; 
    OmegaKH = 0.0_WP; LambdaKH = 0.0_WP; tauKH = 0.0_WP; dKH = 0.0_WP; 
    dKH1 = 0.0_WP; dKH2 = 0.0_WP; dKH21 = 0.0_WP; dKH22 = 0.0_WP;
    OmegaRT = 0.0_WP; LambdaRT = 0.0_WP; tauRT = 0.0_WP; dRT = 0.0_WP; 
    critKH = 0.0_WP; critKH1 = 0.0_WP; critRT = 0.0_WP; 
    tauTHM = 0.0_WP; db = 0.0_WP; diff = 0.0_WP;

    spray%omega_bre1 = 0.0_WP; spray%omega_bre2 = 0.0_WP; spray%omega_bre3 = 0.0_WP

    T_d = spray%T_fuel*spray%Td

    rho => spray%rho; Y_l => spray%Y_l; u_l => spray%u_l; u_g => spray%u_g; Td => spray%Td; Tfuel => spray%T_fuel; Tboil => spray%T_boil
    di => spray%di; h => spray%h; dsd => spray%dsd; Red => spray%Red; CD => spray%CD
    Re => spray%Re; We => spray%We; DRg => spray%DRg; VRg => spray%VRg; 
    SR =>spray%SR; DRl => spray%DRl; VRl => spray%VRl
    Rdot => spray%Rdot; Nbub => spray%Nbub; db => spray%db

    u_rel = abs(spray%Crel*u_g-u_l)

    u_rel2 = u_rel**2

    ! Breakup length
    L_breakup = 0.5_WP*spray%B1*sqrt(spray%DRa)

    do k = spray%kmin,spray%kmax
       
       if ( spray%d2(k) > eps .and. spray%dm(k) > eps .and. u_rel(k) > eps ) then

          ! Factors to scale We_l for local droplet condition
          Wef = 1.0_WP !DRl(k)*SR(k)
          !if(SR(k)==0.0_WP) Wef = MAX_REAL_WP
          ! Factors to scale Re_l for local droplet condition
          Ref = 1.0_WP !DRl(k)/VRl(k)

          ! Kelvin-Helmoltz WAVE model
          !B0 = 0.61; B1 = 0.01_WP
          Weg = 0.5_WP*We*u_rel2(k)*di(:,k)/DRg(k) *Wef
          Wel = 0.5_WP*We*di(:,k)*u_rel2(k) *Wef
          Rel = 0.5_WP*Re*di(:,k)*u_rel(k) *Ref
          Z = sqrt(Wel)/(Rel)
          T = Z*sqrt(Weg)

          OmegaKH = (0.34_WP+0.38_WP*Weg**1.5_WP)/(1.0_WP+Z)/(1.0_WP+1.4_WP*T**0.6_WP)*sqrt(8.0_WP/We /Wef /di(:,k)**3)
          ! Avoid division by zero
          !where(Wel == 0.0_WP) OmegaKH = MAX_REAL_WP 

          LambdaKH = (9.02_WP/2.0_WP)*di(:,k)*(1.0_WP+0.45_WP*sqrt(Z))*(1.0_WP+0.4_WP*T**0.7_WP)/&
               (1.0_WP+0.865_WP*Weg**1.67_WP)**0.6_WP
          tauKH = (3.788_WP*spray%B1/2.0_WP)*di(:,k)/(OmegaKH*LambdaKH)
          dKH1 = 2.0_WP*spray%B0*LambdaKH

          dKH21 = (3.0_WP*Pi*di(:,k)**2*u_rel(k)/OmegaKH)**(1.0_WP/3.0_WP)
          dKH22 = (1.5_WP*di(:,k)**2*LambdaKH)**(1.0_WP/3.0_WP)

          where(dKH21 < dKH22) critKH1 = 1.0_WP
          dKH2 = critKH1*dKH21 + (1.0-critKH1)*dKH22

          where(dKH1 > 1.0_WP) critKH = 1.0_WP
          dKH = (1.0_WP-critKH)*dKH1 + critKH*dKH2

          ! Rayleigh-Taylor model
          !C3 = 1.0_WP
          OmegaRT = 0.5_WP*(1.0_WP/DRg(k))*sqrt((1.0_WP/(DRg(k)+1.0_WP))) &
               * ((DRg(k)-1.0_WP)*CD(:,k)*u_rel2(k)/di(:,k))**(0.75_WP)*(We *Wef )**0.25_WP
          LambdaRT = spray%C3*4.0_WP*Pi*DRg(k)*sqrt(1.0_WP/We /Wef /(DRg(k)-1.0_WP))*sqrt(di(:,k)/CD(:,k)/u_rel2(k))
          tauRT = 1.0_WP/OmegaRT
          dRT = LambdaRT

          ! Corresponding to modified KH-RT breakup model in CONVERGE: First/Primary breakup is always KH breakup
          if(k > 4) then 
             where(di(:,k)>dRT .and. tauRT<tauKH) critRT = 1.0_WP
          end if

          tb = critRT*tauRT + (1.0_WP-critRT)*tauKH;
          dst = critRT*dRT + (1.0_WP-critRT)*dKH;

          diff = di(:,k) - dst
          where(diff < 0.0_WP) diff = 0.0_WP
          K_bre1(k) = sum(diff/tb*dsd(:,k))
          K_bre2(k) = 2.0_WP*sum(di(:,k)*diff/tb*dsd(:,k))
          K_bre3(k) = 3.0_WP*sum(di(:,k)*di(:,k)*diff/tb*dsd(:,k))

          if(K_bre1(k) /= K_bre1(k)) then
             K_bre1(k) = 0.0_WP
          end if

       call computeSatBubbleProperties(spray,Tfuel)
       Rcrit = 2.0*spray%sigma/(spray%P_bub-spray%P_a)

       call computeSatBubbleProperties(spray,T_d(k))
       !Rcrit = 2.0*spray%sigma/(spray%P_bub-spray%P_a)

       if (spray%flash_boiling .and. spray%P_bub-spray%P_a.gt.5.0_WP)then
          Pdif = spray%P_bub-spray%P_a 
          Rdot = sqrt((2.0_WP/3.0_WP)*Pdif/spray%rho_l) 

          Vd = (Pi/6.0_WP)*(di(spray%nd,k)*spray%D_eff)**3.0_WP
          !Nden = 5.757e12_WP*exp(-5.279_WP/(T_d(k) - Tboil))+eps  
          Nden = 5.757e12_WP*exp(-5.279_WP/(spray%T_fuel - Tboil))+eps  
          Nbub = ceiling(Nden*Pi*((di(:,k)*spray%D_eff)**3.0_WP)/6.0_WP)

          eta_thm = (6.0_WP*0.55_WP*Vd/(Pi*Nbub*0.45))**(1.0_WP/3.0_WP)
          omega_thm = 2.0_WP*Rdot/spray%D_eff 

          tauTHM = spray%D0*(di(spray%nd,k)*spray%D_eff)/(eta_thm*omega_thm*spray%tau)
          dst = spray%D1*eta_thm/spray%D_eff

          diff = di(:,k)- dst

          if(diff(1)<0.0_WP)then
           where(diff < 0.0_WP) diff = 0.0_WP
           !write(*,*)'dstb is less than di: dst and dd', dst, di(1,k)*spray%D_eff
           !stop
          end if 

          K_bref1(k) = sum(diff/tauTHM*dsd(:,k))     
          K_bref2(k) = 2.0_WP*sum(di(:,k)*diff/tauTHM*dsd(:,k))    
          K_bref3(k) = 3.0_WP*sum(di(:,k)*di(:,k)*diff/tauTHM*dsd(:,k)) 

          if(isnan(K_bref1(k)).or.isnan(K_bref2(k)).or.isnan(K_bref3(k)))then
           write(*,*)'NaN resulting from thermal breakup model', K_bref1(k), K_bref2(k), K_bref3(k)
           stop
          end if 

        else 
          tauTHM = 9999.0_WP
        end if

        if(tauTHM(1)<tauRT(1).and.tauTHM(1)<tauKH(1))then
          K_bre1(k) = K_bref1(k)   
          K_bre2(k) = K_bref2(k)     
          K_bre3(k) = K_bref3(k)   
        end if

       end if
    end do

    spray%omega_bre1 = K_bre1*rho*Y_l
    spray%omega_bre2 = K_bre2*rho*Y_l
    spray%omega_bre3 = K_bre3*rho*Y_l

  end subroutine breakupModel

  subroutine turbulenceModel(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: k
    real(WP), dimension(spray%nzo) :: factor

    ! Assuming only gas phase turbulence induced by gas phase velocity

    if(spray%ck_type == 'Dynamic') then
       spray%c_k = spray%omega_ent/(spray%VRtg/spray%Re + 1.0_WP/spray%VRg/spray%Re)
    end if
    !spray%c_k = min(spray%c_k,12.0_WP)
    !spray%c_k(spray%kmino:spray%kmin) = 0.0_WP
    !spray%c_k = sum(spray%rho*spray%Y_g*spray%b**2*spray%c_k)/sum(spray%rho*spray%Y_g*spray%b**2)

    !spray%omega_k_g_p = 0.5_WP*spray%c_k*spray%u_g**2*spray%mu_t_g
    spray%omega_k_g_p = spray%c_k*spray%u_g**2/spray%Re*spray%VRtg*spray%Y_g
    spray%omega_k_g_d = spray%rho*spray%Y_g*spray%eps_g*spray%b*spray%b

    factor = spray%eps_g/spray%k_g
    
    !where (spray%k_g == 0.0_WP) factor = 0.0_WP
    where (factor /= factor .or. spray%k_g == 0.0_WP) factor = 0.0_WP
    
    spray%omega_eps_g_p = factor*spray%c_eps1*spray%omega_k_g_p
    spray%omega_eps_g_d = factor*spray%c_eps2*spray%omega_k_g_d

    !spray%omega_zvar_g_p = spray%c_k*spray%zmix_g**2*spray%mu_t_g
    spray%omega_zvar_g_p = spray%c_k*spray%zmix_g**2/spray%Re*spray%VRtg*spray%Y_g
    spray%omega_zvar_g_d = spray%c_zvar*factor*spray%rho*spray%Y_g*spray%zvar_g*spray%b*spray%b

  end subroutine turbulenceModel

  subroutine computeTimeStep(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    real(WP), dimension(spray%nzo) :: u_rel
    
    u_rel = abs(spray%Crel*spray%u_g-spray%u_l) !spray%Crel*abs(spray%u_g-spray%u_l)

    spray%CFL_bre = u_rel/spray%dz*spray%dt

    spray%CFL_conv = max(spray%u_l,spray%u_g)*spray%dt/spray%dz

    spray%CFL = spray%MaxCFL*spray%MaxCFL/max(spray%MaxCFL,maxval(max(spray%CFL_conv,spray%CFL_bre)))

    !spray%CFL = spray%MaxCFL

    spray%dt = spray%CFL*spray%dz

    !if (spray%step < 100) spray%dt = 0.1_WP*spray%dt

    !if (spray%ndtime > 0.15e-3/spray%tau .and. spray%ndtime < 0.25e-3/spray%tau ) spray%dt = 0.1_WP*spray%dt

  end subroutine computeTimeStep

  subroutine advanceTime(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    
    real(WP) :: dt_now

    spray%step = spray%step + 1

    if((spray%ndftime-spray%ndtime)<spray%dt) then
        dt_now = spray%ndftime-spray%ndtime
    else
        dt_now = spray%dt
    end if

    spray%ndtime = spray%ndtime + dt_now

  end subroutine advanceTime

  subroutine read_ROI_from_file(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------  
    real(WP), dimension(:,:), pointer :: scaled
    logical :: exist_file
    character(len=128) :: cmd, line, rnk
    integer :: nlines, ioerr, i, idx
    real(WP) :: max_roi

    if(associated(spray%roi)) return

    inquire(file=trim(spray%roi_file), exist=exist_file)

    i = 1

    if(exist_file) then
       
       write(rnk,"(I0.2)") spray%irank

!!$       cmd = "cat "//trim(spray%roi_file)//" | sed '/^\s*#/d;/^\s*$/d' | wc -l > nlines.txt_"//trim(rnk)
!!$       call system(cmd)
!!$       open(unit=104,file='nlines.txt_'//trim(rnk))
!!$       read(unit=104,fmt='(i6)',iostat=ioerr) nlines
!!$       cmd = 'rm nlines.txt_'//trim(rnk)
!!$       call system(cmd)
!!$       close(unit=104)
       nlines = 0 
       open (unit=104, file =trim(spray%roi_file))
       do
          read(unit=104,fmt=*,iostat=ioerr)
          if (ioerr/=0) exit
          nlines = nlines + 1
       end do
       close(unit=104)

       !nlines = 10000

       allocate(spray%roi(nlines,2)); spray%roi = 0.0_WP
       open(unit=105,file=trim(spray%roi_file),form="formatted",status="old",action="read")
       write(*,*) 'Reading ROI profile...'
       do while (.true.)
          read(unit=105,fmt='(a)',iostat=ioerr) line

          if (ioerr .ne. 0) then
             exit
          end if

          read(line,*)  spray%roi(i,1), spray%roi(i,2)

          i = i + 1
       end do
       close(unit=105)

       allocate(scaled(nlines,2))

       max_roi = maxval(spray%roi(:,2))

       scaled(:,1) = spray%roi(:,1)/spray%tau
       scaled(:,2) = spray%roi(:,2)/max_roi

       spray%roi = scaled
       !spray%roi(1,2) = 0.0_WP

       deallocate(scaled)

    end if

  end subroutine read_ROI_from_file

!!$  ! Liquid Fuel Property Table should be in the following format
!!$  ! Header should be commented using '#' and physical properties
!!$  ! as a function of temperature
!!$  ! Temperature[K] Viscosity[Ns/m^2] SurfTension[N/m] HeatOfVap[J/kg] VaporPres[Pa] ThermConductivity[W/(mK)] Density[kg/m^3] SpecificHeatCapacity[J/(kgK)]
!!$  subroutine readLiquidFuelPropertiesTable(spray)
!!$    implicit none
!!$
!!$    ! ---------------------------------
!!$    type(spray_t), pointer, intent(inout) :: spray
!!$
!!$    ! ---------------------------------  
!!$    logical :: exist_file
!!$    character(len=128) :: cmd, line
!!$    integer :: nlines, ioerr, i
!!$
!!$    inquire(file=trim(spray%LFPTname), exist=exist_file)
!!$
!!$    i = 1
!!$
!!$    if(exist_file) then
!!$       cmd = "cat "//trim(spray%LFPTname)//" | sed '/^\s*#/d;/^\s*$/d' | wc -l > nlfpt.txt"
!!$       call system(cmd)
!!$       open(unit=1040,file='nlfpt.txt')
!!$       read(unit=1040,fmt='(i6)',iostat=ioerr) nlines
!!$       cmd = 'rm nlfpt.txt'
!!$       call system(cmd)
!!$       close(unit=1040)
!!$       allocate(spray%LFPT(nlines,8)); spray%LFPT = 0.0_WP
!!$       open(unit=1050,file=trim(spray%LFPTname),form="formatted",status="old",action="read")
!!$       write(*,*) 'Reading Liquid Fuel Properties from Table...'
!!$       do while (.true.)
!!$          read(unit=1050,fmt='(a)',iostat=ioerr) line
!!$
!!$          if (ioerr .ne. 0) then
!!$             exit
!!$          end if
!!$
!!$          if ( index(line,'!') > 0 .or. index(line,'#') > 0 ) then
!!$             cycle
!!$          end if
!!$          
!!$          read(line,*)  spray%LFPT(i,1), spray%LFPT(i,2), spray%LFPT(i,3), spray%LFPT(i,4),&
!!$                        spray%LFPT(i,5), spray%LFPT(i,6), spray%LFPT(i,7), spray%LFPT(i,8)
!!$
!!$          !write(*,*) spray%LFPT(i,1), spray%LFPT(i,2), spray%LFPT(i,3), spray%LFPT(i,4), spray%LFPT(i,5), spray%LFPT(i,6), spray%LFPT(i,7), spray%LFPT(i,8)
!!$
!!$          i = i + 1
!!$
!!$       end do
!!$       close(unit=1050)
!!$    end if
!!$
!!$  end subroutine readLiquidFuelPropertiesTable
!!$
!!$  ! Vapor Fuel Property Table should be in the following format
!!$  ! Header should be commented using '#' and physical properties
!!$  ! as a function of temperature
!!$  ! Temperature[K] Viscosity[Ns/m^2] ThermConductivity[W/(mK)] SpecificHeatCapacity[J/(kgK)] DiffusionCoefficient[m^2/s]
!!$  subroutine readVaporFuelPropertiesTable(spray)
!!$    implicit none
!!$
!!$    ! ---------------------------------
!!$    type(spray_t), pointer, intent(inout) :: spray
!!$
!!$    ! ---------------------------------  
!!$    logical :: exist_file
!!$    character(len=128) :: cmd, line
!!$    integer :: nlines, ioerr, i
!!$
!!$    inquire(file=trim(spray%VFPTname), exist=exist_file)
!!$
!!$    i = 1
!!$
!!$    if(exist_file) then
!!$       cmd = "cat "//trim(spray%VFPTname)//" | sed '/^\s*#/d;/^\s*$/d' | wc -l > nlvfpt.txt"
!!$       call system(cmd)
!!$       open(unit=1041,file='nlvfpt.txt')
!!$       read(unit=1041,fmt='(i6)',iostat=ioerr) nlines
!!$       cmd = 'rm nlvfpt.txt'
!!$       call system(cmd)
!!$       close(unit=1041)
!!$       allocate(spray%VFPT(nlines,5)); spray%VFPT = 0.0_WP
!!$       open(unit=105,file=trim(spray%VFPTname),form="formatted",status="old",action="read")
!!$       write(*,*) 'Reading Vapor Fuel Properties from Table...'
!!$       do while (.true.)
!!$          read(unit=1051,fmt='(a)',iostat=ioerr) line
!!$
!!$          if (ioerr .ne. 0) then
!!$             exit
!!$          end if
!!$
!!$          if ( index(line,'!') > 0 .or. index(line,'#') > 0 ) then
!!$             cycle
!!$          end if
!!$          
!!$          read(line,*)  spray%VFPT(i,1), spray%VFPT(i,2), spray%VFPT(i,3), spray%VFPT(i,4), spray%VFPT(i,5)
!!$          i = i + 1
!!$       end do
!!$       close(unit=1051)
!!$    end if
!!$
!!$  end subroutine readVaporFuelPropertiesTable

  subroutine applyBC(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: i
    real(WP) :: eps = 1E-16_WP, turb_intensity

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    ! Left boundary (injector)

    if(associated(spray%roi)) then
       !Interpolate
       do i=1,maxval(shape(spray%roi))-1
          if (spray%ndtime .ge. spray%roi(i,1) .and. spray%ndtime .lt. spray%roi(i+1,1)) then
             spray%u_l(kmino:kmin-1) = (spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2)

             spray%Y_l(kmino:kmin-1) = (spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2)

             spray%d3(kmino:kmin-1) = spray%init_d3*((spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2))
             spray%d2(kmino:kmin-1) = spray%init_d2*((spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2))
             spray%dm(kmino:kmin-1) = spray%init_dm*((spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2))
             spray%dvar(kmino:kmin-1) = spray%init_dvar*((spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2))
          end if
       end do
    else
       spray%u_l(kmino:kmin-1) = 1.0_WP
       spray%Y_l(kmino:kmin-1) = 1.0_WP
       spray%d3(kmino:kmin-1) = spray%init_d3
       spray%d2(kmino:kmin-1) = spray%init_d2
       spray%dm(kmino:kmin-1) = spray%init_dm
       spray%dvar(kmino:kmin-1) = spray%init_dvar
    end if

    spray%rho(kmino:kmin-1) = 1.0_WP

    if(spray%dsd_type(kmin-1) > type_delta) then
       spray%Y_l(kmino:kmin-1) = 1.0_WP
       spray%d3(kmino:kmin-1) = spray%init_d3
       spray%d2(kmino:kmin-1) = spray%init_d2
       spray%dm(kmino:kmin-1) = spray%init_dm
       spray%dvar(kmino:kmin-1) = spray%init_dvar
    end if

    spray%Y_v(kmino:kmin-1) = 0.0_WP
    spray%Y_a(kmino:kmin-1) = 1.0_WP-spray%Y_l(kmino:kmin-1)
    spray%Y_g(kmino:kmin-1) = spray%Y_v(kmino:kmin-1) + spray%Y_a(kmino:kmin-1)

    spray%b(kmino:kmin-1)   = 0.5_WP

    spray%u_g(kmino:kmin-1) = 0.0_WP

    spray%Td(kmino:kmin-1) = 1.0_WP 

    ! Boundary conditions for k, eps, and zvar are in the gas phase grid
    ! i.e. first cell next to the nozzle exit
    ! Based on Tamanini (1981)
    turb_intensity = 0.16_WP/(spray%VRg(kmin)/spray%DRg(kmin)*spray%Re*spray%u_g(kmin))**0.125_WP
    if(spray%u_g(kmin) == 0.0_WP) turb_intensity = 0.0_WP
    
    spray%k_g(kmino:kmin) = 0.5_WP*(turb_intensity*spray%u_g(kmin))**2
    !spray%k_g(kmino:kmin) = 1.5_WP*(turb_intensity*spray%u_g(kmin))**2

    !spray%k_g(kmino:kmin) = 1.0E-04_WP

    ! Based on Tamanini (1981): Production = Dissipation
    !spray%eps_g(kmino:kmin) = sqrt(0.5_WP*spray%c_k*spray%c_mu/spray%rho(kmin)/spray%Y_g(kmin))*(1.0_WP/spray%DRa/spray%DRg(kmin))**(0.25_WP)\
    !                        * spray%u_g(kmin)*spray%k_g(kmin)/spray%b(kmin)

    ! For alternative mu_t_g def
    !spray%eps_g(kmino:kmin) = sqrt(0.5_WP*spray%c_k*spray%c_mu)*(spray%Y_a(kmin)/spray%Y_g(kmin))**(0.25_WP)*spray%k_g(kmin)*spray%u_g(kmin)/spray%b(kmin)
    
    !spray%eps_g(kmino:kmin) = sqrt(0.5_WP*spray%c_k*spray%c_mu)*spray%k_g(kmin)*spray%u_g(kmin)/spray%b(kmin)
    spray%eps_g(kmino:kmin) = 0.5_WP*spray%c_mu**0.75_WP*(spray%k_g(kmin)**1.5_WP)/0.07_WP/spray%b(kmin)
    !spray%eps_g(kmino:kmin) = 1.0E-10_WP

    spray%mu_t_g(kmino:kmin) = 0.0_WP

    spray%zvar_g(kmino:kmin) = 0.0_WP !spray%zvar_g(kmin)

    if(spray%nsc > 0) then
       do i = 1,spray%nsc
          spray%Sc(i,kmino:kmin-1) = 0.0_WP
       end do
    end if
    
    ! Right boundary (Neumann)
    spray%rho(kmax+1:kmaxo) = spray%rho(kmax)

    spray%Y_l(kmax+1:kmaxo) = spray%Y_l(kmax)
    spray%Y_v(kmax+1:kmaxo) = spray%Y_v(kmax)
    spray%Y_a(kmax+1:kmaxo) = spray%Y_a(kmax)
    spray%Y_g(kmax+1:kmaxo) = spray%Y_g(kmax)

    spray%b(kmax+1:kmaxo) = spray%b(kmax)

    spray%u_l(kmax+1:kmaxo) = spray%u_l(kmax)
    spray%u_g(kmax+1:kmaxo) = spray%u_g(kmax)

    spray%d3(kmax+1:kmaxo) = spray%d3(kmax)
    spray%d2(kmax+1:kmaxo) = spray%d2(kmax)
    spray%dm(kmax+1:kmaxo) = spray%dm(kmax)
    spray%dvar(kmax+1:kmaxo) = spray%dvar(kmax)

    spray%Td(kmax+1:kmaxo) = spray%Td(kmax)

    spray%k_g(kmax+1:kmaxo) = spray%k_g(kmax)
    spray%eps_g(kmax+1:kmaxo) = spray%eps_g(kmax)
    spray%zvar_g(kmax+1:kmaxo) = spray%zvar_g(kmax)

    if(spray%nsc > 0) then
       do i = 1,spray%nsc
          spray%Sc(i,kmax+1:kmaxo) = spray%Sc(i,kmax)
       end do
    end if

  end subroutine applyBC

  subroutine write_output(spray,step,time)
    implicit none

    ! ---------------------------------
    type(spray_t), intent(inout) :: spray
    integer, intent(in) :: step
    real(WP), intent(in) :: time
    ! ---------------------------------
    character(len=512) :: rowfmt, rowfmth, fname, stp, tm
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, s, i, nsc
    logical :: exist

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    !write(fname,"(A9,I0.6,A1,ES15.5)") 'data.out_', step, '_', time

    if (mod(spray%step,spray%datafreq) == 0 .or. spray%saveDataFile .or. spray%end) then
       write(fname,"(A,A1,I0.6)") trim(spray%datafilename), '_', step

       open(unit=100,file=trim(fname),form="formatted",status="unknown",action="write")

       rowfmth = '(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)'
       rowfmt = "(ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, &
                  ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, &
                  ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, I2, ES15.5E3, ES15.5E3, &
                  ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"

       write(100,FMT=rowfmth) '# z<1> ', 'rho<2> ', 'Y_l<3> ', 'Y_v<4> ', 'Y_a<5> ', 'Y_g<6> ', 'u_l<7> ', 'u_g<8> ', 'dm<9> ', 'dvar<10> ', 'd2<11> ', 'd3<12> ', 'Td<13> ', 'Tg<14> ', 'b<15> ', 'k_g<16> ', 'eps_g<17> ', 'mu_t_g<18> ', 'zvar_g<19> ', 'zmix_g<20> ', 'chi_g<21> ', 'chi_g_stl<22> ', 'dsd_type<23> ', 'Pr_g<24> ', 'Sc_g<25> ', 'c_k<26> ', 'VRg<27> ', 'VRtg<28> ', 'SR<29> ', 'DRl<30> ', 'Y_ref<31> ', 'alpha<32> ', 'T_v<33> ', 'DRg<34> '
       do k = kmin,kmax
          write(100,FMT=rowfmt) spray%z(k), spray%rho(k), spray%Y_l(k), spray%Y_v(k), spray%Y_a(k), &
               spray%Y_g(k), spray%u_l(k), spray%u_g(k), spray%dm(k), spray%dvar(k), spray%d2(k), &
               spray%d3(k), spray%Td(k), spray%Tg(k), spray%b(k), spray%k_g(k), &
               spray%eps_g(k), spray%mu_t_g(k), spray%zvar_g(k), spray%zmix_g(k), spray%chi_g(k), spray%chi_g_stl(k), spray%dsd_type(k), spray%Pr_g(k), spray%Sc_g(k),spray%c_k(k),spray%VRg(k),spray%VRtg(k),spray%SR(k),spray%DRl(k),spray%Y_ref(k),spray%alpha(k),spray%Tv(k),spray%DRg(k)
       end do

       close(unit=100)

       spray%saveDataFile = .false.

    end if

    if (mod(spray%step,spray%outfreq) == 0 .or. spray%end) then
       write(fname,"(A)") trim(spray%plfilename)

       inquire(file=trim(fname), exist=exist)

       if(exist) then
          open(unit=101,file=trim(fname),form="formatted",status="old",position="append",action="write")
       else
          open(unit=101,file=trim(fname),form="formatted",status="new",action="write")

          rowfmth = '(A,A,A,A,A,A,A,A)'

          write(101,fmt=rowfmth) '# step<1> ', 'time<2> ', 'LPL<3> ', 'VPL<4> ', 'U_tip<5> ', 'D_tip<5> ','Chi_st<6> '

       end if

       rowfmt = "(I0.6,ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"

       !do s = 1,step
       s = 1
       write(101,FMT=rowfmt) spray%step, spray%time(s), spray%LPL(s), spray%VPL(s), spray%utip(s), spray%dtip(s), spray%chi_st(s), spray%chi_st1(s)
       !end do

       close(unit=101)

       if(spray%combustion_model == 'MRIF') then
          write(fname,"(A,A1,I0.6)") trim('flamelet_chi.out')

          if(exist) then
             open(unit=202,file=trim(fname),form="formatted",status="old",position="append",action="write")
          else
             open(unit=202,file=trim(fname),form="formatted",status="new",action="write")
          end if

          nsc = spray%nsc
          !write(202,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

          do k = spray%kmin,spray%kmax
             write(202, '(<nsc+2>(ES15.5E3))') spray%step, spray%time(1), (spray%chi_st_m(i), i=1,nsc)
          end do
          close(unit=202)

       end if

    end if

    if (mod(spray%step,spray%datafreq) == 0 .and. maxval(spray%dsd_type) == 7 .or. spray%end) then
       write(fname,"(A,A1,I0.6)") 'dsd.out', '_', step

       open(unit=102,file=trim(fname),form="formatted",status="unknown",action="write")

       rowfmt = "(ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"

       do k = kmin-1,kmax
          write(102,FMT=rowfmt) spray%z(k), spray%dsdlam(1,k), spray%dsdlam(2,k), spray%dsdlam(3,k), &
               spray%dsdlam(4,k), spray%dm(k), spray%d2(k), spray%d3(k)
       end do

       close(unit=102)
    end if

    if (mod(spray%step,spray%datafreq) == 0 .or. spray%end) then
       write(fname,"(A,A1,I0.6)") 'omega.out', '_', step

       open(unit=103,file=trim(fname),form="formatted",status="unknown",action="write")

       rowfmth = '(A,A,A,A,A,A,A,A,A,A,A,A,A)'

       write(103,fmt=rowfmth) '# z<1> ', 'omega_ent<2> ', 'omega_drag<3> ', 'omega_bre1<4> ', 'omega_bre2<5> ', 'omega_vap<6> ', 'omega_T<7> ', 'omega_k_g_p<8> ', 'omega_k_g_d<9> ', 'omega_eps_g_p<10> ', 'omega_eps_g_d<11> ', 'omega_zvar_g_p<12> ', 'omega_zvar_g_d<13> '

       rowfmt = "(ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"

       do k = kmin,kmax
          write(103,FMT=rowfmt) spray%z(k), spray%omega_ent(k), spray%f_drag(k), spray%omega_bre1(k), spray%omega_bre2(k), spray%omega_vap(k), spray%omega_T(k), spray%omega_k_g_p(k), spray%omega_k_g_d(k)/spray%b(k)**2, spray%omega_eps_g_p(k), spray%omega_eps_g_d(k)/spray%b(k)**2, spray%omega_zvar_g_p(k), spray%omega_zvar_g_d(k)/spray%b(k)**2
       end do

       close(unit=103)

#ifdef SOOT
       write(fname,"(A)") trim('Soot.out')

       inquire(file=trim(fname), exist=exist)

       if(exist) then
          open(unit=104,file=trim(fname),form="formatted",status="old",position="append",action="write")
       else
          open(unit=104,file=trim(fname),form="formatted",status="new",action="write")

          rowfmth = '(A,A,A)'

          write(104,fmt=rowfmth) '# step<1> ', 'time<2> ', 'Soot Mass (kg)<3>'

       end if

       rowfmt = "(I0.6,ES15.5E3, ES15.5E3)"

       !do s = 1,step
       s = 1
       write(104,FMT=rowfmt) spray%step, spray%time(s), spray%TotSootMass
       !end do

       close(unit=104)
#endif

    end if

    if(spray%end) then

       write(*,*) 'End of Spray Simulation'

    end if

  end subroutine write_output

  subroutine write_params(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    ! ---------------------------------
    character(len=128) :: rowfmt, rowfmth, fname, stp, tm
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, s

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    write(fname,"(A)") 'simulation_param.out'

    open(unit=99,file=trim(fname),form="formatted",status="replace",action="write")

    rowfmt = "(A,ES15.5E3)"

    write(99,FMT=rowfmt) 'Re :',spray%Re
    write(99,FMT=rowfmt) 'We :',spray%We
    write(99,FMT=rowfmt) 'DRa :',spray%DRa
    write(99,FMT=rowfmt) 'VRa :',spray%VRa
    write(99,FMT=rowfmt) 'DRv :',spray%DRv
    write(99,FMT=rowfmt) 'VRv :',spray%VRv
    write(99,FMT=rowfmt) 'WR :',spray%WR
    write(99,FMT=rowfmt) 'De :',spray%De
    write(99,FMT=*) ''
    write(99,FMT=rowfmt) 'Spray half-cone angle :',spray%theta
    write(99,FMT=rowfmt) 'Spreading coefficient :',spray%beta
    write(99,FMT=*) ''
    write(99,FMT=rowfmt) 'Discharge coefficient of nozzle :',spray%C_d
    write(99,FMT=rowfmt) 'Injection velocity :',spray%U_inj
    write(99,FMT=rowfmt) 'Effective jet diameter :',spray%D_eff
    write(99,FMT=rowfmt) 'Max. Mass flow rate(mg/ms) :',spray%rho_l*pi/4*spray%D_eff**2*spray%U_inj*1.0E+03_WP
    write(99,FMT=rowfmt) 'Max. Momentum flow rate(kgm/s^2) :',spray%rho_l*pi/4*spray%D_eff**2*spray%U_inj**2
    write(99,FMT=*) ''
    write(99,FMT=rowfmt) 'Stoichiometric mixture fraction :',spray%Zmix_st

    close(unit=99)

  end subroutine write_params

  subroutine getPenetration(spray,step)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    integer, intent(in) :: step
    ! ---------------------------------
    integer :: k, istep

    !istep = step + 1
    istep = 1

    !call add_element(spray%time,istep)
    !call add_element(spray%LPL,istep)
    !call add_element(spray%VPL,istep)

    spray%time(istep) = spray%ndtime*spray%tau*1000.0_WP   

!!$    do k = spray%kmino,spray%kmaxo
!!$       if (spray%rho(k)*spray%Y_l(k) .le. 1.0E-03_WP .and. spray%dm(k) == 0.0_WP) then
!!$          spray%LPL(istep) = spray%D_eff*(spray%z(k)-spray%z(spray%kmin)+0.5_WP*spray%dz)*1000.0_WP
!!$          exit
!!$       end if
!!$    end do
    k = spray%kmaxo
    do while(k >= spray%kmino)
       if (spray%rho(k)*spray%Y_l(k) .ge. 0.2E-03_WP) then
          spray%LPL(istep) = spray%D_eff*(spray%z(k)-spray%z(spray%kmin)+0.5_WP*spray%dz)*1000.0_WP
          spray%utip(istep) = spray%U_inj*spray%u_l(k)
          !spray%dtip(istep) = spray%D_eff*spray%dm(k)
          exit
       end if
       k = k - 1
    end do
    
    do k = spray%kmino,spray%kmaxo
      if(spray%z(k).eq.217.5_WP)then
        spray%dtip(istep) = spray%D_eff*spray%dm(k)
        !spray%utip(istep) = spray%U_inj*spray%u_l(k)
      end if
    end do

    do k = spray%kmino,spray%kmaxo
       if (spray%Y_v(k) .ge. 0.2E-03_WP) then
          spray%VPL(istep) = spray%D_eff*(spray%z(k)-spray%z(spray%kmin)+0.5_WP*spray%dz)*1000.0_WP
       end if
    end do
    
  end subroutine getPenetration

  subroutine getScalarDissipationMRIF(spray,step)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    integer, intent(in) :: step
    ! ---------------------------------
    type(srif_t), pointer     :: mrif, next_mrif
    real(WP), dimension(:), pointer :: zz=>null(), bpdf=>null()
    character(len=512) :: fname, rowfmt, rowfmth
    integer, dimension(2) :: bound
    integer :: i, k, istep, col=1, nzz, nflm
    integer, dimension(1) :: izmax
    real(WP) :: mean, var, Zmax, Zmix_max
    real(WP) :: dummy=0.0_WP, nu, z_integral, dzz, Zmix_st, PZmix_st, dz_st
    real(WP), dimension(spray%nsc) :: numer, denom

    !Zmix_st = 0.121205_WP
    !Zmix_st = 0.146836_WP
    !Zmix_st = 0.2_WP
    Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

    Zmax = maxval(spray%Zmix_g)
    izmax = maxloc(spray%Zmix_g)

    ! Warning!!! If you change anything in the following two lines,
    ! remember to change it also in spray_combust.f90
    ! Zmix_max < 1.0 results in too high scalar dissipation rate and 
    ! it becomes hard to ignite. Therefore, Zmix_max = 1.0 is used here.
    Zmix_max = maxval(spray%Zmix_g) + 2.0_WP*sqrt(maxval(spray%zvar_g))
    Zmix_max = 1.0_WP !min(1.0,Zmix_max)
    !write(*,*) izmax(1), Zmax

    numer = 0.0_WP
    denom = 0.0_WP
    do k = spray%kmino,spray%kmaxo
       nflm = 1
       nullify(mrif)
       mrif=>spray%first_rif
       do while (associated(mrif))

          zz=>mrif%zz
          bpdf=>mrif%bpdf
          if(associated(zz)) then
             nzz = size(zz)
             if(.not.associated(bpdf)) allocate(bpdf(nzz))
          else
             nzz = 1001
             allocate(zz(nzz))
             dzz = 1.0_WP/(nzz-1)
             zz = (/ (dzz*real(i,WP),i=0,nzz-1,1) /)
             if(.not.associated(bpdf)) allocate(bpdf(nzz))
          end if

          if(spray%Zmix_g(k) >= 0.0_WP .and. spray%Zvar_g(k) >= 0.0_WP) then
             ! Attention: mrif%bpdf => P(Z)*dZ
             call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),nzz,zz,bpdf,bound)

             call interpolate1val(zz,bpdf,Zmix_st,PZmix_st,dz_st)
             z_integral = sum((zz(bound(1):bound(2))/Zmix_st)**2       &
                  * (log(zz(bound(1):bound(2))/Zmix_max)/log(Zmix_st/Zmix_max)) &
                  * bpdf(bound(1):bound(2)))

             ! Limit maximum scalar dissipation rate to 1000.0 (1/s)
             !spray%chi_g_stl(k) = min(1000_WP*spray%tau,max(0.0_WP,spray%chi_g(k)/z_integral))
             spray%chi_g_stl(k) = max(0.0_WP,spray%chi_g(k)/z_integral)

             ! Integrate over spray volume
             numer(nflm) = numer(nflm) + spray%Wt(nflm,k)*spray%chi_g_stl(k)**1.5_WP/spray%DRg(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
             denom(nflm) = denom(nflm) + spray%Wt(nflm,k)*spray%chi_g_stl(k)**0.5_WP/spray%DRg(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
          end if

          nflm = nflm + 1
          next_mrif => mrif%next
          mrif => next_mrif
       end do
    end do

    ! Initially, scalar dissipation rate is set to high value
    ! which is physically possible.
    if(spray%combustion_model == 'MRIF') then
       do i = 1,spray%nsc
          if(numer(i)==0.0_WP .and. denom(i)==0.0_WP) then
             spray%chi_st_m(i) = 0.0_WP !min(1000.0_WP,1.0_WP/spray%tau)
          end if
          if(denom(i) > 0.0_WP) then
             spray%chi_st_m(i) = numer(i)/denom(i)/spray%tau
          end if
       end do
       where(spray%chi_st_m/=spray%chi_st_m .or. spray%chi_st_m < 0.0_WP) spray%chi_st_m = 0.0_WP
       if(sum(numer)==0.0_WP .and. sum(denom)==0.0_WP) then
          spray%chi_st(1) = 0.0_WP !min(1000.0_WP,1.0_WP/spray%tau)
       end if
       if(sum(denom) > 0.0_WP) then
          spray%chi_st(1) = sum(numer)/sum(denom)/spray%tau
       end if
    end if

    !if(associated(bpdf)) deallocate(bpdf)
    
  end subroutine getScalarDissipationMRIF

  subroutine getScalarDissipation(spray,step)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    integer, intent(in) :: step
    ! ---------------------------------
    real(WP), dimension(:), pointer :: zz=>null(), bpdf=>null()
    character(len=512) :: fname, rowfmt, rowfmth
    integer, dimension(2) :: bound
    integer :: i, k, istep, col=1, nzz
    integer, dimension(1) :: izmax
    real(WP) :: mean, var, Zmax, Zmix_max
    real(WP) :: dummy=0.0_WP, nu, z_integral, dzz, Zmix_st, PZmix_st, dz_st
    real(WP), dimension(:), pointer :: numer, denom

    allocate(numer(1)); numer = 0.0_WP
    allocate(denom(1)); denom = 0.0_WP

    !istep = step + 1
    istep = 1

    !call add_element(spray%chi_st,istep)

    zz=>spray%zz
    bpdf=>spray%bpdf

    if(associated(zz)) then
       nzz = size(zz)
       if(.not.associated(bpdf)) allocate(bpdf(nzz))
    else
       nzz = 1001
       allocate(zz(nzz))
       dzz = 1.0_WP/(nzz-1)
       zz = (/ (dzz*real(i,WP),i=0,nzz-1,1) /)
       if(.not.associated(bpdf)) allocate(bpdf(nzz))
    end if

    !Zmix_st = 0.121205_WP
    !Zmix_st = 0.146836_WP
    !Zmix_st = 0.2_WP
    Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

    Zmax = maxval(spray%Zmix_g)
    izmax = maxloc(spray%Zmix_g)

    ! Warning!!! If you change anything in the following two lines,
    ! remember to change it also in spray_combust.f90
    ! Zmix_max < 1.0 results in too high scalar dissipation rate and 
    ! it becomes hard to ignite. Therefore, Zmix_max = 1.0 is used here.
    Zmix_max = maxval(spray%Zmix_g) + 2.0_WP*sqrt(maxval(spray%zvar_g))
    Zmix_max = min(1.0,Zmix_max) !1.0_WP !min(1.0,Zmix_max)
    !write(*,*) izmax(1), Zmax

    if (mod(spray%step,spray%datafreq) == 0 .or. spray%saveDataFile .or. spray%end) then
       write(fname,"(A,A1,I0.6)") trim('bpdf.out'), '_', step

       open(unit=100,file=trim(fname),form="formatted",status="unknown",action="write")

       rowfmth = '*'
       rowfmt = '(<nzz>(ES15.5E3))'

       !write(100,FMT=rowfmth) ''
    end if

    numer = 0.0_WP
    denom = 0.0_WP
    do k = spray%kmino,spray%kmaxo

       if(spray%Zmix_g(k) >= 0.0_WP .and. spray%Zvar_g(k) >= 0.0_WP) then
          ! Attention: mrif%bpdf => P(Z)*dZ
          call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),nzz,zz,bpdf,bound)

          call interpolate1val(zz,bpdf,Zmix_st,PZmix_st,dz_st)
          z_integral = sum((zz(bound(1):bound(2))/Zmix_st)**2       &
                     * (log(zz(bound(1):bound(2))/Zmix_max)/log(Zmix_st/Zmix_max)) &
                     * bpdf(bound(1):bound(2)))

          if (mod(spray%step,spray%datafreq) == 0 .or. spray%saveDataFile .or. spray%end) then
             write(100, '(ES15.5E3,ES15.5E3,ES15.5E3,ES15.5E3,ES15.5E3)') spray%z(k), PZmix_st/dz_st, sum(bpdf(bound(1):bound(2))),z_integral,Zmix_max
             !write(100, '(<nzz>(ES15.5E3))') bpdf/dzz
          end if

          ! Limit maximum scalar dissipation rate to 1000.0 (1/s)
          !spray%chi_g_stl(k) = min(1000_WP*spray%tau,max(0.0_WP,spray%chi_g(k)/z_integral))
          spray%chi_g_stl(k) = max(0.0_WP,spray%chi_g(k)/z_integral)
          
          ! Integrate over spray volume
          if(spray%combustion_model == 'MRIF') then
             numer(1:spray%nsc) = numer(1:spray%nsc) + spray%Wt(1:spray%nsc,k)*spray%chi_g_stl(k)**1.5_WP/spray%DRg(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
             denom(1:spray%nsc) = denom(1:spray%nsc) + spray%Wt(1:spray%nsc,k)*spray%chi_g_stl(k)**0.5_WP/spray%DRg(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
          else
             numer(1) = numer(1) + spray%chi_g_stl(k)**1.5_WP/spray%DRg(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
             denom(1) = denom(1) + spray%chi_g_stl(k)**0.5_WP/spray%DRg(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
             !numer = numer + spray%chi_g_stl(k)**1.5_WP/spray%DRg(k)*spray%dz*Pi*spray%b(k)**2
             !denom = denom + spray%chi_g_stl(k)**0.5_WP/spray%DRg(k)*spray%dz*Pi*spray%b(k)**2
          end if
       end if

    end do
    if (mod(spray%step,spray%datafreq) == 0 .or. spray%saveDataFile .or. spray%end) then
       close(unit=100)
    end if

    ! Initially, scalar dissipation rate is set to high value
    ! which is physically possible.
    if(numer(1)==0.0_WP .and. denom(1)==0.0_WP) then
       spray%chi_st(istep) = min(1000.0_WP,1.0_WP/spray%tau)
    end if
    if(denom(1) > 0.0_WP) then
       spray%chi_st(istep) = numer(1)/denom(1)/spray%tau
!!$       if(spray%chi_st(istep) > 20000.0_WP) then
!!$          write(*,*) spray%chi_st(istep)
!!$          spray%saveDataFile = .true.
!!$          call write_output(spray,spray%step,spray%ndtime)
!!$          GO TO 1000
!!$       end if
    end if

    !if(associated(bpdf)) deallocate(bpdf)
    
  end subroutine getScalarDissipation

  subroutine getScalarDissipation1(spray,step)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    integer, intent(in) :: step
    ! ---------------------------------
    real(WP), dimension(:), pointer :: zz=>null(), bpdf=>null()
    integer, dimension(2) :: bound
    integer :: i, k, istep, col=1, nzz
    integer, dimension(1) :: izmax
    real(WP) :: mean, var, Zmax
    real(WP) :: dummy=0.0_WP, nu, z_integral, dzz, Zmix_st, PZmix_st, dz_st, numer, denom

    !istep = step + 1
    istep = 1

    !call add_element(spray%chi_st,istep)

    zz=>spray%zz
    bpdf=>spray%bpdf

    if(associated(zz)) then
       nzz = size(zz)
       if(.not.associated(bpdf)) allocate(bpdf(nzz))
    else
       nzz = 1001
       allocate(zz(nzz))
       dzz = 1.0_WP/(nzz-1)
       zz = (/ (dzz*real(i,WP),i=0,nzz-1,1) /)
       if(.not.associated(bpdf)) allocate(bpdf(nzz))
    end if

    Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

    Zmax = maxval(spray%Zmix_g)
    izmax = maxloc(spray%Zmix_g)

    numer = 0.0_WP
    denom = 0.0_WP
    do k = izmax(1),spray%kmaxo !spray%kmino,spray%kmaxo

       if(spray%Zmix_g(k) > 0.0_WP .and. spray%Zvar_g(k) > 0.0_WP) then
          ! Attention: mrif%bpdf => P(Z)*dZ
          call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),nzz,zz,bpdf,bound)
          call interpolate1val(zz,bpdf,Zmix_st,PZmix_st,dz_st)
          z_integral = sum((zz(bound(1):bound(2))/Zmix_st)**2       &
               *(log(zz(bound(1):bound(2)))/log(Zmix_st)) &
               *bpdf(bound(1):bound(2)))
          ! Limit maximum scalar dissipation rate to 1000.0 (1/s)
          spray%chi_g_stl(k) = min(1000_WP*spray%tau,max(0.0_WP,spray%chi_g(k)/z_integral))
          !spray%chi_g_stl(k) = max(0.0_WP,spray%chi_g(k)/z_integral)

          ! Integrate over spray volume
          numer = numer + spray%chi_g_stl(k)**1.5_WP*spray%rho(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
          denom = denom + spray%chi_g_stl(k)**0.5_WP*spray%rho(k)*PZmix_st/dz_st*spray%dz*Pi*spray%b(k)**2
       end if

    end do

    if(denom > 0.0_WP) then
       spray%chi_st1(istep) = numer/denom/spray%tau
!!$       if(spray%chi_st(istep) > 20000.0_WP) then
!!$          write(*,*) spray%chi_st(istep)
!!$          spray%saveDataFile = .true.
!!$          call write_output(spray,spray%step,spray%ndtime)
!!$          GO TO 1000
!!$       end if
    end if

    !if(associated(bpdf)) deallocate(bpdf)
    
  end subroutine getScalarDissipation1

!
! Random Sample from normal (Gaussian) distribution
!
  FUNCTION rand_normal(mean,stdev) RESULT(c)
    DOUBLE PRECISION :: mean,stdev,c,temp(2), r, theta
    DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462
    IF(stdev < 0.0d0) THEN

       WRITE(*,*) "Standard Deviation must be +ve"
    ELSE
       CALL RANDOM_NUMBER(temp)
       r=(-2.0d0*log(temp(1)))**0.5
       theta = 2.0d0*PI*temp(2)
       c= mean+stdev*r*sin(theta)
    END IF
  END FUNCTION rand_normal

  subroutine add_element(array,nelem)
    implicit none

    ! ---------------------------------
    real(WP), dimension(:), pointer, intent(inout) :: array
    integer, intent(in) :: nelem
    ! ---------------------------------

    real(WP), dimension(:), pointer :: array_new

    allocate(array_new(nelem)); array_new = 0.0_WP

    array_new(1:nelem-1) = array

    deallocate(array); array=>null()

    array => array_new

    nullify(array_new)

  end subroutine add_element

end module spray_func
