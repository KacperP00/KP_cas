module spray_func
  use precision
  use math
  use spray_defs
  use pc_defs
  use pc_func
  use solver_defs
  use solver_func

  implicit none

contains

  subroutine init_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k
    real(WP) :: scal, Cm, Cvar, var
    
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    spray%R_gas = 8.3144598_WP   ! J/K/mol
 
    ! Initialize grid and grid variables
    call allocate_spray_grid_vars(spray)

    ! Load liquid fuel properties table if provided
    if ( spray%LFPTname /= 'noname' ) then
       call readLiquidFuelPropertiesTable(spray)
    end if

    ! Load vapor fuel properties table if provided
    if ( spray%VFPTname /= 'noname' ) then
       call readVaporFuelPropertiesTable(spray)
    end if

    if ( associated(spray%LFPT) .and. .not.associated(spray%VFPT) ) then
       write(*,*) 'Warning! Vapor Fuel Properties Table not provided.'
       write(*,*) 'Warning! Vapor Fuel Properties will be used from pc_database if available.'
    else if ( .not.associated(spray%LFPT) .and. associated(spray%VFPT) ) then
       write(*,*) 'Warning! Liquid Fuel Properties Table not provided.'
       write(*,*) 'Warning! Liquid Fuel Properties will be used from pc_database if available.'
    end if

    ! Initialize solvers
    call solver_init(spray)

    do k=kmino,kmaxo
       spray%z(k) = k*spray%dz - 3*spray%dz
    end do
    
    ! Liquid fuel properties
    if (associated(spray%LFPT)) then
       call getLiquidFuelPropertiesFromLFPT(spray)
    else
       call computeLiquidFuelProperties(spray)
    end if

    ! Ambient properties
    call computeAmbientProperties(spray)

    ! Vapor fuel properties
    if (associated(spray%VFPT)) then
       call getVaporFuelPropertiesFromVFPT(spray,spray%T_a)
    else
       call computeVaporFuelProperties(spray,spray%T_a)
    end if

    ! Compute nozzle flow
    call nozzle_flow_model(spray)

    ! Compute non-dimensional parameters
    call compute_constNonDparams(spray)

    ! Time scale
    spray%tau = spray%noz_D/spray%U_inj

    spray%ndtime = 0.0_WP; spray%ndftime = spray%final_time/spray%tau
    
    spray%step = 0

    ! Read Rate of Injection profile if provided
    call read_ROI_from_file(spray)

    ! Compute spray half cone angle and spreading coefficient
    call compute_beta(spray)

    ! Scaling parameters of initial DSD (from Spray A primary breakup simulations)
    scal = 2.2278E+09_WP; Cm = 0.009363081_WP; Cvar = 2.38132E-05_WP

    if(spray%init_dm == -9999.0_WP) then
       spray%init_dm = Cm*scal*spray%DRa/(spray%Re*spray%We)
    end if

    var = Cvar*(scal*spray%DRa/(spray%Re*spray%We))**2

    if(spray%init_d2 == -9999.0_WP) then
       spray%init_d2 = spray%init_dm**2 + var
    end if

    ! Initialize non-dimensionalized flow variables and source terms
    spray%Y_l(kmino:kmin-1) = 1.0_WP; spray%Y_l(kmin:kmaxo) = 0.0_WP
    spray%Y_v = 0.0_WP;               
    spray%Y_a(kmino:kmin-1) = 0.0_WP; spray%Y_a(kmin:kmaxo) = 1.0_WP
    spray%Y_g = spray%Y_v + spray%Y_a

    spray%rho = 1.0_WP/(spray%DRa-spray%Y_l*(spray%DRa-1.0_WP))

    spray%u_l(kmino:kmin-1) = 1.0_WP; spray%u_l(kmin:kmaxo) = 0.0_WP
    spray%u_g = 0.0_WP

    spray%d2(kmino:kmin-1) = spray%init_d2; spray%d2(kmin:kmaxo) = 0.0_WP
    spray%dm(kmino:kmin-1) = spray%init_dm; spray%dm(kmin:kmaxo) = 0.0_WP

    spray%Td = 1.0_WP
    spray%Td(kmino:kmin-1) = 1.0_WP
    spray%Tg = (spray%Y_v*(1.0_WP - spray%De) + spray%Y_a*spray%T_a/spray%T_fuel)/spray%Y_g
    spray%Tg(kmino:kmin-1) = 1.0_WP

    spray%b = 0.5_WP

    ! Reference temperature and mass fraction for evaporation model
    call computeRefTemperature(spray)

    ! Compute initial droplet size distribution
    call compute_DSD(spray)

    call computeTimeStep(spray)

    ! Post processing
    allocate(spray%time(floor(spray%ndftime/(spray%dt))+10)); spray%time = 0.0_WP
    allocate(spray%LPL(floor(spray%ndftime/(spray%dt))+10));  spray%LPL  = 0.0_WP
    allocate(spray%VPL(floor(spray%ndftime/(spray%dt))+10));  spray%VPL  = 0.0_WP

    ! Compute gas mixture properties at reference temperature
    call computeGasMixtureProperties(spray,spray%T_ref)

    ! Compute non-dimensional parameters
    call compute_varNonDparams(spray)

    spray%omega_ent  = 0.0_WP; spray%omega_vap  = 0.0_WP; spray%f_drag  = 0.0_WP;
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

    ! Output important parameters
    call write_params(spray)

  end subroutine init_spray

  subroutine run_spray(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer :: niter
    real(WP), dimension(spray%kmino:spray%kmaxo) :: T_d

    do while (spray%ndtime < spray%ndftime)
       
       write(*,*) 'step: ', spray%step, 'Time: ', spray%ndtime, 'dt: ', spray%dt!, 'End: ', spray%ndftime

       do niter = 1,spray%solver%rk%stage

          spray%solver%rk%niter = niter

          ! Compute time step
          call computeTimeStep(spray)

          ! Apply BCs
          call applyBC(spray)

          ! update models: Entrainment, Drag, Breakup, Evaporation and heating
          call entrainmentTerm(spray)

          call dragModel(spray)

          call breakupModel(spray)

          call evaporationModel(spray)

          ! Update State Vector
          call buildStateVector(spray)

          ! Update Flux Vector
          call buildFluxVector(spray)

          ! Update Source Vector
          call buildSourceVector(spray)

          ! Solve hyperbolic system of equations
          call solver_run(spray)

          ! Update flow variables
          call updateFlowVariables(spray)

          ! Compute DSD
          call compute_DSD(spray)

          ! Update reference temperature
          call computeRefTemperature(spray)

          ! Update liquid fuel properties
          T_d = spray%Td*spray%T_fuel
          if (associated(spray%LFPT)) then
             call updateLiquidFuelPropertiesFromLFPT(spray,T_d)
          else
             call updateLiquidFuelProperties(spray,T_d)
          end if

          ! Update reference properties
          if( spray%fixed_DRa > 0.0_WP .or. spray%fixed_DRv > 0.0_WP .or. &
              spray%fixed_VRa > 0.0_WP .or. spray%fixed_VRv > 0.0_WP ) then
             call computeGasMixtureProperties_fixed_DR_VR(spray,spray%T_ref)
          else
             call computeGasMixtureProperties(spray,spray%T_ref)
          end if

          ! Update varying Non-dimensional parameters
          call compute_varNonDparams(spray)

          ! Advance time 
          call advanceTime(spray)
       
       end do

       spray%step = spray%step + 1

       ! Compute Penetration Lengths
       call getPenetration(spray,spray%step)

       if (mod(spray%step,spray%outfreq) == 0) call write_output(spray,spray%step,spray%ndtime)

    end do

  end subroutine run_spray

  subroutine nozzle_flow_model(spray)
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

  end subroutine nozzle_flow_model

  subroutine compute_constNonDparams(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    ! Fuel Jet exit Reynolds number
    spray%Re = spray%rho_l*spray%U_inj*spray%noz_D/spray%visc_l
    
    ! Fuel Jet Weber number
    spray%We = spray%rho_l*spray%U_inj**2*spray%noz_D/spray%sigma                 

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
   
  end subroutine compute_varNonDparams

  subroutine compute_beta(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP) :: dd, ld, theta, DRa, beta
    
    dd = spray%noz_DoDi; ld = spray%noz_LD
    DRa = spray%DRa

    spray%Cnoz = 41.75_WP*(ld**(-0.22_WP))*(dd)**0.15_WP

    spray%theta = spray%Cnoz*(1.0_WP/DRa)**0.26_WP

    spray%beta = tan(spray%theta*Pi/180.0_WP)

  end subroutine compute_beta

  subroutine compute_DSD(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(spray%nzo):: Dbar, q, ratio, Dm, D2
    real(WP), dimension(spray%nd):: D, f
    real(WP), dimension(spray%nd,spray%nzo) :: dsd
    real(WP) :: L, h, norm, r_min, r_max
    integer :: i, k

    Dbar = 0.0_WP
    Dm = spray%dm
    D2 = spray%d2

    L = 10.0_WP*Dm(1);

    h = L/spray%nd; spray%h = h

    D = (/ (h*real(i,WP),i=0,spray%nd,1) /)

    D(1) = 1E-16

    dsd(:,:) = 0.0_WP
    
    r_min = 1.05_WP; r_max = 1.9_WP
    ratio = D2/Dm**2

    do k=spray%kmino,spray%kmaxo

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
             dsd(:,k) = f/norm;
          end if

       end if

    end do
  
    dsd(1,:) = 0.0_WP;
  
    spray%di = D
    spray%dsd = dsd

  contains
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

  end subroutine compute_DSD

  subroutine updateFlowVariables(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l
    real(WP), dimension(:,:), pointer :: W=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_a=>null(), Y_g=>null(), &
                                       u_l=>null(), u_g=>null(), d2=>null(), dm=>null(), &
                                       Td=>null(), b=>null()
    real(WP), pointer :: DRv=>null(), DRa=>null()
    real(WP), dimension(spray%nzo) :: T_low, T_high
    integer :: k

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_a => spray%Y_a; Y_g => spray%Y_g
    u_l => spray%u_l; u_g => spray%u_g; d2 => spray%d2; dm => spray%dm; Td => spray%Td; b => spray%b
    DRv => spray%DRv; DRa => spray%DRa
    W => spray%solver%W

    T_high = spray%NBP/spray%T_fuel
    T_low = spray%MP/spray%T_fuel

    rho = 0.0_WP; u_l = 0.0_WP; u_g = 0.0_WP; Y_l = 0.0_WP; Y_a = 1.0_WP; Y_v = 0.0_WP;
    dm = 0.0_WP; d2 = 0.0_WP; Td = 1.0_WP

    !u_g(:) = W(3,:)/(W(1,:)+W(2,:))

    do k = spray%kmino,spray%kmaxo

       if ( (W(1,k)+W(2,k)) .gt. 0.0_WP ) then
          u_g(k) = W(3,k)/(W(1,k)+W(2,k))
       end if

       Y_l(k) = max(0.0_WP,W(4,k)/(W(1,k)+W(2,k)+W(4,k)));
       if ( Y_l(k) < 1E-16 ) Y_l(k) = 0.0_WP
       Y_v(k) = max(0.0_WP,W(2,k)/(W(1,k)+W(2,k)+W(4,k)));
       if ( Y_v(k) < 1E-16 ) Y_v(k) = 0.0_WP
       Y_a(k) = max(0.0_WP,1.0_WP - Y_l(k) - Y_v(k));
       rho(k) = 1.0_WP/(Y_l(k) + DRv*Y_v(k) + DRa*Y_a(k))
       b(k) = sqrt((W(1,k)+W(2,k)+W(4,k))/rho(k))

       if(Y_l(k) >= 1E-16_WP) then

          dm(k) = max(0.0_WP,W(6,k)/rho(k)/Y_l(k)/b(k)**2)

          d2(k) = max(0.0_WP,W(7,k)/rho(k)/Y_l(k)/b(k)**2)

          if( dm(k) > 0.0_WP .and. d2(k) > 0.0_WP ) then

             u_l(k) = W(5,k)/rho(k)/Y_l(k)/b(k)**2
             
             Td(k) = max(T_low(k),min(T_high(k),W(8,k)/rho(k)/Y_l(k)/b(k)**2))
             
          end if

       end if

    end do

    Y_g = 1.0_WP - Y_l

    spray%Tg = (spray%Y_v*(1.0_WP - spray%De*spray%CR*spray%LR) &
             +  spray%Y_a*spray%T_a/spray%T_fuel)/spray%Y_g
    spray%Tg(spray%kmino:spray%kmin-1) = 1.0_WP

  end subroutine updateFlowVariables

  subroutine updateFlowVariablesOld(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    type(pc_t), pointer :: pc_l
    real(WP), dimension(:,:), pointer :: W=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_a=>null(), Y_g=>null(), &
                                       u_l=>null(), u_g=>null(), d2=>null(), dm=>null(), &
                                       Td=>null(), b=>null()
    real(WP), pointer :: DRv=>null(), DRa=>null()
    real(WP), dimension(spray%nzo) :: T_low, T_high
    integer :: k

    pc_l => spray%pc_l(1)
    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_a => spray%Y_a; Y_g => spray%Y_g
    u_l => spray%u_l; u_g => spray%u_g; d2 => spray%d2; dm => spray%dm; Td => spray%Td; b => spray%b
    DRv => spray%DRv; DRa => spray%DRa
    W => spray%solver%W

    T_high = pc_l%NormalBoilingPoint/spray%T_fuel
    T_low = pc_l%MeltingPoint/spray%T_fuel

    rho = 0.0_WP; u_l = 0.0_WP; u_g = 0.0_WP; Y_l = 0.0_WP; Y_a = 1.0_WP; Y_v = 0.0_WP;
    dm = 0.0_WP; d2 = 0.0_WP; Td = 1.0_WP

    !u_g(:) = W(3,:)/(W(1,:)+W(2,:))

    do k = spray%kmin,spray%kmax

       if ( (W(1,k)+W(2,k)) .gt. 0.0_WP ) then
          u_g(k) = W(3,k)/(W(1,k)+W(2,k))
       end if

       Y_l(k) = max(0.0_WP,W(4,k)/(W(1,k)+W(2,k)+W(4,k)));      
       Y_a(k) = max(0.0_WP,W(1,k)/(W(1,k)+W(2,k)+W(4,k)));

       if( W(4,k) .gt. 0.0_WP .and. Y_l(k) .gt. 0.0_WP ) then

             u_l(k) = W(5,k)/W(4,k)

             dm(k) = max(0.0_WP,W(6,k)/W(4,k))
             d2(k) = max(0.0_WP,W(7,k)/W(4,k))

             if ( dm(k) > 0.0_WP .and. d2(k) > 0.0_WP ) then
                Td(k) = max(T_low(k),min(T_high(k),W(8,k)/W(4,k)))
             end if

       end if

    end do

    Y_v = 1.0_WP - Y_l - Y_a
    rho = 1.0_WP/(Y_l + DRv*Y_v + DRa*Y_a)
    b = sqrt((W(1,:)+W(2,:)+W(4,:))/rho)

    Y_g = 1.0_WP - Y_l

  end subroutine updateFlowVariablesOld

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
       pc_l%T = min(spray%T_fuel,pc_l%NormalBoilingPoint)
       pc_l%p = spray%P_a

       ! Compute fuel properties: Liquid phase
       call computeLiqDensity(pc_l)
       call computeLiqViscosity(pc_l)
       call computeLiqHeatCapacity(pc_l)
       call computeVapPressure(pc_l)
       call computeHeatOfVap(pc_l)
       call computeSurfaceTension(pc_l)

       ! Set fuel properties to spray
       spray%C_l(k) = pc_l%liqHeatCapacity
       spray%p_vap(k) = pc_l%vapPressure%val
       spray%L_f(k) = pc_l%HeatOfVap

       nullify(pc_l)
    end do
    
    pc_l => spray%pc_l(1)

    ! Set fuel properties to spray (constant throughout spray simulation)
    spray%sigma = pc_l%SurfaceTension%val
    spray%rho_l = pc_l%liqDensity
    spray%visc_l = pc_l%liqViscosity%val
    spray%MW_f = pc_l%MolecularWeight/1000.0_WP
    spray%MP = pc_l%MeltingPoint
    spray%NBP = pc_l%NormalBoilingPoint

  end subroutine computeLiquidFuelProperties

  subroutine getLiquidFuelPropertiesFromLFPT(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:), pointer :: values
    real(WP) :: T
    integer :: col
 
    allocate(values(size(spray%LFPT,2))); values = 0.0_WP

    ! Set temperature of liquid and vapor phases on the grid
    T = max(minval(spray%LFPT(:,1)),min(spray%T_fuel,maxval(spray%LFPT(:,1))))

    ! Linearly interpolate fuel properties: Liquid phase
    ! Use temperature column as a x-value for interpolation
    col = 1
    call interpolate1(spray%LFPT,col,T,values)

    ! Set fuel properties to spray
    spray%C_l = values(8)
    spray%p_vap = values(5)
    spray%L_f = values(4)

    ! Set fuel properties to spray (constant throughout spray simulation)
    spray%sigma = values(3)
    spray%rho_l = values(7)
    spray%visc_l = values(2)

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
       write(*,*) 'Fuel molecular weight : <value>'

       spray%MP = 263.0_WP

       write(*,*) 'Setting melting point of the fuel to ', spray%MP, '...'

    end if
  
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

  subroutine updateLiquidFuelProperties(spray,T)
    implicit none

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

       pc_l => spray%pc_l(k)

       if ( pc_l%name == 'pc' ) then
          call choosePC(pc_l,spray%Fuel)
       end if

       ! Set temperature of liquid and vapor phaseson the grid
       pc_l%T = T(k)
       pc_l%p = spray%P_a

       ! Compute fuel properties: Liquid phase

       call computeLiqHeatCapacity(pc_l)
       call computeVapPressure(pc_l)
       call computeHeatOfVap(pc_l)

       ! Set fuel properties to spray
       spray%C_l(k) = pc_l%liqHeatCapacity
       spray%p_vap(k) = pc_l%vapPressure%val
       spray%L_f(k) = pc_l%HeatOfVap
       
       nullify(pc_l)
    end do
  
  end subroutine updateLiquidFuelProperties

  subroutine updateLiquidFuelPropertiesFromLFPT(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in) :: T

    ! ---------------------------------
    real(WP), dimension(:), pointer :: values
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: col, k
 
    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    allocate(values(size(spray%LFPT,2))); values = 0.0_WP
    
    do k=kmino,kmaxo

       ! Linearly interpolate fuel properties: Liquid phase
       ! Use temperature column as a x-value for interpolation
       col = 1
       call interpolate1(spray%LFPT,col,T(k),values)

       ! Set fuel properties to spray
       spray%C_l(k) = values(8)
       spray%p_vap(k) = values(5)
       spray%L_f(k) = values(4)

    end do
  
  end subroutine updateLiquidFuelPropertiesFromLFPT

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
    call computeDiffusionCoeffcientFuller(pc_v)
    
    ! Set fuel vapor properties in spray object
    spray%rho_v = pc_v%IG_vapDensity
    spray%visc_v = pc_v%vapViscosity%val
    spray%lambda_v = pc_v%vapThermalConductivity%val
    spray%Cp_v = pc_v%IG_HeatCapacity

  end subroutine computeVaporFuelProperties

  subroutine getVaporFuelPropertiesFromVFPT(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), intent(in):: T

    ! ---------------------------------
    real(WP), dimension(:), pointer :: values
    integer :: col
    real(WP) :: R_gas
    
    if (spray%MW_f == -9999.0_WP) then
       write(*,*) 'Error: Molecular weight of the fuel not provide. Please add following line to input file:'
       write(*,*) 'Fuel molecular weight : <value>'
       call abort
    end if

    R_gas = 8.3144598E03_WP   ! J/K/kmol
    
    allocate(values(size(spray%VFPT,2))); values = 0.0_WP

    ! Linearly interpolate fuel properties: Liquid phase
    ! Use temperature column as a x-value for interpolation
    col = 1
    call interpolate1(spray%VFPT,col,T,values)

    spray%rho_v = spray%P_a*spray%MW_f/(R_gas*T)
    spray%visc_v = values(2)
    spray%lambda_v = values(3)
    spray%Cp_v = values(4)
 
  end subroutine getVaporFuelPropertiesFromVFPT

  subroutine computeRefTemperature(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    ! Reference temperature and mass fraction for evaporation model
    spray%T_ref = spray%Tg*spray%T_fuel !spray%T_a !(spray%T_a + 2.0_WP*spray%T_fuel)/3.0_WP !spray%T_a ! (spray%T_a + 2.0_WP*spray%T_fuel)/3.0_WP
    !spray%Y_ref = spray%Y_v ! (spray%Y_v + 2.0_WP*spray%Y_vs)/3.0_WP

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
       call computeDiffusionCoeffcientFuller(pc_v)

       ! Set fuel vapor properties in spray object
       spray%rho_rv(k) = pc_v%IG_vapDensity
       spray%visc_rv(k) = pc_v%vapViscosity%val
       spray%lambda_rv(k) = pc_v%vapThermalConductivity%val
       spray%Cp_rv(k) = pc_v%IG_HeatCapacity
       spray%G_rv(k) = pc_v%DiffusionCoefficientFuller

       nullify(pc_v)
    end do
  end subroutine computeRefVaporFuelProperties

  subroutine getRefVaporFuelPropertiesFromVFPT(spray,T)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    real(WP), dimension(:), intent(in):: T

    ! ---------------------------------
    real(WP), dimension(:), pointer :: values
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: col, k
    real(WP) :: R_gas
    
    if (spray%MW_f == -9999.0_WP) then
       write(*,*) 'Error: Molecular weight of the fuel not provide. Please add following line to input file:'
       write(*,*) 'Fuel molecular weight : <value>'
       call abort
    end if

    R_gas = 8.3144598E03_WP   ! J/K/kmol
    
    allocate(values(size(spray%VFPT,2))); values = 0.0_WP

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo
   
    do k=kmino,kmaxo
       ! Linearly interpolate fuel properties: Liquid phase
       ! Use temperature column as a x-value for interpolation
       col = 1
       call interpolate1(spray%VFPT,col,T(k),values)

       spray%rho_rv(k) = spray%P_a*spray%MW_f/(R_gas*T(k))
       spray%visc_rv(k) = values(2)
       spray%lambda_rv(k) = values(3)
       spray%Cp_rv(k) = values(4)
       spray%G_rv(k) = values(5)
    end do

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
    spray%rho_g(spray%kmino:spray%kmin-1) = spray%rho_g(spray%kmin)

    phi_va = ((1.0_WP+sqrt(spray%visc_rv/spray%visc_ra)*(spray%WR)**0.25_WP)**2)/sqrt(8.0_WP*(1.0_WP+1.0_WP/spray%WR))
    phi_av = ((1.0_WP+sqrt(spray%visc_ra/spray%visc_rv)*(1.0_WP/spray%WR)**0.25_WP)**2)/sqrt(8.0_WP*(1.0_WP+spray%WR))

    spray%visc_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%visc_rv + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%visc_ra
    spray%visc_g(spray%kmino:spray%kmin-1) = spray%visc_g(spray%kmin)

    spray%lambda_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%lambda_rv + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%lambda_ra
    spray%lambda_g(spray%kmino:spray%kmin-1) = spray%lambda_g(spray%kmin)

    spray%Cp_g = (Y_v*spray%Cp_rv + Y_a*spray%Cp_ra)/(Y_v+Y_a)
    spray%Cp_g(spray%kmino:spray%kmin-1) = spray%Cp_g(spray%kmin)

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

    phi_va = ((1.0_WP+sqrt(spray%VRa/spray%VRv)*(spray%WR)**0.25_WP)**2)/sqrt(8.0_WP*(1.0_WP+1.0_WP/spray%WR))
    phi_av = ((1.0_WP+sqrt(spray%VRv/spray%VRa)*(1.0_WP/spray%WR)**0.25_WP)**2)/sqrt(8.0_WP*(1.0_WP+spray%WR))

    spray%visc_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%visc_l/spray%VRv + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%visc_l/spray%VRa
    spray%visc_g(spray%kmino:spray%kmin-1) = spray%visc_g(spray%kmin)

    spray%lambda_g = Y_v/(Y_v+Y_a*phi_va/spray%WR)*spray%lambda_rv + (Y_a/(Y_a+Y_v*phi_av*spray%WR))*spray%lambda_ra
    spray%lambda_g(spray%kmino:spray%kmin-1) = spray%lambda_g(spray%kmin)

    spray%Cp_g = (Y_v*spray%Cp_rv + Y_a*spray%Cp_ra)/(Y_v+Y_a)
    spray%Cp_g(spray%kmino:spray%kmin-1) = spray%Cp_g(spray%kmin)

    spray%Sc_g = spray%visc_g/(spray%rho_g*spray%G_rv)
    spray%Pr_g = spray%visc_g*spray%Cp_g/spray%lambda_g

  end subroutine computeGasMixtureProperties_fixed_DR_VR

  subroutine entrainmentTerm(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), dimension(:), pointer :: u_g=>null(), omega_ent=>null()
    real(WP), pointer :: DRa=>null()
    
    u_g => spray%u_g; omega_ent => spray%omega_ent;  DRa => spray%DRa

    omega_ent(spray%kmin:spray%kmax) = spray%beta*u_g(spray%kmin:spray%kmax)/DRa
    
  end subroutine entrainmentTerm

  subroutine evaporationModel(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: eps = 1E-16_WP
    real(WP), dimension(:,:), pointer :: dsd=>null(), Red=>null(), Shd=>null(), Nud=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), Y_v=>null(), Y_ref=>null(), &
                                       di=>null(), &
                                       DRg=>null(), VRg=>null(), Sc_g=>null(), Pr_g=>null(), &
                                       Td=>null(), CR=>null(), LR=>null()
    real(WP), pointer :: Re=>null(), h=>null(), De=>null()
    real(WP), dimension(spray%nzo) :: Lv, Lk, Xeq, T_d, Ystar_fe, Bdeq, zeta

    real(WP), dimension(spray%nd) :: zetta, D, Xneq, Ystar_f, Bd, coeff, K_vap, Qd
    real(WP) :: TBd, K_T
    integer :: k

    spray%omega_vap = 0.0_WP; spray%omega_vapdm = 0.0_WP; spray%omega_vapd2 = 0.0_WP
    spray%omega_T = 0.0_WP; spray%Shd = 0.0_WP; spray%Nud = 0.0_WP; 

    rho => spray%rho; Y_l => spray%Y_l; Y_v => spray%Y_v; Y_ref => spray%Y_ref; Td => spray%Td
    di => spray%di; h => spray%h; dsd => spray%dsd; Red => spray%Red; Shd => spray%Shd; Nud => spray%Nud
    DRg => spray%DRg; VRg => spray%VRg; Sc_g => spray%Sc_g; Pr_g => spray%Pr_g
    Re => spray%Re; CR => spray%CR; LR => spray%LR; De => spray%De

    Lv = spray%L_f; TBd = spray%NBP; T_d = spray%T_fuel*spray%Td; D = spray%noz_D*di

    Xeq = (101325.0_WP/spray%P_a)*exp((Lv*spray%MW_f/spray%R_gas)*(1.0_WP/TBd-1.0_WP/T_d))

    Lk = spray%visc_ra*sqrt(2*Pi*T_d*spray%R_gas/spray%MW_f)/(Sc_g*spray%P_a)

    !Ystar_fe = 1.0_WP/(1.0_WP + spray%WR*(spray%P_a/spray%p_vap-1.0_WP))

    Ystar_fe = Xeq/(Xeq + (1.0_WP-Xeq)*spray%WR)
    spray%Y_ref = spray%Y_v !(spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP !spray%Y_v ! (spray%Y_v + 2.0_WP*Ystar_fe)/3.0_WP
    where(Ystar_fe < 0.0_WP) Ystar_fe = 0.0_WP

    Bdeq = (Ystar_fe - Y_ref)/(1.0_WP-Ystar_fe)
    where(Bdeq < eps) Bdeq = eps

    zeta = 0.5_WP*(Pr_g/Sc_g)*log(1.0_WP+Bdeq)

    do k = spray%kmin,spray%kmax
       
       if ( spray%d2(k) > eps .and. spray%dm(k) > eps ) then

          Shd(:,k) = (2.0_WP + 0.552_WP*(Red(:,k)**0.5_WP)*(Sc_g(k)**(1.0_WP/3.0_WP)))

          zetta = zeta(k)*Shd(:,k)

          Xneq = Xeq(k) - 2.0_WP*(Lk(k)/D)*zetta
          where(Xneq <0.0_WP) Xneq = 0.0_WP

          Ystar_f = Xneq/(Xneq + (1.0_WP-Xneq)*spray%WR)
          where(Ystar_f < 0.0_WP) Ystar_f = 0.0_WP

          Bd = (Ystar_f - Y_ref(k))/(1.0_WP-Ystar_f)
          where(Bd < eps) Bd = eps

          coeff = log(1.0_WP+Bd)
          zetta =  0.5*(Pr_g(k)/Sc_g(k))*coeff*Shd(:,k);

          K_vap = 8.0_WP*coeff*Shd(:,k)/(VRg(k)*Sc_g(k))

          spray%omega_vap(k) = (1.5_WP*rho(k)*Y_l(k))*sum((h*dsd(:,k)*K_vap/di**2))/Re
          spray%omega_vapd2(k) = (1.5_WP*rho(k)*Y_l(k))*sum((h*dsd(:,k)*K_vap))/Re
          spray%omega_vapdm(k) = (1.5_WP*rho(k)*Y_l(k))*sum((h*dsd(:,k)*K_vap/di))/Re

          Nud(:,k) = (2.0_WP + 0.552_WP*(Red(:,k)**0.5_WP)*(Pr_g(k)**(1.0_WP/3.0_WP)))
          !Qd = (1.0_WP/(Pr_g(k)*VRg(k)))*(spray%Cp_g(k)/spray%C_l(k))*(zetta/(exp(zetta)-1.0_WP+eps))*(spray%T_a/spray%T_fuel-Td(k))*Nud(:,k)/di;
          Qd = (1.0_WP/(Pr_g(k)*VRg(k)))*(spray%Cp_g(k)/spray%C_l(k))*(zetta/(exp(zetta)-1.0_WP+eps))*(spray%Tg(k)-Td(k))*Nud(:,k)/di;
          K_T = sum((6.0_WP*Qd/di-1.5_WP*De*CR(k)*LR(k)*K_vap/di**2)*h*dsd(:,k))/Re

          spray%omega_T(k) = K_T*rho(k)*Y_l(k)
       end if

    end do

  end subroutine evaporationModel

  subroutine dragModel(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    real(WP), parameter :: eps = 1E-16_WP
    real(WP), dimension(:,:), pointer :: dsd=>null(), Red=>null(), CD=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), u_l=>null(), u_g=>null(), &
                                       di=>null(), DRg=>null(), VRg=>null()
    real(WP), pointer :: Re=>null(), h=>null()
    real(WP), dimension(spray%nzo) :: u_rel

    real(WP), dimension(spray%nd) :: Hvs
    real(WP) :: CDd
    integer :: k

    spray%f_drag = 0.0_WP; spray%Red = 0.0_WP; spray%CD = 0.0_WP

    rho => spray%rho; Y_l => spray%Y_l; u_l => spray%u_l; u_g => spray%u_g; 
    di => spray%di; h => spray%h; dsd => spray%dsd; Red => spray%Red
    Re => spray%Re; DRg => spray%DRg; VRg => spray%VRg; CD => spray%CD

    u_rel = abs(u_g-u_l)

    Hvs = 0.0_WP

    do k = spray%kmin,spray%kmax
       
       Red(:,k) = (Re*VRg(k)/DRg(k))*u_rel(k)*di

       Red(1,k) = Red(2,k)

       where(Red(:,k) <= 1000.0_WP) Hvs = 1.0_WP
       
       if ( spray%d2(k) > eps .and. spray%dm(k) > eps ) then

          CD(:,k) = (24.0_WP/(Red(:,k)+eps) + 4.0_WP/((Red(:,k)+eps)**(1.0_WP/3.0_WP)))*Hvs + 0.424_WP*(1.0_WP-Hvs)
          CDd = sum(CD(:,k)/di*dsd(:,k)*h)

          spray%f_drag(k) = 0.75_WP*CDd*rho(k)*Y_l(k)*u_rel(k)*(u_g(k) - u_l(k))/DRg(k)
         
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
    real(WP), dimension(:,:), pointer :: dsd=>null(), Red=>null(), CD=>null()
    real(WP), dimension(:), pointer :: rho=>null(), Y_l=>null(), u_l=>null(), u_g=>null(), &
                                       di=>null(), DRg=>null(), VRg=>null()
    real(WP), pointer :: Re=>null(), We=>null(), h=>null()
    real(WP), dimension(spray%nzo) :: u_rel, u_rel2, K_bre1, K_bre2

    real(WP), dimension(spray%nd) :: Weg, Wel, Rel, Z, T, OmegaKH, LambdaKH, tauKH, dKH, &
                                     OmegaRT, LambdaRT, tauRT, dRT, crit, tb, dst, diff
    real(WP) :: B0, B1, C3 
    integer :: k

    spray%omega_bre1 = 0.0_WP; spray%omega_bre2 = 0.0_WP

    rho => spray%rho; Y_l => spray%Y_l; u_l => spray%u_l; u_g => spray%u_g; 
    di => spray%di; h => spray%h; dsd => spray%dsd; Red => spray%Red; CD => spray%CD
    Re => spray%Re; We => spray%We; DRg => spray%DRg; VRg => spray%VRg

    u_rel = abs(u_g-u_l)

    u_rel2 = u_rel**2

    crit = 0.0_WP

    do k = spray%kmin,spray%kmax
       
       if ( spray%d2(k) > eps .and. spray%dm(k) > eps ) then
          
          if( u_rel(k) > eps ) then
             
             ! Kelvin-Helmoltz WAVE model
             B0 = 0.61; B1 = 10
             Weg = 0.5_WP*We*u_rel2(k)*di/DRg(k)
             Wel = 0.5_WP*We*di*u_rel2(k)
             Rel = 0.5_WP*Re*di*u_rel(k)
             Z = sqrt(Wel)/(Rel)
             T = Z*sqrt(Weg)
             OmegaKH = (0.34_WP+0.38_WP*Weg**1.5_WP)/(1.0_WP+Z)/(1.0_WP+1.4_WP*T**0.6_WP)*sqrt(8.0_WP/We/di**3)
             LambdaKH = (9.02_WP/2.0_WP)*di*(1.0_WP+0.45_WP*sqrt(Z))*(1.0_WP+0.4_WP*T**0.7_WP)/(1.0_WP+0.865_WP*Weg**1.67_WP)**0.6_WP
             tauKH = (3.726_WP*B1/2.0_WP)*di/(OmegaKH*LambdaKH)
             dKH = 2.0_WP*B0*LambdaKH

             ! Rayleigh-Taylor model
             C3 = 2.5_WP
             OmegaRT = 0.5_WP*(We/DRg(k))**0.25_WP*(1.0_WP/(DRg(k)+1.0_WP)**0.5_WP) &
                     * ((DRg(k)-1.0_WP)*CD(:,k))**(0.75_WP)*(u_rel2(k)/di/DRg(k))**(0.75_WP)
             LambdaRT = C3*4.0_WP*Pi*DRg(k)*sqrt(1.0_WP/We/(DRg(k)-1.0_WP))*sqrt(di/CD(:,k)/u_rel(k))
             tauRT = 1.0_WP/OmegaRT
             dRT = LambdaRT

             where(di>dRT .and. tauRT<tauKH) crit = 1.0_WP

             tb = crit*tauRT + (1.0_WP-crit)*tauKH;
             dst = crit*dRT + (1.0_WP-crit)*dKH;

          end if

          diff = di - dst
          where(diff < 0.0_WP) diff = 0.0_WP
          K_bre1(k) = sum(diff/tb*dsd(:,k)*h)
          K_bre2(k) = sum(2*di*diff/tb*dsd(:,k)*h)
         
       end if

    end do

    spray%omega_bre1 = K_bre1*rho*Y_l
    spray%omega_bre2 = K_bre2*rho*Y_l

  end subroutine breakupModel

  subroutine computeTimeStep(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    spray%CFL = spray%MaxCFL

    spray%dt = spray%CFL*spray%dz

    if (spray%step < 100) spray%dt = 0.1_WP*spray%dt

!!$    if (spray%ndtime < 1500) spray%dt = 0.1_WP*spray%dt

  end subroutine computeTimeStep

  subroutine advanceTime(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------

    spray%ndtime = spray%ndtime + spray%dt

  end subroutine advanceTime

  subroutine read_ROI_from_file(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------  
    real(WP), dimension(:,:), pointer :: scaled
    logical :: exist_file
    character(len=128) :: cmd, line
    integer :: nlines, ioerr, i, idx
    real(WP) :: max_roi

    inquire(file=trim(spray%roi_file), exist=exist_file)

    i = 1

    if(exist_file) then
       cmd = "cat "//trim(spray%roi_file)//" | sed '/^\s*#/d;/^\s*$/d' | wc -l > nlines.txt"
       call system(cmd)
       open(unit=104,file='nlines.txt')
       read(unit=104,fmt='(i)',iostat=ioerr) nlines
       cmd = 'rm nlines.txt'
       call system(cmd)
       close(unit=104)
       allocate(spray%roi(nlines,2)); spray%roi = 0.0_WP
       open(unit=105,file=trim(spray%roi_file),form="formatted",status="old",action="read")
       write(*,*) 'Reading ROI profile...'
       do while (.true.)
          read(unit=105,fmt='(a)',iostat=ioerr) line

          if (ioerr .ne. 0) then
             exit
          end if

          if (index(line,',') > 0) then
             idx = index(line,',')
             read(line(1:idx-1),*)  spray%roi(i,1)
             read(line(idx+1:len(line)),*) spray%roi(i,2)
          end if
          i = i + 1
       end do
       close(unit=105)

       allocate(scaled(nlines,2))

       max_roi = maxval(spray%roi(:,2))

       scaled(:,1) = spray%roi(:,1)/spray%tau
       scaled(:,2) = spray%roi(:,2)/max_roi

       spray%roi = scaled

       deallocate(scaled)

    end if

  end subroutine read_ROI_from_file

  ! Liquid Fuel Property Table should be in the following format
  ! Header should be commented using '#' and physical properties
  ! as a function of temperature
  ! Temperature[K] Viscosity[Ns/m^2] SurfTension[N/m] HeatOfVap[J/kg] VaporPres[Pa] ThermConductivity[W/(mK)] Density[kg/m^3] SpecificHeatCapacity[J/(kgK)]
  subroutine readLiquidFuelPropertiesTable(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------  
    logical :: exist_file
    character(len=128) :: cmd, line
    integer :: nlines, ioerr, i

    inquire(file=trim(spray%LFPTname), exist=exist_file)

    i = 1

    if(exist_file) then
       cmd = "cat "//trim(spray%LFPTname)//" | sed '/^\s*#/d;/^\s*$/d' | wc -l > nlfpt.txt"
       call system(cmd)
       open(unit=1040,file='nlfpt.txt')
       read(unit=1040,fmt='(i)',iostat=ioerr) nlines
       cmd = 'rm nlfpt.txt'
       call system(cmd)
       close(unit=1040)
       allocate(spray%LFPT(nlines,8)); spray%LFPT = 0.0_WP
       open(unit=1050,file=trim(spray%LFPTname),form="formatted",status="old",action="read")
       write(*,*) 'Reading Liquid Fuel Properties from Table...'
       do while (.true.)
          read(unit=1050,fmt='(a)',iostat=ioerr) line

          if (ioerr .ne. 0) then
             exit
          end if

          if (index(line,'!') .or. index(line,'#')) then
             cycle
          end if
          
          read(line,*)  spray%LFPT(i,1), spray%LFPT(i,2), spray%LFPT(i,3), spray%LFPT(i,4), spray%LFPT(i,5), spray%LFPT(i,6), spray%LFPT(i,7), spray%LFPT(i,8)

          write(*,*) spray%LFPT(i,1), spray%LFPT(i,2), spray%LFPT(i,3), spray%LFPT(i,4), spray%LFPT(i,5), spray%LFPT(i,6), spray%LFPT(i,7), spray%LFPT(i,8)

          i = i + 1

       end do
       close(unit=1050)
    end if

  end subroutine readLiquidFuelPropertiesTable

  ! Vapor Fuel Property Table should be in the following format
  ! Header should be commented using '#' and physical properties
  ! as a function of temperature
  ! Temperature[K] Viscosity[Ns/m^2] ThermConductivity[W/(mK)] SpecificHeatCapacity[J/(kgK)] DiffusionCoefficient[m^2/s]
  subroutine readVaporFuelPropertiesTable(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------  
    logical :: exist_file
    character(len=128) :: cmd, line
    integer :: nlines, ioerr, i

    inquire(file=trim(spray%VFPTname), exist=exist_file)

    i = 1

    if(exist_file) then
       cmd = "cat "//trim(spray%VFPTname)//" | sed '/^\s*#/d;/^\s*$/d' | wc -l > nlvfpt.txt"
       call system(cmd)
       open(unit=1041,file='nlvfpt.txt')
       read(unit=1041,fmt='(i)',iostat=ioerr) nlines
       cmd = 'rm nlvfpt.txt'
       call system(cmd)
       close(unit=1041)
       allocate(spray%VFPT(nlines,5)); spray%VFPT = 0.0_WP
       open(unit=105,file=trim(spray%VFPTname),form="formatted",status="old",action="read")
       write(*,*) 'Reading Fuel Properties from Table...'
       do while (.true.)
          read(unit=1051,fmt='(a)',iostat=ioerr) line

          if (ioerr .ne. 0) then
             exit
          end if

          if (index(line,'!') .or. index(line,'#')) then
             cycle
          end if
          
          read(line,*)  spray%VFPT(i,1), spray%VFPT(i,2), spray%VFPT(i,3), spray%VFPT(i,4), spray%VFPT(i,5)
          i = i + 1
       end do
       close(unit=1051)
    end if

  end subroutine readVaporFuelPropertiesTable

  subroutine applyBC(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray

    ! ---------------------------------
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: i
    real(WP) :: eps = 1E-16_WP

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    ! Left boundary (injector)
    spray%rho(kmino:kmin-1) = 1.0_WP

    spray%Y_l(kmino:kmin-1) = 1.0_WP
    spray%Y_v(kmino:kmin-1) = 0.0_WP
    spray%Y_a(kmino:kmin-1) = 0.0_WP
    spray%Y_g(kmino:kmin-1) = spray%Y_v(kmino:kmin-1) + spray%Y_a(kmino:kmin-1)

    spray%b(kmino:kmin-1)   = 0.5_WP

    if(allocated(spray%roi)) then
       !Interpolate
       do i=1,maxval(shape(spray%roi))-1
          if (spray%ndtime .ge. spray%roi(i,1) .and. spray%ndtime .lt. spray%roi(i+1,1)) then
             spray%u_l(kmino:kmin-1) = (spray%roi(i+1,2)-spray%roi(i,2))/(spray%roi(i+1,1)-spray%roi(i,1))* &
                     (spray%ndtime - spray%roi(i,1)) + spray%roi(i,2)
          end if
       end do
    else
       spray%u_l(kmino:kmin-1) = 1.0_WP
    end if
    !spray%u_l(kmino:kmin-1) = max((1.0_WP/ramp)*min(ramp,spray%ndtime),1.0_WP/ramp)!%1.0;%(1.0/35)*min(35,t)+(t==0)*(1/35

    spray%u_g(kmino:kmin-1) = 0.0_WP

    spray%d2(kmino:kmin-1) = spray%init_d2
    spray%dm(kmino:kmin-1) = spray%init_dm

    spray%Td(kmino:kmin-1) = 1.0_WP 

    ! Right boundary (Neumann)
    spray%rho(kmax+1:kmaxo) = spray%rho(kmax)

    spray%Y_l(kmax+1:kmaxo) = spray%Y_l(kmax)
    spray%Y_v(kmax+1:kmaxo) = spray%Y_v(kmax)
    spray%Y_a(kmax+1:kmaxo) = spray%Y_a(kmax)
    spray%Y_g(kmax+1:kmaxo) = spray%Y_g(kmax)

    spray%b(kmax+1:kmaxo) = spray%b(kmax)

    spray%u_l(kmax+1:kmaxo) = spray%u_l(kmax)
    spray%u_g(kmax+1:kmaxo) = spray%u_g(kmax)

    spray%d2(kmax+1:kmaxo) = spray%d2(kmax)
    spray%dm(kmax+1:kmaxo) = spray%dm(kmax)

    spray%Td(kmax+1:kmaxo) = spray%Td(kmax)
   
  end subroutine applyBC

  subroutine write_output(spray,step,time)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    integer, intent(in) :: step
    real(WP), intent(in) :: time
    ! ---------------------------------
    character(len=256) :: rowfmt, rowfmth, fname, stp, tm
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, s

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    !write(fname,"(A9,I0.6,A1,ES15.5)") 'data.out_', step, '_', time

    write(fname,"(A,A1,I0.6)") trim(spray%datafilename), '_', step

    open(unit=100,file=trim(fname),form="formatted",status="replace",action="write")

    !rowfmth = '(A,A,A,A,A,A,A,A,A,A,A,A)'
    rowfmt = "(ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"

    !write(100,fmt=rowfmth) 'z,', 'rho,', 'Y_l,', 'Y_v,', 'Y_a,', 'Y_g,', 'u_l,', 'u_g,', 'dm,', 'd2,', 'Td,', 'b'
    do k = kmin-1,kmax
       write(100,FMT=rowfmt) spray%z(k), spray%rho(k), spray%Y_l(k), spray%Y_v(k), spray%Y_a(k), spray%Y_g(k), spray%u_l(k), spray%u_g(k), spray%dm(k), spray%d2(k), spray%Td(k), spray%Tg(k), spray%b(k)
    end do

    close(unit=100)

    write(fname,"(A)") trim(spray%plfilename)

    open(unit=101,file=trim(fname),form="formatted",status="replace",action="write")

    rowfmth = "(ES15.5E3, ES15.5E3, ES15.5E3)"

    do s = 1,step
       write(101,FMT=rowfmth) spray%time(s), spray%LPL(s), spray%VPL(s)
    end do

    close(unit=101)

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
    write(99,FMT=rowfmt) 'Max. Mass flow rate(mg/ms) :',spray%rho_l*pi/4*spray%noz_D**2*spray%U_inj*1.0E+03_WP
    write(99,FMT=rowfmt) 'Max. Momentum flow rate(kgm/s^2) :',spray%rho_l*pi/4*spray%noz_D**2*spray%U_inj**2

    close(unit=99)

  end subroutine write_params

  subroutine getPenetration(spray,step)
    implicit none

    ! ---------------------------------
    type(spray_t), pointer, intent(inout) :: spray
    integer, intent(in) :: step
    ! ---------------------------------
    integer :: k

    spray%time(step) = spray%ndtime*spray%tau*1000.0_WP   

    do k = spray%kmino,spray%kmaxo
        if (spray%rho(k)*spray%Y_l(k) .ge. 1.0E-03_WP) then
            spray%LPL(step) = spray%noz_D*(spray%z(k)-spray%z(spray%kmin-1))*1000.0_WP
        end if
        if (spray%Y_v(k) .ge. 1.0E-03_WP) then
            spray%VPL(step) = spray%noz_D*(spray%z(k)-spray%z(spray%kmin-1))*1000.0_WP
        end if
    end do
  end subroutine getPenetration

end module spray_func
