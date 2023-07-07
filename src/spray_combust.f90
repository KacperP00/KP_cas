!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com                                           !
!--------------------------------------------------------------------------------!

! Adapted from knock-indication routine in CIAO from Marco Davidovic
module spray_combust
   use spray_defs
   use solver_defs
   use combust_mduc_data
   use combust_mduc_interface
   use MPI

   implicit none

   !include 'mpif.h'


   private ! default property
   
   ! public members
   public spray_combust_init
   public spray_combust_poststep_mrif

   type(mduc_t), target                 :: mduc_  ! mduc solver memory
   real(WP), dimension(:), pointer      :: x  ! flamelet grid
   integer                              :: nspec
   integer                              :: nprog = 4
   integer, dimension(:), allocatable   :: iprog_species
   real(WP)                             :: PROG_THRESH 
   logical                              :: USE_KNOCK
   

contains
  
  subroutine spray_combust_init(spray)
    use parser
    implicit none

    type(spray_t), intent(inout) :: spray

    ! ----------------------------------

    character(len=str_medium)         :: mduc_input_file, dir, prefix
    type(mduc_t), pointer             :: mduc
    type(srif_t), pointer             :: mrif, next_mrif
    character(len=str_medium), dimension(:), allocatable :: species,prog_species
    real(WP)                          :: Pavg,Hmin,Hmax,nu
    integer                           :: s,c,i,count
    logical                           :: exist_dir, exist
    ! ----------------------------------

    if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') then
       write(*,*) 'No combustion model specified...'
       write(*,*) 'Computing spray without combustion...'
       return
    end if

    call read_param(spray%inp_fname,'RIF Input File',mduc_input_file,'mduc.inp')

    if(spray%combustion_model == 'MRIF') then
       call init_flamelets(spray)
    else
       allocate(mduc)
    end if

    count = 1
    mrif=>spray%first_rif
    do while (associated(mrif))
       mduc=>mrif%mduc
       !mduc => mduc_

       ! Initialize mduc solver
       call mduc_init(mduc,trim(mduc_input_file),0,spray%nranks,spray%comm)

       ! Set output directiory for flamelet data
       write(dir,"(A,A1,I0.3)") trim('rif'), '_', count
       call mduc_set_data_dir(mduc,dir)
       count = count + 1

       write(prefix,"(A)") trim('flamelet')
       call mduc_set_file_prefix(mduc,prefix)

       if(spray%irank .eq. 0) then
          inquire(directory=trim(dir), exist=exist_dir)

          if(.not. exist_dir) then
             call system('mkdir '//trim(dir))
             call system('mkdir '//trim(dir)//'/results')
          end if
       end if

       ! Check if times of ACAS and MDUC are consistent 
       if (abs(mduc_get_time(mduc,0)-spray%ndtime*spray%tau).gt.0.1_WP*spray%dt*spray%tau) then
          write(*,'(A,E15.6,A,E15.6)') 'Error: RIF and Spray are at different times. RIF: ', &
               mduc_get_time(mduc,0),'; spray: ',spray%time
       end if

       ! Progress variable
       allocate(mrif%prog_species(nprog))
       allocate(mrif%iprog_species(nprog))
       mrif%prog_species(1) = 'CO'
       mrif%prog_species(2) = 'CO2'
       mrif%prog_species(3) = 'OH'
       mrif%prog_species(4) = 'CH2O'

       ! get number of species in flamelet
       nspec = mduc_get_nnons_species(mduc)

       ! get species from flamelet 
       allocate(mrif%species(nspec))
       call mduc_get_species_names(mduc,mrif%species,len(mrif%species))

       ! find the progress variable species
       do c=1,nprog
          do s=1,nspec
             if (trim(mrif%species(s)).eq.trim(mrif%prog_species(c))) then
                mrif%iprog_species(c) = s
                write(*,*) trim(mrif%prog_species(c)), ' has index' , s
             end if
          end do
       end do

       ! allocate flamelet grid
       allocate(x(mduc%nx1))

       ! get the flamelet grid
       call mduc_get_grid_1D(mduc,x)

       !spray%zz => x
       mrif%zz => x
       ! allocate flamelet solution variables
       allocate(mrif%fY_1d(mduc%nx1,nspec))
       allocate(mrif%fT_1d(mduc%nx1))
       allocate(mrif%chi(mduc%nx1))

       if(spray%twoD .and. spray%irank .eq. 0) then
#ifdef SOOT
          allocate(mrif%sootM002D(spray%nzo,spray%nro+1)); mrif%sootM002D= 0.0_WP
          allocate(mrif%sootM102D(spray%nzo,spray%nro+1)); mrif%sootM102D = 0.0_WP
#endif
          ! get species from flamelet
          allocate(mrif%Yspec2D(nprog,spray%nzo,spray%nro+1)); mrif%Yspec2D = 0.0_WP
          allocate(mrif%Temp2D(spray%nzo,spray%nro+1)); mrif%Temp2D = 0.0_WP
       end if

#ifdef SOOT
       allocate(mrif%fS_1d(mduc%nx1,4))
       allocate(mrif%frrho_1d(mduc%nx1)); mrif%frrho_1d = 0.0_WP
       allocate(mrif%sootM00(spray%nzo)); mrif%sootM00 = 0.0_WP
       allocate(mrif%sootM10(spray%nzo)); mrif%sootM10 = 0.0_WP
#endif
       ! get species from flamelet
       allocate(mrif%Yspec(nspec,spray%nzo)); mrif%Yspec = 0.0_WP
       allocate(mrif%bpdf(mduc%nx1))
       allocate(mrif%Temp(spray%nzo)); mrif%Temp = 0.0_WP

       Hmin = 0.0_WP
       Hmax = 1.0_WP

       ! initialize flow data in flamelet solver
       call mduc_init_flow(mduc,(spray%ndtime-spray%dt)*spray%tau,spray%P_a,Hmin,Hmax)

       if (spray%step.eq.1) then
          ! overwrite mduc pressure 
          call mduc_set_pressure(mduc,-1,spray%P_a)
       end if

       ! update the time step in flamelet solver
       call mduc_update_timestep(spray%step,mduc)

       next_mrif => mrif%next
       mrif => next_mrif

    end do

    if(spray%twoD .and. spray%irank .eq. 0) then
       allocate(spray%Yspec2D(nprog,spray%nzo,spray%nro+1))
       allocate(spray%Tf2D(spray%nzo,spray%nro+1))
       allocate(spray%MassPDF(spray%first_rif%mduc%nx1))
#ifdef SOOT
       allocate(spray%SootNumDen2D(spray%nzo,spray%nro+1))
       allocate(spray%SootVolume2D(spray%nzo,spray%nro+1))
       allocate(spray%SootMass2D(spray%nzo,spray%nro+1))
#endif
    end if

    allocate(spray%Yspec(nspec,spray%nzo))
    allocate(spray%Tf(spray%nzo))
    allocate(spray%MassPDF(spray%first_rif%mduc%nx1))
#ifdef SOOT
    allocate(spray%SootNumDen(spray%nzo))
    allocate(spray%SootVolume(spray%nzo))
    allocate(spray%SootMass(spray%nzo))
#endif

  end subroutine spray_combust_init

  subroutine init_flamelets(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), intent(inout) :: spray

    ! ---------------------------------
    type(srif_t), pointer :: mrif, next_mrif
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, ierr, i, nend
    real(WP) :: max_mfr,doi,tot_inj_mass

    max_mfr = spray%rho_l*pi/4*spray%D_eff**2*spray%U_inj
    doi = spray%final_time
    nend = size(spray%roi,1)
    k = 1
    tot_inj_mass = 0.0_WP

    do k=2,nend
       tot_inj_mass = tot_inj_mass + max_mfr &
            * 0.5_WP*(spray%roi(k,2)+spray%roi(k-1,2)) &
            * (spray%roi(k,1)-spray%roi(k-1,1))*spray%tau
    end do

    spray%inj_mass = tot_inj_mass*doi/(spray%roi(nend,1)*spray%tau)
    !write(*,*) tot_inj_mass,spray%inj_mass,doi,(spray%roi(nend,1)*spray%tau)

    ! Allocate scalars for flamelets
    allocate(spray%Sc(spray%nsc,spray%nzo)); spray%Sc = 0.0_WP
    allocate(spray%Wt(spray%nsc,spray%nzo)); spray%Wt = 0.0_WP

    allocate(spray%first_rif)
    mrif => spray%first_rif
    allocate(mrif%mduc)
    do i = 1,spray%nsc-1
       call list_insert(mrif)
       next_mrif => mrif%next
       mrif => next_mrif
    end do

  end subroutine init_flamelets

  subroutine list_insert( elem )
    type(srif_t), pointer  :: elem

    type(srif_t), pointer :: next

    allocate(next)
    allocate(next%mduc)

    next%next => elem%next
    elem%next => next
  end subroutine list_insert

  subroutine spray_combust_poststep_mrif(spray)
    implicit none

    type(spray_t), intent(inout) :: spray

    ! ----------------------------------
    type(srif_t), pointer     :: mrif, next_mrif
    type(mduc_t), pointer     :: mduc
    real(WP)                  :: Pavg,Hmin,Hmax,Nu
    real(WP), dimension(:,:), pointer :: fY,fS
    real(WP), dimension(:), pointer :: fT, chi, frrho
    real(WP)                  :: Zmix_st, Zmix_max
    real(WP)                  :: Tempp, sootM00, SootNumDen, &
                                 sootM10, SootVolume, &
                                 SootMass, gasDen, gasMass, gasVol, cellVol, gasVolFrac
    integer                   :: i,c, ispec, k, nsc, nflm, inert, ierr, ii, iii
    logical                   :: reset, exist
    character(len=str_medium), dimension(:), pointer :: species
    real(WP), dimension(:), pointer :: Yspec, bpdf
    character(len=512) :: fname
    integer, dimension(2) :: bound
    real(WP) :: curr_vap_mass, factor, T_corr, dr, dA, Tfavg, rhoavg, cellVolverify
    real(WP), dimension(spray%first_rif%mduc%nx1) :: dzz 
    ! ----------------------------------

    if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') return

    !if(spray%chi_st(1) == 0.0_WP) return

    Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

    ! Warning!!! If you change anything in the following two lines,
    ! remember to change it also in getScalarDissipationRate spray_func.f90
    ! Zmix_max < 1.0 results in too high scalar dissipation rate and 
    ! it becomes hard to ignite. Therefore, Zmix_max = 1.0 is used here.
    Zmix_max = maxval(spray%Zmix_g) + 2.0_WP*sqrt(maxval(spray%zvar_g))
    Zmix_max = min(1.0,Zmix_max) !1.0_WP !min(1.0,Zmix_max)

    Pavg = spray%P_a

    ! This is for enthalpy flamelets, use Hmin=0 and Hmax=1 here
    Hmin = 0.0_WP
    Hmax = 1.0_WP

    curr_vap_mass = sum(spray%rho*spray%rho_l*spray%Y_v*spray%dz*Pi*spray%b**2*spray%D_eff**3)
    !curr_vap_mass = sum(spray%rho*spray%rho_l*(spray%Y_v+spray%Y_l)*spray%dz*Pi*spray%b**2*spray%D_eff**3)
    factor = int(curr_vap_mass/(spray%inj_mass/spray%nsc))

    !write(*,*) 'Flamelet mass: ',curr_vap_mass,factor
    spray%MassPDF = 0.0_WP
    ! Advance flamelet
    nflm = 1
    nullify(mrif)
    mrif=>spray%first_rif
    do while (associated(mrif))
       mduc=>mrif%mduc
       !mduc => mduc_
       inert = 0
       reset = .false.
       if(factor<nflm-1) then
          inert = 1
          reset = .true.
          !nflm = nflm + 1
          !next_mrif => mrif%next
          !mrif => next_mrif
          !cycle
       end if
       x => mrif%zz
       fY => mrif%fY_1d
       fT => mrif%fT_1d
       chi => mrif%chi

       if(spray%chi_st_m(nflm) == 0.0_WP .and. factor == 0) then
          inert = 1
          reset = .true.
       end if

       chi = spray%C_chi*spray%chi_st_m(nflm)*(x/Zmix_st)**2*(log(x/Zmix_max)/log(Zmix_st/Zmix_max))

       where(chi/=chi) chi = 0.0_WP
       where(chi<0.0_WP) chi = 0.0_WP

       !write(*,*) 'x, chi '

#ifdef SOOT
       fS => mrif%fS_1d
       frrho => mrif%frrho_1d
       call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,inert,reset,fY,fT,frrho,fS)

       ! Set to 1.0 if density from CATS is used
       frrho = 1.0_WP
#else
       call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,inert,reset,fY,fT)
#endif

       do k = spray%kmin,spray%kmax

          ! Attention: mrif%bpdf => P(Z)*dZ
          call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),mduc%nx1,x,mrif%bpdf,bound)

          do ii = 2,mrif%mduc%nx1-1
             dzz(ii) = 0.5_WP*(x(ii+1) - x(ii-1))
          end do
          dzz(1) = 0.5_WP*(x(2) - x(1))
          dzz(mrif%mduc%nx1) = 0.5_WP*(x(mrif%mduc%nx1) - x(mrif%mduc%nx1-1))

          ! Convolute
          do ispec = 1,nspec

             mrif%Yspec(ispec,k) = spray%Wt(nflm,k)*sum(fY(bound(1):bound(2),ispec)*mrif%bpdf(bound(1):bound(2)))

          end do

#ifdef SOOT
          mrif%sootM00(k) = spray%Wt(nflm,k)*sum(fS(bound(1):bound(2),1)*mrif%bpdf(bound(1):bound(2)))
          mrif%sootM10(k) = spray%Wt(nflm,k)*sum(fS(bound(1):bound(2),2)*mrif%bpdf(bound(1):bound(2)))
#endif

          mrif%Temp(k) = spray%Wt(nflm,k)*sum(fT(bound(1):bound(2))*mrif%bpdf(bound(1):bound(2)))

          ! compute mass distribution
          cellVol = spray%dz*Pi*spray%b(k)**2*spray%D_eff**3
          gasMass = spray%rho(k)*spray%rho_l*spray%Y_g(k)*cellVol

          spray%MassPDF = spray%MassPDF + gasMass*spray%Wt(nflm,k)*mrif%bpdf/dzz
#ifdef SOOT
          spray%TotSootMass = spray%TotSootMass + gasMass*mrif%sootM10(k)*24.0_WP
#endif
       end do

       nflm = nflm + 1
       next_mrif => mrif%next
       mrif => next_mrif
    end do

    do k = spray%kmin,spray%kmax
       nullify(mrif)
       mrif=>spray%first_rif
       ! get number of species in flamelet
       nspec = mduc_get_nnons_species(mrif%mduc)
       spray%Yspec(:,k) = 0.0_WP
       !Tempp = 0.0_WP
       spray%Tf(k) = 0.0_WP
#ifdef SOOT
       spray%SootNumDen(k) = 0.0_WP
       spray%SootVolume(k) = 0.0_WP
       spray%SootMass(k) = 0.0_WP
#endif
       do while (associated(mrif))
          mduc=>mrif%mduc
          !mduc => mduc_
          ! Sum
          do ispec = 1,nspec
             spray%Yspec(ispec,k) = spray%Yspec(ispec,k) + mrif%Yspec(ispec,k)
          end do

          !Tempp = Tempp + mrif%Temp(k)
          spray%Tf(k) = spray%Tf(k) + mrif%Temp(k)

#ifdef SOOT
          spray%SootNumDen(k) = spray%SootNumDen(k) + mrif%sootM00(k)
          spray%SootVolume(k) = spray%SootVolume(k) + mrif%sootM10(k)
#endif
          next_mrif => mrif%next
          mrif => next_mrif
       end do

#ifdef SOOT

       ! Temperature correction for gas density at high temperatures due to combustion
       ! This is essential because of incompressible assumption, the density of the mixture is constant, which must be scaled with temperature
       if(spray%Tf(k) > 0.0_WP) then
          T_corr = spray%Tg(k)*spray%T_fuel/spray%Tf(k)
       else
          T_corr = 1.0_WP
       end if

       cellVol = spray%dz*Pi*spray%b(k)**2*spray%D_eff**3

       ! TRY TODO: Use gas density from the flamelet solution for gasDen
       ! To calculate gas volume fraction gas density from CATS should be used 
       ! rather than flamelet solution to be consistent with the formulation
       ! Flamelet gas density should only be used for calculating soot moments
       !gasVolFrac = spray%rho(k)*spray%rho_l*spray%Y_g(k)/spray%rho_g(k)
       !gasVol = gasVolFrac*cellVol

       ! If density from CATS is used for soot moments calculation
       gasMass = spray%rho(k)*spray%rho_l*spray%Y_g(k)*cellVol*T_corr
       gasDen = spray%rho_g(k)*T_corr !spray%rho(k)*spray%rho_l*spray%Y_g(k)
       gasVol = gasMass/gasDen
       spray%SootNumDen(k) = spray%SootNumDen(k)*gasDen
       spray%SootVolume(k) = spray%SootVolume(k)*24.0_WP/1800.0_WP*gasVol*gasDen

       !spray%SootNumDen(k) = spray%SootNumDen(k)
       !spray%SootVolume(k) = spray%SootVolume(k)*24.0_WP/1800.0_WP*gasVol
       spray%SootMass(k) = 1800.0_WP*spray%SootVolume(k)
#endif

    end do

#ifdef SOOT
    spray%TotSootMass = sum(spray%SootMass)
#endif

    if(spray%twoD .and. spray%irank .eq. 0) then

       if (mod(spray%step,spray%datafreq) == 0 .or. spray%saveDataFile .or. spray%end) then

          nflm = 1
          nullify(mrif)
          mrif=>spray%first_rif
          do while (associated(mrif))
             mduc=>mrif%mduc
             x => mrif%zz
             fY => mrif%fY_1d
             fT => mrif%fT_1d
             fS => mrif%fS_1d

             do k = spray%kmin,spray%kmax

                do ii = 1,spray%nro+1
                   call betaPDF(spray%Zmix2D(k,ii),spray%Zvar2D(k,ii),mduc%nx1,x,mrif%bpdf,bound)

                   ! Convolute
                   do ispec = 1,nprog

                      mrif%Yspec2D(ispec,k,ii) = spray%Wt(nflm,k)*sum(fY(bound(1):bound(2),mrif%iprog_species(ispec))*mrif%bpdf(bound(1):bound(2)))

                   end do

#ifdef SOOT
                   mrif%sootM002D(k,ii) = spray%Wt(nflm,k)*sum(fS(bound(1):bound(2),1)*mrif%bpdf(bound(1):bound(2)))
                   mrif%sootM102D(k,ii) = spray%Wt(nflm,k)*sum(fS(bound(1):bound(2),2)*mrif%bpdf(bound(1):bound(2)))
#endif

                   mrif%Temp2D(k,ii) = spray%Wt(nflm,k)*sum(fT(bound(1):bound(2))*mrif%bpdf(bound(1):bound(2)))

                end do
             end do
             !end if

             nflm = nflm + 1
             next_mrif => mrif%next
             mrif => next_mrif
          end do

          do k = spray%kmin,spray%kmax
             Tfavg = 0.0_WP; rhoavg = 0.0_WP
             cellVolverify = 0.0_WP
             do ii = 1,spray%nro+1
                nullify(mrif)
                mrif=>spray%first_rif
                ! get number of species in flamelet
                nspec = mduc_get_nnons_species(mrif%mduc)

                spray%Yspec2D(:,k,ii) = 0.0_WP
                !Tempp = 0.0_WP
                spray%Tf2D(k,ii) = 0.0_WP
#ifdef SOOT
                spray%SootNumDen2D(k,ii) = 0.0_WP
                spray%SootVolume2D(k,ii) = 0.0_WP
                spray%SootMass2D(k,ii) = 0.0_WP
#endif
                ! Collec contribution from different flamelets
                do while (associated(mrif))
                   mduc=>mrif%mduc
                   !mduc => mduc_

                   ! Sum
                   do ispec = 1,nprog
                      spray%Yspec2D(ispec,k,ii) = spray%Yspec2D(ispec,k,ii) + mrif%Yspec2D(ispec,k,ii)
                   end do
                   spray%Tf2D(k,ii) = spray%Tf2D(k,ii) + mrif%Temp2D(k,ii)
#ifdef SOOT
                   spray%SootNumDen2D(k,ii) = spray%SootNumDen2D(k,ii) + mrif%sootM002D(k,ii)
                   spray%SootVolume2D(k,ii) = spray%SootVolume2D(k,ii) + mrif%sootM102D(k,ii)
#endif

                   next_mrif => mrif%next
                   mrif => next_mrif
                end do

                if(ii==1) dr = spray%r(k,ii+1)-spray%r(k,ii)
                dA = 2.0_WP*Pi*(spray%r(k,ii)+0.5_WP*dr)*dr*spray%D_eff**2
                Tfavg = Tfavg + spray%rho2D(k,ii)*spray%Tf2D(k,ii)*dA
                rhoavg = rhoavg + spray%rho2D(k,ii)*dA

#ifdef SOOT
                ! Temperature correction for gas density at high temperatures due to combustion
                ! This is essential because of incompressible assumption, the density of the inert mixture is constant,
                ! which must be scaled with temperature of flamelet
                if(spray%Tf2D(k,ii) > 0.0_WP) then
                   T_corr = spray%Tg(k)*spray%T_fuel/spray%Tf2D(k,ii)
                else
                   T_corr = 1.0_WP
                end if

                cellVol = dA*spray%dz*spray%D_eff
                !write(*,*) spray%r(k,ii), cellVol
                !cellVolverify = cellVolverify + cellVol*spray%rho2D(k,ii)*spray%rho_l

                !!! Warning: This method is not recommended as it is not consistent
                ! TRY TODO: Use gas density from the flamelet solution for gasDen
                ! To calculate gas volume fraction gas density from CATS should be used 
                ! rather than flamelet solution to be consistent with the formulation
                ! Flamelet gas density should only be used for calculating soot moments
                !gasVolFrac = spray%rho2D(k,ii)*spray%rho_l*spray%Y_g(k)/spray%rho_g(k)
                !gasVol = gasVolFrac*cellVol

                ! Scaling of density for temperature correction 
                ! Comment following five lines if density from flamelet solution is used
                !gasMass = spray%rho2D(k,ii)*spray%rho_l*spray%Y_g(k)*cellVol*T_corr
                !gasDen = spray%rho_g(k)*T_corr !spray%rho(k)*spray%rho_l*spray%Y_g(k)
                !gasVol = gasMass/gasDen

                gasVol = spray%gVF2D(k,ii)*cellVol
                gasDen = spray%rho_g2D(k,ii)*T_corr
                !write(*,*) spray%rho_g(k),gasDen,cellVol,gasVol
                !cellVolverify = cellVolverify + cellVol*spray%rho_g2D(k,ii)
                spray%SootNumDen2D(k,ii) = spray%SootNumDen2D(k,ii)*gasDen
                spray%SootVolume2D(k,ii) = spray%SootVolume2D(k,ii)*24.0_WP/1800.0_WP*gasVol*gasDen

                !spray%SootNumDen2D(k,ii) = spray%SootNumDen2D(k,ii)!*gasDen
                !spray%SootVolume2D(k,ii) = spray%SootVolume2D(k,ii)*24.0_WP/1800.0_WP*gasVol!*gasDen
                spray%SootMass2D(k,ii) = 1800.0_WP*spray%SootVolume2D(k,ii)
                spray%SootVolume2D(k,ii) = spray%SootVolume2D(k,ii)/gasVol
#endif
             end do
             !write(*,*) k, spray%b(k), spray%dz, spray%D_eff, cellVolverify, 1.545_WP**2*spray%rho_g(k)*Pi*spray%b(k)**2*spray%dz*spray%D_eff**3
          end do
       end if

    end if

    ! Update gas temperature set by update_flow_variables
    ! This is done to model the effect of high temperature during combustion on evaporation. It sets the gas-phase temperature from flamelet temperature. The blending of the temperature from flamelet solution with the inert gas temperature is performed to mimic the lower temperature experienced by the liquid phase and to avoid erroneous effects of high combustion temperatures on liquid penetration, where they should occur.
    do k = spray%kmin,spray%kmax

       if(spray%Zmix_g(k) > 0.0_WP .and. spray%Tf(k) > 0.0_WP) then
          spray%Tg(k) = spray%Tf(k)/spray%T_fuel !(1.0_WP-spray%Y_l(k))*spray%Tf(k)/spray%T_fuel + spray%Y_l(k)*spray%Tg(k)
          !spray%Tg(k) = (spray%rho(k)*spray%Y_g(k)*spray%Cp_g(k)*spray%Tf(k)/spray%T_fuel + spray%rho(k)*spray%Y_l(k)*spray%C_l(k)*spray%Tg(k))/(spray%rho(k)*spray%Y_g(k)*spray%Cp_g(k) + spray%rho(k)*spray%Y_l(k)*spray%C_l(k))
       end if

    end do

    if(spray%irank .eq. 0) then
       call write_combust_output(spray)
    end if

    return

  end subroutine spray_combust_poststep_mrif
  
  subroutine write_combust_output(spray)
    implicit none

    ! ---------------------------------
    type(spray_t), intent(inout) :: spray
    ! ---------------------------------
    type(srif_t), pointer     :: mrif, next_mrif
    type(mduc_t), pointer     :: mduc
    character(len=512) :: rowfmt, rowfmth, stp, tm
    character(len=512) :: fname, fname1, fname2, fname3, fname4
    integer, pointer :: kmin, kmax, kmino, kmaxo
    integer :: k, s, i, nsc
    logical :: exist

    kmin => spray%kmin; kmax => spray%kmax
    kmino => spray%kmino; kmaxo => spray%kmaxo

    !write(fname,"(A9,I0.6,A1,ES15.5)") 'data.out_', step, '_', time

    if (mod(spray%step,spray%datafreq) == 0 .or. spray%saveDataFile .or. spray%end) then

       write(fname,"(A,A1,I0.6)") trim('species.out'), '_', spray%step

       open(unit=100,file=trim(fname),form="formatted",status="replace",action="write")

#ifdef SOOT
       write(100,'(<nspec+5>(A))') '#z<1> ','Temperature<2> ','SootNumDen<3> ','SootVolume<4> ','SootMass<5> ',(spray%first_rif%species(i), i=1,nspec)
#else
       write(100,'(<nspec+2>(A))') '#z<1> ','Temperature<2> ',(spray%first_rif%species(i), i=1,nspec)
#endif

       do k = spray%kmin,spray%kmax

#ifdef SOOT
          write(100, '(<nspec+5>(ES15.5E3))') spray%z(k),spray%Tf(k),spray%SootNumDen(k),spray%SootVolume(k),spray%SootMass(k),(spray%Yspec(i,k), i=1,nspec)
#else
          write(100, '(<nspec+2>(ES15.5E3))') spray%z(k), spray%Tf(k), (spray%Yspec(i,k), i=1,nspec)
#endif

       end do
       close(unit=100)

       if(spray%combustion_model == 'MRIF') then
          write(fname,"(A,A1,I0.6)") trim('flamelets.out'), '_', spray%step

          open(unit=200,file=trim(fname),form="formatted",status="replace",action="write")

          nsc = spray%nsc
          !write(200,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

          do k = spray%kmin,spray%kmax
             write(200, '(<nsc+1>(ES15.5E3))') spray%z(k), (spray%Sc(i,k), i=1,nsc)
          end do
          close(unit=200)

          write(fname,"(A,A1,I0.6)") trim('weights.out'), '_', spray%step

          open(unit=300,file=trim(fname),form="formatted",status="replace",action="write")

          nsc = spray%nsc
          !write(300,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

          do k = spray%kmin,spray%kmax
             write(300, '(<nsc+1>(ES15.5E3))') spray%z(k), (spray%Wt(i,k), i=1,nsc)
          end do
          close(unit=300)

          write(fname,"(A)") trim('mass_pdf.out')

          inquire(file=trim(fname), exist=exist)

          if(exist) then
             open(unit=400,file=trim(fname),form="formatted",status="old",position="append",action="write")
          else
             open(unit=400,file=trim(fname),form="formatted",status="replace",action="write")
          end if
          nsc = size(spray%MassPDF)
          write(*,*) nsc
          !write(400,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

          write(400, '(<nsc+2>(ES15.5E3))') spray%step, spray%time(1), (spray%MassPDF(i), i=1,nsc)
          close(unit=400)
       end if


       if (spray%twoD) then

          nullify(mrif)
          mrif=>spray%first_rif

          write(fname,"(A,A1,I0.6)") trim('Temp2D.out'), '_', spray%step

          open(unit=1000,file=trim(fname),form="formatted",status="replace",action="write")

          write(1000,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro+1,(spray%xi(i), i=1,spray%nro+1)

          write(fname1,"(A,A1,I0.6)") trim('YOH_2D.out'), '_', spray%step

          open(unit=1001,file=trim(fname1),form="formatted",status="replace",action="write")

          write(1001,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro+1,(spray%xi(i), i=1,spray%nro+1)

          write(fname2,"(A,A1,I0.6)") trim('YCH2O_2D.out'), '_', spray%step

          open(unit=1002,file=trim(fname2),form="formatted",status="replace",action="write")

          write(1002,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro+1,(spray%xi(i), i=1,spray%nro+1)
#ifdef SOOT
          write(fname3,"(A,A1,I0.6)") trim('SootMass_2D.out'), '_', spray%step

          open(unit=1003,file=trim(fname3),form="formatted",status="replace",action="write")

          write(1003,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro+1,(spray%xi(i), i=1,spray%nro+1)

          write(fname4,"(A,A1,I0.6)") trim('SootVolFrac_2D.out'), '_', spray%step

          open(unit=1004,file=trim(fname4),form="formatted",status="replace",action="write")

          write(1004,'(<1>(I0.3),<spray%nro+1>(ES15.5E3))') spray%nro+1,(spray%xi(i), i=1,spray%nro+1)
#endif
          do k = spray%kmin,spray%kmax

             write(1000, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Tf2D(k,i), i=1,spray%nro+1)

!!$             write(1001, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Yspec2D(mrif%iprog_species(3),k,i), i=1,spray%nro)
!!$
!!$             write(1002, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Yspec2D(mrif%iprog_species(4),k,i), i=1,spray%nro)

             write(1001, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Yspec2D(3,k,i), i=1,spray%nro+1)

             write(1002, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%Yspec2D(4,k,i), i=1,spray%nro+1)

#ifdef SOOT
             write(1003, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%SootMass2D(k,i), i=1,spray%nro+1)
             write(1004, '(<spray%nro+2>(ES15.5E3))') spray%z(k), (spray%SootVolume2D(k,i), i=1,spray%nro+1)
#endif
          end do

          close(unit=1000)
          close(unit=1001)
          close(unit=1002)
#ifdef SOOT
          close(unit=1003)
          close(unit=1004)
#endif
       end if

    end if
  end subroutine write_combust_output

end module 
