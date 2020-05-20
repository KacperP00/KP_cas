module fpt_func
  use fpt_defs

  public :: fptInterpLiqViscosity, fptInterpSurfaceTension, fptInterpHeatOfVap, &
            fptInterpVapPressure, fptInterpLiqThermalConductivity, &
            fptInterpLiqDensity, fptInterpLiqHeatCapacity, &
            fptInterpVapViscosity, fptInterpVapThermalConductivity, &
            fptInterpVapDensity,fptInterpVapHeatCapacity, &
            fptInterpDiffusionCoeffcient, fptInterpDiffusionCoeffcientWilkeLee, &
            fptInterpDiffusionCoeffcientFuller

contains

  ! Fuel Property File should be in a specific format (in SI units)
  ! If some property data is not available, set the values to -9999.0
  ! An example for OME4: 

  ! # Fuel Properties
  ! #=======================
  ! FUEL_NAME
  ! OME4
  ! SMILES
  ! COCOCOCOCOC
  ! CHEMICAL_FORMULA
  ! C6H14O5
  ! MOLECULAR_WEIGHT
  ! 166.172E-3
  ! T_CRIT
  ! 616.581
  ! P_CRIT
  ! 2.95007
  ! V_CRIT
  ! 0.533478
  ! MELTING_POINT
  ! 256.499
  ! NORMAL_BOILING_POINT
  ! 424.45
  ! LIQUID
  ! # Temperature liqViscosity SurfaceTension HeatOfVap vapPressure liqThermalConductivity liqDensity liqHeatCapacity
  ! 0.0000000	0.0046588	0.0267094	5.55795E+07	3.58794E+00	0.2874091	1105.5191521	1609.5183970
  ! 10.0000000	0.0046588	0.0267094	5.55795E+07	3.58794E+00	0.2874091	1105.5191521	1609.5183970
  ! ...
  ! 620.0000000	0.0003053	0.0000000	0.00000E+00	4.10684E+06	0.1610231	870.3151730	1972.2997845
  ! VAPOR
  ! # Temperature vapViscosity, vapThermalConductivity, vapDensity, vapHeatCapacity, DiffusionCoefficient
  ! 2.00000E+02	2.46987E-06	5.09586E-03	4.99651E+02	1.44490E+03	-9.99900E+03
  ! ...
  ! 6.30000E+02	1.30086E-05	3.69982E-02	1.58619E+02	1.99107E+03	-9.99900E+03
  ! ...
  ! 9.90000E+02	1.30086E-05	3.69982E-02	1.00940E+02	1.99107E+03	-9.99900E+03
  ! END

  subroutine readFuelProperties(fpt)
    implicit none

    ! ---------------------------------
    type(fpt_t), pointer, intent(inout) :: fpt

    ! ---------------------------------  
    logical :: exist_file
    character(len=1024) :: cmd, line, dummy
    integer :: ioerr, i, j, iC, iH, iO
    
    iC = -1; iH = -1; iO = -1

    inquire(file=trim(fpt%file), exist=exist_file)

    j = 1

    if(exist_file) then

       fpt%nlines = 0 
       fpt%lrows = 0
       fpt%vrows = 0
       open (unit=1040, file =trim(fpt%file))
       do
          read(unit=1040,fmt=*,iostat=ioerr) line

          if (ioerr .ne. 0) then
             exit
          end if

          if(line == 'LIQUID') then
             do
                read(unit=1040,fmt=*,iostat=ioerr) line
                fpt%nlines = fpt%nlines + 1
                if(line == 'VAPOR') then
                   do
                      read(unit=1040,fmt=*,iostat=ioerr) line
                      fpt%nlines = fpt%nlines + 1
                      if(line == 'END') exit
                      fpt%vrows = fpt%vrows + 1
                   end do
                end if
                if(line == 'END') exit
                fpt%lrows = fpt%lrows + 1
             end do
          end if
          if (ioerr/=0) exit
          fpt%nlines = fpt%nlines + 1
       end do
       close(unit=1040)

       open(unit=1050,file=trim(fpt%file),form="formatted",status="old",action="read")
       write(*,*) 'Reading Liquid Fuel Properties from '//trim(fpt%file)//'...'
       do while (.true.)
          read(unit=1050,fmt='(a)',iostat=ioerr) line

          if (ioerr .ne. 0) then
             exit
          end if

          if ( index(line,'!') > 0 .or. index(line,'#') > 0 ) then
             cycle
          end if

          read(line,*) dummy
          if(trim(dummy) == 'FUEL_NAME') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) dummy
             fpt%name = dummy
          end if
          if(trim(dummy) == 'SMILES') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) dummy
             fpt%SMILES = dummy
          end if
          if(trim(dummy) == 'CHEMICAL_FORMULA') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) dummy
             iC = scan(trim(dummy),'C')
             iH = scan(trim(dummy),'H')
             iO = scan(trim(dummy),'O')
             if(iC >=0) read(dummy(iC+1:iH-1),*) fpt%ChemicalFormula%C
             if(iH >=0) read(dummy(iH+1:iO-1),*) fpt%ChemicalFormula%H
             if(iO >=0) read(dummy(iO+1:),*) fpt%ChemicalFormula%O
          end if
          if(trim(dummy) == 'MOLECULAR_WEIGHT') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) fpt%MolecularWeight
          end if
          if(trim(dummy) == 'T_CRIT') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) fpt%Tcrit
          end if
          if(trim(dummy) == 'P_CRIT') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) fpt%pcrit
          end if
          if(trim(dummy) == 'V_CRIT') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) fpt%vcrit
          end if
          if(trim(dummy) == 'MELTING_POINT') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) fpt%MeltingPoint
          end if
          if(trim(dummy) == 'NORMAL_BOILING_POINT') then
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) fpt%NormalBoilingPoint
          end if
          if(trim(dummy) == 'LIQUID') then
             allocate(fpt%lTemp(fpt%lrows)); fpt%lTemp = -9999.0
             allocate(fpt%liqViscosity(fpt%lrows)); fpt%liqViscosity = -9999.0
             allocate(fpt%SurfaceTension(fpt%lrows)); fpt%SurfaceTension = -9999.0
             allocate(fpt%HeatOfVap(fpt%lrows)); fpt%HeatOfVap = -9999.0
             allocate(fpt%vapPressure(fpt%lrows)); fpt%vapPressure = -9999.0
             allocate(fpt%liqThermalConductivity(fpt%lrows)); fpt%liqThermalConductivity = -9999.0
             allocate(fpt%liqDensity(fpt%lrows)); fpt%liqDensity = -9999.0
             allocate(fpt%liqHeatCapacity(fpt%lrows)); fpt%liqHeatCapacity = -9999.0
             i = 1
             do
                read(unit=1050,fmt='(a)',iostat=ioerr) line
                if ( index(line,'!') > 0 .or. index(line,'#') > 0 ) then
                   cycle
                end if
                read(line,*) dummy
                if(trim(dummy) == 'VAPOR') exit
                read(line,*) fpt%lTemp(i), &
                             fpt%liqViscosity(i), &
                             fpt%SurfaceTension(i), &
                             fpt%HeatOfVap(i), &
                             fpt%vapPressure(i) , &
                             fpt%liqThermalConductivity(i), &
                             fpt%liqDensity(i), &
                             fpt%liqHeatCapacity(i)

                i = i + 1
             end do
          end if
          if(trim(dummy) == 'VAPOR') then
             allocate(fpt%vTemp(fpt%vrows)); fpt%vTemp = -9999.0
             allocate(fpt%vapViscosity(fpt%vrows)); fpt%vapViscosity = -9999.0
             allocate(fpt%vapThermalConductivity(fpt%vrows)); fpt%vapThermalConductivity = -9999.0
             allocate(fpt%vapDensity(fpt%vrows)); fpt%vapDensity = -9999.0
             allocate(fpt%vapHeatCapacity(fpt%vrows)); fpt%vapHeatCapacity = -9999.0
             allocate(fpt%DiffusionCoefficient(fpt%vrows)); fpt%DiffusionCoefficient = -9999.0
             i = 1
             do
                read(unit=1050,fmt='(a)',iostat=ioerr) line
                if ( index(line,'!') > 0 .or. index(line,'#') > 0 ) then
                   cycle
                end if
                read(line,*) dummy
                if(trim(dummy) == 'END') exit
                read(line,*) fpt%vTemp(i), &
                             fpt%vapViscosity(i), &
                             fpt%vapThermalConductivity(i), &
                             fpt%vapDensity(i), &
                             fpt%vapHeatCapacity(i), &
                             fpt%DiffusionCoefficient(i)
                i = i + 1
             end do
          end if

          j = j + 1

       end do
       close(unit=1050)
    end if

  end subroutine readFuelProperties

  subroutine computeDiffusionCoeffcientWilkeLeeFPT(fpt,p,T,DiffCoeff)
    implicit none

    ! ---------------------------------------------
    type(fpt_t), pointer, intent(inout) :: fpt
    real(WP), intent(in) :: p, T
    real(WP), intent(out) :: DiffCoeff
    ! ---------------------------------------------
    real(WP) :: M_fuel, M_air, M_AB, V_b_fuel, sigma_air, sigma_fuel, sigma_AB, &
         epsilon_air, epsilon_fuel, epsilon_AB, T_star, A, B, C, D, E, F,&
         G, H, Omega_D

    M_fuel = fpt%MolecularWeight*1000.0_WP;
    M_air = 28.96_WP;
    M_AB = 2.0_WP* (((1/M_air) + (1/M_fuel))**(-1));

    ! V_b = liquid molar volume @ normal boiling point 
    ! tabled values according to Schroeder (Chaper 4):
    V_b_fuel =  7.0_WP * ( fpt%ChemicalFormula%C + &
                           fpt%ChemicalFormula%H + &
                           fpt%ChemicalFormula%O )

    sigma_air = 3.617_WP ! air /Angström
    sigma_fuel = 1.18_WP*V_b_fuel**(1/3)  ! fuel /Angström
    sigma_AB = 0.5_WP* (sigma_air + sigma_fuel)

    epsilon_air = 97.0_WP ! air /K
    epsilon_fuel = 1.15_WP*fpt%NormalBoilingPoint; ! normal boiling point /K
    epsilon_AB = (epsilon_air * epsilon_fuel)**(0.5_WP); ! /K

    T_star = T / epsilon_AB

    A = 1.06036_WP; B = 0.15610_WP; C = 0.19300_WP; D = 0.47635_WP;
    E = 1.03587_WP; F = 1.52996_WP; G = 1.76474_WP; H = 3.89411_WP;

    Omega_D = (A/(T_star**B)) + (C/(exp(D * T_star))) + &
         (E/(exp(F * T_star))) + (G/(exp(H * T_star)))

    DiffCoeff = ((3.03_WP - 0.98_WP/(M_AB**(0.5_WP)))*(1E-3_WP) &
         * T**(1.5_WP))/((p/(100000.0_WP))*(M_AB**(0.5_WP)) &
         * (sigma_AB**2)*Omega_D)/(10000.0_WP)

  end subroutine computeDiffusionCoeffcientWilkeLeeFPT

  subroutine computeDiffusionCoeffcientFullerFPT(fpt,p,T,DiffCoeff)
    implicit none

    ! ---------------------------------------------
    type(fpt_t), pointer, intent(inout) :: fpt
    real(WP), intent(in) :: p, T
    real(WP), intent(out) :: DiffCoeff
    ! ---------------------------------------------
    real(WP) :: M_fuel, M_air, M_AB, Sigma_A, Sigma_B

    M_fuel = fpt%MolecularWeight*1000.0_WP
    M_air = 28.96_WP
    M_AB = 2.0_WP* (((1/M_air) + (1/M_fuel))**(-1))
    Sigma_A = 19.7_WP; ! air, Table 11.1
    Sigma_B = 15.9_WP * fpt%ChemicalFormula%C + &
         2.31_WP * fpt%ChemicalFormula%H + &
         6.11_WP * fpt%ChemicalFormula%O

    ! aromatic rings should be considered differently!!! with -18.3; heterocyclinc rings with also -18.3

    DiffCoeff = (0.00143_WP*T**1.75_WP) / ((p/(100000.0_WP)) &
         * (M_AB**0.5_WP) * (Sigma_A**(1.0_WP/3.0_WP) &
         + Sigma_B**(1.0_WP/3.0_WP))**2)/(10000.0_WP)

  end subroutine computeDiffusionCoeffcientFullerFPT
!!$
!!$! Solid Density - kg/m^3 (kmol/m^3)
!!$subroutine fptInterpSolDensity(fpt)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%liqDensityMol%range(1,1) .and. fpt%T < fpt%liqDensityMol%range(2,1) ) then
!!$
!!$  fpt%solDensityMol%val = fptInterpVal(pc,fpt%solDensityMol%eqn,&
!!$                                       fpt%T,&
!!$                                       fpt%solDensityMol%A,&
!!$                                       fpt%solDensityMol%B,&
!!$                                       fpt%solDensityMol%C,&
!!$                                       fpt%solDensityMol%D,&
!!$                                       fpt%solDensityMol%E)
!!$  fpt%solDensity = fpt%solDensityMol%val*fpt%MolecularWeight
!!$  
!!$end subroutine fptInterpSolDensity
!!$
!!$! Liquid Density - kg/m^3 (kmol/m^3)
!!$subroutine fptInterpLiqDensity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%liqDensityMol%range(1,1) .and. fpt%T < fpt%liqDensityMol%range(2,1) ) then
!!$     fpt%liqDensityMol%val = fptInterpVal(pc,fpt%liqDensityMol%eqn,&
!!$                                       fpt%T,&
!!$                                       fpt%liqDensityMol%A,&
!!$                                       fpt%liqDensityMol%B,&
!!$                                       fpt%liqDensityMol%C,&
!!$                                       fpt%liqDensityMol%D,&
!!$                                       fpt%liqDensityMol%E)
!!$  else if (fpt%T <= fpt%liqDensityMol%range(1,1)) then
!!$     !Write(*,*) 'Warning! liqDensity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%liqDensityMol%range(1,1)
!!$     fpt%liqDensityMol%val = fpt%liqDensityMol%range(1,2)
!!$  else if (fpt%T >= fpt%liqDensityMol%range(2,1)) then
!!$     !Write(*,*) 'Warning! liqDensity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%liqDensityMol%range(2,1)
!!$     fpt%liqDensityMol%val = fpt%liqDensityMol%range(2,2)
!!$  end if
!!$
!!$  fpt%liqDensity = fpt%liqDensityMol%val*fpt%MolecularWeight
!!$
!!$end subroutine fptInterpLiqDensity
!!$
!!$! Vapor Density - kg/m^3 (kmol/m^3) Ideal Gas Law
!!$subroutine fptInterpIG_vapDensity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$  real(WP) :: R_gas
!!$
!!$  R_gas = 8.3144598E03_WP   ! J/K/kmol
!!$
!!$  fpt%IG_vapDensity = fpt%p*fpt%MolecularWeight/(R_gas*fpt%T)
!!$
!!$end subroutine fptInterpIG_vapDensity
!!$
!!$! Vapor Pressure - Pa
!!$subroutine fptInterpVapPressure(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%vapPressure%range(1,1) .and. fpt%T < fpt%vapPressure%range(2,1) ) then
!!$     fpt%vapPressure%val = fptInterpVal(pc,fpt%vapPressure%eqn,&
!!$                                        fpt%T,&
!!$                                        fpt%vapPressure%A,&
!!$                                        fpt%vapPressure%B,&
!!$                                        fpt%vapPressure%C,&
!!$                                        fpt%vapPressure%D,&
!!$                                        fpt%vapPressure%E)
!!$  else if (fpt%T <= fpt%vapPressure%range(1,1)) then
!!$     !Write(*,*) 'Warning! vapPressure out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%vapPressure%range(1,1)
!!$     fpt%vapPressure%val = fpt%vapPressure%range(1,2)
!!$  else if (fpt%T >= fpt%vapPressure%range(2,1)) then
!!$     !Write(*,*) 'Warning! vapPressure out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%vapPressure%range(2,1)
!!$     fpt%vapPressure%val = fpt%vapPressure%range(2,2)
!!$  end if
!!$
!!$end subroutine fptInterpVapPressure
!!$
!!$! Heat of Vaporization - J/kg (J/kmol)
!!$subroutine fptInterpHeatOfVap(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%HeatOfVapMol%range(1,1) .and. fpt%T < fpt%HeatOfVapMol%range(2,1) ) then
!!$     fpt%HeatOfVapMol%val = fptInterpVal(pc,fpt%HeatOfVapMol%eqn,&
!!$                                         fpt%T,&
!!$                                         fpt%HeatOfVapMol%A,&
!!$                                         fpt%HeatOfVapMol%B,&
!!$                                         fpt%HeatOfVapMol%C,&
!!$                                         fpt%HeatOfVapMol%D,&
!!$                                         fpt%HeatOfVapMol%E)
!!$  else if (fpt%T <= fpt%HeatOfVapMol%range(1,1)) then
!!$     !Write(*,*) 'Warning! HeatOfVap out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%HeatOfVapMol%range(1,1)
!!$     fpt%HeatOfVapMol%val = fpt%HeatOfVapMol%range(1,2)
!!$  else if (fpt%T >= fpt%HeatOfVapMol%range(2,1)) then
!!$     !Write(*,*) 'Warning! HeatOfVap out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%HeatOfVapMol%range(2,1)
!!$     fpt%HeatOfVapMol%val = fpt%HeatOfVapMol%range(2,2)
!!$  end if
!!$
!!$  fpt%HeatOfVap = fpt%HeatOfVapMol%val/fpt%MolecularWeight
!!$
!!$end subroutine fptInterpHeatOfVap
!!$
!!$! Solid Heat Capacity - J/kg K (J/kmol K)
!!$subroutine fptInterpSolHeatCapacity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%solHeatCapacityMol%range(1,1) .and. fpt%T < fpt%solHeatCapacityMol%range(2,1) ) then
!!$     fpt%solHeatCapacityMol%val = fptInterpVal(pc,fpt%solHeatCapacityMol%eqn,&
!!$                                               fpt%T,&
!!$                                               fpt%solHeatCapacityMol%A,&
!!$                                               fpt%solHeatCapacityMol%B,&
!!$                                               fpt%solHeatCapacityMol%C,&
!!$                                               fpt%solHeatCapacityMol%D,&
!!$                                               fpt%solHeatCapacityMol%E)
!!$  else if (fpt%T <= fpt%solHeatCapacityMol%range(1,1)) then
!!$     !Write(*,*) 'Warning! HeatCapacity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%solHeatCapacityMol%range(1,1)
!!$     fpt%solHeatCapacityMol%val = fpt%solHeatCapacityMol%range(1,2)
!!$  else if (fpt%T >= fpt%solHeatCapacityMol%range(2,1)) then
!!$     !Write(*,*) 'Warning! HeatCapacity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%solHeatCapacityMol%range(2,1)
!!$     fpt%solHeatCapacityMol%val = fpt%solHeatCapacityMol%range(2,2)
!!$  end if
!!$
!!$  fpt%solHeatCapacity = fpt%solHeatCapacityMol%val*fpt%MolecularWeight
!!$
!!$end subroutine fptInterpSolHeatCapacity
!!$
!!$! Liquid Heat Capacity - J/kg K (J/kmol K)
!!$subroutine fptInterpLiqHeatCapacity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%liqHeatCapacityMol%range(1,1) .and. fpt%T < fpt%liqHeatCapacityMol%range(2,1) ) then
!!$     fpt%liqHeatCapacityMol%val = fptInterpVal(pc,fpt%liqHeatCapacityMol%eqn,&
!!$                                               fpt%T,&
!!$                                               fpt%liqHeatCapacityMol%A,&
!!$                                               fpt%liqHeatCapacityMol%B,&
!!$                                               fpt%liqHeatCapacityMol%C,&
!!$                                               fpt%liqHeatCapacityMol%D,&
!!$                                               fpt%liqHeatCapacityMol%E)
!!$  else if (fpt%T <= fpt%liqHeatCapacityMol%range(1,1)) then
!!$     !Write(*,*) 'Warning! liqHeatCapacity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%liqHeatCapacityMol%range(1,1)
!!$     fpt%liqHeatCapacityMol%val = fpt%liqHeatCapacityMol%range(1,2)
!!$  else if (fpt%T >= fpt%liqHeatCapacityMol%range(2,1)) then
!!$     !Write(*,*) 'Warning! liqHeatCapacity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%liqHeatCapacityMol%range(2,1)
!!$     fpt%liqHeatCapacityMol%val = fpt%liqHeatCapacityMol%range(2,2)
!!$  end if
!!$
!!$  fpt%liqHeatCapacity = fpt%liqHeatCapacityMol%val/fpt%MolecularWeight
!!$
!!$end subroutine fptInterpLiqHeatCapacity
!!$
!!$! Ideal Gas heat Capacity - J/kg K (J/kmol K)
!!$subroutine fptInterpIG_HeatCapacity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%IG_HeatCapacityMol%range(1,1) .and. fpt%T < fpt%IG_HeatCapacityMol%range(2,1) ) then
!!$     fpt%IG_HeatCapacityMol%val = fptInterpVal(pc,fpt%IG_HeatCapacityMol%eqn,&
!!$                                               fpt%T,&
!!$                                               fpt%IG_HeatCapacityMol%A,&
!!$                                               fpt%IG_HeatCapacityMol%B,&
!!$                                               fpt%IG_HeatCapacityMol%C,&
!!$                                               fpt%IG_HeatCapacityMol%D,&
!!$                                               fpt%IG_HeatCapacityMol%E)
!!$  else if (fpt%T <= fpt%IG_HeatCapacityMol%range(1,1)) then
!!$     !Write(*,*) 'Warning! IG_HeatCapacity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%IG_HeatCapacityMol%range(1,1)
!!$     fpt%IG_HeatCapacityMol%val = fpt%IG_HeatCapacityMol%range(1,2)
!!$  else if (fpt%T >= fpt%IG_HeatCapacityMol%range(2,1)) then
!!$     !Write(*,*) 'Warning! IG_HeatCapacity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%IG_HeatCapacityMol%range(2,1)
!!$     fpt%IG_HeatCapacityMol%val = fpt%IG_HeatCapacityMol%range(2,2)
!!$  end if
!!$
!!$  fpt%IG_HeatCapacity = fpt%IG_HeatCapacityMol%val/fpt%MolecularWeight
!!$
!!$end subroutine fptInterpIG_HeatCapacity
!!$
!!$! Second Viral Coefficient - m^3/kg (m^3/kmol)
!!$subroutine fptInterpSeconfViralCoef(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%SecondViralCoefMol%range(1,1) .and. fpt%T < fpt%SecondViralCoefMol%range(2,1) ) then
!!$     fpt%SecondViralCoefMol%val = fptInterpVal(pc,fpt%SecondViralCoefMol%eqn,&
!!$                                               fpt%T,&
!!$                                               fpt%SecondViralCoefMol%A,&
!!$                                               fpt%SecondViralCoefMol%B,&
!!$                                               fpt%SecondViralCoefMol%C,&
!!$                                               fpt%SecondViralCoefMol%D,&
!!$                                               fpt%SecondViralCoefMol%E)
!!$  else if (fpt%T <= fpt%SecondViralCoefMol%range(1,1)) then
!!$     !Write(*,*) 'Warning! SecondViralCoef out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%SecondViralCoefMol%range(1,1)
!!$     fpt%SecondViralCoefMol%val = fpt%SecondViralCoefMol%range(1,2)
!!$  else if (fpt%T >= fpt%SecondViralCoefMol%range(2,1)) then
!!$     !Write(*,*) 'Warning! SecondViralCoef out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%SecondViralCoefMol%range(2,1)
!!$     fpt%SecondViralCoefMol%val = fpt%SecondViralCoefMol%range(2,2)
!!$  end if
!!$
!!$  fpt%SecondViralCoef = fpt%SecondViralCoefMol%val/fpt%MolecularWeight
!!$
!!$end subroutine fptInterpSeconfViralCoef
!!$
!!$! Liquid Viscosity - Pa s
!!$subroutine fptInterpLiqViscosity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%liqViscosity%range(1,1) .and. fpt%T < fpt%liqViscosity%range(2,1) ) then
!!$     fpt%liqViscosity%val = fptInterpVal(pc,fpt%liqViscosity%eqn,&
!!$                                         fpt%T,&
!!$                                         fpt%liqViscosity%A,&
!!$                                         fpt%liqViscosity%B,&
!!$                                         fpt%liqViscosity%C,&
!!$                                         fpt%liqViscosity%D,&
!!$                                         fpt%liqViscosity%E)
!!$  else if (fpt%T <= fpt%liqViscosity%range(1,1)) then
!!$     !Write(*,*) 'Warning! liqViscosity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%liqViscosity%range(1,1)
!!$     fpt%liqViscosity%val = fpt%liqViscosity%range(1,2)
!!$  else if (fpt%T >= fpt%liqViscosity%range(2,1)) then
!!$     !Write(*,*) 'Warning! liqViscosity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%liqViscosity%range(2,1)
!!$     fpt%liqViscosity%val = fpt%liqViscosity%range(2,2)
!!$  end if
!!$
!!$end subroutine fptInterpLiqViscosity
!!$
!!$! Vapor Viscosity - Pa s
!!$subroutine fptInterpVapViscosity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%vapViscosity%range(1,1) .and. fpt%T < fpt%vapViscosity%range(2,1) ) then
!!$     fpt%vapViscosity%val = fptInterpVal(pc,fpt%vapViscosity%eqn,&
!!$                                         fpt%T,&
!!$                                         fpt%vapViscosity%A,&
!!$                                         fpt%vapViscosity%B,&
!!$                                         fpt%vapViscosity%C,&
!!$                                         fpt%vapViscosity%D,&
!!$                                         fpt%vapViscosity%E)
!!$  else if (fpt%T <= fpt%vapViscosity%range(1,1)) then
!!$     !Write(*,*) 'Warning! vapViscosity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%vapViscosity%range(1,1)
!!$     fpt%vapViscosity%val = fpt%vapViscosity%range(1,2)
!!$  else if (fpt%T >= fpt%vapViscosity%range(2,1)) then
!!$     !Write(*,*) 'Warning! vapViscosity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%vapViscosity%range(2,1)
!!$     fpt%vapViscosity%val = fpt%vapViscosity%range(2,2)
!!$  end if
!!$
!!$end subroutine fptInterpVapViscosity
!!$
!!$! Liquid Thermal Conductivity - W/mK
!!$subroutine fptInterpLiqThermalConductivity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%liqThermalConductivity%range(1,1) .and. fpt%T < fpt%liqThermalConductivity%range(2,1) ) then
!!$     fpt%liqThermalConductivity%val = fptInterpVal(pc,fpt%liqThermalConductivity%eqn,&
!!$                                                   fpt%T,&
!!$                                                   fpt%liqThermalConductivity%A,&
!!$                                                   fpt%liqThermalConductivity%B,&
!!$                                                   fpt%liqThermalConductivity%C,&
!!$                                                   fpt%liqThermalConductivity%D,&
!!$                                                   fpt%liqThermalConductivity%E)
!!$  else if (fpt%T <= fpt%liqThermalConductivity%range(1,1)) then
!!$     !Write(*,*) 'Warning! liqThermalConductivity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%liqThermalConductivity%range(1,1)
!!$     fpt%liqThermalConductivity%val = fpt%liqThermalConductivity%range(1,2)
!!$  else if (fpt%T >= fpt%liqThermalConductivity%range(2,1)) then
!!$     !Write(*,*) 'Warning! liqThermalConductivity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%liqThermalConductivity%range(2,1)
!!$     fpt%liqThermalConductivity%val = fpt%liqThermalConductivity%range(2,2)
!!$  end if
!!$
!!$end subroutine fptInterpLiqThermalConductivity
!!$
!!$! Vapor Thermal Conductivity - W/mK
!!$subroutine fptInterpVapThermalConductivity(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%vapThermalConductivity%range(1,1) .and. fpt%T < fpt%vapThermalConductivity%range(2,1) ) then
!!$     fpt%vapThermalConductivity%val = fptInterpVal(pc,fpt%vapThermalConductivity%eqn,&
!!$                                                   fpt%T,&
!!$                                                   fpt%vapThermalConductivity%A,&
!!$                                                   fpt%vapThermalConductivity%B,&
!!$                                                   fpt%vapThermalConductivity%C,&
!!$                                                   fpt%vapThermalConductivity%D,&
!!$                                                   fpt%vapThermalConductivity%E)
!!$  else if (fpt%T <= fpt%vapThermalConductivity%range(1,1)) then
!!$     !Write(*,*) 'Warning! vapThermalConductivity out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%vapThermalConductivity%range(1,1)
!!$     fpt%vapThermalConductivity%val = fpt%vapThermalConductivity%range(1,2)
!!$  else if (fpt%T >= fpt%vapThermalConductivity%range(2,1)) then
!!$     !Write(*,*) 'Warning! vapThermalConductivity out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%vapThermalConductivity%range(2,1)
!!$     fpt%vapThermalConductivity%val = fpt%vapThermalConductivity%range(2,2)
!!$  end if
!!$
!!$end subroutine fptInterpVapThermalConductivity
!!$
!!$! Surface Tension - N/m
!!$subroutine fptInterpSurfaceTension(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$
!!$  if (fpt%T > fpt%SurfaceTension%range(1,1) .and. fpt%T < fpt%SurfaceTension%range(2,1) ) then
!!$     fpt%SurfaceTension%val = fptInterpVal(pc,fpt%SurfaceTension%eqn,&
!!$                                           fpt%T,&
!!$                                           fpt%SurfaceTension%A,&
!!$                                           fpt%SurfaceTension%B,&
!!$                                           fpt%SurfaceTension%C,&
!!$                                           fpt%SurfaceTension%D,&
!!$                                           fpt%SurfaceTension%E)
!!$  else if (fpt%T <= fpt%SurfaceTension%range(1,1)) then
!!$     !Write(*,*) 'Warning! SurfaceTension out of range. Setting it for minimum of the range...'
!!$     !fpt%T = fpt%SurfaceTension%range(1,1)
!!$     fpt%SurfaceTension%val = fpt%SurfaceTension%range(1,2)
!!$  else if (fpt%T >= fpt%SurfaceTension%range(2,1)) then
!!$     !Write(*,*) 'Warning! SurfaceTension out of range. Setting it for maximum of the range...'
!!$     !fpt%T = fpt%SurfaceTension%range(2,1)
!!$     fpt%SurfaceTension%val = fpt%SurfaceTension%range(2,2)
!!$  end if
!!$
!!$end subroutine fptInterpSurfaceTension
!!$
!!$! Diffusion Coefficients - m^2/s
!!$! According to J. Poling, The Properties of Gases and Liquids, 5th edition
!!$
!!$subroutine fptInterpDiffusionCoeffcientWilkeLee(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$  real(WP) :: M_fuel, M_air, M_AB, V_b_fuel, sigma_air, sigma_fuel, sigma_AB, &
!!$              epsilon_air, epsilon_fuel, epsilon_AB, T_star, A, B, C, D, E, F,&
!!$              G, H, Omega_D
!!$
!!$  M_fuel = fpt%MolecularWeight;
!!$  M_air = 28.96_WP;
!!$  M_AB = 2.0_WP* (((1/M_air) + (1/M_fuel))**(-1));
!!$  
!!$  ! V_b = liquid molar volume @ normal boiling point 
!!$  ! tabled values according to Schroeder (Chaper 4):
!!$  V_b_fuel =  7.0_WP * ( fpt%ChemicalFormula%C + &
!!$                         fpt%ChemicalFormula%H + &
!!$                         fpt%ChemicalFormula%O )
!!$
!!$  sigma_air = 3.617_WP ! air /Angström
!!$  sigma_fuel = 1.18_WP*V_b_fuel**(1/3)  ! fuel /Angström
!!$  sigma_AB = 0.5_WP* (sigma_air + sigma_fuel)
!!$  
!!$  epsilon_air = 97.0_WP ! air /K
!!$  epsilon_fuel = 1.15_WP*fpt%NormalBoilingPoint; ! normal boiling point /K
!!$  epsilon_AB = (epsilon_air * epsilon_fuel)**(0.5_WP); ! /K
!!$
!!$  T_star = fpt%T / epsilon_AB
!!$
!!$  A = 1.06036_WP; B = 0.15610_WP; C = 0.19300_WP; D = 0.47635_WP;
!!$  E = 1.03587_WP; F = 1.52996_WP; G = 1.76474_WP; H = 3.89411_WP;
!!$
!!$  Omega_D = (A/(T_star**B)) + (C/(exp(D * T_star))) + &
!!$            (E/(exp(F * T_star))) + (G/(exp(H * T_star)))
!!$
!!$  fpt%DiffusionCoefficientWilkeLee = ((3.03_WP - 0.98_WP/(M_AB**(0.5_WP)))*(1E-3_WP) &
!!$                                  * fpt%T**(1.5_WP))/((fpt%p/(100000.0_WP))*(M_AB**(0.5_WP)) &
!!$                                  * (sigma_AB**2)*Omega_D)/(10000.0_WP)
!!$
!!$end subroutine fptInterpDiffusionCoeffcientWilkeLee
!!$
!!$subroutine fptInterpDiffusionCoeffcientFuller(pc)
!!$  implicit none
!!$
!!$  ! ---------------------------------------------
!!$  type(fpt_t), pointer, intent(inout) :: fpt
!!$  ! ---------------------------------------------
!!$  real(WP) :: M_fuel, M_air, M_AB, Sigma_A, Sigma_B
!!$
!!$  M_fuel = fpt%MolecularWeight
!!$  M_air = 28.96_WP
!!$  M_AB = 2.0_WP* (((1/M_air) + (1/M_fuel))**(-1))
!!$  Sigma_A = 19.7_WP; ! air, Table 11.1
!!$  Sigma_B = 15.9_WP * fpt%ChemicalFormula%C + &
!!$            2.31_WP * fpt%ChemicalFormula%H + &
!!$            6.11_WP * fpt%ChemicalFormula%O
!!$ 
!!$  ! aromatic rings should be considered differently!!! with -18.3; heterocyclinc rings with also -18.3
!!$              
!!$  fpt%DiffusionCoefficientFuller = (0.00143_WP*fpt%T**1.75_WP) / ((fpt%p/(100000.0_WP)) &
!!$                                * (M_AB**0.5_WP) * (Sigma_A**(1.0_WP/3.0_WP) &
!!$                                + Sigma_B**(1.0_WP/3.0_WP))**2)/(10000.0_WP)
!!$ 
!!$end subroutine fptInterpDiffusionCoeffcientFuller
!!$
!!$real(WP) function fptInterpVal(pc,eqn,T,A,B,C,D,E)
!!$ implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  integer, intent(in) :: eqn
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$
!!$  select case (eqn)
!!$  case (100)
!!$     fptInterpVal = eqn100(pc,T,A,B,C,D,E)
!!$  case (101)
!!$     fptInterpVal = eqn101(pc,T,A,B,C,D,E)
!!$  case (102)
!!$     fptInterpVal = eqn102(pc,T,A,B,C,D,E)
!!$  case (104)
!!$     fptInterpVal = eqn104(pc,T,A,B,C,D,E)
!!$  case (105)
!!$     fptInterpVal = eqn105(pc,T,A,B,C,D,E)
!!$  case (106)
!!$     fptInterpVal = eqn106(pc,T,A,B,C,D,E)
!!$  case (107)
!!$     fptInterpVal = eqn107(pc,T,A,B,C,D,E)
!!$  case (114)
!!$     fptInterpVal = eqn114(pc,T,A,B,C,D,E)
!!$  case (115)
!!$     fptInterpVal = eqn115(pc,T,A,B,C,D,E)
!!$  case default 
!!$     write(*,*) 'No equation available to fptInterp property.'
!!$  end select
!!$
!!$end function fptInterpVal
!!$
!!$real(WP) function eqn100(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$
!!$  eqn100 = A + B*T + C*T**2 + D*T**3 + E*T**4
!!$
!!$end function eqn100
!!$
!!$real(WP) function eqn101(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$
!!$  eqn101 = exp(A + B/T + C*log(T) + D*(T**E))
!!$
!!$end function eqn101
!!$
!!$real(WP) function eqn102(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$
!!$  eqn102 = (A * T**B) / (1 + C/T + D/(T**2))
!!$
!!$end function eqn102
!!$
!!$real(WP) function eqn104(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$
!!$  eqn104 = A + (B/T) + (C/((T)**3)) + (D/((T)**8)) + (E/((T)**9))
!!$
!!$end function eqn104
!!$
!!$real(WP) function eqn105(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$
!!$  eqn105 = A/(B**(1+(1-T/C)**D))
!!$
!!$end function eqn105
!!$
!!$real(WP) function eqn106(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$  real(WP) :: Tred
!!$  
!!$
!!$  Tred = T / fpt%Tcrit;
!!$  eqn106 = A * (1 - Tred)**(B + C*Tred + D*Tred**2 + E*Tred**3)
!!$
!!$end function eqn106
!!$
!!$real(WP) function eqn107(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$  real(WP) :: n, m
!!$
!!$  n = C/T
!!$  m = E/T
!!$  eqn107 = A + B*(n/sinh(n))**2 + D*(m/cosh(m))**2;
!!$
!!$end function eqn107
!!$
!!$real(WP) function eqn114(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$  real(WP) :: n, m, Tred, t1
!!$
!!$  Tred = T / fpt%Tcrit
!!$  t1 = (1.0_WP-Tred)
!!$  eqn114 = (A**2)/t1+B-2.0_WP*A*C*t1-2.0_WP*A*D*(t1**2)-(C**2)*(t1**3)/3.0_WP-C*D*(t1**4)/2.0_WP-(D**2)*(t1**5)/5.0_WP
!!$
!!$end function eqn114
!!$
!!$real(WP) function eqn115(pc,T,A,B,C,D,E)
!!$  implicit none
!!$  ! --------------------------------------------
!!$  type(fpt_t), pointer, intent(in) :: pc
!!$  real(WP), intent(in) :: T, A, B, C, D, E
!!$  ! --------------------------------------------
!!$  real(WP) :: t1
!!$
!!$  t1 = T/1000.0_WP
!!$
!!$  eqn115 = A + B*t1 + C*t1**2 + D*t1**3 + E*t1**(-2)
!!$
!!$end function eqn115
!!$
!!$end module pc_func

end module fpt_func
