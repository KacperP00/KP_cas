!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com                                           !
!--------------------------------------------------------------------------------!

module fpt_func
  use fpt_defs

  public :: readFuelProperties, computeDiffusionCoefficientWilkeLeeFPT, &
            computeDiffusionCoefficientFullerFPT, computeDiffusionCoefficientTheoryFPT

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
             fpt%ChemicalFormula%C = 0
             fpt%ChemicalFormula%H = 0
             fpt%ChemicalFormula%O = 0
             read(unit=1050,fmt='(a)',iostat=ioerr) line
             read(line,*) dummy
             iC = scan(trim(dummy),'C')
             iH = scan(trim(dummy),'H')
             iO = scan(trim(dummy),'O')
             if(iC >=0) read(dummy(iC+1:iH-1),*) fpt%ChemicalFormula%C
             if(iH > 0 .and. iO > 0) then
                read(dummy(iH+1:iO-1),*) fpt%ChemicalFormula%H
             else
                read(dummy(iH+1:),*) fpt%ChemicalFormula%H
             end if
             if(iO > 0) read(dummy(iO+1:),*) fpt%ChemicalFormula%O
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

  subroutine computeDiffusionCoefficientWilkeLeeFPT(fpt,p,T,DiffCoeff)
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
    M_AB = 2.0_WP/((1/M_air) + (1/M_fuel));

    ! V_b = liquid molar volume @ normal boiling point 
    ! tabled values according to Schroeder (Chaper 4):
    V_b_fuel =  7.0_WP * ( fpt%ChemicalFormula%C + &
                           fpt%ChemicalFormula%H + &
                           fpt%ChemicalFormula%O )

    sigma_air = 3.617_WP ! air /Angström
    sigma_fuel = 1.18_WP*V_b_fuel**(1.0_WP/3.0_WP)  ! fuel /Angström
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

  end subroutine computeDiffusionCoefficientWilkeLeeFPT

  subroutine computeDiffusionCoefficientFullerFPT(fpt,p,T,DiffCoeff)
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

  end subroutine computeDiffusionCoefficientFullerFPT

  subroutine computeDiffusionCoefficientTheoryFPT(fpt,p,T,DiffCoeff)
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
    M_AB = 2.0_WP/((1/M_air) + (1/M_fuel));

    ! V_b = liquid molar volume @ normal boiling point 
    ! tabled values according to Schroeder (Chaper 4):
    V_b_fuel =  7.0_WP * ( fpt%ChemicalFormula%C + &
                           fpt%ChemicalFormula%H + &
                           fpt%ChemicalFormula%O )

    sigma_air = 3.617_WP ! air /Angström
    sigma_fuel = 1.18_WP*V_b_fuel**(1.0_WP/3.0_WP)  ! fuel /Angström
    sigma_AB = 0.5_WP* (sigma_air + sigma_fuel)

    epsilon_air = 97.0_WP ! air /K
    epsilon_fuel = 1.15_WP*fpt%NormalBoilingPoint; ! normal boiling point /K
    epsilon_AB = (epsilon_air * epsilon_fuel)**(0.5_WP); ! /K

    T_star = T / epsilon_AB

    A = 1.06036_WP; B = 0.15610_WP; C = 0.19300_WP; D = 0.47635_WP;
    E = 1.03587_WP; F = 1.52996_WP; G = 1.76474_WP; H = 3.89411_WP;

    Omega_D = (A/(T_star**B)) + (C/(exp(D * T_star))) + &
         (E/(exp(F * T_star))) + (G/(exp(H * T_star)))

    DiffCoeff = 0.00266_WP * T**(1.5_WP)/((p/(100000.0_WP))*(M_AB**(0.5_WP)) &
         * (sigma_AB**2)*Omega_D)/(10000.0_WP)

  end subroutine computeDiffusionCoefficientTheoryFPT

end module fpt_func
