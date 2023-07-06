module pc_database
  use precision
  use pc_defs

  implicit none

contains

  ! XXXXXXX NOT FULLY UPDATED XXXXXXXXX
  ! Source: Thermodynamic and transport properties of diesel fuel from book 
  ! Multiphase Flow Dynamics 4 by Nikolay Ivanove Kolev
  ! All parameters and coefficents are not available, therefore those not
  ! available are taken from ndodecane
  subroutine pc_diesel(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

    pc%name = 'diesel'
    pc%CAS = ''
    ! components %%%%% Arbitrary %%%%%
    pc%ChemicalFormula%C = 12;         
    pc%ChemicalFormula%H = 26;
    pc%ChemicalFormula%O = 0;


    pc%MolecularWeight = 170.000_WP;   ! kg/kmol  varies between 148 and 200
    pc%Tcrit = 658.40_WP;              ! K    
    pc%pcrit = 3.0000E6_WP;            ! Pa 
    pc%vcrit = 0.71600_WP;             ! m^3/kmol retained same as ndodecane
    pc%CritCompressFactor = 0.238_WP;  ! -        retained same as ndodecane


    pc%MeltingPoint = 263.57_WP;             ! K  retained same as ndodecane
    pc%TripletPointT = 263.57_WP;            ! K  retained same as ndodecane
    pc%TripletPointp = 6.1520E-01_WP;        ! Pa retained same as ndodecane
    pc%NormalBoilingPoint = 553.15_WP;       ! K  varies from 453.15 K - 653.15 
    pc%liqMolarVolume  = 0.228605_WP;        ! m^3/kmol retained same as ndodecane


    pc%IG_HeatOfFormation = -2.9072E08_WP;            ! J/kmol retained same as ndodecane
    pc%IG_GibbsOfFormation = 3.3180E07_WP;            ! J/kmol retained same as ndodecane
    pc%IG_AbsoluteEntropy = 6.2415E5_WP;              ! J/kmolK retained same as ndodecane
    pc%HeatFusionAtMeltPoint = 3.6840E07_WP;          ! J/kmol retained same as ndodecane
    pc%StandardNetHeatOfCombustion = -7.5137E9_WP;    ! J/kmol retained same as ndodecane
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol retained same as ndodecane
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol retained same as ndodecane
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK retained same as ndodecane


    pc%AcentricFactor = 0.5764_WP;            ! - retained same as ndodecane
    pc%RadiusOfGyration = 5.9140E-10_WP;      ! m retained same as ndodecane
    pc%SolubilityParameter = 1.5900E4_WP;     ! (J/m^3)^.5 retained same as ndodecane
    pc%DipoleMoment = 0.0_WP;                   ! C^m retained same as ndodecane
    pc%VanDerWaalsVolume = 0.12960_WP;        ! m^3/kmol retained same as ndodecane


    pc%VanDerWaalsArea = 1.7740E09_WP;        ! m^2/kmol retained same as ndodecane
    pc%RefractiveIndex = 1.4151_WP;           ! - retained same as ndodecane
    pc%FlashPoint = 347.00_WP;                ! K retained same as ndodecane
    pc%FlammabilityLimits = (/ 0.6_WP, 4.7_WP /);     ! vol-! retained same as ndodecane
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K retained same as ndodecane
    pc%AutoignitionTemp = 476.00_WP;          ! K retained same as ndodecane


    !Retained same as ndodecane
    pc%solDensityMol%A = 5.2800E00_WP;
    pc%solDensityMol%B = 0.0_WP; 
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 263.57_WP, MAX_REAL_WP /)
    pc%solDensityMol%range(:,2) = (/ 5.2800E00_WP, MAX_REAL_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 828.59744_WP/pc%MolecularWeight;
    pc%liqDensityMol%B = 0.63993_WP/pc%MolecularWeight;
    pc%liqDensityMol%C = -0.00216_WP/pc%MolecularWeight;
    pc%liqDensityMol%D = 0.0_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ pc%MeltingPoint, pc%Tcrit /);
    pc%liqDensityMol%range(:,2) = (/ 4.98359080656_WP, 1.84463377882_WP /);
    pc%liqDensityMol%eqn = 100;

    !Retained same as ndodecane
    pc%vapPressure%A =  1.3747E02_WP;
    pc%vapPressure%B = -1.1976E04_WP;
    pc%vapPressure%C = -1.6698E01_WP;
    pc%vapPressure%D =  8.0906E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%vapPressure%range(:,2) = (/ 6.1534E-01_WP, 1.8223E06_WP /);
    pc%vapPressure%eqn = 101;  

    !Retained same as ndodecane
    pc%HeatOfVapMol%A = 1.6E05_WP/pc%MolecularWeight;
    pc%HeatOfVapMol%B = -187.629_WP/pc%MolecularWeight;
    pc%HeatOfVapMol%C = 0.365_WP/pc%MolecularWeight;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 23103486.6009_WP, 33097150.536_WP /);
    pc%HeatOfVapMol%eqn = 100;     

    !Retained same as ndodecane
    pc%solHeatCapacityMol%A = -4.0970E04_WP;
    pc%solHeatCapacityMol%B = 2.8390E03_WP;
    pc%solHeatCapacityMol%C = -1.3232E01_WP;
    pc%solHeatCapacityMol%D = 2.7972E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 29.97_WP, 258.41_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.2983E04_WP, 2.9175E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    !Retained same as ndodecane
    pc%liqHeatCapacityMol%A = 5.0821E05_WP;
    pc%liqHeatCapacityMol%B = -1.3687E03_WP;
    pc%liqHeatCapacityMol%C = 3.1015E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 263.57_WP, 330.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.6292E05_WP, 3.9429E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = -38118.9718_WP;
    pc%IG_HeatCapacityMol%B = 1033816.42_WP;
    pc%IG_HeatCapacityMol%C = -601695.64_WP;
    pc%IG_HeatCapacityMol%D = 135329.194_WP;
    pc%IG_HeatCapacityMol%E = 216.8348_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 423.15_WP, 1273.15_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 303067.824251_WP, 582196.957273_WP /);
    pc%IG_HeatCapacityMol%eqn = 115;

    !Retained same as ndodecane
    pc%SecondViralCoefMol%A = 8.800E-01_WP;
    pc%SecondViralCoefMol%B = -1.0910E03_WP;
    pc%SecondViralCoefMol%C = -5.0300E07_WP;
    pc%SecondViralCoefMol%D = -5.4871E21_WP;
    pc%SecondViralCoefMol%E = 1.4959E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 329.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.0701E01_WP, 1.3759E-01_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    !Retained same as ndodecane
    pc%liqViscosity%A = -2.0607E01_WP;
    pc%liqViscosity%B = 1.9430E03_WP;
    pc%liqViscosity%C = 1.3205E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 263.57_WP, 489.47_WP/);
    pc%liqViscosity%range(:,2) = (/ 2.8109E-03_WP, 2.1196E-04_WP /);
    pc%liqViscosity%eqn = 101;

    !Retained same as ndodecane
    pc%vapViscosity%A = 6.3440E-08_WP;
    pc%vapViscosity%B = 8.2870E-01_WP;
    pc%vapViscosity%C = 2.1950E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 263.57_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.5111E-06_WP, 1.5932E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 0.13924_WP;
    pc%liqThermalConductivity%B = 3.78253E-05_WP;
    pc%liqThermalConductivity%C = -2.89732E-07_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 263.57_WP, 489.47_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.2908218003E-01_WP, 8.83401018E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    !Retained same as ndodecane
    pc%vapThermalConductivity%A = 5.7190E-06_WP;
    pc%vapThermalConductivity%B = 1.4699E00_WP;
    pc%vapThermalConductivity%C = 5.7940E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 489.47_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.3537E-02_WP, 9.3010E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    !Retained same as ndodecane
    pc%SurfaceTension%A = 5.5493E-02_WP;
    pc%SurfaceTension%B = 1.3262E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.8150E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;

  end subroutine pc_diesel

  ! Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000
  subroutine pc_ndodecane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

    pc%name = 'ndodecane'
    pc%CAS = '112-40-3'
    ! components
    pc%ChemicalFormula%C = 12;         
    pc%ChemicalFormula%H = 26;
    pc%ChemicalFormula%O = 0;


    pc%MolecularWeight = 170.338_WP;   ! kg/kmol    
    pc%Tcrit = 658.00_WP;              ! K    
    pc%pcrit = 1.8200E6_WP;            ! Pa 
    pc%vcrit = 0.71600_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.238_WP;  ! -


    pc%MeltingPoint = 263.57_WP;             ! K
    pc%TripletPointT = 263.57_WP;            ! K
    pc%TripletPointp = 6.1520E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 489.47_WP;       ! K
    pc%liqMolarVolume  = 0.228605_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.9072E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 3.3180E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 6.2415E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 3.6840E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -7.5137E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.5764_WP;            ! -
    pc%RadiusOfGyration = 5.9140E-10_WP;      ! m
    pc%SolubilityParameter = 1.5900E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.12960_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.7740E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4151_WP;           ! -
    pc%FlashPoint = 347.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.6_WP, 4.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 476.00_WP;          ! K


    pc%solDensityMol%A = 5.2800E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 263.57_WP, MAX_REAL_WP /)
    pc%solDensityMol%range(:,2) = (/ 5.2800E00_WP, MAX_REAL_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.5541E-01_WP;
    pc%liqDensityMol%B = 2.5511E-01_WP;
    pc%liqDensityMol%C = 6.5800E02_WP;
    pc%liqDensityMol%D = 2.9368E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 4.5132E00_WP, 1.3932E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.3747E02_WP;
    pc%vapPressure%B = -1.1976E04_WP;
    pc%vapPressure%C = -1.6698E01_WP;
    pc%vapPressure%D =  8.0906E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%vapPressure%range(:,2) = (/ 6.1534E-01_WP, 1.8223E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.7337E07_WP;
    pc%HeatOfVapMol%B = 4.0681E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.2802E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -4.0970E04_WP;
    pc%solHeatCapacityMol%B = 2.8390E03_WP;
    pc%solHeatCapacityMol%C = -1.3232E01_WP;
    pc%solHeatCapacityMol%D = 2.7972E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 29.97_WP, 258.41_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.2983E04_WP, 2.9175E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 5.0821E05_WP;
    pc%liqHeatCapacityMol%B = -1.3687E03_WP;
    pc%liqHeatCapacityMol%C = 3.1015E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 263.57_WP, 330.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.6292E05_WP, 3.9429E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.1295E05_WP;
    pc%IG_HeatCapacityMol%B = 6.6330E05_WP;
    pc%IG_HeatCapacityMol%C = 1.7155E03_WP;
    pc%IG_HeatCapacityMol%D = 4.5161E05_WP;
    pc%IG_HeatCapacityMol%E = 7.7750E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.2442E05_WP, 7.4325E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 8.800E-01_WP;
    pc%SecondViralCoefMol%B = -1.0910E03_WP;
    pc%SecondViralCoefMol%C = -5.0300E07_WP;
    pc%SecondViralCoefMol%D = -5.4871E21_WP;
    pc%SecondViralCoefMol%E = 1.4959E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 329.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.0701E01_WP, 1.3759E-01_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0607E01_WP;
    pc%liqViscosity%B = 1.9430E03_WP;
    pc%liqViscosity%C = 1.3205E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 263.57_WP, 489.47_WP/);
    pc%liqViscosity%range(:,2) = (/ 2.8109E-03_WP, 2.1196E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 6.3440E-08_WP;
    pc%vapViscosity%B = 8.2870E-01_WP;
    pc%vapViscosity%C = 2.1950E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 263.57_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.5111E-06_WP, 1.5932E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0470E-01_WP;
    pc%liqThermalConductivity%B = -2.3260E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 263.57_WP, 489.47_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.4339E-01_WP, 9.0849E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.7190E-06_WP;
    pc%vapThermalConductivity%B = 1.4699E00_WP;
    pc%vapThermalConductivity%C = 5.7940E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 489.47_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.3537E-02_WP, 9.3010E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5493E-02_WP;
    pc%SurfaceTension%B = 1.3262E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.8150E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;

  end subroutine pc_ndodecane


  subroutine pc_1butanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

    ! all alternative names
    pc%name = '1butanol'; 
    pc%CAS = '71-36-3';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 74.123_WP;   ! kg/kmol    
    pc%Tcrit = 563.05_WP;             ! K    
    pc%pcrit = 4.4230E06_WP;          ! Pa 
    pc%vcrit = 0.27400_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.260_WP; ! -


    pc%MeltingPoint = 183.85_WP;              ! K
    pc%TripletPointT = 184.51_WP;             ! K
    pc%TripletPointp = 5.7220E-04_WP;         ! Pa
    pc%NormalBoilingPoint = 390.81_WP;        ! K
    pc%liqMolarVolume = 0.091991_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.7460E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.5030E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.6148E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  9.3720E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.4560E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.5935_WP;            ! -
    pc%RadiusOfGyration = 3.2250E-10_WP;      ! m
    pc%SolubilityParameter = 2.3350E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05240_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.6200E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3971_WP;           ! -
    pc%FlashPoint = 302.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.4_WP, 11.2_WP/);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 616.00_WP;          ! K


    pc%solDensityMol%A = 1.4100E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 184.51_WP, MAX_REAL_WP /);
    pc%solDensityMol%range(:,2) = (/ 1.4100E01_WP, MAX_REAL_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.6500E-01_WP;
    pc%liqDensityMol%B = 2.6660E-01_WP;
    pc%liqDensityMol%C = 5.6305E02_WP;
    pc%liqDensityMol%D = 2.4419E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 184.51_WP, 563.05_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.2016E01_WP, 3.6197E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  9.3173E01_WP;
    pc%vapPressure%B = -9.1859E03_WP;
    pc%vapPressure%C = -9.7464E00_WP;
    pc%vapPressure%D =  4.7796E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 184.51_WP, 563.05_WP /);
    pc%vapPressure%range(:,2) = (/ 5.7220E-04_WP, 4.3392E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.7390E07_WP;
    pc%HeatOfVapMol%B = 1.7300E-01_WP;
    pc%HeatOfVapMol%C = 2.9150E-01_WP;
    pc%HeatOfVapMol%D = 0;
    pc%HeatOfVapMol%E = 0;
    pc%HeatOfVapMol%range(:,1) = (/ 184.51_WP, 563.05_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.0575E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.1680E04_WP;
    pc%solHeatCapacityMol%B = 1.0262E03_WP;
    pc%solHeatCapacityMol%C = -1.7010E00_WP;
    pc%solHeatCapacityMol%D = -1.5820E-02_WP;
    pc%solHeatCapacityMol%E = 6.6500E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 184.51_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 8.0477E0_WP, 9.7456E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.9120E05_WP;
    pc%liqHeatCapacityMol%B = -7.3040E02_WP;
    pc%liqHeatCapacityMol%C = 2.2998E00_WP;
    pc%liqHeatCapacityMol%D = 0;
    pc%liqHeatCapacityMol%E = 0;
    pc%liqHeatCapacityMol%range(:,1) = (/ 184.51_WP, 390.81_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.3473E05_WP, 2.5701E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 7.4540E04_WP;
    pc%IG_HeatCapacityMol%B = 2.5907E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6073E03_WP;
    pc%IG_HeatCapacityMol%D = 1.7320E05_WP;
    pc%IG_HeatCapacityMol%E = 7.1240E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 8.1616E04_WP, 2.8509E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.8900E-01_WP;
    pc%SecondViralCoefMol%B = -1.8200E02_WP;
    pc%SecondViralCoefMol%C = -4.0500E07_WP;
    pc%SecondViralCoefMol%D = -2.2700E20_WP;
    pc%SecondViralCoefMol%E = 4.5400E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 282.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ 3.9128E00_WP, 5.5659E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -3.5426E01_WP;
    pc%liqViscosity%B = 3.1845E03_WP;
    pc%liqViscosity%C = 3.2965E00_WP;
    pc%liqViscosity%D = -3.0000E-27_WP;
    pc%liqViscosity%E = 1.0000E01_WP;
    pc%liqViscosity%range(:,1) = (/ 190.00_WP, 390.81_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.5440E-01_WP, 3.8873E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.4031E-06_WP;
    pc%vapViscosity%B = 4.6110E-01_wP;
    pc%vapViscosity%C = 5.3700E+02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 184.51_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.9786E-06_WP, 2.2066E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1360E-01_WP;
    pc%liqThermalConductivity%B = -2.0340E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 184.51_WP, 390.81_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7607E-01_WP, 1.3411E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -4.4940E-02_WP;
    pc%vapThermalConductivity%B = 4.4600E-02_WP;
    pc%vapThermalConductivity%C = -1.3552E+03_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 370.00_WP, 800.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1971E-02_WP, 8.7247E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 4.9830E-02_WP;
    pc%SurfaceTension%B = -8.5400E-05_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 273.15_WP, 413.15_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.6503E-02_WP, 1.4547E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_1butanol


  subroutine pc_1decanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1decanol'; 
    pc%CAS = '112-30-1';
    ! components
    pc%ChemicalFormula%C = 10;         
    pc%ChemicalFormula%H = 22;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 158.284_WP;  ! kg/kmol    
    pc%Tcrit = 687.30_WP;             ! K    
    pc%pcrit = 2.3150E06_WP;          ! Pa 
    pc%vcrit = 0.64900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.263_WP; ! -


    pc%MeltingPoint = 280.05_WP;              ! K
    pc%TripletPointT = 280.05_WP;             ! K
    pc%TripletPointp = 1.2019E-01_WP;         ! Pa
    pc%NormalBoilingPoint = 504.07_WP;        ! K
    pc%liqMolarVolume = 0.191854_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.9670E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -9.8790E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  5.9700E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = 3.7656E07_WP;           ! J/kmol
    pc%StandardNetHeatOfCombustion = -6.1159E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.2650E08_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -4.7810E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 4.1690E05_WP;         ! J/kmolK


    pc%AcentricFactor = 0.6219_WP;            ! -
    pc%RadiusOfGyration = 5.4990E-10_WP;      ! m
    pc%SolubilityParameter = 2.02500E04_WP;   ! (J/m^3)^.5 
    pc%DipoleMoment = 5.4000E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.11378_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.5730E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4350_WP;           ! -
    pc%FlashPoint = 355.37_WP;                ! K
    pc%FlammabilityLimits = (/ 0.7_WP, 5.5_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 376.0_WP, 416.0_WP /); ! K
    pc%AutoignitionTemp = 561.15_WP;          ! K


    pc%solDensityMol%A = 6.1900E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 280.05_WP, MAX_REAL_WP /);
    pc%solDensityMol%range(:,2) = (/ 6.1900E00_WP, MAX_REAL_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.7384E-01_WP;
    pc%liqDensityMol%B = 2.4241E-01_WP;
    pc%liqDensityMol%C = 6.8730E02_WP;
    pc%liqDensityMol%D = 2.6646E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 280.05_WP, 687.30_WP /);
    pc%liqDensityMol%range(:,2) = (/ 5.2902E00_WP, 1.5422E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  2.5059E02_WP;
    pc%vapPressure%B = -1.9169E04_WP;
    pc%vapPressure%C = -3.2903E01_WP;
    pc%vapPressure%D =  1.4627E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 280.05_WP, 687.30_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2019E-01_WP, 2.3302E06_WP /);
    pc%vapPressure%eqn = 101;

    pc%HeatOfVapMol%A = 1.1750E08_WP;
    pc%HeatOfVapMol%B = 6.5112E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 280.05_WP, 687.30_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 8.3570E07_WP, 0.0_WP /); 
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 4.11026E05_WP;
    pc%solHeatCapacityMol%B = -2.4283E03_WP;
    pc%solHeatCapacityMol%C = 7.3148E00_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 180.00_WP, 270.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.8017E05_WP, 3.5787E05_WP /);   
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 4.9885E06_WP;
    pc%liqHeatCapacityMol%B = -5.2898E04_WP;
    pc%liqHeatCapacityMol%C = 2.1635E02_WP;
    pc%liqHeatCapacityMol%D = -3.7538E-01_WP;
    pc%liqHeatCapacityMol%E = 2.3674E-04_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 280.05_WP, 504.07_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.5373E05_WP, 5.0206E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.6984E05_WP;
    pc%IG_HeatCapacityMol%B = 5.5920E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5680E03_WP;
    pc%IG_HeatCapacityMol%D = 3.9380E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2050E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.8502E05_WP, 6.2186E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.8836E-01_WP;
    pc%SecondViralCoefMol%B = -5.8491E02_WP;
    pc%SecondViralCoefMol%C = -1.1672E08_WP;
    pc%SecondViralCoefMol%D = 3.8936E21_WP;
    pc%SecondViralCoefMol%E = -1.5303E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 343.65_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.9662E00_WP, 6.3948E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -8.0656E01_WP;
    pc%liqViscosity%B = 6.3255E03_WP;
    pc%liqViscosity%C = 9.6460E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 285.00_WP, 504.07_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.9500E-02_WP, 3.0879E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 5.5065E-08_WP;
    pc%vapViscosity%B = 8.3410E-01_WP;
    pc%vapViscosity%C = 7.9560E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 280.05_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.7154E-06_WP, 1.6215E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2800E-01_WP;
    pc%liqThermalConductivity%B = -2.2300E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 280.05_WP, 503.35_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6555E-01_WP, 1.1575E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -3.0720E-01_WP;
    pc%vapThermalConductivity%B = 4.8900E-01_WP;
    pc%vapThermalConductivity%C = -6.7500E04_WP;
    pc%vapThermalConductivity%D = -2.9400E07_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 504.07_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.5906E-02_WP, 9.3886E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1263E-02_WP;
    pc%SurfaceTension%B = 1.0395E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 280.05_WP, 687.30_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.9754E-02_WP, 0.0000E00_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1decanol


  subroutine pc_1dodecanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1dodecanol'; 
    pc%CAS = '112-53-8';
    ! components
    pc%ChemicalFormula%C = 12;         
    pc%ChemicalFormula%H = 26;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 186.338_WP;  ! kg/kmol    
    pc%Tcrit = 719.4_WP;              ! K    
    pc%pcrit = 1.9940E06_WP;          ! Pa 
    pc%vcrit = 0.716_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.239_WP; ! -


    pc%MeltingPoint = 296.95_WP;              ! K
    pc%TripletPointT = 296.95_WP;             ! K
    pc%TripletPointp = 0.101105_WP;           ! Pa
    pc%NormalBoilingPoint = 536.95_WP;        ! K
    pc%liqMolarVolume = 0.224322_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -4.3670E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -8.0820E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  6.7521E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = 3.1380E07_WP;           ! J/kmol
    pc%StandardNetHeatOfCombustion = -7.3380E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.1500E08_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -5.2850E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 4.8300E05_WP;         ! J/kmolK    


    pc%AcentricFactor = 0.666353_WP;          ! -
    pc%RadiusOfGyration = 6.1190E-10_WP;      ! m
    pc%SolubilityParameter = 1.9750E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.6400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.13424_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.8430E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4413_WP;           ! -
    pc%FlashPoint = 400.0_WP;                 ! K
    pc%FlammabilityLimits = (/ 0.6_WP, 5.1_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 397_WP, 441_WP /); ! K
    pc%AutoignitionTemp = 548.15_WP;          ! K


    pc%solDensityMol%A = 5.5881E00_WP;
    pc%solDensityMol%B = -1.6655E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 78.15_WP, 296.95_WP /);
    pc%solDensityMol%range(:,2) = (/ 5.4579E00_WP, 5.0935E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.5818E-01_WP;
    pc%liqDensityMol%B = 2.5626E-01_WP;
    pc%liqDensityMol%C = 7.1940E02_WP;
    pc%liqDensityMol%D = 2.9963E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 296.95_WP, 719.40_WP /);
    pc%liqDensityMol%range(:,2) = (/ 4.4623E00_WP, 1.3977E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  3.8678E02_WP;
    pc%vapPressure%B = -2.4382E04_WP;
    pc%vapPressure%C = -5.6170E01_WP;
    pc%vapPressure%D =  4.3253E-02_WP;
    pc%vapPressure%E = 1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 296.95_WP, 719.40_WP /);
    pc%vapPressure%range(:,2) = (/ 1.011E-01_WP ,1.9716E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 1.2925E08_WP;
    pc%HeatOfVapMol%B = 6.7609E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 296.95_WP, 719.40_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 9.0182E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 3.5700E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 293.15_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.5700E05_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 5.1422E06_WP;
    pc%liqHeatCapacityMol%B = -5.4970E04_WP;
    pc%liqHeatCapacityMol%C = 2.2872E02_WP;
    pc%liqHeatCapacityMol%D = -4.0331E-01_WP;
    pc%liqHeatCapacityMol%E = 2.5844E-04_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 296.95_WP, 486.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 4.3614E05_WP, 5.7106E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.0150E05_WP;
    pc%IG_HeatCapacityMol%B = 6.3270E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5625E03_WP;
    pc%IG_HeatCapacityMol%D = 4.6733E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2120E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.1943E05_WP, 7.3406E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 6.2391E-01_WP;
    pc%SecondViralCoefMol%B = -7.8675E02_WP;
    pc%SecondViralCoefMol%C = -1.6374E08_WP;
    pc%SecondViralCoefMol%D = 7.1921E21_WP;
    pc%SecondViralCoefMol%E = -2.8712E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 359.70_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.9011E00_WP, 5.1100E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -5.3664E01_WP;
    pc%liqViscosity%B = 5.2265E03_WP;
    pc%liqViscosity%C = 5.6140E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 293.15_WP, 536.95_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.9373E-02_WP, 1.7667E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 5.5812E-08_WP;
    pc%vapViscosity%B = 8.2000E-01_WP;
    pc%vapViscosity%C = 5.4330E01_WP;
    pc%vapViscosity%D = 1.4110E04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 296.95_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.4285E-06_WP, 1.5065E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0444E-01_WP;
    pc%liqThermalConductivity%B = -1.8378E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 296.95_WP, 536.95_WP /);
    pc%liqThermalConductivity%range(:,2) = (/  1.4987E-01_WP, 1.0576E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -2.8300E-03_WP;
    pc%vapThermalConductivity%B = 9.1261E-01_WP;
    pc%vapThermalConductivity%C = -1.9181E04_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 536.95_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.5266E-02_WP, 8.5114E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.3034E-02_WP;
    pc%SurfaceTension%B = 1.1008E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 293.15_WP, 719.40_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.9808E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1dodecanol


  subroutine pc_1heptanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1heptanol'; 
    pc%CAS = '111-70-6';
    ! components
    pc%ChemicalFormula%C = 7;         
    pc%ChemicalFormula%H = 16;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 116.203_WP;  ! kg/kmol    
    pc%Tcrit = 632.60_WP;             ! K    
    pc%pcrit = 3.0580E06_WP;          ! Pa 
    pc%vcrit = 0.43500_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.253_WP; ! -


    pc%MeltingPoint = 239.15_WP;              ! K
    pc%TripletPointT = 239.15_WP;             ! K
    pc%TripletPointp = 1.1773E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 449.45_WP;        ! K
    pc%liqMolarVolume = 0.141801_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.3640E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.2530E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.7919E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.8180E07_WP;           ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.2860E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.4650E08_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -4.0330E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 3.2590E05_WP;         ! J/kmolK


    pc%AcentricFactor = 0.5670_WP;            ! -
    pc%RadiusOfGyration = 4.3800E-10_WP;      ! m
    pc%SolubilityParameter = 2.1540E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.8000E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.08309_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.1680E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4223_WP;           ! -
    pc%FlashPoint = 350.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.0_WP, 7.2_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 344.0_WP, 379.0_WP /); ! K
    pc%AutoignitionTemp = 555.00_WP;          ! K


    pc%solDensityMol%A = 8.6800E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 239.15_WP, 8.6800E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.0531E-01_WP;
    pc%liqDensityMol%B = 2.6329E-01_WP;
    pc%liqDensityMol%C = 6.3260E02_WP;
    pc%liqDensityMol%D = 2.7375E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 239.15_WP, 632.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 7.4209E00_WP, 2.2990E+00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  2.4331E02_WP;
    pc%vapPressure%B = -1.7320E04_WP;
    pc%vapPressure%C = -3.2183E01_WP;
    pc%vapPressure%D =  1.6443E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 239.15_WP, 632.60_WP /);
    pc%vapPressure%range(:,2) = (/ 1.1773E-02_WP, 3.0663E+06_WP /);
    pc%vapPressure%eqn = 101;

    pc%HeatOfVapMol%A = 1.0075E08_WP;
    pc%HeatOfVapMol%B = 6.1052E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 239.15_WP, 632.60_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 7.5393E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 1.9300E04_WP;
    pc%solHeatCapacityMol%B = 9.8150E02_WP;
    pc%solHeatCapacityMol%C = -3.0750E00_WP;
    pc%solHeatCapacityMol%D = 6.9530E-03_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 80.00_WP, 239.15_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 8.1700E04_WP, 1.7326E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;


    pc%liqHeatCapacityMol%A = 2.4168E06_WP;
    pc%liqHeatCapacityMol%B = -2.6105E04_WP;
    pc%liqHeatCapacityMol%C = 1.1003E02_WP;
    pc%liqHeatCapacityMol%D = -1.9172E-01_WP;
    pc%liqHeatCapacityMol%E = 1.1968E-04_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 239.15_WP, 449.45_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.3590E05_WP, 3.8773E+05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.2215E05_WP;
    pc%IG_HeatCapacityMol%B = 3.9910E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5800E03_WP;
    pc%IG_HeatCapacityMol%D = 2.8350E05_WP;
    pc%IG_HeatCapacityMol%E = 7.1770E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3330E05_WP, 4.5346E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.8979E-01_WP;
    pc%SecondViralCoefMol%B = -3.0797E02_WP;
    pc%SecondViralCoefMol%C = -9.8592E07_WP;
    pc%SecondViralCoefMol%D = -2.7797E20_WP;
    pc%SecondViralCoefMol%E = 1.3614E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 316.30_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.1382E00_WP, 5.5313E-02_WP /); 
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -7.5245E01_WP;
    pc%liqViscosity%B = 5.7567E03_WP;
    pc%liqViscosity%C = 8.9146E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 239.15_WP, 449.45_WP /);
    pc%liqViscosity%range(:,2) = (/ 9.5615E-02_WP, 3.4016E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.5720E-07_WP;
    pc%vapViscosity%B = 6.5020E-01_WP;
    pc%vapViscosity%C = 2.4860E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 239.15_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.4397E-06_WP, 1.8384E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2390E-01_WP;
    pc%liqThermalConductivity%B = -2.2600E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 239.15_WP, 573.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6985E-01_WP, 9.4368E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -6.1993E-02_WP;
    pc%vapThermalConductivity%B = 2.7920E-01_WP;
    pc%vapThermalConductivity%C = -3.3360E03_WP;
    pc%vapThermalConductivity%D = -1.6420E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 504.07_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.3447E-02_WP, 1.0722E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.3765E-02_WP;
    pc%SurfaceTension%B = 1.0974E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 239.15_WP, 632.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1928E-02_WP, 0.0000E00_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1heptanol


  subroutine pc_1hexanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1hexanol'; 
    pc%CAS = '111-27-3';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 102.77_WP;   ! kg/kmol    
    pc%Tcrit = 610.30_WP;             ! K    
    pc%pcrit = 3.4170E06_WP;          ! Pa 
    pc%vcrit = 0.38700_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.261_WP; ! -


    pc%MeltingPoint = 228.55_WP;              ! K
    pc%TripletPointT = 228.55_WP;             ! K
    pc%TripletPointp = 2.9550E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 430.55_WP;        ! K
    pc%liqMolarVolume = 0.125200_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.1650E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -1.3440E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 4.4010E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.5400E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.6766E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = - 1.4980E+08_WP;    ! J/kmol
    pc%StandardHeatOfFormation = -3.7750E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.8740E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.5764_WP;            ! -
    pc%RadiusOfGyration = 4.1440E-10_WP;      ! m
    pc%SolubilityParameter = 2.1830E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5000E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.07286_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.0330E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4161_WP;           ! -
    pc%FlashPoint = 356.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 8.2_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 332.0_WP, 366.0_WP /); ! K
    pc%AutoignitionTemp = 558.00_WP;          ! K


    pc%solDensityMol%A = 9.8900E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 228.55_WP, 9.8900E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.8425E-01_WP;
    pc%liqDensityMol%B = 2.6477E-01_WP;
    pc%liqDensityMol%C = 6.1030E02_WP;
    pc%liqDensityMol%D = 2.4395E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 228.55_WP, 610.30_WP /);
    pc%liqDensityMol%range(:,2) = (/ 8.4540E00_WP, 2.5843E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.2508E02_WP;
    pc%vapPressure%B = -1.1692E04_WP;
    pc%vapPressure%C = -1.4258E01_WP;
    pc%vapPressure%D =  1.1102E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 228.55_WP, 610.30_WP /);
    pc%vapPressure%range(:,2) = (/ 2.9550E-02_WP, 3.4221E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 8.9962E07_WP;
    pc%HeatOfVapMol%B = 5.5179E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 228.55_WP, 610.30_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.9442E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.8080E04_WP;
    pc%solHeatCapacityMol%B = 1.4083E03_WP;
    pc%solHeatCapacityMol%C = -4.2360E00_WP;
    pc%solHeatCapacityMol%D = -1.1200E-03_WP;
    pc%solHeatCapacityMol%E = 3.0800E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 200.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.0892E+03_WP, 1.3446E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.6386E06_WP;
    pc%liqHeatCapacityMol%B = -1.7261E04_WP;
    pc%liqHeatCapacityMol%C = 7.1721E01_WP;
    pc%liqHeatCapacityMol%D = -1.2026E-01_WP;
    pc%liqHeatCapacityMol%E = 7.1087E-05_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 228.55_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.9821E05_WP, 3.5197E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0625E05_WP;
    pc%IG_HeatCapacityMol%B = 3.5210E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5835E03_WP;
    pc%IG_HeatCapacityMol%D = 2.4620E05_WP;
    pc%IG_HeatCapacityMol%E = 7.1575E+02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.1607E05_WP, 3.9726E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.5305E-01_WP;
    pc%SecondViralCoefMol%B = -2.6298E02_WP;
    pc%SecondViralCoefMol%C = -7.5678E07_WP;
    pc%SecondViralCoefMol%D = -2.2383E20_WP;
    pc%SecondViralCoefMol%E = 2.0576E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 305.15_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.352E00_WP, 5.5299E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -4.3289E01_WP;
    pc%liqViscosity%B = 4.0412E03_WP;
    pc%liqViscosity%C = 4.2714E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 250.00_WP, 430.55_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.9018E-02_WP, 3.3662E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.5773E-07_WP;
    pc%vapViscosity%B = 7.1890E-01_WP;
    pc%vapViscosity%C = 1.6330E+02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 228.55_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.5672E-06_WP, 1.9450E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1930E-01_WP;
    pc%liqThermalConductivity%B = -2.2000E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 228.55_WP, 575.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6902E-01_WP, 9.2800E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -4.9355E06_WP;
    pc%vapThermalConductivity%B = -1.6530E-01_WP;
    pc%vapThermalConductivity%C = 1.5631E09_WP;
    pc%vapThermalConductivity%D = -1.5752E13_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 430.55_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.2264E-02_WP, 1.1104E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0895E-02_WP;
    pc%SurfaceTension%B = 1.0073E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 228.55_WP, 610.30_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1727E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1hexanol


  subroutine pc_1nonanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1nonanol'; 
    pc%CAS = '143-08-8';
    ! components
    pc%ChemicalFormula%C = 9;         
    pc%ChemicalFormula%H = 20;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 144.257_WP;   ! kg/kmol    
    pc%Tcrit = 670.70_WP;             ! K    
    pc%pcrit = 2.5280E06_WP;           ! Pa 
    pc%vcrit = 0.57200_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.259_WP; ! -


    pc%MeltingPoint = 268.15_WP;              ! K
    pc%TripletPointT = 268.15_WP;             ! K
    pc%TripletPointp = 7.1702E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 486.25_WP;        ! K
    pc%liqMolarVolume = 0.174942_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.8120E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.1220E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  5.5770E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.8800E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -5.5007E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.3300E08_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -4.5360E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 3.8460E05_WP;         ! J/kmolK    


    pc%AcentricFactor = 0.5997_WP;            ! -
    pc%RadiusOfGyration = 5.1080E-10_WP;      ! m
    pc%SolubilityParameter = 2.0720E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.3700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.10355_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.4380E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4319_WP;           ! -
    pc%FlashPoint = 347.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.8_WP, 6.1_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 366.0_WP, 405.0_WP /); ! K
    pc%AutoignitionTemp = 550.00_WP;          ! K


    pc%solDensityMol%A = 6.8500E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 268.15_WP, 6.8500E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 4.4231E-01_WP;
    pc%liqDensityMol%B = 2.5300E-01_WP;
    pc%liqDensityMol%D = 2.5261E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 268.15_WP, 670.70_WP /);
    pc%liqDensityMol%range(:,2) = (/ 5.8515E00_WP, 1.7483E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  2.5234E02_WP;
    pc%vapPressure%B = -1.8775E04_WP;
    pc%vapPressure%C = -3.3278E01_WP;
    pc%vapPressure%D =  1.5535E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 268.15_WP, 670.70_WP /);
    pc%vapPressure%range(:,2) = (/ 7.1702E-02_WP, 2.5489E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 1.1484E08_WP;
    pc%HeatOfVapMol%B = 6.6695E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 268.15_WP, 670.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 8.1701E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 3.8373E05_WP;
    pc%solHeatCapacityMol%B = -1.4621E03_WP;
    pc%solHeatCapacityMol%C = 4.4071E00_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 200.00_WP, 260.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.6759E05_WP, 3.0150E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.0483E07_WP;
    pc%liqHeatCapacityMol%B = -1.1522E05_WP;
    pc%liqHeatCapacityMol%C = 4.7687E02_WP;
    pc%liqHeatCapacityMol%D = -8.5381E-01_WP;
    pc%liqHeatCapacityMol%E = 5.6246E-04_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 310.00_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.5059E05_WP, 4.6494E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.5400E05_WP;
    pc%IG_HeatCapacityMol%B = 4.9360E05_WP;
    pc%IG_HeatCapacityMol%C = 1.580E03_WP;
    pc%IG_HeatCapacityMol%D = 3.5880E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2111E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.6777E05_WP, 5.6606E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.3008E-01_WP;
    pc%SecondViralCoefMol%B = -5.1453E02_WP;
    pc%SecondViralCoefMol%C = -9.3063E07_WP;
    pc%SecondViralCoefMol%D = 2.9028E+21_WP;
    pc%SecondViralCoefMol%E = - 1.1382E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 335.35_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.6432E00_WP, 5.9569E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -7.1348E00_WP;
    pc%liqViscosity%B = 2.7763E03_WP;
    pc%liqViscosity%C = -1.2064E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 280.00_WP, 486.25_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.8002E-02_WP, 1.3792E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.2000E-07_WP;
    pc%vapViscosity%B = 7.4000E-01_WP;
    pc%vapViscosity%C = 1.8000E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 268.15_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.4992E-06_WP, 1.6877E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2920E-01_WP;
    pc%liqThermalConductivity%B = -2.3000E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 268.15_WP, 578.65_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6753E-01_WP, 9.6111E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -3.0715E01_WP;
    pc%vapThermalConductivity%B = -1.0750E-01_WP;
    pc%vapThermalConductivity%C = 8.1070E03_WP;
    pc%vapThermalConductivity%D = -1.5683E08_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 486.25_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.4464E-02_WP, 9.8947E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1752E-02_WP;
    pc%SurfaceTension%B = 1.0562E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 268.15_WP, 670.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.0183E-02_WP, 0.0000E00_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1nonanol


  subroutine pc_1octanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1octanol'; 
    pc%CAS = '111-87-5';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 130.230_WP;  ! kg/kmol    
    pc%Tcrit = 652.50_WP;             ! K    
    pc%pcrit = 2.7770E06_WP;          ! Pa 
    pc%vcrit = 0.49700_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.254_WP; ! -


    pc%MeltingPoint = 257.65_WP;              ! K
    pc%TripletPointT = 257.65_WP;             ! K
    pc%TripletPointp = 8.5564E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 468.35_WP;        ! K
    pc%liqMolarVolume = 0.158231_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.5550E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.1550E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  5.1850E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.2500E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.8983E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.3690E08_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -4.2650E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 3.5200E05_WP;         ! J/kmolK


    pc%AcentricFactor = 0.5829_WP;            ! -
    pc%RadiusOfGyration = 4.7870E-10_WP;      ! m
    pc%SolubilityParameter = 2.1020E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5000E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.09332_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.3030E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4276_WP;           ! -
    pc%FlashPoint = 354.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.9_WP, 6.4_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 355.0_WP, 392.0_WP /); ! K
    pc%AutoignitionTemp = 555.00_WP;          ! K


    pc%solDensityMol%A = 7.6600E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 257.65_WP, 7.6600E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.1836E-01_WP;
    pc%liqDensityMol%B = 2.5771E-01_WP;
    pc%liqDensityMol%C = 6.5250E02_WP;
    pc%liqDensityMol%D = 2.7713E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 257.65_WP, 652.50_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.5440E00_WP, 2.0114E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.5390E02_WP;
    pc%vapPressure%B = -1.4181E+04_WP;
    pc%vapPressure%C = -1.8251E+01_WP;
    pc%vapPressure%D =  1.2261E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 257.65_WP, 652.50_WP /);
    pc%vapPressure%range(:,2) = (/ 8.5564E-02_WP, 2.7625E+06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 1.0730E+08_WP;
    pc%HeatOfVapMol%B = 6.4524E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 257.65_WP, 652.50_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 7.7596E+07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 4.1740E+04_WP;
    pc%solHeatCapacityMol%B = 5.7900E+02_WP;
    pc%solHeatCapacityMol%C = 2.7800E-01_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 101.90_WP, 216.90_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0363E+05_WP, 1.8040E+05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 5.7137E+05_WP;
    pc%liqHeatCapacityMol%B = -4.8490E+03_WP;
    pc%liqHeatCapacityMol%C = 1.9725E+01_WP;
    pc%liqHeatCapacityMol%D = -2.1532E-02_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 250.00_WP, 468.35_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.5550E+05_WP, 4.1500E+05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.3805E+05_WP;
    pc%IG_HeatCapacityMol%B = 4.4590E+05_WP;
    pc%IG_HeatCapacityMol%C = 1.5751E+03_WP;
    pc%IG_HeatCapacityMol%D = 3.2016E+05_WP;
    pc%IG_HeatCapacityMol%E = 7.1880E+02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.5055E+05_WP, 5.0965E+05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.2384E-01_WP;
    pc%SecondViralCoefMol%B = -3.5035E+02_WP;
    pc%SecondViralCoefMol%C = -1.2731E+08_WP;
    pc%SecondViralCoefMol%D = -3.0518E+20_WP;
    pc%SecondViralCoefMol%E = -1.1146E+22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 326.25_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.0600E+00_WP, 5.2540E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -5.0674E+01_WP;
    pc%liqViscosity%B = 4.7251E+03_WP;
    pc%liqViscosity%C = 5.2499E+00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 280.00_WP, 425.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.4749E-02_WP, 4.1665E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.7520E-07_WP;
    pc%vapViscosity%B = 6.9410E-01_WP;
    pc%vapViscosity%C = 2.0680E+02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 257.65_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.5827E-06_WP, 1.7547E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.3160E-01_WP;
    pc%liqThermalConductivity%B = -2.4070E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 257.65_WP, 570.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6958E-01_WP, 9.4365E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -3.0238E-03_WP;
    pc%vapThermalConductivity%B = 8.7450E-01_WP;
    pc%vapThermalConductivity%C = -1.3352E+04_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 468.35_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.3796E-02_WP, 1.0288E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1595E-02_WP;
    pc%SurfaceTension%B = 1.0483E+00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 257.65_WP, 652.50_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.0355E-02_WP, 0.0000E+00_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1octanol


  subroutine pc_1pentanal(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1pentanal'; 
    pc%CAS = '110-62-3';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 86.134_WP;   ! kg/kmol    
    pc%Tcrit = 566.10_WP;             ! K    
    pc%pcrit = 3.9700E06_WP;          ! Pa 
    pc%vcrit = 0.31300_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.264_WP; ! -


    pc%MeltingPoint = 182.00_WP;              ! K
    pc%TripletPointT = 182.00_WP;             ! K
    pc%TripletPointp = 5.2282E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 376.15_WP;        ! K
    pc%liqMolarVolume = 0.107035_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.2780E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.0710E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.8289E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.5000E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.9100E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.3472_WP;            ! -
    pc%RadiusOfGyration = 3.5160E-10_WP;      ! m
    pc%SolubilityParameter = 1.8440E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 8.5700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05950_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.5400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3917_WP;           ! -
    pc%FlashPoint = 285.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 9.5_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 480.00_WP;          ! K


    pc%solDensityMol%A = 1.2300E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 182.00_WP, 1.2300E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.3871E-01_WP;
    pc%liqDensityMol%B = 2.6252E-01_WP;
    pc%liqDensityMol%C = 5.6610E02_WP;
    pc%liqDensityMol%D = 2.9444E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 182.00_WP, 566.10_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0534E01_WP, 3.1948E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.4958E02_WP;
    pc%vapPressure%B = -8.8900E03_WP;
    pc%vapPressure%C = -2.0697E01_WP;
    pc%vapPressure%D =  2.2101E-02_WP;
    pc%vapPressure%E =  1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 182.00_WP, 566.10_WP /);
    pc%vapPressure%range(:,2) = (/ 5.2282E-02_WP, 3.9685E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.1478E07_WP;
    pc%HeatOfVapMol%B = 3.7541E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 182.00_WP, 566.10_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.4502E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 9.3786E03_WP;
    pc%solHeatCapacityMol%B = 7.4571E02_WP;
    pc%solHeatCapacityMol%C = -8.5714E-01_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 80.00_WP, 150.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 6.3550E04_WP, 1.0195E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.1205E05_WP;
    pc%liqHeatCapacityMol%B = 2.5778E02_WP;
    pc%liqHeatCapacityMol%C = 0.0_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 200.00_WP, 376.15_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.6361E05_WP, 2.0901E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0743E05_WP;
    pc%IG_HeatCapacityMol%B = 2.8363E05_WP;
    pc%IG_HeatCapacityMol%C = 1.9549E03_WP;
    pc%IG_HeatCapacityMol%D = 2.0146E05_WP;
    pc%IG_HeatCapacityMol%E = 8.9044E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.0960E05_WP, 3.2404E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.9740E-01_WP;
    pc%SecondViralCoefMol%B = -2.2484E02_WP;
    pc%SecondViralCoefMol%C = -4.0107E07_WP;
    pc%SecondViralCoefMol%D = -1.9523E20_WP;
    pc%SecondViralCoefMol%E = 3.0063E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 283.05_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.5262E00_WP, 3.5616E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.0846E01_WP;
    pc%liqViscosity%B = 9.8001E02_WP;
    pc%liqViscosity%C = -5.4565E-03_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 182.00_WP, 376.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 4.1290E-03_WP, 2.5532E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.2700E-07_WP;
    pc%vapViscosity%B = 6.7670E-01_WP;
    pc%vapViscosity%C = 1.9174E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 182.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.7404E-06_WP, 2.0415E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2697E-01_WP;
    pc%liqThermalConductivity%B = -3.3227E-04_WP;
    pc%liqThermalConductivity%C = 1.1770E-07_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 182.00_WP, 513.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7040E-01_WP, 8.7459E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -4.9187E06_WP;
    pc%vapThermalConductivity%B = -1.0297E-01_WP;
    pc%vapThermalConductivity%C = 2.6911E09_WP;
    pc%vapThermalConductivity%D = -2.3179E13_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 376.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.7048E-02_WP, 1.1788E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.4779E-02_WP;
    pc%SurfaceTension%B = 1.2556E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 182.00_WP, 566.10_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.9804E-02_WP, 0.0000E00_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1pentanal


  subroutine pc_1pentanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1pentanol'; 
    pc%CAS = '71-41-0';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 88.150_WP;   ! kg/kmol    
    pc%Tcrit = 588.10_WP;             ! K    
    pc%pcrit = 3.8970E06_WP;          ! Pa 
    pc%vcrit = 0.32600_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.260_WP; ! -


    pc%MeltingPoint = 195.56_WP;              ! K
    pc%TripletPointT = 195.56_WP;             ! K
    pc%TripletPointp = 3.4525E-4_WP;          ! Pa
    pc%NormalBoilingPoint = 410.95_WP;        ! K
    pc%liqMolarVolume = 0.108538_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.9874E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -1.4602E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy =  4.0250E05_WP;            ! J/kmolK
    pc%HeatFusionAtMeltPoint =  9.7906E06_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.0605E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = -1.5610E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -3.5160E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.5890E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.5731_WP;            ! -
    pc%RadiusOfGyration = 3.6790E-10_WP;      ! m
    pc%SolubilityParameter = 2.2520E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.6706E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06263_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.9800E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4080_WP;           ! -
    pc%FlashPoint = 305.93_WP;                ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 10.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 316.0_WP, 353.0_WP /); ! K
    pc%AutoignitionTemp = 573.15_WP;          ! K


    pc%solDensityMol%A = 1.1800E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 195.56_WP, 1.1800E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.1754E-01_WP;
    pc%liqDensityMol%B = 2.6732E-01_WP;
    pc%liqDensityMol%C = 5.8810E02_WP;
    pc%liqDensityMol%D = 2.5348E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 195.56_WP, 588.10_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0061E01_WP, 3.0583E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.6264E02_WP;
    pc%vapPressure%B = -1.2413E04_WP;
    pc%vapPressure%C = -2.0383E01_WP;
    pc%vapPressure%D =  1.0482E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 195.56_WP, 588.10_WP /);
    pc%vapPressure%range(:,2) = (/ 3.4524E-04_WP, 3.9062E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 8.3972E07_WP;
    pc%HeatOfVapMol%B = 5.3507E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 195.56_WP, 588.10_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.7638E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.8050E04_WP;
    pc%solHeatCapacityMol%B = 1.6164E03_WP;
    pc%solHeatCapacityMol%C = -9.4040E00_WP;
    pc%solHeatCapacityMol%D = 2.4970E-0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 31.56_WP, 178.06_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.4382E04_WP, 1.1258E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.0120E05_WP;
    pc%liqHeatCapacityMol%B = -6.5130E02_WP;
    pc%liqHeatCapacityMol%C = 2.2750E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 200.14_WP, 389.15_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.6198E05_WP, 2.9227E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 9.0600E04_WP;
    pc%IG_HeatCapacityMol%B = 3.0620E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6054E03_WP;
    pc%IG_HeatCapacityMol%D = 2.1150E05_WP;
    pc%IG_HeatCapacityMol%E = -7.1797E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 9.8902E04_WP, 3.4133E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.1875E-01_WP;
    pc%SecondViralCoefMol%B = -2.1461E02_WP;
    pc%SecondViralCoefMol%C = -5.9275E07_WP;
    pc%SecondViralCoefMol%D = -1.3203E20_WP;
    pc%SecondViralCoefMol%E = 9.8814E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 294.05_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.6034E00_WP, 5.8109E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.2758E01_WP;
    pc%liqViscosity%B = 2.9169E03_WP;
    pc%liqViscosity%C = 1.2839E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 253.15_WP, 420.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.6075E-02_WP, 3.1659E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.8903E-07_WP;
    pc%vapViscosity%B = 7.0310E-01_WP;
    pc%vapViscosity%C = 1.7590E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 410.95_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 9.1106E-06_WP, 2.0676E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0060E-01_WP;
    pc%liqThermalConductivity%B = -1.6030E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 273.15_WP, 353.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5681E-01_WP, 1.4399E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.8960E03_WP;
    pc%vapThermalConductivity%B = 8.9850E-01_WP;
    pc%vapThermalConductivity%C = 1.2735E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 410.95_WP, 990.95_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.0849E-02_WP, 1.1087E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0660E-02_WP;
    pc%SurfaceTension%B = -8.3920E-05_WP;
    pc%SurfaceTension%C = -3.8300E-09_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 283.15_WP, 373.15_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.6591E-02_WP, 1.8812E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_1pentanol


  subroutine pc_1propanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1propanol'; 
    pc%CAS = '71-23-8';
    ! components
    pc%ChemicalFormula%C = 3;         
    pc%ChemicalFormula%H = 8;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 60.096_WP;   ! kg/kmol    
    pc%Tcrit = 536.78_WP;             ! K    
    pc%pcrit = 5.1750E+06_WP;         ! Pa 
    pc%vcrit = 0.21900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.254_WP; ! -


    pc%MeltingPoint = 146.95_WP;              ! K
    pc%TripletPointT = 146.95_WP;             ! K
    pc%TripletPointp = 3.0828E-07_WP;         ! Pa
    pc%NormalBoilingPoint = 370.35_WP;        ! K
    pc%liqMolarVolume = 0.075089_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.5520E+08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -1.5990E+08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy =  3.2247E+05_WP;            ! J/kmolK
    pc%HeatFusionAtMeltPoint = 5.3720E+06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.8438E+09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK


    pc%AcentricFactor = 0.6218_WP;            ! -
    pc%RadiusOfGyration = 2.359E-10_WP;       ! m
    pc%SolubilityParameter = 2.4450E+04_WP;   ! (J/m^3)^.5 
    pc%DipoleMoment = 5.6000E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.04217_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 6.2800E+08_WP;       ! m^2/kmol
    pc%RefractiveIndex = 1.3837_WP;           ! -
    pc%FlashPoint = 288.15_WP;                ! K
    pc%FlammabilityLimits = (/ 2.0_WP, 12.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 644.00_WP;          ! K


    pc%solDensityMol%A = 1.7800E+01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 146.95_WP, 1.7800E+01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.2350E+00_WP;
    pc%liqDensityMol%B = 2.7136E-01_WP;
    pc%liqDensityMol%C = 5.3678E+02_WP;
    pc%liqDensityMol%D = 2.4000E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 146.95_WP, 536.78_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.5231E+01_WP, 4.5511E+00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  8.8134E+01_WP;
    pc%vapPressure%B = -8.4986E+03_WP;
    pc%vapPressure%C = -9.0766E+00_WP;
    pc%vapPressure%D =  8.3303E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 146.95_WP, 536.78_WP /);
    pc%vapPressure%range(:,2) = (/ 3.0828E-07_WP, 5.1214E+06_WP /);
    pc%vapPressure%eqn = 101;

    pc%HeatOfVapMol%A = 6.3300E+07_WP;
    pc%HeatOfVapMol%B = 3.5750E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 146.95_WP, 536.78_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.6460E+07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 8.7700E+03_WP;
    pc%solHeatCapacityMol%B = 8.1580E+02_WP;
    pc%solHeatCapacityMol%C = 2.6180E+00_WP;
    pc%solHeatCapacityMol%D = -7.6720E-02_WP;
    pc%solHeatCapacityMol%E = 3.3530E-04_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 140.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.8141E+03_WP, 7.5044E+04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.5876E+05_WP;
    pc%liqHeatCapacityMol%B = -6.3500E+02_WP;
    pc%liqHeatCapacityMol%C = 1.9690E+00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 146.95_WP, 400.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.0797E+05_WP, 2.1980E+05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 6.1900E+04_WP;
    pc%IG_HeatCapacityMol%B = 2.0213E+05_WP;
    pc%IG_HeatCapacityMol%C = 1.6293E+03_WP;
    pc%IG_HeatCapacityMol%D = 1.2956E+05_WP;
    pc%IG_HeatCapacityMol%E = 7.2740E+02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 6.6651E+04_WP, 2.2458E+05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.6200E-02_WP;
    pc%SecondViralCoefMol%B = -5.5600E+01_WP;
    pc%SecondViralCoefMol%C = -3.3370E+07_WP;
    pc%SecondViralCoefMol%D = -3.2000E+18_WP;
    pc%SecondViralCoefMol%E = -2.5610E+22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 295.19_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.9998E+00_WP, -7.5487E-04_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.7528E02_WP;
    pc%liqViscosity%B = 6.1044E03_WP;
    pc%liqViscosity%C = 2.8761E01_WP;
    pc%liqViscosity%D = -5.1323E-02_WP;
    pc%liqViscosity%E = 1.0000E00;
    pc%liqViscosity%range(:,1) = (/ 146.95_WP, 370.35_WP /);
    pc%liqViscosity%range(:,2) = (/ 9.3820E00_WP, 4.5317E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 7.9420E-07_WP;
    pc%vapViscosity%B = 5.4910E-01_WP;
    pc%vapViscosity%C = 4.1580E+02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 200.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.7317E-06_WP, 2.4902E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2030E-01_WP;
    pc%liqThermalConductivity%B = -2.1550E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 200.00_WP, 370.35_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7720E-01_WP, 1.4049E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -6.1384E+02_WP;
    pc%vapThermalConductivity%B = 7.9270E-01_WP;
    pc%vapThermalConductivity%C = -1.1574E+09_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 372.13_WP, 720.25_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1531E-02_WP, 7.0337E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 4.5330E-02_WP;
    pc%SurfaceTension%B = -6.8800E-05_WP;
    pc%SurfaceTension%C = -1.6000E-08_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 283.15_WP, 370.35_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.4566E-02_WP, 1.7655E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_1propanol


  subroutine pc_1tetradecanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1tetradecanol'; 
    pc%CAS = '112-72-1';
    ! components
    pc%ChemicalFormula%C = 14;         
    pc%ChemicalFormula%H = 30;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 214.392_WP;  ! kg/kmol    
    pc%Tcrit = 747.00_WP;             ! K    
    pc%pcrit = 1.8100E06_WP;          ! Pa 
    pc%vcrit = 0.80200_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.234_WP; ! -


    pc%MeltingPoint = 310.65_WP;              ! K
    pc%TripletPointT = 310.65_WP;             ! K
    pc%TripletPointp = 7.9998E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 568.80_WP;        ! K
    pc%liqMolarVolume = 0.260434_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -4.7590E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -6.2200E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  7.5390E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.5100E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -8.4914E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.7432_WP;            ! -
    pc%RadiusOfGyration = 6.7300E-10_WP;      ! m
    pc%SolubilityParameter = 1.9430E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.1702E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.15470_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 2.1130E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4454_WP;           ! -
    pc%FlashPoint = 414.15_WP;                ! K
    pc%FlammabilityLimits = (/ 0.5_WP, 5.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 550.00_WP;          ! K


    pc%solDensityMol%A = 4.4513E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 294.00_WP, 298.15_WP /);
    pc%solDensityMol%range(:,2) = (/ 4.4513E00_WP, 4.4513E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.2737E-01_WP;
    pc%liqDensityMol%B = 2.6236E-01_WP;
    pc%liqDensityMol%C = 7.4700E02_WP;
    pc%liqDensityMol%D = 3.2417E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 310.065_WP, 747.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 3.8397E00_WP, 1.2478E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  4.1146E02_WP;
    pc%vapPressure%B = -2.6293E04_WP;
    pc%vapPressure%C = -5.9876E01_WP;
    pc%vapPressure%D =  4.5908E-02_WP;
    pc%vapPressure%E =  1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 310.65_WP, 747.00_WP /);
    pc%vapPressure%range(:,2) = (/ 7.9998E-02_WP, 1.8151E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 2.6345E08_WP;
    pc%HeatOfVapMol%B = 2.6157E00_WP;
    pc%HeatOfVapMol%C = -1.9979E00_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 310.65_WP, 747.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 1.0092E08_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -4.4453E05_WP;
    pc%solHeatCapacityMol%B = 2.9100E03_WP;
    pc%solHeatCapacityMol%C = -8.0430E00_WP;
    pc%solHeatCapacityMol%D = 2.3700E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 269.00_WP, 310.65_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.3826E05_WP, 4.5946E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = -4.4041E06_WP;
    pc%liqHeatCapacityMol%B = 3.7646E04_WP;
    pc%liqHeatCapacityMol%C = -1.0306E02_WP;
    pc%liqHeatCapacityMol%D = 1.2078E-01_WP;
    pc%liqHeatCapacityMol%E = -4.9723E-05_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 324.00_WP, 568.80_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 5.3443E05_WP, 6.8753E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.3340E05_WP;
    pc%IG_HeatCapacityMol%B = 7.2700E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5623E03_WP;
    pc%IG_HeatCapacityMol%D = 5.4210E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2300E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.3668E05_WP, 8.4672E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 7.2937E-01_WP;
    pc%SecondViralCoefMol%B = - 9.0895E02_WP;
    pc%SecondViralCoefMol%C = -2.3210E08_WP;
    pc%SecondViralCoefMol%D = 1.0921E22_WP;
    pc%SecondViralCoefMol%E = -4.4086E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 373.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -8.4889E00_WP, 5.4944E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1007E01_WP;
    pc%liqViscosity%B = 3.5413E03_WP;
    pc%liqViscosity%C = -8.1269E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 310.65_WP, 568.80_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.3966E-02_WP, 4.8381E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 9.3882E-07_WP;
    pc%vapViscosity%B = 4.7630E-01_WP;
    pc%vapViscosity%C = 8.4022E02_WP;
    pc%vapViscosity%D = -5.3350E04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 310.65_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.5823E-06_WP, 1.4106E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.3150E-01_WP;
    pc%liqThermalConductivity%B = -1.9870E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 310.65_WP, 568.80_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6977E-01_WP, 1.1848E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -3.0040E00_WP;
    pc%vapThermalConductivity%B = -3.4874E-02_WP;
    pc%vapThermalConductivity%C = - 1.7663E03_WP;
    pc%vapThermalConductivity%D = -2.8478E07_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 568.80_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.6716E-02_WP, 8.0730E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.9479E-02_WP;
    pc%SurfaceTension%B = 1.3005E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 310.65_WP, 747.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.9561E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1tetradecanol


  subroutine pc_1tridecanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1tridecanol'; 
    pc%CAS = '112-70-9';
    ! components
    pc%ChemicalFormula%C = 13;         
    pc%ChemicalFormula%H = 28;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 200.365_WP;  ! kg/kmol    
    pc%Tcrit = 734.00_WP;             ! K    
    pc%pcrit = 1.9350E06_WP;          ! Pa 
    pc%vcrit = 0.74900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.237_WP; ! -


    pc%MeltingPoint = 303.75_WP;              ! K
    pc%TripletPointT = 303.75_WP;             ! K
    pc%TripletPointp = 9.4584E-02_WP;          ! Pa
    pc%NormalBoilingPoint = 553.60_WP;        ! K
    pc%liqMolarVolume = 0.241531_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -4.6040E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -7.5660E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  7.1470E05_WP;            ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.3300E07_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -7.9014E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.7124_WP;            ! -
    pc%RadiusOfGyration = 6.4170E-10_WP;      ! m
    pc%SolubilityParameter = 1.9760E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.500E-30_WP;           ! C^m
    pc%VanDerWaalsVolume = 0.14447_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.9780E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4433_WP;           ! -
    pc%FlashPoint = 394.15_WP;                ! K
    pc%FlammabilityLimits = (/ 0.6_WP, 5.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 550.00_WP;          ! K


    pc%solDensityMol%A = 4.7913E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 294.00_WP, 4.7213E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.5093E-01_WP;
    pc%liqDensityMol%B = 2.6286E-01_WP;
    pc%liqDensityMol%C = 7.3400E02_WP;
    pc%liqDensityMol%D = 3.1074E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 303.75_WP, 734.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 4.1402E00_WP, 1.3320E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  2.6599E02_WP;
    pc%vapPressure%B = -2.1388E04_WP;
    pc%vapPressure%C = -3.4855E01_WP;
    pc%vapPressure%D =  1.4125E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 303.75_WP, 734.00_WP /);
    pc%vapPressure%range(:,2) = (/ 9.4584E-02_WP, 1.9240E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 2.3170E08_WP;
    pc%HeatOfVapMol%B = 2.3543E00_WP;
    pc%HeatOfVapMol%C = -1.7420E00_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 303.75_WP, 734.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 9.6831E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 3.8300E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 293.15_WP, 293.15_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.8300E08_WP, 3.8300E08_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.0247E06_WP;
    pc%liqHeatCapacityMol%B = -3.2237E04_WP;
    pc%liqHeatCapacityMol%C = 1.3859E02_WP;
    pc%liqHeatCapacityMol%D = -2.4455E-01_WP;
    pc%liqHeatCapacityMol%E = 1.5407E-04_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 303.75_WP, 553.60_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 4.7757E05_WP, 6.3232E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.1780E05_WP;
    pc%IG_HeatCapacityMol%B = 6.8080E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5702E03_WP;
    pc%IG_HeatCapacityMol%D = 5.0750E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2440E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.2082E05_WP, 7.2076E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.4000E-02_WP;
    pc%SecondViralCoefMol%B = -5.5700E01_WP;
    pc%SecondViralCoefMol%C = -1.2900E07_WP;
    pc%SecondViralCoefMol%D = -6.4900E19_WP;
    pc%SecondViralCoefMol%E = -2.5800E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 367.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -8.2103E00_WP, 5.3664E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.6661E01_WP;
    pc%liqViscosity%B = 3.7541E03_WP;
    pc%liqViscosity%C = 1.6684E-02_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 303.75_WP, 553.60_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.4899E-02_WP, 5.6891E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.2436E-08_WP;
    pc%vapViscosity%B = 8.8190E-01_WP;
    pc%vapViscosity%C = -4.0200E00_WP;
    pc%vapViscosity%D = 2.7480E04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 303.75_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.3038E-06_WP, 1.4531E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1542E-01_WP;
    pc%liqThermalConductivity%B = -2.1522E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 303.75_WP, 553.60_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5005E-01_WP, 9.6274E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 3.4379E03_WP;
    pc%vapThermalConductivity%B = 9.5884E-01_WP;
    pc%vapThermalConductivity%C = 3.0868E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 553.60_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.6319E-02_WP, 8.3812E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.9791E-02_WP;
    pc%SurfaceTension%B = 1.2646E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 303.75_WP, 734.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.0428E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1tridecanol


  subroutine pc_1undecanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '1undecanol';
    pc%CAS = '112-42-5';
    ! components
    pc%ChemicalFormula%C = 11;         
    pc%ChemicalFormula%H = 24;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 172.311_WP;   ! kg/kmol    
    pc%Tcrit = 703.60_WP;              ! K    
    pc%pcrit = 2.1470E06_WP;           ! Pa 
    pc%vcrit = 0.71800_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.264_WP;  ! -


    pc%MeltingPoint = 288.45_WP;              ! K
    pc%TripletPointT = 288.45_WP;             ! K
    pc%TripletPointp = 1.2889E-01_WP          ! Pa
    pc%NormalBoilingPoint = 518.15_WP;        ! K
    pc%liqMolarVolume = 0.207672_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -4.1840E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -9.1530E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  6.3620E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.7110E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -6.7261E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.2220E08_WP;              ! J/kmol
    pc%StandardHeatOfFormation = -5.0480E08_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = 4.4920E05_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.6247_WP;            ! -
    pc%RadiusOfGyration = 5.8080E-10_WP;      ! m
    pc%SolubilityParameter = 1.9920E04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.12401_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.7080E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4386_WP;           ! -
    pc%FlashPoint = 383.15_WP;                ! K
    pc%FlammabilityLimits = (/ 0.7_WP, 5.2_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 387.0_WP, 428.0_WP /); ! K
    pc%AutoignitionTemp = 538.00_WP;          ! K


    pc%solDensityMol%A = 5.6800E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 288.45_WP, 5.6800E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.3013E-01_WP;
    pc%liqDensityMol%B = 2.3703E-01_WP;
    pc%liqDensityMol%C = 7.0360E02_WP;
    pc%liqDensityMol%D = 2.7001E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 288.45_WP, 703.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 4.8537E00_WP, 1.3928E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  2.2521E02_WP;
    pc%vapPressure%B = -1.8532E04_WP;
    pc%vapPressure%C = -2.8939E01_WP;
    pc%vapPressure%D =  1.0981E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 288.45_WP, 703.60_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2889E-01_WP, 2.1419E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 1.2015E08_WP;
    pc%HeatOfVapMol%B = 6.3047E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 288.46_WP, 703.60_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 8.6152E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 3.3000E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 289.05_WP, 289.05_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.3000E05_WP, 3.3000E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.2945E05_WP;
    pc%liqHeatCapacityMol%B = -3.0395E03_WP;
    pc%liqHeatCapacityMol%C = 2.7927E01_WP;
    pc%liqHeatCapacityMol%D = -6.1847E-02_WP;
    pc%liqHeatCapacityMol%E = 4.3042E-05_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 289.05_WP, 520.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.9103E05_WP, 5.5125E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.8590E05_WP;
    pc%IG_HeatCapacityMol%B = 5.8690E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5718E03_WP;
    pc%IG_HeatCapacityMol%D = 4.3260E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2270E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.0232E05_WP, 6.7834E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 5.5925E-01_WP;
    pc%SecondViralCoefMol%B = -7.1346E02_WP;
    pc%SecondViralCoefMol%C = -1.2861E08_WP;
    pc%SecondViralCoefMol%D = 5.6156E21_WP;
    pc%SecondViralCoefMol%E = -2.1970E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 351.80_WP, -1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.1054E00_WP, 4.5665E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -7.3477E01_WP;
    pc%liqViscosity%B = 5.9928E03_WP;
    pc%liqViscosity%C = 8.6131E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 288.45_WP, 590.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.0000E-02_WP, 2.3230E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 5.9537E-08_WP;
    pc%vapViscosity%B = 8.1842E-01_WP;
    pc%vapViscosity%C = 9.0245E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 288.45_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.6767E-06_WP, 1.5578E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1211E-01_WP;
    pc%liqThermalConductivity%B = -2.1815E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 288.45_WP, 561.20_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.4918E-01_WP, 8.9684E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.4988E03_WP;
    pc%vapThermalConductivity%B = 9.5209E-01_WP;
    pc%vapThermalConductivity%C = 2.0167E10_WP;
    pc%vapThermalConductivity%D = -2.6800E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 518.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.4658E-02_WP, 8.8994E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.3432E-02_WP;
    pc%SurfaceTension%B = 1.1218E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 288.45_WP, 703.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.9565E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_1undecanol


  subroutine pc_2butanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2butanol';
    pc%CAS = '78-92-2';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 74.123_WP;   ! kg/kmol    
    pc%Tcrit = 536.05_WP;             ! K    
    pc%pcrit = 4.1790E06_WP;           ! Pa 
    pc%vcrit = 0.26900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.252_WP; ! -


    pc%MeltingPoint = 158.45_WP;              ! K
    pc%TripletPointT = 158.45_WP;             ! K
    pc%TripletPointp = 1.1323E-06_WP;          ! Pa
    pc%NormalBoilingPoint = 372.70_WP;        ! K
    pc%liqMolarVolume = 0.092017_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.9290E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.6960E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.6469E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  5.9710E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.4408E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.5722_WP;            ! -
    pc%RadiusOfGyration = 3.2030E-10_WP;      ! m
    pc%SolubilityParameter = 2.2540E04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05239_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.6200E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3949_WP;           ! -
    pc%FlashPoint = 297.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.7_WP, 9.8_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 679.00_WP;          ! K


    pc%solDensityMol%A = 1.4700E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 158.45_WP, 1.4700E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.6600E-01_WP;
    pc%liqDensityMol%B = 2.6064E-01_WP;
    pc%liqDensityMol%C = 5.3605E02_WP;
    pc%liqDensityMol%D = 2.7460E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 158.45_WP, 536.05_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.2570E01_WP, 3.7063E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A = 1.5254E02_WP;
    pc%vapPressure%B = -1.1111E04_WP;
    pc%vapPressure%C = -1.9025E01_WP;
    pc%vapPressure%D = 1.0426E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 158.45_WP, 536.05_WP /);
    pc%vapPressure%range(:,2) = (/ 1.1323E-06_WP, 4.2014E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.2560E07_WP;
    pc%HeatOfVapMol%B = 4.7740E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 158.45_WP, 536.05_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.1383E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.0530E04_WP;
    pc%solHeatCapacityMol%B = 1.2313E03_WP;
    pc%solHeatCapacityMol%C = -6.8730E00_WP;
    pc%solHeatCapacityMol%D = 1.9800E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 150.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 6.4599E03_WP, 8.6347E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.0670E05_WP;
    pc%liqHeatCapacityMol%B = -1.0204E03_WP;
    pc%liqHeatCapacityMol%C = 3.2900E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 158.45_WP, 372.70_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.2762E05_WP, 2.8340E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 8.2020E04_WP;
    pc%IG_HeatCapacityMol%B = 2.5220E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6010E03_WP;
    pc%IG_HeatCapacityMol%D = 1.5864E05_WP;
    pc%IG_HeatCapacityMol%E = -7.0415E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 8.8896E04_WP, 2.8513E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.8777E-01_WP;
    pc%SecondViralCoefMol%B = -1.7264E02_WP;
    pc%SecondViralCoefMol%C = -3.6937E07_WP;
    pc%SecondViralCoefMol%D = -5.7699E19_WP;
    pc%SecondViralCoefMol%E = 4.7178E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 268.03_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.8800E00_WP, 6.1730E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.0638E02_WP;
    pc%liqViscosity%B = 7.4343E03_WP;
    pc%liqViscosity%C = 1.3285E01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 288.15_WP, 373.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 4.7932E-03_WP, 4.1659E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.2534E-07_WP;
    pc%vapViscosity%B = 7.6560E-01_WP;
    pc%vapViscosity%C = 9.2847E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 158.45_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.8197E-06_WP, 2.2716E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0460E-01_WP;
    pc%liqThermalConductivity%B = -2.3542E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 200.00_WP, 398.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5752E-01_WP, 1.1087E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.0086E03_WP;
    pc%vapThermalConductivity%B = 9.7450E-01_WP;
    pc%vapThermalConductivity%C = 1.5008E10_WP;
    pc%vapThermalConductivity%D = -9.2508E11_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 372.70_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.9153E-02_WP, 1.1959E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 4.6450E-02_WP;
    pc%SurfaceTension%B = -7.4730E-05_WP;
    pc%SurfaceTension%C = -1.3000E-08_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 273.15_WP, 372.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.5068E-02_WP, 1.6792E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_2butanol


  subroutine pc_2butanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2butanone';
    pc%CAS = '78-93-3';
    ! components
    pc%ChemicalFormula%C = 2;         
    pc%ChemicalFormula%H = 6;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 72.107_WP;   ! kg/kmol    
    pc%Tcrit = 535.50_WP;             ! K    
    pc%pcrit = 4.1543E06_WP;           ! Pa 
    pc%vcrit = 0.26700_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.249_WP; ! -


    pc%MeltingPoint = 186.48_WP;              ! K
    pc%TripletPointT = 186.48_WP;             ! K
    pc%TripletPointp = 1.3941E00_WP;          ! Pa
    pc%NormalBoilingPoint = 352.79_WP;        ! K
    pc%liqMolarVolume = 0.090204_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.3495E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.4606E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.3811E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  8.4391E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.2616E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.3241_WP;            ! -
    pc%RadiusOfGyration = 3.1350E-10_WP;      ! m
    pc%SolubilityParameter = 1.8796E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 9.2064E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.04927_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.1900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3764_WP;           ! -
    pc%FlashPoint = 267.04_WP;                ! K
    pc%FlammabilityLimits = (/ 1.8_WP, 10.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 788.71_WP;          ! K


    pc%solDensityMol%A = 1.4700E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 186.48_WP, 1.4700E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.3300E-01_WP;
    pc%liqDensityMol%B = 2.514E-01_WP;
    pc%liqDensityMol%C = 5.3550E02_WP;
    pc%liqDensityMol%D = 2.8570E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 186.48_WP, 535.50_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.2593E01_WP, 3.7112E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.1474E02_WP;
    pc%vapPressure%B = -7.1300E03_WP;
    pc%vapPressure%C = -1.5184E01_WP;
    pc%vapPressure%D =  1.7234E-02_WP;
    pc%vapPressure%E =  1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 186.48_WP, 535.50_WP /);
    pc%vapPressure%range(:,2) = (/ 1.3941E00_WP, 4.1976E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.1400E07_WP;
    pc%HeatOfVapMol%B = 7.5450E-01_WP;
    pc%HeatOfVapMol%C = -5.9500E-01_WP;
    pc%HeatOfVapMol%D = 2.3300E-01_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 186.48_WP, 535.50_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.0176E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -9.712E03_WP;
    pc%solHeatCapacityMol%B = 1.3640E03_WP;
    pc%solHeatCapacityMol%C = -7.1303E00_WP;
    pc%solHeatCapacityMol%D = 1.6570E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 10.00_WP, 186.48_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.2315E03_WP, 1.0415E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.3230E05_WP;
    pc%liqHeatCapacityMol%B = 2.0087E02_WP;
    pc%liqHeatCapacityMol%C = -9.5970E-01_WP;
    pc%liqHeatCapacityMol%D = 1.9533E-03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 186.48_WP, 373.15_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.4905E05_WP, 1.7511E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 7.8400E04_WP;
    pc%IG_HeatCapacityMol%B = 2.1032E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5488E03_WP;
    pc%IG_HeatCapacityMol%D = 1.1855E05_WP;
    pc%IG_HeatCapacityMol%E = 6.9300E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 8.3967E04_WP, 2.4816E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 9.0100E-02_WP;
    pc%SecondViralCoefMol%B = -8.6700E01_WP;
    pc%SecondViralCoefMol%C = -6.0500E07_WP;
    pc%SecondViralCoefMol%D = 3.4650E20_WP;
    pc%SecondViralCoefMol%E = -1.0330E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 277.75_WP, 2677.50_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.7632E00_WP, 5.4567E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -8.1340E-01_WP;
    pc%liqViscosity%B = 5.0978E02_WP;
    pc%liqViscosity%C = -1.5324E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 188.71_WP, 533.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.1504E-03_WP, 7.6450E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.9759E-08_WP;
    pc%vapViscosity%B = 1.0252E00_WP;
    pc%vapViscosity%C = -1.6480E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 273.15_WP, 573.15_WP /);
    pc%vapViscosity%range(:,2) = (/ 6.6159E-06_WP, 1.3684E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1920E-01_WP;
    pc%liqThermalConductivity%B = -2.4840E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 186.48_WP, 352.79_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7288E-01_WP, 1.3157E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -4.2027E06_WP;
    pc%vapThermalConductivity%B = -1.5240E-01_WP;
    pc%vapThermalConductivity%C = 2.0846E09_WP;
    pc%vapThermalConductivity%D = -1.4577E13_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 352.79_WP, 992.79_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.5457E-02_WP, 1.1571E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.9441E-02_WP;
    pc%SurfaceTension%B = 1.1165E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 186.48_WP, 355.50_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.4101E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2butanone


  subroutine pc_2heptanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2heptanone';
    pc%CAS = '110-43-0';
    ! components
    pc%ChemicalFormula%C = 7;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 114.188_WP;   ! kg/kmol    
    pc%Tcrit = 611.40_WP;             ! K    
    pc%pcrit = 2.9400E06_WP;           ! Pa 
    pc%vcrit = 0.43400_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.251_WP; ! -


    pc%MeltingPoint = 238.15_WP;              ! K
    pc%TripletPointT = 238.15_WP;             ! K
    pc%TripletPointp = 3.5448E00_WP;          ! Pa
    pc%NormalBoilingPoint = 424.18_WP;        ! K
    pc%liqMolarVolume = 0.140675_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.0045E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.1960E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  2.8064E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.8400E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.0995E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.3480E08_WP;              ! J/kmol
    pc%StandardHeatOfFormation = -3.4820E08_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = 3.4100E05_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.4190_WP;            ! -
    pc%RadiusOfGyration = 4.4980E-10_WP;      ! m
    pc%SolubilityParameter = 1.7800E04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 8.7100E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.07996_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.1240E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4066_WP;           ! -
    pc%FlashPoint = 312.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.1_WP, 7.9_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 311.0_WP, 350.0_WP /); ! K
    pc%AutoignitionTemp = 666.15_WP;          ! K


    pc%solDensityMol%A = 8.8400E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 238.15_WP, 8.8400E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.8247E-01_WP;
    pc%liqDensityMol%B = 2.5279E-01_WP;
    pc%liqDensityMol%C = 6.1140E02_WP;
    pc%liqDensityMol%D = 2.9818E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 238.15_WP, 611.40_WP /);
    pc%liqDensityMol%range(:,2) = (/ 7.5514E00_WP, 2.3042E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.5494E01_WP;
    pc%vapPressure%B = -7.8965E03_WP;
    pc%vapPressure%C = -7.5047E00_WP;
    pc%vapPressure%D =  8.9130E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 238.15_WP, 611.40_WP /);
    pc%vapPressure%range(:,2) = (/ 3.5448E00_WP, 2.9456E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.1425E07_WP;
    pc%HeatOfVapMol%B = 3.9802E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 238.15_WP, 611.40_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.0471E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 2.0390E05;
    pc%solHeatCapacityMol%B = 0;
    pc%solHeatCapacityMol%C = 0;
    pc%solHeatCapacityMol%D = 0;
    pc%solHeatCapacityMol%E = 0;
    pc%solHeatCapacityMol%range(:,1) = (/ 238.15_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.0390E05_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.6504E05_WP;
    pc%liqHeatCapacityMol%B = -3.7568E02_WP;
    pc%liqHeatCapacityMol%C = 1.0024E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 238.15_WP, 490.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.3242E05_WP, 3.2163E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.2507E05_WP;
    pc%IG_HeatCapacityMol%B = 2.1480E05_WP;
    pc%IG_HeatCapacityMol%C = 6.9120E02_WP;
    pc%IG_HeatCapacityMol%D = 1.6190E05_WP;
    pc%IG_HeatCapacityMol%E = 1.7593E03_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1200.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.2688E05_WP, 3.8446E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.9183E-01_WP;
    pc%SecondViralCoefMol%B = -3.4109E02_WP;
    pc%SecondViralCoefMol%C = -7.8913E07_WP;
    pc%SecondViralCoefMol%D = -4.0026E20_WP;
    pc%SecondViralCoefMol%E = 5.6762E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 305.70_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.3996E00_WP, 4.1041E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.3929E01_WP;
    pc%liqViscosity%B = 1.3219E03_WP;
    pc%liqViscosity%C = 4.0382E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 250.00_WP, 424.18_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.6422E-03_WP, 2.3184E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 8.8629E-08_WP;
    pc%vapViscosity%B = 7.8376E-01_WP;
    pc%vapViscosity%C = 1.0018E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 238.15_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.5495E-06_WP, 1.8088E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1080E-01_WP;
    pc%liqThermalConductivity%B = -2.4600E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 236.15_WP, 424.05_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5222E-01_WP, 1.0648E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.0493E03_WP;
    pc%vapThermalConductivity%B = 1.0323E00_WP;
    pc%vapThermalConductivity%C = 2.2983E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 424.18_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.9506E-02_WP, 1.1145E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.1492E-02_WP;
    pc%SurfaceTension%B = 1.2802E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 236.15_WP, 611.40_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.2692E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2heptanone


  subroutine pc_2hexanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2hexanol';
    pc%CAS = '626-93-7';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 102.77_WP;   ! kg/kmol    
    pc%Tcrit = 586.20_WP;             ! K    
    pc%pcrit = 3.4000E06_WP;          ! Pa 
    pc%vcrit = 0.38000_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.265_WP; ! -


    pc%MeltingPoint = 223.00_WP;              ! K
    pc%TripletPointT = 223.00_WP;             ! K
    pc%TripletPointp = 9.4925E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 413.04_WP;        ! K
    pc%liqMolarVolume = 0.126073_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.3350E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -1.5020E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 4.3600E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.6660E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.5655_WP;            ! -
    pc%RadiusOfGyration = 3.7870E-10_WP;      ! m
    pc%SolubilityParameter = 2.1115E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = MAX_REAL_WP;                 ! C^m
    pc%VanDerWaalsVolume = 0.07285_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.0320E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4128_WP;           ! -
    pc%FlashPoint = 314.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 8.3_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = MAX_REAL_WP;             ! K


    pc%solDensityMol%A = 9.9600E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 223.00_WP, 9.9600E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.9120E-01_WP;
    pc%liqDensityMol%B = 2.6275E-01_WP;
    pc%liqDensityMol%C = 5.8620E02_WP;
    pc%liqDensityMol%D = 2.6945E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 223.00_WP, 586.20_WP /);
    pc%liqDensityMol%range(:,2) = (/ 8.5168E00_WP, 2.6306E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.2842E02_WP;
    pc%vapPressure%B = -1.0914E04_WP;
    pc%vapPressure%C = -1.5189E01_WP;
    pc%vapPressure%D =  5.9636E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 223.00_WP, 586.20_WP /);
    pc%vapPressure%range(:,2) = (/ 9.4925E-02_WP, 3.4076E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 9.1260E07_WP;
    pc%HeatOfVapMol%B = 6.2140E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 223.00_WP, 586.20_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.7778E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 1.9600E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 223.00_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.9600E05_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.6063E05_WP;
    pc%liqHeatCapacityMol%B = -1.9047E03_WP;
    pc%liqHeatCapacityMol%C = 5.2540E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 230.00_WP, 330.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.0049E05_WP, 3.0424E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.2575E05_WP;
    pc%IG_HeatCapacityMol%B = 3.5890E05_WP;
    pc%IG_HeatCapacityMol%C = 1.8215E03_WP;
    pc%IG_HeatCapacityMol%D = 2.5557E05_WP;
    pc%IG_HeatCapacityMol%E = 8.1300E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 300.00_WP, 1200.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.5898E05_WP, 3.7756E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.1711E-01_WP;
    pc%SecondViralCoefMol%B = -7.5358E01_WP;
    pc%SecondViralCoefMol%C = -1.0028E08_WP;
    pc%SecondViralCoefMol%D = 9.8894E20_WP;
    pc%SecondViralCoefMol%E = -4.2519E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 293.10_WP, 2923.10_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.2600E01_WP, 8.7315E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.2444E02_WP;
    pc%liqViscosity%B = 8.3425E03_WP;
    pc%liqViscosity%C = 1.5972E01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 230.00_WP, 413.04_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.6948E-01_WP, 3.2420E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.0200E-07_WP;
    pc%vapViscosity%B = 7.7710E-01_WP;
    pc%vapViscosity%C = 1.0050E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 223.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.6978E-06_WP, 1.9875E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1560E-01_WP;
    pc%liqThermalConductivity%B = -2.6510E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 223.00_WP, 463.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5648E-01_WP, 9.2859E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -7.9554E06_WP;
    pc%vapThermalConductivity%B = -1.6390E-01_WP;
    pc%vapThermalConductivity%C = 2.4484E09_WP;
    pc%vapThermalConductivity%D = -2.4038E13_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 413.04_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1961E-02_WP, 1.1877E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1290E-02_WP;
    pc%SurfaceTension%B = 1.0435E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 223.00_WP, 586.20_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1124E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2hexanol


  subroutine pc_2hexanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2hexanone';
    pc%CAS = '591-78-6';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 100.161_WP;   ! kg/kmol    
    pc%Tcrit = 587.61_WP;             ! K    
    pc%pcrit = 3.2870E06_WP;          ! Pa 
    pc%vcrit = 0.37800_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.254_WP; ! -


    pc%MeltingPoint = 217.35_WP;              ! K
    pc%TripletPointT = 217.35_WP;             ! K
    pc%TripletPointp = 1.4456E00_WP;          ! Pa
    pc%NormalBoilingPoint = 400.70_WP;        ! K
    pc%liqMolarVolume = 0.124128_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.7983E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -1.3008E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 4.1786E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.4900E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.4900E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = -1.3960E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -3.2201E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 3.0811E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.3846_WP;            ! -
    pc%RadiusOfGyration = 4.0930E-10_WP;      ! m
    pc%SolubilityParameter = 1.8050E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 8.9400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06973_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 9.8900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3987_WP;           ! -
    pc%FlashPoint = 296.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 8.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 294.0_WP, 330.0_WP /); ! K
    pc%AutoignitionTemp = 696.15_WP;          ! K


    pc%solDensityMol%A = 1.0200E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 217.35_WP, 1.0200E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.7816E-01_WP;
    pc%liqDensityMol%B = 2.5634E-01_WP;
    pc%liqDensityMol%C = 5.8761E02_WP;
    pc%liqDensityMol%D = 2.8365E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 217.35_WP, 587.61_WP /);
    pc%liqDensityMol%range(:,2) = (/ 8.7319E00_WP, 2.6455E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.0744E02_WP;
    pc%vapPressure%B = -8.5286E03_WP;
    pc%vapPressure%C = -1.2679E01_WP;
    pc%vapPressure%D =  8.4606E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 217.35_WP, 587.61_WP /);
    pc%vapPressure%range(:,2) = (/ 1.4456E00_WP, 3.2859E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.6231E07_WP;
    pc%HeatOfVapMol%B = 3.8207E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 217.35_WP, 587.61_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.7135E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 5.4305E03_WP;
    pc%solHeatCapacityMol%B = 8.2129E02_WP;
    pc%solHeatCapacityMol%C = 5.0967E00_WP;
    pc%solHeatCapacityMol%D = -7.0876E-02_WP;
    pc%solHeatCapacityMol%E = 2.1969E-04_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 39.76_WP, 217.69_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 4.2236E04_WP, 1.8794E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.0825E05_WP;
    pc%liqHeatCapacityMol%B = -1.0747E02_WP;
    pc%liqHeatCapacityMol%C = 2.0620E-01_WP;
    pc%liqHeatCapacityMol%D = 7.0293E-04_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 217.35_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.0185E05_WP, 2.7087E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0940E05_WP;
    pc%IG_HeatCapacityMol%B = 1.8070E05_WP;
    pc%IG_HeatCapacityMol%C = 6.8900E02_WP;
    pc%IG_HeatCapacityMol%D = 1.4740E05_WP;
    pc%IG_HeatCapacityMol%E = 1.7720E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1200.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.1815E05_WP, 3.3207E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.4931E-01_WP;
    pc%SecondViralCoefMol%B = -2.8704E02_WP;
    pc%SecondViralCoefMol%C = -5.8288E07_WP;
    pc%SecondViralCoefMol%D = -2.8661E20_WP;
    pc%SecondViralCoefMol%E = 4.2633E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 293.80_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.5749E00_WP, 4.0669E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1445E01_WP;
    pc%liqViscosity%B = 1.1872E03_WP;
    pc%liqViscosity%C = 2.9076E-03_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 217.35_WP, 400.70_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.5614E-03_WP, 2.1077E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 9.7820E-08_WP;
    pc%vapViscosity%B = 7.7720E-01_WP;
    pc%vapViscosity%C = 9.9530E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 217.35_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.3968E-06_WP, 1.9091E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1076E-01_WP;
    pc%liqThermalConductivity%B = -2.4000E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 217.35_WP, 400.85_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5860E-01_WP, 1.1456E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -1.2158E00_WP;
    pc%vapThermalConductivity%B = 2.6637E-02_WP;
    pc%vapThermalConductivity%C = -1.7116E03_WP;
    pc%vapThermalConductivity%D = -1.3176E07_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.00_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 7.7543E-03_WP, 1.0523E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.2446E-02_WP;
    pc%SurfaceTension%B = 1.2674E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 217.35_WP, 587.61_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.4777E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2hexanone


  subroutine pc_2methyl1butanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2methyl1butanol';
    pc%CAS = '137-32-6';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 88.150_WP;   ! kg/kmol    
    pc%Tcrit = 575.40_WP;             ! K    
    pc%pcrit = 3.9400E06_WP;           ! Pa 
    pc%vcrit = 0.32700_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.269_WP; ! -


    pc%MeltingPoint = 195.00_WP;              ! K
    pc%TripletPointT = 195.00_WP;             ! K
    pc%TripletPointp = 1.2845E-03_WP;         ! Pa
    pc%NormalBoilingPoint = 401.85_WP;        ! K
    pc%liqMolarVolume = 0.108257_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -5.0209E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.4671E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.9351E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  8.2700E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.0620E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.5940E08_WP;              ! J/kmol
    pc%StandardHeatOfFormation = -3.5660E08_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = 2.5330E05_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.5736_WP;            ! -
    pc%RadiusOfGyration = 3.6120E-10_WP;      ! m
    pc%SolubilityParameter = 2.2090E+04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 6.2700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06262_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.9700E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4086_WP;           ! -
    pc%FlashPoint = 316.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.4_WP, 9.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 313.0_WP, 345.0_WP /); ! K
    pc%AutoignitionTemp = 658.15_WP;          ! K


    pc%solDensityMol%A = 1.1800E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 195.0_WP, 1.1800E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.1662E-01_WP;
    pc%liqDensityMol%B = 2.6703E-01_WP;
    pc%liqDensityMol%C = 5.7540E02_WP;
    pc%liqDensityMol%D = 2.4336E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 195.00_WP, 575.40_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0092E01_WP, 3.0582E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  9.6196E01_WP;
    pc%vapPressure%B = -9.6589E03_WP;
    pc%vapPressure%C = -1.0112E01_WP;
    pc%vapPressure%D =  8.4773E-19_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 195.00_WP, 575.40_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2846E-03_WP, 3.9057E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.7920E07_WP;
    pc%HeatOfVapMol%B = 4.6991E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 195.00_WP, 575.40_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.4149E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    !keine Daten in Stoffdatenblatt
    pc%solHeatCapacityMol%A = 0.0_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 8.3025E04_WP;
    pc%liqHeatCapacityMol%B = 4.5971E02_WP;
    pc%liqHeatCapacityMol%C = 0.0_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 195.00_WP, 401.85_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.7267E05_WP, 2.6776E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0890E05_WP;
    pc%IG_HeatCapacityMol%B = 2.1850E05_WP;
    pc%IG_HeatCapacityMol%C = 8.5300E02_WP;
    pc%IG_HeatCapacityMol%D = 1.4000E05_WP;
    pc%IG_HeatCapacityMol%E = 2.9060E03_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 298.15_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3247E05_WP, 3.4718E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.1048E-01_WP;
    pc%SecondViralCoefMol%B = -2.0521E02_WP;
    pc%SecondViralCoefMol%C = -5.4378E07_WP;
    pc%SecondViralCoefMol%D = -2.0676E20_WP;
    pc%SecondViralCoefMol%E = 2.8200E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 287.70_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.1030E00_WP, 5.7554E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -6.4568E01_WP;
    pc%liqViscosity%B = 5.2770E03_WP;
    pc%liqViscosity%C = 7.2659E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 273.15_WP, 403.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.1253E-02_WP, 3.7514E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.0470E-07_WP;
    pc%vapViscosity%B = 7.8178E-01_WP;
    pc%vapViscosity%C = 9.2644E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 195.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.3795E-06_WP, 2.1223E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2553E-01_WP;
    pc%liqThermalConductivity%B = -2.8271E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 195.00_WP, 449.50_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7040E-01_WP, 9.8452E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 7.5670E00_WP;
    pc%vapThermalConductivity%B = -2.7146E-01_WP;
    pc%vapThermalConductivity%C = - 3.8945E03_WP;
    pc%vapThermalConductivity%D = 1.2574E07_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 401.85_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1482E-02_WP, 1.1986E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5638E-02_WP;
    pc%SurfaceTension%B = 1.0827E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 195.00_WP, 575.40_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.5545E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2methyl1butanol


  subroutine pc_2methyl1propanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2methyl1propanol';
    pc%CAS = '78-83-1';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 74.123_WP;   ! kg/kmol    
    pc%Tcrit = 547.78_WP;             ! K    
    pc%pcrit = 4.3000E06_WP;          ! Pa 
    pc%vcrit = 0.27300_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.258_WP; ! -


    pc%MeltingPoint = 165.15_WP;              ! K
    pc%TripletPointT = 165.15_WP;             ! K
    pc%TripletPointp = 6.2545E-07_WP;         ! Pa
    pc%NormalBoilingPoint = 380.81_WP;        ! K
    pc%liqMolarVolume = 0.092913_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.8320E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.5490E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.4800E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  6.3220E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.4490E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.5848_WP;            ! -
    pc%RadiusOfGyration = 3.3320E-10_WP;      ! m
    pc%SolubilityParameter = 2.2910E+04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.4700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05239_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.6200E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3938_WP;           ! -
    pc%FlashPoint = 301.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.7_WP, 10.9_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 681.15_WP;          ! K


    pc%solDensityMol%A = 1.4100E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 165.15_WP, 1.4100E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.9680E-01_WP;
    pc%liqDensityMol%B = 2.7226E-01_WP;
    pc%liqDensityMol%C = 5.4778E02_WP;
    pc%liqDensityMol%D = 2.3890E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 165.15_WP, 547.78_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.2084E01_WP, 3.6612E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.8779E02_WP;
    pc%vapPressure%B = -1.2955E04_WP;
    pc%vapPressure%C = -2.4285E01_WP;
    pc%vapPressure%D =  1.4262E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 165.15_WP, 547.78_WP /);
    pc%vapPressure%range(:,2) = (/ 6.2544E-07_WP, 4.3312E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.6190E07_WP;
    pc%HeatOfVapMol%B = 5.0470E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 165.15_WP, 547.78_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.3570E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.0520E04_WP;
    pc%solHeatCapacityMol%B = 1.1474E03_WP;
    pc%solHeatCapacityMol%C = -5.7430E00_WP;
    pc%solHeatCapacityMol%D = 1.5260E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 160.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 5.4503E+03_WP, 8.8548E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.8020E05_WP;
    pc%liqHeatCapacityMol%B = -7.3800E02_WP;
    pc%liqHeatCapacityMol%C = 2.4930E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 165.15_WP, 380.81_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.2631E05_WP, 2.6069E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 8.7940E04_WP;
    pc%IG_HeatCapacityMol%B = 2.4160E05_WP;
    pc%IG_HeatCapacityMol%C = 1.7180E03_WP;
    pc%IG_HeatCapacityMol%D = 1.6540E05_WP;
    pc%IG_HeatCapacityMol%E = 7.9870E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 300.00_WP, 1200.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.1089E05_WP, 2.6355E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.4800E-01_WP;
    pc%SecondViralCoefMol%B = -1.3400E02_WP;
    pc%SecondViralCoefMol%C = -2.8400E07_WP;
    pc%SecondViralCoefMol%D = 1.7000E18_WP;
    pc%SecondViralCoefMol%E = -7.7520E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 301.25_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.1043E00_WP, 5.0250E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -4.8035E01_WP;
    pc%liqViscosity%B = 4.3067E03_WP;
    pc%liqViscosity%C = 4.8948E00_WP;
    pc%liqViscosity%D = -3.5000E-28_WP;
    pc%liqViscosity%E = 1.0E01_WP;
    pc%liqViscosity%range(:,1) = (/ 210.93_WP, 544.26_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.4105E-01_WP, 4.1665E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.1945E-07_WP;
    pc%vapViscosity%B = 6.4700E-01_WP;
    pc%vapViscosity%C = 2.2314E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 165.15_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.6991E-06_WP, 2.2800E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.8503E-01_WP;
    pc%liqThermalConductivity%B = -1.7855E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 165.15_WP, 380.81_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5554E-01_WP, 1.1704E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 1.8980E03_WP;
    pc%vapThermalConductivity%B = 8.6980E-01_WP;
    pc%vapThermalConductivity%C = 6.9300E09_WP;
    pc%vapThermalConductivity%D = 1.2574E07_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 380.81_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.8322E-02_WP, 1.1142E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 4.7506E-02_WP;
    pc%SurfaceTension%B = -8.6704E-05_WP;
    pc%SurfaceTension%C = 1.0150E-08_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 283.15_WP, 380.81_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.3770E-02_WP, 1.5960E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_2methyl1propanol


  subroutine pc_2methyl2propanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2methyl2propanol';
    pc%CAS = '75-65-0';
    ! components
    pc%ChemicalFormula%C = 2;         
    pc%ChemicalFormula%H = 6;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 74.123_WP;   ! kg/kmol    
    pc%Tcrit = 506.21_WP;             ! K    
    pc%pcrit = 3.9730E6_WP;           ! Pa 
    pc%vcrit = 0.27500_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.260_WP; ! -


    pc%MeltingPoint = 298.97_WP;              ! K
    pc%TripletPointT = 298.97_WP;             ! K
    pc%TripletPointp = 5.9356E03_WP;          ! Pa
    pc%NormalBoilingPoint = 355.57_WP;        ! K
    pc%liqMolarVolume = 0.094744_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.1240E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.7760E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.2630E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  6.7030E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.4239E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.6115_WP;            ! -
    pc%RadiusOfGyration = 3.0670E-10_WP;      ! m
    pc%SolubilityParameter = 2.1600E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05238_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.6200E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3852_WP;           ! -
    pc%FlashPoint = 284.26_WP;                ! K
    pc%FlammabilityLimits = (/ 2.4_WP, 8.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 751.00_WP;          ! K


    pc%solDensityMol%A = 1.3000E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 298.97_WP, 1.3000E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.2120E-01_WP;
    pc%liqDensityMol%B = 2.5440E-01_WP;
    pc%liqDensityMol%C = 5.0621E02_WP;
    pc%liqDensityMol%D = 2.7600E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 298.97_WP, 506.21_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0555E01_WP, 3.6211E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.7231E02_WP;
    pc%vapPressure%B = -1.1590E04_WP;
    pc%vapPressure%C = -2.2118E01_WP;
    pc%vapPressure%D =  1.3709E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 298.97_WP, 506.21_WP /);
    pc%vapPressure%range(:,2) = (/ 5.9356E03_WP, 3.9910E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.7320E07_WP;
    pc%HeatOfVapMol%B = 5.6450E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 298.97_WP, 506.21_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.6703E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 4.6000E03_WP;
    pc%solHeatCapacityMol%B = 4.1200E02_WP;
    pc%solHeatCapacityMol%C = 4.1750E00_WP;
    pc%solHeatCapacityMol%D = -3.4270E-02_WP;
    pc%solHeatCapacityMol%E = 8.0080E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 25.00_WP, 286.13_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.7005E04_WP, 1.9826E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = -9.2546E05_WP;
    pc%liqHeatCapacityMol%B = 7.8949E03_WP;
    pc%liqHeatCapacityMol%C = -1.7661E01_WP;
    pc%liqHeatCapacityMol%D = 1.3617E-02_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 298.96_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.2016E05_WP, 2.9455E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 7.7040E04_WP;
    pc%IG_HeatCapacityMol%B = 2.5390E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5502E03_WP;
    pc%IG_HeatCapacityMol%D = 1.6690E05_WP;
    pc%IG_HeatCapacityMol%E = -6.7930E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 8.5670E04_WP, 2.8508E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.8698E-01_WP;
    pc%SecondViralCoefMol%B = -1.5686E02_WP;
    pc%SecondViralCoefMol%C = -3.2431E07_WP;
    pc%SecondViralCoefMol%D = -3.7866E19_WP;
    pc%SecondViralCoefMol%E =-2.8168E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 253.11_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.0200E00_WP, 7.2796E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.1640E02_WP;
    pc%liqViscosity%B = 1.3205E04_WP;
    pc%liqViscosity%C = 2.9254E01_WP;
    pc%liqViscosity%D = -2.4616E-27_WP;
    pc%liqViscosity%E = 1.0000E01_WP;
    pc%liqViscosity%range(:,1) = (/ 298.15_WP, 451.21_WP /);
    pc%liqViscosity%range(:,2) = (/ 4.3119E-03_WP, 1.0137E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 9.6050E-07_WP;
    pc%vapViscosity%B = 4.8560E-01_WP;
    pc%vapViscosity%C = 3.8100E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 298.97_WP, 600.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 6.7267E-06_WP, 1.3124E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.5960E-01_WP;
    pc%liqThermalConductivity%B = -1.4700E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 293.15_WP, 479.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.1651E-01_WP, 8.9165E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -1.4000E-02_WP;
    pc%vapThermalConductivity%B = 3.6250E-01_WP;
    pc%vapThermalConductivity%C = -2.3760E03_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 355.57_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.0717E-02_WP, 1.2445E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 4.0190E-02_WP;
    pc%SurfaceTension%B = -5.3100E-05_WP;
    pc%SurfaceTension%C = -5.1200E-08_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 298.15_WP, 355.57_WP /);
    pc%SurfaceTension%range(:,2) = (/ 1.9807E-02_WP, 1.4836E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_2methyl2propanol


  subroutine pc_2methylnonane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2methylnonane';
    pc%CAS = '871-83-0';
    ! components
    pc%ChemicalFormula%C = 10;         
    pc%ChemicalFormula%H = 22;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 142.265_WP;  ! kg/kmol    
    pc%Tcrit = 610.00_WP;             ! K    
    pc%pcrit = 2.1200E06_WP;          ! Pa 
    pc%vcrit = 0.58300_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.244_WP; ! -


    pc%MeltingPoint = 198.50_WP;              ! K
    pc%TripletPointT = 198.50_WP;             ! K
    pc%TripletPointp = 2.4868E-03_WP;          ! Pa
    pc%NormalBoilingPoint = 440.15_WP;        ! K
    pc%liqMolarVolume = 0.196780_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.5652E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = 2.8400E07_WP;             ! J/kmol
    pc%IG_AbsoluteEntropy =  5.3802E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.7500E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -6.2888E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.4723_WP;            ! -
    pc%RadiusOfGyration = 5.4290E-10_WP;      ! m
    pc%SolubilityParameter = 1.5440E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = MAX_REAL_WP;                 ! C^m
    pc%VanDerWaalsVolume = 0.10917_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.5030E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4075_WP;           ! -
    pc%FlashPoint = 314.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.7_WP, 5.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 487.15_WP;          ! K


    pc%solDensityMol%A = 6.5300E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 198.50_WP, 6.5300E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 4.5413E-01_WP;
    pc%liqDensityMol%B = 2.6475E-01_WP;
    pc%liqDensityMol%C = 6.100E02_WP;
    pc%liqDensityMol%D = 3.0083E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 198.50_WP, 610.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 5.5854E00_WP, 1.7153E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.3571E02_WP;
    pc%vapPressure%B = -1.0597E04_WP;
    pc%vapPressure%C = -1.6769E01_WP;
    pc%vapPressure%D =  1.0148E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 198.50_WP, 610.00_WP /);
    pc%vapPressure%range(:,2) = (/ 2.4867E-03_WP, 2.1194E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.5700E07_WP;
    pc%HeatOfVapMol%B = 4.2540E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 198.50_WP, 610.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.5570E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 1.6222E06_WP;
    pc%solHeatCapacityMol%B = -5.4504E04_WP;
    pc%solHeatCapacityMol%C = 7.0650E02_WP;
    pc%solHeatCapacityMol%D = -3.9062E00_WP;
    pc%solHeatCapacityMol%E = 7.8986E-03_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 80.0_WP, 198.50_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0703E05_WP, 3.5197E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.0520E05_WP;
    pc%liqHeatCapacityMol%B = -4.7130E02_WP;
    pc%liqHeatCapacityMol%C = 1.6330E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 198.50_WP, 440.15_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.7599E05_WP, 4.1412E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.5390E05_WP;
    pc%IG_HeatCapacityMol%B = 5.6550E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5750E03_WP;
    pc%IG_HeatCapacityMol%D = 3.8870E05_WP;
    pc%IG_HeatCapacityMol%E = 7.0400E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.5677E05_WP, 6.1969E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.4066E-01_WP;
    pc%SecondViralCoefMol%B = -4.3379E02_WP;
    pc%SecondViralCoefMol%C = -9.7744E07_WP;
    pc%SecondViralCoefMol%D = -3.4420E20_WP;
    pc%SecondViralCoefMol%E = 5.9453E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 305.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.5199E00_WP, 2.2494E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0613E01_WP;
    pc%liqViscosity%B = 1.7419E03_WP;
    pc%liqViscosity%C = 1.3450E00_WP;
    pc%liqViscosity%D = -1.2703E-28_WP;
    pc%liqViscosity%E = 1E01_WP;
    pc%liqViscosity%range(:,1) = (/ 198.50_WP, 593.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 8.9011E-03_WP, 5.6989E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 7.7800E-07_WP;
    pc%vapViscosity%B = 4.6360E-01_WP;
    pc%vapViscosity%C = 3.1607E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 198.50_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.4877E-06_WP, 1.4538E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0262E-01_WP;
    pc%liqThermalConductivity%B = -2.4000E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 198.5_WP, 488.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5498E-01_WP, 8.5500E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 3.6874E-05_WP;
    pc%vapThermalConductivity%B = 1.1200E00_WP;
    pc%vapThermalConductivity%C = -1.5392E02_WP;
    pc%vapThermalConductivity%D = 1.7277E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 440.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1849E-02_WP, 8.2911E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1634E-02_WP;
    pc%SurfaceTension%B = 1.2440E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 198.50_WP, 610.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1642E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2methylnonane


  subroutine pc_2methyloctane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2methyloctane';
    pc%CAS = '3221-61 -2';
    ! components
    pc%ChemicalFormula%C = 9;         
    pc%ChemicalFormula%H = 20;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 128.258_WP;  ! kg/kmol    
    pc%Tcrit = 582.80_WP;             ! K    
    pc%pcrit = 2.3100E06_WP;          ! Pa 
    pc%vcrit = 0.54100_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.258_WP; ! -


    pc%MeltingPoint = 192.78_WP;              ! K
    pc%TripletPointT = 192.78_WP;             ! K
    pc%TripletPointp = 1.8181E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 416.43_WP;        ! K
    pc%liqMolarVolume = 0.180734_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.3590E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = 2.0080E7_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.9880E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.8000E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -5.6791E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = 8.2190E06_WP;              ! J/kmol
    pc%StandardHeatOfFormation = -2.8082E08_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = 3.8790E05_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.4583_WP;            ! -
    pc%RadiusOfGyration = 4.9730E-10_WP;      ! m
    pc%SolubilityParameter = 1.53540E4_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.09894_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.3680E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4008_WP;           ! -
    pc%FlashPoint = 296.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.9_WP, 5.4_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 298.0_WP, 333.0_WP /); ! K
    pc%AutoignitionTemp = 493.15_WP;          ! K


    pc%solDensityMol%A = 7.1600E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 192.78_WP, 7.1600E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 4.8369E-01_WP;
    pc%liqDensityMol%B = 2.6168E-01_WP;
    pc%liqDensityMol%C = 5.8280E02_WP;
    pc%liqDensityMol%D = 2.8064E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 192.78_WP, 582.80_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.1230E00_WP, 1.8484E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.8095E01_WP;
    pc%vapPressure%B = -7.7507E03_WP;
    pc%vapPressure%C = -7.9631E00_WP;
    pc%vapPressure%D =  1.4327E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 192.78_WP, 582.80_WP /);
    pc%vapPressure%range(:,2) = (/ 1.8180E-02_WP, 2.3026E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.9830E07_WP;
    pc%HeatOfVapMol%B = 3.9749E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 192.78_WP, 582.80_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.1002E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 3.4830E04_WP;
    pc%solHeatCapacityMol%B = 7.6538E02_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 80.0_WP, 160.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 9.6060E04_WP, 1.5729E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.0720E05_WP;
    pc%liqHeatCapacityMol%B = 2.0160E01_WP;
    pc%liqHeatCapacityMol%C = 7.5925E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 192.78_WP, 500.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.3930E05_WP, 4.0709E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.3743E05_WP;
    pc%IG_HeatCapacityMol%B = 5.1693E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5604E03_WP;
    pc%IG_HeatCapacityMol%D = 3.4726E05_WP;
    pc%IG_HeatCapacityMol%E = 6.9593E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.4022E05_WP, 5.6295E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.0640E-01_WP;
    pc%SecondViralCoefMol%B = -3.8576E02_WP;
    pc%SecondViralCoefMol%C = -7.2988E07_WP;
    pc%SecondViralCoefMol%D = -1.6901E20_WP;
    pc%SecondViralCoefMol%E = 2.6260E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 291.40_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.4846E00_WP, 2.7595E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1283E01_WP;
    pc%liqViscosity%B = 1.1519E03_WP;
    pc%liqViscosity%C = 1.7447E-03_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 192.78_WP, 416.43_WP /);
    pc%liqViscosity%range(:,2) = (/ 4.9985E-03_WP, 2.0218E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 5.3451E-07_WP;
    pc%vapViscosity%B = 5.2741E-01_WP;
    pc%vapViscosity%C = 2.7431E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 192.78_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.5382E-06_WP, 1.6029E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0102E-01_WP;
    pc%liqThermalConductivity%B = -2.4765E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 192.78_WP, 469.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5328E-01_WP, 8.4872E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 4.6511E-05_WP;
    pc%vapThermalConductivity%B = 1.0987E00_WP;
    pc%vapThermalConductivity%C = -1.1797E02_WP;
    pc%vapThermalConductivity%D = 1.6726E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 416.43_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.0894E-02_WP, 8.7652E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.3154E-02_WP;
    pc%SurfaceTension%B = 1.2669E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 192.78_WP, 582.80_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1956E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2methyloctane


  subroutine pc_2nonanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2nonanone';
    pc%CAS = '821-55-6';
    ! components
    pc%ChemicalFormula%C = 9;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 142.241_WP;  ! kg/kmol    
    pc%Tcrit = 652.50_WP;             ! K    
    pc%pcrit = 2.4100E6_WP;           ! Pa 
    pc%vcrit = 0.56000_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.249_WP; ! -


    pc%MeltingPoint = 265.65_WP;              ! K
    pc%TripletPointT = 265.65_WP;             ! K
    pc%TripletPointp = 4.7226E00_WP;          ! Pa
    pc%NormalBoilingPoint = 467.15_WP;        ! K
    pc%liqMolarVolume = 0.173947_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.4070E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.0410E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  5.3560E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.5600E7_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -5.3212E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = - 1.2270E8_WP;              ! J/kmol
    pc%StandardHeatOfFormation = -3.9680E8_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = 4.0830E5_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.4979_WP;            ! -
    pc%RadiusOfGyration = 5.3210E-10_WP;      ! m
    pc%SolubilityParameter = 1.7560E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 9.1400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.10042_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.3940E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4187_WP;           ! -
    pc%FlashPoint = 337.15_WP;                ! K
    pc%FlammabilityLimits = (/ 0.9_WP, 5.9_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 340.0_WP, 380.0_WP /); ! K
    pc%AutoignitionTemp = 633.15_WP;          ! K


    pc%solDensityMol%A = 6.9300E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 265.65_WP, 6.9300E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 4.4319E-01_WP;
    pc%liqDensityMol%B = 2.4818E-01_WP;
    pc%liqDensityMol%C = 6.5250E02_WP;
    pc%liqDensityMol%D = 2.8763E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 265.65_WP, 652.50_WP /);
    pc%liqDensityMol%range(:,2) = (/ 5.9233E00_WP, 1.7858E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.9160E01_WP;
    pc%vapPressure%B = -7.8739E03_WP;
    pc%vapPressure%C = -5.0102E00_WP;
    pc%vapPressure%D =  9.1592E-19_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 265.65_WP, 652.50_WP /);
    pc%vapPressure%range(:,2) = (/ 4.7226E00_WP, 2.4054E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.1729E07_WP;
    pc%HeatOfVapMol%B = 4.0270E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 265.65_WP, 652.50_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.8112E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 2.5700E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 265.65_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.5700E05_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.5400E05_WP;
    pc%liqHeatCapacityMol%B = -5.5822E02_WP;
    pc%liqHeatCapacityMol%C = 1.3135E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 265.65_WP, 510.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.9840E05_WP, 4.1095E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.5347E05_WP;
    pc%IG_HeatCapacityMol%B = 4.2140E05_WP;
    pc%IG_HeatCapacityMol%C = 1.3160E03_WP;
    pc%IG_HeatCapacityMol%D = 2.4060E05_WP;
    pc%IG_HeatCapacityMol%E = 6.4420E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.6947E05_WP, 5.1866E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.8301E-01_WP;
    pc%SecondViralCoefMol%B = -4.4922E02_WP;
    pc%SecondViralCoefMol%C = -1.3676E08_WP;
    pc%SecondViralCoefMol%D = -8.4430E20_WP;
    pc%SecondViralCoefMol%E = 1.1643E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 326.25_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -8.7298E00_WP, 4.2979E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1678E01_WP;
    pc%liqViscosity%B = 1.5184E03_WP;
    pc%liqViscosity%C = -2.9833E-02_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 265.65_WP, 467.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.1793E-03_WP, 1.8209E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 7.3861E-08_WP;
    pc%vapViscosity%B = 7.9750E-01_WP;
    pc%vapViscosity%C = 9.7893E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 265.65_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.6297E-06_WP, 1.6609E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.9050E-01_WP;
    pc%liqThermalConductivity%B = -1.9230E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 265.65_WP, 593.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.3942E-01_WP, 7.6437E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -2.2439E-03_WP;
    pc%vapThermalConductivity%B = 1.0196E00_WP;
    pc%vapThermalConductivity%C = -2.6430E04_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 467.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1276E-02_WP, 1.0103E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0312E-02_WP;
    pc%SurfaceTension%B = 1.0822E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 265.65_WP, 652.50_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.8574E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2nonanone


  subroutine pc_2octanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2octanone';
    pc%CAS = '111-13-7';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 16;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 128.214_WP;  ! kg/kmol    
    pc%Tcrit = 632.70_WP;             ! K    
    pc%pcrit = 2.6400E06_WP;          ! Pa 
    pc%vcrit = 0.49700_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.249_WP; ! -


    pc%MeltingPoint = 252.85_WP;              ! K
    pc%TripletPointT = 252.85_WP;             ! K
    pc%TripletPointp = 4.6773E00_WP;          ! Pa
    pc%NormalBoilingPoint = 446.15_WP;        ! K
    pc%liqMolarVolume = 0.157287_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.2160E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.1380E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.9620E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.4419E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.6984E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.2917E08_WP;              ! J/kmol
    pc%StandardHeatOfFormation = -3.7270E08_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = 3.7384E05_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.4549_WP;            ! -
    pc%RadiusOfGyration = 4.8920E-10_WP;      ! m
    pc%SolubilityParameter = 1.7580E4_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 8.2100E-30_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.09019_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.2590E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4133_WP;           ! -
    pc%FlashPoint = 324.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.0_WP, 6.1_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 326.0_WP, 363.0_WP /); ! K
    pc%AutoignitionTemp = 653.00_WP;          ! K


    pc%solDensityMol%A = 7.7800E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 252.85_WP, 7.7800E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.0006E-01_WP;
    pc%liqDensityMol%B = 2.48511E-01_WP;
    pc%liqDensityMol%C = 6.3270E02_WP;
    pc%liqDensityMol%D = 2.9942E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 192.78_WP, 582.80_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.1230E00_WP, 1.8484E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  6.3775E01_WP;
    pc%vapPressure%B = -7.7113E03_WP;
    pc%vapPressure%C = -5.7359E00_WP;
    pc%vapPressure%D =  3.0902E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 252.85_WP, 632.70_WP /);
    pc%vapPressure%range(:,2) = (/ 4.6773E00_WP, 2.6466E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.5363E07_WP;
    pc%HeatOfVapMol%B = 3.8718E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 252.85_WP, 632.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.3646E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.9270E04_WP;
    pc%solHeatCapacityMol%B = 1.9170E03_WP;
    pc%solHeatCapacityMol%C = -8.2019E00_WP;
    pc%solHeatCapacityMol%D = 1.6844E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 13.00_WP, 252.85_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 4.3019E03_WP, 2.1336E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.0040E05_WP;
    pc%liqHeatCapacityMol%B = -4.2620E02_WP;
    pc%liqHeatCapacityMol%C = 1.1172E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 252.86_WP, 500.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.6406E05_WP, 3.6660E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.3901E05_WP;
    pc%IG_HeatCapacityMol%B = 3.8060E05_WP;
    pc%IG_HeatCapacityMol%C = 1.3717E03_WP;
    pc%IG_HeatCapacityMol%D = 2.2573E05_WP;
    pc%IG_HeatCapacityMol%E = 6.6096E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.4162E05_WP, 4.6547E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.3845E-01_WP;
    pc%SecondViralCoefMol%B = -3.9878E02_WP;
    pc%SecondViralCoefMol%C = -1.0424E08_WP;
    pc%SecondViralCoefMol%D = -5.3291E20_WP;
    pc%SecondViralCoefMol%E = 7.0168E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 316.35_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.3161E00_WP, 4.1692E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1736E01_WP;
    pc%liqViscosity%B = 1.4152E03_WP;
    pc%liqViscosity%C = 3.6180E-04_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 252.85_WP, 446.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.1614E-03_WP, 1.9128E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 8.0901E-08_WP;
    pc%vapViscosity%B = 7.9062E-01_WP;
    pc%vapViscosity%C = 9.9338E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 252.85_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.6110E-06_WP, 1.7325E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1320E-01_WP;
    pc%liqThermalConductivity%B = -2.4940E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 252.85_WP, 499.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5014E-01_WP, 8.8749E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -2.0184E-03_WP;
    pc%vapThermalConductivity%B = 1.0027E00_WP;
    pc%vapThermalConductivity%C = -2.0406E04_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 446.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.0463E-02_WP, 1.0597E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.9062E-02_WP;
    pc%SurfaceTension%B = 1.2967E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 252.85_WP, 632.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.0488E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2octanone


  subroutine pc_2pentanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2pentanol';
    pc%CAS = '6032-29-7';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 88.1497_WP;  ! kg/kmol    
    pc%Tcrit = 560.3_WP;              ! K    
    pc%pcrit = 3.6750E6_WP;           ! Pa 
    pc%vcrit = 0.329_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.26_WP;  ! -


    pc%MeltingPoint = 200.0_WP;              ! K
    pc%TripletPointT = 200.0_WP;             ! K
    pc%TripletPointp = 5.9659E-3_WP;         ! Pa
    pc%NormalBoilingPoint = 392.15_WP;       ! K
    pc%liqMolarVolume = 0.10957_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.1380E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.5930E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.9640E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  8.4800E06_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.0515E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = -1.6850E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -3.6520E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.5500E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.558707_WP;          ! -
    pc%RadiusOfGyration = 3.6190E-10_WP;      ! m
    pc%SolubilityParameter = 2.1700E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.3900E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06262_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.9700E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4044_WP;           ! -
    pc%FlashPoint = 307.0_WP;                 ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 9.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 307.0_WP, 338.0_WP /); ! K
    pc%AutoignitionTemp = 616.0_WP;           ! K


    pc%solDensityMol%A = 1.1700E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.0126E-01_WP;
    pc%liqDensityMol%B = 2.6371E-01_WP;
    pc%liqDensityMol%C = 5.6030E02_WP;
    pc%liqDensityMol%D = 2.5302E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 200.00_WP, 560.30_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0008E01_WP, 3.0384E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.1227E02_WP;
    pc%vapPressure%B = -1.0221E04_WP;
    pc%vapPressure%C = -1.2511E01_WP;
    pc%vapPressure%D =  8.5670E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 200.00_WP, 560.30_WP /);
    pc%vapPressure%range(:,2) = (/ 5.9659E-03_WP, 3.6767E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 8.1306E07_WP;
    pc%HeatOfVapMol%B = 5.3685E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 200.00_WP, 560.30_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.4147E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 0.0_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.7400E05_WP;
    pc%liqHeatCapacityMol%B = -6.5134E01_WP;
    pc%liqHeatCapacityMol%C = 9.5351E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 280.15_WP, 420.22_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.3059E05_WP, 3.1500E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0385E05_WP;
    pc%IG_HeatCapacityMol%B = 2.9140E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6383E03_WP;
    pc%IG_HeatCapacityMol%D = 1.9435E05_WP;
    pc%IG_HeatCapacityMol%E = 7.5094E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 298.15_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3604E05_WP, 3.4084E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.3739E-01_WP;
    pc%SecondViralCoefMol%B = -2.4456E02_WP;
    pc%SecondViralCoefMol%C = -3.2955E07_WP;
    pc%SecondViralCoefMol%D = 3.8515E20_WP;
    pc%SecondViralCoefMol%E = -1.6083E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 280.15_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.1139E00_WP, 6.4596E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.0809E02_WP;
    pc%liqViscosity%B = 7.5182E03_WP;
    pc%liqViscosity%C = 1.3552E01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 250.00_WP, 392.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 4.1157E-02_WP, 3.3865E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.1030E-07_WP;
    pc%vapViscosity%B = 6.4060E-01_WP;
    pc%vapViscosity%C = 2.3282E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 200.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.2711E-06_WP, 2.1022E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2037E-01_WP;
    pc%liqThermalConductivity%B = -2.8791E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 200.00_WP, 435.60_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6279E-01_WP, 9.4956E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.8508E03_WP;
    pc%vapThermalConductivity%B = 9.0265E-01_WP;
    pc%vapThermalConductivity%C = 1.2060E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 392.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.0326E-02_WP, 1.2066E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5779E-02_WP;
    pc%SurfaceTension%B = 1.1409E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 200.00_WP, 560.30_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.3705E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2pentanol


  subroutine pc_2pentanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------



    ! all alternative names
    pc%name = '2pentanone';
    pc%CAS = '107-87-9';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 86.1338_WP;  ! kg/kmol    
    pc%Tcrit = 561.08_WP;             ! K    
    pc%pcrit = 3.6940E06_WP;          ! Pa 
    pc%vcrit = 0.301_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.238_WP; ! -


    pc%MeltingPoint = 196.29_WP;              ! K
    pc%TripletPointT = 196.29_WP;             ! K
    pc%TripletPointp = 0.752347_WP;           ! Pa
    pc%NormalBoilingPoint = 375.46_WP;        ! K
    pc%liqMolarVolume = 0.10725_WP;           ! m^3/kmol


    pc%IG_HeatOfFormation = -2.5920E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.3830E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.7860E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.0620E07_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.8796E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = -1.4520E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -2.9729E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.7410E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.343288_WP;          ! -
    pc%RadiusOfGyration = 3.6180E-10_WP;      ! m
    pc%SolubilityParameter = 1.8290E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 9.2400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.0595_WP;         ! m^3/kmol


    pc%VanDerWaalsArea = 8.5400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.388_WP;            ! -
    pc%FlashPoint = 280.0_WP;                 ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 8.2_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ 278.0_WP, 309.0_WP /); ! K
    pc%AutoignitionTemp = 725.15_WP;          ! K


    pc%solDensityMol%A = 1.2200E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 196.29_WP, 1.2200E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.0411E-01_WP;
    pc%liqDensityMol%B = 2.7207E-01_WP;
    pc%liqDensityMol%C = 5.6108E02_WP;
    pc%liqDensityMol%D = 3.0669E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 196.29_WP, 561.08_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0398E01_WP, 3.3231E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  8.4635E01_WP;
    pc%vapPressure%B = -7.0784E03_WP;
    pc%vapPressure%C = -9.3000E00_WP;
    pc%vapPressure%D =  6.2702E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 196.29_WP, 561.08_WP /);
    pc%vapPressure%range(:,2) = (/ 7.5235E-01_WP, 3.7062E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.1740E07_WP;
    pc%HeatOfVapMol%B = 3.9422E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 196.29_WP, 561.08_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.3663E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -2.5097E04_WP;
    pc%solHeatCapacityMol%B = 2.0373E03_WP;
    pc%solHeatCapacityMol%C = -1.4289E01_WP;
    pc%solHeatCapacityMol%D = 4.4874E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 30.00_WP, 110.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.4373E04_WP, 8.5836E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.9459E05_WP;
    pc%liqHeatCapacityMol%B = -2.6386E02_WP;
    pc%liqHeatCapacityMol%C = 7.6808E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 196.29_WP, 375.46_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.7239E05_WP, 2.0380E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 9.0053E04_WP;
    pc%IG_HeatCapacityMol%B = 2.7085E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6592E03_WP;
    pc%IG_HeatCapacityMol%D = 1.8012E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4396E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 9.5908E04_WP, 3.0797E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.2325E-01_WP;
    pc%SecondViralCoefMol%B = -2.5789E02_WP;
    pc%SecondViralCoefMol%C = -3.5849E07_WP;
    pc%SecondViralCoefMol%D = -3.0892E20_WP;
    pc%SecondViralCoefMol%E = 5.9158E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 280.54_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.8752E00_WP, 4.0691E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1055E01_WP;
    pc%liqViscosity%B = 1.0053E03_WP;
    pc%liqViscosity%C = 3.9301E-03_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 250.00_WP, 375.46_WP /);
    pc%liqViscosity%range(:,2) = (/ 9.0091E-04_WP, 2.3541E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.4630E-07_WP;
    pc%vapViscosity%B = 6.6530E-01_WP;
    pc%vapViscosity%C = 2.0870E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 196.29_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.0030E-06_WP, 2.0186E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1610E-01_WP;
    pc%liqThermalConductivity%B = -2.4866E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 196.29_WP, 375.46_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6729E-01_WP, 1.2274E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -1.7190E-02_WP;
    pc%vapThermalConductivity%B = 4.8320E-01_WP;
    pc%vapThermalConductivity%C = -3.7980E03_WP;
    pc%vapThermalConductivity%D = -1.2350E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.00_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 8.7672E-03_WP, 1.2002E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.3802E-02_WP;
    pc%SurfaceTension%B = 1.0733E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 196.29_WP, 561.08_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.3893E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2pentanone


  subroutine pc_2propanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '2propanol';
    pc%CAS = '67-63-0';
    ! components
    pc%ChemicalFormula%C = 3;         
    pc%ChemicalFormula%H = 8;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 60.096_WP;   ! kg/kmol    
    pc%Tcrit = 508.31_WP;             ! K    
    pc%pcrit = 4.7643E6_WP;           ! Pa 
    pc%vcrit = 0.22013_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.248_WP; ! -


    pc%MeltingPoint = 185.28_WP;              ! K
    pc%TripletPointT = 185.28_WP;             ! K
    pc%TripletPointp = 3.2000E-2_WP;          ! Pa
    pc%NormalBoilingPoint = 355.41_WP;        ! K
    pc%liqMolarVolume = 0.076784_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.7242E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.7339E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.0991E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  5.4099E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.8300E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.6689_WP;            ! -
    pc%RadiusOfGyration = 2.8070E-10_WP;      ! m
    pc%SolubilityParameter = 2.3575E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.5372E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.04216_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 6.2700E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3752_WP;           ! -
    pc%FlashPoint = 284.82_WP;                ! K
    pc%FlammabilityLimits = (/ 2.0_WP, 12.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 672.04_WP;          ! K


    pc%solDensityMol%A = 1.71E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 185.28_WP, 1.7100E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.18E00_WP;
    pc%liqDensityMol%B = 2.6475E-01_WP;
    pc%liqDensityMol%C = 5.0831E02_WP;
    pc%liqDensityMol%D = 2.43E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 185.28_WP, 508.31_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.4656E01_WP, 4.4570E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  9.2935E01_WP;
    pc%vapPressure%B = -8.1771E03_WP;
    pc%vapPressure%C = -1.0031E01_WP;
    pc%vapPressure%D =  3.9988E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 185.28_WP, 508.31_WP /);
    pc%vapPressure%range(:,2) = (/ 3.2000E-02_WP, 4.7673E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.6980E07_WP;
    pc%HeatOfVapMol%B = 8.7E-02_WP;
    pc%HeatOfVapMol%C = 3.0070E-01_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 185.28_WP, 508.31_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.2121E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -8.5E03_WP;
    pc%solHeatCapacityMol%B = 8.7600E02_WP;
    pc%solHeatCapacityMol%C = -8.100E-01_WP;
    pc%solHeatCapacityMol%D = -3.025E-02_WP;
    pc%solHeatCapacityMol%E = 1.4130E-04_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 12.00_WP, 176.66_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.8460E03_WP, 9.1821E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 4.6640E05_WP;
    pc%liqHeatCapacityMol%B = -4.1086E03_WP;
    pc%liqHeatCapacityMol%C = 1.4506E01_WP;
    pc%liqHeatCapacityMol%D = -1.4126E-02_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 188.32_WP, 463.15_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.1277E05_WP, 2.7174E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 4.7460E04_WP;
    pc%IG_HeatCapacityMol%B = 1.9350E05_WP;
    pc%IG_HeatCapacityMol%C = 1.1240E03_WP;
    pc%IG_HeatCapacityMol%D = 9.3800E04_WP;
    pc%IG_HeatCapacityMol%E = 4.6000E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 100.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 4.8262E04_WP, 2.1652E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.0200E-02_WP;
    pc%SecondViralCoefMol%B = -6.2900E00_WP;
    pc%SecondViralCoefMol%C = -4.3100E07_WP;
    pc%SecondViralCoefMol%D = 1.9860E20_WP;
    pc%SecondViralCoefMol%E = -8.6100E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 279.58_WP, 2541.60_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.8942E00_WP, 2.5100E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -8.2300E00_WP;
    pc%liqViscosity%B = 2.2822E03_WP;
    pc%liqViscosity%C = -9.8495E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 187.35_WP, 354.45_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.0037E-01_WP, 5.1389E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.9930E-07_WP;
    pc%vapViscosity%B = 7.2330E-01_WP;
    pc%vapViscosity%C = 1.7800E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 200.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.8683E-06_WP, 2.5019E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0290E-01_WP;
    pc%liqThermalConductivity%B = -2.2780E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 185.28_WP, 410.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6069E-01_WP, 1.0950E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -8.0642E01_WP;
    pc%vapThermalConductivity%B = -1.4549E00_WP;
    pc%vapThermalConductivity%C = -6.0442E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 355.41_WP, 450.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.2388E-02_WP, 3.2427E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 3.818E-02_WP;
    pc%SurfaceTension%B = -3.8180E-05_WP;
    pc%SurfaceTension%C = -6.5100E-08_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 273.15_WP, 355.41_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.2894E-02_WP, 1.6387E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_2propanol


  subroutine pc_3heptanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '3heptanone';
    pc%CAS = '106-35-4';
    ! components
    pc%ChemicalFormula%C = 7;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 114.188_WP;  ! kg/kmol    
    pc%Tcrit = 606.60_WP;             ! K    
    pc%pcrit = 2.9200E6_WP;           ! Pa 
    pc%vcrit = 0.43300_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.251_WP; ! -


    pc%MeltingPoint = 234.15_WP;              ! K
    pc%TripletPointT = 234.15_WP;             ! K
    pc%TripletPointp = 2.2992E00_WP;          ! Pa
    pc%NormalBoilingPoint = 420.55_WP;        ! K
    pc%liqMolarVolume = 0.140257_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.0100E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.2250E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.5800E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.8100E7_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.0980E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.3580E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -3.4860E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 3.4140E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.4076_WP;            ! -
    pc%RadiusOfGyration = 4.5390E-10_WP;      ! m
    pc%SolubilityParameter = 1.7940E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 9.3700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.07996_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.1240E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4066_WP;           ! -
    pc%FlashPoint = 310.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.1_WP, 6.8_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 309.0_WP, 344.0_WP /); ! K
    pc%AutoignitionTemp = 683.00_WP;          ! K


    pc%solDensityMol%A = 8.8600E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 234.15_WP, 8.8600E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.9268E-01_WP;
    pc%liqDensityMol%B = 2.5663E-01_WP;
    pc%liqDensityMol%C = 6.0660E02_WP;
    pc%liqDensityMol%D = 2.7766E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 234.15_WP, 606.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 7.5751E00_WP, 2.3095E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.8463E01_WP;
    pc%vapPressure%B = -8.0772E03_WP;
    pc%vapPressure%C = -7.9062E00_WP;
    pc%vapPressure%D =  8.0521E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 234.15_WP, 606.60_WP /);
    pc%vapPressure%range(:,2) = (/ 2.2992E-04_WP, 2.9190E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.3357E07_WP;
    pc%HeatOfVapMol%B = 4.2167E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 234.15_WP, 606.60_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.1579E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 2.0400E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 234.15_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.0400E05_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.7073E05_WP;
    pc%liqHeatCapacityMol%B = -3.99893E02_WP;
    pc%liqHeatCapacityMol%C = 1.0601E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 234.15_WP, 480.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.3522E05_WP, 3.2303E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.2768E05_WP;
    pc%IG_HeatCapacityMol%B = 3.3810E05_WP;
    pc%IG_HeatCapacityMol%C = 1.3831E03_WP;
    pc%IG_HeatCapacityMol%D = 1.8880E05_WP;
    pc%IG_HeatCapacityMol%E = 6.5030E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3968E05_WP, 4.1386E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.9056E-01_WP;
    pc%SecondViralCoefMol%B = -3.3949E02_WP;
    pc%SecondViralCoefMol%C = -7.6495E07_WP;
    pc%SecondViralCoefMol%D = -4.1164E20_WP;
    pc%SecondViralCoefMol%E = 6.0973E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 303.3_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.5114E00_WP, 4.1554E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -9.3877E00_WP;
    pc%liqViscosity%B = 1.2049E03_WP;
    pc%liqViscosity%C = -3.2618E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 234.15_WP, 421.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.4268E-03_WP, 2.0397E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 8.9656E-08_WP;
    pc%vapViscosity%B = 7.8236E-01_WP;
    pc%vapViscosity%C = 1.0014E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 234.15_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.4848E-06_WP, 1.8122E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0260E-01_WP;
    pc%liqThermalConductivity%B = -2.2340E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 234.15_WP, 553.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5029E-01_WP, 7.9026E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 1.3486E03_WP;
    pc%vapThermalConductivity%B = 1.0313E00_WP;
    pc%vapThermalConductivity%C = 1.4832E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 420.55_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.9429E-02_WP, 1.1287E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.8308E-02_WP;
    pc%SurfaceTension%B = 1.2123E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 234.15_WP, 606.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.2279E-02_WP, 8.3816E-04_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_3heptanone


  subroutine pc_3hexanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '3hexanone';
    pc%CAS = '589-38-8';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 100.161;  ! kg/kmol    
    pc%Tcrit = 582.82;             ! K    
    pc%pcrit = 3.3200E6;           ! Pa 
    pc%vcrit = 0.37800;            ! m^3/kmol 
    pc%CritCompressFactor = 0.259; ! -


    pc%MeltingPoint = 217.15_WP;              ! K
    pc%TripletPointT = 217.15_WP;             ! K
    pc%TripletPointp = 2.2202E00_WP;          ! Pa
    pc%NormalBoilingPoint = 396.65_WP;        ! K
    pc%liqMolarVolume = 0.123630_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.7760E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.2600E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.0920E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.3490E7_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.4920E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.3670E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -3.2020E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 3.0531E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.3801_WP;            ! -
    pc%RadiusOfGyration = 4.0090E-10_WP;      ! m
    pc%SolubilityParameter = 1.7930E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 9.5700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06973_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 9.8900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3980_WP;           ! -
    pc%FlashPoint = 287.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.0_WP, 8.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 288.0_WP, 327.0_WP /); ! K
    pc%AutoignitionTemp = 712.00_WP;          ! K


    pc%solDensityMol%A = 1.0200E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 217.50_WP, 1.0200E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.7666E-01_WP;
    pc%liqDensityMol%B = 2.5578E-01_WP;
    pc%liqDensityMol%C = 5.8282E02_WP;
    pc%liqDensityMol%D = 2.7746E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 217.50_WP, 582.82_WP /);
    pc%liqDensityMol%range(:,2) = (/ 8.7631E00_WP, 2.6455E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.3155E01_WP;
    pc%vapPressure%B = -7.2429E03_WP;
    pc%vapPressure%C = -7.2569E00_WP;
    pc%vapPressure%D =  1.2741E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 217.50_WP, 582.82_WP /);
    pc%vapPressure%range(:,2) = (/ 2.2202E00_WP, 3.3216E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.6232E07_WP;
    pc%HeatOfVapMol%B = 3.9972E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 217.50_WP, 582.82_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.6655E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -4.1000E03_WP;
    pc%solHeatCapacityMol%B = 4.2800E02_WP;
    pc%solHeatCapacityMol%C = 2.6860E01_WP;
    pc%solHeatCapacityMol%D = -3.5340E-01_WP;
    pc%solHeatCapacityMol%E = 1.3463E-03_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 12.52_WP, 140.13_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 4.8084E03_WP, 1.2999E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.3596E05_WP;
    pc%liqHeatCapacityMol%B = -3.4594E02_WP;
    pc%liqHeatCapacityMol%C = 9.4278E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 217.50_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.0532E05_WP, 2.7632E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.1237E05_WP;
    pc%IG_HeatCapacityMol%B = 2.9360E05_WP;
    pc%IG_HeatCapacityMol%C = 1.4010E03_WP;
    pc%IG_HeatCapacityMol%D = 1.6010E05_WP;
    pc%IG_HeatCapacityMol%E = 6.5050E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.1443E05_WP, 3.5874E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.4188E-01_WP;
    pc%SecondViralCoefMol%B = -2.7496E02_WP;
    pc%SecondViralCoefMol%C = -5.8174E07_WP;
    pc%SecondViralCoefMol%D = -2.1114E20_WP;
    pc%SecondViralCoefMol%E = 2.4535E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 297.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.0275E00_WP, 4.1329E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.3682E01_WP;
    pc%liqViscosity%B = 1.2834E03_WP;
    pc%liqViscosity%C = 3.3755E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 217.50_WP, 396.65_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.5631E-03_WP, 2.1848E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 9.8882E-08_WP;
    pc%vapViscosity%B = 7.7550E-01_WP;
    pc%vapViscosity%C = 9.9825E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 217.50_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.4032E-06_WP, 1.9067E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.3493E-01_WP;
    pc%liqThermalConductivity%B = -2.9120E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 217.50_WP, 466.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7159E-01_WP, 9.9231E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -3.3262E-01_WP;
    pc%vapThermalConductivity%B = 1.2054E-01_WP;
    pc%vapThermalConductivity%C = -2.4726E03_WP;
    pc%vapThermalConductivity%D = -5.4934E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.00_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 7.9990E-03_WP, 1.0980E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.7740E-02_WP;
    pc%SurfaceTension%B = 1.1689E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 217.50_WP, 582.82_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.3447E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_3hexanone


  subroutine pc_3pentanone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '3pentanone';
    pc%CAS = '96-22-0';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 86.134_WP;   ! kg/kmol    
    pc%Tcrit = 560.95_WP;             ! K    
    pc%pcrit = 3.7400E6_WP;           ! Pa 
    pc%vcrit = 0.33600_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.269_WP; ! -


    pc%MeltingPoint = 234.18_WP;              ! K
    pc%TripletPointT = 234.18_WP;             ! K
    pc%TripletPointp = 7.3422E01_WP;          ! Pa
    pc%NormalBoilingPoint = 375.14_WP;        ! K
    pc%liqMolarVolume = 0.106410_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.5790E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.3440E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.7000E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.1590E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.8804E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.4200E08_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -2.9651E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 2.6600E05_WP;         ! J/kmolK    


    pc%AcentricFactor = 0.3448_WP;            ! -
    pc%RadiusOfGyration = 3.5830E-10_WP;      ! m
    pc%SolubilityParameter = 1.8410E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 9.4100E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05950_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.5400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3900_WP;           ! -
    pc%FlashPoint = 286.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 8.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 277.0_WP, 308.0_WP /); ! K
    pc%AutoignitionTemp = 725.00_WP;          ! K


    pc%solDensityMol%A = 1.1800E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 234.18_WP, 1.1800E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 7.1811E-01_WP;
    pc%liqDensityMol%B = 2.4129E-01_WP;
    pc%liqDensityMol%C = 5.6095E02_WP;
    pc%liqDensityMol%D = 2.7996E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 234.18_WP, 560.95_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0102E01_WP, 2.9761E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A = 4.4286E01_WP;
    pc%vapPressure%B = -5.4151E03_WP;
    pc%vapPressure%C = -3.0913E00_WP;
    pc%vapPressure%D =  1.8580E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 234.18_WP, 560.95_WP /);
    pc%vapPressure%range(:,2) = (/ 7.3422E01_WP, 3.6993E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.2359E07_WP;
    pc%HeatOfVapMol%B = 4.0465E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 234.18_WP, 560.95_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.2075E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -9.2023E03_WP;
    pc%solHeatCapacityMol%B = 1.1108E03_WP;
    pc%solHeatCapacityMol%C = 5.6357E00_WP;
    pc%solHeatCapacityMol%D = -1.4318E-01_WP;
    pc%solHeatCapacityMol%E = 6.8742E-04_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 11.16_WP, 110.39_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.7078E03_WP, 9.1569E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.9302E05_WP;
    pc%liqHeatCapacityMol%B = -1.7643E02_WP;
    pc%liqHeatCapacityMol%C = 5.6690E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 234.18_WP, 375.14_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.8279E05_WP, 2.0661E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 9.6896E04_WP;
    pc%IG_HeatCapacityMol%B = 2.4907E05_WP;
    pc%IG_HeatCapacityMol%C = 1.4177E03_WP;
    pc%IG_HeatCapacityMol%D = 1.3010E05_WP;
    pc%IG_HeatCapacityMol%E = 6.4670E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.0536E05_WP, 3.0358E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.9680E-01_WP;
    pc%SecondViralCoefMol%B = -2.1873E02_WP;
    pc%SecondViralCoefMol%C = -4.5060E07_WP;
    pc%SecondViralCoefMol%D = -1.5690E20_WP;
    pc%SecondViralCoefMol%E = 1.9653E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 280.47_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.8930E00_WP, 3.7623E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.8695E00_WP;
    pc%liqViscosity%B = 5.9632E02_WP;
    pc%liqViscosity%C = -1.2025E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 234.18_WP, 375.14_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.0240E-03_WP, 2.2319E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.1640E-07_WP;
    pc%vapViscosity%B = 7.6150E-01_WP;
    pc%vapViscosity%C = 1.0794E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 234.18_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 5.0785E-06_WP, 2.0227E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1569E-01_WP;
    pc%liqThermalConductivity%B = -2.4081E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 234.18_WP, 375.14_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5930E-01_WP, 1.2535E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.2775E01_WP;
    pc%vapThermalConductivity%B = 1.0019E00_WP;
    pc%vapThermalConductivity%C = 1.9100E08_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.00_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 8.9821E-03_WP, 1.2082E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.0021E-02_WP;
    pc%SurfaceTension%B = 1.1704E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 234.18_WP, 560.95_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1888E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_3pentanone


  subroutine pc_3petanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = '3petanol';
    pc%CAS = '584-02-1';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 88.150_WP;  ! kg/kmol    
    pc%Tcrit = 559.60_WP;             ! K    
    pc%pcrit = 3.8800E6_WP;           ! Pa 
    pc%vcrit = 0.32500_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.271_WP; ! -


    pc%MeltingPoint = 204.15_WP;              ! K
    pc%TripletPointT = 204.15_WP;             ! K
    pc%TripletPointp = 1.4819E-02_WP;         ! Pa
    pc%NormalBoilingPoint = 388.45_WP;        ! K
    pc%liqMolarVolume = 0.108015_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.1673E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.5816E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.8200E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  9.0800E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.0483E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.7050E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -3.6890E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.4910E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.5383_WP;            ! -
    pc%RadiusOfGyration = 3.2530E-10_WP;      ! m
    pc%SolubilityParameter = 2.1780E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.4700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06262_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.9700E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4079_WP;           ! -
    pc%FlashPoint = 313.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 9.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 300.0_WP, 332.0_WP /); ! K
    pc%AutoignitionTemp = 708.15_WP;          ! K


    pc%solDensityMol%A = 1.1900E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 204.15_WP, 1.1900E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.0198E-01_WP;
    pc%liqDensityMol%B = 2.6069E-01_WP;
    pc%liqDensityMol%C = 5.5960E02_WP;
    pc%liqDensityMol%D = 2.6161E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 204.15_WP, 559.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0152E01_WP, 3.0764E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.3775E02_WP;
    pc%vapPressure%B = -1.1249E04_WP;
    pc%vapPressure%C = -1.6331E01_WP;
    pc%vapPressure%D =  2.8016E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 204.15_WP, 559.60_WP /);
    pc%vapPressure%range(:,2) = (/ 1.4819E-02_WP, 3.9083E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 8.4954E07_WP;
    pc%HeatOfVapMol%B = 6.0261E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 204.15_WP, 559.60_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.4626E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 0.0_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP  /);
    pc%solHeatCapacityMol%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP  /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = -1.6787E06_WP;
    pc%liqHeatCapacityMol%B = 1.1243E04_WP;
    pc%liqHeatCapacityMol%C = -1.6088E01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 303.02_WP, 326.16_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.5094E05_WP, 2.7687E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0751E05_WP;
    pc%IG_HeatCapacityMol%B = 3.0130E05_WP;
    pc%IG_HeatCapacityMol%C = 1.8100E03_WP;
    pc%IG_HeatCapacityMol%D = 2.1755E05_WP;
    pc%IG_HeatCapacityMol%E = 8.1170E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 300.00_WP, 1200.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3596E05_WP, 3.2093E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.2585E-01_WP;
    pc%SecondViralCoefMol%B = -2.3629E02_WP;
    pc%SecondViralCoefMol%C = -3.0124E07_WP;
    pc%SecondViralCoefMol%D = 3.5131E20_WP;
    pc%SecondViralCoefMol%E = -1.5281E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 279.80_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.1804E00_WP, 5.9407E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.5290E02_WP;
    pc%liqViscosity%B = 1.4984E04_WP;
    pc%liqViscosity%C = 3.4620E01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 250.00_WP, 388.45_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.6349E-01_WP, 3.6484E-04_WP /);
    pc%liqViscosity%eqn = 101;


    pc%vapViscosity%A = 1.1891E-07_WP;
    pc%vapViscosity%B = 7.6401E-01_WP;
    pc%vapViscosity%C = 1.0337E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%vapViscosity%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1707E-01_WP;
    pc%liqThermalConductivity%B = -2.7901E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 3.1596E03_WP;
    pc%vapThermalConductivity%B = 9.0260E-01_WP;
    pc%vapThermalConductivity%C = 1.2930E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.8608E-02_WP;
    pc%SurfaceTension%B = 1.1579E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%SurfaceTension%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_3petanol


  subroutine pc_acetone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'acetone';
    pc%CAS = '67-64-1';
    ! components
    pc%ChemicalFormula%C = 3;         
    pc%ChemicalFormula%H = 6;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 58.080_WP;   ! kg/kmol    
    pc%Tcrit = 508.20_WP;             ! K    
    pc%pcrit = 4.7015E6_WP;          ! Pa 
    pc%vcrit = 0.20900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.233_WP; ! -


    pc%MeltingPoint = 178.45_WP;              ! K
    pc%TripletPointT = 178.45_WP;             ! K
    pc%TripletPointp = 2.5938E00_WP;          ! Pa
    pc%NormalBoilingPoint = 329.44_WP;        ! K
    pc%liqMolarVolume = 0.073931_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.1715E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.5272E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  2.9535E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  5.6912E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.6592E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.3064_WP;            ! -
    pc%RadiusOfGyration = 2.7460E-10_WP;      ! m
    pc%SolubilityParameter = 1.9774E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 9.6066E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.03904_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 5.8400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3560_WP;           ! -
    pc%FlashPoint = 255.37_WP;                ! K
    pc%FlammabilityLimits = (/ 2.6_WP, 12.8_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 810.93_WP;          ! K


    pc%solDensityMol%A = 1.6677E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 174.15_WP, 1.6677E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.2298E00_WP;
    pc%liqDensityMol%B = 2.5760E-01_WP;
    pc%liqDensityMol%C = 5.0820E02_WP;
    pc%liqDensityMol%D = 2.9903E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 178.45_WP, 508.20_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.5721E01_WP, 4.7741E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.0720E01_WP;
    pc%vapPressure%B = -5.6850E03_WP;
    pc%vapPressure%C = -7.3510E00_WP;
    pc%vapPressure%D =  6.300E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 178.45_WP, 508.20_WP /);
    pc%vapPressure%range(:,2) = (/ 2.5938E00_WP, 4.6733E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.9170E07_WP;
    pc%HeatOfVapMol%B = 1.0360E00_WP;
    pc%HeatOfVapMol%C = -1.2940E00_WP;
    pc%HeatOfVapMol%D = 6.7200E-01_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 178.45_WP, 508.2_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.6887E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.2000E04_WP;
    pc%solHeatCapacityMol%B = 1.2760E03_WP;
    pc%solHeatCapacityMol%C = -6.1700E00_WP;
    pc%solHeatCapacityMol%D = 1.3200E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 30.00_WP, 178.45_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.1083E04_WP, 9.4233E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.3560E05_WP;
    pc%liqHeatCapacityMol%B = -1.7700E02_WP;
    pc%liqHeatCapacityMol%C = 2.8370E-01_WP;
    pc%liqHeatCapacityMol%D = 6.8900E-04_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 178.45_WP, 329.44_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.169E05_WP, 1.3271E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 4.8120E04_WP;
    pc%IG_HeatCapacityMol%B = 1.6440E05_WP;
    pc%IG_HeatCapacityMol%C = 1.2500E03_WP;
    pc%IG_HeatCapacityMol%D = 7.1700E04_WP;
    pc%IG_HeatCapacityMol%E = -5.244E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 100.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 4.8340E04_WP, 1.8699E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.1070E-01_WP;
    pc%SecondViralCoefMol%B = -9.9200E01_WP;
    pc%SecondViralCoefMol%C = -3.0540E07_WP;
    pc%SecondViralCoefMol%D = 9.8900E18_WP;
    pc%SecondViralCoefMol%E = -1.6454E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 254.10_WP, 2540.50_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.298E00_WP, 6.9790E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.4918E01_WP;
    pc%liqViscosity%B = 1.0234E03_WP;
    pc%liqViscosity%C = 5.9610E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 223.15_WP, 329.44_WP /);
    pc%liqViscosity%range(:,2) = (/ 8.1834E-04_WP, 2.3506E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.1005E-08_WP;
    pc%vapViscosity%B = 9.7620E-01_WP;
    pc%vapViscosity%C = 2.3139E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 300.00_WP, 650.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 7.5393E-06_WP, 1.6680E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.5020E-01_WP;
    pc%liqThermalConductivity%B = -2.9800E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 270.00_WP, 315.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6974E-01_WP, 1.5633E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -2.6800E01_WP;
    pc%vapThermalConductivity%B = 9.0980E-01_WP;
    pc%vapThermalConductivity%C = -1.2650E08_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.15_WP, 571.70_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 9.5298E-03_WP, 3.9056E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.2200E-02_WP;
    pc%SurfaceTension%B = 1.1240E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 178.45_WP, 508.20_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.8251E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_acetone


  subroutine pc_diethylether(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'diethylether';
    pc%CAS = '60-29-7';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 74.123_WP;   ! kg/kmol    
    pc%Tcrit = 466.70_WP;             ! K    
    pc%pcrit = 3.6400E6_WP;           ! Pa 
    pc%vcrit = 0.28000_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.263_WP; ! -


    pc%MeltingPoint = 156.85_WP;              ! K
    pc%TripletPointT = 156.85_WP;             ! K
    pc%TripletPointp = 3.9545E-01_WP;         ! Pa
    pc%NormalBoilingPoint = 307.58_WP;        ! K
    pc%liqMolarVolume = 0.104692_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.5210E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.2210E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.4230E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  7.1900E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.2032E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.2811_WP;            ! -
    pc%RadiusOfGyration = 3.1770E-10_WP;      ! m
    pc%SolubilityParameter = 1.5420E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 3.8400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05150_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.5400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3495_WP;           ! -
    pc%FlashPoint = 228.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.9_WP, 48.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp =433.15_WP;           ! K


    pc%solDensityMol%A = 1.3400E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 156.85_WP, 1.3400E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.5540E-01_WP;
    pc%liqDensityMol%B = 2.6847E-01_WP;
    pc%liqDensityMol%C = 4.6670E02_WP;
    pc%liqDensityMol%D = 2.8140E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 156.85_WP, 466.70_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.1487E01_WP, 3.5587E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A = 1.3690E02_WP;
    pc%vapPressure%B = -6.9543E03_WP;
    pc%vapPressure%C = -1.9254E01_WP;
    pc%vapPressure%D =  2.4508E-02_WP;
    pc%vapPressure%E =  1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 156.85_WP, 466.70_WP /);
    pc%vapPressure%range(:,2) = (/ 3.9545E-01_WP, 3.6412E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.0600E07_WP;
    pc%HeatOfVapMol%B = 3.8680E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 156.85_WP, 466.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.4651E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -2.3760E04_WP;
    pc%solHeatCapacityMol%B = 1.8611E03_WP;
    pc%solHeatCapacityMol%C = -1.2638E01_WP;
    pc%solHeatCapacityMol%D = 3.7650E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 30.00_WP, 156.85_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.1715E04_WP, 1.0252E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 4.4400E04_WP;
    pc%liqHeatCapacityMol%B = 1.3010E03_WP;
    pc%liqHeatCapacityMol%C = -5.5000E00_WP;
    pc%liqHeatCapacityMol%D = 8.7630E-03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 156.92_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.4698E05_WP, 3.3202E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 8.6210E04_WP;
    pc%IG_HeatCapacityMol%B = 2.5510E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5413E03_WP;
    pc%IG_HeatCapacityMol%D = 1.4370E05_WP;
    pc%IG_HeatCapacityMol%E = -6.8890E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 9.3157E04_WP, 2.9244E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.6330E-01_WP;
    pc%SecondViralCoefMol%B = -1.3680E02_WP;
    pc%SecondViralCoefMol%C = -2.0080E07_WP;
    pc%SecondViralCoefMol%D = -1.7500E19_WP;
    pc%SecondViralCoefMol%E = 2.6570E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 233.35_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.6987E00_WP, 6.6150E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = 1.0197E01_WP;
    pc%liqViscosity%B = -6.3800E01_WP;
    pc%liqViscosity%C = -3.2260E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 200.00_WP, 373.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 7.3593E-00_WP, 1.1412E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.9480E-06_WP;
    pc%vapViscosity%B = 4.1000E-01_WP;
    pc%vapViscosity%C = 4.9580E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 156.85_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.7200E-06_WP, 2.2116E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.4950E-01_WP;
    pc%liqThermalConductivity%B = -4.0700E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 156.85_WP, 433.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.8566E-01_WP, 7.3208E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -4.4894E-03_WP;
    pc%vapThermalConductivity%B = 6.1550E-01_WP;
    pc%vapThermalConductivity%C = -3.2663E03_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 200.00_WP, 600.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 7.6364E-03_WP, 5.1806E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.7928E-02_WP;
    pc%SurfaceTension%B = 1.2376E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 156.85_WP, 466.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.4893E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_diethylether

  subroutine pc_oxymethylenether4(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'OME4';
    pc%CAS = '';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 5;

    pc%MolecularWeight = 166.173_WP;   ! kg/kmol    
    pc%Tcrit = 615.0_WP;             ! K    
    pc%pcrit = 1.48691E06_WP;          ! Pa XXXX NOT SURE
    pc%vcrit = 0.21300_WP;            ! m^3/kmol XXXX NOT SURE
    pc%CritCompressFactor = 0.211_WP; ! -


    pc%MeltingPoint = 268.0_WP;              ! K FROM EXCEL SHEET
    pc%TripletPointT = 168.35_WP;             ! K DO NOT KNOW
    pc%TripletPointp = 1.9634E00_WP;         ! Pa DO NOT KNOW
    pc%NormalBoilingPoint = 456.66_WP;        ! K FROM EXCEL SHEET
    pc%liqMolarVolume = 0.089127_WP;          ! m^3/kmol DO NOT KNOW

    pc%IG_HeatOfFormation = -3.4820E08_WP;             ! J/kmol DO NOT KNOW
    pc%IG_GibbsOfFormation = -2.2631E08_WP;            ! J/kmol DO NOT KNOW
    pc%IG_AbsoluteEntropy =  3.3572E05_WP;             ! J/kmolK DO NOT KNOW
    pc%HeatFusionAtMeltPoint =  8.3303E06_WP;          ! J/kmol DO NOT KNOW
    pc%StandardNetHeatOfCombustion = -1.7998E09;    ! J/kmol DO NOT KNOW
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol DO NOT KNOW
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol DO NOT KNOW
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK DO NOT KNOW


    pc%AcentricFactor = 0.2904_WP;            ! - DO NOT KNOW
    pc%RadiusOfGyration = 2.8910E-10_WP;      ! m DO NOT KNOW
    pc%SolubilityParameter = 1.7189E04_WP;    ! (J/m^3)^.5 DO NOT KNOW
    pc%DipoleMoment = 2.4684E-30_WP;          ! C^m DO NOT KNOW
    pc%VanDerWaalsVolume = 0.04497_WP;        ! m^3/kmol DO NOT KNOW


    pc%VanDerWaalsArea = 6.7900E08_WP;        ! m^2/kmol DO NOT KNOW
    pc%RefractiveIndex = 1.3504_WP;           ! - DO NOT KNOW
    pc%FlashPoint = 255.37_WP;                ! K DO NOT KNOW
    pc%FlammabilityLimits = (/ 1.6_WP, 17.6_WP /);    ! vol-! DO NOT KNOW
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K DO NOT KNOW
    pc%AutoignitionTemp = 510.35_WP;          ! K DO NOT KNOW

    ! XXXX NOT KNOWN
    pc%solDensityMol%A = 1.5600E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 168.35_WP, 1.5600E01_WP /);
    pc%solDensityMol%eqn = 100;

    ! XXXX FROM EXCEL SHEET
    pc%liqDensityMol%A = 1.38435695329566E+01_WP;
    pc%liqDensityMol%B = -6.64940622122727E-02_WP;
    pc%liqDensityMol%C = 2.34788283896903E-04_WP;
    pc%liqDensityMol%D = -3.85301058535382E-07_WP;
    pc%liqDensityMol%E = 2.20853685616797E-10_WP;
    pc%liqDensityMol%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.652820567_WP, 3.665080563_WP /);
    pc%liqDensityMol%eqn = 100;

    ! XXXX FROM EXCEL SHEET
    pc%vapPressure%A = 73.4273462767746_WP;
    pc%vapPressure%B = -8497.77116527637_WP;
    pc%vapPressure%C = -7.09835749350697_WP;
    pc%vapPressure%D = 0.00000000000000001_WP;
    pc%vapPressure%E = 6.10934905086933_WP;
    pc%vapPressure%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%vapPressure%range(:,2) = (/ 7.68484E+00_WP, 3.72994E+06_WP /);
    pc%vapPressure%eqn = 101;  

    ! XXXX FROM EXCEL SHEET
    pc%HeatOfVapMol%A = 8.42909520204301E+07_WP;
    pc%HeatOfVapMol%B = 1.200766275140770E+00_WP;
    pc%HeatOfVapMol%C = -1.047590592324820E+00_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.505535408532550E+07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    ! XXXX NOT AVAILABLE
    pc%solHeatCapacityMol%A = -1.6050E04_WP;
    pc%solHeatCapacityMol%B = 1.3984E03_WP;
    pc%solHeatCapacityMol%C = -5.8400E00_WP;
    pc%solHeatCapacityMol%D = 8.5900E-03_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 168.35_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.6410E03_WP, 9.4841E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    ! FROM EXCEL SHEET
    pc%liqHeatCapacityMol%A = 3.04940766054870E+05_WP;
    pc%liqHeatCapacityMol%B = -2.74831030094651E+02_WP;
    pc%liqHeatCapacityMol%C = 5.02855394544221E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.67403135847247E+05_WP, 3.26318553839206E+05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    ! FROM EXCEL SHEET
    pc%IG_HeatCapacityMol%A = 1.71002033600000E+05_WP;
    pc%IG_HeatCapacityMol%B = 2.57840000000000E+02_WP;
    pc%IG_HeatCapacityMol%C = 0.0_WP;
    pc%IG_HeatCapacityMol%D = 0.0_WP;
    pc%IG_HeatCapacityMol%E = 0.0_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 268.0_WP, 620.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.40103153600000E+05_WP, 3.30862833600000E+05_WP /);
    pc%IG_HeatCapacityMol%eqn = 100;

    ! XXXX NOT AVAILABLE
    pc%SecondViralCoefMol%A = 1.6669E-01_WP;
    pc%SecondViralCoefMol%B = -1.6805E02_WP;
    pc%SecondViralCoefMol%C = -1.9195E07_WP;
    pc%SecondViralCoefMol%D = 2.9849E19_WP;
    pc%SecondViralCoefMol%E = -8.7864E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 240.30_WP, 2400.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.5200_WP, 9.5290E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    ! XXXX FROM EXCEL SHEET
    pc%liqViscosity%A = -4.73279672E+01_WP;
    pc%liqViscosity%B = 3.27503801E+03_WP;
    pc%liqViscosity%C = 5.28039807E+00_WP;
    pc%liqViscosity%D = 3.28403307E-06_WP;
    pc%liqViscosity%E = 4.86541006E-01_WP;
    pc%liqViscosity%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.7537908E-03_WP, 3.0538220E-04_WP /);
    pc%liqViscosity%eqn = 101;

    ! XXXX FROM EXCEL SHEET
    pc%vapViscosity%A = 1.72298102E-08_WP;
    pc%vapViscosity%B = 1.23762536E+00_WP;
    pc%vapViscosity%C = 1.80430750E+03_WP;
    pc%vapViscosity%D = -4.83910313E+04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 268.0_WP, 620.0_WP /);
    pc%vapViscosity%range(:,2) = (/ 2.4698702112E-06_WP, 1.3008621821E-05_WP /);
    pc%vapViscosity%eqn = 102;

    ! XXXX FROM EXCEL SHEET
    pc%liqThermalConductivity%A = 6.81596746E-01_WP;
    pc%liqThermalConductivity%B = -2.66777745E-03_WP;
    pc%liqThermalConductivity%C = 5.98408357E-06_WP;
    pc%liqThermalConductivity%D = -6.79546584E-09_WP;
    pc%liqThermalConductivity%E = 3.06374344E-12_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.8143332E-01_WP, 1.6174228E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    ! XXXX FROM EXCEL SHEET assuming Pr=0.7
    pc%vapThermalConductivity%A = 9.7323914200E-02_WP;
    pc%vapThermalConductivity%B = 1.22192393E-01_WP;
    pc%vapThermalConductivity%C = -2.30249695E+03_WP;
    pc%vapThermalConductivity%D = 3.26150483E+06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.15_WP, 620.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 5.0958598408E-03_WP, 3.6998199098E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    ! XXXX FROM EXCEL SHEET
    pc%SurfaceTension%A = 4.28E-02_WP;
    pc%SurfaceTension%B = 0.858893924_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 293.35_WP, 615.6_WP /);
    pc%SurfaceTension%range(:,2) = (/ 0.0246_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106; 

  end subroutine pc_oxymethylenether4

  subroutine pc_oxymethylenether4_pcpsaft(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'OME4';
    pc%CAS = '';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 5;

    pc%MolecularWeight = 166.172_WP;   ! kg/kmol    
    pc%Tcrit = 616.581_WP;             ! K    
    pc%pcrit = 2.95007_WP;             ! Pa
    pc%vcrit = 0.533478_WP;            ! m^3/kmol
    pc%CritCompressFactor = 0.211_WP;  ! -


    pc%MeltingPoint = 256.499_WP;              ! K 
    pc%TripletPointT = MAX_REAL_WP;            ! K DO NOT KNOW
    pc%TripletPointp = MAX_REAL_WP;            ! Pa DO NOT KNOW
    pc%NormalBoilingPoint = 424.45_WP;         ! K
    pc%liqMolarVolume = 0.201301_WP;           ! m^3/kmol

    pc%IG_HeatOfFormation = MAX_REAL_WP;             ! J/kmol DO NOT KNOW
    pc%IG_GibbsOfFormation = MAX_REAL_WP;            ! J/kmol DO NOT KNOW
    pc%IG_AbsoluteEntropy = MAX_REAL_WP;             ! J/kmolK DO NOT KNOW
    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;          ! J/kmol DO NOT KNOW
    pc%StandardNetHeatOfCombustion = MAX_REAL_WP;    ! J/kmol DO NOT KNOW
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol DO NOT KNOW
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol DO NOT KNOW
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK DO NOT KNOW


    pc%AcentricFactor = 0.3836_WP;            ! - PCP-SAFT calculations
    pc%RadiusOfGyration = MAX_REAL_WP;        ! m DO NOT KNOW
    pc%SolubilityParameter = MAX_REAL_WP;     ! (J/m^3)^.5 DO NOT KNOW
    pc%DipoleMoment = MAX_REAL_WP;            ! C^m DO NOT KNOW
    pc%VanDerWaalsVolume = MAX_REAL_WP;       ! m^3/kmol DO NOT KNOW


    pc%VanDerWaalsArea = MAX_REAL_WP;         ! m^2/kmol DO NOT KNOW
    pc%RefractiveIndex = MAX_REAL_WP;         ! - DO NOT KNOW
    pc%FlashPoint = MAX_REAL_WP;              ! K DO NOT KNOW
    pc%FlammabilityLimits = (/ MAX_REAL_WP, MAX_REAL_WP /);    ! vol-! DO NOT KNOW
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K DO NOT KNOW
    pc%AutoignitionTemp = MAX_REAL_WP;          ! K DO NOT KNOW

    ! XXXX NOT KNOWN
    pc%solDensityMol%A = 0.0_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solDensityMol%eqn = 100;

    ! XXXX
    pc%liqDensityMol%A = 1.38435695329566E+01_WP;
    pc%liqDensityMol%B = -6.64940622122727E-02_WP;
    pc%liqDensityMol%C = 2.34788283896903E-04_WP;
    pc%liqDensityMol%D = -3.85301058535382E-07_WP;
    pc%liqDensityMol%E = 2.20853685616797E-10_WP;
    pc%liqDensityMol%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.652820567_WP, 3.665080563_WP /);
    pc%liqDensityMol%eqn = 100;

    ! XXXX FROM EXCEL SHEET
    pc%vapPressure%A = 73.4273462767746_WP;
    pc%vapPressure%B = -8497.77116527637_WP;
    pc%vapPressure%C = -7.09835749350697_WP;
    pc%vapPressure%D = 0.00000000000000001_WP;
    pc%vapPressure%E = 6.10934905086933_WP;
    pc%vapPressure%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%vapPressure%range(:,2) = (/ 7.68484E+00_WP, 3.72994E+06_WP /);
    pc%vapPressure%eqn = 101;  

    ! XXXX FROM EXCEL SHEET
    pc%HeatOfVapMol%A = 8.42909520204301E+07_WP;
    pc%HeatOfVapMol%B = 1.200766275140770E+00_WP;
    pc%HeatOfVapMol%C = -1.047590592324820E+00_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.505535408532550E+07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    ! XXXX NOT AVAILABLE
    pc%solHeatCapacityMol%A = -1.6050E04_WP;
    pc%solHeatCapacityMol%B = 1.3984E03_WP;
    pc%solHeatCapacityMol%C = -5.8400E00_WP;
    pc%solHeatCapacityMol%D = 8.5900E-03_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 168.35_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.6410E03_WP, 9.4841E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    ! FROM EXCEL SHEET
    pc%liqHeatCapacityMol%A = 3.04940766054870E+05_WP;
    pc%liqHeatCapacityMol%B = -2.74831030094651E+02_WP;
    pc%liqHeatCapacityMol%C = 5.02855394544221E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.67403135847247E+05_WP, 3.26318553839206E+05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    ! FROM EXCEL SHEET
    pc%IG_HeatCapacityMol%A = 1.71002033600000E+05_WP;
    pc%IG_HeatCapacityMol%B = 2.57840000000000E+02_WP;
    pc%IG_HeatCapacityMol%C = 0.0_WP;
    pc%IG_HeatCapacityMol%D = 0.0_WP;
    pc%IG_HeatCapacityMol%E = 0.0_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 268.0_WP, 620.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.40103153600000E+05_WP, 3.30862833600000E+05_WP /);
    pc%IG_HeatCapacityMol%eqn = 100;

    ! XXXX NOT AVAILABLE
    pc%SecondViralCoefMol%A = 1.6669E-01_WP;
    pc%SecondViralCoefMol%B = -1.6805E02_WP;
    pc%SecondViralCoefMol%C = -1.9195E07_WP;
    pc%SecondViralCoefMol%D = 2.9849E19_WP;
    pc%SecondViralCoefMol%E = -8.7864E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 240.30_WP, 2400.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.5200_WP, 9.5290E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    ! XXXX FROM EXCEL SHEET
    pc%liqViscosity%A = -4.73279672E+01_WP;
    pc%liqViscosity%B = 3.27503801E+03_WP;
    pc%liqViscosity%C = 5.28039807E+00_WP;
    pc%liqViscosity%D = 3.28403307E-06_WP;
    pc%liqViscosity%E = 4.86541006E-01_WP;
    pc%liqViscosity%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.7537908E-03_WP, 3.0538220E-04_WP /);
    pc%liqViscosity%eqn = 101;

    ! XXXX FROM EXCEL SHEET
    pc%vapViscosity%A = 1.72298102E-08_WP;
    pc%vapViscosity%B = 1.23762536E+00_WP;
    pc%vapViscosity%C = 1.80430750E+03_WP;
    pc%vapViscosity%D = -4.83910313E+04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 268.0_WP, 620.0_WP /);
    pc%vapViscosity%range(:,2) = (/ 2.4698702112E-06_WP, 1.3008621821E-05_WP /);
    pc%vapViscosity%eqn = 102;

    ! XXXX FROM EXCEL SHEET
    pc%liqThermalConductivity%A = 6.81596746E-01_WP;
    pc%liqThermalConductivity%B = -2.66777745E-03_WP;
    pc%liqThermalConductivity%C = 5.98408357E-06_WP;
    pc%liqThermalConductivity%D = -6.79546584E-09_WP;
    pc%liqThermalConductivity%E = 3.06374344E-12_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 268.0_WP, 615.6_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.8143332E-01_WP, 1.6174228E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    ! XXXX FROM EXCEL SHEET assuming Pr=0.7
    pc%vapThermalConductivity%A = 9.7323914200E-02_WP;
    pc%vapThermalConductivity%B = 1.22192393E-01_WP;
    pc%vapThermalConductivity%C = -2.30249695E+03_WP;
    pc%vapThermalConductivity%D = 3.26150483E+06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.15_WP, 620.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 5.0958598408E-03_WP, 3.6998199098E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    ! XXXX FROM EXCEL SHEET
    pc%SurfaceTension%A = 4.28E-02_WP;
    pc%SurfaceTension%B = 0.858893924_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 293.35_WP, 615.6_WP /);
    pc%SurfaceTension%range(:,2) = (/ 0.0246_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106; 

  end subroutine pc_oxymethylenether4_pcpsaft

  subroutine pc_dimethoxymethane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'OME1';
    pc%CAS = '	109-87-5';
    ! components
    pc%ChemicalFormula%C = 3;         
    pc%ChemicalFormula%H = 8;
    pc%ChemicalFormula%O = 2;

    pc%MolecularWeight = 76.095_WP;   ! kg/kmol    
    pc%Tcrit = 480.6_WP;             ! K    
    pc%pcrit = 3.9517E06_WP;          ! Pa 
    pc%vcrit = 0.21300_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.211_WP; ! -


    pc%MeltingPoint = 168.35_WP;              ! K
    pc%TripletPointT = 168.35_WP;             ! K
    pc%TripletPointp = 1.9634E00_WP;         ! Pa
    pc%NormalBoilingPoint = 315.00_WP;        ! K
    pc%liqMolarVolume = 0.089127_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.4820E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -2.2631E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.3572E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  8.3303E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.7998E09;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.2904_WP;            ! -
    pc%RadiusOfGyration = 2.8910E-10_WP;      ! m
    pc%SolubilityParameter = 1.7189E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 2.4684E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.04497_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 6.7900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3504_WP;           ! -
    pc%FlashPoint = 255.37_WP;                ! K
    pc%FlammabilityLimits = (/ 1.6_WP, 17.6_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 510.35_WP;          ! K


    pc%solDensityMol%A = 1.5600E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 168.35_WP, 1.5600E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.4355E00_WP;
    pc%liqDensityMol%B = 3.0576E-01_WP;
    pc%liqDensityMol%C = 4.8060E02_WP;
    pc%liqDensityMol%D = 3.17755E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 168.35_WP, 480.6_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.3194E01_WP, 4.6949E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.6361E01_WP;
    pc%vapPressure%B = -5.0483E03_WP;
    pc%vapPressure%C = -5.0136E00_WP;
    pc%vapPressure%D =  2.3659E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 168.35_WP, 480.6_WP /);
    pc%vapPressure%range(:,2) = (/ 1.9634E00_WP, 3.9500E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.1300E07_WP;
    pc%HeatOfVapMol%B = 1.7170E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 168.35_WP, 480.6_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.5183E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.6050E04_WP;
    pc%solHeatCapacityMol%B = 1.3984E03_WP;
    pc%solHeatCapacityMol%C = -5.8400E00_WP;
    pc%solHeatCapacityMol%D = 8.5900E-03_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.00_WP, 168.35_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.6410E03_WP, 9.4841E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.6050E05_WP;
    pc%liqHeatCapacityMol%B = -2.2930E02_WP;
    pc%liqHeatCapacityMol%C = 7.86600E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    !pc%liqHeatCapacityMol%range(:,1) = (/ 168.35_WP, 310.00_WP /);
    !pc%liqHeatCapacityMol%range(:,2) = (/ 1.4417E05_WP, 1.6495E05_WP /);
    pc%liqHeatCapacityMol%range(:,1) = (/ 168.35_WP, 480.60_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.4417E05_WP, 231984.424776_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 7.4976E04_WP;
    pc%IG_HeatCapacityMol%B = 1.6166E05_WP;
    pc%IG_HeatCapacityMol%C = 8.6287E02_WP;
    pc%IG_HeatCapacityMol%D = 7.8964E05_WP;
    pc%IG_HeatCapacityMol%E = 4.6718E03_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 298.15_WP, 1000.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 9.1671E04_WP, 2.0789E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.6669E-01_WP;
    pc%SecondViralCoefMol%B = -1.6805E02_WP;
    pc%SecondViralCoefMol%C = -1.9195E07_WP;
    pc%SecondViralCoefMol%D = 2.9849E19_WP;
    pc%SecondViralCoefMol%E = -8.7864E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 240.30_WP, 2400.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.5200_WP, 9.5290E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -6.0611E00_WP;
    pc%liqViscosity%B = 6.3620E02_WP;
    pc%liqViscosity%C = -7.1926E-1_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    !pc%liqViscosity%range(:,1) = (/ 168.35_WP, 378.50_WP /);
    !pc%liqViscosity%range(:,2) = (/ 2.5567E-03_WP, 1.7514E-04_WP /);
    pc%liqViscosity%range(:,1) = (/ 168.35_WP, 480.60_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.5567E-03_WP, 0.00010320492920665706_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 4.3034E-07_WP;
    pc%vapViscosity%B = 6.1720E-01_WP;
    pc%vapViscosity%C = 2.2554E+02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 315.0_WP, 975.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 8.7348E-06_WP, 2.4449E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.5239E-01_WP;
    pc%liqThermalConductivity%B = -3.6877E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    !pc%liqThermalConductivity%range(:,1) = (/ 168.35_WP, 315.0_WP /);
    !pc%liqThermalConductivity%range(:,2) = (/ 1.9031E-01_WP, 1.3623E-01_WP /);
    pc%liqThermalConductivity%range(:,1) = (/ 168.35_WP, 480.6_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.9031E-01_WP, 0.07515913800000001_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 1.1366E-03_WP;
    pc%vapThermalConductivity%B = 6.9110E-01_WP;
    pc%vapThermalConductivity%C = 4.1879E+02_WP;
    pc%vapThermalConductivity%D = 1.8761E+05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 315.00_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.4350E-02_WP, 8.3763E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 7.0967E-02_WP;
    pc%SurfaceTension%B = 1.2676_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 168.35_WP, 480.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 4.1083E-02_WP, 0.00E00_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_dimethoxymethane


  subroutine pc_dimethylether(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'dimethylether';
    pc%CAS = '115-10-6';
    ! components
    pc%ChemicalFormula%C = 2;         
    pc%ChemicalFormula%H = 6;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 46.069_WP;   ! kg/kmol    
    pc%Tcrit = 400.10_WP;             ! K    
    pc%pcrit = 5.3700E6_WP;           ! Pa 
    pc%vcrit = 0.17000_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.274_WP; ! -


    pc%MeltingPoint = 131.66_WP;              ! K
    pc%TripletPointT = 131.65_WP;             ! K
    pc%TripletPointp = 3.0496E00_WP;          ! Pa
    pc%NormalBoilingPoint = 248.31_WP;        ! K
    pc%liqMolarVolume = 0.070231_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -1.8410E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.1280E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  2.6670E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  4.9370E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.3284E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.2002_WP;            ! -
    pc%RadiusOfGyration = 2.1540E-10_WP;      ! m
    pc%SolubilityParameter = 1.5120E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 4.3400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.03104_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 4.8400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.2984_WP;           ! -
    pc%FlashPoint = 232.00_WP;                ! K
    pc%FlammabilityLimits = (/ 3.3_WP, 27.3_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp =623.15_WP;           ! K


    pc%solDensityMol%A = 2.2200E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 131.65_WP, 2.2200E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.5693E00_WP;
    pc%liqDensityMol%B = 2.6790E-01_WP;
    pc%liqDensityMol%C = 4.0010E02_WP;
    pc%liqDensityMol%D = 2.8820E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.8950E01_WP, 5.8578E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A = 4.4704E01_WP;
    pc%vapPressure%B = -3.5256E03_WP;
    pc%vapPressure%C = -3.4444E00_WP;
    pc%vapPressure%D =  5.4574E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%vapPressure%range(:,2) = (/ 3.0496E00_WP, 5.2735E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 2.9940E07_WP;
    pc%HeatOfVapMol%B = 3.5050E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 2.6032E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.3880E04_WP;
    pc%solHeatCapacityMol%B = 1.2580E03_WP;
    pc%solHeatCapacityMol%C = -8.2200E00_WP;
    pc%solHeatCapacityMol%D = 2.6752E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 24.26_WP, 131.64_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.2183E04_WP, 7.0305E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.1010E05_WP;
    pc%liqHeatCapacityMol%B = -1.5747E02_WP;
    pc%liqHeatCapacityMol%C = 5.1853E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    !pc%liqHeatCapacityMol%range(:,1) = (/ 131.65_WP, 250.00_WP /);
    !pc%liqHeatCapacityMol%range(:,2) = (/ 9.8356E04_WP, 1.0314E05_WP /);
    pc%liqHeatCapacityMol%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 9.8356E04_WP, 130102.54058530001_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 5.1480E04_WP;
    pc%IG_HeatCapacityMol%B = 1.4420E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6034E03_WP;
    pc%IG_HeatCapacityMol%D = 7.7470E04_WP;
    pc%IG_HeatCapacityMol%E = 7.2540E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 5.4363E04_WP, 1.6581E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.1570E-01_WP;
    pc%SecondViralCoefMol%B = -1.1010E02_WP;
    pc%SecondViralCoefMol%C = -5.2550E06_WP;
    pc%SecondViralCoefMol%D = -4.6140E18_WP;
    pc%SecondViralCoefMol%E = 7.6350E20_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.4028E00_WP, 4.0743E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.0620E01_WP;
    pc%liqViscosity%B = 4.4899E02_WP;
    pc%liqViscosity%C = 8.3967E-05_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    !pc%liqViscosity%range(:,1) = (/ 131.65_WP, 248.31_WP /);
    !pc%liqViscosity%range(:,2) = (/ 7.3982E-04_WP, 1.4903E-04_WP /);
    pc%liqViscosity%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%liqViscosity%range(:,2) = (/ 7.3982E-04_WP, 7.505401960287075E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.6800E-06_WP;
    pc%vapViscosity%B = 3.9750E-01_WP;
    pc%vapViscosity%C = 5.3400E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 131.65_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.6879E-06_WP, 2.7215E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 3.1174E-01_WP;
    pc%liqThermalConductivity%B = -5.6380E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    !pc%liqThermalConductivity%range(:,1) = (/ 131.65_WP, 320.03_WP /);
    !pc%liqThermalConductivity%range(:,2) = (/ 2.3752E-01_WP, 1.3131E-01_WP /);
    pc%liqThermalConductivity%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.3752E-01_WP, 0.08616362_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.9975E-02_WP;
    pc%vapThermalConductivity%B = 2.6670E-01_WP;
    pc%vapThermalConductivity%C = 1.0186E03_WP;
    pc%vapThermalConductivity%D = 1.0988E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 248.31_WP, 1500.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.1388E-02_WP, 1.9458E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.0960E-02_WP;
    pc%SurfaceTension%B = 1.2286E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 131.65_WP, 400.10_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.7336E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_dimethylether


  subroutine pc_dinbutylether(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'dinbutylether';
    pc%CAS = '142-96-1';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 130.23_WP;   ! kg/kmol    
    pc%Tcrit = 584.1_WP;              ! K    
    pc%pcrit = 2.4600E6_WP;           ! Pa 
    pc%vcrit = 0.487_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.247_WP; ! -


    pc%MeltingPoint = 175.3_WP;               ! K
    pc%TripletPointT = 175.3_WP;              ! K
    pc%TripletPointp = 7.1413E-04_WP;         ! Pa
    pc%NormalBoilingPoint = 414.15_WP;        ! K
    pc%liqMolarVolume = 0.170405_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.3340E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -8.8270E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  5.0140E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.6300E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.9469E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -9.9490E07_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -3.7790E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 3.8980E05_WP;         ! J/kmolK    


    pc%AcentricFactor = 0.447646_WP;           ! -
    pc%RadiusOfGyration = 4.7970E-10_WP;      ! m
    pc%SolubilityParameter = 1.5830E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 3.9000E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.09192_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.2880E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3968_WP;           ! -
    pc%FlashPoint = 298.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 7.6_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 308.0_WP, 341.0_WP /); ! K
    pc%AutoignitionTemp =467.1_WP;           ! K


    pc%solDensityMol%A = 7.7300E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 175.30_WP, 7.730E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.5941E-01_WP;
    pc%liqDensityMol%B = 2.7243E-01_WP;
    pc%liqDensityMol%C = 5.8410E02_WP;
    pc%liqDensityMol%D = 2.9932E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 175.30_WP, 584.10_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.6071E00_WP, 2.0534E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A = 7.2227E01_WP;
    pc%vapPressure%B = -7.5376E03_WP;
    pc%vapPressure%C = -7.0596E00_WP;
    pc%vapPressure%D =  9.1442E-18_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 175.30_WP, 584.10_WP /);
    pc%vapPressure%range(:,2) = (/ 7.1414E-04_WP, 2.4592E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.9616E07_WP;
    pc%HeatOfVapMol%B = 3.8833E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 175.30_WP, 584.10_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.1902E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 0.0_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.7072E05_WP;
    pc%liqHeatCapacityMol%B = -2.5983E02_WP;
    pc%liqHeatCapacityMol%C = 9.5427E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 175.30_WP, 450.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.5450E05_WP, 3.4704E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.6122E05_WP;
    pc%IG_HeatCapacityMol%B = 4.4777E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6831E03_WP;
    pc%IG_HeatCapacityMol%D = 2.9180E05_WP;
    pc%IG_HeatCapacityMol%E = 7.8160E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.6841E05_WP, 5.2145E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.3911E-01_WP;
    pc%SecondViralCoefMol%B = -3.7344E02_WP;
    pc%SecondViralCoefMol%C = -7.7353E06_WP;
    pc%SecondViralCoefMol%D = -1.3080E20_WP;
    pc%SecondViralCoefMol%E = 4.4188E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 292.05_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.2304E00_WP, 6.7226E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = 1.0027E01_WP;
    pc%liqViscosity%B = 2.0600E02_WP;
    pc%liqViscosity%C = -3.1607E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 175.30_WP, 414.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 5.9306E-03_WP, 1.9890E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 7.7147E-08_WP;
    pc%vapViscosity%B = 7.9906E-01_WP;
    pc%vapViscosity%C = 8.0765E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 175.30_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.2784E-06_WP, 1.7814E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.9418E-01_WP;
    pc%liqThermalConductivity%B = -2.2246E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 175.30_WP, 523.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5518E-01_WP, 7.7800E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 3.2694E-03_WP;
    pc%vapThermalConductivity%B = 5.8633E-01_WP;
    pc%vapThermalConductivity%C = 1.2599E03_WP;
    pc%vapThermalConductivity%D = 3.0089E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 323.15_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.2440E-02_WP, 7.3297E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1346E-02_WP;
    pc%SurfaceTension%B = 1.1604E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 175.30_WP, 584.10_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.3937E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_dinbutylether


  subroutine pc_dinpropylether(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'dinpropylether';
    pc%CAS = '1 11-43-3';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 102.177_WP;  ! kg/kmol    
    pc%Tcrit = 530.6_WP;              ! K    
    pc%pcrit = 3.0280E6_WP;           ! Pa 
    pc%vcrit = 0.382_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.262_WP; ! -


    pc%MeltingPoint = 149.95_WP;              ! K
    pc%TripletPointT = 149.95_WP;             ! K
    pc%TripletPointp = 7.6265E-04_WP;         ! Pa
    pc%NormalBoilingPoint = 363.23_WP;        ! K
    pc%liqMolarVolume = 0.137649_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.9290E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.0550E08_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.2250E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.0770E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.7250E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.1204E07_WP;       ! J/kmol
    pc%StandardHeatOfFormation = -3.2882E08_WP;        ! J/kmol
    pc%StandardAbsoluteEntropy = 3.2388E05_WP;         ! J/kmolK    


    pc%AcentricFactor = 0.368756_WP;          ! -
    pc%RadiusOfGyration = 3.9140E-10_WP;      ! m
    pc%SolubilityParameter = 1.5550E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 4.0400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.07196_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.0220E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.378_WP;            ! -
    pc%FlashPoint = 294.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 9.5_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 264.0_WP, 301.0_WP /); ! K
    pc%AutoignitionTemp =488.15_WP;           ! K


    pc%solDensityMol%A = 9.8800E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 149.95_WP, 9.8800E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 7.2390E-01_WP;
    pc%liqDensityMol%B = 2.7653E-01_WP;
    pc%liqDensityMol%C = 5.3060E02_WP;
    pc%liqDensityMol%D = 2.7940E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 149.95_WP, 530.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 8.4474E00_WP, 2.6178E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A = 6.5935E01_WP;
    pc%vapPressure%B = -6.1985E03_WP;
    pc%vapPressure%C = -6.3423E00_WP;
    pc%vapPressure%D =  2.0909E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 149.95_WP, 530.60_WP /);
    pc%vapPressure%range(:,2) = (/ 7.6265E-04_WP, 3.0416E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.8430E07_WP;
    pc%HeatOfVapMol%B = 3.6670E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 149.95_WP, 530.6_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.2877E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.8140E04_WP;
    pc%solHeatCapacityMol%B = 1.6964E03_WP;
    pc%solHeatCapacityMol%C = -8.0340E00_WP;
    pc%solHeatCapacityMol%D = 1.8100E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 30.00_WP, 149.95_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.6010E04_WP, 1.1662E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.9374E05_WP;
    pc%liqHeatCapacityMol%B = -1.2470E02_WP;
    pc%liqHeatCapacityMol%C = 7.3500E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 149.95_WP, 363.23_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.9157E05_WP, 2.4542E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.2030E05_WP;
    pc%IG_HeatCapacityMol%B = 3.4123E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5455E03_WP;
    pc%IG_HeatCapacityMol%D = 2.0480E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2320E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.2805E05_WP, 4.0072E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.4620E-01_WP;
    pc%SecondViralCoefMol%B = -2.6016E02_WP;
    pc%SecondViralCoefMol%C = -3.9596E07_WP;
    pc%SecondViralCoefMol%D = 1.0688E19_WP;
    pc%SecondViralCoefMol%E = -1.1007E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 265.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.0957E00_WP, 6.1028E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -8.2119E00_WP;
    pc%liqViscosity%B = 8.9352E02_WP;
    pc%liqViscosity%C = -4.5940E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 149.95_WP, 363.23_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.0516E-02_WP, 2.1174E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 4.0850E-08_WP;
    pc%vapViscosity%B = 8.9530E-01_WP;
    pc%vapViscosity%C = 0.0_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 149.95_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.6251E-06_WP, 1.9819E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1370E-01_WP;
    pc%liqThermalConductivity%B = -2.9200E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 149.95_WP, 420.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6991E-01_WP, 9.1060E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 7.9334E-02_WP;
    pc%vapThermalConductivity%B = 2.8330E-01_WP;
    pc%vapThermalConductivity%C = 4.5613E03_WP;
    pc%vapThermalConductivity%D = 1.3546E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 363.23_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.7690E-02_WP, 8.1192E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.4950E-02_WP;
    pc%SurfaceTension%B = 1.2254E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 149.95_WP, 530.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.6578E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_dinpropylether


  subroutine pc_ethane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'ethane';
    pc%CAS = '74-64- 0';
    ! components
    pc%ChemicalFormula%C = 2;         
    pc%ChemicalFormula%H = 6;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 30.070_WP;   ! kg/kmol    
    pc%Tcrit = 305.32_WP;             ! K    
    pc%pcrit = 4.8720E6_WP;           ! Pa 
    pc%vcrit = 0.14550_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.279_WP; ! -


    pc%MeltingPoint = 90.35_WP;              ! K
    pc%TripletPointT = 90.35_WP;             ! K
    pc%TripletPointp = 1.1300E00_WP;          ! Pa
    pc%NormalBoilingPoint = 184.55_WP;        ! K
    pc%liqMolarVolume = 0.055229_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -8.3820E07_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -3.1920E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  2.2912E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  2.6590E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.4286E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.0995_WP;            ! -
    pc%RadiusOfGyration = 1.8260E-10_WP;      ! m
    pc%SolubilityParameter = 1.2400E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                 ! C^m
    pc%VanDerWaalsVolume = 0.02734_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 4.2400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.1849_WP;           ! -
    pc%FlashPoint = MAX_REAL_WP;                   ! K
    pc%FlammabilityLimits = (/ 2.9_WP, 13.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 745.00_WP;          ! K


    pc%solDensityMol%A = 2.5263E01_WP;
    pc%solDensityMol%B = -1.0950E-02_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 23.15_WP, 90.35_WP /);
    pc%solDensityMol%range(:,2) = (/ 2.5010E01_WP, 2.4274E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.9122E00_WP;
    pc%liqDensityMol%B = 2.7937E-01_WP;
    pc%liqDensityMol%C = 3.0532E02_WP;
    pc%liqDensityMol%D = 2.9187E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 90.35_WP, 305.32_WP /);
    pc%liqDensityMol%range(:,2) = (/ 2.1640E01_WP, 6.8447E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.1857E01_WP;
    pc%vapPressure%B = -2.5987E03_WP;
    pc%vapPressure%C = -5.1283E00_WP;
    pc%vapPressure%D =  1.4913E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 90.35_WP, 305.32_WP /);
    pc%vapPressure%range(:,2) = (/ 1.1273E00_WP, 4.8522E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 2.1091E07_WP;
    pc%HeatOfVapMol%B = 6.0646E-01_WP;
    pc%HeatOfVapMol%C = -5.5492E-01_WP;
    pc%HeatOfVapMol%D = 3.2799E-01_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 90.35_WP, 305.32_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 1.7879E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 8.1584E03_WP;
    pc%solHeatCapacityMol%B = -1.4778E03_WP;
    pc%solHeatCapacityMol%C = 9.7183E01_WP;
    pc%solHeatCapacityMol%D = -1.6114E00_WP;
    pc%solHeatCapacityMol%E = 9.1233E-03_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 15.53_WP, 63.03_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.1420E03_WP, 4.1591E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 4.4009E01_WP;
    pc%liqHeatCapacityMol%B = 8.9718E04_WP;
    pc%liqHeatCapacityMol%C = 9.1877E02_WP;
    pc%liqHeatCapacityMol%D = -1.8860E03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 92.00_WP, 290.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 6.8554E04_WP, 1.2444E05_WP /);
    pc%liqHeatCapacityMol%eqn = 114;    

    pc%IG_HeatCapacityMol%A = 4.0326E04_WP;
    pc%IG_HeatCapacityMol%B = 1.3422E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6555E03_WP;
    pc%IG_HeatCapacityMol%D = 7.3223E04_WP;
    pc%IG_HeatCapacityMol%E = 7.5287E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 4.2557E04_WP, 1.4562E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 8.0950E-02_WP;
    pc%SecondViralCoefMol%B = -6.1710E01_WP;
    pc%SecondViralCoefMol%C = -1.4350E06_WP;
    pc%SecondViralCoefMol%D = 6.7600E16_WP;
    pc%SecondViralCoefMol%E = -9.7400E18_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 152.71_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.1318E-01_WP, 3.9385E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -3.4134E00_WP;
    pc%liqViscosity%B = 1.9705E02_WP;
    pc%liqViscosity%C = -1.2193E00_WP;
    pc%liqViscosity%D = -9.2023E-26_WP;
    pc%liqViscosity%E = 1E01_WP;
    pc%liqViscosity%range(:,1) = (/ 90.35_WP, 300.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.2020E-03_WP, 3.5196E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.1005E-08_WP;
    pc%vapViscosity%B = 9.7620E-01_WP;
    pc%vapViscosity%C = 2.3139E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) =(/ 90.35_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) =(/ 2.6429E-06_WP, 2.5827E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 3.5758E-01_WP;
    pc%liqThermalConductivity%B = -1.1458E-03_WP;
    pc%liqThermalConductivity%C = 6.1866E-07_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 90.35_WP, 300.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.5911E-01_WP, 6.9519E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 7.3869E-05_WP;
    pc%vapThermalConductivity%B = 1.1689E00_WP;
    pc%vapThermalConductivity%C = 5.0073E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 164.55_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 8.8628E-03_WP, 1.5807E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 4.8643E-02_WP;
    pc%SurfaceTension%B = 1.1981E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 90.35_WP, 305.32_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1949E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_ethane


  subroutine pc_ethanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'ethanol';
    pc%CAS = '64-17-5';
    ! components
    pc%ChemicalFormula%C = 2;         
    pc%ChemicalFormula%H = 6;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 46.069_WP;   ! kg/kmol    
    pc%Tcrit = 513.92_WP;             ! K    
    pc%pcrit = 6.1480E6_WP;           ! Pa 
    pc%vcrit = 0.16700_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.240_WP; ! -


    pc%MeltingPoint = 159.05_WP;              ! K
    pc%TripletPointT = 159.05_WP;             ! K
    pc%TripletPointp = 4.8459E-4_WP;          ! Pa
    pc%NormalBoilingPoint = 351.44_WP;        ! K
    pc%liqMolarVolume = 0.058618_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.3495E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.6785E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  2.8064E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  4.9310E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -1.2350E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.6452_WP;            ! -
    pc%RadiusOfGyration = 2.2590E-10_WP;      ! m
    pc%SolubilityParameter = 2.6130E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.6400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.03194_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 4.9300E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3594_WP;           ! -
    pc%FlashPoint = 286.00_WP;                ! K
    pc%FlammabilityLimits = (/ 4.3_WP, 19.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 696.00_WP;          ! K


    pc%solDensityMol%A = 2.29001E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 159.05_WP, 2.29E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.6480E00_WP;
    pc%liqDensityMol%B = 2.7627E-01_WP;
    pc%liqDensityMol%C = 5.1391E02_WP;
    pc%liqDensityMol%D = 2.3310E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 159.05_WP, 513.92_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.9413E01_WP, 5.9652E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.4475E01_WP;
    pc%vapPressure%B = -7.1643E03_WP;
    pc%vapPressure%C = -7.3270E00_WP;
    pc%vapPressure%D =  3.1340E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 159.05_WP, 513.92_WP /);
    pc%vapPressure%range(:,2) = (/ 4.8459E-04_WP, 6.1171E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.6900E07_WP;
    pc%HeatOfVapMol%B = 3.3590E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 159.05_WP, 513.92_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.0245E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.3500E04_WP;
    pc%solHeatCapacityMol%B = 1.1755E03_WP;
    pc%solHeatCapacityMol%C = -8.0430E00_WP;
    pc%solHeatCapacityMol%D = 2.3700E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 25.00_WP, 159.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.1231E04_WP, 6.5336E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.0264E05_WP;
    pc%liqHeatCapacityMol%B = -1.3963E02_WP;
    pc%liqHeatCapacityMol%C = -3.0341E-02_WP;
    pc%liqHeatCapacityMol%D = 2.0386E-03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 159.05_WP, 390.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 8.7867E04_WP, 1.6450E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 4.8200E04_WP;
    pc%IG_HeatCapacityMol%B = 1.4577E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6628E03_WP;
    pc%IG_HeatCapacityMol%D = 9.3900E04_WP;
    pc%IG_HeatCapacityMol%E = 7.4470E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 5.2236E04_WP, 1.6576E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.4000E-02_WP;
    pc%SecondViralCoefMol%B = -5.5700E01_WP;
    pc%SecondViralCoefMol%C = -1.2900E07_WP;
    pc%SecondViralCoefMol%D = -6.4900E19_WP;
    pc%SecondViralCoefMol%E = -2.5800E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 309.72_WP, 2581.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.3203E00_WP, 2.1669E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = 7.8750E00_WP;
    pc%liqViscosity%B = 7.8198E02_WP;
    pc%liqViscosity%C = -3.0418E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 200.00_WP, 440.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.3148E-02_WP, 1.4160E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.0613E-07_WP;
    pc%vapViscosity%B = 8.0660E-01_WP;
    pc%vapViscosity%C = 5.2700E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 200.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 6.0295E-06_WP, 2.6505E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.4680E-01_WP;
    pc%liqThermalConductivity%B = -2.6400E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 159.05_WP, 353.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.0481E-01_WP, 1.5357E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -1.0109E-02_WP;
    pc%vapThermalConductivity%B = 6.4750E-01_WP;
    pc%vapThermalConductivity%C = -7.3320E03_WP;
    pc%vapThermalConductivity%D = -2.6800E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 293.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.4747E-02_WP, 1.3417E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 3.7640E-02_WP;
    pc%SurfaceTension%B = -2.1570E-05_WP;
    pc%SurfaceTension%C = -1.0250E-07_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 273.15_WP, 503.15_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.4101E-02_WP, 8.3816E-04_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_ethanol


  subroutine pc_gammavalerolactone(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'gammavalerolactone';
    pc%CAS = '108-29-2';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 8;
    pc%ChemicalFormula%O = 2;

    pc%MolecularWeight = 100.117_WP;  ! kg/kmol    
    pc%Tcrit = 727.00_WP;             ! K    
    pc%pcrit = 4.8400E6_WP;           ! Pa 
    pc%vcrit = 0.21900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.223_WP; ! -


    pc%MeltingPoint = 242.15_WP;              ! K
    pc%TripletPointT = 242.15_WP;             ! K
    pc%TripletPointp = 3.1332E-01_WP;         ! Pa
    pc%NormalBoilingPoint = 480.65_WP;        ! K
    pc%liqMolarVolume = 0.095430_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -4.1800E08_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -2.9700E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.5000E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                 ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.4929E09_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    


    pc%AcentricFactor = 0.4030_WP;            ! -
    pc%RadiusOfGyration = 3.4319E-10_WP;      ! m
    pc%SolubilityParameter = 2.287E04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 1.4343E-29_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05631_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.5900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4330_WP;           ! -
    pc%FlashPoint = 354.00_WP;                ! K
    pc%FlammabilityLimits = (/ MAX_REAL_WP, 10.6_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = MAX_REAL_WP;             ! K


    pc%solDensityMol%A = 1.2900E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 242.15_WP, 1.2900E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.8680E-01_WP;
    pc%liqDensityMol%B = 2.7532E-01_WP;
    pc%liqDensityMol%C = 7.2700E02_WP;
    pc%liqDensityMol%D = 3.4900E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 242.15_WP, 727.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0983E01_WP, 3.5842E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  8.3640E01_WP;
    pc%vapPressure%B = -8.7843E03_WP;
    pc%vapPressure%C = -8.8878E00_WP;
    pc%vapPressure%D =  4.5338E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 242.15_WP, 727.00_WP /);
    pc%vapPressure%range(:,2) = (/ 3.1332E-01_WP, 4.8391E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.1320E07_WP;
    pc%HeatOfVapMol%B = 2.9680E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 242.15_WP, 727.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.4373E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 1.4830E05_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 242.15_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.4830E05_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 9.5800E04_WP;
    pc%liqHeatCapacityMol%B = 2.7750E02_WP;
    pc%liqHeatCapacityMol%C = 0.0_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 363.50_WP, 543.50_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.9667E05_WP, 2.4662E05_WP /);
    pc%liqHeatCapacityMol%eqn = 107;    

    pc%IG_HeatCapacityMol%A = 8.9700E04_WP;
    pc%IG_HeatCapacityMol%B = 1.8340E05_WP;
    pc%IG_HeatCapacityMol%C = 6.7500E02_WP;
    pc%IG_HeatCapacityMol%D = 9.3000E04_WP;
    pc%IG_HeatCapacityMol%E = 2.2500E03_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 300.00_WP , 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3190E05_WP, 2.9902E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.1157E-01_WP;
    pc%SecondViralCoefMol%B = -2.9963E02_WP;
    pc%SecondViralCoefMol%C = -8.9927E07_WP;
    pc%SecondViralCoefMol%D = -2.9311E20_WP;
    pc%SecondViralCoefMol%E = -2.6964E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 363.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.6900E00_WP, -1.4841E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.8539E00_WP;
    pc%liqViscosity%B = 5.6334E02_WP;
    pc%liqViscosity%C = -1.1445E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 242.15_WP, 572.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.1024E-03_WP, 1.0770E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 8.0120E-08_WP;
    pc%vapViscosity%B = 8.5930E-01_WP;
    pc%vapViscosity%C = 6.1350E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) =(/ 242.15_WP, 1000.0_WP /);
    pc%vapViscosity%range(:,2) =(/ 7.1501E-06_WP, 2.8562E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2568E-01_WP;
    pc%liqThermalConductivity%B = -2.1770E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 242.15_WP, 480.65_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7296E-01_WP, 1.2104E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 3.8300E-04_WP;
    pc%vapThermalConductivity%B = 8.7240E-01_WP;
    pc%vapThermalConductivity%C = 7.0550E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 480.65_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 3.3925E-02_WP, 1.5807E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.6815E-02_WP;
    pc%SurfaceTension%B = 1.3153E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 242.15_WP, 727.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.9217E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_gammavalerolactone


  subroutine pc_isobutane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'isobutane';
    pc%CAS = '75-28-5';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 58.123_WP;   ! kg/kmol    
    pc%Tcrit = 407.80_WP;             ! K    
    pc%pcrit = 3.6400E6_WP;           ! Pa 
    pc%vcrit = 0.25900_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.278_WP; ! -


    pc%MeltingPoint = 113.54_WP;              ! K
    pc%TripletPointT = 113.54_WP;             ! K
    pc%TripletPointp = 1.2062E-2_WP;          ! Pa
    pc%NormalBoilingPoint = 261.43_WP;        ! K
    pc%liqMolarVolume = 0.105350_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -1.3499E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -2.1440E7_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  2.9550E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  4.5400E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.6490E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -2.0760E07_WP;      ! J/kmol
    pc%StandardHeatOfFormation = - 1.3418E08_WP;      ! J/kmol
    pc%StandardAbsoluteEntropy = 2.9539E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.1835_WP;            ! -
    pc%RadiusOfGyration = 2.9480E-10_WP;      ! m
    pc%SolubilityParameter = 1.2560E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 4.4000E-31_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.04779_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 6.9300E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3175_WP;           ! -
    pc%FlashPoint = MAX_REAL_WP;                   ! K
    pc%FlammabilityLimits = (/ 1.8_WP, 8.4_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 191.0_WP, 212.0_WP /); ! K
    pc%AutoignitionTemp = 733.15_WP;          ! K


    pc%solDensityMol%A = 1.4700E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 113.54_WP, 1.4700E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.0631E00_WP;
    pc%liqDensityMol%B = 2.7506E-01_WP;
    pc%liqDensityMol%C = 4.0780E02_WP;
    pc%liqDensityMol%D = 2.7580E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 113.54_WP, 407.80_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.2574E01_WP, 3.8650E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.0843E02_WP;
    pc%vapPressure%B = -5.0399E03_WP;
    pc%vapPressure%C = -1.5012E01_WP;
    pc%vapPressure%D =  2.2725E-02_WP;
    pc%vapPressure%E =  1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 113.54_WP, 407.80_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2062E-02_WP, 3.6299E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 3.1880E07_WP;
    pc%HeatOfVapMol%B = 3.9006E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 113.54_WP, 407.80_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 2.8070E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 1.2546E05_WP;
    pc%solHeatCapacityMol%B = -2.1580E03_WP;
    pc%solHeatCapacityMol%C = 1.5616E01_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 77.00_WP, 107.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 5.1881E04_WP, 7.3342E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.7237E05_WP;
    pc%liqHeatCapacityMol%B = -1.7839E03_WP;
    pc%liqHeatCapacityMol%C = 1.4759E01_WP;
    pc%liqHeatCapacityMol%D = -4.7909E-02_WP;
    pc%liqHeatCapacityMol%E = 5.8050E-05_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 113.54_WP, 380.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 9.9613E04_WP, 2.0725E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 6.5490E04_WP;
    pc%IG_HeatCapacityMol%B = 2.4776E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5870E03_WP;
    pc%IG_HeatCapacityMol%D = 1.5750E05_WP;
    pc%IG_HeatCapacityMol%E = -7.0699E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 7.2181E04_WP, 2.6656E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.4990E-01_WP;
    pc%SecondViralCoefMol%B = -1.3056E02_WP;
    pc%SecondViralCoefMol%C = -8.7980E06_WP;
    pc%SecondViralCoefMol%D = -9.9500E17_WP;
    pc%SecondViralCoefMol%E = 5.2500E19_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 204.70_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.7532E00_WP, 6.0253E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.3912E01_WP;
    pc%liqViscosity%B = 7.9709E02_WP;
    pc%liqViscosity%C = 4.5508E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 110.00_WP, 310.95_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.0716E-02_WP, 1.5877E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.0871E-07_WP;
    pc%vapViscosity%B = 7.8135E-01_WP;
    pc%vapViscosity%C = 7.0639E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 150.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.7065E-06_WP, 2.2422E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0455E-01_WP;
    pc%liqThermalConductivity%B = -3.6589E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 113.54_WP, 400.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6301E-01_WP, 5.8194E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 8.9772E-02_WP;
    pc%vapThermalConductivity%B = 1.8501E-01_WP;
    pc%vapThermalConductivity%C = 6.3923E02_WP;
    pc%vapThermalConductivity%D = 1.1147E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 261.43_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.2726E-02_WP, 1.1701E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1359E-02_WP;
    pc%SurfaceTension%B = 1.2532E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 113.54_WP , 407.80_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.4121E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_isobutane


  subroutine pc_isoheptane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'isoheptane';
    pc%CAS = '591-76-4';
    ! components
    pc%ChemicalFormula%C = 7;         
    pc%ChemicalFormula%H = 16;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 100.204_WP;  ! kg/kmol    
    pc%Tcrit = 530.4_WP;              ! K    
    pc%pcrit = 2.7400E6_WP;           ! Pa 
    pc%vcrit = 0.421_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.262_WP; ! -


    pc%MeltingPoint = 154.9_WP;              ! K
    pc%TripletPointT = 154.9_WP;             ! K
    pc%TripletPointp = 4.2989E-03_WP;         ! Pa
    pc%NormalBoilingPoint = 363.199_WP;      ! K
    pc%liqMolarVolume = 0.14869_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -1.9460E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 3.4700E06_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  4.2041E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  9.1840E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.4596E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = -2.4762E06_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -2.2949E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 3.2334E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.328006_WP;          ! -
    pc%RadiusOfGyration = 4.1670E-10_WP;      ! m
    pc%SolubilityParameter = 1.4730E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.07848_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.0980E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.38228_WP;          ! -
    pc%FlashPoint = 250.0_WP;                 ! K
    pc%FlammabilityLimits = (/ 1.0_WP, 6.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 260.0_WP, 291.0_WP /); ! K
    pc%AutoignitionTemp = 553.15_WP;          ! K


    pc%solDensityMol%A = 9.1400E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 154.90_WP, 9.1400E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.3993E-01_WP;
    pc%liqDensityMol%B = 2.6950E-01_WP;
    pc%liqDensityMol%C = 5.3040E02_WP;
    pc%liqDensityMol%D = 2.7930E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 154.90_WP, 530.40_WP /);
    pc%liqDensityMol%range(:,2) = (/ 7.8100E00_WP, 2.3745E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.9531E01_WP;
    pc%vapPressure%B = -5.8264E03_WP;
    pc%vapPressure%C = -5.4269E00_WP;
    pc%vapPressure%D =  1.4542E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 154.90_WP, 530.40_WP /);
    pc%vapPressure%range(:,2) = (/ 4.2989E-03_WP, 2.7401E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.7126E07_WP;
    pc%HeatOfVapMol%B = 3.6975E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 154.90_WP, 530.40_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.1476E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 3.7436E04_WP;
    pc%solHeatCapacityMol%B = 4.1437E02_WP;
    pc%solHeatCapacityMol%C = 1.3235E00_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 80.00_WP, 154.90_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 7.9068E04_WP, 1.3342E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.7401E05_WP;
    pc%liqHeatCapacityMol%B = -1.0578E02_WP;
    pc%liqHeatCapacityMol%C = 9.0525E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 154.90_WP, 363.20_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.7935E05_WP, 2.5501E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0600E05_WP;
    pc%IG_HeatCapacityMol%B = 4.2507E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5770E03_WP;
    pc%IG_HeatCapacityMol%D = 2.7580E05_WP;
    pc%IG_HeatCapacityMol%E = 6.9080E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.1915E05_WP, 4.5148E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.4123E-01_WP;
    pc%SecondViralCoefMol%B = -2.9400E02_WP;
    pc%SecondViralCoefMol%C = -3.8010E07_WP;
    pc%SecondViralCoefMol%D = -4.3597E19_WP;
    pc%SecondViralCoefMol%E = 6.3133E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 265.20_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.7141E00_WP, 3.3966E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.2224E01_WP;
    pc%liqViscosity%B = 1.0213E03_WP;
    pc%liqViscosity%C = 1.5187E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 154.90_WP, 363.20_WP /);
    pc%liqViscosity%range(:,2) = (/ 7.7129E-03_WP, 2.0000E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.0132E-06_WP;
    pc%vapViscosity%B = 4.5610E-01_WP;
    pc%vapViscosity%C = 3.5980E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 154.90_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.0414E-06_WP, 1.7399E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.9690E-01_WP;
    pc%liqThermalConductivity%B = -2.6920E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 154.90_WP, 420.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5520E-01_WP, 8.3836E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.9926E-05_WP;
    pc%vapThermalConductivity%B = 1.0785E00_WP;
    pc%vapThermalConductivity%C = -8.5924E01_WP;
    pc%vapThermalConductivity%D = 1.5868E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 363.20_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.7582E-02_WP, 9.6076E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0288E-02_WP;
    pc%SurfaceTension%B = 1.1923E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 154.90_WP, 530.40_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.3314E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_isoheptane


  subroutine pc_isohexane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'isohexane';
    pc%CAS = '107-83-5';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 86.1772_WP;  ! kg/kmol    
    pc%Tcrit = 497.7_WP;              ! K    
    pc%pcrit = 3.0400E6_WP;           ! Pa 
    pc%vcrit = 0.368_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.27_WP;  ! -


    pc%MeltingPoint = 119.55_WP;              ! K
    pc%TripletPointT = 119.55_WP;             ! K
    pc%TripletPointp = 2.0656E-5_WP;          ! Pa
    pc%NormalBoilingPoint = 333.41_WP;        ! K
    pc%liqMolarVolume = 0.132925_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -1.7455E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -5.3380E6_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.8089E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  6.2680E6_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion =-3.8492E9_WP;     ! J/kmol
    pc%StandardGibbsOfFormation = -8.4999E06_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -2.0464E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.9058E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.279149_WP;          ! -
    pc%RadiusOfGyration = 3.7840E-10_WP;      ! m
    pc%SolubilityParameter = 1.4400E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.06825_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 9.6300E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.36873_WP;          ! -
    pc%FlashPoint = 238.0_WP;                 ! K
    pc%FlammabilityLimits = (/ 1.2_WP, 7.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 239.0_WP, 269.0_WP /); ! K
    pc%AutoignitionTemp = 579.0_WP;           ! K


    pc%solDensityMol%A = 1.0800E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 119.55_WP,1.0800E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 7.2701E-01_WP;
    pc%liqDensityMol%B = 2.6754E-01_WP;
    pc%liqDensityMol%C = 4.9770E02_WP;
    pc%liqDensityMol%D = 2.8268E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 119.55_WP, 497.70_WP /);
    pc%liqDensityMol%range(:,2) = (/ 9.2041E00_WP, 2.7174E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.3579E01_WP;
    pc%vapPressure%B = -5.0412E03_WP;
    pc%vapPressure%C = -4.6404E00_WP;
    pc%vapPressure%D =  1.9443E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 119.55_WP, 497.70_WP /);
    pc%vapPressure%range(:,2) = (/ 2.0656E-05_WP, 3.0439E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.2522E07_WP;
    pc%HeatOfVapMol%B = 3.8070E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 119.55_WP, 497.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.8300E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.6127E04_WP;
    pc%solHeatCapacityMol%B = 1.6109E03_WP;
    pc%solHeatCapacityMol%C = -9.0041E00_WP;
    pc%solHeatCapacityMol%D = 2.6825E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 8.00_WP, 119.53_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0108E04_WP, 9.3590E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.4222E05_WP;
    pc%liqHeatCapacityMol%B = -4.7830E01_WP;
    pc%liqHeatCapacityMol%C = 7.3900E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 119.55_WP, 333.41_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.4706E05_WP, 2.0842E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 9.0300E04_WP;
    pc%IG_HeatCapacityMol%B = 3.8010E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6020E03_WP;
    pc%IG_HeatCapacityMol%D = 2.4530E05_WP;
    pc%IG_HeatCapacityMol%E = -6.9160E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.0192E05_WP, 3.9617E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.0487E-01_WP;
    pc%SecondViralCoefMol%B = -2.0373E02_WP;
    pc%SecondViralCoefMol%C = -3.1708E07_WP;
    pc%SecondViralCoefMol%D = -4.4141E18_WP;
    pc%SecondViralCoefMol%E = -3.7253E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 248.73_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.3958E00_WP, 5.9655E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.2860E01_WP;
    pc%liqViscosity%B = 9.4691E02_WP;
    pc%liqViscosity%C = 2.6191E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 119.55_WP, 333.41_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.5060E-02_WP, 2.0380E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.1164E-06_WP;
    pc%vapViscosity%B = 4.5370E-01_WP;
    pc%vapViscosity%C = 3.7474E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 119.55_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 2.3658E-06_WP, 1.8651E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.9334E-01_WP;
    pc%liqThermalConductivity%B = -2.8038E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 119.55_WP, 389.25_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5982E-01_WP, 8.4202E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 6.1119E-05_WP;
    pc%vapThermalConductivity%B = 1.0861E00_WP;
    pc%vapThermalConductivity%C = -5.9592E01_WP;
    pc%vapThermalConductivity%D = 1.4126E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 333.41_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.6063E-02_WP, 1.0242E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0048E-02_WP;
    pc%SurfaceTension%B = 1.1899E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 119.55_WP, 497.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.6093E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_isohexane

  subroutine pc_isooctane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

    ! all alternative names
    pc%name = 'isooctane';
    pc%CAS = '540-84-1';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 114.231_WP;   ! kg/kmol    
    pc%Tcrit = 543.96_WP;              ! K    
    pc%pcrit = 2.5680E6_WP;            ! Pa 
    pc%vcrit = 0.46800_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.266_WP;  ! -

    pc%MeltingPoint = 166.85_WP;              ! K
    pc%TripletPointT = 165.76_WP;              ! K
    pc%TripletPointp = 0.000E-0_WP;           ! Pa
    pc%NormalBoilingPoint = 372.39_WP;        ! K
    pc%liqMolarVolume = 0.1651_WP;            ! m^3/kmol

    pc%IG_HeatOfFormation = -0.0000E0_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = 0.0000E0_WP;             ! J/kmol
    pc%IG_AbsoluteEntropy =  0.0000E0_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  0.0000E0_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -0.0000E0_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;        ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;         ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;         ! J/kmolK   

    pc%AcentricFactor = 0.3030_WP;            ! -
    pc%RadiusOfGyration = 4.091E-10_WP;       ! m
    pc%SolubilityParameter = 0.0000E0_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                 ! C^m
    pc%VanDerWaalsVolume = 0.00000_WP;        ! m^3/kmol

    pc%VanDerWaalsArea = 0.0000E0_WP;                           ! m^2/kmol
    pc%RefractiveIndex = 0.0000_WP;                             ! -
    pc%FlashPoint = 261.15_WP;                                  ! K
    pc%FlammabilityLimits = (/ 0.0_WP, 0.0_WP /);               ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 690.928_WP;                            ! K

    pc%liqDensityMol%A = 0.24563_WP;
    pc%liqDensityMol%B = 0.27373_WP;
    pc%liqDensityMol%C = 0.28460_WP;
    pc%liqDensityMol%D = 0.0_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 165.78_WP, 543.96_WP /);
    pc%liqDensityMol%eqn = 116;

    pc%vapPressure%A =  50.3422_WP;
    pc%vapPressure%B = -3.2789E03_WP;
    pc%vapPressure%C = -1.6111E01_WP;
    pc%vapPressure%D =  7.4260E-03_WP;
    pc%vapPressure%E = -9.1804E-14_WP;
    pc%vapPressure%range(:,1) = (/ 165.78_WP, 543.96_WP /);
    pc%vapPressure%eqn = 117;

    pc%HeatOfVapMol%A = 42.901_WP;
    pc%HeatOfVapMol%B = 0.281_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 165.78_WP, 543.96_WP /);
    pc%HeatOfVapMol%eqn = 118; 

    pc%liqHeatCapacityMol%A = 122.772_WP;
    pc%liqHeatCapacityMol%B = 7.9485E-01_WP;
    pc%liqHeatCapacityMol%C = -2.4977E-03_WP;
    pc%liqHeatCapacityMol%D = 3.5652E-06_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 167.0_WP, 490.0_WP /);
    pc%liqHeatCapacityMol%eqn = 119;

    pc%vapViscosity%A = 4.636_WP;
    pc%vapViscosity%B = 2.0622E-01_WP;
    pc%vapViscosity%C = -3.2545E-05_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 200.0_WP, 992.0_WP /);    ! Ideally the given viscosity range started from 355 K
    pc%vapViscosity%eqn = 120;

    pc%liqViscosity%A = -15.0420_WP;
    pc%liqViscosity%B = 2.0021E+03_WP;
    pc%liqViscosity%C = 3.7102E-02_WP; 
    pc%liqViscosity%D = -3.4486E-05_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 167.0_WP, 544.0_WP /);   ! Ideally the given viscosity range started from 273 K
    pc%liqViscosity%eqn = 121;

    pc%vapThermalConductivity%A = -0.00622_WP;
    pc%vapThermalConductivity%B = 4.1978E-05_WP;
    pc%vapThermalConductivity%C = 7.7251E-08_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 200.0_WP, 580.0_WP /);  ! Ideally the given thermal conductivity range started from 355 K
    pc%vapThermalConductivity%eqn = 122; 

    pc%liqThermalConductivity%A = -1.6144_WP;
    pc%liqThermalConductivity%B = 0.7698_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 166.0_WP, 517.0_WP /);
    pc%liqThermalConductivity%eqn = 123;

    pc%SurfaceTension%A = 47.434_WP;
    pc%SurfaceTension%B = 1.1975_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 165.78_WP, 543.96_WP /);
    pc%SurfaceTension%eqn = 124;    

    pc%IG_HeatCapacityMol%A = -21.703_WP;
    pc%IG_HeatCapacityMol%B = 8.5849E-01_WP;
    pc%IG_HeatCapacityMol%C = -5.5323E-04_WP;
    pc%IG_HeatCapacityMol%D = 2.0892E-07_WP;
    pc%IG_HeatCapacityMol%E = -3.7285E-11_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%eqn = 125;

  end subroutine pc_isooctane

  subroutine pc_2methylheptane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

    ! all alternative names
    pc%name = '2methylheptane';
    pc%CAS = '107-83-5';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 114.231_WP;   ! kg/kmol    
    pc%Tcrit = 559.64_WP;              ! K    
    pc%pcrit = 2.4840E6_WP;            ! Pa 
    pc%vcrit = 0.48800_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.261_WP;  ! -


    pc%MeltingPoint = 164.16_WP;              ! K
    pc%TripletPointT = 164.16_WP;             ! K
    pc%TripletPointp = 1.0606E-3_WP;          ! Pa
    pc%NormalBoilingPoint = 390.80_WP;        ! K
    pc%liqMolarVolume = 0.164409_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -2.1535E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = 1.1690E7_WP;             ! J/kmol
    pc%IG_AbsoluteEntropy =  4.5957E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  1.1878E7_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -5.0694E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.3772_WP;            ! -
    pc%RadiusOfGyration = 4.5700E-10_WP;      ! m
    pc%SolubilityParameter = 1.5050E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.088710_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.2330E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3926_WP;           ! -
    pc%FlashPoint = 277.15_WP;                ! K
    pc%FlammabilityLimits = (/ 0.9_WP, 5.8_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 520.00_WP;          ! K


    pc%solDensityMol%A = 8.1000E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 164.16_WP, 8.1000E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.2966E-01_WP;
    pc%liqDensityMol%B = 2.6042E-01_WP;
    pc%liqDensityMol%C = 5.5964E02_WP;
    pc%liqDensityMol%D = 2.7016E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 164.16_WP, 559.64_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.9236E00_WP, 2.0339E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  9.0727E01_WP;
    pc%vapPressure%B = -7.4715E03_WP;
    pc%vapPressure%C = -1.0243E01_WP;
    pc%vapPressure%D =  6.8840E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 164.16_WP, 559.64_WP /);
    pc%vapPressure%range(:,2) = (/ 1.0606E-03_WP, 2.4757E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.3089E07_WP;
    pc%HeatOfVapMol%B = 3.8112E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 164.16_WP, 559.64_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.6509E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -2.0180E04_WP;
    pc%solHeatCapacityMol%B = 1.8980E03_WP;
    pc%solHeatCapacityMol%C = -9.1850E00_WP;
    pc%solHeatCapacityMol%D = 2.3100E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 164.16_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.4291E04_WP, 1.4606E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.0489E05_WP;
    pc%liqHeatCapacityMol%B = -1.2073E02_WP;
    pc%liqHeatCapacityMol%C = 9.2846E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 164.16_WP, 500.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.1009E05_WP, 3.7664E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.2215E05_WP;
    pc%IG_HeatCapacityMol%B = 4.7270E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5866E03_WP;
    pc%IG_HeatCapacityMol%D = 3.1650E05_WP;
    pc%IG_HeatCapacityMol%E = 6.9960E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 150.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.2460E05_WP, 5.0777E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.7362E-01_WP;
    pc%SecondViralCoefMol%B = -3.4177E02_WP;
    pc%SecondViralCoefMol%C = -5.4179E07_WP;
    pc%SecondViralCoefMol%D = -9.1933E19_WP;
    pc%SecondViralCoefMol%E = 1.3873E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 279.82_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.5475E00_WP, 2.9717E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.1342E01_WP;
    pc%liqViscosity%B = 1.0740E03_WP;
    pc%liqViscosity%C = 1.3051E-02_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 164.16_WP, 390.80_WP /);
    pc%liqViscosity%range(:,2) = (/ 8.7998E-03_WP, 2.0026E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 4.4600E-07_WP;
    pc%vapViscosity%B = 5.5355E-01_WP;
    pc%vapViscosity%C = 2.2223E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 164.16_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.1903E-06_WP, 1.6704E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.9945E-01_WP;
    pc%liqThermalConductivity%B = -2.5820E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 164.16_WP, 448.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5706E-01_WP, 8.3776E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.4211E-05_WP;
    pc%vapThermalConductivity%B = 1.0809E00_WP;
    pc%vapThermalConductivity%C = -1.3266E02_WP;
    pc%vapThermalConductivity%D = 1.7358E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 390.80_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.9106E-02_WP, 9.1069E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0322E-02_WP;
    pc%SurfaceTension%B = 1.2008E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 164.16_WP, 559.64_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.3166E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_2methylheptane

  subroutine pc_isopentane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'isopentane';
    pc%CAS = '107-83-5';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 72.1503_WP;  ! kg/kmol    
    pc%Tcrit = 460.4_WP;              ! K    
    pc%pcrit = 3.3800E6_WP;           ! Pa 
    pc%vcrit = 0.306_WP;              ! m^3/kmol 
    pc%CritCompressFactor = 0.27_WP;  ! -


    pc%MeltingPoint = 113.25_WP;              ! K
    pc%TripletPointT = 113.25_WP;             ! K
    pc%TripletPointp = 1.2144E-4_WP;          ! Pa
    pc%NormalBoilingPoint = 300.994_WP;       ! K
    pc%liqMolarVolume  = 0.117065_WP;         ! m^3/kmol


    pc%IG_HeatOfFormation = -1.5370E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.4050E7_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.4374E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  5.1500E06_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.2395E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -1.4420E07_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -1.7891E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.6041E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.227875_WP;            ! -
    pc%RadiusOfGyration = 3.3240E-10_WP;      ! m
    pc%SolubilityParameter = 1.3870E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 4.3400E-30_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.05802_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.2800E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.35088_WP;          ! -
    pc%FlashPoint = 217.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.3_WP, 8.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 216.0_WP, 244.0_WP /); ! K
    pc%AutoignitionTemp = 693.15_WP;          ! K


    pc%solDensityMol%A = 1.2600E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 113.25_WP, 1.2600E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.1991E-01_WP;
    pc%liqDensityMol%B = 2.7815E-01_WP;
    pc%liqDensityMol%C = 4.6040E02_WP;
    pc%liqDensityMol%D = 2.8667E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 113.25_WP, 460.40_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0764E01_WP, 3.3072E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.1308E01_WP;
    pc%vapPressure%B = -4.9760E03_WP;
    pc%vapPressure%C = -7.7169E00_WP;
    pc%vapPressure%D =  8.7271E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 113.25_WP, 460.40_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2144E-04_WP, 3.3656E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 3.7593E07_WP;
    pc%HeatOfVapMol%B = 3.9173E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 113.25_WP, 460.40_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.3657E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -8.1500E03_WP;
    pc%solHeatCapacityMol%B = 8.5500E02_WP;
    pc%solHeatCapacityMol%C = 1.4470E01_WP;
    pc%solHeatCapacityMol%D = -2.7730E-01_WP;
    pc%solHeatCapacityMol%E = 1.3350E-03_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 113.36_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.2733E04_WP, 9.1223E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.0830E05_WP;
    pc%liqHeatCapacityMol%B = 1.4600E02_WP;
    pc%liqHeatCapacityMol%C = -2.9200E-01_WP;
    pc%liqHeatCapacityMol%D = 1.5100E-03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 113.25_WP, 310.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.2328E05_WP, 1.7048E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 7.4600E04_WP;
    pc%IG_HeatCapacityMol%B = 3.2650E05_WP;
    pc%IG_HeatCapacityMol%C = 1.5450E03_WP;
    pc%IG_HeatCapacityMol%D = 1.9230E05_WP;
    pc%IG_HeatCapacityMol%E = 6.6670E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 8.5462E04_WP, 3.3792E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.8170E-01_WP;
    pc%SecondViralCoefMol%B = -1.9553E02_WP;
    pc%SecondViralCoefMol%C = -1.5973E07_WP;
    pc%SecondViralCoefMol%D = -4.6950E18_WP;
    pc%SecondViralCoefMol%E =3.6600E20_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 230.20_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.3708E00_WP, 4.6614E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.2596E01_WP;
    pc%liqViscosity%B = 8.8911E02_WP;
    pc%liqViscosity%C = 2.0469E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 150.00_WP, 310.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.5423E-03_WP, 1.9283E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.4344E-08_WP;
    pc%vapViscosity%B = 9.7376E-01_WP;
    pc%vapViscosity%C = -9.1597E01_WP;
    pc%vapViscosity%D = 1.8720E04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 150.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 2.6214E-06_WP, 2.1905E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1246E-01_WP;
    pc%liqThermalConductivity%B = -3.3581E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 113.25_WP, 368.13_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7443E-01_WP, 8.8838E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 8.9680E-04_WP;
    pc%vapThermalConductivity%B = 7.7420E-01_WP;
    pc%vapThermalConductivity%C = 4.5600E02_WP;
    pc%vapThermalConductivity%D = 2.3064E05_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.15_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.1981E-02_WP, 1.1176E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0876E-02_WP;
    pc%SurfaceTension%B = 1.2066E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 113.25_WP, 460.40_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.6188E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_isopentane


  subroutine pc_methane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'methane';
    pc%CAS = '74-82-8';
    ! components
    pc%ChemicalFormula%C = 1;         
    pc%ChemicalFormula%H = 4;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 16.043_WP;    ! kg/kmol    
    pc%Tcrit = 190.56_WP;              ! K    
    pc%pcrit = 4.5990E06_WP;           ! Pa 
    pc%vcrit = 0.09860_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.286_WP;  ! -


    pc%MeltingPoint = 90.69_WP;              ! K
    pc%TripletPointT = 90.69_WP;             ! K
    pc%TripletPointp = 1.1696E04_WP;         ! Pa
    pc%NormalBoilingPoint = 111.66_WP;       ! K
    pc%liqMolarVolume  = 0.037969_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -7.4520E07_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -5.0490E07_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy =  1.8627E05_WP;            ! J/kmolK
    pc%HeatFusionAtMeltPoint =  9.4140E05_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -8.0262E08_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.0115_WP;            ! -
    pc%RadiusOfGyration = 1.1180E-10_WP;      ! m
    pc%SolubilityParameter = 1.1600E04_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.01705_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 2.8800E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.0004_WP;           ! -
    pc%FlashPoint = MAX_REAL_WP;                   ! K
    pc%FlammabilityLimits = (/ 5.0_WP, 15.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 810.00_WP;          ! K


    pc%solDensityMol%A = 3.3022E01_WP;
    pc%solDensityMol%B = -1.5870E-02_WP;
    pc%solDensityMol%C = -1.5500E-04_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 23.15_WP, 90.67_WP /);
    pc%solDensityMol%range(:,2) = (/ 3.2572E01_WP, 3.0309E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 2.9214E00_WP;
    pc%liqDensityMol%B = 2.8976E-01_WP;
    pc%liqDensityMol%C = 1.9056E02_WP;
    pc%liqDensityMol%D = 2.8818E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 90.69_WP, 190.56_WP /);
    pc%liqDensityMol%range(:,2) = (/ 2.8180E01_WP, 1.0082E01_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  3.9205E01_WP;
    pc%vapPressure%B = -1.3244E03_WP;
    pc%vapPressure%C = -3.4366E00_WP;
    pc%vapPressure%D =  3.1019E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 90.69_WP, 190.56_WP /);
    pc%vapPressure%range(:,2) = (/ 1.1687E04_WP, 4.5897E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 1.0194E07_WP;
    pc%HeatOfVapMol%B = 2.6087E-01_WP;
    pc%HeatOfVapMol%C = -1.4694E-01_WP;
    pc%HeatOfVapMol%D = 2.2154E-01_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 90.69_WP, 190.56_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 8.7235E06_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -2.5500E03_WP;
    pc%solHeatCapacityMol%B = 1.2499E03_WP;
    pc%solHeatCapacityMol%C = -1.4170E01_WP;
    pc%solHeatCapacityMol%D = 6.6480E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 22.85_WP, 90.67_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.9405E04_WP, 4.3841E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 6.5708E01_WP;
    pc%liqHeatCapacityMol%B = 3.8883E04_WP;
    pc%liqHeatCapacityMol%C = -2.5795E02_WP;
    pc%liqHeatCapacityMol%D = 6.1407E02_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 90.69_WP, 190.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 5.3605E04_WP, 1.4978E06_WP /);
    pc%liqHeatCapacityMol%eqn = 114;    

    pc%IG_HeatCapacityMol%A = 3.3298E04_WP;
    pc%IG_HeatCapacityMol%B = 7.9933E04_WP;
    pc%IG_HeatCapacityMol%C = 2.0869E03_WP;
    pc%IG_HeatCapacityMol%D = 4.1602E04_WP;
    pc%IG_HeatCapacityMol%E = 9.9196E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 50.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 3.3298E04_WP, 8.8904E04_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 5.1075E-02_WP;
    pc%SecondViralCoefMol%B = -2.5181E01_WP;
    pc%SecondViralCoefMol%C = -2.5601E05_WP;
    pc%SecondViralCoefMol%D = 5.9777E15_WP;
    pc%SecondViralCoefMol%E =-5.7697E17_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 110.83_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.3028E-01_WP, 3.4212E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -6.1572E00_WP;
    pc%liqViscosity%B = 1.7815E02_WP;
    pc%liqViscosity%C = -9.5239E-01_WP;
    pc%liqViscosity%D = -9.0606E-24_WP;
    pc%liqViscosity%E = 1.0E01_WP;
    pc%liqViscosity%range(:,1) = (/ 90.69_WP, 188.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.0634E-04_WP, 2.2625E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 5.2546E-07_WP;
    pc%vapViscosity%B = 5.9006E-01_WP;
    pc%vapViscosity%C = 1.0567E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 90.69_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.4684E-06_WP, 2.7996E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 4.1768E-01_WP;
    pc%liqThermalConductivity%B = -2.4528E-03_WP;
    pc%liqThermalConductivity%C = 3.5588E-06_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 90.69_WP, 180.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.2451E-01_WP, 9.1481E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 6.3252E03_WP;
    pc%vapThermalConductivity%B = 4.3041E-01_WP;
    pc%vapThermalConductivity%C = 7.7040E08_WP;
    pc%vapThermalConductivity%D = -3.8725E10_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 90.69_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.1625E-02_WP, 1.6904E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 3.6557E-02_WP;
    pc%SurfaceTension%B = 1.1466E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 90.69_WP, 190.56_WP /);
    pc%SurfaceTension%range(:,2) = (/ 1.7428E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_methane


  subroutine pc_methanol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'methanol';
    pc%CAS = '67-56-1';
    ! components
    pc%ChemicalFormula%C = 1;         
    pc%ChemicalFormula%H = 4;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 32.042_WP;    ! kg/kmol    
    pc%Tcrit = 512.64_WP;              ! K    
    pc%pcrit = 8.0970E06_WP;           ! Pa 
    pc%vcrit = 0.11800_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.224_WP;  ! -


    pc%MeltingPoint = 175.47_WP;             ! K
    pc%TripletPointT = 175.47_WP;            ! K
    pc%TripletPointp = 1.1147E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 337.85_WP;       ! K
    pc%liqMolarVolume  = 0.040581_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.0094E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -1.6232E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy =  2.3988E05_WP;            ! J/kmolK
    pc%HeatFusionAtMeltPoint =  3.2150E06_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -6.38200E8_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.5640_WP;            ! -
    pc%RadiusOfGyration = 1.5520E-10_WP;      ! m
    pc%SolubilityParameter = 2.9590E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 5.6700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.02171_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 3.5800E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3265_WP;           ! -
    pc%FlashPoint = 284.00_WP;                ! K
    pc%FlammabilityLimits = (/ 7.3_WP, 36.0_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 737.00_WP;          ! K


    pc%solDensityMol%A = 3.0585E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 163.15_WP, 90.67_WP /);
    pc%solDensityMol%range(:,2) = (/ 3.0585E01_WP, 3.0309E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 2.2880E00_WP;
    pc%liqDensityMol%B = 2.6850E-01_WP;
    pc%liqDensityMol%C = 5.1264E02_WP;
    pc%liqDensityMol%D = 2.4530E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 175.47_WP, 512.64_WP /);
    pc%liqDensityMol%range(:,2) = (/ 2.7912E01_WP, 8.5214E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  8.1768E01_WP;
    pc%vapPressure%B = -6.8760E03_WP;
    pc%vapPressure%C = -8.7078E00_WP;
    pc%vapPressure%D =  7.1926E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 175.47_WP, 512.64_WP /);
    pc%vapPressure%range(:,2) = (/ 1.1147E-01_WP, 8.1402E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.2390E07_WP;
    pc%HeatOfVapMol%B = 3.6820E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 175.47_WP, 512.64_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.4900E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.2706E04_WP;
    pc%solHeatCapacityMol%B = 9.5878E02_WP;
    pc%solHeatCapacityMol%C = -5.2332E00_WP;
    pc%solHeatCapacityMol%D = 1.3152E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 25.00_WP, 150.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 8.1982E03_WP, 5.7752E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.0580E05_WP;
    pc%liqHeatCapacityMol%B = -3.6223E02_WP;
    pc%liqHeatCapacityMol%C = 9.3790E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 175.47_WP, 400.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 7.1117E04_WP, 1.1097E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 3.9252E04_WP;
    pc%IG_HeatCapacityMol%B = 8.7900E04_WP;
    pc%IG_HeatCapacityMol%C = 1.9165E03_WP;
    pc%IG_HeatCapacityMol%D = 5.3654E03_WP;
    pc%IG_HeatCapacityMol%E = 8.9670E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 3.9802E04_WP, 1.0533E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.3000E-02_WP;
    pc%SecondViralCoefMol%B = -1.4400E01_WP;
    pc%SecondViralCoefMol%C = -2.0600E07_WP;
    pc%SecondViralCoefMol%D = 7.6000E19_WP;
    pc%SecondViralCoefMol%E =-4.6200E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 298.15_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.0770E00_WP, -2.7019E-03_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.5317E01_WP;
    pc%liqViscosity%B = 1.7892E03_WP;
    pc%liqViscosity%C = 2.0690E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 175.47_WP, 337.85_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.1928E-02_WP, 3.4423E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.0663E-07_WP;
    pc%vapViscosity%B = 6.9655E-01_WP;
    pc%vapViscosity%C = 2.0500E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 240.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 7.5232E-06_WP, 3.1281E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.8370E-01_WP;
    pc%liqThermalConductivity%B = -2.8100E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 175.47_WP, 337.85_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.3439E-01_WP, 1.8876E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -7.7630E00_WP;
    pc%vapThermalConductivity%B = 1.0279E00_WP;
    pc%vapThermalConductivity%C = -7.4360E07_WP;
    pc%vapThermalConductivity%D = 6.7700E09_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.00_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.3652E-02_WP, 1.3927E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 3.5130E-02_WP;
    pc%SurfaceTension%B = -7.0400E-06_WP;
    pc%SurfaceTension%C = -1.2160E-07_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 273.10_WP, 503.15_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.4138E-02_WP, 8.0358E-04_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_methanol


  subroutine pc_methyltertbutylether(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'methyltertbutylether';
    pc%CAS = '1634-04-4';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 88.150_WP;    ! kg/kmol    
    pc%Tcrit = 497.10_WP;              ! K    
    pc%pcrit = 3.4300E6_WP;            ! Pa 
    pc%vcrit = 0.32900_WP;             ! m^3/kmol 
    pc%CritCompressFactor =0.273_WP;  ! -


    pc%MeltingPoint = 164.55_WP;             ! K
    pc%TripletPointT = 164.55_WP;            ! K
    pc%TripletPointp = 5.3566E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 328.35_WP;       ! K
    pc%liqMolarVolume  = 0.119887_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.8350E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.1750E8_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy =  3.5780E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  7.6000E06_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.1049E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.2661_WP;            ! -
    pc%RadiusOfGyration = 3.1790E-10_WP;      ! m
    pc%SolubilityParameter = 1.5070E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 4.5400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.06171_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 9.0800E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3663_WP;           ! -
    pc%FlashPoint = 245.00_WP;                ! K
    pc%FlammabilityLimits = (/ 2.0_WP, 15.1_WP /);    ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 733.15_WP;          ! K


    pc%solDensityMol%A = 1.1400E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 164.55_WP, 90.67_WP /);
    pc%solDensityMol%range(:,2) = (/ 1.1400E01_WP, 3.0309E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.2157E-01_WP;
    pc%liqDensityMol%B = 2.7032E-01_WP;
    pc%liqDensityMol%C = 4.9710E02_WP;
    pc%liqDensityMol%D = 2.8290E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 164.55_WP, 497.10_WP /);
    pc%liqDensityMol%range(:,2) = (/ 9.7682E00_WP, 3.0392E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.5875E01_WP;
    pc%vapPressure%B = -5.1316E03_WP;
    pc%vapPressure%C = -4.9604E00_WP;
    pc%vapPressure%D =  1.9123E-17_WP;
    pc%vapPressure%E =  6.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 164.55_WP, 497.10_WP /);
    pc%vapPressure%range(:,2) = (/ 5.3566E-01_WP, 3.4106E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.2024E07_WP;
    pc%HeatOfVapMol%B = 3.7826E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 164.55_WP, 497.10_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.6096E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.2579E04_WP;
    pc%solHeatCapacityMol%B = 1.4658E03_WP;
    pc%solHeatCapacityMol%C = -7.7409E00_WP;
    pc%solHeatCapacityMol%D = 9.6626E-03_WP;
    pc%solHeatCapacityMol%E = 9.9006E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 11.50_WP, 164.55_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.2704E03_WP, 1.3466E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.4012E05_WP;
    pc%liqHeatCapacityMol%B = -9.0000E00_WP;
    pc%liqHeatCapacityMol%C = 5.6300E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 164.55_WP, 328.35_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.5388E05_WP, 1.9786E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 9.9328E04_WP;
    pc%IG_HeatCapacityMol%B = 3.0667E05_WP;
    pc%IG_HeatCapacityMol%C = 1.7426E03_WP;
    pc%IG_HeatCapacityMol%D = 2.0764E05_WP;
    pc%IG_HeatCapacityMol%E = 7.9559E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.0394E05_WP, 3.4321E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.8083E-01_WP;
    pc%SecondViralCoefMol%B = -1.8999E02_WP;
    pc%SecondViralCoefMol%C = -2.2566E07_WP;
    pc%SecondViralCoefMol%D = -9.3110E18_WP;
    pc%SecondViralCoefMol%E = 9.3054E18_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 248.60_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.8896E00_WP, 4.6565E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -8.4217E00_WP;
    pc%liqViscosity%B = 8.5730E02_WP;
    pc%liqViscosity%C = -4.3178E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 180.00_WP, 449.93_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.7363E-03_WP, 1.0579E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.5439E-07_WP;
    pc%vapViscosity%B = 7.3601E-01_WP;
    pc%vapViscosity%C = 1.0824E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 164.55_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.9839E-06_WP, 2.2491E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2555E-01_WP;
    pc%liqThermalConductivity%B = -3.2780E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 164.55_WP, 328.35_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7161E-01_WP, 1.1792E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 3.1710E-04_WP;
    pc%vapThermalConductivity%B = 8.8303E-01_WP;
    pc%vapThermalConductivity%C = 4.6940E02_WP;
    pc%vapThermalConductivity%D = 8.6314E04_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 273.00_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.1584E-02_WP, 9.0857E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.9810E-02_WP;
    pc%SurfaceTension%B = 1.2352E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 164.55_WP, 497.10_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.6402E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_methyltertbutylether


  subroutine pc_nbutane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'nbutane';
    pc%CAS = '106-97-8';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 58.123_WP;    ! kg/kmol    
    pc%Tcrit = 425.12_WP;              ! K    
    pc%pcrit = 3.7960E6_WP;            ! Pa 
    pc%vcrit = 0.25500_WP;             ! m^3/kmol 
    pc%CritCompressFactor =0.274_WP;   ! -


    pc%MeltingPoint = 134.86_WP;             ! K
    pc%TripletPointT = 134.86_WP;            ! K
    pc%TripletPointp = 6.7358E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 272.65_WP;       ! K
    pc%liqMolarVolume  = 0.096484_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -1.2579E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -1.6700E07_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy =  3.0991E5_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint =  4.6610E06_WP;         ! J/kmol
    pc%StandardNetHeatOfCombustion = - 2.6573E9_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.2002_WP;            ! -
    pc%RadiusOfGyration = 2.8860E-10_WP;      ! m
    pc%SolubilityParameter = 1.3700E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.04780_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 6.9400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3292_WP;           ! -
    pc%FlashPoint = MAX_REAL_WP;                   ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 9.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 561.00_WP;          ! K


    pc%solDensityMol%A = 1.5194E01_WP;
    pc%solDensityMol%B = -6.8446E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 23.15_WP, 134.86_WP /);
    pc%solDensityMol%range(:,2) = (/ 1.5036E01_WP, 1.4271E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.0677E00_WP;
    pc%liqDensityMol%B = 2.7188E-01_WP;
    pc%liqDensityMol%C = 4.2512E02_WP;
    pc%liqDensityMol%D = 2.8688E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 134.86_WP, 425.12_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.2620E01_WP, 3.9271E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  6.6343E01_WP;
    pc%vapPressure%B = -4.3632E03_WP;
    pc%vapPressure%C = -7.0460E00_WP;
    pc%vapPressure%D =  9.4509E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 134.86_WP, 425.12_WP /);
    pc%vapPressure%range(:,2) = (/ 6.7441E-01_WP, 3.7699E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 3.6238E07_WP;
    pc%HeatOfVapMol%B = 8.3370E-01_WP;
    pc%HeatOfVapMol%C = -8.2274E-01_WP;
    pc%HeatOfVapMol%D = 3.9613E-01_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 134.86_WP, 425.12_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 2.8684E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.8858E04_WP;
    pc%solHeatCapacityMol%B = 1.5065E03_WP;
    pc%solHeatCapacityMol%C = -9.8874E00_WP;
    pc%solHeatCapacityMol%D = 3.2073E-02_WP;
    pc%solHeatCapacityMol%E = 1.6746E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 22.35_WP, 103.54_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0235E04_WP, 6.8653E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 6.4730E01_WP;
    pc%liqHeatCapacityMol%B = 1.6184E05_WP;
    pc%liqHeatCapacityMol%C = 9.8341E02_WP;
    pc%liqHeatCapacityMol%D = -1.4315E03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 134.86_WP, 420.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.1380E-05_WP, 5.0822E05_WP /);
    pc%liqHeatCapacityMol%eqn = 114;    

    pc%IG_HeatCapacityMol%A = 7.1340E04_WP;
    pc%IG_HeatCapacityMol%B = 2.4300E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6300E03_WP;
    pc%IG_HeatCapacityMol%D = 1.5033E05_WP;
    pc%IG_HeatCapacityMol%E = 7.3042E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 7.6733E04_WP, 2.6602E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.4900E-01_WP;
    pc%SecondViralCoefMol%B = -1.4300E02_WP;
    pc%SecondViralCoefMol%C = -1.0010E07_WP;
    pc%SecondViralCoefMol%D = -3.0300E18_WP;
    pc%SecondViralCoefMol%E = 4.0300E20_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 212.58_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.8377E00_WP, 5.0701E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -7.2471E00_WP;
    pc%liqViscosity%B = 5.3482E02_WP;
    pc%liqViscosity%C = -5.7469E-01_WP;
    pc%liqViscosity%D = -4.6625E-27_WP;
    pc%liqViscosity%E = 1E01_WP;
    pc%liqViscosity%range(:,1) = (/ 134.86_WP, 420.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.2434E-03_WP, 3.5664E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.2980E-07_WP;
    pc%vapViscosity%B = 6.9442E-01_WP;
    pc%vapViscosity%C = 2.2766E02_WP;
    pc%vapViscosity%D = -1.4610E04_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 134.86_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.6738E-06_WP, 2.2947E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.7349E-01_WP;
    pc%liqThermalConductivity%B = -7.1267E-04_WP;
    pc%liqThermalConductivity%C = 5.1555E-07_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 134.86_WP, 400.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.8676E-01_WP, 7.0910E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.1094E-02_WP;
    pc%vapThermalConductivity%B = 4.5253E-01_WP;
    pc%vapThermalConductivity%C = 5.4555E03_WP;
    pc%vapThermalConductivity%D = 1.9798E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 272.65_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.3570E-02_WP, 1.3799E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.1960E-02_WP;
    pc%SurfaceTension%B = 1.2181E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 134.86_WP, 425.12_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.2644E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_nbutane


  subroutine pc_nbutylformate(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'nbutylformate';
    pc%CAS = '592-84-7';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 2;

    pc%MolecularWeight = 102.133_WP;   ! kg/kmol    
    pc%Tcrit = 559.0_WP;               ! K    
    pc%pcrit = 3.5100E6_WP;            ! Pa 
    pc%vcrit = 0.336_WP;               ! m^3/kmol 
    pc%CritCompressFactor =0.254_WP;   ! -


    pc%MeltingPoint = 181.25_WP;             ! K
    pc%TripletPointT = 181.25_WP;            ! K
    pc%TripletPointp = 0.043691_WP;          ! Pa
    pc%NormalBoilingPoint = 379.25_WP;       ! K
    pc%liqMolarVolume  = 0.115158_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -4.2710E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = -2.8440E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 4.0800E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.3900E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.7099E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = -2.9390E08_WP;      ! J/kmol
    pc%StandardHeatOfFormation = -4.6920E08_WP;       ! J/kmol
    pc%StandardAbsoluteEntropy = 2.9860E05_WP;        ! J/kmolK    


    pc%AcentricFactor = 0.387793_WP;          ! -
    pc%RadiusOfGyration = 3.7850E-10_WP;      ! m
    pc%SolubilityParameter = 1.8230E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 6.7700E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.0632_WP;         ! m^3/kmol


    pc%VanDerWaalsArea = 9.1400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3874_WP;           ! -
    pc%FlashPoint = 291.0_WP;                 ! K
    pc%FlammabilityLimits = (/ 1.7_WP, 8.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 284.0_WP, 313.0_WP /); ! K
    pc%AutoignitionTemp = 595.0_WP;           ! K


    pc%solDensityMol%A = 1.1400E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 181.25,1.1400E01 /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 7.9185E-01_WP;
    pc%liqDensityMol%B = 2.6606E-01_WP;
    pc%liqDensityMol%C = 5.5900E02_WP;
    pc%liqDensityMol%D = 2.7850E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 181.25_WP, 559.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 9.7549E00_WP, 2.9762E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.1554E02_WP;
    pc%vapPressure%B = -7.9070E03_WP;
    pc%vapPressure%C = -1.4945E01_WP;
    pc%vapPressure%D =  1.4711E-02_WP;
    pc%vapPressure%E =  1.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 181.25_WP, 559.00_WP /);
    pc%vapPressure%range(:,2) = (/ 4.3690E-02_WP, 3.5210E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.5800E07_WP;
    pc%HeatOfVapMol%B = 4.1300E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 181.25_WP, 559.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.7461E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 0.0_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.0210E05_WP;
    pc%liqHeatCapacityMol%B = 3.2795E02_WP;
    pc%liqHeatCapacityMol%C = 0.0_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 290.15_WP, 379.25_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.9725E05_WP, 2.2648E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0596E05_WP;
    pc%IG_HeatCapacityMol%B = 2.9480E05_WP;
    pc%IG_HeatCapacityMol%C = 1.8990E03_WP;
    pc%IG_HeatCapacityMol%D = 2.2542E05_WP;
    pc%IG_HeatCapacityMol%E = -8.2219E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 298.15_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.3348E05_WP, 3.3415E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.0700E-01_WP;
    pc%SecondViralCoefMol%B = -2.3200E02_WP;
    pc%SecondViralCoefMol%C = -4.0100E07_WP;
    pc%SecondViralCoefMol%D = -1.7314E20_WP;
    pc%SecondViralCoefMol%E = 2.9240E20_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 279.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -4.2995E00_WP, 4.0446E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.7801E01_WP;
    pc%liqViscosity%B = 1.4978E03_WP;
    pc%liqViscosity%C = 9.5330E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 250.00_WP, 379.25_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.4356E-03_WP, 2.7721E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 6.2067E-07_WP;
    pc%vapViscosity%B = 5.7800E-01_WP;
    pc%vapViscosity%C = 3.6965E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 181.25_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.1243E-06_WP, 2.4561E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1050E-01_WP;
    pc%liqThermalConductivity%B = -2.3200E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 181.25_WP, 503.15_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6845E-01_WP, 9.3769E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 1.3460E03_WP;
    pc%vapThermalConductivity%B = 9.9860E-01_WP;
    pc%vapThermalConductivity%C = 1.2465E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 379.25_WP, 999.25_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.5403E-02_WP, 1.0678E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.7032E-02_WP;
    pc%SurfaceTension%B = 1.1073E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 181.25_WP, 559.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.6953E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_nbutylformate


  subroutine pc_ndecane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'ndecane';
    pc%CAS = '124-18-5';
    ! components
    pc%ChemicalFormula%C = 10;         
    pc%ChemicalFormula%H = 22;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 142.285_WP;   ! kg/kmol    
    pc%Tcrit = 617.70_WP;              ! K    
    pc%pcrit = 2.1100E6_WP;            ! Pa 
    pc%vcrit = 0.60000_WP;             ! m^3/kmol 
    pc%CritCompressFactor =0.247_WP;   ! -


    pc%MeltingPoint = 243.51_WP;             ! K
    pc%TripletPointT = 243.51_WP;            ! K
    pc%TripletPointp = 1.3930_WP;            ! Pa
    pc%NormalBoilingPoint = 447.30_WP;       ! K
    pc%liqMolarVolume  = 0.195827_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.4946E8_WP;             ! J/kmol
    pc%IG_GibbsOfFormation = 3.3180E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 5.4570E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 2.8710E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -6.2942E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.4923_WP;            ! -
    pc%RadiusOfGyration = 5.1480E-10_WP;      ! m
    pc%SolubilityParameter = 1.5700E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.10920_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.5040E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4097_WP;           ! -
    pc%FlashPoint = 319.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.7_WP, 5.4_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ 284.0_WP, 313.0_WP /); ! K
    pc%AutoignitionTemp = 474.00_WP;          ! K


    pc%solDensityMol%A = 6.8978E00_WP;
    pc%solDensityMol%B = -2.4600E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 193.15_WP, 243.15_WP /);
    pc%solDensityMol%range(:,2) = (/ 6.4227_WP, 6.2988E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 4.2831E-01_WP;
    pc%liqDensityMol%B = 2.5745E-01_WP;
    pc%liqDensityMol%C = 6.1770E02_WP;
    pc%liqDensityMol%D = 2.8912E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 243.15_WP, 617.70_WP /);
    pc%liqDensityMol%range(:,2) = (/ 5.3811E00_WP, 1.6637E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.1273E02_WP;
    pc%vapPressure%B = -9.7496E03_WP;
    pc%vapPressure%C = -1.3245E01_WP;
    pc%vapPressure%D =  7.1266E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 243.51_WP, 617.70_WP /);
    pc%vapPressure%range(:,2) = (/ 1.3930E00_WP, 2.0908E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.6126E07_WP;
    pc%HeatOfVapMol%B = 3.9797E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 243.51_WP, 617.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.4168E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -3.4900E04_WP;
    pc%solHeatCapacityMol%B = 2.5280E03_WP;
    pc%solHeatCapacityMol%C = -1.2650E01_WP;
    pc%solHeatCapacityMol%D = 2.8480E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 19.93_WP, 240.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0684E04_WP, 2.3689E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.7862E05_WP;
    pc%liqHeatCapacityMol%B = -1.9791E02_WP;
    pc%liqHeatCapacityMol%C = 1.0737E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 243.51_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.9409E05_WP, 4.1478E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.6720E05_WP;
    pc%IG_HeatCapacityMol%B = 5.3530E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6141E03_WP;
    pc%IG_HeatCapacityMol%D = 3.7820E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4200E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.7967E05_WP, 6.0932E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.800E-01_WP;
    pc%SecondViralCoefMol%B = -5.8850E02_WP;
    pc%SecondViralCoefMol%C = -7.6060E07_WP;
    pc%SecondViralCoefMol%D = -1.6130E21_WP;
    pc%SecondViralCoefMol%E = 3.9940E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 309.23_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -7.8395E00_WP, 6.5078E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.6468E01_WP;
    pc%liqViscosity%B = 1.5335E03_WP;
    pc%liqViscosity%C = 7.5110E-01_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 243.51_WP, 448.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.3742E-03_WP, 2.1165E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.6400E-08_WP;
    pc%vapViscosity%B = 9.4870E-01_WP;
    pc%vapViscosity%C = 7.1000E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 243.51_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.7547E-06_WP, 1.7295E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0630E-01_WP;
    pc%liqThermalConductivity%B = -2.5000E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 243.51_WP, 447.30_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.4542E-01_WP, 9.4475E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -6.6840E02_WP;
    pc%vapThermalConductivity%B = 9.3230E-01_WP;
    pc%vapThermalConductivity%C = -4.0710E09_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 447.30_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1731E-02_WP, 1.0286E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5435E-02_WP;
    pc%SurfaceTension%B = 1.3095E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 243.51_WP, 617.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.8756E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_ndecane


  subroutine pc_ndodecane_1(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'ndodecane';
    pc%CAS = '112-40-3';
    ! components
    pc%ChemicalFormula%C = 12;         
    pc%ChemicalFormula%H = 26;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 170.338_WP;   ! kg/kmol    
    pc%Tcrit = 658.00_WP;              ! K    
    pc%pcrit = 1.8200E6_WP;            ! Pa 
    pc%vcrit = 0.71600_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.238_WP;  ! -


    pc%MeltingPoint = 263.57_WP;             ! K
    pc%TripletPointT = 263.57_WP;            ! K
    pc%TripletPointp = 6.1520E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 489.47_WP;       ! K
    pc%liqMolarVolume  = 0.228605_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.9072E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 3.3180E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 6.2415E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 3.6840E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -7.5137E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.5764_WP;            ! -
    pc%RadiusOfGyration = 5.9140E-10_WP;      ! m
    pc%SolubilityParameter = 1.5900E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.12960_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.7740E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4151_WP;           ! -
    pc%FlashPoint = 347.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.6_WP, 4.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 476.00_WP;          ! K


    pc%solDensityMol%A = 5.2800E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 263.57_WP, 5.2800E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.5541E-01_WP;
    pc%liqDensityMol%B = 2.5511E-01_WP;
    pc%liqDensityMol%C = 6.5800E02_WP;
    pc%liqDensityMol%D = 2.9368E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 4.5132E00_WP, 1.3932E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.3747E02_WP;
    pc%vapPressure%B = -1.1976E04_WP;
    pc%vapPressure%C = -1.6698E01_WP;
    pc%vapPressure%D =  8.0906E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%vapPressure%range(:,2) = (/ 6.1534E-01_WP, 1.8223E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 7.7337E07_WP;
    pc%HeatOfVapMol%B = 4.0681E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.2802E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -4.0970E04_WP;
    pc%solHeatCapacityMol%B = 2.8390E03_WP;
    pc%solHeatCapacityMol%C = -1.3232E01_WP;
    pc%solHeatCapacityMol%D = 2.7972E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 29.97_WP, 258.41_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.2983E04_WP, 2.9175E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 5.0821E05_WP;
    pc%liqHeatCapacityMol%B = -1.3687E03_WP;
    pc%liqHeatCapacityMol%C = 3.1015E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 263.57_WP, 330.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.6292E05_WP, 3.9429E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.1295E05_WP;
    pc%IG_HeatCapacityMol%B = 6.6330E05_WP;
    pc%IG_HeatCapacityMol%C = 1.7155E03_WP;
    pc%IG_HeatCapacityMol%D = 4.5161E05_WP;
    pc%IG_HeatCapacityMol%E = 7.7750E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.2442E05_WP, 7.4325E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 8.800E-01_WP;
    pc%SecondViralCoefMol%B = -1.0910E03_WP;
    pc%SecondViralCoefMol%C = -5.0300E07_WP;
    pc%SecondViralCoefMol%D = -5.4871E21_WP;
    pc%SecondViralCoefMol%E = 1.4959E24_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 329.00_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.0701E01_WP, 1.3759E-01_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0607E01_WP;
    pc%liqViscosity%B = 1.9430E03_WP;
    pc%liqViscosity%C = 1.3205E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 263.57_WP, 489.47_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.8109E-03_WP, 2.1196E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 6.3440E-08_WP;
    pc%vapViscosity%B = 8.2870E-01_WP;
    pc%vapViscosity%C = 2.1950E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 263.57_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.5111E-06_WP, 1.5932E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0470E-01_WP;
    pc%liqThermalConductivity%B = -2.3260E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 263.57_WP, 489.47_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.4339E-01_WP, 9.0849E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.7190E-06_WP;
    pc%vapThermalConductivity%B = 1.4699E00_WP;
    pc%vapThermalConductivity%C = 5.7940E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 489.47_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.3537E-02_WP, 9.3010E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5493E-02_WP;
    pc%SurfaceTension%B = 1.3262E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 263.57_WP, 658.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.8150E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_ndodecane_1


  subroutine pc_nheptane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'nheptane';
    pc%CAS = '142-82-5';
    ! components
    pc%ChemicalFormula%C = 7;         
    pc%ChemicalFormula%H = 16;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 100.204_WP;   ! kg/kmol    
    pc%Tcrit = 540.20_WP;              ! K    
    pc%pcrit = 2.7400E6_WP;            ! Pa 
    pc%vcrit = 0.42800_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.261_WP;  ! -


    pc%MeltingPoint = 182.57_WP;             ! K
    pc%TripletPointT = 182.57_WP;            ! K
    pc%TripletPointp = 1.8269E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 371.58_WP;       ! K
    pc%liqMolarVolume  = 0.147024_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -1.8765E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 8.1650E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 4.2798E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.4050E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -4.4647E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.3495_WP;            ! -
    pc%RadiusOfGyration = 4.1730E-10_WP;      ! m
    pc%SolubilityParameter = 1.5200E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.07849_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.0990E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3851_WP;           ! -
    pc%FlashPoint = 269.00_WP;                ! K
    pc%FlammabilityLimits = (/ 1.0_WP, 7.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 477.00_WP;          ! K


    pc%solDensityMol%A = 9.1435E00_WP;
    pc%solDensityMol%B = -2.8756E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 133.15_WP, 182.57_WP /);
    pc%solDensityMol%range(:,2) = (/ 8.7606E00_WP, 8.6185E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 6.1259E-01_WP;
    pc%liqDensityMol%B = 2.6211E-01_WP;
    pc%liqDensityMol%C = 5.4020E02_WP;
    pc%liqDensityMol%D = 2.8141E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 182.57_WP, 540.20_WP /);
    pc%liqDensityMol%range(:,2) = (/ 7.6998E00_WP, 2.3371E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  8.7829E01_WP;
    pc%vapPressure%B = -6.9964E03_WP;
    pc%vapPressure%C = -9.8802E00_WP;
    pc%vapPressure%D =  7.2099E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 182.57_WP, 540.20_WP /);
    pc%vapPressure%range(:,2) = (/ 1.8269E-01_WP, 2.7192E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.0014E07_WP;
    pc%HeatOfVapMol%B = 3.8795E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 182.57_WP, 540.20_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.2619E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -2.3900E04_WP;
    pc%solHeatCapacityMol%B = 1.9408E03_WP;
    pc%solHeatCapacityMol%C = -1.0480E01_WP;
    pc%solHeatCapacityMol%D = 2.7190E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 182.58_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0942E04_WP, 1.4658E+05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 6.1260E01_WP;
    pc%liqHeatCapacityMol%B = 3.1441E05_WP;
    pc%liqHeatCapacityMol%C = 1.8246E03_WP;
    pc%liqHeatCapacityMol%D = -2.5479E03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 182.57_WP, 520.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.9989E05_WP, 4.0657E05_WP /);
    pc%liqHeatCapacityMol%eqn = 114;    

    pc%IG_HeatCapacityMol%A = 1.2015E05_WP;
    pc%IG_HeatCapacityMol%B = 4.0010E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6766E03_WP;
    pc%IG_HeatCapacityMol%D = 2.7400E05_WP;
    pc%IG_HeatCapacityMol%E = 7.5640E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.2828E05_WP, 4.4283E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.7460E-01_WP;
    pc%SecondViralCoefMol%B = -2.9100E02_WP;
    pc%SecondViralCoefMol%C = -4.4180E07_WP;
    pc%SecondViralCoefMol%D = -8.8000E19_WP;
    pc%SecondViralCoefMol%E = 1.2850E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 270.10_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ - 4.4720E00_WP, 6.7507E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.4451E01_WP;
    pc%liqViscosity%B = 1.5331E03_WP;
    pc%liqViscosity%C = 2.0087E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 182.57_WP, 373.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.7197E-03_WP, 2.1455E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 6.6720E-08_WP;
    pc%vapViscosity%B = 8.2837E-01_WP;
    pc%vapViscosity%C = 8.5752E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 182.57_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.3910E-06_WP, 1.8777E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1500E-01_WP;
    pc%liqThermalConductivity%B = -3.0300E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 182.57_WP, 371.58_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5968E-01_WP, 1.0241E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -7.0028E-02_WP;
    pc%vapThermalConductivity%B = 3.8068E-01_WP;
    pc%vapThermalConductivity%C = -7.0499E03_WP;
    pc%vapThermalConductivity%D = -2.4005E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 339.15_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.5827E-02_WP, 1.1493E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.4143E-02_WP;
    pc%SurfaceTension%B = 1.2512E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 182.57_WP, 540.20_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.2317E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_nheptane


  subroutine pc_nhexane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'nhexane';
    pc%CAS = '110-54-3';
    ! components
    pc%ChemicalFormula%C = 6;         
    pc%ChemicalFormula%H = 14;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 86.177_WP;    ! kg/kmol    
    pc%Tcrit = 507.60_WP;              ! K    
    pc%pcrit = 3.0250E6_WP;            ! Pa 
    pc%vcrit = 0.37100_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.266_WP;  ! -


    pc%MeltingPoint = 177.83_WP;             ! K
    pc%TripletPointT = 177.83_WP;            ! K
    pc%TripletPointp = 9.0170E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 341.88_WP;       ! K
    pc%liqMolarVolume  = 0.131362_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -1.6694E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -6.6340E04_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 3.8874E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.3080E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.8551E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.3013_WP;            ! -
    pc%RadiusOfGyration = 3.7690E-10_WP;      ! m
    pc%SolubilityParameter = 1.4900E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.06826_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 9.6400E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3723_WP;           ! -
    pc%FlashPoint = 251.50_WP;                ! K
    pc%FlammabilityLimits = (/ 1.1_WP, 7.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 498.00_WP;          ! K


    pc%solDensityMol%A = 1.0657E01_WP;
    pc%solDensityMol%B = -4.6530E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 133.15_WP, 177.83_WP /);
    pc%solDensityMol%range(:,2) = (/ 1.0037E01_WP, 9.8296E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 7.0824E-01_WP;
    pc%liqDensityMol%B = 2.6411E-01_WP;
    pc%liqDensityMol%C = 5.0760E02_WP;
    pc%liqDensityMol%D = 2.7537E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 177.83_WP, 507.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 8.7470E00_WP, 2.6816E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.0465E02_WP;
    pc%vapPressure%B = -6.9955E03_WP;
    pc%vapPressure%C = -1.2702E01_WP;
    pc%vapPressure%D =  1.2381E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 177.83_WP, 507.60_WP /);
    pc%vapPressure%range(:,2) = (/ 9.0169E-01_WP, 3.0449E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.4544E07_WP;
    pc%HeatOfVapMol%B = 3.9002E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 177.83_WP, 507.60_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.7647E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.9370E04_WP;
    pc%solHeatCapacityMol%B = 1.6560E03_WP;
    pc%solHeatCapacityMol%C = -8.3720E00_WP;
    pc%solHeatCapacityMol%D = 2.0200E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 177.86_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0563E04_WP, 1.2398E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.7212E05_WP;
    pc%liqHeatCapacityMol%B = -1.8378E02_WP;
    pc%liqHeatCapacityMol%C = 8.8734E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 177.83_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.6750E05_WP, 2.7534E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.0440E05_WP;
    pc%IG_HeatCapacityMol%B = 3.5230E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6946E03_WP;
    pc%IG_HeatCapacityMol%D = 2.3690E05_WP;
    pc%IG_HeatCapacityMol%E = 7.6160E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.1117E05_WP, 3.8620E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.9700E-01_WP;
    pc%SecondViralCoefMol%B = -2.0000E02_WP;
    pc%SecondViralCoefMol%C = -3.7100E07_WP;
    pc%SecondViralCoefMol%D = 1.6700E19_WP;
    pc%SecondViralCoefMol%E = -6.1800E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 253.70_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -3.3097E00_WP, 5.2675E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0715E01_WP;
    pc%liqViscosity%B = 1.2075E03_WP;
    pc%liqViscosity%C = 1.4993E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 177.83_WP, 343.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.1182E-03_WP, 2.1541E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.7514E-07_WP;
    pc%vapViscosity%B = 7.0737E-01_WP;
    pc%vapViscosity%C = 1.5714E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 177.83_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.6305E-06_WP, 2.0050E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.2492E-01_WP;
    pc%liqThermalConductivity%B = -3.5300E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 177.83_WP, 370.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.6209E-01_WP, 9.4199E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -6.5050E02_WP;
    pc%vapThermalConductivity%B = 8.0530E-01_WP;
    pc%vapThermalConductivity%C = -1.4121E09_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 339.09_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.7036E-02_WP, 1.2003E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5003E-02_WP;
    pc%SurfaceTension%B = 1.2674E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 177.83_WP, 507.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.1841E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_nhexane


  subroutine pc_nnonane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'nnonane';
    pc%CAS = '111-84-2';
    ! components
    pc%ChemicalFormula%C = 9;         
    pc%ChemicalFormula%H = 20;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 128.258_WP;   ! kg/kmol    
    pc%Tcrit = 594.60_WP;              ! K    
    pc%pcrit = 2.2900E6_WP;            ! Pa 
    pc%vcrit = 0.54400_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.252_WP;  ! -


    pc%MeltingPoint = 219.66_WP;             ! K
    pc%TripletPointT = 219.66_WP;            ! K
    pc%TripletPointp = 4.3058E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 423.97_WP;       ! K
    pc%liqMolarVolume  = 0.179559_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.2874E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 2.4980E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 5.0640E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 1.5470E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -5.6846E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.4435_WP;            ! -
    pc%RadiusOfGyration = 4.8140E-10_WP;      ! m
    pc%SolubilityParameter = 1.5600E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.09895_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.3690E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4031_WP;           ! -
    pc%FlashPoint = 304.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.7_WP, 5.6_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 478.00_WP;          ! K


    pc%solDensityMol%A = 7.4283E00_WP;
    pc%solDensityMol%B = -2.3699E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 193.15, 219.66 /);
    pc%solDensityMol%range(:,2) = (/ 6.9706E00, 6.9077E00 /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 4.8387E-01_WP;
    pc%liqDensityMol%B = 2.6147E-01_WP;
    pc%liqDensityMol%C = 5.9460E02_WP;
    pc%liqDensityMol%D = 2.8281E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 219.66_WP, 594.60_WP /);
    pc%liqDensityMol%range(:,2) = (/ 6.0070E00_WP, 1.8506E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.0935E02_WP;
    pc%vapPressure%B = -9.0304E03_WP;
    pc%vapPressure%C = -1.2882E01_WP;
    pc%vapPressure%D =  7.8544E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 219.66_WP, 594.60_WP /);
    pc%vapPressure%range(:,2) = (/ 4.3058E-01_WP, 2.3054E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.0370E07_WP;
    pc%HeatOfVapMol%B = 3.8522E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 219.66_WP, 594.60_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.0545E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -3.1900E04_WP;
    pc%solHeatCapacityMol%B = 2.3720E03_WP;
    pc%solHeatCapacityMol%C = -1.2440E01_WP;
    pc%solHeatCapacityMol%D = 3.0180E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 217.21_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.0805E04_WP, 2.0569E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.8308E05_WP;
    pc%liqHeatCapacityMol%B = -1.1398E03_WP;
    pc%liqHeatCapacityMol%C = 2.7101E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 219.66_WP, 325.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.6346E05_WP, 2.9890E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.5175E05_WP;
    pc%IG_HeatCapacityMol%B = 4.9150E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6448E03_WP;
    pc%IG_HeatCapacityMol%D = 3.4700E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4960E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.6257E05_WP, 5.5407E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 3.9060E-01_WP;
    pc%SecondViralCoefMol%B = -4.6860E02_WP;
    pc%SecondViralCoefMol%C = -6.6800E07_WP;
    pc%SecondViralCoefMol%D = -7.8930E20_WP;
    pc%SecondViralCoefMol%E = 1.8200E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 297.83_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -6.5900E00_WP, 5.8381E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.1149E01_WP;
    pc%liqViscosity%B = 1.6580E03_WP;
    pc%liqViscosity%C = 1.4540E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 219.66_WP, 423.97_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.1483E-03_WP, 2.1559E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.0344E-07_WP;
    pc%vapViscosity%B = 7.7301E-01_WP;
    pc%vapViscosity%C = 2.2047E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 219.66_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.3347E-06_WP, 1.7668E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0900E-01_WP;
    pc%liqThermalConductivity%B = -2.6400E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 219.66_WP, 423.97_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5101E-01_WP, 9.7072E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -6.5771E-02_WP;
    pc%vapThermalConductivity%B = 2.7198E-01_WP;
    pc%vapThermalConductivity%C = -3.4823E03_WP;
    pc%vapThermalConductivity%D = -1.5803E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 423.97_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.1299E-02_WP, 1.0597E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.4975E-02_WP;
    pc%SurfaceTension%B = 1.2897E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 219.66_WP, 594.60_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.0331E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_nnonane


  subroutine pc_noctane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'noctane';
    pc%CAS = '111-65-9';
    ! components
    pc%ChemicalFormula%C = 8;         
    pc%ChemicalFormula%H = 18;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 114.231_WP;   ! kg/kmol    
    pc%Tcrit = 568.70_WP;              ! K    
    pc%pcrit = 2.4900E6_WP;            ! Pa 
    pc%vcrit = 0.48600_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.256_WP;  ! -


    pc%MeltingPoint = 216.38_WP;             ! K
    pc%TripletPointT = 216.38_WP;            ! K
    pc%TripletPointp = 2.1083E00_WP;         ! Pa
    pc%NormalBoilingPoint = 398.83_WP;       ! K
    pc%liqMolarVolume  = 0.163374_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -2.0875E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 1.6000E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 4.6723E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 2.0740E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -5.0741E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.3996_WP;            ! -
    pc%RadiusOfGyration = 4.5460E-10_WP;      ! m
    pc%SolubilityParameter = 1.5400E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.08872_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.3690E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3951_WP;           ! -
    pc%FlashPoint = 286.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.8_WP, 6.5_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 479.00_WP;          ! K


    pc%solDensityMol%A = 8.3409E00_WP;
    pc%solDensityMol%B = -3.1515E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 133.15_WP, 216.38_WP /);
    pc%solDensityMol%range(:,2) = (/ 7.9213E00_WP, 7.6590E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 5.3731E-01_WP;
    pc%liqDensityMol%B = 2.6115E-01_WP;
    pc%liqDensityMol%C = 5.6870E02_WP;
    pc%liqDensityMol%D = 2.8034E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 216.38, 568.70 /);
    pc%liqDensityMol%range(:,2) = (/ 6.6558E00, 2.0575E00 /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  9.6084E01_WP;
    pc%vapPressure%B = -7.9002E03_WP;
    pc%vapPressure%C = -1.1003E01_WP;
    pc%vapPressure%D =  7.1802E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 216.38_WP, 568.70_WP /);
    pc%vapPressure%range(:,2) = (/ 2.1083E00_WP, 2.4673E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 5.5180E07_WP;
    pc%HeatOfVapMol%B = 3.8467E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 216.38_WP, 568.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 4.5898E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -2.4000E04_WP;
    pc%solHeatCapacityMol%B = 1.9472E03_WP;
    pc%solHeatCapacityMol%C = -8.5360E00_WP;
    pc%solHeatCapacityMol%D = 1.3400E-02_WP;
    pc%solHeatCapacityMol%E = 2.0940E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 216.37_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.1640E04_WP, 1.7933E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 2.2483E05_WP;
    pc%liqHeatCapacityMol%B = -1.8663E02_WP;
    pc%liqHeatCapacityMol%C = 9.5891E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 216.38_WP, 460.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 2.2934E05_WP, 3.4189E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 1.3554E05_WP;
    pc%IG_HeatCapacityMol%B = 4.4310E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6356E03_WP;
    pc%IG_HeatCapacityMol%D = 3.0540E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4640E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 1.4529E05_WP, 4.9764E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 2.7390E-01_WP;
    pc%SecondViralCoefMol%B = -3.2144E02_WP;
    pc%SecondViralCoefMol%C = -6.6830E07_WP;
    pc%SecondViralCoefMol%D = -1.2710E20_WP;
    pc%SecondViralCoefMol%E = 1.6103E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 284.38_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.4098E00_WP, 3.9801E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0463E01_WP;
    pc%liqViscosity%B = 1.4974E03_WP;
    pc%liqViscosity%C = 1.3790E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 216.38_WP, 398.83_WP /);
    pc%liqViscosity%range(:,2) = (/ 2.1813E-03_WP, 2.1383E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.1191E-08_WP;
    pc%vapViscosity%B = 9.2925E-01_WP;
    pc%vapViscosity%C = 5.5092E01_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 216.38_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.6772E-06_WP, 1.8134E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1560E-01_WP;
    pc%liqThermalConductivity%B = -2.9483E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 216.38_WP, 398.83_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5180E-01_WP, 9.8013E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -8.7580E03_WP;
    pc%vapThermalConductivity%B = 8.4480E-01_WP;
    pc%vapThermalConductivity%C = -2.7121E10_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 339.00_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.5025E-02_WP, 1.1053E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.2789E-02_WP;
    pc%SurfaceTension%B = 1.2323E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 216.38_WP, 568.70_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.9261E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_noctane


  subroutine pc_npentane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

    ! all alternative names
    pc%name = 'npentane';
    pc%CAS = '109-66-0';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 12;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 72.150_WP;    ! kg/kmol    
    pc%Tcrit = 469.70_WP;              ! K    
    pc%pcrit = 3.3700E6_WP;            ! Pa 
    pc%vcrit = 0.31300_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.270_WP;  ! -


    pc%MeltingPoint = 143.42_WP;             ! K
    pc%TripletPointT = 143.42_WP;            ! K
    pc%TripletPointp = 6.8642E-02_WP;        ! Pa
    pc%NormalBoilingPoint = 309.22_WP;       ! K
    pc%liqMolarVolume  = 0.116045_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -1.4676E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 8.8130E06_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 3.4945E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 8.4015E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -3.2449E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.2515_WP;            ! -
    pc%RadiusOfGyration = 3.3370E-10_WP;      ! m
    pc%SolubilityParameter = 1.4400E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.05803_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 8.2900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.3547_WP;           ! -
    pc%FlashPoint = 233.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.3_WP, 8.0_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 516.00_WP;          ! K


    pc%solDensityMol%A = 1.3533E01_WP;
    pc%solDensityMol%B = -6.0608E-03_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 23.15_WP, 143.42_WP /);
    pc%solDensityMol%range(:,2) = (/ 1.3393E01_WP, 1.2664E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 8.4947E-01_WP;
    pc%liqDensityMol%B = 2.6726E-01_WP;
    pc%liqDensityMol%C = 4.6970E02_WP;
    pc%liqDensityMol%D = 2.7789E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 143.42_WP, 469.70_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.0474E01_WP, 3.1784E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  7.8741E01_WP;
    pc%vapPressure%B = -5.4203E03_WP;
    pc%vapPressure%C = -8.8253E00_WP;
    pc%vapPressure%D =  9.6171E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 143.42_WP, 469.70_WP /);
    pc%vapPressure%range(:,2) = (/ 6.8642E-02_WP, 3.3642E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 3.9109E07_WP;
    pc%HeatOfVapMol%B = 3.8681E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 143.42_WP, 469.70_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.3968E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.0050E04_WP;
    pc%solHeatCapacityMol%B = 1.0103E03_WP;
    pc%solHeatCapacityMol%C = 4.2700E00_WP;
    pc%solHeatCapacityMol%D = -1.0169E-01_WP;
    pc%solHeatCapacityMol%E = 4.2200E-04_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 12.29_WP, 134.60_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 2.8324E03_WP, 9.3832E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.5908E05_WP;
    pc%liqHeatCapacityMol%B = -2.7050E02_WP;
    pc%liqHeatCapacityMol%C = 9.9537E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 143.42_WP, 390.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.4076E05_WP, 2.0498E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 8.8050E04_WP;
    pc%IG_HeatCapacityMol%B = 3.0110E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6502E03_WP;
    pc%IG_HeatCapacityMol%D = 1.8920E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4760E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 9.4039E04_WP, 3.2927E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.9160E-01_WP;
    pc%SecondViralCoefMol%B = -1.9800E02_WP;
    pc%SecondViralCoefMol%C = -1.6400E07_WP;
    pc%SecondViralCoefMol%D = -1.7920E19_WP;
    pc%SecondViralCoefMol%E = 2.6850E21_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 234.80_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -2.6205E00_WP, 5.4740E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0383E01;
    pc%liqViscosity%B = 1.0504E03;
    pc%liqViscosity%C = 1.4874E00;
    pc%liqViscosity%D = -2.0167E-27;
    pc%liqViscosity%E = 1E01;
    pc%liqViscosity%range(:,1) = (/ 143.42_WP, 465.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.4377E-03_WP, 4.7971E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 6.3412E-08;
    pc%vapViscosity%B = 8.4758E-01;
    pc%vapViscosity%C = 4.1718E01;
    pc%vapViscosity%D = 0;
    pc%vapViscosity%E = 0;
    pc%vapViscosity%range(:,1) = (/ 143.42, 1000.00 /);
    pc%vapViscosity%range(:,2) = (/ 3.3051E-06, 2.1240E-05 /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.5370E-01;
    pc%liqThermalConductivity%B = -5.7600E-04;
    pc%liqThermalConductivity%C = 3.4400E-07;
    pc%liqThermalConductivity%D = 0;
    pc%liqThermalConductivity%E = 0;
    pc%liqThermalConductivity%range(:,1) = (/ 143.42, 445.00 /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7817E-01, 6.5501E-02 /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -6.8440E02;
    pc%vapThermalConductivity%B = 7.6400E-01;
    pc%vapThermalConductivity%C = -1.0550E09;
    pc%vapThermalConductivity%D = 0;
    pc%vapThermalConductivity%E = 0;
    pc%vapThermalConductivity%range(:,1) = (/ 273.15, 1000.00 /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.2879E-02, 1.2707E-01 /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.2020E-02;
    pc%SurfaceTension%B = 1.2041E00;
    pc%SurfaceTension%C = 0;
    pc%SurfaceTension%D = 0;
    pc%SurfaceTension%E = 0;
    pc%SurfaceTension%range(:,1) = (/ 143.42, 469.70 /);
    pc%SurfaceTension%range(:,2) = (/ 3.3546E-02, 0.0 /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_npentane

  subroutine pc_ntetradecane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'ntetradecane';
    pc%CAS = '629-59-4';
    ! components
    pc%ChemicalFormula%C = 14;         
    pc%ChemicalFormula%H = 30;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 198.392_WP;  ! kg/kmol    
    pc%Tcrit = 693.00_WP;             ! K    
    pc%pcrit = 1.5700E06_WP;          ! Pa 
    pc%vcrit = 0.83000_WP;            ! m^3/kmol 
    pc%CritCompressFactor = 0.226_WP; ! -


    pc%MeltingPoint = 279.01_WP;              ! K
    pc%TripletPointT = 279.01_WP;             ! K
    pc%TripletPointp = 2.5269E-01_WP;         ! Pa
    pc%NormalBoilingPoint = 526.73_WP;        ! K
    pc%liqMolarVolume = 0.261271_WP;          ! m^3/kmol


    pc%IG_HeatOfFormation = -3.3244E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 6.5990E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 7.0259E05_WP;             ! J/kmolK
    pc%HeatFusionAtMeltPoint = 4.5070E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -8.7328E09_WP;   ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.6430_WP;            ! -
    pc%RadiusOfGyration = 6.4270E-10_WP;      ! m
    pc%SolubilityParameter = 1.6100E04_WP;    ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.15010_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 2.0440E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4269_WP;           ! -
    pc%FlashPoint = 373.15_WP;                ! K
    pc%FlammabilityLimits = (/ 0.5_WP, 4.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 473.00_WP;          ! K


    pc%solDensityMol%A = 4.6071E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 0.0_WP, 4.6071E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.0545E-01_WP;
    pc%liqDensityMol%B = 2.5350E-01_WP;
    pc%liqDensityMol%C = 6.9300E02_WP;
    pc%liqDensityMol%D = 3.0538E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 279.01_WP, 693.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 3.8924E00_WP, 1.2049E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.4047E02_WP;
    pc%vapPressure%B = -1.3231E04_WP;
    pc%vapPressure%C = -1.6859E01_WP;
    pc%vapPressure%D =  6.5877E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 279.01_WP, 693.00_WP /);
    pc%vapPressure%range(:,2) = (/ 2.5268E-01_WP, 1.5693E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 9.0539E07_WP;
    pc%HeatOfVapMol%B = 4.4467E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 279.01_WP, 693.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 7.2002E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -4.9782E04_WP;
    pc%solHeatCapacityMol%B = 3.3010E03_WP;
    pc%solHeatCapacityMol%C = -1.5807E01_WP;
    pc%solHeatCapacityMol%D = 3.6919E-02_WP;
    pc%solHeatCapacityMol%E = -1.3660E-05_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 29.98_WP, 279.01_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 3.5958E04_WP, 3.5981E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.5314E05_WP;
    pc%liqHeatCapacityMol%B = 2.9130E01_WP;
    pc%liqHeatCapacityMol%C = 8.6116E-01_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 279.01_WP, 526.73_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 4.2831E05_WP, 6.0741E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.3082E05_WP;
    pc%IG_HeatCapacityMol%B = 7.8678E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6823E03_WP;
    pc%IG_HeatCapacityMol%D = 5.4486E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4310E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.4864E05_WP, 8.6225E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.8209E-01_WP;
    pc%SecondViralCoefMol%B = -5.0132E02_WP;
    pc%SecondViralCoefMol%C = -2.8594E08_WP;
    pc%SecondViralCoefMol%D = -2.5720E21_WP;
    pc%SecondViralCoefMol%E = 4.5523E23_WP
    pc%SecondViralCoefMol%range(:,1) = (/ 346.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.3893E01_WP, 6.3065E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.0486E01_WP;
    pc%liqViscosity%B = 2.0884E03_WP;
    pc%liqViscosity%C = 1.2852E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 279.01_WP, 528.15_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.1395E-03_WP, 2.0874E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 4.4565E-08_WP;
    pc%vapViscosity%B = 8.6840E-01_WP;
    pc%vapViscosity%C = 2.2816E02_WP;
    pc%vapViscosity%D = -4.3472E03_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 279.01_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.3635E-06_WP, 1.4672E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.02936E-01_WP;
    pc%liqThermalConductivity%B = -2.1798E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 279.01_WP, 526.73_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.4211E-01_WP, 8.8113E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -1.6362E02_WP;
    pc%vapThermalConductivity%B = 9.1930E-01_WP;
    pc%vapThermalConductivity%C = -1.0876E09_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 526.73_WP, 1000.0_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.5172E-02_WP, 8.6152E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.6276E-02_WP;
    pc%SurfaceTension%B = 1.3631E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 279.01_WP, 693.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.7883E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_ntetradecane


  subroutine pc_ntridecane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'ntridecane';
    pc%CAS = '629-50-5';
    ! components
    pc%ChemicalFormula%C = 13;         
    pc%ChemicalFormula%H = 28;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 184.365_WP;   ! kg/kmol    
    pc%Tcrit = 675.00_WP;              ! K    
    pc%pcrit = 1.6800E6_WP;            ! Pa 
    pc%vcrit = 0.77500_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.232_WP;  ! -


    pc%MeltingPoint = 267.76_WP;             ! K
    pc%TripletPointT = 267.76_WP;            ! K
    pc%TripletPointp = 2.5096E-02_WP;        ! Pa
    pc%NormalBoilingPoint = 508.62_WP;       ! K
    pc%liqMolarVolume  = 0.244631_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -3.1177E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = 5.7710E07_WP;            ! J/kmol
    pc%IG_AbsoluteEntropy = 6.6337E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 2.8500E07_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -8.1229E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.6174_WP;            ! -
    pc%RadiusOfGyration = 6.1960E-10_WP;      ! m
    pc%SolubilityParameter = 1.6100E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.13960_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 1.9090E09_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4235_WP;           ! -
    pc%FlashPoint = 352.00_WP;                ! K
    pc%FlammabilityLimits = (/ 0.6_WP, 4.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 475.00_WP;          ! K


    pc%solDensityMol%A = 5.0311E00_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 0.0_WP, 5.0311E00_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 3.2160E-01_WP;
    pc%liqDensityMol%B = 2.5040E-01_WP;
    pc%liqDensityMol%C = 6.7500E02_WP;
    pc%liqDensityMol%D = 3.0710E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 267.76_WP, 675.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 4.2035E00_WP, 1.2843E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.3745E02_WP;
    pc%vapPressure%B = -1.2549E04_WP;
    pc%vapPressure%C = -1.6543E01_WP;
    pc%vapPressure%D =  7.1275E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 267.76_WP, 675.00_WP /);
    pc%vapPressure%range(:,2) = (/ 2.5096E-01_WP, 1.6786E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 8.4339E07_WP;
    pc%HeatOfVapMol%B = 4.2570E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 267.76_WP, 675.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 6.8015E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -3.0600E04_WP;
    pc%solHeatCapacityMol%B = 2.3520E03_WP;
    pc%solHeatCapacityMol%C = -3.2440E00_WP;
    pc%solHeatCapacityMol%D = -3.5173E-02_WP;
    pc%solHeatCapacityMol%E = 1.3218E-04_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 255.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.4882E04_WP, 3.3389E05_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 3.5018E05_WP;
    pc%liqHeatCapacityMol%B = -1.0470E02_WP;
    pc%liqHeatCapacityMol%C = 1.0022E00_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 267.76_WP, 508.62_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 3.9400E05_WP, 5.5619E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.1496E05_WP;
    pc%IG_HeatCapacityMol%B = 7.3045E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6695E03_WP;
    pc%IG_HeatCapacityMol%D = 4.9998E05_WP;
    pc%IG_HeatCapacityMol%E = 7.4102E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 2.3156E05_WP, 8.0251E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 4.4720E-01_WP;
    pc%SecondViralCoefMol%B = -4.7897E02_WP;
    pc%SecondViralCoefMol%C = -2.3434E08_WP;
    pc%SecondViralCoefMol%D = -1.7856E21_WP;
    pc%SecondViralCoefMol%E = 3.0980E23_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 337.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.2222E01_WP, 5.8391E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -2.1008E01_WP;
    pc%liqViscosity%B = 2.0425E03_WP;
    pc%liqViscosity%C = 1.3692E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 267.76_WP, 508.62_WP /);
    pc%liqViscosity%range(:,2) = (/ 3.2603E-03_WP, 2.1182E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.5585E-08_WP;
    pc%vapViscosity%B = 8.9870E-01_WP;
    pc%vapViscosity%C = 1.6530E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 267.76_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 3.3441E-06_WP, 1.5168E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.0447E-01_WP;
    pc%liqThermalConductivity%B = -2.2612E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 267.76_WP, 508.62_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.4392E-01_WP, 8.9461E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 5.3701E-06_WP;
    pc%vapThermalConductivity%B = 1.4751E00_WP;
    pc%vapThermalConductivity%C = 5.9909E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 508.62_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.4219E-02_WP, 8.9415E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.5322E-02_WP;
    pc%SurfaceTension%B = 1.3257E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 267.76_WP, 675.00_WP /);
    pc%SurfaceTension%range(:,2) = (/ 2.8312E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_ntridecane


  subroutine pc_propane(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'propane';
    pc%CAS = '74-98-6';
    ! components
    pc%ChemicalFormula%C = 3;         
    pc%ChemicalFormula%H = 8;
    pc%ChemicalFormula%O = 0;

    pc%MolecularWeight = 44.096_WP;    ! kg/kmol    
    pc%Tcrit = 369.83_WP;              ! K    
    pc%pcrit = 4.2480E6_WP;            ! Pa 
    pc%vcrit = 0.20000_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.276_WP;  ! -


    pc%MeltingPoint = 85.47_WP;              ! K
    pc%TripletPointT = 85.47_WP;             ! K
    pc%TripletPointp = 1.6850E-04_WP;        ! Pa
    pc%NormalBoilingPoint = 231.11_WP;       ! K
    pc%liqMolarVolume  = 0.075700_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -1.0468E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -2.4390E07_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 2.7020E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 3.5240E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.0431E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.1523_WP;            ! -
    pc%RadiusOfGyration = 2.4310E-10_WP;      ! m
    pc%SolubilityParameter = 1.3100E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 0.0_WP;                   ! C^m
    pc%VanDerWaalsVolume = 0.03757_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 5.5900E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.2861_WP;           ! -
    pc%FlashPoint = MAX_REAL_WP;                   ! K
    pc%FlammabilityLimits = (/ 2.0_WP, 9.5_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 723.00_WP;          ! K


    pc%solDensityMol%A = 1.8861E01_WP;
    pc%solDensityMol%B = -2.0332E-02_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 23.15_WP, 85.47_WP /);
    pc%solDensityMol%range(:,2) = (/ 1.8390E01_WP, 1.7123E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.3757E00_WP;
    pc%liqDensityMol%B = 2.7453E-01_WP;
    pc%liqDensityMol%C = 3.6983E02_WP;
    pc%liqDensityMol%D = 2.9359E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 85.47_WP, 369.83_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.6583E01_WP, 5.0111E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  5.9078E01_WP;
    pc%vapPressure%B = -3.4926E03_WP;
    pc%vapPressure%C = -6.0669E00_WP;
    pc%vapPressure%D =  1.0919E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 85.47_WP, 369.83_WP /);
    pc%vapPressure%range(:,2) = (/ 1.6788E-04_WP, 4.2135E06_WP /);
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 2.9209E07_WP;
    pc%HeatOfVapMol%B = 7.8237E-01_WP;
    pc%HeatOfVapMol%C = -7.7319E-01_WP;
    pc%HeatOfVapMol%D = 3.9246E-01_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 85.47_WP, 369.83_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 2.4787E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.1230E04_WP;
    pc%solHeatCapacityMol%B = 1.0590E03_WP;
    pc%solHeatCapacityMol%C = -3.6000E00_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 30.00_WP, 84.00_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 1.7300E04_WP, 5.2324E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 6.2983E01_WP;
    pc%liqHeatCapacityMol%B = 1.1363E05_WP;
    pc%liqHeatCapacityMol%C = 6.3321E02_WP;
    pc%liqHeatCapacityMol%D = -8.7346E02_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 85.47_WP, 360.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 8.4879E04_WP, 2.6079E05_WP /);
    pc%liqHeatCapacityMol%eqn = 114;    

    pc%IG_HeatCapacityMol%A = 5.1920E04_WP;
    pc%IG_HeatCapacityMol%B = 1.9245E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6265E03_WP;
    pc%IG_HeatCapacityMol%D = 1.1680E05_WP;
    pc%IG_HeatCapacityMol%E = 7.2360E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 200.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 5.6322E04_WP, 2.0556E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.1270E-01_WP;
    pc%SecondViralCoefMol%B = -9.9200E01_WP;
    pc%SecondViralCoefMol%C = -4.5100E06_WP;
    pc%SecondViralCoefMol%D = 3.0900E17_WP;
    pc%SecondViralCoefMol%E = -7.0500E19_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 184.91_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.1900E00_WP, 4.5230E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -6.9281E00_WP;
    pc%liqViscosity%B = 4.2076E02_WP;
    pc%liqViscosity%C = -6.3276E-01_WP;
    pc%liqViscosity%D = -1.7130E-26_WP;
    pc%liqViscosity%E = 1E01_WP;
    pc%liqViscosity%range(:,1) = (/ 85.47_WP, 360.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 8.0682E-03_WP, 4.0666E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 2.4993E-07_WP;
    pc%vapViscosity%B = 6.8612E-01_WP;
    pc%vapViscosity%C = 1.7934E02_WP;
    pc%vapViscosity%D = -8.2546E03_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 85.47_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 2.6865E-06_WP, 2.4411E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.6755E-01_WP;
    pc%liqThermalConductivity%B = -6.6457E-04_WP;
    pc%liqThermalConductivity%C = 2.7740E-07_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 85.47_WP, 350.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 2.1278E-01_WP, 6.8932E-02_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = -1.1200E00_WP;
    pc%vapThermalConductivity%B = 1.0972E-01_WP;
    pc%vapThermalConductivity%C = -9.8346E03_WP;
    pc%vapThermalConductivity%D = -7.5358E06_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 231.11_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.1142E-02_WP, 1.4599E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 5.0920E-02_WP;
    pc%SurfaceTension%B = 1.2197E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 85.47_WP, 369.83_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.6956E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_propane


  subroutine pc_tetrahydrofuran(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'tetrahydrofuran';
    pc%CAS = '109-99-9';
    ! components
    pc%ChemicalFormula%C = 4;         
    pc%ChemicalFormula%H = 8;
    pc%ChemicalFormula%O = 1;

    pc%MolecularWeight = 72.107_WP;    ! kg/kmol    
    pc%Tcrit = 540.15_WP;              ! K    
    pc%pcrit = 5.1678E6_WP;            ! Pa 
    pc%vcrit = 0.22393_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.259_WP;  ! -


    pc%MeltingPoint = 164.65_WP;             ! K
    pc%TripletPointT = 164.65_WP;            ! K
    pc%TripletPointp = 1.2753E-01_WP;        ! Pa
    pc%NormalBoilingPoint = 338.00_WP;       ! K
    pc%liqMolarVolume  = 0.081942_WP;        ! m^3/kmol


    pc%IG_HeatOfFormation = -1.8418E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -7.9680E07_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 2.9729E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = 8.5400E06_WP;          ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.3250E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.2264_WP;            ! -
    pc%RadiusOfGyration = 2.6940E-10_WP;      ! m
    pc%SolubilityParameter = 1.9129E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 5.4371E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.04462_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 6.0000E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4050_WP;           ! -
    pc%FlashPoint = 258.11_WP;                ! K
    pc%FlammabilityLimits = (/ 2.0, 11.8 /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 594.26_WP;          ! K


    pc%solDensityMol%A = 1.6400E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 164.65_WP, 1.6400E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 1.2543E00_WP;
    pc%liqDensityMol%B = 2.8084E-01_WP;
    pc%liqDensityMol%C = 5.4015E02_WP;
    pc%liqDensityMol%D = 2.9120E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 164.65_WP, 540.15_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.3998E01_WP, 4.4662E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  8.5184E01_WP;
    pc%vapPressure%B = -6.3373E03_WP;
    pc%vapPressure%C = -9.5958E00_WP;
    pc%vapPressure%D =  8.1661E-06_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 164.65_WP, 540.15_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2753E-01_WP, 5.1720E06_WP /)
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 4.4430E07_WP;
    pc%HeatOfVapMol%B = 3.9080E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 164.65_WP, 540.15_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 3.8545E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = -1.8200E04_WP;
    pc%solHeatCapacityMol%B = 1.4470E03_WP;
    pc%solHeatCapacityMol%C = -9.4700E00_WP;
    pc%solHeatCapacityMol%D = 2.6620E-02_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ 20.00_WP, 164.65_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ 7.1650E03_WP, 8.2141E04_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 1.7165E05_WP;
    pc%liqHeatCapacityMol%B = -7.9835E02_WP;
    pc%liqHeatCapacityMol%C = 2.8812E00_WP;
    pc%liqHeatCapacityMol%D = -2.4813E-03_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 164.65_WP, 320.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.0723E05_WP, 1.2991E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 4.0430E04_WP;
    pc%IG_HeatCapacityMol%B = 2.5410E05_WP;
    pc%IG_HeatCapacityMol%C = 1.3866E03_WP;
    pc%IG_HeatCapacityMol%D = 1.4200E05_WP;
    pc%IG_HeatCapacityMol%E = -6.3880E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 50.00_WP, 1500.00_WP /);
    pc%IG_HeatCapacityMol%range(:,2) = (/ 4.0430E04_WP, 2.5462E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.2580E-01_WP;
    pc%SecondViralCoefMol%B = -1.4320E02_WP;
    pc%SecondViralCoefMol%C = -2.0790E07_WP;
    pc%SecondViralCoefMol%D = 6.3140E19_WP;
    pc%SecondViralCoefMol%E = -1.9480E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 304.58_WP, 2701.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -1.0912E00_WP, 7.1728E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -1.0335E01_WP;
    pc%liqViscosity%B = 8.8360E02_WP;
    pc%liqViscosity%C = -5.2650E-02_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 164.65_WP, 350.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 5.3152E-03_WP, 2.9787E-04_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 3.7780E-07_WP;
    pc%vapViscosity%B = 6.5330E-01_WP;
    pc%vapViscosity%C = 2.7101E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 338.00_WP, 998.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 9.4123E-06_WP, 2.7056E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 1.9428E-01_WP;
    pc%liqThermalConductivity%B = -2.4900E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 164.65_WP, 338.00_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.5328E-01_WP, 1.1012E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 9.5521E-06_WP;
    pc%vapThermalConductivity%B = 1.4561E00_WP;
    pc%vapThermalConductivity%C = 6.6222E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 338.00_WP, 998.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 1.5534E-02_WP, 1.3369E-01_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.7130E-02_WP;
    pc%SurfaceTension%B = 1.2310E00_WP;
    pc%SurfaceTension%C = 0.0_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 164.65_WP, 540.15_WP /);
    pc%SurfaceTension%range(:,2) = (/ 4.2908E-02_WP, 0.0_WP /);
    pc%SurfaceTension%eqn = 106;    

  end subroutine pc_tetrahydrofuran


  subroutine pc_tetrahydrofurfurylalcohol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------




    ! all alternative names
    pc%name = 'tetrahydrofurfurylalcohol';
    pc%CAS = '97-99-4';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 2;

    pc%MolecularWeight = 102.133_WP;   ! kg/kmol    
    pc%Tcrit = 639.00_WP;              ! K    
    pc%pcrit = 4.6600E6_WP;            ! Pa 
    pc%vcrit = 0.29000_WP;             ! m^3/kmol 
    pc%CritCompressFactor = 0.254_WP;  ! -


    pc%MeltingPoint = 193.00_WP;             ! K
    pc%TripletPointT = MAX_REAL_WP;            ! K
    pc%TripletPointp = MAX_REAL_WP;            ! Pa
    pc%NormalBoilingPoint = 450.80_WP;    ! K
    pc%liqMolarVolume  = 0.097319_WP;     ! m^3/kmol


    pc%IG_HeatOfFormation = -3.6920E08_WP;            ! J/kmol
    pc%IG_GibbsOfFormation = -2.1600E08_WP;           ! J/kmol
    pc%IG_AbsoluteEntropy = 3.7400E5_WP;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.7412E9_WP;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    


    pc%AcentricFactor = 0.7044_WP;            ! -
    pc%RadiusOfGyration = 3.4810E-10_WP;      ! m
    pc%SolubilityParameter = 2.2920E4_WP;     ! (J/m^3)^.5 
    pc%DipoleMoment = 7.3400E-30_WP;          ! C^m
    pc%VanDerWaalsVolume = 0.05830_WP;        ! m^3/kmol


    pc%VanDerWaalsArea = 7.4600E08_WP;        ! m^2/kmol
    pc%RefractiveIndex = 1.4499_WP;           ! -
    pc%FlashPoint = 343.15_WP;                ! K
    pc%FlammabilityLimits = (/ 1.5_WP, 9.7_WP /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP, MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 553.15_WP;          ! K


    pc%solDensityMol%A = 1.2900E01_WP;
    pc%solDensityMol%B = 0.0_WP;
    pc%solDensityMol%C = 0.0_WP;
    pc%solDensityMol%D = 0.0_WP;
    pc%solDensityMol%E = 0.0_WP;
    pc%solDensityMol%range(:,1) = (/ 0.0_WP, 1.2900E01_WP /);
    pc%solDensityMol%eqn = 100;

    pc%liqDensityMol%A = 9.7001E-01_WP;
    pc%liqDensityMol%B = 2.8130E-01_WP;
    pc%liqDensityMol%C = 6.3900E02_WP;
    pc%liqDensityMol%D = 2.3837E-01_WP;
    pc%liqDensityMol%E = 0.0_WP;
    pc%liqDensityMol%range(:,1) = (/ 193.00_WP, 639.00_WP /);
    pc%liqDensityMol%range(:,2) = (/ 1.1046E01_WP, 3.4483E00_WP /);
    pc%liqDensityMol%eqn = 105;

    pc%vapPressure%A =  1.5059E02_WP;
    pc%vapPressure%B = -1.1574E04_WP;
    pc%vapPressure%C = -1.9025E01_WP;
    pc%vapPressure%D =  1.4141E-05_WP;
    pc%vapPressure%E =  2.0000E00_WP;
    pc%vapPressure%range(:,1) = (/ 193.00_WP, 639.00_WP /);
    pc%vapPressure%range(:,2) = (/ 1.2655E-04_WP, 4.6469E06_WP /)
    pc%vapPressure%eqn = 101;  

    pc%HeatOfVapMol%A = 6.4109E07_WP;
    pc%HeatOfVapMol%B = 2.8538E-01_WP;
    pc%HeatOfVapMol%C = 0.0_WP;
    pc%HeatOfVapMol%D = 0.0_WP;
    pc%HeatOfVapMol%E = 0.0_WP;
    pc%HeatOfVapMol%range(:,1) = (/ 193.00_WP, 639.00_WP /);
    pc%HeatOfVapMol%range(:,2) = (/ 5.7857E07_WP, 0.0_WP /);
    pc%HeatOfVapMol%eqn = 106;     

    pc%solHeatCapacityMol%A = 0.0_WP;
    pc%solHeatCapacityMol%B = 0.0_WP;
    pc%solHeatCapacityMol%C = 0.0_WP;
    pc%solHeatCapacityMol%D = 0.0_WP;
    pc%solHeatCapacityMol%E = 0.0_WP;
    pc%solHeatCapacityMol%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%range(:,2) = (/ MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacityMol%eqn = 100;

    pc%liqHeatCapacityMol%A = 5.2700E04_WP;
    pc%liqHeatCapacityMol%B = 4.3580E02_WP;
    pc%liqHeatCapacityMol%C = 0.0_WP;
    pc%liqHeatCapacityMol%D = 0.0_WP;
    pc%liqHeatCapacityMol%E = 0.0_WP;
    pc%liqHeatCapacityMol%range(:,1) = (/ 296.65_WP, 469.50_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.8198E05_WP, 2.5731E05_WP /);
    pc%liqHeatCapacityMol%eqn = 100;    

    pc%IG_HeatCapacityMol%A = 2.6823E05_WP;
    pc%IG_HeatCapacityMol%B = 1.1670E05_WP;
    pc%IG_HeatCapacityMol%C = 1.6860E03_WP;
    pc%IG_HeatCapacityMol%D = -3.5200E05_WP;
    pc%IG_HeatCapacityMol%E = 2.4520E02_WP;
    pc%IG_HeatCapacityMol%range(:,1) = (/ 300.00_WP, 1500.00_WP /);
    pc%liqHeatCapacityMol%range(:,2) = (/ 1.8198E05_WP, 2.5731E05_WP /);
    pc%IG_HeatCapacityMol%eqn = 107;

    pc%SecondViralCoefMol%A = 1.9766E-01_WP;
    pc%SecondViralCoefMol%B = -1.8908E02_WP;
    pc%SecondViralCoefMol%C = -8.1632E07_WP;
    pc%SecondViralCoefMol%D = -5.5857E20_WP;
    pc%SecondViralCoefMol%E = 8.4850E22_WP;
    pc%SecondViralCoefMol%range(:,1) = (/ 319.50_WP, 1500.00_WP /);
    pc%SecondViralCoefMol%range(:,2) = (/ -5.5954E00_WP, 4.7400E-02_WP /);
    pc%SecondViralCoefMol%eqn = 104;  

    pc%liqViscosity%A = -7.9742E00_WP;
    pc%liqViscosity%B = 2.7454E03_WP;
    pc%liqViscosity%C = -1.1468E00_WP;
    pc%liqViscosity%D = 0.0_WP;
    pc%liqViscosity%E = 0.0_WP;
    pc%liqViscosity%range(:,1) = (/ 220.00_WP, 503.00_WP /);
    pc%liqViscosity%range(:,2) = (/ 1.8628E-01_WP, 6.4431E-05_WP /);
    pc%liqViscosity%eqn = 101;

    pc%vapViscosity%A = 1.6196E-07_WP;
    pc%vapViscosity%B = 7.4453E-01_WP;
    pc%vapViscosity%C = 1.4579E02_WP;
    pc%vapViscosity%D = 0.0_WP;
    pc%vapViscosity%E = 0.0_WP;
    pc%vapViscosity%range(:,1) = (/ 193.00_WP, 1000.00_WP /);
    pc%vapViscosity%range(:,2) = (/ 4.6420E-06_WP, 2.4204E-05_WP /);
    pc%vapViscosity%eqn = 102;

    pc%liqThermalConductivity%A = 2.1782E-01_WP;
    pc%liqThermalConductivity%B = -2.4050E-04_WP;
    pc%liqThermalConductivity%C = 0.0_WP;
    pc%liqThermalConductivity%D = 0.0_WP;
    pc%liqThermalConductivity%E = 0.0_WP;
    pc%liqThermalConductivity%range(:,1) = (/ 193.00_WP, 450.80_WP /);
    pc%liqThermalConductivity%range(:,2) = (/ 1.7140E-01_WP, 1.0940E-01_WP /);
    pc%liqThermalConductivity%eqn = 100;

    pc%vapThermalConductivity%A = 2.8570E-04_WP;
    pc%vapThermalConductivity%B = 9.0229E-01_WP;
    pc%vapThermalConductivity%C = 7.0294E02_WP;
    pc%vapThermalConductivity%D = 0.0_WP;
    pc%vapThermalConductivity%E = 0.0_WP;
    pc%vapThermalConductivity%range(:,1) = (/ 450.80_WP, 1000.00_WP /);
    pc%vapThermalConductivity%range(:,2) = (/ 2.7698E-02_WP, 8.5424E-02_WP /);
    pc%vapThermalConductivity%eqn = 102; 

    pc%SurfaceTension%A = 6.5073E-02_WP;
    pc%SurfaceTension%B = -8.8258E-05_WP;
    pc%SurfaceTension%C = -1.6138E-08_WP;
    pc%SurfaceTension%D = 0.0_WP;
    pc%SurfaceTension%E = 0.0_WP;
    pc%SurfaceTension%range(:,1) = (/ 294.75_WP, 450.80_WP /);
    pc%SurfaceTension%range(:,2) = (/ 3.7657E-02_WP, 2.2007E-02_WP /);
    pc%SurfaceTension%eqn = 100;    

  end subroutine pc_tetrahydrofurfurylalcohol

end module pc_database
