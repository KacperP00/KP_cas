module pc_database
  use precision
  use pc_defs

  implicit none

contains

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

!!$
!!$  subroutine pc_1butanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$    ! all alternative names
!!$    pc%name = {'1butanol','1-butanol','1-Butanol','n-butanol','n-Butanol, Butan-1-ol'}; 
!!$    pc%CAS = '71-36-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 74.123_WP;   ! kg/kmol    
!!$    pc%Tcrit = 563.05_WP;             ! K    
!!$    pc%pcrit = 4.4230E06_WP;          ! Pa 
!!$    pc%vcrit = 0.27400_WP;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.260_WP; ! -
!!$
!!$
!!$    pc%MeltingPoint = 183.85_WP;              ! K
!!$    pc%TripletPointT = 184.51_WP;             ! K
!!$    pc%TripletPointp = 5.7220E-04_WP;         ! Pa
!!$    pc%NormalBoilingPoint = 390.81_WP;        ! K
!!$    pc%liqMolarVolume = 0.091991_WP;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.7460E08_WP;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5030E08_WP;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.6148E05_WP;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  9.3720E06_WP;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.4560E09_WP;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5935_WP;            ! -
!!$    pc%RadiusOfGyration = 3.2250E-10_WP;      ! m
!!$    pc%SolubilityParameter = 2.3350E04_WP;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5400E-30_WP;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05240_WP;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.6200E08_WP;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3971_WP;           ! -
!!$    pc%FlashPoint = 302.00_WP;                ! K
!!$    pc%FlammabilityLimits = (/ 1.4_WP, 11.2_WP/);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 616.00_WP;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.4100E01_WP;
!!$    pc%solDensity%B = 0.0_WP;
!!$    pc%solDensity%C = 0.0_WP;
!!$    pc%solDensity%D = 0.0_WP;
!!$    pc%solDensity%E = 0.0_WP;
!!$    pc%solDensity%range(:,1) = (/ 184.51_WP, MAX_REAL_WP /);
!!$    pc%solDensity%range(:,2) = (/ 1.4100E01_WP, MAX_REAL_WP /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.6500E-01_WP;
!!$    pc%liqDensity%B = 2.6660E-01_WP;
!!$    pc%liqDensity%C = 5.6305E02_WP;
!!$    pc%liqDensity%D = 2.4419E-01_WP;
!!$    pc%liqDensity%E = 0.0_WP;
!!$    pc%liqDensity%range(:,1) = (/ 184.51, 563.05 /);
!!$    pc%liqDensity%range(:,2) = (/ 1.2016E01, 3.6197E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  9.3173E01;
!!$    pc%vapPressure%B = -9.1859E03;
!!$    pc%vapPressure%C = -9.7464E00;
!!$    pc%vapPressure%D =  4.7796E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 184.51, 563.05 /);
!!$    pc%vapPressure%range(:,2) = (/ 5.7220E-04, 4.3392E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.7390E07;
!!$    pc%HeatOfVap%B = 1.7300E-01;
!!$    pc%HeatOfVap%C = 2.9150E-01;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 184.51, 563.05 /);
!!$    pc%HeatOfVap%range(:,2) = (/ 6.0575E07, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.1680E04;
!!$    pc%solHeatCapacity%B = 1.0262E03;
!!$    pc%solHeatCapacity%C = -1.7010E00;
!!$    pc%solHeatCapacity%D = -1.5820E-02;
!!$    pc%solHeatCapacity%E = 6.6500E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 184.51 /);
!!$    pc%solHeatCapacity%range(:,2) = (/ 8.0477E0, 9.7456E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.9120E05;
!!$    pc%liqHeatCapacity%B = -7.3040E02;
!!$    pc%liqHeatCapacity%C = 2.2998E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 184.51, 390.81 /);
!!$    pc%liqHeatCapacity%range(:,2) = (/ 1.3473E05, 2.5701E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 7.4540E04;
!!$    pc%IG_HeatCapacity%B = 2.5907E05;
!!$    pc%IG_HeatCapacity%C = 1.6073E03;
!!$    pc%IG_HeatCapacity%D = 1.7320E05;
!!$    pc%IG_HeatCapacity%E = 7.1240E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1500.00 /);
!!$    pc%IG_HeatCapacity%range(:,2) = (/ 8.1616E04, 2.8509E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.8900E-01;
!!$    pc%SecondViralCoef%B = -1.8200E02;
!!$    pc%SecondViralCoef%C = -4.0500E07;
!!$    pc%SecondViralCoef%D = -2.2700E20;
!!$    pc%SecondViralCoef%E = 4.5400E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 282.00, 1500.00 /);
!!$    pc%SecondViralCoef%range(:,2) = (/ 3.9128E00, 5.5659E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -3.5426E01;
!!$    pc%liqViscosity%B = 3.1845E03;
!!$    pc%liqViscosity%C = 3.2965E00;
!!$    pc%liqViscosity%D = -3.0000E-27;
!!$    pc%liqViscosity%E = 1.0000E01;
!!$    pc%liqViscosity%range(:,1) = (/ 190.00, 390.81 /);
!!$    pc%liqViscosity%range(:,2) = (/ 2.5440E-01, 3.8873E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.4031E-06;
!!$    pc%vapViscosity%B = 4.6110E-01;
!!$    pc%vapViscosity%C = 5.3700E+02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 184.51, 1000.00 /);
!!$    pc%vapViscosity%range(:,2) = (/ 3.9786E-06, 2.2066E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1360E-01;
!!$    pc%liqThermalConductivity%B = -2.0340E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 184.51, 390.81 /);
!!$    pc%liqThermalConductivity%range(:,2) = (/ 1.7607E-01, 1.3411E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -4.4940E-02;
!!$    pc%vapThermalConductivity%B = 4.4600E-02;
!!$    pc%vapThermalConductivity%C = -1.3552E+03;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 370.00, 800.00 /);
!!$    pc%vapThermalConductivity%range(:,2) = (/ 2.1971E-02, 8.7247E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 4.9830E-02;
!!$    pc%SurfaceTension%B = -8.5400E-05;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 273.15, 413.15 /);
!!$    pc%SurfaceTension%range(:,2) = (/ 2.6503E-02, 1.4547E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_1butanol
!!$
!!$
!!$  subroutine pc_1decanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1decanol','1-decanol','1-Decanol','n-decanol','n-Decanol', 'decyl alcohol'}; 
!!$    pc%CAS = '112-30-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 10;         
!!$    pc%ChemicalFormula%H = 22;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 158.284;  ! kg/kmol    
!!$    pc%Tcrit = 687.30;             ! K    
!!$    pc%pcrit = 2.3150E06;          ! Pa 
!!$    pc%vcrit = 0.64900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.263; ! -
!!$
!!$
!!$    pc%MeltingPoint = 280.05;              ! K
!!$    pc%TripletPointT = 280.05;             ! K
!!$    pc%TripletPointp = 1.2019E-01;         ! Pa
!!$    pc%NormalBoilingPoint = 504.07;        ! K
!!$    pc%liqMolarVolume = 0.191854;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.9670E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -9.8790E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  5.9700E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 3.7656E07;           ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -6.1159E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.2650E08;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -4.7810E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 4.1690E05;         ! J/kmolK
!!$
!!$
!!$    pc%AcentricFactor = 0.6219;            ! -
!!$    pc%RadiusOfGyration = 5.4990E-10;      ! m
!!$    pc%SolubilityParameter = 2.02500E04;   ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.4000E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.11378;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.5730E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4350;           ! -
!!$    pc%FlashPoint = 355.37;                ! K
!!$    pc%FlammabilityLimits = (/ 0.7 5.5 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 376 416 /); ! K
!!$    pc%AutoignitionTemp = 561.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 6.1900E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 280.05, MAX_REAL_WP /);
!!$    pc%solDensity%range(:,2) = (/ 6.1900E00, MAX_REAL_WP /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.7384E-01;
!!$    pc%liqDensity%B = 2.4241E-01;
!!$    pc%liqDensity%C = 6.8730E02;
!!$    pc%liqDensity%D = 2.6646E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 280.05, 5.2902E00; 687.30, 1.5422E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  2.5059E02;
!!$    pc%vapPressure%B = -1.9169E04;
!!$    pc%vapPressure%C = -3.2903E01;
!!$    pc%vapPressure%D =  1.4627E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 280.05, 1.2019E-01; 687.30, 2.3302E06 /);
!!$    pc%vapPressure%eqn = 101;
!!$
!!$    pc%HeatOfVap%A = 1.1750E08;
!!$    pc%HeatOfVap%B = 6.5112E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 280.05, 8.3570E07; 687.30, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 4.11026E05;
!!$    pc%solHeatCapacity%B = -2.4283E03;
!!$    pc%solHeatCapacity%C = 7.3148E00;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 180.00, 2.8017E05; 270.00, 3.5787E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 4.9885E06;
!!$    pc%liqHeatCapacity%B = -5.2898E04;
!!$    pc%liqHeatCapacity%C = 2.1635E02;
!!$    pc%liqHeatCapacity%D = -3.7538E-01;
!!$    pc%liqHeatCapacity%E = 2.3674E-04;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 280.05, 3.5373E05; 504.07, 5.0206E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.6984E05;
!!$    pc%IG_HeatCapacity%B = 5.5920E05;
!!$    pc%IG_HeatCapacity%C = 1.5680E03;
!!$    pc%IG_HeatCapacity%D = 3.9380E05;
!!$    pc%IG_HeatCapacity%E = 7.2050E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.8502E05; 1500.00, 6.2186E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.8836E-01;
!!$    pc%SecondViralCoef%B = -5.8491E02;
!!$    pc%SecondViralCoef%C = -1.1672E08;
!!$    pc%SecondViralCoef%D = 3.8936E21;
!!$    pc%SecondViralCoef%E = -1.5303E24;
!!$    pc%SecondViralCoef%range(:,1) = (/ 343.65, -6.9662E00; 1500.00, 6.3948E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -8.0656E01;
!!$    pc%liqViscosity%B = 6.3255E03;
!!$    pc%liqViscosity%C = 9.6460E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 285.00, 1.9500E-02; 504.07, 3.0879E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 5.5065E-08;
!!$    pc%vapViscosity%B = 8.3410E-01;
!!$    pc%vapViscosity%C = 7.9560E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 280.05, 4.7154E-06; 1000.00, 1.6215E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2800E-01;
!!$    pc%liqThermalConductivity%B = -2.2300E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 280.05, 1.6555E-01; 503.35, 1.1575E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -3.0720E-01;
!!$    pc%vapThermalConductivity%B = 4.8900E-01;
!!$    pc%vapThermalConductivity%C = -6.7500E04;
!!$    pc%vapThermalConductivity%D = -2.9400E07;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 504.07, 2.5906E-02; 1000.00, 9.3886E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1263E-02;
!!$    pc%SurfaceTension%B = 1.0395E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 280.05, 2.9754E-02; 687.30, 0.0000E00 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1decanol
!!$
!!$
!!$  subroutine pc_1dodecanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1dodecanol','1-dodecanol','1-Dodecanol','n-dodecanol','Lauryl 24'}; 
!!$    pc%CAS = '112-53-8';
!!$    ! components
!!$    pc%ChemicalFormula%C = 12;         
!!$    pc%ChemicalFormula%H = 26;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 186.338;  ! kg/kmol    
!!$    pc%Tcrit = 719.4;              ! K    
!!$    pc%pcrit = 1.9940E06;          ! Pa 
!!$    pc%vcrit = 0.716;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.239; ! -
!!$
!!$
!!$    pc%MeltingPoint = 296.95;              ! K
!!$    pc%TripletPointT = 296.95;             ! K
!!$    pc%TripletPointp = 0.101105;           ! Pa
!!$    pc%NormalBoilingPoint = 536.95;        ! K
!!$    pc%liqMolarVolume = 0.224322;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -4.3670E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -8.0820E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  6.7521E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 3.1380E07;           ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -7.3380E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.1500E08;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -5.2850E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 4.8300E05;         ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.666353;          ! -
!!$    pc%RadiusOfGyration = 6.1190E-10;      ! m
!!$    pc%SolubilityParameter = 1.9750E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.6400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.13424;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.8430E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4413;           ! -
!!$    pc%FlashPoint = 400.0;                 ! K
!!$    pc%FlammabilityLimits = (/ 0.6 5.1 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 397 441 /); ! K
!!$    pc%AutoignitionTemp = 548.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 5.5881E00;
!!$    pc%solDensity%B = -1.6655E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 78.15, 5.4579E00; 296.95, 5.0935E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.5818E-01;
!!$    pc%liqDensity%B = 2.5626E-01;
!!$    pc%liqDensity%C = 7.1940E02;
!!$    pc%liqDensity%D = 2.9963E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 296.95, 4.4623E00; 719.40, 1.3977E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  3.8678E02;
!!$    pc%vapPressure%B = -2.4382E04;
!!$    pc%vapPressure%C = -5.6170E01;
!!$    pc%vapPressure%D =  4.3253E-02;
!!$    pc%vapPressure%E = 1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 296.95,1.011E-01; 719.40,1.9716E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 1.2925E08;
!!$    pc%HeatOfVap%B = 6.7609E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 296.95, 9.0182E07; 719.40, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 3.5700E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 293.15, 3.5700E05; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 5.1422E06;
!!$    pc%liqHeatCapacity%B = -5.4970E04;
!!$    pc%liqHeatCapacity%C = 2.2872E02;
!!$    pc%liqHeatCapacity%D = -4.0331E-01;
!!$    pc%liqHeatCapacity%E = 2.5844E-04;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 296.95,4.3614E05; 486.00, 5.7106E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.0150E05;
!!$    pc%IG_HeatCapacity%B = 6.3270E05;
!!$    pc%IG_HeatCapacity%C = 1.5625E03;
!!$    pc%IG_HeatCapacity%D = 4.6733E05;
!!$    pc%IG_HeatCapacity%E = 7.2120E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 2.1943E05; 1500.00, 7.3406E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 6.2391E-01;
!!$    pc%SecondViralCoef%B = -7.8675E02;
!!$    pc%SecondViralCoef%C = -1.6374E08;
!!$    pc%SecondViralCoef%D = 7.1921E21;
!!$    pc%SecondViralCoef%E = -2.8712E24;
!!$    pc%SecondViralCoef%range(:,1) = (/ 359.70,-7.9011E00; 1500.00, 5.1100E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -5.3664E01;
!!$    pc%liqViscosity%B = 5.2265E03;
!!$    pc%liqViscosity%C = 5.6140E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 293.15, 1.9373E-02; 536.95, 1.7667E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 5.5812E-08;
!!$    pc%vapViscosity%B = 8.2000E-01;
!!$    pc%vapViscosity%C = 5.4330E01;
!!$    pc%vapViscosity%D = 1.4110E04;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 296.95, 4.4285E-06; 1000.00, 1.5065E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0444E-01;
!!$    pc%liqThermalConductivity%B = -1.8378E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 296.95, 1.4987E-01; 536.95, 1.0576E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -2.8300E-03;
!!$    pc%vapThermalConductivity%B = 9.1261E-01;
!!$    pc%vapThermalConductivity%C = -1.9181E04;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 536.95, 2.5266E-02; 1000.00, 8.5114E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.3034E-02;
!!$    pc%SurfaceTension%B = 1.1008E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 293.15, 2.9808E-02; 719.40, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1dodecanol
!!$
!!$
!!$  subroutine pc_1heptanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1heptanol','1-heptanol','1-Heptanol','n-heptanol','n-Heptanol', 'heptyl alcohol'}; 
!!$    pc%CAS = '111-70-6';
!!$    ! components
!!$    pc%ChemicalFormula%C = 7;         
!!$    pc%ChemicalFormula%H = 16;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 116.203;  ! kg/kmol    
!!$    pc%Tcrit = 632.60;             ! K    
!!$    pc%pcrit = 3.0580E06;          ! Pa 
!!$    pc%vcrit = 0.43500;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.253; ! -
!!$
!!$
!!$    pc%MeltingPoint = 239.15;              ! K
!!$    pc%TripletPointT = 239.15;             ! K
!!$    pc%TripletPointp = 1.1773E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 449.45;        ! K
!!$    pc%liqMolarVolume = 0.141801;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.3640E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.2530E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.7919E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.8180E07;           ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.2860E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.4650E08;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -4.0330E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.2590E05;         ! J/kmolK
!!$
!!$
!!$    pc%AcentricFactor = 0.5670;            ! -
!!$    pc%RadiusOfGyration = 4.3800E-10;      ! m
!!$    pc%SolubilityParameter = 2.1540E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.8000E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.08309;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.1680E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4223;           ! -
!!$    pc%FlashPoint = 350.00;                ! K
!!$    pc%FlammabilityLimits = (/ 1.0 7.2 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 344 379 /); ! K
!!$    pc%AutoignitionTemp = 555.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 8.6800E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 239.15, 8.6800E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.0531E-01;
!!$    pc%liqDensity%B = 2.6329E-01;
!!$    pc%liqDensity%C = 6.3260E02;
!!$    pc%liqDensity%D = 2.7375E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 239.15, 7.4209E00; 632.60, 2.2990E+00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  2.4331E02;
!!$    pc%vapPressure%B = -1.7320E04;
!!$    pc%vapPressure%C = -3.2183E01;
!!$    pc%vapPressure%D =  1.6443E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 239.15, 1.1773E-02; 632.60, 3.0663E+06 /);
!!$    pc%vapPressure%eqn = 101;
!!$
!!$    pc%HeatOfVap%A = 1.0075E08;
!!$    pc%HeatOfVap%B = 6.1052E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 239.15, 7.5393E07; 632.60, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 1.9300E04;
!!$    pc%solHeatCapacity%B = 9.8150E02;
!!$    pc%solHeatCapacity%C = -3.0750E00;
!!$    pc%solHeatCapacity%D = 6.9530E-03;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 80.00, 8.1700E04; 239.15, 1.7326E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$
!!$    pc%liqHeatCapacity%A = 2.4168E06;
!!$    pc%liqHeatCapacity%B = -2.6105E04;
!!$    pc%liqHeatCapacity%C = 1.1003E02;
!!$    pc%liqHeatCapacity%D = -1.9172E-01;
!!$    pc%liqHeatCapacity%E = 1.1968E-04;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 239.15, 2.3590E05; 449.45, 3.8773E+05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2215E05;
!!$    pc%IG_HeatCapacity%B = 3.9910E05;
!!$    pc%IG_HeatCapacity%C = 1.5800E03;
!!$    pc%IG_HeatCapacity%D = 2.8350E05;
!!$    pc%IG_HeatCapacity%E = 7.1770E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.3330E05; 1500.00, 4.5346E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.8979E-01;
!!$    pc%SecondViralCoef%B = -3.0797E02;
!!$    pc%SecondViralCoef%C = -9.8592E07;
!!$    pc%SecondViralCoef%D = -2.7797E20;
!!$    pc%SecondViralCoef%E = 1.3614E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 316.30, -6.1382E00; 1500.00, 5.5313E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -7.5245E01;
!!$    pc%liqViscosity%B = 5.7567E03;
!!$    pc%liqViscosity%C = 8.9146E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 239.15, 9.5615E-02; 449.45, 3.4016E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.5720E-07;
!!$    pc%vapViscosity%B = 6.5020E-01;
!!$    pc%vapViscosity%C = 2.4860E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 239.15, 4.4397E-06; 1000.00, 1.8384E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2390E-01;
!!$    pc%liqThermalConductivity%B = -2.2600E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 239.15, 1.6985E-01; 573.15, 9.4368E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -6.1993E-02;
!!$    pc%vapThermalConductivity%B = 2.7920E-01;
!!$    pc%vapThermalConductivity%C = -3.3360E03;
!!$    pc%vapThermalConductivity%D = -1.6420E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 504.07, 2.3447E-02; 1000.00, 1.0722E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.3765E-02;
!!$    pc%SurfaceTension%B = 1.0974E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 239.15, 3.1928E-02; 632.60, 0.0000E00 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1heptanol
!!$
!!$
!!$  subroutine pc_1hexanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1hexanol','1-hexanol','1-Hexanol','1Hexanol','n-Hexanol','n-hexanol',}; 
!!$    pc%CAS = '111-27-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 102.77;   ! kg/kmol    
!!$    pc%Tcrit = 610.30;             ! K    
!!$    pc%pcrit = 3.4170E06;          ! Pa 
!!$    pc%vcrit = 0.38700;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.261; ! -
!!$
!!$
!!$    pc%MeltingPoint = 228.55;              ! K
!!$    pc%TripletPointT = 228.55;             ! K
!!$    pc%TripletPointp = 2.9550E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 430.55;        ! K
!!$    pc%liqMolarVolume = 0.125200;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.1650E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.3440E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 4.4010E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.5400E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.6766E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = - 1.4980E+08;    ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.7750E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.8740E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5764;            ! -
!!$    pc%RadiusOfGyration = 4.1440E-10;      ! m
!!$    pc%SolubilityParameter = 2.1830E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5000E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.07286;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.0330E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4161;           ! -
!!$    pc%FlashPoint = 356.00;                ! K
!!$    pc%FlammabilityLimits = (/ 1.2 8.2 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 332 366 /); ! K
!!$    pc%AutoignitionTemp = 558.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 9.8900E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 228.55, 9.8900E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.8425E-01;
!!$    pc%liqDensity%B = 2.6477E-01;
!!$    pc%liqDensity%C = 6.1030E02;
!!$    pc%liqDensity%D = 2.4395E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 228.55, 8.4540E00; 610.30, 2.5843E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.2508E02;
!!$    pc%vapPressure%B = -1.1692E04;
!!$    pc%vapPressure%C = -1.4258E01;
!!$    pc%vapPressure%D =  1.1102E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 228.55, 2.9550E-02; 610.30, 3.4221E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 8.9962E07;
!!$    pc%HeatOfVap%B = 5.5179E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 228.55, 6.9442E07; 610.30, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.8080E04;
!!$    pc%solHeatCapacity%B = 1.4083E03;
!!$    pc%solHeatCapacity%C = -4.2360E00;
!!$    pc%solHeatCapacity%D = -1.1200E-03;
!!$    pc%solHeatCapacity%E = 3.0800E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 15.00, 2.0892E+03; 200.00, 1.3446E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.6386E06;
!!$    pc%liqHeatCapacity%B = -1.7261E04;
!!$    pc%liqHeatCapacity%C = 7.1721E01;
!!$    pc%liqHeatCapacity%D = -1.2026E-01;
!!$    pc%liqHeatCapacity%E = 7.1087E-05;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 228.55, 1.9821E05; 460.00, 3.5197E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0625E05;
!!$    pc%IG_HeatCapacity%B = 3.5210E05;
!!$    pc%IG_HeatCapacity%C = 1.5835E03;
!!$    pc%IG_HeatCapacity%D = 2.4620E05;
!!$    pc%IG_HeatCapacity%E = 7.1575E+02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.1607E05; 1500.00, 3.9726E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.5305E-01;
!!$    pc%SecondViralCoef%B = -2.6298E02;
!!$    pc%SecondViralCoef%C = -7.5678E07;
!!$    pc%SecondViralCoef%D = -2.2383E20;
!!$    pc%SecondViralCoef%E = 2.0576E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 305.15, -5.352E00; 1500.00, 5.5299E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -4.3289E01;
!!$    pc%liqViscosity%B = 4.0412E03;
!!$    pc%liqViscosity%C = 4.2714E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 250.00, 2.9018E-02; 430.55, 3.3662E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.5773E-07;
!!$    pc%vapViscosity%B = 7.1890E-01;
!!$    pc%vapViscosity%C = 1.6330E+02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 228.55, 4.5672E-06; 1000.00, 1.9450E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1930E-01;
!!$    pc%liqThermalConductivity%B = -2.2000E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 228.55, 1.6902E-01; 575.00, 9.2800E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -4.9355E06;
!!$    pc%vapThermalConductivity%B = -1.6530E-01;
!!$    pc%vapThermalConductivity%C = 1.5631E09;
!!$    pc%vapThermalConductivity%D = -1.5752E13;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 430.55, 2.2264E-02; 1000.0, 1.1104E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0895E-02;
!!$    pc%SurfaceTension%B = 1.0073E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 228.55, 3.1727E-02; 610.30, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1hexanol
!!$
!!$
!!$  subroutine pc_1nonanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1nonanol','1-nonanol','1Nonanol','1-Nonanol','Alcohol C-9','n-Nonyl Alcohol'}; 
!!$    pc%CAS = '143-08-8';
!!$    ! components
!!$    pc%ChemicalFormula%C = 9;         
!!$    pc%ChemicalFormula%H = 20;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 144.257;   ! kg/kmol    
!!$    pc%Tcrit = 670.70;             ! K    
!!$    pc%pcrit = 2.5280E06;           ! Pa 
!!$    pc%vcrit = 0.57200;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.259; ! -
!!$
!!$
!!$    pc%MeltingPoint = 268.15;              ! K
!!$    pc%TripletPointT = 268.15;             ! K
!!$    pc%TripletPointp = 7.1702E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 486.25;        ! K
!!$    pc%liqMolarVolume = 0.174942;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.8120E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.1220E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  5.5770E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.8800E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -5.5007E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.3300E08;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -4.5360E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.8460E05;         ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5997;            ! -
!!$    pc%RadiusOfGyration = 5.1080E-10;      ! m
!!$    pc%SolubilityParameter = 2.0720E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.3700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.10355;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.4380E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4319;           ! -
!!$    pc%FlashPoint = 347.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.8 6.1 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 366 405 /); ! K
!!$    pc%AutoignitionTemp = 550.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 6.8500E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 268.15, 6.8500E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 4.4231E-01;
!!$    pc%liqDensity%B = 2.5300E-01;
!!$    pc%liqDensity%C = 6.7070E+02;
!!$    pc%liqDensity%D = 2.5261E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 268.15, 5.8515E00; 670.70, 1.7483E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  2.5234E02;
!!$    pc%vapPressure%B = -1.8775E04;
!!$    pc%vapPressure%C = -3.3278E01;
!!$    pc%vapPressure%D =  1.5535E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 268.15, 7.1702E-02; 670.70, 2.5489E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 1.1484E08;
!!$    pc%HeatOfVap%B = 6.6695E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 268.15, 8.1701E07; 670.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 3.8373E05;
!!$    pc%solHeatCapacity%B = -1.4621E03;
!!$    pc%solHeatCapacity%C = 4.4071E00;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 200.00, 2.6759E05; 260.00, 3.0150E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.0483E07;
!!$    pc%liqHeatCapacity%B = -1.1522E05;
!!$    pc%liqHeatCapacity%C = 4.7687E02;
!!$    pc%liqHeatCapacity%D = -8.5381E-01;
!!$    pc%liqHeatCapacity%E = 5.6246E-04;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 310.00, 3.5059E05; 460.00, 4.6494E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.5400E05;
!!$    pc%IG_HeatCapacity%B = 4.9360E05;
!!$    pc%IG_HeatCapacity%C = 1.580E03;
!!$    pc%IG_HeatCapacity%D = 3.5880E05;
!!$    pc%IG_HeatCapacity%E = 7.2111E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.6777E05; 1500.00, 5.6606E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.3008E-01;
!!$    pc%SecondViralCoef%B = -5.1453E02;
!!$    pc%SecondViralCoef%C = -9.3063E07;
!!$    pc%SecondViralCoef%D = 2.9028E+21;
!!$    pc%SecondViralCoef%E = - 1.1382E24;
!!$    pc%SecondViralCoef%range(:,1) = (/ 335.35, -6.6432E00; 1500.00, 5.9569E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -7.1348E00;
!!$    pc%liqViscosity%B = 2.7763E03;
!!$    pc%liqViscosity%C = -1.2064E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 280.00, 1.8002E-02; 486.25, 1.3792E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.2000E-07;
!!$    pc%vapViscosity%B = 7.4000E-01;
!!$    pc%vapViscosity%C = 1.8000E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 268.15, 4.4992E-06; 1000.00, 1.6877E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2920E-01;
!!$    pc%liqThermalConductivity%B = -2.3000E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 268.15, 1.6753E-01; 578.65, 9.6111E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -3.0715E01;
!!$    pc%vapThermalConductivity%B = -1.0750E-01;
!!$    pc%vapThermalConductivity%C = 8.1070E03;
!!$    pc%vapThermalConductivity%D = -1.5683E08;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 486.25, 2.4464E-02; 1000.0, 9.8947E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1752E-02;
!!$    pc%SurfaceTension%B = 1.0562E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 268.15, 3.0183E-02; 670.70, 0.0000E00 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1nonanol
!!$
!!$
!!$  subroutine pc_1octanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1octanol','1-octanol','1-Octanol','n-octanol','n-Octanol, octanol'}; 
!!$    pc%CAS = '111-87-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 8;         
!!$    pc%ChemicalFormula%H = 18;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 130.230;  ! kg/kmol    
!!$    pc%Tcrit = 652.50;             ! K    
!!$    pc%pcrit = 2.7770E06;          ! Pa 
!!$    pc%vcrit = 0.49700;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.254; ! -
!!$
!!$
!!$    pc%MeltingPoint = 257.65;              ! K
!!$    pc%TripletPointT = 257.65;             ! K
!!$    pc%TripletPointp = 8.5564E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 468.35;        ! K
!!$    pc%liqMolarVolume = 0.158231;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.5550E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.1550E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  5.1850E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.2500E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.8983E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.3690E08;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -4.2650E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.5200E05;         ! J/kmolK
!!$
!!$
!!$    pc%AcentricFactor = 0.5829;            ! -
!!$    pc%RadiusOfGyration = 4.7870E-10;      ! m
!!$    pc%SolubilityParameter = 2.1020E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5000E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.09332;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.3030E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4276;           ! -
!!$    pc%FlashPoint = 354.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.9 6.4 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 355 392 /); ! K
!!$    pc%AutoignitionTemp = 555.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 7.6600E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 257.65, 7.6600E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.1836E-01;
!!$    pc%liqDensity%B = 2.5771E-01;
!!$    pc%liqDensity%C = 6.5250E02;
!!$    pc%liqDensity%D = 2.7713E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 257.65, 6.5440E00; 652.50, 2.0114E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.5390E02;
!!$    pc%vapPressure%B = -1.4181E+04;
!!$    pc%vapPressure%C = -1.8251E+01;
!!$    pc%vapPressure%D =  1.2261E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 257.65, 8.5564E-02; 652.50, 2.7625E+06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 1.0730E+08;
!!$    pc%HeatOfVap%B = 6.4524E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 257.65, 7.7596E+07; 652.50, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 4.1740E+04;
!!$    pc%solHeatCapacity%B = 5.7900E+02;
!!$    pc%solHeatCapacity%C = 2.7800E-01;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 101.90, 1.0363E+05; 216.90 1.8040E+05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 5.7137E+05;
!!$    pc%liqHeatCapacity%B = -4.8490E+03;
!!$    pc%liqHeatCapacity%C = 1.9725E+01;
!!$    pc%liqHeatCapacity%D = -2.1532E-02;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 250.00, 2.5550E+05; 468.35, 4.1500E+05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.3805E+05;
!!$    pc%IG_HeatCapacity%B = 4.4590E+05;
!!$    pc%IG_HeatCapacity%C = 1.5751E+03;
!!$    pc%IG_HeatCapacity%D = 3.2016E+05;
!!$    pc%IG_HeatCapacity%E = 7.1880E+02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.5055E+05; 1500.00, 5.0965E+05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.2384E-01;
!!$    pc%SecondViralCoef%B = -3.5035E+02;
!!$    pc%SecondViralCoef%C = -1.2731E+08;
!!$    pc%SecondViralCoef%D = -3.0518E+20;
!!$    pc%SecondViralCoef%E = -1.1146E+22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 326.25, -7.0600E+00; 1500.00, 5.2540E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -5.0674E+01;
!!$    pc%liqViscosity%B = 4.7251E+03;
!!$    pc%liqViscosity%C = 5.2499E+00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 280.00, 1.4749E-02; 425.00, 4.1665E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.7520E-07;
!!$    pc%vapViscosity%B = 6.9410E-01;
!!$    pc%vapViscosity%C = 2.0680E+02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 257.65, 4.5827E-06; 1000.00, 1.7547E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.3160E-01;
!!$    pc%liqThermalConductivity%B = -2.4070E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 257.65, 1.6958E-01; 570.15, 9.4365E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -3.0238E-03;
!!$    pc%vapThermalConductivity%B = 8.7450E-01;
!!$    pc%vapThermalConductivity%C = -1.3352E+04;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 468.35, 2.3796E-02; 1000.00, 1.0288E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1595E-02;
!!$    pc%SurfaceTension%B = 1.0483E+00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 257.65, 3.0355E-02; 652.50, 0.0000E+00 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1octanol
!!$
!!$
!!$  subroutine pc_1pentanal(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1pentanal','1-pentanal','1-Pentanal','1Pentanal','n-Pentanal','n-pentanal'}; 
!!$    pc%CAS = '110-62-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 86.134;   ! kg/kmol    
!!$    pc%Tcrit = 566.10;             ! K    
!!$    pc%pcrit = 3.9700E06;          ! Pa 
!!$    pc%vcrit = 0.31300;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.264; ! -
!!$
!!$
!!$    pc%MeltingPoint = 182.00;              ! K
!!$    pc%TripletPointT = 182.00;             ! K
!!$    pc%TripletPointp = 5.2282E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 376.15;        ! K
!!$    pc%liqMolarVolume = 0.107035;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.2780E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.0710E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.8289E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.5000E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.9100E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3472;            ! -
!!$    pc%RadiusOfGyration = 3.5160E-10;      ! m
!!$    pc%SolubilityParameter = 1.8440E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 8.5700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05950;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.5400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3917;           ! -
!!$    pc%FlashPoint = 285.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.5 9.5 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 480.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.2300E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 182.00, 1.2300E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.3871E-01;
!!$    pc%liqDensity%B = 2.6252E-01;
!!$    pc%liqDensity%C = 5.6610E02;
!!$    pc%liqDensity%D = 2.9444E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 182.00, 1.0534E01; 566.10, 3.1948E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.4958E02;
!!$    pc%vapPressure%B = -8.8900E03;
!!$    pc%vapPressure%C = -2.0697E01;
!!$    pc%vapPressure%D =  2.2101E-02;
!!$    pc%vapPressure%E =  1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 182.00, 5.2282E-02; 566.10, 3.9685E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.1478E07;
!!$    pc%HeatOfVap%B = 3.7541E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 182.00, 4.4502E07; 566.10, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 9.3786E03;
!!$    pc%solHeatCapacity%B = 7.4571E02;
!!$    pc%solHeatCapacity%C = -8.5714E-01;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 80.00, 6.3550E04; 150.00, 1.0195E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.1205E05;
!!$    pc%liqHeatCapacity%B = 2.5778E02;
!!$    pc%liqHeatCapacity%C = 0;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 200.00, 1.6361E05; 376.15, 2.0901E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0743E05;
!!$    pc%IG_HeatCapacity%B = 2.8363E05;
!!$    pc%IG_HeatCapacity%C = 1.9549E03;
!!$    pc%IG_HeatCapacity%D = 2.0146E05;
!!$    pc%IG_HeatCapacity%E = 8.9044E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.0960E05; 1500.00, 3.2404E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.9740E-01;
!!$    pc%SecondViralCoef%B = -2.2484E02;
!!$    pc%SecondViralCoef%C = -4.0107E07;
!!$    pc%SecondViralCoef%D = -1.9523E20;
!!$    pc%SecondViralCoef%E = 3.0063E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 283.05, -4.5262E00; 1500.00, 3.5616E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.0846E01;
!!$    pc%liqViscosity%B = 9.8001E02;
!!$    pc%liqViscosity%C = -5.4565E-03;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 182.00, 4.1290E-03; 376.15, 2.5532E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.2700E-07;
!!$    pc%vapViscosity%B = 6.7670E-01;
!!$    pc%vapViscosity%C = 1.9174E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 182.00, 3.7404E-06; 1000.00, 2.0415E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2697E-01;
!!$    pc%liqThermalConductivity%B = -3.3227E-04;
!!$    pc%liqThermalConductivity%C = 1.1770E-07;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 182.00, 1.7040E-01; 513.15, 8.7459E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -4.9187E06;
!!$    pc%vapThermalConductivity%B = -1.0297E-01;
!!$    pc%vapThermalConductivity%C = 2.6911E09;
!!$    pc%vapThermalConductivity%D = -2.3179E13;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 376.15, 1.7048E-02; 1000.0, 1.1788E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.4779E-02;
!!$    pc%SurfaceTension%B = 1.2556E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 182.00, 3.9804E-02; 566.10, 0.0000E00 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1pentanal
!!$
!!$
!!$  subroutine pc_1pentanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1pentanol','1-pentanol','1-Pentanol','1Pentanol','n-Pentanol','n-pentanol'}; 
!!$    pc%CAS = '71-41-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 88.150;   ! kg/kmol    
!!$    pc%Tcrit = 588.10;             ! K    
!!$    pc%pcrit = 3.8970E06;          ! Pa 
!!$    pc%vcrit = 0.32600;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.260; ! -
!!$
!!$
!!$    pc%MeltingPoint = 195.56;              ! K
!!$    pc%TripletPointT = 195.56;             ! K
!!$    pc%TripletPointp = 3.4525E-4;          ! Pa
!!$    pc%NormalBoilingPoint = 410.95;        ! K
!!$    pc%liqMolarVolume = 0.108538;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.9874E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.4602E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.0250E05;            ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  9.7906E06;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.0605E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.5610E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.5160E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.5890E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5731;            ! -
!!$    pc%RadiusOfGyration = 3.6790E-10;      ! m
!!$    pc%SolubilityParameter = 2.2520E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.6706E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06263;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.9800E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4080;           ! -
!!$    pc%FlashPoint = 305.93;                ! K
!!$    pc%FlammabilityLimits = (/ 1.2 10.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 316 353 /); ! K
!!$    pc%AutoignitionTemp = 573.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1800E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 195.56, 1.1800E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.1754E-01;
!!$    pc%liqDensity%B = 2.6732E-01;
!!$    pc%liqDensity%C = 5.8810E02;
!!$    pc%liqDensity%D = 2.5348E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 195.56, 1.0061E01; 588.10, 3.0583E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.6264E02;
!!$    pc%vapPressure%B = -1.2413E04;
!!$    pc%vapPressure%C = -2.0383E01;
!!$    pc%vapPressure%D =  1.0482E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 195.56, 3.4524E-04; 588.10, 3.9062E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 8.3972E07;
!!$    pc%HeatOfVap%B = 5.3507E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 195.56, 6.7638E07; 588.10, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.8050E04;
!!$    pc%solHeatCapacity%B = 1.6164E03;
!!$    pc%solHeatCapacity%C = -9.4040E00;
!!$    pc%solHeatCapacity%D = 2.4970E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 31.56, 2.4382E04; 178.06, 1.1258E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.0120E05;
!!$    pc%liqHeatCapacity%B = -6.5130E02;
!!$    pc%liqHeatCapacity%C = 2.2750E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 200.14, 1.6198E05; 389.15, 2.9227E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 9.0600E04;
!!$    pc%IG_HeatCapacity%B = 3.0620E05;
!!$    pc%IG_HeatCapacity%C = 1.6054E03;
!!$    pc%IG_HeatCapacity%D = 2.1150E05;
!!$    pc%IG_HeatCapacity%E = -7.1797E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 9.8902E04; 1500.00, 3.4133E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.1875E-01;
!!$    pc%SecondViralCoef%B = -2.1461E02;
!!$    pc%SecondViralCoef%C = -5.9275E07;
!!$    pc%SecondViralCoef%D = -1.3203E20;
!!$    pc%SecondViralCoef%E = 9.8814E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 294.05, -4.6034E00; 1500.00, 5.8109E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.2758E01;
!!$    pc%liqViscosity%B = 2.9169E03;
!!$    pc%liqViscosity%C = 1.2839E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 253.15, 1.6075E-02; 420.00, 3.1659E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.8903E-07;
!!$    pc%vapViscosity%B = 7.0310E-01;
!!$    pc%vapViscosity%C = 1.7590E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 410.95, 9.1106E-06; 1000.00, 2.0676E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0060E-01;
!!$    pc%liqThermalConductivity%B = -1.6030E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 273.15, 1.5681E-01; 353.15, 1.4399E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.8960E03;
!!$    pc%vapThermalConductivity%B = 8.9850E-01;
!!$    pc%vapThermalConductivity%C = 1.2735E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 410.95, 2.0849E-02; 990.95, 1.1087E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0660E-02;
!!$    pc%SurfaceTension%B = -8.3920E-05;
!!$    pc%SurfaceTension%C = -3.8300E-09;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 283.15, 2.6591E-02; 373.15, 1.8812E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_1pentanol
!!$
!!$
!!$  subroutine pc_1propanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1propanol','1-propanol','1Ppropanol','n-propanol','n-Propanol', 'propanol-1'}; 
!!$    pc%CAS = '71-23-8';
!!$    ! components
!!$    pc%ChemicalFormula%C = 3;         
!!$    pc%ChemicalFormula%H = 8;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 60.096;   ! kg/kmol    
!!$    pc%Tcrit = 536.78;             ! K    
!!$    pc%pcrit = 5.1750E+06;         ! Pa 
!!$    pc%vcrit = 0.21900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.254; ! -
!!$
!!$
!!$    pc%MeltingPoint = 146.95;              ! K
!!$    pc%TripletPointT = 146.95;             ! K
!!$    pc%TripletPointp = 3.0828E-07;         ! Pa
!!$    pc%NormalBoilingPoint = 370.35;        ! K
!!$    pc%liqMolarVolume = 0.075089;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.5520E+08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5990E+08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.2247E+05;            ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 5.3720E+06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.8438E+09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK
!!$
!!$
!!$    pc%AcentricFactor = 0.6218;            ! -
!!$    pc%RadiusOfGyration = 2.359E-10;       ! m
!!$    pc%SolubilityParameter = 2.4450E+04;   ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.6000E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.04217;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 6.2800E+08;       ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3837;           ! -
!!$    pc%FlashPoint = 288.15;                ! K
!!$    pc%FlammabilityLimits = (/ 2.0 12.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 644.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.7800E+01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 146.95, 1.7800E+01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.2350E+00;
!!$    pc%liqDensity%B = 2.7136E-01;
!!$    pc%liqDensity%C = 5.3678E+02;
!!$    pc%liqDensity%D = 2.4000E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 146.95, 1.5231E+01; 536.78, 4.5511E+00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  8.8134E+01;
!!$    pc%vapPressure%B = -8.4986E+03;
!!$    pc%vapPressure%C = -9.0766E+00;
!!$    pc%vapPressure%D =  8.3303E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 146.95, 3.0828E-07; 536.78, 5.1214E+06 /);
!!$    pc%vapPressure%eqn = 101;
!!$
!!$    pc%HeatOfVap%A = 6.3300E+07;
!!$    pc%HeatOfVap%B = 3.5750E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 146.95, 5.6460E+07; 536.78, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 8.7700E+03;
!!$    pc%solHeatCapacity%B = 8.1580E+02;
!!$    pc%solHeatCapacity%C = 2.6180E+00;
!!$    pc%solHeatCapacity%D = -7.6720E-02;
!!$    pc%solHeatCapacity%E = 3.3530E-04;
!!$    pc%solHeatCapacity%range(:,1) = (/ 15.00, 3.8141E+03; 140.00, 7.5044E+04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.5876E+05;
!!$    pc%liqHeatCapacity%B = -6.3500E+02;
!!$    pc%liqHeatCapacity%C = 1.9690E+00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 146.95, 1.0797E+05; 400.00 2.1980E+05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 6.1900E+04;
!!$    pc%IG_HeatCapacity%B = 2.0213E+05;
!!$    pc%IG_HeatCapacity%C = 1.6293E+03;
!!$    pc%IG_HeatCapacity%D = 1.2956E+05;
!!$    pc%IG_HeatCapacity%E = 7.2740E+02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 6.6651E+04; 1500.00, 2.2458E+05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.6200E-02;
!!$    pc%SecondViralCoef%B = -5.5600E+01;
!!$    pc%SecondViralCoef%C = -3.3370E+07;
!!$    pc%SecondViralCoef%D = -3.2000E+18;
!!$    pc%SecondViralCoef%E = -2.5610E+22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 295.19, -2.9998E+00; 1500.00, -7.5487E-04 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.7528E02;
!!$    pc%liqViscosity%B = 6.1044E03;
!!$    pc%liqViscosity%C = 2.8761E01;
!!$    pc%liqViscosity%D = -5.1323E-02;
!!$    pc%liqViscosity%E = 1.0000E00;
!!$    pc%liqViscosity%range(:,1) = (/ 146.95, 9.3820E00; 370.35, 4.5317E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 7.9420E-07;
!!$    pc%vapViscosity%B = 5.4910E-01;
!!$    pc%vapViscosity%C = 4.1580E+02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 200.00, 4.7317E-06; 1000.00, 2.4902E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2030E-01;
!!$    pc%liqThermalConductivity%B = -2.1550E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 200.00, 1.7720E-01; 370.35, 1.4049E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -6.1384E+02;
!!$    pc%vapThermalConductivity%B = 7.9270E-01;
!!$    pc%vapThermalConductivity%C = -1.1574E+09;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 372.13, 2.1531E-02; 720.25, 7.0337E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 4.5330E-02;
!!$    pc%SurfaceTension%B = -6.8800E-05;
!!$    pc%SurfaceTension%C = -1.6000E-08;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 283.15, 2.4566E-02; 370.35, 1.7655E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_1propanol
!!$
!!$
!!$  subroutine pc_1tetradecanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1tetradecanol','n-tetradecanol','n-Tetradecanol','1Tetradecanol','1-tetradecanol','1-Tetradecanol'}; 
!!$    pc%CAS = '112-72-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 14;         
!!$    pc%ChemicalFormula%H = 30;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 214.392;  ! kg/kmol    
!!$    pc%Tcrit = 747.00;             ! K    
!!$    pc%pcrit = 1.8100E06;          ! Pa 
!!$    pc%vcrit = 0.80200;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.234; ! -
!!$
!!$
!!$    pc%MeltingPoint = 310.65;              ! K
!!$    pc%TripletPointT = 310.65;             ! K
!!$    pc%TripletPointp = 7.9998E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 568.80;        ! K
!!$    pc%liqMolarVolume = 0.260434;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -4.7590E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -6.2200E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  7.5390E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.5100E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -8.4914E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.7432;            ! -
!!$    pc%RadiusOfGyration = 6.7300E-10;      ! m
!!$    pc%SolubilityParameter = 1.9430E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.1702E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.15470;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 2.1130E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4454;           ! -
!!$    pc%FlashPoint = 414.15;                ! K
!!$    pc%FlammabilityLimits = (/ 0.5 5.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 550.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 4.4513E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 294.00, 4.4513E00; 298.15, 4.4513E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.2737E-01;
!!$    pc%liqDensity%B = 2.6236E-01;
!!$    pc%liqDensity%C = 7.4700E02;
!!$    pc%liqDensity%D = 3.2417E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 310.065, 3.8397E00; 747.00, 1.2478E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  4.1146E02;
!!$    pc%vapPressure%B = -2.6293E04;
!!$    pc%vapPressure%C = -5.9876E01;
!!$    pc%vapPressure%D =  4.5908E-02;
!!$    pc%vapPressure%E =  1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 310.65, 7.9998E-02; 747.00, 1.8151E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 2.6345E08;
!!$    pc%HeatOfVap%B = 2.6157E00;
!!$    pc%HeatOfVap%C = -1.9979E00;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 310.65, 1.0092E08; 747.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -4.4453E05;
!!$    pc%solHeatCapacity%B = 2.9100E03;
!!$    pc%solHeatCapacity%C = -8.0430E00;
!!$    pc%solHeatCapacity%D = 2.3700E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 269.00, 3.3826E05; 310.65, 4.5946E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = -4.4041E06;
!!$    pc%liqHeatCapacity%B = 3.7646E04;
!!$    pc%liqHeatCapacity%C = -1.0306E02;
!!$    pc%liqHeatCapacity%D = 1.2078E-01;
!!$    pc%liqHeatCapacity%E = -4.9723E-05;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 324.00, 5.3443E05; 568.80, 6.8753E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.3340E05;
!!$    pc%IG_HeatCapacity%B = 7.2700E05;
!!$    pc%IG_HeatCapacity%C = 1.5623E03;
!!$    pc%IG_HeatCapacity%D = 5.4210E05;
!!$    pc%IG_HeatCapacity%E = 7.2300E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 2.3668E05; 1500.00, 8.4672E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 7.2937E-01;
!!$    pc%SecondViralCoef%B = - 9.0895E02;
!!$    pc%SecondViralCoef%C = -2.3210E08;
!!$    pc%SecondViralCoef%D = 1.0921E22;
!!$    pc%SecondViralCoef%E = -4.4086E24;
!!$    pc%SecondViralCoef%range(:,1) = (/ 373.50, -8.4889E00; 1500.00, 5.4944E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1007E01;
!!$    pc%liqViscosity%B = 3.5413E03;
!!$    pc%liqViscosity%C = -8.1269E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 310.65, 1.3966E-02; 568.80, 4.8381E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 9.3882E-07;
!!$    pc%vapViscosity%B = 4.7630E-01;
!!$    pc%vapViscosity%C = 8.4022E02;
!!$    pc%vapViscosity%D = -5.3350E04;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 310.65, 4.5823E-06; 1000.00, 1.4106E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.3150E-01;
!!$    pc%liqThermalConductivity%B = -1.9870E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 310.65, 1.6977E-01; 568.80, 1.1848E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -3.0040E00;
!!$    pc%vapThermalConductivity%B = -3.4874E-02;
!!$    pc%vapThermalConductivity%C = - 1.7663E03;
!!$    pc%vapThermalConductivity%D = -2.8478E07;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 568.80, 2.6716E-02; 1000.0, 8.0730E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.9479E-02;
!!$    pc%SurfaceTension%B = 1.3005E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 310.65, 2.9561E-02; 747.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1tetradecanol
!!$
!!$
!!$  subroutine pc_1tridecanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1tridecanol','n-tridecanol','n-Tridecanol','1-Tridecanol'}; 
!!$    pc%CAS = '112-70-9';
!!$    ! components
!!$    pc%ChemicalFormula%C = 13;         
!!$    pc%ChemicalFormula%H = 28;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 200.365;  ! kg/kmol    
!!$    pc%Tcrit = 734.00;             ! K    
!!$    pc%pcrit = 1.9350E06;          ! Pa 
!!$    pc%vcrit = 0.74900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.237; ! -
!!$
!!$
!!$    pc%MeltingPoint = 303.75;              ! K
!!$    pc%TripletPointT = 303.75;             ! K
!!$    pc%TripletPointp = 9.4584E-02;          ! Pa
!!$    pc%NormalBoilingPoint = 553.60;        ! K
!!$    pc%liqMolarVolume = 0.241531;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -4.6040E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -7.5660E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  7.1470E05;            ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.3300E07;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -7.9014E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.7124;            ! -
!!$    pc%RadiusOfGyration = 6.4170E-10;      ! m
!!$    pc%SolubilityParameter = 1.9760E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.500E-30;           ! C^m
!!$    pc%VanDerWaalsVolume = 0.14447;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.9780E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4433;           ! -
!!$    pc%FlashPoint = 394.15;                ! K
!!$    pc%FlammabilityLimits = (/ 0.6 5.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 550.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 4.7913E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 294.00, 4.7213E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.5093E-01;
!!$    pc%liqDensity%B = 2.6286E-01;
!!$    pc%liqDensity%C = 7.3400E02;
!!$    pc%liqDensity%D = 3.1074E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 303.75, 4.1402E00; 734.00, 1.3320E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  2.6599E02;
!!$    pc%vapPressure%B = -2.1388E04;
!!$    pc%vapPressure%C = -3.4855E01;
!!$    pc%vapPressure%D =  1.4125E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 303.75, 9.4584E-02; 734.00, 1.9240E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 2.3170E08;
!!$    pc%HeatOfVap%B = 2.3543E00;
!!$    pc%HeatOfVap%C = -1.7420E00;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 303.75, 9.6831E07; 734.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 3.8300E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 293.15, 3.8300E08; 293.15, 3.8300E08 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.0247E06;
!!$    pc%liqHeatCapacity%B = -3.2237E04;
!!$    pc%liqHeatCapacity%C = 1.3859E02;
!!$    pc%liqHeatCapacity%D = -2.4455E-01;
!!$    pc%liqHeatCapacity%E = 1.5407E-04;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 303.75, 4.7757E05; 553.60, 6.3232E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.1780E05;
!!$    pc%IG_HeatCapacity%B = 6.8080E05;
!!$    pc%IG_HeatCapacity%C = 1.5702E03;
!!$    pc%IG_HeatCapacity%D = 5.0750E05;
!!$    pc%IG_HeatCapacity%E = 7.2440E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 2.2082E05; 1500.00, 7.2076E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.4000E-02;
!!$    pc%SecondViralCoef%B = -5.5700E01;
!!$    pc%SecondViralCoef%C = -1.2900E07;
!!$    pc%SecondViralCoef%D = -6.4900E19;
!!$    pc%SecondViralCoef%E = -2.5800E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 367.00, -8.2103E00; 1500.00, 5.3664E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.6661E01;
!!$    pc%liqViscosity%B = 3.7541E03;
!!$    pc%liqViscosity%C = 1.6684E-02;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 303.75, 1.4899E-02; 553.60, 5.6891E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.2436E-08;
!!$    pc%vapViscosity%B = 8.8190E-01;
!!$    pc%vapViscosity%C = -4.0200E00;
!!$    pc%vapViscosity%D = 2.7480E04;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 303.75, 4.3038E-06; 1000.00, 1.4531E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1542E-01;
!!$    pc%liqThermalConductivity%B = -2.1522E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 303.75, 1.5005E-01; 553.60, 9.6274E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 3.4379E03;
!!$    pc%vapThermalConductivity%B = 9.5884E-01;
!!$    pc%vapThermalConductivity%C = 3.0868E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 553.60, 2.6319E-02; 1000.0, 8.3812E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.9791E-02;
!!$    pc%SurfaceTension%B = 1.2646E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 303.75, 3.0428E-02; 734.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1tridecanol
!!$
!!$
!!$  subroutine pc_1undecanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'1undecanol','1-undecanol','1-Undecanol','1Undecanol','n-Undecanol','n-undecanol'}; 
!!$    pc%CAS = '112-42-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 11;         
!!$    pc%ChemicalFormula%H = 24;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 172.311;   ! kg/kmol    
!!$    pc%Tcrit = 703.60;             ! K    
!!$    pc%pcrit = 2.1470E06;           ! Pa 
!!$    pc%vcrit = 0.71800;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.264; ! -
!!$
!!$
!!$    pc%MeltingPoint = 288.45;              ! K
!!$    pc%TripletPointT = 288.45;             ! K
!!$    pc%TripletPointp = 1.2889E-01;          ! Pa
!!$    pc%NormalBoilingPoint = 518.15;        ! K
!!$    pc%liqMolarVolume = 0.207672;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -4.1840E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -9.1530E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  6.3620E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.7110E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -6.7261E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.2220E08;              ! J/kmol
!!$    pc%StandardHeatOfFormation = -5.0480E08;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 4.4920E05;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.6247;            ! -
!!$    pc%RadiusOfGyration = 5.8080E-10;      ! m
!!$    pc%SolubilityParameter = 1.9920E04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.12401;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.7080E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4386;           ! -
!!$    pc%FlashPoint = 383.15;                ! K
!!$    pc%FlammabilityLimits = (/ 0.7 5.2 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 387 428 /); ! K
!!$    pc%AutoignitionTemp = 538.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 5.6800E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 288.45, 5.6800E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.3013E-01;
!!$    pc%liqDensity%B = 2.3703E-01;
!!$    pc%liqDensity%C = 7.0360E02;
!!$    pc%liqDensity%D = 2.7001E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 288.45, 4.8537E00; 703.60, 1.3928E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  2.2521E02;
!!$    pc%vapPressure%B = -1.8532E04;
!!$    pc%vapPressure%C = -2.8939E01;
!!$    pc%vapPressure%D =  1.0981E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 288.45, 1.2889E-01; 703.60, 2.1419E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 1.2015E08;
!!$    pc%HeatOfVap%B = 6.3047E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 288.46, 8.6152E07; 703.60, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 3.3000E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 289.05, 3.3000E05; 289.05, 3.3000E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.2945E05;
!!$    pc%liqHeatCapacity%B = -3.0395E03;
!!$    pc%liqHeatCapacity%C = 2.7927E01;
!!$    pc%liqHeatCapacity%D = -6.1847E-02;
!!$    pc%liqHeatCapacity%E = 4.3042E-05;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 289.05, 3.9103E05; 520.00, 5.5125E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.8590E05;
!!$    pc%IG_HeatCapacity%B = 5.8690E05;
!!$    pc%IG_HeatCapacity%C = 1.5718E03;
!!$    pc%IG_HeatCapacity%D = 4.3260E05;
!!$    pc%IG_HeatCapacity%E = 7.2270E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 2.0232E05; 1500.00, 6.7834E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 5.5925E-01;
!!$    pc%SecondViralCoef%B = -7.1346E02;
!!$    pc%SecondViralCoef%C = -1.2861E08;
!!$    pc%SecondViralCoef%D = 5.6156E21;
!!$    pc%SecondViralCoef%E = -2.1970E24;
!!$    pc%SecondViralCoef%range(:,1) = (/ 351.80, -7.1054E00; 1500.00, 4.5665E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -7.3477E01;
!!$    pc%liqViscosity%B = 5.9928E03;
!!$    pc%liqViscosity%C = 8.6131E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 288.45, 2.0000E-02; 590.15, 2.3230E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 5.9537E-08;
!!$    pc%vapViscosity%B = 8.1842E-01;
!!$    pc%vapViscosity%C = 9.0245E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 288.45, 4.6767E-06; 1000.00, 1.5578E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1211E-01;
!!$    pc%liqThermalConductivity%B = -2.1815E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 288.45, 1.4918E-01; 561.20, 8.9684E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.4988E03;
!!$    pc%vapThermalConductivity%B = 9.5209E-01;
!!$    pc%vapThermalConductivity%C = 2.0167E10;
!!$    pc%vapThermalConductivity%D = -2.6800E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 518.15, 2.4658E-02; 1000.0, 8.8994E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.3432E-02;
!!$    pc%SurfaceTension%B = 1.1218E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 288.45, 2.9565E-02; 703.60, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_1undecanol
!!$
!!$
!!$  subroutine pc_2butanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2butanol','2-butanol','2-Butanol','2Butanol','Butane'}; 
!!$    pc%CAS = '78-92-2';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 74.123;   ! kg/kmol    
!!$    pc%Tcrit = 536.05;             ! K    
!!$    pc%pcrit = 4.1790E06;           ! Pa 
!!$    pc%vcrit = 0.26900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.252; ! -
!!$
!!$
!!$    pc%MeltingPoint = 158.45;              ! K
!!$    pc%TripletPointT = 158.45;             ! K
!!$    pc%TripletPointp = 1.1323E-06;          ! Pa
!!$    pc%NormalBoilingPoint = 372.70;        ! K
!!$    pc%liqMolarVolume = 0.092017;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.9290E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.6960E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.6469E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  5.9710E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.4408E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5722;            ! -
!!$    pc%RadiusOfGyration = 3.2030E-10;      ! m
!!$    pc%SolubilityParameter = 2.2540E04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05239;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.6200E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3949;           ! -
!!$    pc%FlashPoint = 297.00;                ! K
!!$    pc%FlammabilityLimits = (/ 1.7 9.8 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 679.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.4700E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 158.45, 1.4700E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.6600E-01;
!!$    pc%liqDensity%B = 2.6064E-01;
!!$    pc%liqDensity%C = 5.3605E02;
!!$    pc%liqDensity%D = 2.7460E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 158.45, 1.2570E01; 536.05, 3.7063E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A = 1.5254E02;
!!$    pc%vapPressure%B = -1.1111E04;
!!$    pc%vapPressure%C = -1.9025E01;
!!$    pc%vapPressure%D = 1.0426E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 158.45, 1.1323E-06; 536.05, 4.2014E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 7.2560E07;
!!$    pc%HeatOfVap%B = 4.7740E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 158.45, 6.1383E07; 536.05, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.0530E04;
!!$    pc%solHeatCapacity%B = 1.2313E03;
!!$    pc%solHeatCapacity%C = -6.8730E00;
!!$    pc%solHeatCapacity%D = 1.9800E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 15.00, 6.4599E03; 150.00, 8.6347E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.0670E05;
!!$    pc%liqHeatCapacity%B = -1.0204E03;
!!$    pc%liqHeatCapacity%C = 3.2900E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 158.45, 1.2762E05; 372.70, 2.8340E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 8.2020E04;
!!$    pc%IG_HeatCapacity%B = 2.5220E05;
!!$    pc%IG_HeatCapacity%C = 1.6010E03;
!!$    pc%IG_HeatCapacity%D = 1.5864E05;
!!$    pc%IG_HeatCapacity%E = -7.0415E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 8.8896E04; 1500.00, 2.8513E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.8777E-01;
!!$    pc%SecondViralCoef%B = -1.7264E02;
!!$    pc%SecondViralCoef%C = -3.6937E07;
!!$    pc%SecondViralCoef%D = -5.7699E19;
!!$    pc%SecondViralCoef%E = 4.7178E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 268.03, -3.8800E00; 1500.00, 6.1730E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.0638E02;
!!$    pc%liqViscosity%B = 7.4343E03;
!!$    pc%liqViscosity%C = 1.3285E01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 288.15, 4.7932E-03; 373.15, 4.1659E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.2534E-07;
!!$    pc%vapViscosity%B = 7.6560E-01;
!!$    pc%vapViscosity%C = 9.2847E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 158.45, 3.8197E-06; 1000.00, 2.2716E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0460E-01;
!!$    pc%liqThermalConductivity%B = -2.3542E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 200.00, 1.5752E-01; 398.15, 1.1087E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.0086E03;
!!$    pc%vapThermalConductivity%B = 9.7450E-01;
!!$    pc%vapThermalConductivity%C = 1.5008E10;
!!$    pc%vapThermalConductivity%D = -9.2508E11;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 372.70, 1.9153E-02; 1000.0, 1.1959E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 4.6450E-02;
!!$    pc%SurfaceTension%B = -7.4730E-05;
!!$    pc%SurfaceTension%C = -1.3000E-08;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 273.15, 2.5068E-02; 372.70, 1.6792E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_2butanol
!!$
!!$
!!$  subroutine pc_2butanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2butanone','2-butanone','2-Butanone','2-butanone','butanone'}; 
!!$    pc%CAS = '78-93-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 2;         
!!$    pc%ChemicalFormula%H = 6;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 72.107;   ! kg/kmol    
!!$    pc%Tcrit = 535.50;             ! K    
!!$    pc%pcrit = 4.1543E06;           ! Pa 
!!$    pc%vcrit = 0.26700;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.249; ! -
!!$
!!$
!!$    pc%MeltingPoint = 186.48;              ! K
!!$    pc%TripletPointT = 186.48;             ! K
!!$    pc%TripletPointp = 1.3941E00;          ! Pa
!!$    pc%NormalBoilingPoint = 352.79;        ! K
!!$    pc%liqMolarVolume = 0.090204;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.3495E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.4606E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.3811E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  8.4391E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.2616E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3241;            ! -
!!$    pc%RadiusOfGyration = 3.1350E-10;      ! m
!!$    pc%SolubilityParameter = 1.8796E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.2064E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.04927;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.1900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3764;           ! -
!!$    pc%FlashPoint = 267.04;                ! K
!!$    pc%FlammabilityLimits = (/ 1.8 10.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 788.71;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.4700E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 186.48, 1.4700E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.3300E-01;
!!$    pc%liqDensity%B = 2.514E-01;
!!$    pc%liqDensity%C = 5.3550E02;
!!$    pc%liqDensity%D = 2.8570E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 186.48, 1.2593E01; 535.50, 3.7112E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.1474E02;
!!$    pc%vapPressure%B = -7.1300E03;
!!$    pc%vapPressure%C = -1.5184E01;
!!$    pc%vapPressure%D =  1.7234E-02;
!!$    pc%vapPressure%E =  1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 186.48, 1.3941E00; 535.50, 4.1976E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.1400E07;
!!$    pc%HeatOfVap%B = 7.5450E-01;
!!$    pc%HeatOfVap%C = -5.9500E-01;
!!$    pc%HeatOfVap%D = 2.3300E-01;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 186.48, 4.0176E07; 535.50, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -9.712E03;
!!$    pc%solHeatCapacity%B = 1.3640E03;
!!$    pc%solHeatCapacity%C = -7.1303E00;
!!$    pc%solHeatCapacity%D = 1.6570E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 10.00, 3.2315E03; 186.48, 1.0415E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.3230E05;
!!$    pc%liqHeatCapacity%B = 2.0087E02;
!!$    pc%liqHeatCapacity%C = -9.5970E-01;
!!$    pc%liqHeatCapacity%D = 1.9533E-03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 186.48, 1.4905E05; 373.15, 1.7511E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 7.8400E04;
!!$    pc%IG_HeatCapacity%B = 2.1032E05;
!!$    pc%IG_HeatCapacity%C = 1.5488E03;
!!$    pc%IG_HeatCapacity%D = 1.1855E05;
!!$    pc%IG_HeatCapacity%E = 6.9300E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 8.3967E04; 1500.00, 2.4816E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 9.0100E-02;
!!$    pc%SecondViralCoef%B = -8.6700E01;
!!$    pc%SecondViralCoef%C = -6.0500E07;
!!$    pc%SecondViralCoef%D = 3.4650E20;
!!$    pc%SecondViralCoef%E = -1.0330E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 277.75, -3.7632E00; 2677.50, 5.4567E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -8.1340E-01;
!!$    pc%liqViscosity%B = 5.0978E02;
!!$    pc%liqViscosity%C = -1.5324E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 188.71, 2.1504E-03; 533.15, 7.6450E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.9759E-08;
!!$    pc%vapViscosity%B = 1.0252E00;
!!$    pc%vapViscosity%C = -1.6480E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 273.15, 6.6159E-06; 573.15, 1.3684E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1920E-01;
!!$    pc%liqThermalConductivity%B = -2.4840E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 186.48, 1.7288E-01; 352.79, 1.3157E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -4.2027E06;
!!$    pc%vapThermalConductivity%B = -1.5240E-01;
!!$    pc%vapThermalConductivity%C = 2.0846E09;
!!$    pc%vapThermalConductivity%D = -1.4577E13;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 352.79, 1.5457E-02; 992.79, 1.1571E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.9441E-02;
!!$    pc%SurfaceTension%B = 1.1165E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 186.48, 2.4101E-02; 355.50, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2butanone
!!$
!!$
!!$  subroutine pc_2heptanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2heptanone','2-heptanone','2Heptanone','2-Heptanone'}; 
!!$    pc%CAS = '110-43-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 7;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 114.188;   ! kg/kmol    
!!$    pc%Tcrit = 611.40;             ! K    
!!$    pc%pcrit = 2.9400E06;           ! Pa 
!!$    pc%vcrit = 0.43400;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.251; ! -
!!$
!!$
!!$    pc%MeltingPoint = 238.15;              ! K
!!$    pc%TripletPointT = 238.15;             ! K
!!$    pc%TripletPointp = 3.5448E00;          ! Pa
!!$    pc%NormalBoilingPoint = 424.18;        ! K
!!$    pc%liqMolarVolume = 0.140675;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.0045E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.1960E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.8064E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.8400E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.0995E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.3480E08;              ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.4820E08;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.4100E05;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4190;            ! -
!!$    pc%RadiusOfGyration = 4.4980E-10;      ! m
!!$    pc%SolubilityParameter = 1.7800E04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 8.7100E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.07996;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.1240E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4066;           ! -
!!$    pc%FlashPoint = 312.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.1 7.9 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 311 350 /); ! K
!!$    pc%AutoignitionTemp = 666.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 8.8400E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 238.15, 8.8400E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.8247E-01;
!!$    pc%liqDensity%B = 2.5279E-01;
!!$    pc%liqDensity%C = 6.1140E02;
!!$    pc%liqDensity%D = 2.9818E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 238.15, 7.5514E00; 611.40, 2.3042E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.5494E01;
!!$    pc%vapPressure%B = -7.8965E03;
!!$    pc%vapPressure%C = -7.5047E00;
!!$    pc%vapPressure%D =  8.9130E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 238.15, 3.5448E00; 611.40, 2.9456E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.1425E07;
!!$    pc%HeatOfVap%B = 3.9802E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 238.15, 5.0471E07; 611.40, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 2.0390E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 238.15, 2.0390E05; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.6504E05;
!!$    pc%liqHeatCapacity%B = -3.7568E02;
!!$    pc%liqHeatCapacity%C = 1.0024E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 238.15, 2.3242E05; 490.00, 3.2163E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2507E05;
!!$    pc%IG_HeatCapacity%B = 2.1480E05;
!!$    pc%IG_HeatCapacity%C = 6.9120E02;
!!$    pc%IG_HeatCapacity%D = 1.6190E05;
!!$    pc%IG_HeatCapacity%E = 1.7593E03;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 1.2688E05; 1200.00, 3.8446E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.9183E-01;
!!$    pc%SecondViralCoef%B = -3.4109E02;
!!$    pc%SecondViralCoef%C = -7.8913E07;
!!$    pc%SecondViralCoef%D = -4.0026E20;
!!$    pc%SecondViralCoef%E = 5.6762E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 305.70, -6.3996E00; 1500.00, 4.1041E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.3929E01;
!!$    pc%liqViscosity%B = 1.3219E03;
!!$    pc%liqViscosity%C = 4.0382E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 250.00, 1.6422E-03; 424.18, 2.3184E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 8.8629E-08;
!!$    pc%vapViscosity%B = 7.8376E-01;
!!$    pc%vapViscosity%C = 1.0018E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 238.15, 4.5495E-06; 1000.00, 1.8088E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1080E-01;
!!$    pc%liqThermalConductivity%B = -2.4600E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 236.15, 1.5222E-01; 424.05, 1.0648E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.0493E03;
!!$    pc%vapThermalConductivity%B = 1.0323E00;
!!$    pc%vapThermalConductivity%C = 2.2983E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 424.18, 1.9506E-02; 1000.0, 1.1145E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.1492E-02;
!!$    pc%SurfaceTension%B = 1.2802E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 236.15, 3.2692E-02; 611.40, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2heptanone
!!$
!!$
!!$  subroutine pc_2hexanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2hexanol','2-hexanol','2-Hexanol','2Hexanol'}; 
!!$    pc%CAS = '626-93-7';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 102.77;   ! kg/kmol    
!!$    pc%Tcrit = 586.20;             ! K    
!!$    pc%pcrit = 3.4000E06;          ! Pa 
!!$    pc%vcrit = 0.38000;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.265; ! -
!!$
!!$
!!$    pc%MeltingPoint = 223.00;              ! K
!!$    pc%TripletPointT = 223.00;             ! K
!!$    pc%TripletPointp = 9.4925E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 413.04;        ! K
!!$    pc%liqMolarVolume = 0.126073;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.3350E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5020E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 4.3600E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.6660E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5655;            ! -
!!$    pc%RadiusOfGyration = 3.7870E-10;      ! m
!!$    pc%SolubilityParameter = 2.1115E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = MAX_REAL_WP;                 ! C^m
!!$    pc%VanDerWaalsVolume = 0.07285;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.0320E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4128;           ! -
!!$    pc%FlashPoint = 314.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.2 8.3 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = MAX_REAL_WP;             ! K
!!$
!!$
!!$    pc%solDensity%A = 9.9600E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 223.00, 9.9600E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.9120E-01;
!!$    pc%liqDensity%B = 2.6275E-01;
!!$    pc%liqDensity%C = 5.8620E02;
!!$    pc%liqDensity%D = 2.6945E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 223.00, 8.5168E00; 586.20, 2.6306E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.2842E02;
!!$    pc%vapPressure%B = -1.0914E04;
!!$    pc%vapPressure%C = -1.5189E01;
!!$    pc%vapPressure%D =  5.9636E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 223.00, 9.4925E-02; 586.20, 3.4076E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 9.1260E07;
!!$    pc%HeatOfVap%B = 6.2140E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 223.00, 6.7778E07; 586.20, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 1.9600E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 223.00, 1.9600E05; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.6063E05;
!!$    pc%liqHeatCapacity%B = -1.9047E03;
!!$    pc%liqHeatCapacity%C = 5.2540E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 230.00, 2.0049E05; 330.00, 3.0424E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2575E05;
!!$    pc%IG_HeatCapacity%B = 3.5890E05;
!!$    pc%IG_HeatCapacity%C = 1.8215E03;
!!$    pc%IG_HeatCapacity%D = 2.5557E05;
!!$    pc%IG_HeatCapacity%E = 8.1300E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 300.00, 1.5898E05; 1200.00, 3.7756E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.1711E-01;
!!$    pc%SecondViralCoef%B = -7.5358E01;
!!$    pc%SecondViralCoef%C = -1.0028E08;
!!$    pc%SecondViralCoef%D = 9.8894E20;
!!$    pc%SecondViralCoef%E = -4.2519E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 293.10, -1.2600E01; 2923.10, 8.7315E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.2444E02;
!!$    pc%liqViscosity%B = 8.3425E03;
!!$    pc%liqViscosity%C = 1.5972E01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 230.00, 2.6948E-01; 413.04, 3.2420E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.0200E-07;
!!$    pc%vapViscosity%B = 7.7710E-01;
!!$    pc%vapViscosity%C = 1.0050E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 223.00, 4.6978E-06; 1000.00, 1.9875E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1560E-01;
!!$    pc%liqThermalConductivity%B = -2.6510E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 223.00, 1.5648E-01; 463.00, 9.2859E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -7.9554E06;
!!$    pc%vapThermalConductivity%B = -1.6390E-01;
!!$    pc%vapThermalConductivity%C = 2.4484E09;
!!$    pc%vapThermalConductivity%D = -2.4038E13;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 413.04, 2.1961E-02; 1000.0, 1.1877E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1290E-02;
!!$    pc%SurfaceTension%B = 1.0435E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 223.00, 3.1124E-02; 586.20, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2hexanol
!!$
!!$
!!$  subroutine pc_2hexanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2hexanone','2-hexanone','2-Hexanone','2Hexanone'}; 
!!$    pc%CAS = '591-78-6';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 100.161;   ! kg/kmol    
!!$    pc%Tcrit = 587.61;             ! K    
!!$    pc%pcrit = 3.2870E06;          ! Pa 
!!$    pc%vcrit = 0.37800;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.254; ! -
!!$
!!$
!!$    pc%MeltingPoint = 217.35;              ! K
!!$    pc%TripletPointT = 217.35;             ! K
!!$    pc%TripletPointp = 1.4456E00;          ! Pa
!!$    pc%NormalBoilingPoint = 400.70;        ! K
!!$    pc%liqMolarVolume = 0.124128;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.7983E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.3008E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 4.1786E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.4900E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.4900E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.3960E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.2201E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.0811E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3846;            ! -
!!$    pc%RadiusOfGyration = 4.0930E-10;      ! m
!!$    pc%SolubilityParameter = 1.8050E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 8.9400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06973;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 9.8900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3987;           ! -
!!$    pc%FlashPoint = 296.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.2 8.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 294 330 /); ! K
!!$    pc%AutoignitionTemp = 696.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.0200E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 217.35, 1.0200E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.7816E-01;
!!$    pc%liqDensity%B = 2.5634E-01;
!!$    pc%liqDensity%C = 5.8761E02;
!!$    pc%liqDensity%D = 2.8365E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 217.35, 8.7319E00; 587.61, 2.6455E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.0744E02;
!!$    pc%vapPressure%B = -8.5286E03;
!!$    pc%vapPressure%C = -1.2679E01;
!!$    pc%vapPressure%D =  8.4606E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 217.35, 1.4456E00; 587.61, 3.2859E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.6231E07;
!!$    pc%HeatOfVap%B = 3.8207E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 217.35, 4.7135E07; 587.61, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 5.4305E03;
!!$    pc%solHeatCapacity%B = 8.2129E02;
!!$    pc%solHeatCapacity%C = 5.0967E00;
!!$    pc%solHeatCapacity%D = -7.0876E-02;
!!$    pc%solHeatCapacity%E = 2.1969E-04;
!!$    pc%solHeatCapacity%range(:,1) = (/ 39.76, 4.2236E04; 217.69, 1.8794E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.0825E05;
!!$    pc%liqHeatCapacity%B = -1.0747E02;
!!$    pc%liqHeatCapacity%C = 2.0620E-01;
!!$    pc%liqHeatCapacity%D = 7.0293E-04;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 217.35, 2.0185E05; 460.00, 2.7087E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0940E05;
!!$    pc%IG_HeatCapacity%B = 1.8070E05;
!!$    pc%IG_HeatCapacity%C = 6.8900E02;
!!$    pc%IG_HeatCapacity%D = 1.4740E05;
!!$    pc%IG_HeatCapacity%E = 1.7720E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.1815E05; 1200.00, 3.3207E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.4931E-01;
!!$    pc%SecondViralCoef%B = -2.8704E02;
!!$    pc%SecondViralCoef%C = -5.8288E07;
!!$    pc%SecondViralCoef%D = -2.8661E20;
!!$    pc%SecondViralCoef%E = 4.2633E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 293.80, -5.5749E00; 1500.00, 4.0669E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1445E01;
!!$    pc%liqViscosity%B = 1.1872E03;
!!$    pc%liqViscosity%C = 2.9076E-03;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 217.35, 2.5614E-03; 400.70, 2.1077E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 9.7820E-08;
!!$    pc%vapViscosity%B = 7.7720E-01;
!!$    pc%vapViscosity%C = 9.9530E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 217.35, 4.3968E-06; 1000.00, 1.9091E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1076E-01;
!!$    pc%liqThermalConductivity%B = -2.4000E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 217.35, 1.5860E-01; 400.85, 1.1456E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -1.2158E00;
!!$    pc%vapThermalConductivity%B = 2.6637E-02;
!!$    pc%vapThermalConductivity%C = -1.7116E03;
!!$    pc%vapThermalConductivity%D = -1.3176E07;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.00, 7.7543E-03; 1000.0, 1.0523E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.2446E-02;
!!$    pc%SurfaceTension%B = 1.2674E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 217.35, 3.4777E-02; 587.61, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2hexanone
!!$
!!$
!!$  subroutine pc_2methyl1butanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2methyl1butanol','2methyl-1butanol','amylalcohol'}; 
!!$    pc%CAS = '137-32-6';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 88.150;   ! kg/kmol    
!!$    pc%Tcrit = 575.40;             ! K    
!!$    pc%pcrit = 3.9400E06;           ! Pa 
!!$    pc%vcrit = 0.32700;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.269; ! -
!!$
!!$
!!$    pc%MeltingPoint = 195.00;              ! K
!!$    pc%TripletPointT = 195.00;             ! K
!!$    pc%TripletPointp = 1.2845E-03;         ! Pa
!!$    pc%NormalBoilingPoint = 401.85;        ! K
!!$    pc%liqMolarVolume = 0.108257;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -5.0209E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.4671E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.9351E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  8.2700E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.0620E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.5940E08;              ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.5660E08;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.5330E05;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5736;            ! -
!!$    pc%RadiusOfGyration = 3.6120E-10;      ! m
!!$    pc%SolubilityParameter = 2.2090E+04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 6.2700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06262;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.9700E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4086;           ! -
!!$    pc%FlashPoint = 316.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.4 9.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 313 345 /); ! K
!!$    pc%AutoignitionTemp = 658.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1800E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 195.0, 1.1800E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.1662E-01;
!!$    pc%liqDensity%B = 2.6703E-01;
!!$    pc%liqDensity%C = 5.7540E02;
!!$    pc%liqDensity%D = 2.4336E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 195.00, 1.0092E01; 575.40, 3.0582E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  9.6196E01;
!!$    pc%vapPressure%B = -9.6589E03;
!!$    pc%vapPressure%C = -1.0112E01;
!!$    pc%vapPressure%D =  8.4773E-19;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 195.00, 1.2846E-03; 575.40, 3.9057E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 7.7920E07;
!!$    pc%HeatOfVap%B = 4.6991E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 195.00, 6.4149E07; 575.40, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    !keine Daten in Stoffdatenblatt
!!$    pc%solHeatCapacity%A = 0;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP,MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 8.3025E04;
!!$    pc%liqHeatCapacity%B = 4.5971E02;
!!$    pc%liqHeatCapacity%C = 0;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 195.00, 1.7267E05; 401.85, 2.6776E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0890E05;
!!$    pc%IG_HeatCapacity%B = 2.1850E05;
!!$    pc%IG_HeatCapacity%C = 8.5300E02;
!!$    pc%IG_HeatCapacity%D = 1.4000E05;
!!$    pc%IG_HeatCapacity%E = 2.9060E03;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 298.15, 1.3247E05; 1500.00, 3.4718E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.1048E-01;
!!$    pc%SecondViralCoef%B = -2.0521E02;
!!$    pc%SecondViralCoef%C = -5.4378E07;
!!$    pc%SecondViralCoef%D = -2.0676E20;
!!$    pc%SecondViralCoef%E = 2.8200E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 287.70, -5.1030E00; 1500.00, 5.7554E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -6.4568E01;
!!$    pc%liqViscosity%B = 5.2770E03;
!!$    pc%liqViscosity%C = 7.2659E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 273.15, 1.1253E-02; 403.15, 3.7514E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.0470E-07;
!!$    pc%vapViscosity%B = 7.8178E-01;
!!$    pc%vapViscosity%C = 9.2644E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 195.00, 4.3795E-06; 1000.00, 2.1223E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2553E-01;
!!$    pc%liqThermalConductivity%B = -2.8271E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 195.00, 1.7040E-01; 449.50, 9.8452E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 7.5670E00;
!!$    pc%vapThermalConductivity%B = -2.7146E-01;
!!$    pc%vapThermalConductivity%C = - 3.8945E03;
!!$    pc%vapThermalConductivity%D = 1.2574E07;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 401.85, 2.1482E-02; 1000.0, 1.1986E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.5638E-02;
!!$    pc%SurfaceTension%B = 1.0827E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 195.00, 3.5545E-02; 575.40, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2methyl1butanol
!!$
!!$
!!$  subroutine pc_2methyl1propanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2methyl1propanol','2methyl-1propanol','Fermation Butyl Alcohol','Isobutyl Alcohol'}; 
!!$    pc%CAS = '78-83-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 74.123;   ! kg/kmol    
!!$    pc%Tcrit = 547.78;             ! K    
!!$    pc%pcrit = 4.3000E06;          ! Pa 
!!$    pc%vcrit = 0.27300;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.258; ! -
!!$
!!$
!!$    pc%MeltingPoint = 165.15;              ! K
!!$    pc%TripletPointT = 165.15;             ! K
!!$    pc%TripletPointp = 6.2545E-07;         ! Pa
!!$    pc%NormalBoilingPoint = 380.81;        ! K
!!$    pc%liqMolarVolume = 0.092913;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.8320E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5490E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.4800E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  6.3220E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.4490E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5848;            ! -
!!$    pc%RadiusOfGyration = 3.3320E-10;      ! m
!!$    pc%SolubilityParameter = 2.2910E+04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.4700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05239;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.6200E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3938;           ! -
!!$    pc%FlashPoint = 301.00;                ! K
!!$    pc%FlammabilityLimits = (/ 1.7 10.9 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 681.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.4100E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 165.15, 1.4100E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.9680E-01;
!!$    pc%liqDensity%B = 2.7226E-01;
!!$    pc%liqDensity%C = 5.4778E02;
!!$    pc%liqDensity%D = 2.3890E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 165.15, 1.2084E01; 547.78, 3.6612E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.8779E02;
!!$    pc%vapPressure%B = -1.2955E04;
!!$    pc%vapPressure%C = -2.4285E01;
!!$    pc%vapPressure%D =  1.4262E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 165.15, 6.2544E-07; 547.78, 4.3312E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 7.6190E07;
!!$    pc%HeatOfVap%B = 5.0470E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 165.15, 6.3570E07; 547.78, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.0520E04;
!!$    pc%solHeatCapacity%B = 1.1474E03;
!!$    pc%solHeatCapacity%C = -5.7430E00;
!!$    pc%solHeatCapacity%D = 1.5260E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 15.00, 5.4503E+03; 160.00, 8.8548E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.8020E05;
!!$    pc%liqHeatCapacity%B = -7.3800E02;
!!$    pc%liqHeatCapacity%C = 2.4930E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 165.15, 1.2631E05; 380.81, 2.6069E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 8.7940E04;
!!$    pc%IG_HeatCapacity%B = 2.4160E05;
!!$    pc%IG_HeatCapacity%C = 1.7180E03;
!!$    pc%IG_HeatCapacity%D = 1.6540E05;
!!$    pc%IG_HeatCapacity%E = 7.9870E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 300.00, 1.1089E05; 1200.00, 2.6355E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.4800E-01;
!!$    pc%SecondViralCoef%B = -1.3400E02;
!!$    pc%SecondViralCoef%C = -2.8400E07;
!!$    pc%SecondViralCoef%D = 1.7000E18;
!!$    pc%SecondViralCoef%E = -7.7520E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 301.25, -5.1043E00; 1500.00, 5.0250E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -4.8035E01;
!!$    pc%liqViscosity%B = 4.3067E03;
!!$    pc%liqViscosity%C = 4.8948E00;
!!$    pc%liqViscosity%D = -3.5000E-28;
!!$    pc%liqViscosity%E = 1.0E01;
!!$    pc%liqViscosity%range(:,1) = (/ 210.93, 2.4105E-01; 544.26, 4.1665E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.1945E-07;
!!$    pc%vapViscosity%B = 6.4700E-01;
!!$    pc%vapViscosity%C = 2.2314E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 165.15, 3.6991E-06; 1000.00, 2.2800E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.8503E-01;
!!$    pc%liqThermalConductivity%B = -1.7855E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 165.15, 1.5554E-01; 380.81, 1.1704E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 1.8980E03;
!!$    pc%vapThermalConductivity%B = 8.6980E-01;
!!$    pc%vapThermalConductivity%C = 6.9300E09;
!!$    pc%vapThermalConductivity%D = 1.2574E07;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 380.81, 1.8322E-02; 1000.0, 1.1142E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 4.7506E-02;
!!$    pc%SurfaceTension%B = -8.6704E-05;
!!$    pc%SurfaceTension%C = 1.0150E-08;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 283.15, 2.3770E-02; 380.81, 1.5960E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_2methyl1propanol
!!$
!!$
!!$  subroutine pc_2methyl2propanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2methyl2propanol','2methyl-2propanol','t-Butanol'}; 
!!$    pc%CAS = '75-65-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 2;         
!!$    pc%ChemicalFormula%H = 6;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 74.123;   ! kg/kmol    
!!$    pc%Tcrit = 506.21;             ! K    
!!$    pc%pcrit = 3.9730E6;           ! Pa 
!!$    pc%vcrit = 0.27500;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.260; ! -
!!$
!!$
!!$    pc%MeltingPoint = 298.97;              ! K
!!$    pc%TripletPointT = 298.97;             ! K
!!$    pc%TripletPointp = 5.9356E03;          ! Pa
!!$    pc%NormalBoilingPoint = 355.57;        ! K
!!$    pc%liqMolarVolume = 0.094744;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.1240E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.7760E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.2630E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  6.7030E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.4239E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.6115;            ! -
!!$    pc%RadiusOfGyration = 3.0670E-10;      ! m
!!$    pc%SolubilityParameter = 2.1600E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05238;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.6200E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3852;           ! -
!!$    pc%FlashPoint = 284.26;                ! K
!!$    pc%FlammabilityLimits = (/ 2.4 8.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 751.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.3000E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 298.97, 1.3000E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.2120E-01;
!!$    pc%liqDensity%B = 2.5440E-01;
!!$    pc%liqDensity%C = 5.0621E02;
!!$    pc%liqDensity%D = 2.7600E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 298.97, 1.0555E01; 506.21, 3.6211E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.7231E02;
!!$    pc%vapPressure%B = -1.1590E04;
!!$    pc%vapPressure%C = -2.2118E01;
!!$    pc%vapPressure%D =  1.3709E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 298.97, 5.9356E03; 506.21, 3.9910E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 7.7320E07;
!!$    pc%HeatOfVap%B = 5.6450E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 298.97, 4.6703E07; 506.21, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 4.6000E03;
!!$    pc%solHeatCapacity%B = 4.1200E02;
!!$    pc%solHeatCapacity%C = 4.1750E00;
!!$    pc%solHeatCapacity%D = -3.4270E-02;
!!$    pc%solHeatCapacity%E = 8.0080E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 25.00, 1.7005E04; 286.13, 1.9826E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = -9.2546E05;
!!$    pc%liqHeatCapacity%B = 7.8949E03;
!!$    pc%liqHeatCapacity%C = -1.7661E01;
!!$    pc%liqHeatCapacity%D = 1.3617E-02;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 298.96, 2.2016E05; 460.00, 2.9455E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 7.7040E04;
!!$    pc%IG_HeatCapacity%B = 2.5390E05;
!!$    pc%IG_HeatCapacity%C = 1.5502E03;
!!$    pc%IG_HeatCapacity%D = 1.6690E05;
!!$    pc%IG_HeatCapacity%E = -6.7930E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 8.5670E04; 1500.00, 2.8508E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.8698E-01;
!!$    pc%SecondViralCoef%B = -1.5686E02;
!!$    pc%SecondViralCoef%C = -3.2431E07;
!!$    pc%SecondViralCoef%D = -3.7866E19;
!!$    pc%SecondViralCoef%E =-2.8168E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 253.11, -4.0200E00; 1500.00, 7.2796E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.1640E02;
!!$    pc%liqViscosity%B = 1.3205E04;
!!$    pc%liqViscosity%C = 2.9254E01;
!!$    pc%liqViscosity%D = -2.4616E-27;
!!$    pc%liqViscosity%E = 1.0000E01;
!!$    pc%liqViscosity%range(:,1) = (/ 298.15, 4.3119E-03; 451.21, 1.0137E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 9.6050E-07;
!!$    pc%vapViscosity%B = 4.8560E-01;
!!$    pc%vapViscosity%C = 3.8100E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 298.97, 6.7267E-06; 600.00, 1.3124E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.5960E-01;
!!$    pc%liqThermalConductivity%B = -1.4700E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 293.15, 1.1651E-01; 479.15, 8.9165E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -1.4000E-02;
!!$    pc%vapThermalConductivity%B = 3.6250E-01;
!!$    pc%vapThermalConductivity%C = -2.3760E03;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 355.57, 2.0717E-02; 1000.0, 1.2445E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 4.0190E-02;
!!$    pc%SurfaceTension%B = -5.3100E-05;
!!$    pc%SurfaceTension%C = -5.1200E-08;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 298.15, 1.9807E-02; 355.57, 1.4836E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_2methyl2propanol
!!$
!!$
!!$  subroutine pc_2methylnonane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2methylnonane','2-methylnonane'}; 
!!$    pc%CAS = '871-83-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 10;         
!!$    pc%ChemicalFormula%H = 22;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 142.265;  ! kg/kmol    
!!$    pc%Tcrit = 610.00;             ! K    
!!$    pc%pcrit = 2.1200E06;          ! Pa 
!!$    pc%vcrit = 0.58300;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.244; ! -
!!$
!!$
!!$    pc%MeltingPoint = 198.50;              ! K
!!$    pc%TripletPointT = 198.50;             ! K
!!$    pc%TripletPointp = 2.4868E-03;          ! Pa
!!$    pc%NormalBoilingPoint = 440.15;        ! K
!!$    pc%liqMolarVolume = 0.196780;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.5652E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = 2.8400E07;             ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  5.3802E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.7500E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -6.2888E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4723;            ! -
!!$    pc%RadiusOfGyration = 5.4290E-10;      ! m
!!$    pc%SolubilityParameter = 1.5440E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = MAX_REAL_WP;                 ! C^m
!!$    pc%VanDerWaalsVolume = 0.10917;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.5030E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4075;           ! -
!!$    pc%FlashPoint = 314.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.7 5.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 487.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 6.5300E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 198.50, 6.5300E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 4.5413E-01;
!!$    pc%liqDensity%B = 2.6475E-01;
!!$    pc%liqDensity%C = 6.100E02;
!!$    pc%liqDensity%D = 3.0083E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 198.50, 5.5854E00; 610.00, 1.7153E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.3571E02;
!!$    pc%vapPressure%B = -1.0597E04;
!!$    pc%vapPressure%C = -1.6769E01;
!!$    pc%vapPressure%D =  1.0148E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 198.50, 2.4867E-03; 610.00, 2.1194E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.5700E07;
!!$    pc%HeatOfVap%B = 4.2540E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 198.50, 5.5570E07; 610.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 1.6222E06;
!!$    pc%solHeatCapacity%B = -5.4504E04;
!!$    pc%solHeatCapacity%C = 7.0650E02;
!!$    pc%solHeatCapacity%D = -3.9062E00;
!!$    pc%solHeatCapacity%E = 7.8986E-03;
!!$    pc%solHeatCapacity%range(:,1) = (/ 80.0, 1.0703E05; 198.50, 3.5197E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.0520E05;
!!$    pc%liqHeatCapacity%B = -4.7130E02;
!!$    pc%liqHeatCapacity%C = 1.6330E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 198.50, 2.7599E05; 440.15, 4.1412E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.5390E05;
!!$    pc%IG_HeatCapacity%B = 5.6550E05;
!!$    pc%IG_HeatCapacity%C = 1.5750E03;
!!$    pc%IG_HeatCapacity%D = 3.8870E05;
!!$    pc%IG_HeatCapacity%E = 7.0400E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 1.5677E05; 1500.00, 6.1969E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.4066E-01;
!!$    pc%SecondViralCoef%B = -4.3379E02;
!!$    pc%SecondViralCoef%C = -9.7744E07;
!!$    pc%SecondViralCoef%D = -3.4420E20;
!!$    pc%SecondViralCoef%E = 5.9453E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 305.00, -6.5199E00; 1500.00, 2.2494E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.0613E01;
!!$    pc%liqViscosity%B = 1.7419E03;
!!$    pc%liqViscosity%C = 1.3450E00;
!!$    pc%liqViscosity%D = -1.2703E-28;
!!$    pc%liqViscosity%E = 1E01;
!!$    pc%liqViscosity%range(:,1) = (/ 198.50, 8.9011E-03; 593.15, 5.6989E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 7.7800E-07;
!!$    pc%vapViscosity%B = 4.6360E-01;
!!$    pc%vapViscosity%C = 3.1607E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 198.50, 3.4877E-06; 1000.00, 1.4538E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0262E-01;
!!$    pc%liqThermalConductivity%B = -2.4000E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 198.5, 1.5498E-01; 488.00, 8.5500E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 3.6874E-05;
!!$    pc%vapThermalConductivity%B = 1.1200E00;
!!$    pc%vapThermalConductivity%C = -1.5392E02;
!!$    pc%vapThermalConductivity%D = 1.7277E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 440.15, 2.1849E-02; 1000.0, 8.2911E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1634E-02;
!!$    pc%SurfaceTension%B = 1.2440E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 198.50, 3.1642E-02; 610.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2methylnonane
!!$
!!$
!!$  subroutine pc_2methyloctane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2methyloctane','2-methyloctane'}; 
!!$    pc%CAS = '3221-61 -2';
!!$    ! components
!!$    pc%ChemicalFormula%C = 9;         
!!$    pc%ChemicalFormula%H = 20;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 128.258;  ! kg/kmol    
!!$    pc%Tcrit = 582.80;             ! K    
!!$    pc%pcrit = 2.3100E06;          ! Pa 
!!$    pc%vcrit = 0.54100;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.258; ! -
!!$
!!$
!!$    pc%MeltingPoint = 192.78;              ! K
!!$    pc%TripletPointT = 192.78;             ! K
!!$    pc%TripletPointp = 1.8181E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 416.43;        ! K
!!$    pc%liqMolarVolume = 0.180734;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.3590E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = 2.0080E7;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.9880E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.8000E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -5.6791E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = 8.2190E06;              ! J/kmol
!!$    pc%StandardHeatOfFormation = -2.8082E08;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.8790E05;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4583;            ! -
!!$    pc%RadiusOfGyration = 4.9730E-10;      ! m
!!$    pc%SolubilityParameter = 1.53540E4;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.09894;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.3680E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4008;           ! -
!!$    pc%FlashPoint = 296.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.9 5.4 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 298 333 /); ! K
!!$    pc%AutoignitionTemp = 493.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 7.1600E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 192.78, 7.1600E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 4.8369E-01;
!!$    pc%liqDensity%B = 2.6168E-01;
!!$    pc%liqDensity%C = 5.8280E02;
!!$    pc%liqDensity%D = 2.8064E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 192.78, 6.1230E00; 582.80, 1.8484E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.8095E01;
!!$    pc%vapPressure%B = -7.7507E03;
!!$    pc%vapPressure%C = -7.9631E00;
!!$    pc%vapPressure%D =  1.4327E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 192.78, 1.8180E-02; 582.80, 2.3026E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.9830E07;
!!$    pc%HeatOfVap%B = 3.9749E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 192.78, 5.1002E07; 582.80, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 3.4830E04;
!!$    pc%solHeatCapacity%B = 7.6538E02;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 80.0, 9.6060E04; 160.00, 1.5729E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.0720E05;
!!$    pc%liqHeatCapacity%B = 2.0160E01;
!!$    pc%liqHeatCapacity%C = 7.5925E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 192.78, 2.3930E05; 500.00, 4.0709E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.3743E05;
!!$    pc%IG_HeatCapacity%B = 5.1693E05;
!!$    pc%IG_HeatCapacity%C = 1.5604E03;
!!$    pc%IG_HeatCapacity%D = 3.4726E05;
!!$    pc%IG_HeatCapacity%E = 6.9593E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 1.4022E05; 1500.00, 5.6295E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.0640E-01;
!!$    pc%SecondViralCoef%B = -3.8576E02;
!!$    pc%SecondViralCoef%C = -7.2988E07;
!!$    pc%SecondViralCoef%D = -1.6901E20;
!!$    pc%SecondViralCoef%E = 2.6260E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 291.40, -5.4846E00; 1500.00, 2.7595E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1283E01;
!!$    pc%liqViscosity%B = 1.1519E03;
!!$    pc%liqViscosity%C = 1.7447E-03;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 192.78, 4.9985E-03; 416.43, 2.0218E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 5.3451E-07;
!!$    pc%vapViscosity%B = 5.2741E-01;
!!$    pc%vapViscosity%C = 2.7431E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 192.78, 3.5382E-06; 1000.00, 1.6029E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0102E-01;
!!$    pc%liqThermalConductivity%B = -2.4765E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 192.78, 1.5328E-01; 469.00, 8.4872E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 4.6511E-05;
!!$    pc%vapThermalConductivity%B = 1.0987E00;
!!$    pc%vapThermalConductivity%C = -1.1797E02;
!!$    pc%vapThermalConductivity%D = 1.6726E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 416.43, 2.0894E-02; 1000.0, 8.7652E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.3154E-02;
!!$    pc%SurfaceTension%B = 1.2669E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 192.78, 3.1956E-02; 582.80, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2methyloctane
!!$
!!$
!!$  subroutine pc_2nonanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2nonanone','2-nonanone','2-Nonanone'}; 
!!$    pc%CAS = '821-55-6';
!!$    ! components
!!$    pc%ChemicalFormula%C = 9;         
!!$    pc%ChemicalFormula%H = 18;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 142.241;  ! kg/kmol    
!!$    pc%Tcrit = 652.50;             ! K    
!!$    pc%pcrit = 2.4100E6;           ! Pa 
!!$    pc%vcrit = 0.56000;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.249; ! -
!!$
!!$
!!$    pc%MeltingPoint = 265.65;              ! K
!!$    pc%TripletPointT = 265.65;             ! K
!!$    pc%TripletPointp = 4.7226E00;          ! Pa
!!$    pc%NormalBoilingPoint = 467.15;        ! K
!!$    pc%liqMolarVolume = 0.173947;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.4070E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.0410E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  5.3560E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.5600E7;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -5.3212E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = - 1.2270E8;              ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.9680E8;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 4.0830E5;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4979;            ! -
!!$    pc%RadiusOfGyration = 5.3210E-10;      ! m
!!$    pc%SolubilityParameter = 1.7560E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.1400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.10042;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.3940E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4187;           ! -
!!$    pc%FlashPoint = 337.15;                ! K
!!$    pc%FlammabilityLimits = (/ 0.9 5.9 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 340 380 /); ! K
!!$    pc%AutoignitionTemp = 633.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 6.9300E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 265.65, 6.9300E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 4.4319E-01;
!!$    pc%liqDensity%B = 2.4818E-01;
!!$    pc%liqDensity%C = 6.5250E02;
!!$    pc%liqDensity%D = 2.8763E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 265.65, 5.9233E00; 652.50, 1.7858E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.9160E01;
!!$    pc%vapPressure%B = -7.8739E03;
!!$    pc%vapPressure%C = -5.0102E00;
!!$    pc%vapPressure%D =  9.1592E-19;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 265.65, 4.7226E00; 652.50, 2.4054E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 7.1729E07;
!!$    pc%HeatOfVap%B = 4.0270E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 265.65, 5.8112E07; 652.50, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 2.5700E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 265.65, 2.5700E05; MAX_REAL_WP,MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.5400E05;
!!$    pc%liqHeatCapacity%B = -5.5822E02;
!!$    pc%liqHeatCapacity%C = 1.3135E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 265.65, 2.9840E05; 510.00, 4.1095E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.5347E05;
!!$    pc%IG_HeatCapacity%B = 4.2140E05;
!!$    pc%IG_HeatCapacity%C = 1.3160E03;
!!$    pc%IG_HeatCapacity%D = 2.4060E05;
!!$    pc%IG_HeatCapacity%E = 6.4420E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.6947E05; 1500.00, 5.1866E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.8301E-01;
!!$    pc%SecondViralCoef%B = -4.4922E02;
!!$    pc%SecondViralCoef%C = -1.3676E08;
!!$    pc%SecondViralCoef%D = -8.4430E20;
!!$    pc%SecondViralCoef%E = 1.1643E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 326.25, -8.7298E00; 1500.00, 4.2979E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1678E01;
!!$    pc%liqViscosity%B = 1.5184E03;
!!$    pc%liqViscosity%C = -2.9833E-02;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 265.65, 2.1793E-03; 467.15, 1.8209E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 7.3861E-08;
!!$    pc%vapViscosity%B = 7.9750E-01;
!!$    pc%vapViscosity%C = 9.7893E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 265.65, 4.6297E-06; 1000.00, 1.6609E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.9050E-01;
!!$    pc%liqThermalConductivity%B = -1.9230E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 265.65, 1.3942E-01; 593.15, 7.6437E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -2.2439E-03;
!!$    pc%vapThermalConductivity%B = 1.0196E00;
!!$    pc%vapThermalConductivity%C = -2.6430E04;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 467.15, 2.1276E-02; 1000.0, 1.0103E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0312E-02;
!!$    pc%SurfaceTension%B = 1.0822E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 265.65, 2.8574E-02; 652.50, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2nonanone
!!$
!!$
!!$  subroutine pc_2octanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2octanone','2-octanone','methyl hexyl ketone'}; 
!!$    pc%CAS = '111-13-7';
!!$    ! components
!!$    pc%ChemicalFormula%C = 8;         
!!$    pc%ChemicalFormula%H = 16;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 128.214;  ! kg/kmol    
!!$    pc%Tcrit = 632.70;             ! K    
!!$    pc%pcrit = 2.6400E06;          ! Pa 
!!$    pc%vcrit = 0.49700;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.249; ! -
!!$
!!$
!!$    pc%MeltingPoint = 252.85;              ! K
!!$    pc%TripletPointT = 252.85;             ! K
!!$    pc%TripletPointp = 4.6773E00;          ! Pa
!!$    pc%NormalBoilingPoint = 446.15;        ! K
!!$    pc%liqMolarVolume = 0.157287;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.2160E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.1380E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.9620E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.4419E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.6984E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.2917E08;              ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.7270E08;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.7384E05;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4549;            ! -
!!$    pc%RadiusOfGyration = 4.8920E-10;      ! m
!!$    pc%SolubilityParameter = 1.7580E4;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 8.2100E-30;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.09019;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.2590E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4133;           ! -
!!$    pc%FlashPoint = 324.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.0 6.1 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 326 363 /); ! K
!!$    pc%AutoignitionTemp = 653.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 7.7800E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 252.85, 7.7800E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.0006E-01;
!!$    pc%liqDensity%B = 2.48511E-01;
!!$    pc%liqDensity%C = 6.3270E02;
!!$    pc%liqDensity%D = 2.9942E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 192.78, 6.1230E00; 582.80, 1.8484E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  6.3775E01;
!!$    pc%vapPressure%B = -7.7113E03;
!!$    pc%vapPressure%C = -5.7359E00;
!!$    pc%vapPressure%D =  3.0902E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 252.85, 4.6773E00; 632.70, 2.6466E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.5363E07;
!!$    pc%HeatOfVap%B = 3.8718E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 252.85, 5.3646E07; 632.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.9270E04;
!!$    pc%solHeatCapacity%B = 1.9170E03;
!!$    pc%solHeatCapacity%C = -8.2019E00;
!!$    pc%solHeatCapacity%D = 1.6844E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 13.00, 4.3019E03; 252.85, 2.1336E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.0040E05;
!!$    pc%liqHeatCapacity%B = -4.2620E02;
!!$    pc%liqHeatCapacity%C = 1.1172E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 252.86, 2.6406E05; 500.00, 3.6660E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.3901E05;
!!$    pc%IG_HeatCapacity%B = 3.8060E05;
!!$    pc%IG_HeatCapacity%C = 1.3717E03;
!!$    pc%IG_HeatCapacity%D = 2.2573E05;
!!$    pc%IG_HeatCapacity%E = 6.6096E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 1.4162E05; 1500.00, 4.6547E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.3845E-01;
!!$    pc%SecondViralCoef%B = -3.9878E02;
!!$    pc%SecondViralCoef%C = -1.0424E08;
!!$    pc%SecondViralCoef%D = -5.3291E20;
!!$    pc%SecondViralCoef%E = 7.0168E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 316.35, -7.3161E00; 1500.00, 4.1692E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1736E01;
!!$    pc%liqViscosity%B = 1.4152E03;
!!$    pc%liqViscosity%C = 3.6180E-04;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 252.85, 2.1614E-03; 446.15, 1.9128E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 8.0901E-08;
!!$    pc%vapViscosity%B = 7.9062E-01;
!!$    pc%vapViscosity%C = 9.9338E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 252.85, 4.6110E-06; 1000.00, 1.7325E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1320E-01;
!!$    pc%liqThermalConductivity%B = -2.4940E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 252.85, 1.5014E-01; 499.00, 8.8749E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -2.0184E-03;
!!$    pc%vapThermalConductivity%B = 1.0027E00;
!!$    pc%vapThermalConductivity%C = -2.0406E04;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 446.15, 2.0463E-02; 1000.0, 1.0597E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.9062E-02;
!!$    pc%SurfaceTension%B = 1.2967E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 252.85, 3.0488E-02; 632.70, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2octanone
!!$
!!$
!!$  subroutine pc_2pentanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2pentanol','2-pentanol','pentanol-2'}; 
!!$    pc%CAS = '6032-29-7';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 88.1497;  ! kg/kmol    
!!$    pc%Tcrit = 560.3;              ! K    
!!$    pc%pcrit = 3.6750E6;           ! Pa 
!!$    pc%vcrit = 0.329;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.26;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 200.0;              ! K
!!$    pc%TripletPointT = 200.0;             ! K
!!$    pc%TripletPointp = 5.9659E-3;         ! Pa
!!$    pc%NormalBoilingPoint = 392.15;       ! K
!!$    pc%liqMolarVolume = 0.10957;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.1380E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5930E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.9640E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  8.4800E06;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.0515E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.6850E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.6520E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.5500E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.558707;          ! -
!!$    pc%RadiusOfGyration = 3.6190E-10;      ! m
!!$    pc%SolubilityParameter = 2.1700E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.3900E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06262;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.9700E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4044;           ! -
!!$    pc%FlashPoint = 307.0;                 ! K
!!$    pc%FlammabilityLimits = (/ 1.5 9.7 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 307 338 /); ! K
!!$    pc%AutoignitionTemp = 616.0;           ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1700E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.0126E-01;
!!$    pc%liqDensity%B = 2.6371E-01;
!!$    pc%liqDensity%C = 5.6030E02;
!!$    pc%liqDensity%D = 2.5302E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 200.00, 1.0008E01; 560.30, 3.0384E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.1227E02;
!!$    pc%vapPressure%B = -1.0221E04;
!!$    pc%vapPressure%C = -1.2511E01;
!!$    pc%vapPressure%D =  8.5670E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 200.00, 5.9659E-03; 560.30, 3.6767E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 8.1306E07;
!!$    pc%HeatOfVap%B = 5.3685E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 200.00, 6.4147E07; 560.30, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 0;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP,MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.7400E05;
!!$    pc%liqHeatCapacity%B = -6.5134E01;
!!$    pc%liqHeatCapacity%C = 9.5351E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 280.15,2.3059E05; 420.22, 3.1500E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0385E05;
!!$    pc%IG_HeatCapacity%B = 2.9140E05;
!!$    pc%IG_HeatCapacity%C = 1.6383E03;
!!$    pc%IG_HeatCapacity%D = 1.9435E05;
!!$    pc%IG_HeatCapacity%E = 7.5094E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 298.15, 1.3604E05; 1500.00, 3.4084E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.3739E-01;
!!$    pc%SecondViralCoef%B = -2.4456E02;
!!$    pc%SecondViralCoef%C = -3.2955E07;
!!$    pc%SecondViralCoef%D = 3.8515E20;
!!$    pc%SecondViralCoef%E = -1.6083E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 280.15, -7.1139E00; 1500.00, 6.4596E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.0809E02;
!!$    pc%liqViscosity%B = 7.5182E03;
!!$    pc%liqViscosity%C = 1.3552E01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 250.00, 4.1157E-02; 392.15, 3.3865E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.1030E-07;
!!$    pc%vapViscosity%B = 6.4060E-01;
!!$    pc%vapViscosity%C = 2.3282E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 200.00, 4.2711E-06; 1000.00, 2.1022E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2037E-01;
!!$    pc%liqThermalConductivity%B = -2.8791E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 200.00, 1.6279E-01; 435.60, 9.4956E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.8508E03;
!!$    pc%vapThermalConductivity%B = 9.0265E-01;
!!$    pc%vapThermalConductivity%C = 1.2060E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 392.15, 2.0326E-02; 1000.0, 1.2066E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.5779E-02;
!!$    pc%SurfaceTension%B = 1.1409E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 200.00, 3.3705E-02; 560.30, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2pentanol
!!$
!!$
!!$  subroutine pc_2pentanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2pentanone','2-pentanone','Ethylacetone'}; 
!!$    pc%CAS = '107-87-9';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 86.1338;  ! kg/kmol    
!!$    pc%Tcrit = 561.08;             ! K    
!!$    pc%pcrit = 3.6940E06;          ! Pa 
!!$    pc%vcrit = 0.301;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.238; ! -
!!$
!!$
!!$    pc%MeltingPoint = 196.29;              ! K
!!$    pc%TripletPointT = 196.29;             ! K
!!$    pc%TripletPointp = 0.752347;           ! Pa
!!$    pc%NormalBoilingPoint = 375.46;        ! K
!!$    pc%liqMolarVolume = 0.10725;           ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.5920E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.3830E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.7860E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.0620E07;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.8796E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.4520E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -2.9729E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.7410E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.343288;          ! -
!!$    pc%RadiusOfGyration = 3.6180E-10;      ! m
!!$    pc%SolubilityParameter = 1.8290E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.2400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.0595;         ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.5400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.388;            ! -
!!$    pc%FlashPoint = 280.0;                 ! K
!!$    pc%FlammabilityLimits = (/ 1.5 8.2 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 278 309 /); ! K
!!$    pc%AutoignitionTemp = 725.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.2200E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 196.29, 1.2200E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.0411E-01;
!!$    pc%liqDensity%B = 2.7207E-01;
!!$    pc%liqDensity%C = 5.6108E02;
!!$    pc%liqDensity%D = 3.0669E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 196.29, 1.0398E01; 561.08, 3.3231E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  8.4635E01;
!!$    pc%vapPressure%B = -7.0784E03;
!!$    pc%vapPressure%C = -9.3000E00;
!!$    pc%vapPressure%D =  6.2702E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 196.29, 7.5235E-01; 561.08, 3.7062E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.1740E07;
!!$    pc%HeatOfVap%B = 3.9422E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 196.29, 4.3663E07; 561.08, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -2.5097E04;
!!$    pc%solHeatCapacity%B = 2.0373E03;
!!$    pc%solHeatCapacity%C = -1.4289E01;
!!$    pc%solHeatCapacity%D = 4.4874E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 30.00, 2.4373E04; 110.00, 8.5836E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.9459E05;
!!$    pc%liqHeatCapacity%B = -2.6386E02;
!!$    pc%liqHeatCapacity%C = 7.6808E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 196.29, 1.7239E05; 375.46, 2.0380E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 9.0053E04;
!!$    pc%IG_HeatCapacity%B = 2.7085E05;
!!$    pc%IG_HeatCapacity%C = 1.6592E03;
!!$    pc%IG_HeatCapacity%D = 1.8012E05;
!!$    pc%IG_HeatCapacity%E = 7.4396E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 9.5908E04; 1500.00, 3.0797E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.2325E-01;
!!$    pc%SecondViralCoef%B = -2.5789E02;
!!$    pc%SecondViralCoef%C = -3.5849E07;
!!$    pc%SecondViralCoef%D = -3.0892E20;
!!$    pc%SecondViralCoef%E = 5.9158E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 280.54, -4.8752E00; 1500.00, 4.0691E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1055E01;
!!$    pc%liqViscosity%B = 1.0053E03;
!!$    pc%liqViscosity%C = 3.9301E-03;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 250.00, 9.0091E-04; 375.46, 2.3541E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.4630E-07;
!!$    pc%vapViscosity%B = 6.6530E-01;
!!$    pc%vapViscosity%C = 2.0870E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 196.29, 4.0030E-06; 1000.00, 2.0186E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1610E-01;
!!$    pc%liqThermalConductivity%B = -2.4866E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 196.29, 1.6729E-01; 375.46, 1.2274E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -1.7190E-02;
!!$    pc%vapThermalConductivity%B = 4.8320E-01;
!!$    pc%vapThermalConductivity%C = -3.7980E03;
!!$    pc%vapThermalConductivity%D = -1.2350E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.00, 8.7672E-03; 1000.0, 1.2002E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.3802E-02;
!!$    pc%SurfaceTension%B = 1.0733E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 196.29, 3.3893E-02; 561.08, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_2pentanone
!!$
!!$
!!$  subroutine pc_2propanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'2propanol','2-propanol','dimethyl carbinol','isohol'}; 
!!$    pc%CAS = '67-63-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 3;         
!!$    pc%ChemicalFormula%H = 8;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 60.096;   ! kg/kmol    
!!$    pc%Tcrit = 508.31;             ! K    
!!$    pc%pcrit = 4.7643E6;           ! Pa 
!!$    pc%vcrit = 0.22013;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.248; ! -
!!$
!!$
!!$    pc%MeltingPoint = 185.28;              ! K
!!$    pc%TripletPointT = 185.28;             ! K
!!$    pc%TripletPointp = 3.2000E-2;          ! Pa
!!$    pc%NormalBoilingPoint = 355.41;        ! K
!!$    pc%liqMolarVolume = 0.076784;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.7242E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.7339E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.0991E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  5.4099E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.8300E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.6689;            ! -
!!$    pc%RadiusOfGyration = 2.8070E-10;      ! m
!!$    pc%SolubilityParameter = 2.3575E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.5372E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.04216;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 6.2700E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3752;           ! -
!!$    pc%FlashPoint = 284.82;                ! K
!!$    pc%FlammabilityLimits = (/ 2 12.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 672.04;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.71E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 185.28, 1.7100E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.18E00;
!!$    pc%liqDensity%B = 2.6475E-01;
!!$    pc%liqDensity%C = 5.0831E02;
!!$    pc%liqDensity%D = 2.43E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 185.28, 1.4656E01; 508.31, 4.4570E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  9.2935E01;
!!$    pc%vapPressure%B = -8.1771E03;
!!$    pc%vapPressure%C = -1.0031E01;
!!$    pc%vapPressure%D =  3.9988E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 185.28, 3.2000E-02; 508.31, 4.7673E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.6980E07;
!!$    pc%HeatOfVap%B = 8.7E-02;
!!$    pc%HeatOfVap%C = 3.0070E-01;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 185.28, 5.2121E07; 508.31, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -8.5E03;
!!$    pc%solHeatCapacity%B = 8.7600E02;
!!$    pc%solHeatCapacity%C = -8.100E-01;
!!$    pc%solHeatCapacity%D = -3.025E-02;
!!$    pc%solHeatCapacity%E = 1.4130E-04;
!!$    pc%solHeatCapacity%range(:,1) = (/ 12.00, 1.8460E03; 176.66, 9.1821E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 4.6640E05;
!!$    pc%liqHeatCapacity%B = -4.1086E03;
!!$    pc%liqHeatCapacity%C = 1.4506E01;
!!$    pc%liqHeatCapacity%D = -1.4126E-02;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 188.32, 1.1277E05; 463.15, 2.7174E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 4.7460E04;
!!$    pc%IG_HeatCapacity%B = 1.9350E05;
!!$    pc%IG_HeatCapacity%C = 1.1240E03;
!!$    pc%IG_HeatCapacity%D = 9.3800E04;
!!$    pc%IG_HeatCapacity%E = 4.6000E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 100.00, 4.8262E04; 1500.00, 2.1652E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.0200E-02;
!!$    pc%SecondViralCoef%B = -6.2900E00;
!!$    pc%SecondViralCoef%C = -4.3100E07;
!!$    pc%SecondViralCoef%D = 1.9860E20;
!!$    pc%SecondViralCoef%E = -8.6100E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 279.58, -4.8942E00; 2541.60, 2.5100E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -8.2300E00;
!!$    pc%liqViscosity%B = 2.2822E03;
!!$    pc%liqViscosity%C = -9.8495E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 187.35, 3.0037E-01; 354.45, 5.1389E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.9930E-07;
!!$    pc%vapViscosity%B = 7.2330E-01;
!!$    pc%vapViscosity%C = 1.7800E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 200.00, 4.8683E-06; 1000.00, 2.5019E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0290E-01;
!!$    pc%liqThermalConductivity%B = -2.2780E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 185.28, 1.6069E-01; 410.00, 1.0950E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -8.0642E01;
!!$    pc%vapThermalConductivity%B = -1.4549E00;
!!$    pc%vapThermalConductivity%C = -6.0442E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 355.41, 2.2388E-02; 450.00, 3.2427E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 3.818E-02;
!!$    pc%SurfaceTension%B = -3.8180E-05;
!!$    pc%SurfaceTension%C = -6.5100E-08;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 273.15, 2.2894E-02; 355.41, 1.6387E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_2propanol
!!$
!!$
!!$  subroutine pc_3heptanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'3heptanone','3-heptanone','butyl ethyl ketone','heptane-3-one'}; 
!!$    pc%CAS = '106-35-4';
!!$    ! components
!!$    pc%ChemicalFormula%C = 7;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 114.188;  ! kg/kmol    
!!$    pc%Tcrit = 606.60;             ! K    
!!$    pc%pcrit = 2.9200E6;           ! Pa 
!!$    pc%vcrit = 0.43300;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.251; ! -
!!$
!!$
!!$    pc%MeltingPoint = 234.15;              ! K
!!$    pc%TripletPointT = 234.15;             ! K
!!$    pc%TripletPointp = 2.2992E00;          ! Pa
!!$    pc%NormalBoilingPoint = 420.55;        ! K
!!$    pc%liqMolarVolume = 0.140257;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.0100E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.2250E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.5800E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.8100E7;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.0980E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.3580E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.4860E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.4140E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4076;            ! -
!!$    pc%RadiusOfGyration = 4.5390E-10;      ! m
!!$    pc%SolubilityParameter = 1.7940E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.3700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.07996;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.1240E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4066;           ! -
!!$    pc%FlashPoint = 310.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.1 6.8 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 309 344 /); ! K
!!$    pc%AutoignitionTemp = 683.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 8.8600E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 234.15, 8.8600E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.9268E-01;
!!$    pc%liqDensity%B = 2.5663E-01;
!!$    pc%liqDensity%C = 6.0660E02;
!!$    pc%liqDensity%D = 2.7766E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 234.15, 7.5751E00; 606.60, 2.3095E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.8463E01;
!!$    pc%vapPressure%B = -8.0772E03;
!!$    pc%vapPressure%C = -7.9062E00;
!!$    pc%vapPressure%D =  8.0521E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 234.15, 2.2992E-04; 606.60, 2.9190E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.3357E07;
!!$    pc%HeatOfVap%B = 4.2167E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 234.15, 5.1579E07; 606.60, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 2.0400E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 234.15, 2.0400E05; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.7073E05;
!!$    pc%liqHeatCapacity%B = -3.99893E02;
!!$    pc%liqHeatCapacity%C = 1.0601E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 234.15, 2.3522E05; 480.00, 3.2303E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2768E05;
!!$    pc%IG_HeatCapacity%B = 3.3810E05;
!!$    pc%IG_HeatCapacity%C = 1.3831E03;
!!$    pc%IG_HeatCapacity%D = 1.8880E05;
!!$    pc%IG_HeatCapacity%E = 6.5030E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.3968E05; 1500.00, 4.1386E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.9056E-01;
!!$    pc%SecondViralCoef%B = -3.3949E02;
!!$    pc%SecondViralCoef%C = -7.6495E07;
!!$    pc%SecondViralCoef%D = -4.1164E20;
!!$    pc%SecondViralCoef%E = 6.0973E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 303.3, -6.5114E00; 1500.00, 4.1554E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -9.3877E00;
!!$    pc%liqViscosity%B = 1.2049E03;
!!$    pc%liqViscosity%C = -3.2618E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 234.15, 2.4268E-03; 421.15, 2.0397E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 8.9656E-08;
!!$    pc%vapViscosity%B = 7.8236E-01;
!!$    pc%vapViscosity%C = 1.0014E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 234.15, 4.4848E-06; 1000.00, 1.8122E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0260E-01;
!!$    pc%liqThermalConductivity%B = -2.2340E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 234.15, 1.5029E-01; 553.15, 7.9026E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 1.3486E03;
!!$    pc%vapThermalConductivity%B = 1.0313E00;
!!$    pc%vapThermalConductivity%C = 1.4832E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 420.55, 1.9429E-02; 1000.0, 1.1287E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.8308E-02;
!!$    pc%SurfaceTension%B = 1.2123E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 234.15, 3.2279E-02; 606.60, 8.3816E-04 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_3heptanone
!!$
!!$
!!$  subroutine pc_3hexanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'3hexanone','3-hexanone','ethyl-n-propylketone'}; 
!!$    pc%CAS = '589-38-8';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 100.161;  ! kg/kmol    
!!$    pc%Tcrit = 582.82;             ! K    
!!$    pc%pcrit = 3.3200E6;           ! Pa 
!!$    pc%vcrit = 0.37800;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.259; ! -
!!$
!!$
!!$    pc%MeltingPoint = 217.15;              ! K
!!$    pc%TripletPointT = 217.15;             ! K
!!$    pc%TripletPointp = 2.2202E00;          ! Pa
!!$    pc%NormalBoilingPoint = 396.65;        ! K
!!$    pc%liqMolarVolume = 0.123630;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.7760E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.2600E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.0920E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.3490E7;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.4920E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.3670E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.2020E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.0531E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3801;            ! -
!!$    pc%RadiusOfGyration = 4.0090E-10;      ! m
!!$    pc%SolubilityParameter = 1.7930E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.5700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06973;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 9.8900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3980;           ! -
!!$    pc%FlashPoint = 287.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.0 8.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 288 327 /); ! K
!!$    pc%AutoignitionTemp = 712.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.0200E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 217.50, 1.0200E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.7666E-01;
!!$    pc%liqDensity%B = 2.5578E-01;
!!$    pc%liqDensity%C = 5.8282E02;
!!$    pc%liqDensity%D = 2.7746E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 217.50, 8.7631E00; 582.82, 2.6455E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.3155E01;
!!$    pc%vapPressure%B = -7.2429E03;
!!$    pc%vapPressure%C = -7.2569E00;
!!$    pc%vapPressure%D =  1.2741E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 217.50, 2.2202E00; 582.82, 3.3216E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.6232E07;
!!$    pc%HeatOfVap%B = 3.9972E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 217.50, 4.6655E07; 582.82, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -4.1000E03;
!!$    pc%solHeatCapacity%B = 4.2800E02;
!!$    pc%solHeatCapacity%C = 2.6860E01;
!!$    pc%solHeatCapacity%D = -3.5340E-01;
!!$    pc%solHeatCapacity%E = 1.3463E-03;
!!$    pc%solHeatCapacity%range(:,1) = (/ 12.52, 4.8084E03; 140.13, 1.2999E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.3596E05;
!!$    pc%liqHeatCapacity%B = -3.4594E02;
!!$    pc%liqHeatCapacity%C = 9.4278E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 217.50, 2.0532E05; 460.00, 2.7632E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.1237E05;
!!$    pc%IG_HeatCapacity%B = 2.9360E05;
!!$    pc%IG_HeatCapacity%C = 1.4010E03;
!!$    pc%IG_HeatCapacity%D = 1.6010E05;
!!$    pc%IG_HeatCapacity%E = 6.5050E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 1.1443E05; 1500.00, 3.5874E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.4188E-01;
!!$    pc%SecondViralCoef%B = -2.7496E02;
!!$    pc%SecondViralCoef%C = -5.8174E07;
!!$    pc%SecondViralCoef%D = -2.1114E20;
!!$    pc%SecondViralCoef%E = 2.4535E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 297.00, -5.0275E00; 1500.00, 4.1329E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.3682E01;
!!$    pc%liqViscosity%B = 1.2834E03;
!!$    pc%liqViscosity%C = 3.3755E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 217.50, 2.5631E-03; 396.65, 2.1848E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 9.8882E-08;
!!$    pc%vapViscosity%B = 7.7550E-01;
!!$    pc%vapViscosity%C = 9.9825E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 217.50, 4.4032E-06; 1000.00, 1.9067E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.3493E-01;
!!$    pc%liqThermalConductivity%B = -2.9120E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 217.50, 1.7159E-01; 466.00, 9.9231E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -3.3262E-01;
!!$    pc%vapThermalConductivity%B = 1.2054E-01;
!!$    pc%vapThermalConductivity%C = -2.4726E03;
!!$    pc%vapThermalConductivity%D = -5.4934E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.00, 7.9990E-03; 1000.0, 1.0980E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.7740E-02;
!!$    pc%SurfaceTension%B = 1.1689E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 217.50, 3.3447E-02; 582.82, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_3hexanone
!!$
!!$
!!$  subroutine pc_3pentanone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'3pentanone','3-pentanone','metacetone','ethyl ketone'}; 
!!$    pc%CAS = '96-22-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 86.134;   ! kg/kmol    
!!$    pc%Tcrit = 560.95;             ! K    
!!$    pc%pcrit = 3.7400E6;           ! Pa 
!!$    pc%vcrit = 0.33600;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.269; ! -
!!$
!!$
!!$    pc%MeltingPoint = 234.18;              ! K
!!$    pc%TripletPointT = 234.18;             ! K
!!$    pc%TripletPointp = 7.3422E01;          ! Pa
!!$    pc%NormalBoilingPoint = 375.14;        ! K
!!$    pc%liqMolarVolume = 0.106410;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.5790E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.3440E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.7000E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.1590E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.8804E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.4200E08;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -2.9651E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.6600E05;         ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3448;            ! -
!!$    pc%RadiusOfGyration = 3.5830E-10;      ! m
!!$    pc%SolubilityParameter = 1.8410E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.4100E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05950;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.5400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3900;           ! -
!!$    pc%FlashPoint = 286.00;                ! K
!!$    pc%FlammabilityLimits = (/ 1.5 8.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 277 308 /); ! K
!!$    pc%AutoignitionTemp = 725.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1800E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 234.18, 1.1800E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 7.1811E-01;
!!$    pc%liqDensity%B = 2.4129E-01;
!!$    pc%liqDensity%C = 5.6095E02;
!!$    pc%liqDensity%D = 2.7996E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 234.18, 1.0102E01; 560.95, 2.9761E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A = 4.4286E01;
!!$    pc%vapPressure%B = -5.4151E03;
!!$    pc%vapPressure%C = -3.0913E00;
!!$    pc%vapPressure%D =  1.8580E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 234.18, 7.3422E01; 560.95, 3.6993E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.2359E07;
!!$    pc%HeatOfVap%B = 4.0465E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 234.18, 4.2075E07; 560.95, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -9.2023E03;
!!$    pc%solHeatCapacity%B = 1.1108E03;
!!$    pc%solHeatCapacity%C = 5.6357E00;
!!$    pc%solHeatCapacity%D = -1.4318E-01;
!!$    pc%solHeatCapacity%E = 6.8742E-04;
!!$    pc%solHeatCapacity%range(:,1) = (/ 11.16, 3.7078E03; 110.39, 9.1569E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.9302E05;
!!$    pc%liqHeatCapacity%B = -1.7643E02;
!!$    pc%liqHeatCapacity%C = 5.6690E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 234.18, 1.8279E05; 375.14, 2.0661E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 9.6896E04;
!!$    pc%IG_HeatCapacity%B = 2.4907E05;
!!$    pc%IG_HeatCapacity%C = 1.4177E03;
!!$    pc%IG_HeatCapacity%D = 1.3010E05;
!!$    pc%IG_HeatCapacity%E = 6.4670E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.0536E05; 1500.00, 3.0358E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.9680E-01;
!!$    pc%SecondViralCoef%B = -2.1873E02;
!!$    pc%SecondViralCoef%C = -4.5060E07;
!!$    pc%SecondViralCoef%D = -1.5690E20;
!!$    pc%SecondViralCoef%E = 1.9653E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 280.47,-4.8930E00; 1500.00, 3.7623E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.8695E00;
!!$    pc%liqViscosity%B = 5.9632E02;
!!$    pc%liqViscosity%C = -1.2025E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 234.18,1.0240E-03; 375.14, 2.2319E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.1640E-07;
!!$    pc%vapViscosity%B = 7.6150E-01;
!!$    pc%vapViscosity%C = 1.0794E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 234.18, 5.0785E-06; 1000.00, 2.0227E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1569E-01;
!!$    pc%liqThermalConductivity%B = -2.4081E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 234.18, 1.5930E-01; 375.14, 1.2535E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.2775E01;
!!$    pc%vapThermalConductivity%B = 1.0019E00;
!!$    pc%vapThermalConductivity%C = 1.9100E08;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.00, 8.9821E-03; 1000.00, 1.2082E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.0021E-02;
!!$    pc%SurfaceTension%B = 1.1704E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 234.18, 3.1888E-02; 560.95, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_3pentanone
!!$
!!$
!!$  subroutine pc_3petanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'3petanol','3-petanol','pentanol-3','diethyl carbinol'}; 
!!$    pc%CAS = '584-02-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 88.150;  ! kg/kmol    
!!$    pc%Tcrit = 559.60;             ! K    
!!$    pc%pcrit = 3.8800E6;           ! Pa 
!!$    pc%vcrit = 0.32500;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.271; ! -
!!$
!!$
!!$    pc%MeltingPoint = 204.15;              ! K
!!$    pc%TripletPointT = 204.15;             ! K
!!$    pc%TripletPointp = 1.4819E-02;         ! Pa
!!$    pc%NormalBoilingPoint = 388.45;        ! K
!!$    pc%liqMolarVolume = 0.108015;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.1673E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5816E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.8200E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  9.0800E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.0483E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.7050E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.6890E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.4910E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5383;            ! -
!!$    pc%RadiusOfGyration = 3.2530E-10;      ! m
!!$    pc%SolubilityParameter = 2.1780E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.4700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06262;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.9700E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4079;           ! -
!!$    pc%FlashPoint = 313.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.2 9.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 300 332 /); ! K
!!$    pc%AutoignitionTemp = 708.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1900E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 204.15, 1.1900E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.0198E-01;
!!$    pc%liqDensity%B = 2.6069E-01;
!!$    pc%liqDensity%C = 5.5960E02;
!!$    pc%liqDensity%D = 2.6161E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 204.15, 1.0152E01; 559.60, 3.0764E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.3775E02;
!!$    pc%vapPressure%B = -1.1249E04;
!!$    pc%vapPressure%C = -1.6331E01;
!!$    pc%vapPressure%D =  2.8016E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 204.15, 1.4819E-02; 559.60, 3.9083E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 8.4954E07;
!!$    pc%HeatOfVap%B = 6.0261E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 204.15, 6.4626E07; 559.60, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 0;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP,MAX_REAL_WP; MAX_REAL_WP,MAX_REAL_WP  /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = -1.6787E06;
!!$    pc%liqHeatCapacity%B = 1.1243E04;
!!$    pc%liqHeatCapacity%C = -1.6088E01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 303.02, 2.5094E05; 326.16, 2.7687E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0751E05;
!!$    pc%IG_HeatCapacity%B = 3.0130E05;
!!$    pc%IG_HeatCapacity%C = 1.8100E03;
!!$    pc%IG_HeatCapacity%D = 2.1755E05;
!!$    pc%IG_HeatCapacity%E = 8.1170E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 300, 1.3596E05; 1200.00, 3.2093E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.2585E-01;
!!$    pc%SecondViralCoef%B = -2.3629E02;
!!$    pc%SecondViralCoef%C = -3.0124E07;
!!$    pc%SecondViralCoef%D = 3.5131E20;
!!$    pc%SecondViralCoef%E = -1.5281E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 279.80, -7.1804E00; 1500.00, 5.9407E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.5290E02;
!!$    pc%liqViscosity%B = 1.4984E04;
!!$    pc%liqViscosity%C = 3.4620E01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 250.00, 1.6349E-01; 388.45, 3.6484E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$
!!$    pc%vapViscosity%A = 1.1891E-07;
!!$    pc%vapViscosity%B = 7.6401E-01;
!!$    pc%vapViscosity%C = 1.0337E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1707E-01;
!!$    pc%liqThermalConductivity%B = -2.7901E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 3.1596E03;
!!$    pc%vapThermalConductivity%B = 9.0260E-01;
!!$    pc%vapThermalConductivity%C = 1.2930E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.8608E-02;
!!$    pc%SurfaceTension%B = 1.1579E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_3petanol
!!$
!!$
!!$  subroutine pc_acetone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'acetone','Acetone','Allylic Alcohol','Dimethylketal','3-Hydroxypropene','Propanone','beta-Ketopropane'}; 
!!$    pc%CAS = '67-64-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 3;         
!!$    pc%ChemicalFormula%H = 6;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 58.080;   ! kg/kmol    
!!$    pc%Tcrit = 508.20;             ! K    
!!$    pc%pcrit = 4.7015E6;          ! Pa 
!!$    pc%vcrit = 0.20900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.233; ! -
!!$
!!$
!!$    pc%MeltingPoint = 178.45;              ! K
!!$    pc%TripletPointT = 178.45;             ! K
!!$    pc%TripletPointp = 2.5938E00;          ! Pa
!!$    pc%NormalBoilingPoint = 329.44;        ! K
!!$    pc%liqMolarVolume = 0.073931;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.1715E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.5272E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.9535E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  5.6912E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.6592E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3064;            ! -
!!$    pc%RadiusOfGyration = 2.7460E-10;      ! m
!!$    pc%SolubilityParameter = 1.9774E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 9.6066E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.03904;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 5.8400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3560;           ! -
!!$    pc%FlashPoint = 255.37;                ! K
!!$    pc%FlammabilityLimits = (/ 2.6 12.8 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 810.93;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.6677E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 174.15, 1.6677E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.2298E00;
!!$    pc%liqDensity%B = 2.5760E-01;
!!$    pc%liqDensity%C = 5.0820E02;
!!$    pc%liqDensity%D = 2.9903E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 178.45, 1.5721E01; 508.20, 4.7741E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.0720E01;
!!$    pc%vapPressure%B = -5.6850E03;
!!$    pc%vapPressure%C = -7.3510E00;
!!$    pc%vapPressure%D =  6.300E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 178.45, 2.5938E00; 508.20, 4.6733E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.9170E07;
!!$    pc%HeatOfVap%B = 1.0360E00;
!!$    pc%HeatOfVap%C = -1.2940E00;
!!$    pc%HeatOfVap%D = 6.7200E-01;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 178.45, 3.6887E07; 508.2, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.2000E04;
!!$    pc%solHeatCapacity%B = 1.2760E03;
!!$    pc%solHeatCapacity%C = -6.1700E00;
!!$    pc%solHeatCapacity%D = 1.3200E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 30.00, 2.1083E04; 178.45, 9.4233E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.3560E05;
!!$    pc%liqHeatCapacity%B = -1.7700E02;
!!$    pc%liqHeatCapacity%C = 2.8370E-01;
!!$    pc%liqHeatCapacity%D = 6.8900E-04;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 178.45,1.169E05; 329.44, 1.3271E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 4.8120E04;
!!$    pc%IG_HeatCapacity%B = 1.6440E05;
!!$    pc%IG_HeatCapacity%C = 1.2500E03;
!!$    pc%IG_HeatCapacity%D = 7.1700E04;
!!$    pc%IG_HeatCapacity%E = -5.244E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 100.00, 4.8340E04; 1500.00, 1.8699E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.1070E-01;
!!$    pc%SecondViralCoef%B = -9.9200E01;
!!$    pc%SecondViralCoef%C = -3.0540E07;
!!$    pc%SecondViralCoef%D = 9.8900E18;
!!$    pc%SecondViralCoef%E = -1.6454E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 254.10,-5.298E00; 2540.50, 6.9790E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.4918E01;
!!$    pc%liqViscosity%B = 1.0234E03;
!!$    pc%liqViscosity%C = 5.9610E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 223.15, 8.1834E-04; 329.44, 2.3506E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.1005E-08;
!!$    pc%vapViscosity%B = 9.7620E-01;
!!$    pc%vapViscosity%C = 2.3139E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 300.00, 7.5393E-06; 650.00, 1.6680E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.5020E-01;
!!$    pc%liqThermalConductivity%B = -2.9800E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 270.00, 1.6974E-01; 315.00, 1.5633E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -2.6800E01;
!!$    pc%vapThermalConductivity%B = 9.0980E-01;
!!$    pc%vapThermalConductivity%C = -1.2650E08;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.15, 9.5298E-03; 571.70, 3.9056E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.2200E-02;
!!$    pc%SurfaceTension%B = 1.1240E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 178.45, 3.8251E-02; 508.20, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_acetone
!!$
!!$
!!$  subroutine pc_diethylether(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'diethylether','aether','ether','ethyl ether'}; 
!!$    pc%CAS = '60-29-7';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 74.123;   ! kg/kmol    
!!$    pc%Tcrit = 466.70;             ! K    
!!$    pc%pcrit = 3.6400E6;           ! Pa 
!!$    pc%vcrit = 0.28000;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.263; ! -
!!$
!!$
!!$    pc%MeltingPoint = 156.85;              ! K
!!$    pc%TripletPointT = 156.85;             ! K
!!$    pc%TripletPointp = 3.9545E-01;         ! Pa
!!$    pc%NormalBoilingPoint = 307.58;        ! K
!!$    pc%liqMolarVolume = 0.104692;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.5210E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.2210E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.4230E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  7.1900E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.2032E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2811;            ! -
!!$    pc%RadiusOfGyration = 3.1770E-10;      ! m
!!$    pc%SolubilityParameter = 1.5420E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 3.8400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05150;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.5400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3495;           ! -
!!$    pc%FlashPoint = 228.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.9 48.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp =433.15;           ! K
!!$
!!$
!!$    pc%solDensity%A = 1.3400E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 156.85, 1.3400E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.5540E-01;
!!$    pc%liqDensity%B = 2.6847E-01;
!!$    pc%liqDensity%C = 4.6670E02;
!!$    pc%liqDensity%D = 2.8140E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 156.85, 1.1487E01; 466.70, 3.5587E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A = 1.3690E02;
!!$    pc%vapPressure%B = -6.9543E03;
!!$    pc%vapPressure%C = -1.9254E01;
!!$    pc%vapPressure%D =  2.4508E-02;
!!$    pc%vapPressure%E =  1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 156.85, 3.9545E-01; 466.70, 3.6412E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.0600E07;
!!$    pc%HeatOfVap%B = 3.8680E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 156.85, 3.4651E07; 466.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -2.3760E04;
!!$    pc%solHeatCapacity%B = 1.8611E03;
!!$    pc%solHeatCapacity%C = -1.2638E01;
!!$    pc%solHeatCapacity%D = 3.7650E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 30.00, 2.1715E04; 156.85 , 1.0252E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 4.4400E04;
!!$    pc%liqHeatCapacity%B = 1.3010E03;
!!$    pc%liqHeatCapacity%C = -5.5000E00;
!!$    pc%liqHeatCapacity%D = 8.7630E-03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 156.92, 1.4698E05; 460.00, 3.3202E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 8.6210E04;
!!$    pc%IG_HeatCapacity%B = 2.5510E05;
!!$    pc%IG_HeatCapacity%C = 1.5413E03;
!!$    pc%IG_HeatCapacity%D = 1.4370E05;
!!$    pc%IG_HeatCapacity%E = -6.8890E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 9.3157E04; 1500.00, 2.9244E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.6330E-01;
!!$    pc%SecondViralCoef%B = -1.3680E02;
!!$    pc%SecondViralCoef%C = -2.0080E07;
!!$    pc%SecondViralCoef%D = -1.7500E19;
!!$    pc%SecondViralCoef%E = 2.6570E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 233.35,-2.6987E00; 1500.00, 6.6150E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = 1.0197E01;
!!$    pc%liqViscosity%B = -6.3800E01;
!!$    pc%liqViscosity%C = -3.2260E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 200.00,7.3593E-00; 373.15, 1.1412E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.9480E-06;
!!$    pc%vapViscosity%B = 4.1000E-01;
!!$    pc%vapViscosity%C = 4.9580E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 156.85, 3.7200E-06; 1000.00, 2.2116E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.4950E-01;
!!$    pc%liqThermalConductivity%B = -4.0700E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 156.85, 1.8566E-01; 433.15, 7.3208E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -4.4894E-03;
!!$    pc%vapThermalConductivity%B = 6.1550E-01;
!!$    pc%vapThermalConductivity%C = -3.2663E03;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 200.00, 7.6364E-03; 600.00, 5.1806E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.7928E-02;
!!$    pc%SurfaceTension%B = 1.2376E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 156.85, 3.4893E-02; 466.70, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_diethylether
!!$
!!$
!!$  subroutine pc_dimethoxymethane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'OME1','Methylal','ome1','methylal'}; 
!!$    pc%CAS = '	109-87-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 3;         
!!$    pc%ChemicalFormula%H = 8;
!!$    pc%ChemicalFormula%O = 2;
!!$
!!$    pc%MolecularWeight = 76.095;   ! kg/kmol    
!!$    pc%Tcrit = 480.6;             ! K    
!!$    pc%pcrit = 3.9517E06;          ! Pa 
!!$    pc%vcrit = 0.21300;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.211; ! -
!!$
!!$
!!$    pc%MeltingPoint = 168.35;              ! K
!!$    pc%TripletPointT = 168.35;             ! K
!!$    pc%TripletPointp = 1.9634E00;         ! Pa
!!$    pc%NormalBoilingPoint = 315.00;        ! K
!!$    pc%liqMolarVolume = 0.089127;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.4820E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -2.2631E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.3572E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  8.3303E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.7998E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2904;            ! -
!!$    pc%RadiusOfGyration = 2.8910E-10;      ! m
!!$    pc%SolubilityParameter = 1.7189E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 2.4684E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.04497;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 6.7900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3504;           ! -
!!$    pc%FlashPoint = 255.37;                ! K
!!$    pc%FlammabilityLimits = (/ 1.6 17.6 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 510.35;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.5600E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 168.35, 1.5600E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.4355E00;
!!$    pc%liqDensity%B = 3.0576E-01;
!!$    pc%liqDensity%C = 4.8060E02;
!!$    pc%liqDensity%D = 3.17755E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 168.35, 1.3194E01; 480.6, 4.6949E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.6361E01;
!!$    pc%vapPressure%B = -5.0483E03;
!!$    pc%vapPressure%C = -5.0136E00;
!!$    pc%vapPressure%D =  2.3659E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 168.35, 1.9634E00; 480.6, 3.9500E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.1300E07;
!!$    pc%HeatOfVap%B = 1.7170E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 168.35, 3.5183E07; 480.6, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.6050E04;
!!$    pc%solHeatCapacity%B = 1.3984E03;
!!$    pc%solHeatCapacity%C = -5.8400E00;
!!$    pc%solHeatCapacity%D = 8.5900E-03;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 15.00, 3.6410E03; 168.35, 9.4841E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.6050E05;
!!$    pc%liqHeatCapacity%B = -2.2930E02;
!!$    pc%liqHeatCapacity%C = 7.86600E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 168.35, 1.4417E05; 310.00, 1.6495E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 7.4976E04;
!!$    pc%IG_HeatCapacity%B = 1.6166E05;
!!$    pc%IG_HeatCapacity%C = 8.6287E02;
!!$    pc%IG_HeatCapacity%D = 7.8964E05;
!!$    pc%IG_HeatCapacity%E = 4.6718E03;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 298.15, 9.1671E04; 1000.00, 2.0789E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.6669E-01;
!!$    pc%SecondViralCoef%B = -1.6805E02;
!!$    pc%SecondViralCoef%C = -1.9195E07;
!!$    pc%SecondViralCoef%D = 2.9849E19;
!!$    pc%SecondViralCoef%E = -8.7864E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 240.30, -2.5200; 2400.00, 9.5290E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -6.0611E00;
!!$    pc%liqViscosity%B = 6.3620E02;
!!$    pc%liqViscosity%C = -7.1926E-1;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 168.35, 2.5567E-03; 378.50, 1.7514E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 4.3034E-07;
!!$    pc%vapViscosity%B = 6.1720E-01;
!!$    pc%vapViscosity%C = 2.2554E+02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 315, 8.7348E-06; 975.00, 2.4449E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.5239E-01;
!!$    pc%liqThermalConductivity%B = -3.6877E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 168.35, 1.9031E-01; 315, 1.3623E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 1.1366E-03;
!!$    pc%vapThermalConductivity%B = 6.9110E-01;
!!$    pc%vapThermalConductivity%C = 4.1879E+02;
!!$    pc%vapThermalConductivity%D = 1.8761E+05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 315.00, 1.4350E-02; 1000.00, 8.3763E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 7.0967E-02;
!!$    pc%SurfaceTension%B = 1.2676;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 168.35, 4.1083E-02; 480.60, 0.00E00 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_dimethoxymethane
!!$
!!$
!!$  subroutine pc_dimethylether(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'dimethylether','methyl ether','ether','ethyl ether'}; 
!!$    pc%CAS = '115-10-6';
!!$    ! components
!!$    pc%ChemicalFormula%C = 2;         
!!$    pc%ChemicalFormula%H = 6;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 46.069;   ! kg/kmol    
!!$    pc%Tcrit = 400.10;             ! K    
!!$    pc%pcrit = 5.3700E6;           ! Pa 
!!$    pc%vcrit = 0.17000;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.274; ! -
!!$
!!$
!!$    pc%MeltingPoint = 131.66;              ! K
!!$    pc%TripletPointT = 131.65;             ! K
!!$    pc%TripletPointp = 3.0496E00;          ! Pa
!!$    pc%NormalBoilingPoint = 248.31;        ! K
!!$    pc%liqMolarVolume = 0.070231;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.8410E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.1280E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.6670E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  4.9370E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.3284E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2002;            ! -
!!$    pc%RadiusOfGyration = 2.1540E-10;      ! m
!!$    pc%SolubilityParameter = 1.5120E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 4.3400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.03104;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 4.8400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.2984;           ! -
!!$    pc%FlashPoint = 232.00;                ! K
!!$    pc%FlammabilityLimits = (/ 3.3 27.3 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp =623.15;           ! K
!!$
!!$
!!$    pc%solDensity%A = 2.2200E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 131.65, 2.2200E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.5693E00;
!!$    pc%liqDensity%B = 2.6790E-01;
!!$    pc%liqDensity%C = 4.0010E02;
!!$    pc%liqDensity%D = 2.8820E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 131.65, 1.8950E01; 400.10, 5.8578E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A = 4.4704E01;
!!$    pc%vapPressure%B = -3.5256E03;
!!$    pc%vapPressure%C = -3.4444E00;
!!$    pc%vapPressure%D =  5.4574E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 131.65, 3.0496E00; 400.10, 5.2735E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 2.9940E07;
!!$    pc%HeatOfVap%B = 3.5050E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 131.65, 2.6032E07; 400.10, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.3880E04;
!!$    pc%solHeatCapacity%B = 1.2580E03;
!!$    pc%solHeatCapacity%C = -8.2200E00;
!!$    pc%solHeatCapacity%D = 2.6752E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 24.26, 1.2183E04; 131.64 , 7.0305E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.1010E05;
!!$    pc%liqHeatCapacity%B = -1.5747E02;
!!$    pc%liqHeatCapacity%C = 5.1853E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 131.65, 9.8356E04; 250.00, 1.0314E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 5.1480E04;
!!$    pc%IG_HeatCapacity%B = 1.4420E05;
!!$    pc%IG_HeatCapacity%C = 1.6034E03;
!!$    pc%IG_HeatCapacity%D = 7.7470E04;
!!$    pc%IG_HeatCapacity%E = 7.2540E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 5.4363E04; 1500.00, 1.6581E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.1570E-01;
!!$    pc%SecondViralCoef%B = -1.1010E02;
!!$    pc%SecondViralCoef%C = -5.2550E06;
!!$    pc%SecondViralCoef%D = -4.6140E18;
!!$    pc%SecondViralCoef%E = 7.6350E20;
!!$    pc%SecondViralCoef%range(:,1) = (/ 200.00, -1.4028E00; 1500.00, 4.0743E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.0620E01;
!!$    pc%liqViscosity%B = 4.4899E02;
!!$    pc%liqViscosity%C = 8.3967E-05;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 131.65, 7.3982E-04; 248.31, 1.4903E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.6800E-06;
!!$    pc%vapViscosity%B = 3.9750E-01;
!!$    pc%vapViscosity%C = 5.3400E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 131.65, 3.6879E-06; 1000.00, 2.7215E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 3.1174E-01;
!!$    pc%liqThermalConductivity%B = -5.6380E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 131.65, 2.3752E-01; 320.03, 1.3131E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 5.9975E-02;
!!$    pc%vapThermalConductivity%B = 2.6670E-01;
!!$    pc%vapThermalConductivity%C = 1.0186E03;
!!$    pc%vapThermalConductivity%D = 1.0988E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 248.31, 1.1388E-02; 1500.00, 1.9458E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.0960E-02;
!!$    pc%SurfaceTension%B = 1.2286E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 131.65, 3.7336E-02; 400.10, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_dimethylether
!!$
!!$
!!$  subroutine pc_dinbutylether(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'dinbutylether','di-n-butyl-ether','butyl oxide','1,1-oxybisbutane'}; 
!!$    pc%CAS = '142-96-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 8;         
!!$    pc%ChemicalFormula%H = 18;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 130.23;   ! kg/kmol    
!!$    pc%Tcrit = 584.1;              ! K    
!!$    pc%pcrit = 2.4600E6;           ! Pa 
!!$    pc%vcrit = 0.487;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.247; ! -
!!$
!!$
!!$    pc%MeltingPoint = 175.3;               ! K
!!$    pc%TripletPointT = 175.3;              ! K
!!$    pc%TripletPointp = 7.1413E-04;         ! Pa
!!$    pc%NormalBoilingPoint = 414.15;        ! K
!!$    pc%liqMolarVolume = 0.170405;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.3340E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -8.8270E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  5.0140E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.6300E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.9469E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -9.9490E07;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.7790E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.8980E05;         ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.447646;           ! -
!!$    pc%RadiusOfGyration = 4.7970E-10;      ! m
!!$    pc%SolubilityParameter = 1.5830E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 3.9000E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.09192;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.2880E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3968;           ! -
!!$    pc%FlashPoint = 298.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.5 7.6 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 308 341 /); ! K
!!$    pc%AutoignitionTemp =467.15;           ! K
!!$
!!$
!!$    pc%solDensity%A = 7.7300E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 175.30, 7.730E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.5941E-01;
!!$    pc%liqDensity%B = 2.7243E-01;
!!$    pc%liqDensity%C = 5.8410E02;
!!$    pc%liqDensity%D = 2.9932E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 175.30, 6.6071E00; 584.10,2.0534E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A = 7.2227E01;
!!$    pc%vapPressure%B = -7.5376E03;
!!$    pc%vapPressure%C = -7.0596E00;
!!$    pc%vapPressure%D =  9.1442E-18;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 175.30, 7.1414E-04; 584.10, 2.4592E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.9616E07;
!!$    pc%HeatOfVap%B = 3.8833E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 175.30, 5.1902E07; 584.10, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 0;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.7072E05;
!!$    pc%liqHeatCapacity%B = -2.5983E02;
!!$    pc%liqHeatCapacity%C = 9.5427E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 175.30, 2.5450E05; 450.00, 3.4704E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.6122E05;
!!$    pc%IG_HeatCapacity%B = 4.4777E05;
!!$    pc%IG_HeatCapacity%C = 1.6831E03;
!!$    pc%IG_HeatCapacity%D = 2.9180E05;
!!$    pc%IG_HeatCapacity%E = 7.8160E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.6841E05; 1500.00, 5.2145E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.3911E-01;
!!$    pc%SecondViralCoef%B = -3.7344E02;
!!$    pc%SecondViralCoef%C = -7.7353E06;
!!$    pc%SecondViralCoef%D = -1.3080E20;
!!$    pc%SecondViralCoef%E = 4.4188E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 292.05, -6.2304E00; 1500.00, 6.7226E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = 1.0027E01;
!!$    pc%liqViscosity%B = 2.0600E02;
!!$    pc%liqViscosity%C = -3.1607E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 175.30, 5.9306E-03; 414.15, 1.9890E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 7.7147E-08;
!!$    pc%vapViscosity%B = 7.9906E-01;
!!$    pc%vapViscosity%C = 8.0765E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 175.30, 3.2784E-06; 1000.00, 1.7814E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.9418E-01;
!!$    pc%liqThermalConductivity%B = -2.2246E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 175.30, 1.5518E-01; 523.15, 7.7800E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 3.2694E-03;
!!$    pc%vapThermalConductivity%B = 5.8633E-01;
!!$    pc%vapThermalConductivity%C = 1.2599E03;
!!$    pc%vapThermalConductivity%D = 3.0089E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 323.15, 1.2440E-02; 1000.00, 7.3297E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1346E-02;
!!$    pc%SurfaceTension%B = 1.1604E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 175.30, 3.3937E-02; 584.10, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_dinbutylether
!!$
!!$
!!$  subroutine pc_dinpropylether(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'dinpropylether','di-n-propyl-ether','n-propyl ether'}; 
!!$    pc%CAS = '1 11-43-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 8;         
!!$    pc%ChemicalFormula%H = 18;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 102.177;  ! kg/kmol    
!!$    pc%Tcrit = 530.6;              ! K    
!!$    pc%pcrit = 3.0280E6;           ! Pa 
!!$    pc%vcrit = 0.382;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.262; ! -
!!$
!!$
!!$    pc%MeltingPoint = 149.95;              ! K
!!$    pc%TripletPointT = 149.95;             ! K
!!$    pc%TripletPointp = 7.6265E-04;         ! Pa
!!$    pc%NormalBoilingPoint = 363.23;        ! K
!!$    pc%liqMolarVolume = 0.137649;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.9290E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.0550E08;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.2250E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.0770E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.7250E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.1204E07;       ! J/kmol
!!$    pc%StandardHeatOfFormation = -3.2882E08;        ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.2388E05;         ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.368756;          ! -
!!$    pc%RadiusOfGyration = 3.9140E-10;      ! m
!!$    pc%SolubilityParameter = 1.5550E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 4.0400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.07196;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.0220E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.378;            ! -
!!$    pc%FlashPoint = 294.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.2 9.5 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 264 301 /); ! K
!!$    pc%AutoignitionTemp =488.15;           ! K
!!$
!!$
!!$    pc%solDensity%A = 9.8800E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 149.95, 9.8800E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 7.2390E-01;
!!$    pc%liqDensity%B = 2.7653E-01;
!!$    pc%liqDensity%C = 5.3060E02;
!!$    pc%liqDensity%D = 2.7940E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 149.95, 8.4474E00; 530.60, 2.6178E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A = 6.5935E01;
!!$    pc%vapPressure%B = -6.1985E03;
!!$    pc%vapPressure%C = -6.3423E00;
!!$    pc%vapPressure%D =  2.0909E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 149.95, 7.6265E-04; 530.60, 3.0416E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.8430E07;
!!$    pc%HeatOfVap%B = 3.6670E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 149.95, 4.2877E07; 530.6, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.8140E04;
!!$    pc%solHeatCapacity%B = 1.6964E03;
!!$    pc%solHeatCapacity%C = -8.0340E00;
!!$    pc%solHeatCapacity%D = 1.8100E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 30.00, 2.6010E04; 149.95, 1.1662E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.9374E05;
!!$    pc%liqHeatCapacity%B = -1.2470E02;
!!$    pc%liqHeatCapacity%C = 7.3500E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 149.95, 1.9157E05; 363.23, 2.4542E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2030E05;
!!$    pc%IG_HeatCapacity%B = 3.4123E05;
!!$    pc%IG_HeatCapacity%C = 1.5455E03;
!!$    pc%IG_HeatCapacity%D = 2.0480E05;
!!$    pc%IG_HeatCapacity%E = 7.2320E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.2805E05; 1500.00, 4.0072E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.4620E-01;
!!$    pc%SecondViralCoef%B = -2.6016E02;
!!$    pc%SecondViralCoef%C = -3.9596E07;
!!$    pc%SecondViralCoef%D = 1.0688E19;
!!$    pc%SecondViralCoef%E = -1.1007E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 265.50, -4.0957E00; 1500.00, 6.1028E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -8.2119E00;
!!$    pc%liqViscosity%B = 8.9352E02;
!!$    pc%liqViscosity%C = -4.5940E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 149.95, 1.0516E-02; 363.23, 2.1174E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 4.0850E-08;
!!$    pc%vapViscosity%B = 8.9530E-01;
!!$    pc%vapViscosity%C = 0;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 149.95,3.6251E-06; 1000.00, 1.9819E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1370E-01;
!!$    pc%liqThermalConductivity%B = -2.9200E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 149.95, 1.6991E-01; 420.00, 9.1060E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 7.9334E-02;
!!$    pc%vapThermalConductivity%B = 2.8330E-01;
!!$    pc%vapThermalConductivity%C = 4.5613E03;
!!$    pc%vapThermalConductivity%D = 1.3546E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 363.23, 1.7690E-02; 1000.00, 8.1192E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.4950E-02;
!!$    pc%SurfaceTension%B = 1.2254E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 149.95, 3.6578E-02; 530.60, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_dinpropylether
!!$
!!$
!!$  subroutine pc_ethane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'ethane','bimethyl','ethyl hydride','methylmethane'}; 
!!$    pc%CAS = '74-64- 0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 2;         
!!$    pc%ChemicalFormula%H = 6;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 30.070;   ! kg/kmol    
!!$    pc%Tcrit = 305.32;             ! K    
!!$    pc%pcrit = 4.8720E6;           ! Pa 
!!$    pc%vcrit = 0.14550;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.279; ! -
!!$
!!$
!!$    pc%MeltingPoint = 90.35;              ! K
!!$    pc%TripletPointT = 90.35;             ! K
!!$    pc%TripletPointp = 1.1300E00;          ! Pa
!!$    pc%NormalBoilingPoint = 184.55;        ! K
!!$    pc%liqMolarVolume = 0.055229;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -8.3820E07;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -3.1920E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.2912E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  2.6590E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.4286E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.0995;            ! -
!!$    pc%RadiusOfGyration = 1.8260E-10;      ! m
!!$    pc%SolubilityParameter = 1.2400E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0.0;                 ! C^m
!!$    pc%VanDerWaalsVolume = 0.02734;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 4.2400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.1849;           ! -
!!$    pc%FlashPoint = MAX_REAL_WP;                   ! K
!!$    pc%FlammabilityLimits = (/ 2.9 13.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 745.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 2.5263E01;
!!$    pc%solDensity%B = -1.0950E-02;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 23.15, 2.5010E01; 90.35, 2.4274E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.9122E00;
!!$    pc%liqDensity%B = 2.7937E-01;
!!$    pc%liqDensity%C = 3.0532E02;
!!$    pc%liqDensity%D = 2.9187E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 90.35, 2.1640E01; 305.32, 6.8447E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.1857E01;
!!$    pc%vapPressure%B = -2.5987E03;
!!$    pc%vapPressure%C = -5.1283E00;
!!$    pc%vapPressure%D =  1.4913E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 90.35, 1.1273E00; 305.32, 4.8522E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 2.1091E07;
!!$    pc%HeatOfVap%B = 6.0646E-01;
!!$    pc%HeatOfVap%C = -5.5492E-01;
!!$    pc%HeatOfVap%D = 3.2799E-01;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 90.35, 1.7879E07; 305.32, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 8.1584E03;
!!$    pc%solHeatCapacity%B = -1.4778E03;
!!$    pc%solHeatCapacity%C = 9.7183E01;
!!$    pc%solHeatCapacity%D = -1.6114E00;
!!$    pc%solHeatCapacity%E = 9.1233E-03;
!!$    pc%solHeatCapacity%range(:,1) = (/ 15.53, 3.1420E03; 63.03, 4.1591E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 4.4009E01;
!!$    pc%liqHeatCapacity%B = 8.9718E04;
!!$    pc%liqHeatCapacity%C = 9.1877E02;
!!$    pc%liqHeatCapacity%D = -1.8860E03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 92.00, 6.8554E04; 290.00, 1.2444E05 /);
!!$    pc%liqHeatCapacity%eqn = 114;    
!!$
!!$    pc%IG_HeatCapacity%A = 4.0326E04;
!!$    pc%IG_HeatCapacity%B = 1.3422E05;
!!$    pc%IG_HeatCapacity%C = 1.6555E03;
!!$    pc%IG_HeatCapacity%D = 7.3223E04;
!!$    pc%IG_HeatCapacity%E = 7.5287E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 4.2557E04; 1500.00, 1.4562E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 8.0950E-02;
!!$    pc%SecondViralCoef%B = -6.1710E01;
!!$    pc%SecondViralCoef%C = -1.4350E06;
!!$    pc%SecondViralCoef%D = 6.7600E16;
!!$    pc%SecondViralCoef%E = -9.7400E18;
!!$    pc%SecondViralCoef%range(:,1) = (/ 152.71, -7.1318E-01; 1500.00, 3.9385E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -3.4134E00;
!!$    pc%liqViscosity%B = 1.9705E02;
!!$    pc%liqViscosity%C = -1.2193E00;
!!$    pc%liqViscosity%D = -9.2023E-26;
!!$    pc%liqViscosity%E = 1E01;
!!$    pc%liqViscosity%range(:,1) = (/ 90.35, 1.2020E-03; 300.00, 3.5196E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.1005E-08;
!!$    pc%vapViscosity%B = 9.7620E-01;
!!$    pc%vapViscosity%C = 2.3139E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) =(/ 90.35, 2.6429E-06; 1000.00, 2.5827E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 3.5758E-01;
!!$    pc%liqThermalConductivity%B = -1.1458E-03;
!!$    pc%liqThermalConductivity%C = 6.1866E-07;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 90.35, 2.5911E-01; 300.00, 6.9519E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 7.3869E-05;
!!$    pc%vapThermalConductivity%B = 1.1689E00;
!!$    pc%vapThermalConductivity%C = 5.0073E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 164.55, 8.8628E-03; 1000.00, 1.5807E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 4.8643E-02;
!!$    pc%SurfaceTension%B = 1.1981E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 90.35, 3.1949E-02; 305.32, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_ethane
!!$
!!$
!!$  subroutine pc_ethanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'ethanol','Ethanol'}; 
!!$    pc%CAS = '64-17-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 2;         
!!$    pc%ChemicalFormula%H = 6;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 46.069;   ! kg/kmol    
!!$    pc%Tcrit = 513.92;             ! K    
!!$    pc%pcrit = 6.1480E6;           ! Pa 
!!$    pc%vcrit = 0.16700;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.240; ! -
!!$
!!$
!!$    pc%MeltingPoint = 159.05;              ! K
!!$    pc%TripletPointT = 159.05;             ! K
!!$    pc%TripletPointp = 4.8459E-4;          ! Pa
!!$    pc%NormalBoilingPoint = 351.44;        ! K
!!$    pc%liqMolarVolume = 0.058618;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.3495E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.6785E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.8064E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  4.9310E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -1.2350E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.6452;            ! -
!!$    pc%RadiusOfGyration = 2.2590E-10;      ! m
!!$    pc%SolubilityParameter = 2.6130E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.6400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.03194;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 4.9300E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3594;           ! -
!!$    pc%FlashPoint = 286.00;                ! K
!!$    pc%FlammabilityLimits = (/ 4.3 19.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 696.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 2.29001E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 159.05, 2.29E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.6480E00;
!!$    pc%liqDensity%B = 2.7627E-01;
!!$    pc%liqDensity%C = 5.1391E02;
!!$    pc%liqDensity%D = 2.3310E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 159.05, 1.9413E01; 513.92, 5.9652E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.4475E01;
!!$    pc%vapPressure%B = -7.1643E03;
!!$    pc%vapPressure%C = -7.3270E00;
!!$    pc%vapPressure%D =  3.1340E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 159.05, 4.8459E-04; 513.92, 6.1171E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.6900E07;
!!$    pc%HeatOfVap%B = 3.3590E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 159.05, 5.0245E07; 513.92, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.3500E04;
!!$    pc%solHeatCapacity%B = 1.1755E03;
!!$    pc%solHeatCapacity%C = -8.0430E00;
!!$    pc%solHeatCapacity%D = 2.3700E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 25.00, 1.1231E04; 159.00, 6.5336E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.0264E05;
!!$    pc%liqHeatCapacity%B = -1.3963E02;
!!$    pc%liqHeatCapacity%C = -3.0341E-02;
!!$    pc%liqHeatCapacity%D = 2.0386E-03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 159.05, 8.7867E04; 390.00, 1.6450E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 4.8200E04;
!!$    pc%IG_HeatCapacity%B = 1.4577E05;
!!$    pc%IG_HeatCapacity%C = 1.6628E03;
!!$    pc%IG_HeatCapacity%D = 9.3900E04;
!!$    pc%IG_HeatCapacity%E = 7.4470E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 5.2236E04; 1500.00, 1.6576E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.4000E-02;
!!$    pc%SecondViralCoef%B = -5.5700E01;
!!$    pc%SecondViralCoef%C = -1.2900E07;
!!$    pc%SecondViralCoef%D = -6.4900E19;
!!$    pc%SecondViralCoef%E = -2.5800E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 309.72, -2.3203E00; 2581.00, 2.1669E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = 7.8750E00;
!!$    pc%liqViscosity%B = 7.8198E02;
!!$    pc%liqViscosity%C = -3.0418E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 200.00, 1.3148E-02; 440.00, 1.4160E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.0613E-07;
!!$    pc%vapViscosity%B = 8.0660E-01;
!!$    pc%vapViscosity%C = 5.2700E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 200.00, 6.0295E-06; 1000.00, 2.6505E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.4680E-01;
!!$    pc%liqThermalConductivity%B = -2.6400E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 159.05, 2.0481E-01; 353.15, 1.5357E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -1.0109E-02;
!!$    pc%vapThermalConductivity%B = 6.4750E-01;
!!$    pc%vapThermalConductivity%C = -7.3320E03;
!!$    pc%vapThermalConductivity%D = -2.6800E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 293.15, 1.4747E-02; 1000.0, 1.3417E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 3.7640E-02;
!!$    pc%SurfaceTension%B = -2.1570E-05;
!!$    pc%SurfaceTension%C = -1.0250E-07;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 273.15, 2.4101E-02; 503.15, 8.3816E-04 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_ethanol
!!$
!!$
!!$  subroutine pc_gammavalerolactone(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'gammavalerolactone','gamma-pentalactone','4-pentanolide'}; 
!!$    pc%CAS = '108-29-2';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 8;
!!$    pc%ChemicalFormula%O = 2;
!!$
!!$    pc%MolecularWeight = 100.117;  ! kg/kmol    
!!$    pc%Tcrit = 727.00;             ! K    
!!$    pc%pcrit = 4.8400E6;           ! Pa 
!!$    pc%vcrit = 0.21900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.223; ! -
!!$
!!$
!!$    pc%MeltingPoint = 242.15;              ! K
!!$    pc%TripletPointT = 242.15;             ! K
!!$    pc%TripletPointp = 3.1332E-01;         ! Pa
!!$    pc%NormalBoilingPoint = 480.65;        ! K
!!$    pc%liqMolarVolume = 0.095430;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -4.1800E08;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -2.9700E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.5000E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                 ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.4929E09;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;               ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;               ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4030;            ! -
!!$    pc%RadiusOfGyration = 3.4319E-10;      ! m
!!$    pc%SolubilityParameter = 2.287E04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 1.4343E-29;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05631;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.5900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4330;           ! -
!!$    pc%FlashPoint = 354.00;                ! K
!!$    pc%FlammabilityLimits = (/ MAX_REAL_WP 10.6 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = MAX_REAL_WP;             ! K
!!$
!!$
!!$    pc%solDensity%A = 1.2900E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 242.15, 1.2900E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.8680E-01;
!!$    pc%liqDensity%B = 2.7532E-01;
!!$    pc%liqDensity%C = 7.2700E02;
!!$    pc%liqDensity%D = 3.4900E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 242.15, 1.0983E01; 727.00, 3.5842E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  8.3640E01;
!!$    pc%vapPressure%B = -8.7843E03;
!!$    pc%vapPressure%C = -8.8878E00;
!!$    pc%vapPressure%D =  4.5338E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 242.15, 3.1332E-01; 727.00, 4.8391E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.1320E07;
!!$    pc%HeatOfVap%B = 2.9680E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 242.15, 5.4373E07; 727.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 1.4830E05;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 242.15, 1.4830E05; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 9.5800E04;
!!$    pc%liqHeatCapacity%B = 2.7750E02;
!!$    pc%liqHeatCapacity%C = 0;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 363.50, 1.9667E05; 543.50, 2.4662E05 /);
!!$    pc%liqHeatCapacity%eqn = 107;    
!!$
!!$    pc%IG_HeatCapacity%A = 8.9700E04;
!!$    pc%IG_HeatCapacity%B = 1.8340E05;
!!$    pc%IG_HeatCapacity%C = 6.7500E02;
!!$    pc%IG_HeatCapacity%D = 9.3000E04;
!!$    pc%IG_HeatCapacity%E = 2.2500E03;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 300.00 , 1.3190E05; 1500.00, 2.9902E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.1157E-01;
!!$    pc%SecondViralCoef%B = -2.9963E02;
!!$    pc%SecondViralCoef%C = -8.9927E07;
!!$    pc%SecondViralCoef%D = -2.9311E20;
!!$    pc%SecondViralCoef%E = -2.6964E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 363.50, -3.6900E00; 1500.00, -1.4841E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.8539E00;
!!$    pc%liqViscosity%B = 5.6334E02;
!!$    pc%liqViscosity%C = -1.1445E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 242.15, 1.1024E-03; 572.15, 1.0770E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 8.0120E-08;
!!$    pc%vapViscosity%B = 8.5930E-01;
!!$    pc%vapViscosity%C = 6.1350E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) =(/ 242.15, 7.1501E-06; 1000, 2.8562E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2568E-01;
!!$    pc%liqThermalConductivity%B = -2.1770E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 242.15, 1.7296E-01; 480.65, 1.2104E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 3.8300E-04;
!!$    pc%vapThermalConductivity%B = 8.7240E-01;
!!$    pc%vapThermalConductivity%C = 7.0550E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 480.65, 3.3925E-02; 1000.00, 1.5807E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.6815E-02;
!!$    pc%SurfaceTension%B = 1.3153E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 242.15, 3.9217E-02; 727.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_gammavalerolactone
!!$
!!$
!!$  subroutine pc_isobutane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'isobutane','iso-butane','i-butane'}; 
!!$    pc%CAS = '75-28-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 58.123;   ! kg/kmol    
!!$    pc%Tcrit = 407.80;             ! K    
!!$    pc%pcrit = 3.6400E6;           ! Pa 
!!$    pc%vcrit = 0.25900;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.278; ! -
!!$
!!$
!!$    pc%MeltingPoint = 113.54;              ! K
!!$    pc%TripletPointT = 113.54;             ! K
!!$    pc%TripletPointp = 1.2062E-2;          ! Pa
!!$    pc%NormalBoilingPoint = 261.43;        ! K
!!$    pc%liqMolarVolume = 0.105350;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.3499E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -2.1440E7;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.9550E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  4.5400E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.6490E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -2.0760E07;      ! J/kmol
!!$    pc%StandardHeatOfFormation = - 1.3418E08;      ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.9539E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.1835;            ! -
!!$    pc%RadiusOfGyration = 2.9480E-10;      ! m
!!$    pc%SolubilityParameter = 1.2560E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 4.4000E-31;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.04779;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 6.9300E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3175;           ! -
!!$    pc%FlashPoint = MAX_REAL_WP;                   ! K
!!$    pc%FlammabilityLimits = (/ 1.8 8.4 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 191 212 /); ! K
!!$    pc%AutoignitionTemp = 733.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.4700E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 113.54, 1.4700E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.0631E00;
!!$    pc%liqDensity%B = 2.7506E-01;
!!$    pc%liqDensity%C = 4.0780E02;
!!$    pc%liqDensity%D = 2.7580E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 113.54, 1.2574E01; 407.80, 3.8650E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.0843E02;
!!$    pc%vapPressure%B = -5.0399E03;
!!$    pc%vapPressure%C = -1.5012E01;
!!$    pc%vapPressure%D =  2.2725E-02;
!!$    pc%vapPressure%E =  1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 113.54, 1.2062E-02; 407.80, 3.6299E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 3.1880E07;
!!$    pc%HeatOfVap%B = 3.9006E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 113.54, 2.8070E07; 407.80, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 1.2546E05;
!!$    pc%solHeatCapacity%B = -2.1580E03;
!!$    pc%solHeatCapacity%C = 1.5616E01;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 77.00, 5.1881E04; 107.00, 7.3342E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.7237E05;
!!$    pc%liqHeatCapacity%B = -1.7839E03;
!!$    pc%liqHeatCapacity%C = 1.4759E01;
!!$    pc%liqHeatCapacity%D = -4.7909E-02;
!!$    pc%liqHeatCapacity%E = 5.8050E-05;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 113.54, 9.9613E04; 380.00, 2.0725E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 6.5490E04;
!!$    pc%IG_HeatCapacity%B = 2.4776E05;
!!$    pc%IG_HeatCapacity%C = 1.5870E03;
!!$    pc%IG_HeatCapacity%D = 1.5750E05;
!!$    pc%IG_HeatCapacity%E = -7.0699E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 7.2181E04; 1500.00, 2.6656E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.4990E-01;
!!$    pc%SecondViralCoef%B = -1.3056E02;
!!$    pc%SecondViralCoef%C = -8.7980E06;
!!$    pc%SecondViralCoef%D = -9.9500E17;
!!$    pc%SecondViralCoef%E = 5.2500E19;
!!$    pc%SecondViralCoef%range(:,1) = (/ 204.70, -1.7532E00; 1500.00, 6.0253E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.3912E01;
!!$    pc%liqViscosity%B = 7.9709E02;
!!$    pc%liqViscosity%C = 4.5508E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 110.00 , 1.0716E-02; 310.95, 1.5877E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.0871E-07;
!!$    pc%vapViscosity%B = 7.8135E-01;
!!$    pc%vapViscosity%C = 7.0639E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 150.00, 3.7065E-06; 1000.00, 2.2422E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0455E-01;
!!$    pc%liqThermalConductivity%B = -3.6589E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 113.54, 1.6301E-01; 400.00, 5.8194E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 8.9772E-02;
!!$    pc%vapThermalConductivity%B = 1.8501E-01;
!!$    pc%vapThermalConductivity%C = 6.3923E02;
!!$    pc%vapThermalConductivity%D = 1.1147E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 261.43, 1.2726E-02; 1000.0, 1.1701E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1359E-02;
!!$    pc%SurfaceTension%B = 1.2532E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 113.54 , 3.4121E-02; 407.80, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_isobutane
!!$
!!$
!!$  subroutine pc_isoheptane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'isoheptane','iso-heptane'}; 
!!$    pc%CAS = '591-76-4';
!!$    ! components
!!$    pc%ChemicalFormula%C = 7;         
!!$    pc%ChemicalFormula%H = 16;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 100.204;  ! kg/kmol    
!!$    pc%Tcrit = 530.4;              ! K    
!!$    pc%pcrit = 2.7400E6;           ! Pa 
!!$    pc%vcrit = 0.421;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.262; ! -
!!$
!!$
!!$    pc%MeltingPoint = 154.9;              ! K
!!$    pc%TripletPointT = 154.9;             ! K
!!$    pc%TripletPointp = 4.2989E-03;         ! Pa
!!$    pc%NormalBoilingPoint = 363.199;      ! K
!!$    pc%liqMolarVolume = 0.14869;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.9460E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 3.4700E06;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.2041E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  9.1840E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.4596E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = -2.4762E06;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -2.2949E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 3.2334E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.328006;          ! -
!!$    pc%RadiusOfGyration = 4.1670E-10;      ! m
!!$    pc%SolubilityParameter = 1.4730E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.07848;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.0980E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.38228;          ! -
!!$    pc%FlashPoint = 250.0;                 ! K
!!$    pc%FlammabilityLimits = (/ 1.0 6.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 260 291 /); ! K
!!$    pc%AutoignitionTemp = 553.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 9.1400E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 154.90, 9.1400E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.3993E-01;
!!$    pc%liqDensity%B = 2.6950E-01;
!!$    pc%liqDensity%C = 5.3040E02;
!!$    pc%liqDensity%D = 2.7930E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 154.90, 7.8100E00; 530.40, 2.3745E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.9531E01;
!!$    pc%vapPressure%B = -5.8264E03;
!!$    pc%vapPressure%C = -5.4269E00;
!!$    pc%vapPressure%D =  1.4542E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 154.90, 4.2989E-03; 530.40, 2.7401E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.7126E07;
!!$    pc%HeatOfVap%B = 3.6975E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 154.90, 4.1476E07; 530.40, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 3.7436E04;
!!$    pc%solHeatCapacity%B = 4.1437E02;
!!$    pc%solHeatCapacity%C = 1.3235E00;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 80.00, 7.9068E04; 154.90, 1.3342E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.7401E05;
!!$    pc%liqHeatCapacity%B = -1.0578E02;
!!$    pc%liqHeatCapacity%C = 9.0525E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 154.90, 1.7935E05; 363.20, 2.5501E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0600E05;
!!$    pc%IG_HeatCapacity%B = 4.2507E05;
!!$    pc%IG_HeatCapacity%C = 1.5770E03;
!!$    pc%IG_HeatCapacity%D = 2.7580E05;
!!$    pc%IG_HeatCapacity%E = 6.9080E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.1915E05; 1500.00, 4.5148E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.4123E-01;
!!$    pc%SecondViralCoef%B = -2.9400E02;
!!$    pc%SecondViralCoef%C = -3.8010E07;
!!$    pc%SecondViralCoef%D = -4.3597E19;
!!$    pc%SecondViralCoef%E = 6.3133E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 265.20, -3.7141E00; 1500.00, 3.3966E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.2224E01;
!!$    pc%liqViscosity%B = 1.0213E03;
!!$    pc%liqViscosity%C = 1.5187E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 154.90, 7.7129E-03; 363.20, 2.0000E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.0132E-06;
!!$    pc%vapViscosity%B = 4.5610E-01;
!!$    pc%vapViscosity%C = 3.5980E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 154.90, 3.0414E-06; 1000.00, 1.7399E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.9690E-01;
!!$    pc%liqThermalConductivity%B = -2.6920E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 154.90, 1.5520E-01; 420.00, 8.3836E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 5.9926E-05;
!!$    pc%vapThermalConductivity%B = 1.0785E00;
!!$    pc%vapThermalConductivity%C = -8.5924E01;
!!$    pc%vapThermalConductivity%D = 1.5868E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 363.20, 1.7582E-02; 1000.0, 9.6076E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0288E-02;
!!$    pc%SurfaceTension%B = 1.1923E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 154.90, 3.3314E-02; 530.40, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_isoheptane
!!$
!!$
!!$  subroutine pc_isohexane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'isohexane','iso-hexane','2-methylpentane'}; 
!!$    pc%CAS = '107-83-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 86.1772;  ! kg/kmol    
!!$    pc%Tcrit = 497.7;              ! K    
!!$    pc%pcrit = 3.0400E6;           ! Pa 
!!$    pc%vcrit = 0.368;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.27;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 119.55;              ! K
!!$    pc%TripletPointT = 119.55;             ! K
!!$    pc%TripletPointp = 2.0656E-5;          ! Pa
!!$    pc%NormalBoilingPoint = 333.41;        ! K
!!$    pc%liqMolarVolume = 0.132925;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.7455E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -5.3380E6;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.8089E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  6.2680E6;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion =-3.8492E9;     ! J/kmol
!!$    pc%StandardGibbsOfFormation = -8.4999E06;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -2.0464E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.9058E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.279149;          ! -
!!$    pc%RadiusOfGyration = 3.7840E-10;      ! m
!!$    pc%SolubilityParameter = 1.4400E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.06825;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 9.6300E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.36873;          ! -
!!$    pc%FlashPoint = 238.0;                 ! K
!!$    pc%FlammabilityLimits = (/ 1.2 7.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 239 269 /); ! K
!!$    pc%AutoignitionTemp = 579.0;           ! K
!!$
!!$
!!$    pc%solDensity%A = 1.0800E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 119.55,1.0800E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 7.2701E-01;
!!$    pc%liqDensity%B = 2.6754E-01;
!!$    pc%liqDensity%C = 4.9770E02;
!!$    pc%liqDensity%D = 2.8268E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 119.55, 9.2041E00; 497.70, 2.7174E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.3579E01;
!!$    pc%vapPressure%B = -5.0412E03;
!!$    pc%vapPressure%C = -4.6404E00;
!!$    pc%vapPressure%D =  1.9443E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 119.55, 2.0656E-05; 497.70, 3.0439E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.2522E07;
!!$    pc%HeatOfVap%B = 3.8070E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 119.55, 3.8300E07; 497.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.6127E04;
!!$    pc%solHeatCapacity%B = 1.6109E03;
!!$    pc%solHeatCapacity%C = -9.0041E00;
!!$    pc%solHeatCapacity%D = 2.6825E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 8.00,1.0108E04; 119.53,9.3590E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.4222E05;
!!$    pc%liqHeatCapacity%B = -4.7830E01;
!!$    pc%liqHeatCapacity%C = 7.3900E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 119.55, 1.4706E05; 333.41,2.0842E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 9.0300E04;
!!$    pc%IG_HeatCapacity%B = 3.8010E05;
!!$    pc%IG_HeatCapacity%C = 1.6020E03;
!!$    pc%IG_HeatCapacity%D = 2.4530E05;
!!$    pc%IG_HeatCapacity%E = -6.9160E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00,1.0192E05; 1500.00, 3.9617E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.0487E-01;
!!$    pc%SecondViralCoef%B = -2.0373E02;
!!$    pc%SecondViralCoef%C = -3.1708E07;
!!$    pc%SecondViralCoef%D = -4.4141E18;
!!$    pc%SecondViralCoef%E = -3.7253E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 248.73,-3.3958E00; 1500.00, 5.9655E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.2860E01;
!!$    pc%liqViscosity%B = 9.4691E02;
!!$    pc%liqViscosity%C = 2.6191E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 119.55, 2.5060E-02; 333.41, 2.0380E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.1164E-06;
!!$    pc%vapViscosity%B = 4.5370E-01;
!!$    pc%vapViscosity%C = 3.7474E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 119.55, 2.3658E-06; 1000.00, 1.8651E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.9334E-01;
!!$    pc%liqThermalConductivity%B = -2.8038E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 119.55, 1.5982E-01; 389.25, 8.4202E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 6.1119E-05;
!!$    pc%vapThermalConductivity%B = 1.0861E00;
!!$    pc%vapThermalConductivity%C = -5.9592E01;
!!$    pc%vapThermalConductivity%D = 1.4126E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 333.41, 1.6063E-02; 1000.0, 1.0242E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0048E-02;
!!$    pc%SurfaceTension%B = 1.1899E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 119.55, 3.6093E-02; 497.70, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_isohexane
!!$
!!$
!!$  subroutine pc_isooctane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'isooctane','iso-octane','heptane'}; 
!!$    pc%CAS = '107-83-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 8;         
!!$    pc%ChemicalFormula%H = 18;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 114.231;   ! kg/kmol    
!!$    pc%Tcrit = 559.64;              ! K    
!!$    pc%pcrit = 2.4840E6;            ! Pa 
!!$    pc%vcrit = 0.48800;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.261;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 164.16;              ! K
!!$    pc%TripletPointT = 164.16;             ! K
!!$    pc%TripletPointp = 1.0606E-3;          ! Pa
!!$    pc%NormalBoilingPoint = 390.80;        ! K
!!$    pc%liqMolarVolume = 0.164409;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.1535E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = 1.1690E7;             ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  4.5957E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  1.1878E7;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -5.0694E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3772;            ! -
!!$    pc%RadiusOfGyration = 4.5700E-10;      ! m
!!$    pc%SolubilityParameter = 1.5050E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.08871;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.2330E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3926;           ! -
!!$    pc%FlashPoint = 277.15;                ! K
!!$    pc%FlammabilityLimits = (/ 0.9 5.8 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 520.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 8.1000E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 164.16,8.1000E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.2966E-01;
!!$    pc%liqDensity%B = 2.6042E-01;
!!$    pc%liqDensity%C = 5.5964E02;
!!$    pc%liqDensity%D = 2.7016E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 164.16, 6.9236E00; 559.64, 2.0339E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  9.0727E01;
!!$    pc%vapPressure%B = -7.4715E03;
!!$    pc%vapPressure%C = -1.0243E01;
!!$    pc%vapPressure%D =  6.8840E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 164.16, 1.0606E-03; 559.64, 2.4757E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.3089E07;
!!$    pc%HeatOfVap%B = 3.8112E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 164.16, 4.6509E07; 559.64, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -2.0180E04;
!!$    pc%solHeatCapacity%B = 1.8980E03;
!!$    pc%solHeatCapacity%C = -9.1850E00;
!!$    pc%solHeatCapacity%D = 2.3100E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.4291E04; 164.16, 1.4606E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.0489E05;
!!$    pc%liqHeatCapacity%B = -1.2073E02;
!!$    pc%liqHeatCapacity%C = 9.2846E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 164.16, 2.1009E05; 500.00, 3.7664E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2215E05;
!!$    pc%IG_HeatCapacity%B = 4.7270E05;
!!$    pc%IG_HeatCapacity%C = 1.5866E03;
!!$    pc%IG_HeatCapacity%D = 3.1650E05;
!!$    pc%IG_HeatCapacity%E = 6.9960E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 150.00, 1.2460E05; 1500.00, 5.0777E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.7362E-01;
!!$    pc%SecondViralCoef%B = -3.4177E02;
!!$    pc%SecondViralCoef%C = -5.4179E07;
!!$    pc%SecondViralCoef%D = -9.1933E19;
!!$    pc%SecondViralCoef%E = 1.3873E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 279.82, -4.5475E00; 1500.00, 2.9717E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.1342E01;
!!$    pc%liqViscosity%B = 1.0740E03;
!!$    pc%liqViscosity%C = 1.3051E-02;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 164.16, 8.7998E-03; 390.80, 2.0026E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 4.4600E-07;
!!$    pc%vapViscosity%B = 5.5355E-01;
!!$    pc%vapViscosity%C = 2.2223E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 164.16, 3.1903E-06; 1000.00, 1.6704E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.9945E-01;
!!$    pc%liqThermalConductivity%B = -2.5820E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 164.16, 1.5706E-01; 448.00 8.3776E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 5.4211E-05;
!!$    pc%vapThermalConductivity%B = 1.0809E00;
!!$    pc%vapThermalConductivity%C = -1.3266E02;
!!$    pc%vapThermalConductivity%D = 1.7358E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 390.80, 1.9106E-02; 1000.0, 9.1069E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0322E-02;
!!$    pc%SurfaceTension%B = 1.2008E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 164.16, 3.3166E-02; 559.64, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_isooctane
!!$
!!$
!!$  subroutine pc_isopentane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'isopentane','iso-pentane','butane','2-methylbutane'}; 
!!$    pc%CAS = '107-83-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 72.1503;  ! kg/kmol    
!!$    pc%Tcrit = 460.4;              ! K    
!!$    pc%pcrit = 3.3800E6;           ! Pa 
!!$    pc%vcrit = 0.306;              ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.27;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 113.25;              ! K
!!$    pc%TripletPointT = 113.25;             ! K
!!$    pc%TripletPointp = 1.2144E-4;          ! Pa
!!$    pc%NormalBoilingPoint = 300.994;       ! K
!!$    pc%liqMolarVolume  = 0.117065;         ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.5370E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.4050E7;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.4374E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  5.1500E06;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.2395E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -1.4420E07;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -1.7891E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.6041E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.227875;            ! -
!!$    pc%RadiusOfGyration = 3.3240E-10;      ! m
!!$    pc%SolubilityParameter = 1.3870E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 4.3400E-30;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.05802;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.2800E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.35088;          ! -
!!$    pc%FlashPoint = 217.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.3 8.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 216 244 /); ! K
!!$    pc%AutoignitionTemp = 693.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.2600E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 113.25, 1.2600E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.1991E-01;
!!$    pc%liqDensity%B = 2.7815E-01;
!!$    pc%liqDensity%C = 4.6040E02;
!!$    pc%liqDensity%D = 2.8667E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 113.25, 1.0764E01; 460.40, 3.3072E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.1308E01;
!!$    pc%vapPressure%B = -4.9760E03;
!!$    pc%vapPressure%C = -7.7169E00;
!!$    pc%vapPressure%D =  8.7271E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 113.25, 1.2144E-04; 460.40,3.3656E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 3.7593E07;
!!$    pc%HeatOfVap%B = 3.9173E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 113.25, 3.3657E07; 460.40, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -8.1500E03;
!!$    pc%solHeatCapacity%B = 8.5500E02;
!!$    pc%solHeatCapacity%C = 1.4470E01;
!!$    pc%solHeatCapacity%D = -2.7730E-01;
!!$    pc%solHeatCapacity%E = 1.3350E-03;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.2733E04; 113.36, 9.1223E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.0830E05;
!!$    pc%liqHeatCapacity%B = 1.4600E02;
!!$    pc%liqHeatCapacity%C = -2.9200E-01;
!!$    pc%liqHeatCapacity%D = 1.5100E-03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 113.25, 1.2328E05; 310.00, 1.7048E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 7.4600E04;
!!$    pc%IG_HeatCapacity%B = 3.2650E05;
!!$    pc%IG_HeatCapacity%C = 1.5450E03;
!!$    pc%IG_HeatCapacity%D = 1.9230E05;
!!$    pc%IG_HeatCapacity%E = 6.6670E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 8.5462E04; 1500.00, 3.3792E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.8170E-01;
!!$    pc%SecondViralCoef%B = -1.9553E02;
!!$    pc%SecondViralCoef%C = -1.5973E07;
!!$    pc%SecondViralCoef%D = -4.6950E18;
!!$    pc%SecondViralCoef%E =3.6600E20;
!!$    pc%SecondViralCoef%range(:,1) = (/ 230.20, -2.3708E00; 1500.00, 4.6614E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.2596E01;
!!$    pc%liqViscosity%B = 8.8911E02;
!!$    pc%liqViscosity%C = 2.0469E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 150.00, 3.5423E-03; 310.00, 1.9283E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.4344E-08;
!!$    pc%vapViscosity%B = 9.7376E-01;
!!$    pc%vapViscosity%C = -9.1597E01;
!!$    pc%vapViscosity%D = 1.8720E04;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 150.00, 2.6214E-06; 1000.00, 2.1905E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1246E-01;
!!$    pc%liqThermalConductivity%B = -3.3581E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 113.25,1.7443E-01; 368.13, 8.8838E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 8.9680E-04;
!!$    pc%vapThermalConductivity%B = 7.7420E-01;
!!$    pc%vapThermalConductivity%C = 4.5600E02;
!!$    pc%vapThermalConductivity%D = 2.3064E05;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.15, 1.1981E-02; 1000.0, 1.1176E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0876E-02;
!!$    pc%SurfaceTension%B = 1.2066E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 113.25, 3.6188E-02; 460.40, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_isopentane
!!$
!!$
!!$  subroutine pc_methane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'methane','Methane'}; 
!!$    pc%CAS = '74-82-8';
!!$    ! components
!!$    pc%ChemicalFormula%C = 1;         
!!$    pc%ChemicalFormula%H = 4;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 16.043;    ! kg/kmol    
!!$    pc%Tcrit = 190.56;              ! K    
!!$    pc%pcrit = 4.5990E06;           ! Pa 
!!$    pc%vcrit = 0.09860;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.286;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 90.69;              ! K
!!$    pc%TripletPointT = 90.69;             ! K
!!$    pc%TripletPointp = 1.1696E04;         ! Pa
!!$    pc%NormalBoilingPoint = 111.66;       ! K
!!$    pc%liqMolarVolume  = 0.037969;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -7.4520E07;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -5.0490E07;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  1.8627E05;            ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  9.4140E05;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -8.0262E08;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.0115;            ! -
!!$    pc%RadiusOfGyration = 1.1180E-10;      ! m
!!$    pc%SolubilityParameter = 1.1600E04;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.01705;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 2.8800E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.0004;           ! -
!!$    pc%FlashPoint = MAX_REAL_WP;                   ! K
!!$    pc%FlammabilityLimits = (/ 5.0 15.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 810.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 3.3022E01;
!!$    pc%solDensity%B = -1.5870E-02;
!!$    pc%solDensity%C = -1.5500E-04;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 23.15, 3.2572E01; 90.67, 3.0309E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 2.9214E00;
!!$    pc%liqDensity%B = 2.8976E-01;
!!$    pc%liqDensity%C = 1.9056E02;
!!$    pc%liqDensity%D = 2.8818E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 90.69, 2.8180E01; 190.56, 1.0082E01 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  3.9205E01;
!!$    pc%vapPressure%B = -1.3244E03;
!!$    pc%vapPressure%C = -3.4366E00;
!!$    pc%vapPressure%D =  3.1019E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 90.69, 1.1687E04; 190.56, 4.5897E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 1.0194E07;
!!$    pc%HeatOfVap%B = 2.6087E-01;
!!$    pc%HeatOfVap%C = -1.4694E-01;
!!$    pc%HeatOfVap%D = 2.2154E-01;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 90.69, 8.7235E06; 190.56, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -2.5500E03;
!!$    pc%solHeatCapacity%B = 1.2499E03;
!!$    pc%solHeatCapacity%C = -1.4170E01;
!!$    pc%solHeatCapacity%D = 6.6480E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 22.85, 1.9405E04; 90.67, 4.3841E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 6.5708E01;
!!$    pc%liqHeatCapacity%B = 3.8883E04;
!!$    pc%liqHeatCapacity%C = -2.5795E02;
!!$    pc%liqHeatCapacity%D = 6.1407E02;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 90.69, 5.3605E04; 190.00, 1.4978E06 /);
!!$    pc%liqHeatCapacity%eqn = 114;    
!!$
!!$    pc%IG_HeatCapacity%A = 3.3298E04;
!!$    pc%IG_HeatCapacity%B = 7.9933E04;
!!$    pc%IG_HeatCapacity%C = 2.0869E03;
!!$    pc%IG_HeatCapacity%D = 4.1602E04;
!!$    pc%IG_HeatCapacity%E = 9.9196E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 50.00, 3.3298E04; 1500.00, 8.8904E04 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 5.1075E-02;
!!$    pc%SecondViralCoef%B = -2.5181E01;
!!$    pc%SecondViralCoef%C = -2.5601E05;
!!$    pc%SecondViralCoef%D = 5.9777E15;
!!$    pc%SecondViralCoef%E =-5.7697E17;
!!$    pc%SecondViralCoef%range(:,1) = (/ 110.83, -3.3028E-01; 1500.00, 3.4212E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -6.1572E00;
!!$    pc%liqViscosity%B = 1.7815E02;
!!$    pc%liqViscosity%C = -9.5239E-01;
!!$    pc%liqViscosity%D = -9.0606E-24;
!!$    pc%liqViscosity%E = 1.0E01;
!!$    pc%liqViscosity%range(:,1) = (/ 90.69, 2.0634E-04; 188.00, 2.2625E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 5.2546E-07;
!!$    pc%vapViscosity%B = 5.9006E-01;
!!$    pc%vapViscosity%C = 1.0567E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 90.69, 3.4684E-06; 1000.00, 2.7996E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 4.1768E-01;
!!$    pc%liqThermalConductivity%B = -2.4528E-03;
!!$    pc%liqThermalConductivity%C = 3.5588E-06;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 90.69,  2.2451E-01; 180.00, 9.1481E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 6.3252E03;
!!$    pc%vapThermalConductivity%B = 4.3041E-01;
!!$    pc%vapThermalConductivity%C = 7.7040E08;
!!$    pc%vapThermalConductivity%D = -3.8725E10;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 90.69, 1.1625E-02; 1000.0, 1.6904E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 3.6557E-02;
!!$    pc%SurfaceTension%B = 1.1466E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 90.69, 1.7428E-02; 190.56, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_methane
!!$
!!$
!!$  subroutine pc_methanol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'methanol','Methanol','methyl alcohol'}; 
!!$    pc%CAS = '67-56-1';
!!$    ! components
!!$    pc%ChemicalFormula%C = 1;         
!!$    pc%ChemicalFormula%H = 4;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 32.042;    ! kg/kmol    
!!$    pc%Tcrit = 512.64;              ! K    
!!$    pc%pcrit = 8.0970E06;           ! Pa 
!!$    pc%vcrit = 0.11800;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.224;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 175.47;             ! K
!!$    pc%TripletPointT = 175.47;            ! K
!!$    pc%TripletPointp = 1.1147E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 337.85;       ! K
!!$    pc%liqMolarVolume  = 0.040581;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.0094E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.6232E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  2.3988E05;            ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  3.2150E06;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -6.38200E8;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5640;            ! -
!!$    pc%RadiusOfGyration = 1.5520E-10;      ! m
!!$    pc%SolubilityParameter = 2.9590E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.6700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.02171;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 3.5800E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3265;           ! -
!!$    pc%FlashPoint = 284.00;                ! K
!!$    pc%FlammabilityLimits = (/ 7.3 36.0 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 737.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 3.0585E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 163.15, 3.0585E01; 90.67, 3.0309E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 2.2880E00;
!!$    pc%liqDensity%B = 2.6850E-01;
!!$    pc%liqDensity%C = 5.1264E02;
!!$    pc%liqDensity%D = 2.4530E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 175.47, 2.7912E01; 512.64, 8.5214E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  8.1768E01;
!!$    pc%vapPressure%B = -6.8760E03;
!!$    pc%vapPressure%C = -8.7078E00;
!!$    pc%vapPressure%D =  7.1926E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 175.47, 1.1147E-01; 512.64, 8.1402E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.2390E07;
!!$    pc%HeatOfVap%B = 3.6820E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 175.47, 4.4900E07; 512.64, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.2706E04;
!!$    pc%solHeatCapacity%B = 9.5878E02;
!!$    pc%solHeatCapacity%C = -5.2332E00;
!!$    pc%solHeatCapacity%D = 1.3152E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 25.00, 8.1982E03; 150.00, 5.7752E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.0580E05;
!!$    pc%liqHeatCapacity%B = -3.6223E02;
!!$    pc%liqHeatCapacity%C = 9.3790E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 175.47, 7.1117E04; 400.00, 1.1097E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 3.9252E04;
!!$    pc%IG_HeatCapacity%B = 8.7900E04;
!!$    pc%IG_HeatCapacity%C = 1.9165E03;
!!$    pc%IG_HeatCapacity%D = 5.3654E03;
!!$    pc%IG_HeatCapacity%E = 8.9670E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 3.9802E04; 1500.00, 1.0533E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.3000E-02;
!!$    pc%SecondViralCoef%B = -1.4400E01;
!!$    pc%SecondViralCoef%C = -2.0600E07;
!!$    pc%SecondViralCoef%D = 7.6000E19;
!!$    pc%SecondViralCoef%E =-4.6200E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 298.15, -2.0770E00; 1500.00, -2.7019E-03 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.5317E01;
!!$    pc%liqViscosity%B = 1.7892E03;
!!$    pc%liqViscosity%C = 2.0690E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 175.47, 1.1928E-02; 337.85, 3.4423E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.0663E-07;
!!$    pc%vapViscosity%B = 6.9655E-01;
!!$    pc%vapViscosity%C = 2.0500E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 240.00, 7.5232E-06; 1000.00, 3.1281E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.8370E-01;
!!$    pc%liqThermalConductivity%B = -2.8100E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 175.47, 2.3439E-01; 337.85, 1.8876E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -7.7630E00;
!!$    pc%vapThermalConductivity%B = 1.0279E00;
!!$    pc%vapThermalConductivity%C = -7.4360E07;
!!$    pc%vapThermalConductivity%D = 6.7700E09;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.00, 1.3652E-02; 1000.0, 1.3927E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 3.5130E-02;
!!$    pc%SurfaceTension%B = -7.0400E-06;
!!$    pc%SurfaceTension%C = -1.2160E-07;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 273.10, 2.4138E-02; 503.15, 8.0358E-04 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_methanol
!!$
!!$
!!$  subroutine pc_methyltertbutylether(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'methyltertbutylether','methyl tert-buty lether','propane'}; 
!!$    pc%CAS = '1634-04-4';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 88.150;    ! kg/kmol    
!!$    pc%Tcrit = 497.10;              ! K    
!!$    pc%pcrit = 3.4300E6;            ! Pa 
!!$    pc%vcrit = 0.32900;             ! m^3/kmol 
!!$    pc%CritCompressFactor =0.273;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 164.55;             ! K
!!$    pc%TripletPointT = 164.55;            ! K
!!$    pc%TripletPointp = 5.3566E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 328.35;       ! K
!!$    pc%liqMolarVolume  = 0.119887;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.8350E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.1750E8;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.5780E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  7.6000E06;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.1049E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2661;            ! -
!!$    pc%RadiusOfGyration = 3.1790E-10;      ! m
!!$    pc%SolubilityParameter = 1.5070E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 4.5400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.06171;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 9.0800E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3663;           ! -
!!$    pc%FlashPoint = 245.00;                ! K
!!$    pc%FlammabilityLimits = (/ 2.0 15.1 /);    ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 733.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1400E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 164.55, 1.1400E01; 90.67, 3.0309E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.2157E-01;
!!$    pc%liqDensity%B = 2.7032E-01;
!!$    pc%liqDensity%C = 4.9710E02;
!!$    pc%liqDensity%D = 2.8290E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 164.55, 9.7682E00; 497.10, 3.0392E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.5875E01;
!!$    pc%vapPressure%B = -5.1316E03;
!!$    pc%vapPressure%C = -4.9604E00;
!!$    pc%vapPressure%D =  1.9123E-17;
!!$    pc%vapPressure%E =  6.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 164.55, 5.3566E-01; 497.10 , 3.4106E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.2024E07;
!!$    pc%HeatOfVap%B = 3.7826E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 164.55, 3.6096E07; 497.10, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.2579E04;
!!$    pc%solHeatCapacity%B = 1.4658E03;
!!$    pc%solHeatCapacity%C = -7.7409E00;
!!$    pc%solHeatCapacity%D = 9.6626E-03;
!!$    pc%solHeatCapacity%E = 9.9006E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 11.50, 3.2704E03; 164.55, 1.3466E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.4012E05;
!!$    pc%liqHeatCapacity%B = -9.0000E00;
!!$    pc%liqHeatCapacity%C = 5.6300E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 164.55, 1.5388E05; 328.35, 1.9786E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 9.9328E04;
!!$    pc%IG_HeatCapacity%B = 3.0667E05;
!!$    pc%IG_HeatCapacity%C = 1.7426E03;
!!$    pc%IG_HeatCapacity%D = 2.0764E05;
!!$    pc%IG_HeatCapacity%E = 7.9559E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.0394E05; 1500.00, 3.4321E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.8083E-01;
!!$    pc%SecondViralCoef%B = -1.8999E02;
!!$    pc%SecondViralCoef%C = -2.2566E07;
!!$    pc%SecondViralCoef%D = -9.3110E18;
!!$    pc%SecondViralCoef%E = 9.3054E18;
!!$    pc%SecondViralCoef%range(:,1) = (/ 248.60, -2.8896E00; 1500.00, 4.6565E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -8.4217E00;
!!$    pc%liqViscosity%B = 8.5730E02;
!!$    pc%liqViscosity%C = -4.3178E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 180.00, 2.7363E-03; 449.93, 1.0579E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.5439E-07;
!!$    pc%vapViscosity%B = 7.3601E-01;
!!$    pc%vapViscosity%C = 1.0824E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 164.55, 3.9839E-06; 1000.00, 2.2491E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2555E-01;
!!$    pc%liqThermalConductivity%B = -3.2780E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 164.55, 1.7161E-01; 328.35, 1.1792E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 3.1710E-04;
!!$    pc%vapThermalConductivity%B = 8.8303E-01;
!!$    pc%vapThermalConductivity%C = 4.6940E02;
!!$    pc%vapThermalConductivity%D = 8.6314E04;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.00, 1.1584E-02; 1000.00, 9.0857E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.9810E-02;
!!$    pc%SurfaceTension%B = 1.2352E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 164.55, 3.6402E-02; 497.10, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_methyltertbutylether
!!$
!!$
!!$  subroutine pc_nbutane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'nbutane','n-butane','methylethylmethane'}; 
!!$    pc%CAS = '106-97-8';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 58.123;    ! kg/kmol    
!!$    pc%Tcrit = 425.12;              ! K    
!!$    pc%pcrit = 3.7960E6;            ! Pa 
!!$    pc%vcrit = 0.25500;             ! m^3/kmol 
!!$    pc%CritCompressFactor =0.274;   ! -
!!$
!!$
!!$    pc%MeltingPoint = 134.86;             ! K
!!$    pc%TripletPointT = 134.86;            ! K
!!$    pc%TripletPointp = 6.7358E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 272.65;       ! K
!!$    pc%liqMolarVolume  = 0.096484;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.2579E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -1.6700E07;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy =  3.0991E5;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint =  4.6610E06;         ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = - 2.6573E9;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2002;            ! -
!!$    pc%RadiusOfGyration = 2.8860E-10;      ! m
!!$    pc%SolubilityParameter = 1.3700E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.04780;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 6.9400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3292;           ! -
!!$    pc%FlashPoint = MAX_REAL_WP;                   ! K
!!$    pc%FlammabilityLimits = (/ 1.5 9.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 561.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.5194E01;
!!$    pc%solDensity%B = -6.8446E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 23.15, 1.5036E01; 134.86, 1.4271E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.0677E00;
!!$    pc%liqDensity%B = 2.7188E-01;
!!$    pc%liqDensity%C = 4.2512E02;
!!$    pc%liqDensity%D = 2.8688E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 134.86, 1.2620E01; 425.12, 3.9271E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  6.6343E01;
!!$    pc%vapPressure%B = -4.3632E03;
!!$    pc%vapPressure%C = -7.0460E00;
!!$    pc%vapPressure%D =  9.4509E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 134.86, 6.7441E-01; 425.12, 3.7699E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 3.6238E07;
!!$    pc%HeatOfVap%B = 8.3370E-01;
!!$    pc%HeatOfVap%C = -8.2274E-01;
!!$    pc%HeatOfVap%D = 3.9613E-01;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 134.86, 2.8684E07; 425.12, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.8858E04;
!!$    pc%solHeatCapacity%B = 1.5065E03;
!!$    pc%solHeatCapacity%C = -9.8874E00;
!!$    pc%solHeatCapacity%D = 3.2073E-02;
!!$    pc%solHeatCapacity%E = 1.6746E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 22.35, 1.0235E04; 103.54, 6.8653E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 6.4730E01;
!!$    pc%liqHeatCapacity%B = 1.6184E05;
!!$    pc%liqHeatCapacity%C = 9.8341E02;
!!$    pc%liqHeatCapacity%D = -1.4315E03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 134.86, 1.1380E-05; 420.00, 5.0822E05 /);
!!$    pc%liqHeatCapacity%eqn = 114;    
!!$
!!$    pc%IG_HeatCapacity%A = 7.1340E04;
!!$    pc%IG_HeatCapacity%B = 2.4300E05;
!!$    pc%IG_HeatCapacity%C = 1.6300E03;
!!$    pc%IG_HeatCapacity%D = 1.5033E05;
!!$    pc%IG_HeatCapacity%E = 7.3042E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 7.6733E04; 1500.00, 2.6602E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.4900E-01;
!!$    pc%SecondViralCoef%B = -1.4300E02;
!!$    pc%SecondViralCoef%C = -1.0010E07;
!!$    pc%SecondViralCoef%D = -3.0300E18;
!!$    pc%SecondViralCoef%E = 4.0300E20;
!!$    pc%SecondViralCoef%range(:,1) = (/ 212.58,-1.8377E00; 1500.00, 5.0701E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -7.2471E00;
!!$    pc%liqViscosity%B = 5.3482E02;
!!$    pc%liqViscosity%C = -5.7469E-01;
!!$    pc%liqViscosity%D = -4.6625E-27;
!!$    pc%liqViscosity%E = 1E01;
!!$    pc%liqViscosity%range(:,1) = (/ 134.86, 2.2434E-03; 420.00, 3.5664E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.2980E-07;
!!$    pc%vapViscosity%B = 6.9442E-01;
!!$    pc%vapViscosity%C = 2.2766E02;
!!$    pc%vapViscosity%D = -1.4610E04;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 134.86, 3.6738E-06; 1000.00, 2.2947E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.7349E-01;
!!$    pc%liqThermalConductivity%B = -7.1267E-04;
!!$    pc%liqThermalConductivity%C = 5.1555E-07;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 134.86, 1.8676E-01; 400.00, 7.0910E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 5.1094E-02;
!!$    pc%vapThermalConductivity%B = 4.5253E-01;
!!$    pc%vapThermalConductivity%C = 5.4555E03;
!!$    pc%vapThermalConductivity%D = 1.9798E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 272.65, 1.3570E-02; 1000.00, 1.3799E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.1960E-02;
!!$    pc%SurfaceTension%B = 1.2181E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 134.86, 3.2644E-02; 425.12, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_nbutane
!!$
!!$
!!$  subroutine pc_nbutylformate(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'nbutylformate','n-butyl formate','butyl methanoate'}; 
!!$    pc%CAS = '592-84-7';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 2;
!!$
!!$    pc%MolecularWeight = 102.133;   ! kg/kmol    
!!$    pc%Tcrit = 559.0;               ! K    
!!$    pc%pcrit = 3.5100E6;            ! Pa 
!!$    pc%vcrit = 0.336;               ! m^3/kmol 
!!$    pc%CritCompressFactor =0.254;   ! -
!!$
!!$
!!$    pc%MeltingPoint = 181.25;             ! K
!!$    pc%TripletPointT = 181.25;            ! K
!!$    pc%TripletPointp = 0.043691;          ! Pa
!!$    pc%NormalBoilingPoint = 379.25;       ! K
!!$    pc%liqMolarVolume  = 0.115158;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -4.2710E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = -2.8440E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 4.0800E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.3900E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.7099E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = -2.9390E08;      ! J/kmol
!!$    pc%StandardHeatOfFormation = -4.6920E08;       ! J/kmol
!!$    pc%StandardAbsoluteEntropy = 2.9860E05;        ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.387793;          ! -
!!$    pc%RadiusOfGyration = 3.7850E-10;      ! m
!!$    pc%SolubilityParameter = 1.8230E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 6.7700E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.0632;         ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 9.1400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3874;           ! -
!!$    pc%FlashPoint = 291.0;                 ! K
!!$    pc%FlammabilityLimits = (/ 1.7 8.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 284 313 /); ! K
!!$    pc%AutoignitionTemp = 595.0;           ! K
!!$
!!$
!!$    pc%solDensity%A = 1.1400E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 181.25,1.1400E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 7.9185E-01;
!!$    pc%liqDensity%B = 2.6606E-01;
!!$    pc%liqDensity%C = 5.5900E02;
!!$    pc%liqDensity%D = 2.7850E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 181.25, 9.7549E00; 559.00, 2.9762E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.1554E02;
!!$    pc%vapPressure%B = -7.9070E03;
!!$    pc%vapPressure%C = -1.4945E01;
!!$    pc%vapPressure%D =  1.4711E-02;
!!$    pc%vapPressure%E =  1.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 181.25, 4.3690E-02; 559.00, 3.5210E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.5800E07;
!!$    pc%HeatOfVap%B = 4.1300E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 181.25, 4.7461E07; 559.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 0;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.0210E05;
!!$    pc%liqHeatCapacity%B = 3.2795E02;
!!$    pc%liqHeatCapacity%C = 0;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 290.15, 1.9725E05; 379.25, 2.2648E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0596E05;
!!$    pc%IG_HeatCapacity%B = 2.9480E05;
!!$    pc%IG_HeatCapacity%C = 1.8990E03;
!!$    pc%IG_HeatCapacity%D = 2.2542E05;
!!$    pc%IG_HeatCapacity%E = -8.2219E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 298.15, 1.3348E05; 1500.00, 3.3415E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.0700E-01;
!!$    pc%SecondViralCoef%B = -2.3200E02;
!!$    pc%SecondViralCoef%C = -4.0100E07;
!!$    pc%SecondViralCoef%D = -1.7314E20;
!!$    pc%SecondViralCoef%E = 2.9240E20;
!!$    pc%SecondViralCoef%range(:,1) = (/ 279.50, -4.2995E00; 1500.00, 4.0446E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.7801E01;
!!$    pc%liqViscosity%B = 1.4978E03;
!!$    pc%liqViscosity%C = 9.5330E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 250.00, 1.4356E-03; 379.25,2.7721E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 6.2067E-07;
!!$    pc%vapViscosity%B = 5.7800E-01;
!!$    pc%vapViscosity%C = 3.6965E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 181.25, 4.1243E-06; 1000.00, 2.4561E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1050E-01;
!!$    pc%liqThermalConductivity%B = -2.3200E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 181.25, 1.6845E-01; 503.15, 9.3769E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 1.3460E03;
!!$    pc%vapThermalConductivity%B = 9.9860E-01;
!!$    pc%vapThermalConductivity%C = 1.2465E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 379.25, 1.5403E-02; 999.25, 1.0678E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.7032E-02;
!!$    pc%SurfaceTension%B = 1.1073E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 181.25, 3.6953E-02; 559.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_nbutylformate
!!$
!!$
!!$  subroutine pc_ndecane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'ndecane','n-decane','decyl hydride'}; 
!!$    pc%CAS = '124-18-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 10;         
!!$    pc%ChemicalFormula%H = 22;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 142.285;   ! kg/kmol    
!!$    pc%Tcrit = 617.70;              ! K    
!!$    pc%pcrit = 2.1100E6;            ! Pa 
!!$    pc%vcrit = 0.60000;             ! m^3/kmol 
!!$    pc%CritCompressFactor =0.247;   ! -
!!$
!!$
!!$    pc%MeltingPoint = 243.51;             ! K
!!$    pc%TripletPointT = 243.51;            ! K
!!$    pc%TripletPointp = 1.3930;            ! Pa
!!$    pc%NormalBoilingPoint = 447.30;       ! K
!!$    pc%liqMolarVolume  = 0.195827;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.4946E8;             ! J/kmol
!!$    pc%IG_GibbsOfFormation = 3.3180E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 5.4570E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 2.8710E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -6.2942E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4923;            ! -
!!$    pc%RadiusOfGyration = 5.1480E-10;      ! m
!!$    pc%SolubilityParameter = 1.5700E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.10920;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.5040E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4097;           ! -
!!$    pc%FlashPoint = 319.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.7 5.4 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ 284 313 /); ! K
!!$    pc%AutoignitionTemp = 474.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 6.8978E00;
!!$    pc%solDensity%B = -2.4600E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 193.15, 6.4227; 243.15, 6.2988E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 4.2831E-01;
!!$    pc%liqDensity%B = 2.5745E-01;
!!$    pc%liqDensity%C = 6.1770E02;
!!$    pc%liqDensity%D = 2.8912E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 243.15, 5.3811E00; 617.70, 1.6637E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.1273E02;
!!$    pc%vapPressure%B = -9.7496E03;
!!$    pc%vapPressure%C = -1.3245E01;
!!$    pc%vapPressure%D =  7.1266E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 243.51, 1.3930E00; 617.70, 2.0908E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.6126E07;
!!$    pc%HeatOfVap%B = 3.9797E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 243.51, 5.4168E07; 617.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -3.4900E04;
!!$    pc%solHeatCapacity%B = 2.5280E03;
!!$    pc%solHeatCapacity%C = -1.2650E01;
!!$    pc%solHeatCapacity%D = 2.8480E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 19.93, 1.0684E04; 240.00, 2.3689E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.7862E05;
!!$    pc%liqHeatCapacity%B = -1.9791E02;
!!$    pc%liqHeatCapacity%C = 1.0737E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 243.51, 2.9409E05; 460.00, 4.1478E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.6720E05;
!!$    pc%IG_HeatCapacity%B = 5.3530E05;
!!$    pc%IG_HeatCapacity%C = 1.6141E03;
!!$    pc%IG_HeatCapacity%D = 3.7820E05;
!!$    pc%IG_HeatCapacity%E = 7.4200E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.7967E05; 1500.00, 6.0932E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.800E-01;
!!$    pc%SecondViralCoef%B = -5.8850E02;
!!$    pc%SecondViralCoef%C = -7.6060E07;
!!$    pc%SecondViralCoef%D = -1.6130E21;
!!$    pc%SecondViralCoef%E = 3.9940E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 309.23, -7.8395E00; 1500.00, 6.5078E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.6468E01;
!!$    pc%liqViscosity%B = 1.5335E03;
!!$    pc%liqViscosity%C = 7.5110E-01;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 243.51, 2.3742E-03; 448.15, 2.1165E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.6400E-08;
!!$    pc%vapViscosity%B = 9.4870E-01;
!!$    pc%vapViscosity%C = 7.1000E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 243.51, 3.7547E-06; 1000.00, 1.7295E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0630E-01;
!!$    pc%liqThermalConductivity%B = -2.5000E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 243.51, 1.4542E-01; 447.30, 9.4475E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -6.6840E02;
!!$    pc%vapThermalConductivity%B = 9.3230E-01;
!!$    pc%vapThermalConductivity%C = -4.0710E09;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 447.30, 2.1731E-02; 1000.00, 1.0286E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.5435E-02;
!!$    pc%SurfaceTension%B = 1.3095E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 243.51, 2.8756E-02; 617.70, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_ndecane
!!$
!!$
!!$  subroutine pc_ndodecane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'ndodecane','n-dodecane','duodecane','adakane 12'}; 
!!$    pc%CAS = '112-40-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 12;         
!!$    pc%ChemicalFormula%H = 26;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 170.338;   ! kg/kmol    
!!$    pc%Tcrit = 658.00;              ! K    
!!$    pc%pcrit = 1.8200E6;            ! Pa 
!!$    pc%vcrit = 0.71600;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.238;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 263.57;             ! K
!!$    pc%TripletPointT = 263.57;            ! K
!!$    pc%TripletPointp = 6.1520E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 489.47;       ! K
!!$    pc%liqMolarVolume  = 0.228605;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.9072E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 3.3180E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 6.2415E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 3.6840E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -7.5137E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.5764;            ! -
!!$    pc%RadiusOfGyration = 5.9140E-10;      ! m
!!$    pc%SolubilityParameter = 1.5900E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.12960;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.7740E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4151;           ! -
!!$    pc%FlashPoint = 347.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.6 4.7 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 476.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 5.2800E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 263.57, 5.2800E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.5541E-01;
!!$    pc%liqDensity%B = 2.5511E-01;
!!$    pc%liqDensity%C = 6.5800E02;
!!$    pc%liqDensity%D = 2.9368E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 263.57, 4.5132E00; 658.00, 1.3932E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.3747E02;
!!$    pc%vapPressure%B = -1.1976E04;
!!$    pc%vapPressure%C = -1.6698E01;
!!$    pc%vapPressure%D =  8.0906E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 263.57, 6.1534E-01; 658.00, 1.8223E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 7.7337E07;
!!$    pc%HeatOfVap%B = 4.0681E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 263.57, 6.2802E07; 658.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -4.0970E04;
!!$    pc%solHeatCapacity%B = 2.8390E03;
!!$    pc%solHeatCapacity%C = -1.3232E01;
!!$    pc%solHeatCapacity%D = 2.7972E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 29.97, 3.2983E04; 258.41, 2.9175E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 5.0821E05;
!!$    pc%liqHeatCapacity%B = -1.3687E03;
!!$    pc%liqHeatCapacity%C = 3.1015E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 263.57, 3.6292E05; 330.00, 3.9429E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.1295E05;
!!$    pc%IG_HeatCapacity%B = 6.6330E05;
!!$    pc%IG_HeatCapacity%C = 1.7155E03;
!!$    pc%IG_HeatCapacity%D = 4.5161E05;
!!$    pc%IG_HeatCapacity%E = 7.7750E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 2.2442E05; 1500.00, 7.4325E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 8.800E-01;
!!$    pc%SecondViralCoef%B = -1.0910E03;
!!$    pc%SecondViralCoef%C = -5.0300E07;
!!$    pc%SecondViralCoef%D = -5.4871E21;
!!$    pc%SecondViralCoef%E = 1.4959E24;
!!$    pc%SecondViralCoef%range(:,1) = (/ 329.00, -1.0701E01; 1500.00, 1.3759E-01 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.0607E01;
!!$    pc%liqViscosity%B = 1.9430E03;
!!$    pc%liqViscosity%C = 1.3205E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 263.57, 2.8109E-03; 489.47, 2.1196E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 6.3440E-08;
!!$    pc%vapViscosity%B = 8.2870E-01;
!!$    pc%vapViscosity%C = 2.1950E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 263.57, 3.5111E-06; 1000.00, 1.5932E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0470E-01;
!!$    pc%liqThermalConductivity%B = -2.3260E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 263.57, 1.4339E-01; 489.47, 9.0849E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 5.7190E-06;
!!$    pc%vapThermalConductivity%B = 1.4699E00;
!!$    pc%vapThermalConductivity%C = 5.7940E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 489.47, 2.3537E-02; 1000.00, 9.3010E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.5493E-02;
!!$    pc%SurfaceTension%B = 1.3262E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 263.57, 2.8150E-02; 658.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_ndodecane
!!$
!!$
!!$  subroutine pc_nheptane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'nheptane','n-heptane','dipropymethane'}; 
!!$    pc%CAS = '142-82-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 7;         
!!$    pc%ChemicalFormula%H = 16;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 100.204;   ! kg/kmol    
!!$    pc%Tcrit = 540.20;              ! K    
!!$    pc%pcrit = 2.7400E6;            ! Pa 
!!$    pc%vcrit = 0.42800;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.261;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 182.57;             ! K
!!$    pc%TripletPointT = 182.57;            ! K
!!$    pc%TripletPointp = 1.8269E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 371.58;       ! K
!!$    pc%liqMolarVolume  = 0.147024;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.8765E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 8.1650E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 4.2798E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.4050E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -4.4647E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3495;            ! -
!!$    pc%RadiusOfGyration = 4.1730E-10;      ! m
!!$    pc%SolubilityParameter = 1.5200E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.07849;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.0990E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3851;           ! -
!!$    pc%FlashPoint = 269.00;                ! K
!!$    pc%FlammabilityLimits = (/ 1.0 7.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 477.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 9.1435E00;
!!$    pc%solDensity%B = -2.8756E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 133.15, 8.7606E00; 182.57, 8.6185E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 6.1259E-01;
!!$    pc%liqDensity%B = 2.6211E-01;
!!$    pc%liqDensity%C = 5.4020E02;
!!$    pc%liqDensity%D = 2.8141E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 182.57, 7.6998E00; 540.20, 2.3371E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  8.7829E01;
!!$    pc%vapPressure%B = -6.9964E03;
!!$    pc%vapPressure%C = -9.8802E00;
!!$    pc%vapPressure%D =  7.2099E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 182.57, 1.8269E-01; 540.20, 2.7192E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.0014E07;
!!$    pc%HeatOfVap%B = 3.8795E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 182.57, 4.2619E07; 540.20, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -2.3900E04;
!!$    pc%solHeatCapacity%B = 1.9408E03;
!!$    pc%solHeatCapacity%C = -1.0480E01;
!!$    pc%solHeatCapacity%D = 2.7190E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.0942E04; 182.58, 1.4658E+05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 6.1260E01;
!!$    pc%liqHeatCapacity%B = 3.1441E05;
!!$    pc%liqHeatCapacity%C = 1.8246E03;
!!$    pc%liqHeatCapacity%D = -2.5479E03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 182.57, 1.9989E05; 520.00, 4.0657E05 /);
!!$    pc%liqHeatCapacity%eqn = 114;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.2015E05;
!!$    pc%IG_HeatCapacity%B = 4.0010E05;
!!$    pc%IG_HeatCapacity%C = 1.6766E03;
!!$    pc%IG_HeatCapacity%D = 2.7400E05;
!!$    pc%IG_HeatCapacity%E = 7.5640E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.2828E05; 1500.00, 4.4283E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.7460E-01;
!!$    pc%SecondViralCoef%B = -2.9100E02;
!!$    pc%SecondViralCoef%C = -4.4180E07;
!!$    pc%SecondViralCoef%D = -8.8000E19;
!!$    pc%SecondViralCoef%E = 1.2850E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 270.10, - 4.4720E00; 1500.00, 6.7507E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.4451E01;
!!$    pc%liqViscosity%B = 1.5331E03;
!!$    pc%liqViscosity%C = 2.0087E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 182.57, 3.7197E-03; 373.15, 2.1455E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 6.6720E-08;
!!$    pc%vapViscosity%B = 8.2837E-01;
!!$    pc%vapViscosity%C = 8.5752E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 182.57, 3.3910E-06; 1000.00, 1.8777E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1500E-01;
!!$    pc%liqThermalConductivity%B = -3.0300E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 182.57, 1.5968E-01; 371.58, 1.0241E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -7.0028E-02;
!!$    pc%vapThermalConductivity%B = 3.8068E-01;
!!$    pc%vapThermalConductivity%C = -7.0499E03;
!!$    pc%vapThermalConductivity%D = -2.4005E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 339.15, 1.5827E-02; 1000.00, 1.1493E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.4143E-02;
!!$    pc%SurfaceTension%B = 1.2512E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 182.57, 3.2317E-02; 540.20, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_nheptane
!!$
!!$
!!$  subroutine pc_nhexane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'nhexane','n-hexane','hexyl hydride'}; 
!!$    pc%CAS = '110-54-3';
!!$    ! components
!!$    pc%ChemicalFormula%C = 6;         
!!$    pc%ChemicalFormula%H = 14;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 86.177;    ! kg/kmol    
!!$    pc%Tcrit = 507.60;              ! K    
!!$    pc%pcrit = 3.0250E6;            ! Pa 
!!$    pc%vcrit = 0.37100;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.266;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 177.83;             ! K
!!$    pc%TripletPointT = 177.83;            ! K
!!$    pc%TripletPointp = 9.0170E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 341.88;       ! K
!!$    pc%liqMolarVolume  = 0.131362;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.6694E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -6.6340E04;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 3.8874E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.3080E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.8551E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3013;            ! -
!!$    pc%RadiusOfGyration = 3.7690E-10;      ! m
!!$    pc%SolubilityParameter = 1.4900E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.06826;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 9.6400E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3723;           ! -
!!$    pc%FlashPoint = 251.50;                ! K
!!$    pc%FlammabilityLimits = (/ 1.1 7.7 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 498.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.0657E01;
!!$    pc%solDensity%B = -4.6530E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 133.15, 1.0037E01; 177.83, 9.8296E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 7.0824E-01;
!!$    pc%liqDensity%B = 2.6411E-01;
!!$    pc%liqDensity%C = 5.0760E02;
!!$    pc%liqDensity%D = 2.7537E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 177.83, 8.7470E00; 507.60, 2.6816E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.0465E02;
!!$    pc%vapPressure%B = -6.9955E03;
!!$    pc%vapPressure%C = -1.2702E01;
!!$    pc%vapPressure%D =  1.2381E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 177.83, 9.0169E-01; 507.60, 3.0449E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.4544E07;
!!$    pc%HeatOfVap%B = 3.9002E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 177.83, 3.7647E07; 507.60, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.9370E04;
!!$    pc%solHeatCapacity%B = 1.6560E03;
!!$    pc%solHeatCapacity%C = -8.3720E00;
!!$    pc%solHeatCapacity%D = 2.0200E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.0563E04; 177.86, 1.2398E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.7212E05;
!!$    pc%liqHeatCapacity%B = -1.8378E02;
!!$    pc%liqHeatCapacity%C = 8.8734E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 177.83, 1.6750E05; 460.00, 2.7534E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.0440E05;
!!$    pc%IG_HeatCapacity%B = 3.5230E05;
!!$    pc%IG_HeatCapacity%C = 1.6946E03;
!!$    pc%IG_HeatCapacity%D = 2.3690E05;
!!$    pc%IG_HeatCapacity%E = 7.6160E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.1117E05; 1500.00, 3.8620E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.9700E-01;
!!$    pc%SecondViralCoef%B = -2.0000E02;
!!$    pc%SecondViralCoef%C = -3.7100E07;
!!$    pc%SecondViralCoef%D = 1.6700E19;
!!$    pc%SecondViralCoef%E = -6.1800E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 253.70, -3.3097E00; 1500.00, 5.2675E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.0715E01;
!!$    pc%liqViscosity%B = 1.2075E03;
!!$    pc%liqViscosity%C = 1.4993E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 177.83, 2.1182E-03; 343.15, 2.1541E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.7514E-07;
!!$    pc%vapViscosity%B = 7.0737E-01;
!!$    pc%vapViscosity%C = 1.5714E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 177.83, 3.6305E-06; 1000.00, 2.0050E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.2492E-01;
!!$    pc%liqThermalConductivity%B = -3.5300E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 177.83, 1.6209E-01; 370.00, 9.4199E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -6.5050E02;
!!$    pc%vapThermalConductivity%B = 8.0530E-01;
!!$    pc%vapThermalConductivity%C = -1.4121E09;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 339.09, 1.7036E-02; 1000.00, 1.2003E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.5003E-02;
!!$    pc%SurfaceTension%B = 1.2674E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 177.83, 3.1841E-02; 507.60, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_nhexane
!!$
!!$
!!$  subroutine pc_nnonane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'nnonane','n-nonane','nonyl hydride'}; 
!!$    pc%CAS = '111-84-2';
!!$    ! components
!!$    pc%ChemicalFormula%C = 9;         
!!$    pc%ChemicalFormula%H = 20;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 128.258;   ! kg/kmol    
!!$    pc%Tcrit = 594.60;              ! K    
!!$    pc%pcrit = 2.2900E6;            ! Pa 
!!$    pc%vcrit = 0.54400;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.252;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 219.66;             ! K
!!$    pc%TripletPointT = 219.66;            ! K
!!$    pc%TripletPointp = 4.3058E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 423.97;       ! K
!!$    pc%liqMolarVolume  = 0.179559;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.2874E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 2.4980E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 5.0640E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 1.5470E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -5.6846E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.4435;            ! -
!!$    pc%RadiusOfGyration = 4.8140E-10;      ! m
!!$    pc%SolubilityParameter = 1.5600E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.09895;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.3690E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4031;           ! -
!!$    pc%FlashPoint = 304.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.7 5.6 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 478.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 7.4283E00;
!!$    pc%solDensity%B = -2.3699E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 193.15, 6.9706E00; 219.66, 6.9077E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 4.8387E-01;
!!$    pc%liqDensity%B = 2.6147E-01;
!!$    pc%liqDensity%C = 5.9460E02;
!!$    pc%liqDensity%D = 2.8281E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 219.66, 6.0070E00; 594.60, 1.8506E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.0935E02;
!!$    pc%vapPressure%B = -9.0304E03;
!!$    pc%vapPressure%C = -1.2882E01;
!!$    pc%vapPressure%D =  7.8544E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 219.66, 4.3058E-01; 594.60, 2.3054E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.0370E07;
!!$    pc%HeatOfVap%B = 3.8522E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 219.66, 5.0545E07; 594.60, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -3.1900E04;
!!$    pc%solHeatCapacity%B = 2.3720E03;
!!$    pc%solHeatCapacity%C = -1.2440E01;
!!$    pc%solHeatCapacity%D = 3.0180E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.0805E04; 217.21, 2.0569E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.8308E05;
!!$    pc%liqHeatCapacity%B = -1.1398E03;
!!$    pc%liqHeatCapacity%C = 2.7101E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 219.66, 2.6346E05; 325.00, 2.9890E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.5175E05;
!!$    pc%IG_HeatCapacity%B = 4.9150E05;
!!$    pc%IG_HeatCapacity%C = 1.6448E03;
!!$    pc%IG_HeatCapacity%D = 3.4700E05;
!!$    pc%IG_HeatCapacity%E = 7.4960E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.6257E05; 1500.00, 5.5407E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 3.9060E-01;
!!$    pc%SecondViralCoef%B = -4.6860E02;
!!$    pc%SecondViralCoef%C = -6.6800E07;
!!$    pc%SecondViralCoef%D = -7.8930E20;
!!$    pc%SecondViralCoef%E = 1.8200E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 297.83, -6.5900E00; 1500.00, 5.8381E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.1149E01;
!!$    pc%liqViscosity%B = 1.6580E03;
!!$    pc%liqViscosity%C = 1.4540E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 219.66, 3.1483E-03; 423.97, 2.1559E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.0344E-07;
!!$    pc%vapViscosity%B = 7.7301E-01;
!!$    pc%vapViscosity%C = 2.2047E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 219.66, 3.3347E-06; 1000.00, 1.7668E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0900E-01;
!!$    pc%liqThermalConductivity%B = -2.6400E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 219.66, 1.5101E-01; 423.97, 9.7072E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -6.5771E-02;
!!$    pc%vapThermalConductivity%B = 2.7198E-01;
!!$    pc%vapThermalConductivity%C = -3.4823E03;
!!$    pc%vapThermalConductivity%D = -1.5803E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 423.97, 2.1299E-02; 1000.00, 1.0597E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.4975E-02;
!!$    pc%SurfaceTension%B = 1.2897E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 219.66, 3.0331E-02; 594.60, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_nnonane
!!$
!!$
!!$  subroutine pc_noctane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'noctane','n-octane','octyl hydride'}; 
!!$    pc%CAS = '111-65-9';
!!$    ! components
!!$    pc%ChemicalFormula%C = 8;         
!!$    pc%ChemicalFormula%H = 18;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 114.231;   ! kg/kmol    
!!$    pc%Tcrit = 568.70;              ! K    
!!$    pc%pcrit = 2.4900E6;            ! Pa 
!!$    pc%vcrit = 0.48600;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.256;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 216.38;             ! K
!!$    pc%TripletPointT = 216.38;            ! K
!!$    pc%TripletPointp = 2.1083E00;         ! Pa
!!$    pc%NormalBoilingPoint = 398.83;       ! K
!!$    pc%liqMolarVolume  = 0.163374;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -2.0875E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 1.6000E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 4.6723E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 2.0740E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -5.0741E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.3996;            ! -
!!$    pc%RadiusOfGyration = 4.5460E-10;      ! m
!!$    pc%SolubilityParameter = 1.5400E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.08872;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.3690E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3951;           ! -
!!$    pc%FlashPoint = 286.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.8 6.5 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 479.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 8.3409E00;
!!$    pc%solDensity%B = -3.1515E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 133.15, 7.9213E00; 216.38, 7.6590E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 5.3731E-01;
!!$    pc%liqDensity%B = 2.6115E-01;
!!$    pc%liqDensity%C = 5.6870E02;
!!$    pc%liqDensity%D = 2.8034E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 216.38, 6.6558E00; 568.70, 2.0575E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  9.6084E01;
!!$    pc%vapPressure%B = -7.9002E03;
!!$    pc%vapPressure%C = -1.1003E01;
!!$    pc%vapPressure%D =  7.1802E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 216.38, 2.1083E00; 568.70, 2.4673E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 5.5180E07;
!!$    pc%HeatOfVap%B = 3.8467E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 216.38, 4.5898E07; 568.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -2.4000E04;
!!$    pc%solHeatCapacity%B = 1.9472E03;
!!$    pc%solHeatCapacity%C = -8.5360E00;
!!$    pc%solHeatCapacity%D = 1.3400E-02;
!!$    pc%solHeatCapacity%E = 2.0940E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.1640E04; 216.37, 1.7933E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 2.2483E05;
!!$    pc%liqHeatCapacity%B = -1.8663E02;
!!$    pc%liqHeatCapacity%C = 9.5891E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 216.38, 2.2934E05; 460.00, 3.4189E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 1.3554E05;
!!$    pc%IG_HeatCapacity%B = 4.4310E05;
!!$    pc%IG_HeatCapacity%C = 1.6356E03;
!!$    pc%IG_HeatCapacity%D = 3.0540E05;
!!$    pc%IG_HeatCapacity%E = 7.4640E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 1.4529E05; 1500.00, 4.9764E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 2.7390E-01;
!!$    pc%SecondViralCoef%B = -3.2144E02;
!!$    pc%SecondViralCoef%C = -6.6830E07;
!!$    pc%SecondViralCoef%D = -1.2710E20;
!!$    pc%SecondViralCoef%E = 1.6103E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 284.38, -5.4098E00; 1500.00, 3.9801E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.0463E01;
!!$    pc%liqViscosity%B = 1.4974E03;
!!$    pc%liqViscosity%C = 1.3790E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 216.38, 2.1813E-03; 398.83, 2.1383E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.1191E-08;
!!$    pc%vapViscosity%B = 9.2925E-01;
!!$    pc%vapViscosity%C = 5.5092E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 216.38, 3.6772E-06; 1000.00, 1.8134E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1560E-01;
!!$    pc%liqThermalConductivity%B = -2.9483E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 216.38, 1.5180E-01; 398.83, 9.8013E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -8.7580E03;
!!$    pc%vapThermalConductivity%B = 8.4480E-01;
!!$    pc%vapThermalConductivity%C = -2.7121E10;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 339.00, 1.5025E-02; 1000.00, 1.1053E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.2789E-02;
!!$    pc%SurfaceTension%B = 1.2323E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 216.38, 2.9261E-02; 568.70, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_noctane
!!$
!!$
!!$  subroutine pc_npentane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'npentane','n-pentane','pental hydride'}; 
!!$    pc%CAS = '109-66-0';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 12;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 72.150;    ! kg/kmol    
!!$    pc%Tcrit = 469.10;              ! K    
!!$    pc%pcrit = 3.3700E6;            ! Pa 
!!$    pc%vcrit = 0.31300;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.270;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 143.42;             ! K
!!$    pc%TripletPointT = 143.42;            ! K
!!$    pc%TripletPointp = 6.8642E-02;        ! Pa
!!$    pc%NormalBoilingPoint = 309.22;       ! K
!!$    pc%liqMolarVolume  = 0.116045;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.4676E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 8.8130E06;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 3.4945E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 8.4015E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -3.2449E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2515;            ! -
!!$    pc%RadiusOfGyration = 3.3370E-10;      ! m
!!$    pc%SolubilityParameter = 1.4400E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.05803;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 8.2900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.3547;           ! -
!!$    pc%FlashPoint = 233.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.3 8.0 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 516.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.3533E01;
!!$    pc%solDensity%B = -6.0608E-03;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 23.15, 1.3393E01; 143.42, 1.2664E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 8.4947E-01;
!!$    pc%liqDensity%B = 2.6726E-01;
!!$    pc%liqDensity%C = 4.6970E02;
!!$    pc%liqDensity%D = 2.7789E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 143.42, 1.0474E01; 469.70, 3.1784E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  7.8741E01;
!!$    pc%vapPressure%B = -5.4203E03;
!!$    pc%vapPressure%C = -8.8253E00;
!!$    pc%vapPressure%D =  9.6171E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 143.42, 6.8642E-02; 469.70, 3.3642E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 3.9109E07;
!!$    pc%HeatOfVap%B = 3.8681E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 143.42, 3.3968E07; 469.70, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.0050E04;
!!$    pc%solHeatCapacity%B = 1.0103E03;
!!$    pc%solHeatCapacity%C = 4.2700E00;
!!$    pc%solHeatCapacity%D = -1.0169E-01;
!!$    pc%solHeatCapacity%E = 4.2200E-04;
!!$    pc%solHeatCapacity%range(:,1) = (/ 12.29, 2.8324E03; 134.60, 9.3832E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.5908E05;
!!$    pc%liqHeatCapacity%B = -2.7050E02;
!!$    pc%liqHeatCapacity%C = 9.9537E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 143.42, 1.4076E05; 390.00, 2.0498E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 8.8050E04;
!!$    pc%IG_HeatCapacity%B = 3.0110E05;
!!$    pc%IG_HeatCapacity%C = 1.6502E03;
!!$    pc%IG_HeatCapacity%D = 1.8920E05;
!!$    pc%IG_HeatCapacity%E = 7.4760E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 9.4039E04; 1500.00, 3.2927E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.9160E-01;
!!$    pc%SecondViralCoef%B = -1.9800E02;
!!$    pc%SecondViralCoef%C = -1.6400E07;
!!$    pc%SecondViralCoef%D = -1.7920E19;
!!$    pc%SecondViralCoef%E = 2.6850E21;
!!$    pc%SecondViralCoef%range(:,1) = (/ 234.80,-2.6205E00; 1500.00, 5.4740E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.0383E01;
!!$    pc%liqViscosity%B = 1.0504E03;
!!$    pc%liqViscosity%C = 1.4874E00;
!!$    pc%liqViscosity%D = -2.0167E-27;
!!$    pc%liqViscosity%E = 1E01;
!!$    pc%liqViscosity%range(:,1) = (/ 143.42, 3.4377E-03; 465.15, 4.7971E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 6.3412E-08;
!!$    pc%vapViscosity%B = 8.4758E-01;
!!$    pc%vapViscosity%C = 4.1718E01;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 143.42, 3.3051E-06; 1000.00, 2.1240E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.5370E-01;
!!$    pc%liqThermalConductivity%B = -5.7600E-04;
!!$    pc%liqThermalConductivity%C = 3.4400E-07;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 143.42, 1.7817E-01; 445.00, 6.5501E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -6.8440E02;
!!$    pc%vapThermalConductivity%B = 7.6400E-01;
!!$    pc%vapThermalConductivity%C = -1.0550E09;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 273.15, 1.2879E-02; 1000.00, 1.2707E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.2020E-02;
!!$    pc%SurfaceTension%B = 1.2041E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 143.42, 3.3546E-02; 469.70, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_npentane
!!$
!!$
!!$  subroutine pc_ntetradecane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'ntetradecane','n-tetradecane','tetradecyl hydride'}; 
!!$    pc%CAS = '629-59-4';
!!$    ! components
!!$    pc%ChemicalFormula%C = 14;         
!!$    pc%ChemicalFormula%H = 30;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 198.392;  ! kg/kmol    
!!$    pc%Tcrit = 693.00;             ! K    
!!$    pc%pcrit = 1.5700E06;          ! Pa 
!!$    pc%vcrit = 0.83000;            ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.226; ! -
!!$
!!$
!!$    pc%MeltingPoint = 279.01;              ! K
!!$    pc%TripletPointT = 279.01;             ! K
!!$    pc%TripletPointp = 2.5269E-01;         ! Pa
!!$    pc%NormalBoilingPoint = 526.73;        ! K
!!$    pc%liqMolarVolume = 0.261271;          ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.3244E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 6.5990E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 7.0259E05;             ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 4.5070E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -8.7328E09;   ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.6430;            ! -
!!$    pc%RadiusOfGyration = 6.4270E-10;      ! m
!!$    pc%SolubilityParameter = 1.6100E04;    ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.15010;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 2.0440E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4269;           ! -
!!$    pc%FlashPoint = 373.15;                ! K
!!$    pc%FlammabilityLimits = (/ 0.5 4.7 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 473.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 4.6071E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 0, 4.6071E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.0545E-01;
!!$    pc%liqDensity%B = 2.5350E-01;
!!$    pc%liqDensity%C = 6.9300E02;
!!$    pc%liqDensity%D = 3.0538E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 279.01, 3.8924E00; 693.00, 1.2049E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.4047E02;
!!$    pc%vapPressure%B = -1.3231E04;
!!$    pc%vapPressure%C = -1.6859E01;
!!$    pc%vapPressure%D =  6.5877E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 279.01, 2.5268E-01; 693.00, 1.5693E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 9.0539E07;
!!$    pc%HeatOfVap%B = 4.4467E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 279.01, 7.2002E07; 693.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -4.9782E04;
!!$    pc%solHeatCapacity%B = 3.3010E03;
!!$    pc%solHeatCapacity%C = -1.5807E01;
!!$    pc%solHeatCapacity%D = 3.6919E-02;
!!$    pc%solHeatCapacity%E = -1.3660E-05;
!!$    pc%solHeatCapacity%range(:,1) = (/ 29.98, 3.5958E04; 279.01, 3.5981E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.5314E05;
!!$    pc%liqHeatCapacity%B = 2.9130E01;
!!$    pc%liqHeatCapacity%C = 8.6116E-01;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 279.01, 4.2831E05; 526.73, 6.0741E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.3082E05;
!!$    pc%IG_HeatCapacity%B = 7.8678E05;
!!$    pc%IG_HeatCapacity%C = 1.6823E03;
!!$    pc%IG_HeatCapacity%D = 5.4486E05;
!!$    pc%IG_HeatCapacity%E = 7.4310E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 2.4864E05; 1500.00, 8.6225E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.8209E-01;
!!$    pc%SecondViralCoef%B = -5.0132E02;
!!$    pc%SecondViralCoef%C = -2.8594E08;
!!$    pc%SecondViralCoef%D = -2.5720E21;
!!$    pc%SecondViralCoef%E = 4.5523E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 346.50, -1.3893E01; 1500.00, 6.3065E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.0486E01;
!!$    pc%liqViscosity%B = 2.0884E03;
!!$    pc%liqViscosity%C = 1.2852E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 279.01, 3.1395E-03; 528.15, 2.0874E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 4.4565E-08;
!!$    pc%vapViscosity%B = 8.6840E-01;
!!$    pc%vapViscosity%C = 2.2816E02;
!!$    pc%vapViscosity%D = -4.3472E03;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 279.01, 3.3635E-06; 1000.00, 1.4672E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.02936E-01;
!!$    pc%liqThermalConductivity%B = -2.1798E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 279.01, 1.4211E-01; 526.73, 8.8113E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -1.6362E02;
!!$    pc%vapThermalConductivity%B = 9.1930E-01;
!!$    pc%vapThermalConductivity%C = -1.0876E09;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 526.73, 2.5172E-02; 1000.0, 8.6152E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.6276E-02;
!!$    pc%SurfaceTension%B = 1.3631E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 279.01, 2.7883E-02; 693.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_ntetradecane
!!$
!!$
!!$  subroutine pc_ntridecane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'ntridecane','n-tridecane','tridecyl hydride'}; 
!!$    pc%CAS = '629-50-5';
!!$    ! components
!!$    pc%ChemicalFormula%C = 13;         
!!$    pc%ChemicalFormula%H = 28;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 184.365;   ! kg/kmol    
!!$    pc%Tcrit = 675.00;              ! K    
!!$    pc%pcrit = 1.6800E6;            ! Pa 
!!$    pc%vcrit = 0.77500;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.232;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 267.76;             ! K
!!$    pc%TripletPointT = 267.76;            ! K
!!$    pc%TripletPointp = 2.5096E-02;        ! Pa
!!$    pc%NormalBoilingPoint = 508.62;       ! K
!!$    pc%liqMolarVolume  = 0.244631;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.1177E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = 5.7710E07;            ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 6.6337E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 2.8500E07;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -8.1229E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.6174;            ! -
!!$    pc%RadiusOfGyration = 6.1960E-10;      ! m
!!$    pc%SolubilityParameter = 1.6100E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.13960;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 1.9090E09;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4235;           ! -
!!$    pc%FlashPoint = 352.00;                ! K
!!$    pc%FlammabilityLimits = (/ 0.6 4.7 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 475.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 5.0311E00;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 0, 5.0311E00 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 3.2160E-01;
!!$    pc%liqDensity%B = 2.5040E-01;
!!$    pc%liqDensity%C = 6.7500E02;
!!$    pc%liqDensity%D = 3.0710E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 267.76, 4.2035E00; 675.00, 1.2843E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.3745E02;
!!$    pc%vapPressure%B = -1.2549E04;
!!$    pc%vapPressure%C = -1.6543E01;
!!$    pc%vapPressure%D =  7.1275E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 267.76, 2.5096E-01; 675.00, 1.6786E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 8.4339E07;
!!$    pc%HeatOfVap%B = 4.2570E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 267.76, 6.8015E07; 675.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -3.0600E04;
!!$    pc%solHeatCapacity%B = 2.3520E03;
!!$    pc%solHeatCapacity%C = -3.2440E00;
!!$    pc%solHeatCapacity%D = -3.5173E-02;
!!$    pc%solHeatCapacity%E = 1.3218E-04;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 1.4882E04; 255.00, 3.3389E05 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 3.5018E05;
!!$    pc%liqHeatCapacity%B = -1.0470E02;
!!$    pc%liqHeatCapacity%C = 1.0022E00;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 267.76, 3.9400E05; 508.62, 5.5619E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.1496E05;
!!$    pc%IG_HeatCapacity%B = 7.3045E05;
!!$    pc%IG_HeatCapacity%C = 1.6695E03;
!!$    pc%IG_HeatCapacity%D = 4.9998E05;
!!$    pc%IG_HeatCapacity%E = 7.4102E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 2.3156E05; 1500.00, 8.0251E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 4.4720E-01;
!!$    pc%SecondViralCoef%B = -4.7897E02;
!!$    pc%SecondViralCoef%C = -2.3434E08;
!!$    pc%SecondViralCoef%D = -1.7856E21;
!!$    pc%SecondViralCoef%E = 3.0980E23;
!!$    pc%SecondViralCoef%range(:,1) = (/ 337.50, -1.2222E01; 1500.00, 5.8391E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -2.1008E01;
!!$    pc%liqViscosity%B = 2.0425E03;
!!$    pc%liqViscosity%C = 1.3692E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 267.76, 3.2603E-03; 508.62, 2.1182E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.5585E-08;
!!$    pc%vapViscosity%B = 8.9870E-01;
!!$    pc%vapViscosity%C = 1.6530E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 267.76, 3.3441E-06; 1000.00, 1.5168E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.0447E-01;
!!$    pc%liqThermalConductivity%B = -2.2612E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 267.76, 1.4392E-01; 508.62, 8.9461E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 5.3701E-06;
!!$    pc%vapThermalConductivity%B = 1.4751E00;
!!$    pc%vapThermalConductivity%C = 5.9909E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 508.62, 2.4219E-02; 1000.00, 8.9415E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.5322E-02;
!!$    pc%SurfaceTension%B = 1.3257E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 267.76, 2.8312E-02; 675.00, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_ntridecane
!!$
!!$
!!$  subroutine pc_propane(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'propane','dimethylmethane','propyl hydride'}; 
!!$    pc%CAS = '74-98-6';
!!$    ! components
!!$    pc%ChemicalFormula%C = 3;         
!!$    pc%ChemicalFormula%H = 8;
!!$    pc%ChemicalFormula%O = 0;
!!$
!!$    pc%MolecularWeight = 44.096;    ! kg/kmol    
!!$    pc%Tcrit = 369.83;              ! K    
!!$    pc%pcrit = 4.2480E6;            ! Pa 
!!$    pc%vcrit = 0.20000;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.276;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 85.47;              ! K
!!$    pc%TripletPointT = 85.47;             ! K
!!$    pc%TripletPointp = 1.6850E-04;        ! Pa
!!$    pc%NormalBoilingPoint = 231.11;       ! K
!!$    pc%liqMolarVolume  = 0.075700;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.0468E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -2.4390E07;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 2.7020E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 3.5240E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.0431E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.1523;            ! -
!!$    pc%RadiusOfGyration = 2.4310E-10;      ! m
!!$    pc%SolubilityParameter = 1.3100E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 0;                   ! C^m
!!$    pc%VanDerWaalsVolume = 0.03757;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 5.5900E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.2861;           ! -
!!$    pc%FlashPoint = MAX_REAL_WP;                   ! K
!!$    pc%FlammabilityLimits = (/ 2.0 9.5 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 723.00;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.8861E01;
!!$    pc%solDensity%B = -2.0332E-02;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 23.15, 1.8390E01; 85.47, 1.7123E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.3757E00;
!!$    pc%liqDensity%B = 2.7453E-01;
!!$    pc%liqDensity%C = 3.6983E02;
!!$    pc%liqDensity%D = 2.9359E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 85.47, 1.6583E01; 369.83, 5.0111E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  5.9078E01;
!!$    pc%vapPressure%B = -3.4926E03;
!!$    pc%vapPressure%C = -6.0669E00;
!!$    pc%vapPressure%D =  1.0919E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 85.47, 1.6788E-04; 369.83, 4.2135E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 2.9209E07;
!!$    pc%HeatOfVap%B = 7.8237E-01;
!!$    pc%HeatOfVap%C = -7.7319E-01;
!!$    pc%HeatOfVap%D = 3.9246E-01;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 85.47, 2.4787E07; 369.83, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.1230E04;
!!$    pc%solHeatCapacity%B = 1.0590E03;
!!$    pc%solHeatCapacity%C = -3.6000E00;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 30.00, 1.7300E04; 84.00, 5.2324E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 6.2983E01;
!!$    pc%liqHeatCapacity%B = 1.1363E05;
!!$    pc%liqHeatCapacity%C = 6.3321E02;
!!$    pc%liqHeatCapacity%D = -8.7346E02;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 85.47, 8.4879E04; 360.00, 2.6079E05 /);
!!$    pc%liqHeatCapacity%eqn = 114;    
!!$
!!$    pc%IG_HeatCapacity%A = 5.1920E04;
!!$    pc%IG_HeatCapacity%B = 1.9245E05;
!!$    pc%IG_HeatCapacity%C = 1.6265E03;
!!$    pc%IG_HeatCapacity%D = 1.1680E05;
!!$    pc%IG_HeatCapacity%E = 7.2360E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 200.00, 5.6322E04; 1500.00, 2.0556E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.1270E-01;
!!$    pc%SecondViralCoef%B = -9.9200E01;
!!$    pc%SecondViralCoef%C = -4.5100E06;
!!$    pc%SecondViralCoef%D = 3.0900E17;
!!$    pc%SecondViralCoef%E = -7.0500E19;
!!$    pc%SecondViralCoef%range(:,1) = (/ 184.91 ,-1.1900E00; 1500.00, 4.5230E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -6.9281E00;
!!$    pc%liqViscosity%B = 4.2076E02;
!!$    pc%liqViscosity%C = -6.3276E-01;
!!$    pc%liqViscosity%D = -1.7130E-26;
!!$    pc%liqViscosity%E = 1E01;
!!$    pc%liqViscosity%range(:,1) = (/ 85.47, 8.0682E-03; 360.00, 4.0666E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 2.4993E-07;
!!$    pc%vapViscosity%B = 6.8612E-01;
!!$    pc%vapViscosity%C = 1.7934E02;
!!$    pc%vapViscosity%D = -8.2546E03;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 85.47, 2.6865E-06; 1000.00, 2.4411E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.6755E-01;
!!$    pc%liqThermalConductivity%B = -6.6457E-04;
!!$    pc%liqThermalConductivity%C = 2.7740E-07;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 85.47, 2.1278E-01; 350.00, 6.8932E-02 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = -1.1200E00;
!!$    pc%vapThermalConductivity%B = 1.0972E-01;
!!$    pc%vapThermalConductivity%C = -9.8346E03;
!!$    pc%vapThermalConductivity%D = -7.5358E06;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 231.11, 1.1142E-02; 1000.00, 1.4599E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 5.0920E-02;
!!$    pc%SurfaceTension%B = 1.2197E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 85.47, 3.6956E-02; 369.83, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_propane
!!$
!!$
!!$  subroutine pc_tetrahydrofuran(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'tetrahydrofuran','hydrofuran','tetramethylene oxide'}; 
!!$    pc%CAS = '109-99-9';
!!$    ! components
!!$    pc%ChemicalFormula%C = 4;         
!!$    pc%ChemicalFormula%H = 8;
!!$    pc%ChemicalFormula%O = 1;
!!$
!!$    pc%MolecularWeight = 72.107;    ! kg/kmol    
!!$    pc%Tcrit = 540.15;              ! K    
!!$    pc%pcrit = 5.1678E6;            ! Pa 
!!$    pc%vcrit = 0.22393;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.259;  ! -
!!$
!!$
!!$    pc%MeltingPoint = 164.65;             ! K
!!$    pc%TripletPointT = 164.65;            ! K
!!$    pc%TripletPointp = 1.2753E-01;        ! Pa
!!$    pc%NormalBoilingPoint = 338.00;       ! K
!!$    pc%liqMolarVolume  = 0.081942;        ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -1.8418E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -7.9680E07;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 2.9729E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = 8.5400E06;          ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.3250E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.2264;            ! -
!!$    pc%RadiusOfGyration = 2.6940E-10;      ! m
!!$    pc%SolubilityParameter = 1.9129E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 5.4371E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.04462;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 6.0000E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4050;           ! -
!!$    pc%FlashPoint = 258.11;                ! K
!!$    pc%FlammabilityLimits = (/ 2.0 11.8 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 594.26;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.6400E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 164.65, 1.6400E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 1.2543E00;
!!$    pc%liqDensity%B = 2.8084E-01;
!!$    pc%liqDensity%C = 5.4015E02;
!!$    pc%liqDensity%D = 2.9120E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 164.65, 1.3998E01; 540.15, 4.4662E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  8.5184E01;
!!$    pc%vapPressure%B = -6.3373E03;
!!$    pc%vapPressure%C = -9.5958E00;
!!$    pc%vapPressure%D =  8.1661E-06;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 164.65, 1.2753E-01; 540.15, 5.1720E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 4.4430E07;
!!$    pc%HeatOfVap%B = 3.9080E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 164.65, 3.8545E07; 540.15, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = -1.8200E04;
!!$    pc%solHeatCapacity%B = 1.4470E03;
!!$    pc%solHeatCapacity%C = -9.4700E00;
!!$    pc%solHeatCapacity%D = 2.6620E-02;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ 20.00, 7.1650E03; 164.65, 8.2141E04 /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 1.7165E05;
!!$    pc%liqHeatCapacity%B = -7.9835E02;
!!$    pc%liqHeatCapacity%C = 2.8812E00;
!!$    pc%liqHeatCapacity%D = -2.4813E-03;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 164.65, 1.0723E05; 320.00, 1.2991E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 4.0430E04;
!!$    pc%IG_HeatCapacity%B = 2.5410E05;
!!$    pc%IG_HeatCapacity%C = 1.3866E03;
!!$    pc%IG_HeatCapacity%D = 1.4200E05;
!!$    pc%IG_HeatCapacity%E = -6.3880E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 50.00, 4.0430E04; 1500.00, 2.5462E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.2580E-01;
!!$    pc%SecondViralCoef%B = -1.4320E02;
!!$    pc%SecondViralCoef%C = -2.0790E07;
!!$    pc%SecondViralCoef%D = 6.3140E19;
!!$    pc%SecondViralCoef%E = -1.9480E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 304.58, -1.0912E00; 2701.00, 7.1728E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -1.0335E01;
!!$    pc%liqViscosity%B = 8.8360E02;
!!$    pc%liqViscosity%C = -5.2650E-02;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 164.65, 5.3152E-03; 350.00, 2.9787E-04 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 3.7780E-07;
!!$    pc%vapViscosity%B = 6.5330E-01;
!!$    pc%vapViscosity%C = 2.7101E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 338.00, 9.4123E-06; 998.00, 2.7056E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 1.9428E-01;
!!$    pc%liqThermalConductivity%B = -2.4900E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 164.65, 1.5328E-01; 338.00, 1.1012E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 9.5521E-06;
!!$    pc%vapThermalConductivity%B = 1.4561E00;
!!$    pc%vapThermalConductivity%C = 6.6222E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 338.00, 1.5534E-02; 998.00, 1.3369E-01 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.7130E-02;
!!$    pc%SurfaceTension%B = 1.2310E00;
!!$    pc%SurfaceTension%C = 0;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 164.65, 4.2908E-02; 540.15, 0 /);
!!$    pc%SurfaceTension%eqn = 106;    
!!$
!!$  end subroutine pc_tetrahydrofuran
!!$
!!$
!!$  subroutine pc_tetrahydrofurfurylalcohol(pc)
!!$    implicit none
!!$
!!$    ! ---------------------------------------------
!!$    type(pc_t), pointer, intent(inout) :: pc
!!$    ! ---------------------------------------------
!!$
!!$
!!$
!!$
!!$    ! all alternative names
!!$    pc%name = {'tetrahydrofurfurylalcohol','tetrahydrofuryl carbinol','oxolan-2-methanol'}; 
!!$    pc%CAS = '97-99-4';
!!$    ! components
!!$    pc%ChemicalFormula%C = 5;         
!!$    pc%ChemicalFormula%H = 10;
!!$    pc%ChemicalFormula%O = 2;
!!$
!!$    pc%MolecularWeight = 102.133;   ! kg/kmol    
!!$    pc%Tcrit = 639.00;              ! K    
!!$    pc%pcrit = 4.6600E6;            ! Pa 
!!$    pc%vcrit = 0.29000;             ! m^3/kmol 
!!$    pc%CritCompressFactor = 0.254;  ! -
!!$
!!$
!!$    pc%MeltingPoint = MAX_REAL_WP;             ! K
!!$    pc%TripletPointT = MAX_REAL_WP;            ! K
!!$    pc%TripletPointp = MAX_REAL_WP;            ! Pa
!!$    pc%NormalBoilingPoint = 450.80;    ! K
!!$    pc%liqMolarVolume  = 0.097319;     ! m^3/kmol
!!$
!!$
!!$    pc%IG_HeatOfFormation = -3.6920E08;            ! J/kmol
!!$    pc%IG_GibbsOfFormation = -2.1600E08;           ! J/kmol
!!$    pc%IG_AbsoluteEntropy = 3.7400E5;              ! J/kmolK
!!$    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                ! J/kmol
!!$    pc%StandardNetHeatOfCombustion = -2.7412E9;    ! J/kmol
!!$    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
!!$    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
!!$    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
!!$
!!$
!!$    pc%AcentricFactor = 0.7044;            ! -
!!$    pc%RadiusOfGyration = 3.4810E-10;      ! m
!!$    pc%SolubilityParameter = 2.2920E4;     ! (J/m^3)^.5 
!!$    pc%DipoleMoment = 7.3400E-30;          ! C^m
!!$    pc%VanDerWaalsVolume = 0.05830;        ! m^3/kmol
!!$
!!$
!!$    pc%VanDerWaalsArea = 7.4600E08;        ! m^2/kmol
!!$    pc%RefractiveIndex = 1.4499;           ! -
!!$    pc%FlashPoint = 343.15;                ! K
!!$    pc%FlammabilityLimits = (/ 1.5 9.7 /);     ! vol-!
!!$    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
!!$    pc%AutoignitionTemp = 553.15;          ! K
!!$
!!$
!!$    pc%solDensity%A = 1.2900E01;
!!$    pc%solDensity%B = 0;
!!$    pc%solDensity%C = 0;
!!$    pc%solDensity%D = 0;
!!$    pc%solDensity%E = 0;
!!$    pc%solDensity%range(:,1) = (/ 0.0, 1.2900E01 /);
!!$    pc%solDensity%eqn = 100;
!!$
!!$    pc%liqDensity%A = 9.7001E-01;
!!$    pc%liqDensity%B = 2.8130E-01;
!!$    pc%liqDensity%C = 6.3900E02;
!!$    pc%liqDensity%D = 2.3837E-01;
!!$    pc%liqDensity%E = 0;
!!$    pc%liqDensity%range(:,1) = (/ 193.00, 1.1046E01; 639.00, 3.4483E00 /);
!!$    pc%liqDensity%eqn = 105;
!!$
!!$    pc%vapPressure%A =  1.5059E02;
!!$    pc%vapPressure%B = -1.1574E04;
!!$    pc%vapPressure%C = -1.9025E01;
!!$    pc%vapPressure%D =  1.4141E-05;
!!$    pc%vapPressure%E =  2.0000E00;
!!$    pc%vapPressure%range(:,1) = (/ 193.00, 1.2655E-04; 639.00, 4.6469E06 /);
!!$    pc%vapPressure%eqn = 101;  
!!$
!!$    pc%HeatOfVap%A = 6.4109E07;
!!$    pc%HeatOfVap%B = 2.8538E-01;
!!$    pc%HeatOfVap%C = 0;
!!$    pc%HeatOfVap%D = 0;
!!$    pc%HeatOfVap%E = 0;
!!$    pc%HeatOfVap%range(:,1) = (/ 193.00, 5.7857E07; 639.00, 0 /);
!!$    pc%HeatOfVap%eqn = 106;     
!!$
!!$    pc%solHeatCapacity%A = 0;
!!$    pc%solHeatCapacity%B = 0;
!!$    pc%solHeatCapacity%C = 0;
!!$    pc%solHeatCapacity%D = 0;
!!$    pc%solHeatCapacity%E = 0;
!!$    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
!!$    pc%solHeatCapacity%eqn = 100;
!!$
!!$    pc%liqHeatCapacity%A = 5.2700E04;
!!$    pc%liqHeatCapacity%B = 4.3580E02;
!!$    pc%liqHeatCapacity%C = 0;
!!$    pc%liqHeatCapacity%D = 0;
!!$    pc%liqHeatCapacity%E = 0;
!!$    pc%liqHeatCapacity%range(:,1) = (/ 296.65, 1.8198E05; 469.50, 2.5731E05 /);
!!$    pc%liqHeatCapacity%eqn = 100;    
!!$
!!$    pc%IG_HeatCapacity%A = 2.6823E05;
!!$    pc%IG_HeatCapacity%B = 1.1670E05;
!!$    pc%IG_HeatCapacity%C = 1.6860E03;
!!$    pc%IG_HeatCapacity%D = -3.5200E05;
!!$    pc%IG_HeatCapacity%E = 2.4520E02;
!!$    pc%IG_HeatCapacity%range(:,1) = (/ 300.00, 1.3998E05; 1500.00, 3.3693E05 /);
!!$    pc%IG_HeatCapacity%eqn = 107;
!!$
!!$    pc%SecondViralCoef%A = 1.9766E-01;
!!$    pc%SecondViralCoef%B = -1.8908E02;
!!$    pc%SecondViralCoef%C = -8.1632E07;
!!$    pc%SecondViralCoef%D = -5.5857E20;
!!$    pc%SecondViralCoef%E = 8.4850E22;
!!$    pc%SecondViralCoef%range(:,1) = (/ 319.50, -5.5954E00; 1500.00, 4.7400E-02 /);
!!$    pc%SecondViralCoef%eqn = 104;  
!!$
!!$    pc%liqViscosity%A = -7.9742E00;
!!$    pc%liqViscosity%B = 2.7454E03;
!!$    pc%liqViscosity%C = -1.1468E00;
!!$    pc%liqViscosity%D = 0;
!!$    pc%liqViscosity%E = 0;
!!$    pc%liqViscosity%range(:,1) = (/ 220.00, 1.8628E-01; 503.00, 6.4431E-05 /);
!!$    pc%liqViscosity%eqn = 101;
!!$
!!$    pc%vapViscosity%A = 1.6196E-07;
!!$    pc%vapViscosity%B = 7.4453E-01;
!!$    pc%vapViscosity%C = 1.4579E02;
!!$    pc%vapViscosity%D = 0;
!!$    pc%vapViscosity%E = 0;
!!$    pc%vapViscosity%range(:,1) = (/ 193.00, 4.6420E-06; 1000.00, 2.4204E-05 /);
!!$    pc%vapViscosity%eqn = 102;
!!$
!!$    pc%liqThermalConductivity%A = 2.1782E-01;
!!$    pc%liqThermalConductivity%B = -2.4050E-04;
!!$    pc%liqThermalConductivity%C = 0;
!!$    pc%liqThermalConductivity%D = 0;
!!$    pc%liqThermalConductivity%E = 0;
!!$    pc%liqThermalConductivity%range(:,1) = (/ 193.00, 1.7140E-01; 450.80, 1.0940E-01 /);
!!$    pc%liqThermalConductivity%eqn = 100;
!!$
!!$    pc%vapThermalConductivity%A = 2.8570E-04;
!!$    pc%vapThermalConductivity%B = 9.0229E-01;
!!$    pc%vapThermalConductivity%C = 7.0294E02;
!!$    pc%vapThermalConductivity%D = 0;
!!$    pc%vapThermalConductivity%E = 0;
!!$    pc%vapThermalConductivity%range(:,1) = (/ 450.80, 2.7698E-02; 1000.00, 8.5424E-02 /);
!!$    pc%vapThermalConductivity%eqn = 102; 
!!$
!!$    pc%SurfaceTension%A = 6.5073E-02;
!!$    pc%SurfaceTension%B = -8.8258E-05;
!!$    pc%SurfaceTension%C = -1.6138E-08;
!!$    pc%SurfaceTension%D = 0;
!!$    pc%SurfaceTension%E = 0;
!!$    pc%SurfaceTension%range(:,1) = (/ 294.75, 3.7657E-02; 450.80, 2.2007E-02 /);
!!$    pc%SurfaceTension%eqn = 100;    
!!$
!!$  end subroutine pc_tetrahydrofurfurylalcohol

end module pc_database
