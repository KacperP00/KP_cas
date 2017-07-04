subroutine pc_tetrahydrofurfurylalcohol(pc)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------------------

 

 
    ! all alternative names
    pc%name = {'tetrahydrofurfurylalcohol','tetrahydrofuryl carbinol','oxolan-2-methanol'}; 
    pc%CAS = '97-99-4';
    ! components
    pc%ChemicalFormula%C = 5;         
    pc%ChemicalFormula%H = 10;
    pc%ChemicalFormula%O = 2;
 
    pc%MolecularWeight = 102.133;   ! kg/kmol    
    pc%Tcrit = 639.00;              ! K    
    pc%pcrit = 4.6600E6;            ! Pa 
    pc%vcrit = 0.29000;             ! m^3/kmol 
    pc%CritCompressFactor = 0.254;  ! -
    
 
    pc%MeltingPoint = MAX_REAL_WP;             ! K
    pc%TripletPointT = MAX_REAL_WP;            ! K
    pc%TripletPointp = MAX_REAL_WP;            ! Pa
    pc%NormalBoilingPoint = 450.80;    ! K
    pc%liqMolarVolume  = 0.097319;     ! m^3/kmol
    
 
    pc%IG_HeatOfFormation = -3.6920E08;            ! J/kmol
    pc%IG_GibbsOfFormation = -2.1600E08;           ! J/kmol
    pc%IG_AbsoluteEntropy = 3.7400E5;              ! J/kmolK
    pc%HeatFusionAtMeltPoint = MAX_REAL_WP;                ! J/kmol
    pc%StandardNetHeatOfCombustion = -2.7412E9;    ! J/kmol
    pc%StandardGibbsOfFormation = MAX_REAL_WP;             ! J/kmol
    pc%StandardHeatOfFormation = MAX_REAL_WP;              ! J/kmol
    pc%StandardAbsoluteEntropy = MAX_REAL_WP;              ! J/kmolK    
    
 
    pc%AcentricFactor = 0.7044;            ! -
    pc%RadiusOfGyration = 3.4810E-10;      ! m
    pc%SolubilityParameter = 2.2920E4;     ! (J/m^3)^.5 
    pc%DipoleMoment = 7.3400E-30;          ! C^m
    pc%VanDerWaalsVolume = 0.05830;        ! m^3/kmol
    
 
    pc%VanDerWaalsArea = 7.4600E08;        ! m^2/kmol
    pc%RefractiveIndex = 1.4499;           ! -
    pc%FlashPoint = 343.15;                ! K
    pc%FlammabilityLimits = (/ 1.5 9.7 /);     ! vol-!
    pc%FlammabilityLimitTemps = (/ MAX_REAL_WP MAX_REAL_WP /); ! K
    pc%AutoignitionTemp = 553.15;          ! K
    
 
    pc%solDensity%A = 1.2900E01;
    pc%solDensity%B = 0;
    pc%solDensity%C = 0;
    pc%solDensity%D = 0;
    pc%solDensity%E = 0;
    pc%solDensity%range(:,1) = (/ 0.0, 1.2900E01 /);
    pc%solDensity%eqn = 100;
    
    pc%liqDensity%A = 9.7001E-01;
    pc%liqDensity%B = 2.8130E-01;
    pc%liqDensity%C = 6.3900E02;
    pc%liqDensity%D = 2.3837E-01;
    pc%liqDensity%E = 0;
    pc%liqDensity%range(:,1) = (/ 193.00, 1.1046E01; 639.00, 3.4483E00 /);
    pc%liqDensity%eqn = 105;
    
    pc%vapPressure%A =  1.5059E02;
    pc%vapPressure%B = -1.1574E04;
    pc%vapPressure%C = -1.9025E01;
    pc%vapPressure%D =  1.4141E-05;
    pc%vapPressure%E =  2.0000E00;
    pc%vapPressure%range(:,1) = (/ 193.00, 1.2655E-04; 639.00, 4.6469E06 /);
    pc%vapPressure%eqn = 101;  
    
    pc%HeatOfVap%A = 6.4109E07;
    pc%HeatOfVap%B = 2.8538E-01;
    pc%HeatOfVap%C = 0;
    pc%HeatOfVap%D = 0;
    pc%HeatOfVap%E = 0;
    pc%HeatOfVap%range(:,1) = (/ 193.00, 5.7857E07; 639.00, 0 /);
    pc%HeatOfVap%eqn = 106;     

    pc%solHeatCapacity%A = 0;
    pc%solHeatCapacity%B = 0;
    pc%solHeatCapacity%C = 0;
    pc%solHeatCapacity%D = 0;
    pc%solHeatCapacity%E = 0;
    pc%solHeatCapacity%range(:,1) = (/ MAX_REAL_WP, MAX_REAL_WP; MAX_REAL_WP, MAX_REAL_WP /);
    pc%solHeatCapacity%eqn = 100;
    
    pc%liqHeatCapacity%A = 5.2700E04;
    pc%liqHeatCapacity%B = 4.3580E02;
    pc%liqHeatCapacity%C = 0;
    pc%liqHeatCapacity%D = 0;
    pc%liqHeatCapacity%E = 0;
    pc%liqHeatCapacity%range(:,1) = (/ 296.65, 1.8198E05; 469.50, 2.5731E05 /);
    pc%liqHeatCapacity%eqn = 100;    
    
    pc%IG_HeatCapacity%A = 2.6823E05;
    pc%IG_HeatCapacity%B = 1.1670E05;
    pc%IG_HeatCapacity%C = 1.6860E03;
    pc%IG_HeatCapacity%D = -3.5200E05;
    pc%IG_HeatCapacity%E = 2.4520E02;
    pc%IG_HeatCapacity%range(:,1) = (/ 300.00, 1.3998E05; 1500.00, 3.3693E05 /);
    pc%IG_HeatCapacity%eqn = 107;
    
    pc%SecondViralCoef%A = 1.9766E-01;
    pc%SecondViralCoef%B = -1.8908E02;
    pc%SecondViralCoef%C = -8.1632E07;
    pc%SecondViralCoef%D = -5.5857E20;
    pc%SecondViralCoef%E = 8.4850E22;
    pc%SecondViralCoef%range(:,1) = (/ 319.50, -5.5954E00; 1500.00, 4.7400E-02 /);
    pc%SecondViralCoef%eqn = 104;  
    
    pc%liqViscosity%A = -7.9742E00;
    pc%liqViscosity%B = 2.7454E03;
    pc%liqViscosity%C = -1.1468E00;
    pc%liqViscosity%D = 0;
    pc%liqViscosity%E = 0;
    pc%liqViscosity%range(:,1) = (/ 220.00, 1.8628E-01; 503.00, 6.4431E-05 /);
    pc%liqViscosity%eqn = 101;
    
    pc%vapViscosity%A = 1.6196E-07;
    pc%vapViscosity%B = 7.4453E-01;
    pc%vapViscosity%C = 1.4579E02;
    pc%vapViscosity%D = 0;
    pc%vapViscosity%E = 0;
    pc%vapViscosity%range(:,1) = (/ 193.00, 4.6420E-06; 1000.00, 2.4204E-05 /);
    pc%vapViscosity%eqn = 102;
    
    pc%liqThermalConductivity%A = 2.1782E-01;
    pc%liqThermalConductivity%B = -2.4050E-04;
    pc%liqThermalConductivity%C = 0;
    pc%liqThermalConductivity%D = 0;
    pc%liqThermalConductivity%E = 0;
    pc%liqThermalConductivity%range(:,1) = (/ 193.00, 1.7140E-01; 450.80, 1.0940E-01 /);
    pc%liqThermalConductivity%eqn = 100;
       
    pc%vapThermalConductivity%A = 2.8570E-04;
    pc%vapThermalConductivity%B = 9.0229E-01;
    pc%vapThermalConductivity%C = 7.0294E02;
    pc%vapThermalConductivity%D = 0;
    pc%vapThermalConductivity%E = 0;
    pc%vapThermalConductivity%range(:,1) = (/ 450.80, 2.7698E-02; 1000.00, 8.5424E-02 /);
    pc%vapThermalConductivity%eqn = 102; 
    
    pc%SurfaceTension%A = 6.5073E-02;
    pc%SurfaceTension%B = -8.8258E-05;
    pc%SurfaceTension%C = -1.6138E-08;
    pc%SurfaceTension%D = 0;
    pc%SurfaceTension%E = 0;
    pc%SurfaceTension%range(:,1) = (/ 294.75, 3.7657E-02; 450.80, 2.2007E-02 /);
    pc%SurfaceTension%eqn = 100;    
    
end subroutine pc_tetrahydrofurfurylalcohol


