function out = pc_2methylnonane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2methylnonane','2-methylnonane'}; 
    out.CAS = '871-83-0';
    % components
    out.ChemicalFormular.C = 10;         
    out.ChemicalFormular.H = 22;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 142.265;  % kg/kmol    
    out.Tcrit = 610.00;             % K    
    out.pcrit = 2.1200E06;          % Pa 
    out.vcrit = 0.58300;            % m^3/kmol 
    out.CritCompressFactor = 0.244; % -
    
%% 2 
    out.MeltingPoint = 198.50;              % K
    out.TripletPointT = 198.50;             % K
    out.TripletPointp = 2.4868E-03;          % Pa
    out.NormalBoilingPoint = 440.15;        % K
    out.liqMolarVolume = 0.196780;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.5652E08;             % J/kmol
    out.IG_GibbsOfFormation = 2.8400E07;             % J/kmol
    out.IG_AbsoluteEntropy =  5.3802E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.7500E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -6.2888E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4723;            % -
    out.RadiusOfGyration = 5.4290E-10;      % m
    out.SolubilityParameter = 1.5440E4;     % (J/m^3)^.5 
    out.DipoleMoment = NaN;                 % C^m
    out.VanDerWaalsVolume = 0.10917;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.5030E09;        % m^2/kmol
    out.RefractiveIndex = 1.4075;           % -
    out.FlashPoint = 314.00;                % K
    out.FlammabilityLimits = [0.7 5.0];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 487.15;          % K
    
%% 6 
    out.solDensity.A = 6.5300E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [198.50, 6.5300E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 4.5413E-01;
    out.liqDensity.B = 2.6475E-01;
    out.liqDensity.C = 6.100E02;
    out.liqDensity.D = 3.0083E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [198.50, 5.5854E00; 610.00, 1.7153E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.3571E02;
    out.vapPressure.B = -1.0597E04;
    out.vapPressure.C = -1.6769E01;
    out.vapPressure.D =  1.0148E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [198.50, 2.4867E-03; 610.00, 2.1194E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.5700E07;
    out.HeatOfVap.B = 4.2540E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [198.50, 5.5570E07; 610.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 1.6222E06;
    out.solHeatCapacity.B = -5.4504E04;
    out.solHeatCapacity.C = 7.0650E02;
    out.solHeatCapacity.D = -3.9062E00;
    out.solHeatCapacity.E = 7.8986E-03;
    out.solHeatCapacity.range = [80.0, 1.0703E05; 198.50, 3.5197E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.0520E05;
    out.liqHeatCapacity.B = -4.7130E02;
    out.liqHeatCapacity.C = 1.6330E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [198.50, 2.7599E05; 440.15, 4.1412E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.5390E05;
    out.IG_HeatCapacity.B = 5.6550E05;
    out.IG_HeatCapacity.C = 1.5750E03;
    out.IG_HeatCapacity.D = 3.8870E05;
    out.IG_HeatCapacity.E = 7.0400E02;
    out.IG_HeatCapacity.range = [150.00, 1.5677E05; 1500.00, 6.1969E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.4066E-01;
    out.SecondViralCoef.B = -4.3379E02;
    out.SecondViralCoef.C = -9.7744E07;
    out.SecondViralCoef.D = -3.4420E20;
    out.SecondViralCoef.E = 5.9453E22;
    out.SecondViralCoef.range = [305.00, -6.5199E00; 1500.00, 2.2494E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.0613E01;
    out.liqViscosity.B = 1.7419E03;
    out.liqViscosity.C = 1.3450E00;
    out.liqViscosity.D = -1.2703E-28;
    out.liqViscosity.E = 1E01;
    out.liqViscosity.range = [198.50, 8.9011E-03; 593.15, 5.6989E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 7.7800E-07;
    out.vapViscosity.B = 4.6360E-01;
    out.vapViscosity.C = 3.1607E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [198.50, 3.4877E-06; 1000.00, 1.4538E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0262E-01;
    out.liqThermalConductivity.B = -2.4000E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [198.5, 1.5498E-01; 488.00, 8.5500E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 3.6874E-05;
    out.vapThermalConductivity.B = 1.1200E00;
    out.vapThermalConductivity.C = -1.5392E02;
    out.vapThermalConductivity.D = 1.7277E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [440.15, 2.1849E-02; 1000.0, 8.2911E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1634E-02;
    out.SurfaceTension.B = 1.2440E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [198.50, 3.1642E-02; 610.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


