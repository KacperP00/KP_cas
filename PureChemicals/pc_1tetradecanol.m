function out = pc_1tetradecanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1tetradecanol','n-tetradecanol','n-Tetradecanol','1Tetradecanol','1-tetradecanol','1-Tetradecanol'}; 
    out.CAS = '112-72-1';
    % components
    out.ChemicalFormular.C = 14;         
    out.ChemicalFormular.H = 30;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 214.392;  % kg/kmol    
    out.Tcrit = 747.00;             % K    
    out.pcrit = 1.8100E06;          % Pa 
    out.vcrit = 0.80200;            % m^3/kmol 
    out.CritCompressFactor = 0.234; % -
    
%% 2 
    out.MeltingPoint = 310.65;              % K
    out.TripletPointT = 310.65;             % K
    out.TripletPointp = 7.9998E-02;         % Pa
    out.NormalBoilingPoint = 568.80;        % K
    out.liqMolarVolume = 0.260434;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -4.7590E08;             % J/kmol
    out.IG_GibbsOfFormation = -6.2200E07;            % J/kmol
    out.IG_AbsoluteEntropy =  7.5390E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.5100E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -8.4914E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.7432;            % -
    out.RadiusOfGyration = 6.7300E-10;      % m
    out.SolubilityParameter = 1.9430E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.1702E-30;          % C^m
    out.VanDerWaalsVolume = 0.15470;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 2.1130E09;        % m^2/kmol
    out.RefractiveIndex = 1.4454;           % -
    out.FlashPoint = 414.15;                % K
    out.FlammabilityLimits = [0.5 5.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 550.00;          % K
    
%% 6 
    out.solDensity.A = 4.4513E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [294.00, 4.4513E00; 298.15, 4.4513E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.2737E-01;
    out.liqDensity.B = 2.6236E-01;
    out.liqDensity.C = 7.4700E02;
    out.liqDensity.D = 3.2417E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [310.065, 3.8397E00; 747.00, 1.2478E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  4.1146E02;
    out.vapPressure.B = -2.6293E04;
    out.vapPressure.C = -5.9876E01;
    out.vapPressure.D =  4.5908E-02;
    out.vapPressure.E =  1.0000E00;
    out.vapPressure.range = [310.65, 7.9998E-02; 747.00, 1.8151E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 2.6345E08;
    out.HeatOfVap.B = 2.6157E00;
    out.HeatOfVap.C = -1.9979E00;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [310.65, 1.0092E08; 747.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -4.4453E05;
    out.solHeatCapacity.B = 2.9100E03;
    out.solHeatCapacity.C = -8.0430E00;
    out.solHeatCapacity.D = 2.3700E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [269.00, 3.3826E05; 310.65, 4.5946E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = -4.4041E06;
    out.liqHeatCapacity.B = 3.7646E04;
    out.liqHeatCapacity.C = -1.0306E02;
    out.liqHeatCapacity.D = 1.2078E-01;
    out.liqHeatCapacity.E = -4.9723E-05;
    out.liqHeatCapacity.range = [324.00, 5.3443E05; 568.80, 6.8753E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.3340E05;
    out.IG_HeatCapacity.B = 7.2700E05;
    out.IG_HeatCapacity.C = 1.5623E03;
    out.IG_HeatCapacity.D = 5.4210E05;
    out.IG_HeatCapacity.E = 7.2300E02;
    out.IG_HeatCapacity.range = [150.00, 2.3668E05; 1500.00, 8.4672E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 7.2937E-01;
    out.SecondViralCoef.B = - 9.0895E02;
    out.SecondViralCoef.C = -2.3210E08;
    out.SecondViralCoef.D = 1.0921E22;
    out.SecondViralCoef.E = -4.4086E24;
    out.SecondViralCoef.range = [373.50, -8.4889E00; 1500.00, 5.4944E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1007E01;
    out.liqViscosity.B = 3.5413E03;
    out.liqViscosity.C = -8.1269E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [310.65, 1.3966E-02; 568.80, 4.8381E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 9.3882E-07;
    out.vapViscosity.B = 4.7630E-01;
    out.vapViscosity.C = 8.4022E02;
    out.vapViscosity.D = -5.3350E04;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [310.65, 4.5823E-06; 1000.00, 1.4106E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.3150E-01;
    out.liqThermalConductivity.B = -1.9870E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [310.65, 1.6977E-01; 568.80, 1.1848E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -3.0040E00;
    out.vapThermalConductivity.B = -3.4874E-02;
    out.vapThermalConductivity.C = - 1.7663E03;
    out.vapThermalConductivity.D = -2.8478E07;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [568.80, 2.6716E-02; 1000.0, 8.0730E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.9479E-02;
    out.SurfaceTension.B = 1.3005E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [310.65, 2.9561E-02; 747.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


