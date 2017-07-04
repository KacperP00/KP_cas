function out = pc_ethanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'ethanol','Ethanol'}; 
    out.CAS = '64-17-5';
    % components
    out.ChemicalFormular.C = 2;         
    out.ChemicalFormular.H = 6;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 46.069;   % kg/kmol    
    out.Tcrit = 513.92;             % K    
    out.pcrit = 6.1480E6;           % Pa 
    out.vcrit = 0.16700;            % m^3/kmol 
    out.CritCompressFactor = 0.240; % -
    
%% 2 
    out.MeltingPoint = 159.05;              % K
    out.TripletPointT = 159.05;             % K
    out.TripletPointp = 4.8459E-4;          % Pa
    out.NormalBoilingPoint = 351.44;        % K
    out.liqMolarVolume = 0.058618;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.3495E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.6785E8;            % J/kmol
    out.IG_AbsoluteEntropy =  2.8064E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  4.9310E6;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.2350E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.6452;            % -
    out.RadiusOfGyration = 2.2590E-10;      % m
    out.SolubilityParameter = 2.6130E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.6400E-30;          % C^m
    out.VanDerWaalsVolume = 0.03194;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 4.9300E08;        % m^2/kmol
    out.RefractiveIndex = 1.3594;           % -
    out.FlashPoint = 286.00;                % K
    out.FlammabilityLimits = [4.3 19.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 696.00;          % K
    
%% 6 
    out.solDensity.A = 2.29001E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [159.05, 2.29E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.6480E00;
    out.liqDensity.B = 2.7627E-01;
    out.liqDensity.C = 5.1391E02;
    out.liqDensity.D = 2.3310E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [159.05, 1.9413E01; 513.92, 5.9652E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.4475E01;
    out.vapPressure.B = -7.1643E03;
    out.vapPressure.C = -7.3270E00;
    out.vapPressure.D =  3.1340E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [159.05, 4.8459E-04; 513.92, 6.1171E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.6900E07;
    out.HeatOfVap.B = 3.3590E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [159.05, 5.0245E07; 513.92, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.3500E04;
    out.solHeatCapacity.B = 1.1755E03;
    out.solHeatCapacity.C = -8.0430E00;
    out.solHeatCapacity.D = 2.3700E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [25.00, 1.1231E04; 159.00, 6.5336E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.0264E05;
    out.liqHeatCapacity.B = -1.3963E02;
    out.liqHeatCapacity.C = -3.0341E-02;
    out.liqHeatCapacity.D = 2.0386E-03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [159.05, 8.7867E04; 390.00, 1.6450E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 4.8200E04;
    out.IG_HeatCapacity.B = 1.4577E05;
    out.IG_HeatCapacity.C = 1.6628E03;
    out.IG_HeatCapacity.D = 9.3900E04;
    out.IG_HeatCapacity.E = 7.4470E02;
    out.IG_HeatCapacity.range = [200.00, 5.2236E04; 1500.00, 1.6576E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.4000E-02;
    out.SecondViralCoef.B = -5.5700E01;
    out.SecondViralCoef.C = -1.2900E07;
    out.SecondViralCoef.D = -6.4900E19;
    out.SecondViralCoef.E = -2.5800E22;
    out.SecondViralCoef.range = [309.72, -2.3203E00; 2581.00, 2.1669E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = 7.8750E00;
    out.liqViscosity.B = 7.8198E02;
    out.liqViscosity.C = -3.0418E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [200.00, 1.3148E-02; 440.00, 1.4160E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.0613E-07;
    out.vapViscosity.B = 8.0660E-01;
    out.vapViscosity.C = 5.2700E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [200.00, 6.0295E-06; 1000.00, 2.6505E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.4680E-01;
    out.liqThermalConductivity.B = -2.6400E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [159.05, 2.0481E-01; 353.15, 1.5357E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -1.0109E-02;
    out.vapThermalConductivity.B = 6.4750E-01;
    out.vapThermalConductivity.C = -7.3320E03;
    out.vapThermalConductivity.D = -2.6800E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [293.15, 1.4747E-02; 1000.0, 1.3417E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 3.7640E-02;
    out.SurfaceTension.B = -2.1570E-05;
    out.SurfaceTension.C = -1.0250E-07;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [273.15, 2.4101E-02; 503.15, 8.3816E-04];
    out.SurfaceTension.eqn = 100;    
    
end


