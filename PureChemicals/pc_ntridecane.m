function out = pc_ntridecane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'ntridecane','n-tridecane','tridecyl hydride'}; 
    out.CAS = '629-50-5';
    % components
    out.ChemicalFormular.C = 13;         
    out.ChemicalFormular.H = 28;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 184.365;   % kg/kmol    
    out.Tcrit = 675.00;              % K    
    out.pcrit = 1.6800E6;            % Pa 
    out.vcrit = 0.77500;             % m^3/kmol 
    out.CritCompressFactor = 0.232;  % -
    
%% 2 
    out.MeltingPoint = 267.76;             % K
    out.TripletPointT = 267.76;            % K
    out.TripletPointp = 2.5096E-02;        % Pa
    out.NormalBoilingPoint = 508.62;       % K
    out.liqMolarVolume  = 0.244631;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.1177E08;            % J/kmol
    out.IG_GibbsOfFormation = 5.7710E07;            % J/kmol
    out.IG_AbsoluteEntropy = 6.6337E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 2.8500E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -8.1229E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.6174;            % -
    out.RadiusOfGyration = 6.1960E-10;      % m
    out.SolubilityParameter = 1.6100E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.13960;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.9090E09;        % m^2/kmol
    out.RefractiveIndex = 1.4235;           % -
    out.FlashPoint = 352.00;                % K
    out.FlammabilityLimits = [0.6 4.7];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 475.00;          % K
    
%% 6 
    out.solDensity.A = 5.0311E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [0, 5.0311E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.2160E-01;
    out.liqDensity.B = 2.5040E-01;
    out.liqDensity.C = 6.7500E02;
    out.liqDensity.D = 3.0710E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [267.76, 4.2035E00; 675.00, 1.2843E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.3745E02;
    out.vapPressure.B = -1.2549E04;
    out.vapPressure.C = -1.6543E01;
    out.vapPressure.D =  7.1275E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [267.76, 2.5096E-01; 675.00, 1.6786E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 8.4339E07;
    out.HeatOfVap.B = 4.2570E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [267.76, 6.8015E07; 675.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -3.0600E04;
    out.solHeatCapacity.B = 2.3520E03;
    out.solHeatCapacity.C = -3.2440E00;
    out.solHeatCapacity.D = -3.5173E-02;
    out.solHeatCapacity.E = 1.3218E-04;
    out.solHeatCapacity.range = [20.00, 1.4882E04; 255.00, 3.3389E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.5018E05;
    out.liqHeatCapacity.B = -1.0470E02;
    out.liqHeatCapacity.C = 1.0022E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [267.76, 3.9400E05; 508.62, 5.5619E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.1496E05;
    out.IG_HeatCapacity.B = 7.3045E05;
    out.IG_HeatCapacity.C = 1.6695E03;
    out.IG_HeatCapacity.D = 4.9998E05;
    out.IG_HeatCapacity.E = 7.4102E02;
    out.IG_HeatCapacity.range = [200.00, 2.3156E05; 1500.00, 8.0251E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.4720E-01;
    out.SecondViralCoef.B = -4.7897E02;
    out.SecondViralCoef.C = -2.3434E08;
    out.SecondViralCoef.D = -1.7856E21;
    out.SecondViralCoef.E = 3.0980E23;
    out.SecondViralCoef.range = [337.50, -1.2222E01; 1500.00, 5.8391E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.1008E01;
    out.liqViscosity.B = 2.0425E03;
    out.liqViscosity.C = 1.3692E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [267.76, 3.2603E-03; 508.62, 2.1182E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.5585E-08;
    out.vapViscosity.B = 8.9870E-01;
    out.vapViscosity.C = 1.6530E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [267.76, 3.3441E-06; 1000.00, 1.5168E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0447E-01;
    out.liqThermalConductivity.B = -2.2612E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [267.76, 1.4392E-01; 508.62, 8.9461E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 5.3701E-06;
    out.vapThermalConductivity.B = 1.4751E00;
    out.vapThermalConductivity.C = 5.9909E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [508.62, 2.4219E-02; 1000.00, 8.9415E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.5322E-02;
    out.SurfaceTension.B = 1.3257E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [267.76, 2.8312E-02; 675.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


