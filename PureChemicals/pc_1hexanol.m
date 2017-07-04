function out = pc_1hexanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1hexanol','1-hexanol','1-Hexanol','1Hexanol','n-Hexanol','n-hexanol',}; 
    out.CAS = '111-27-3';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 102.77;   % kg/kmol    
    out.Tcrit = 610.30;             % K    
    out.pcrit = 3.4170E06;          % Pa 
    out.vcrit = 0.38700;            % m^3/kmol 
    out.CritCompressFactor = 0.261; % -
    
%% 2 
    out.MeltingPoint = 228.55;              % K
    out.TripletPointT = 228.55;             % K
    out.TripletPointp = 2.9550E-02;         % Pa
    out.NormalBoilingPoint = 430.55;        % K
    out.liqMolarVolume = 0.125200;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.1650E08;            % J/kmol
    out.IG_GibbsOfFormation = -1.3440E08;           % J/kmol
    out.IG_AbsoluteEntropy = 4.4010E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = 1.5400E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.6766E09;   % J/kmol
    out.StandardGibbsOfFormation = - 1.4980E+08;    % J/kmol
    out.StandardHeatOfFormation = -3.7750E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.8740E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5764;            % -
    out.RadiusOfGyration = 4.1440E-10;      % m
    out.SolubilityParameter = 2.1830E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.5000E-30;          % C^m
    out.VanDerWaalsVolume = 0.07286;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.0330E09;        % m^2/kmol
    out.RefractiveIndex = 1.4161;           % -
    out.FlashPoint = 356.00;                % K
    out.FlammabilityLimits = [1.2 8.2];     % vol-%
    out.FlammabilityLimitTemps = [332 366]; % K
    out.AutoignitionTemp = 558.00;          % K
    
%% 6 
    out.solDensity.A = 9.8900E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [228.55, 9.8900E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.8425E-01;
    out.liqDensity.B = 2.6477E-01;
    out.liqDensity.C = 6.1030E02;
    out.liqDensity.D = 2.4395E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [228.55, 8.4540E00; 610.30, 2.5843E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.2508E02;
    out.vapPressure.B = -1.1692E04;
    out.vapPressure.C = -1.4258E01;
    out.vapPressure.D =  1.1102E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [228.55, 2.9550E-02; 610.30, 3.4221E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 8.9962E07;
    out.HeatOfVap.B = 5.5179E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [228.55, 6.9442E07; 610.30, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.8080E04;
    out.solHeatCapacity.B = 1.4083E03;
    out.solHeatCapacity.C = -4.2360E00;
    out.solHeatCapacity.D = -1.1200E-03;
    out.solHeatCapacity.E = 3.0800E-05;
    out.solHeatCapacity.range = [15.00, 2.0892E+03; 200.00, 1.3446E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.6386E06;
    out.liqHeatCapacity.B = -1.7261E04;
    out.liqHeatCapacity.C = 7.1721E01;
    out.liqHeatCapacity.D = -1.2026E-01;
    out.liqHeatCapacity.E = 7.1087E-05;
    out.liqHeatCapacity.range = [228.55, 1.9821E05; 460.00, 3.5197E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0625E05;
    out.IG_HeatCapacity.B = 3.5210E05;
    out.IG_HeatCapacity.C = 1.5835E03;
    out.IG_HeatCapacity.D = 2.4620E05;
    out.IG_HeatCapacity.E = 7.1575E+02;
    out.IG_HeatCapacity.range = [200.00, 1.1607E05; 1500.00, 3.9726E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.5305E-01;
    out.SecondViralCoef.B = -2.6298E02;
    out.SecondViralCoef.C = -7.5678E07;
    out.SecondViralCoef.D = -2.2383E20;
    out.SecondViralCoef.E = 2.0576E22;
    out.SecondViralCoef.range = [305.15, -5.352E00; 1500.00, 5.5299E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -4.3289E01;
    out.liqViscosity.B = 4.0412E03;
    out.liqViscosity.C = 4.2714E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [250.00, 2.9018E-02; 430.55, 3.3662E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.5773E-07;
    out.vapViscosity.B = 7.1890E-01;
    out.vapViscosity.C = 1.6330E+02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [228.55, 4.5672E-06; 1000.00, 1.9450E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1930E-01;
    out.liqThermalConductivity.B = -2.2000E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [228.55, 1.6902E-01; 575.00, 9.2800E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -4.9355E06;
    out.vapThermalConductivity.B = -1.6530E-01;
    out.vapThermalConductivity.C = 1.5631E09;
    out.vapThermalConductivity.D = -1.5752E13;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [430.55, 2.2264E-02; 1000.0, 1.1104E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0895E-02;
    out.SurfaceTension.B = 1.0073E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [228.55, 3.1727E-02; 610.30, 0];
    out.SurfaceTension.eqn = 106;    
    
end


