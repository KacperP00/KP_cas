function out = pc_isopentane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'isopentane','iso-pentane','butane','2-methylbutane'}; 
    out.CAS = '107-83-5';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 72.1503;  % kg/kmol    
    out.Tcrit = 460.4;              % K    
    out.pcrit = 3.3800E6;           % Pa 
    out.vcrit = 0.306;              % m^3/kmol 
    out.CritCompressFactor = 0.27;  % -
    
%% 2 
    out.MeltingPoint = 113.25;              % K
    out.TripletPointT = 113.25;             % K
    out.TripletPointp = 1.2144E-4;          % Pa
    out.NormalBoilingPoint = 300.994;       % K
    out.liqMolarVolume  = 0.117065;         % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.5370E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.4050E7;            % J/kmol
    out.IG_AbsoluteEntropy =  3.4374E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  5.1500E06;         % J/kmol
    out.StandardNetHeatOfCombustion = -3.2395E9;    % J/kmol
    out.StandardGibbsOfFormation = -1.4420E07;      % J/kmol
    out.StandardHeatOfFormation = -1.7891E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.6041E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.227875;            % -
    out.RadiusOfGyration = 3.3240E-10;      % m
    out.SolubilityParameter = 1.3870E4;     % (J/m^3)^.5 
    out.DipoleMoment = 4.3400E-30;                   % C^m
    out.VanDerWaalsVolume = 0.05802;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.2800E08;        % m^2/kmol
    out.RefractiveIndex = 1.35088;          % -
    out.FlashPoint = 217.15;                % K
    out.FlammabilityLimits = [1.3 8.0];     % vol-%
    out.FlammabilityLimitTemps = [216 244]; % K
    out.AutoignitionTemp = 693.15;          % K
    
%% 6 
    out.solDensity.A = 1.2600E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [113.25, 1.2600E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.1991E-01;
    out.liqDensity.B = 2.7815E-01;
    out.liqDensity.C = 4.6040E02;
    out.liqDensity.D = 2.8667E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [113.25, 1.0764E01; 460.40, 3.3072E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.1308E01;
    out.vapPressure.B = -4.9760E03;
    out.vapPressure.C = -7.7169E00;
    out.vapPressure.D =  8.7271E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [113.25, 1.2144E-04; 460.40,3.3656E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 3.7593E07;
    out.HeatOfVap.B = 3.9173E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [113.25, 3.3657E07; 460.40, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -8.1500E03;
    out.solHeatCapacity.B = 8.5500E02;
    out.solHeatCapacity.C = 1.4470E01;
    out.solHeatCapacity.D = -2.7730E-01;
    out.solHeatCapacity.E = 1.3350E-03;
    out.solHeatCapacity.range = [20.00, 1.2733E04; 113.36, 9.1223E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.0830E05;
    out.liqHeatCapacity.B = 1.4600E02;
    out.liqHeatCapacity.C = -2.9200E-01;
    out.liqHeatCapacity.D = 1.5100E-03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [113.25, 1.2328E05; 310.00, 1.7048E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 7.4600E04;
    out.IG_HeatCapacity.B = 3.2650E05;
    out.IG_HeatCapacity.C = 1.5450E03;
    out.IG_HeatCapacity.D = 1.9230E05;
    out.IG_HeatCapacity.E = 6.6670E02;
    out.IG_HeatCapacity.range = [200.00, 8.5462E04; 1500.00, 3.3792E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.8170E-01;
    out.SecondViralCoef.B = -1.9553E02;
    out.SecondViralCoef.C = -1.5973E07;
    out.SecondViralCoef.D = -4.6950E18;
    out.SecondViralCoef.E =3.6600E20;
    out.SecondViralCoef.range = [230.20, -2.3708E00; 1500.00, 4.6614E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.2596E01;
    out.liqViscosity.B = 8.8911E02;
    out.liqViscosity.C = 2.0469E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [150.00, 3.5423E-03; 310.00, 1.9283E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.4344E-08;
    out.vapViscosity.B = 9.7376E-01;
    out.vapViscosity.C = -9.1597E01;
    out.vapViscosity.D = 1.8720E04;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [150.00, 2.6214E-06; 1000.00, 2.1905E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1246E-01;
    out.liqThermalConductivity.B = -3.3581E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [113.25,1.7443E-01; 368.13, 8.8838E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 8.9680E-04;
    out.vapThermalConductivity.B = 7.7420E-01;
    out.vapThermalConductivity.C = 4.5600E02;
    out.vapThermalConductivity.D = 2.3064E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.15, 1.1981E-02; 1000.0, 1.1176E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0876E-02;
    out.SurfaceTension.B = 1.2066E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [113.25, 3.6188E-02; 460.40, 0];
    out.SurfaceTension.eqn = 106;    
    
end


