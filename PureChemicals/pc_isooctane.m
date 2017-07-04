function out = pc_isooctane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'isooctane','iso-octane','heptane'}; 
    out.CAS = '107-83-5';
    % components
    out.ChemicalFormular.C = 8;         
    out.ChemicalFormular.H = 18;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 114.231;   % kg/kmol    
    out.Tcrit = 559.64;              % K    
    out.pcrit = 2.4840E6;            % Pa 
    out.vcrit = 0.48800;             % m^3/kmol 
    out.CritCompressFactor = 0.261;  % -
    
%% 2 
    out.MeltingPoint = 164.16;              % K
    out.TripletPointT = 164.16;             % K
    out.TripletPointp = 1.0606E-3;          % Pa
    out.NormalBoilingPoint = 390.80;        % K
    out.liqMolarVolume = 0.164409;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.1535E8;             % J/kmol
    out.IG_GibbsOfFormation = 1.1690E7;             % J/kmol
    out.IG_AbsoluteEntropy =  4.5957E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.1878E7;          % J/kmol
    out.StandardNetHeatOfCombustion = -5.0694E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3772;            % -
    out.RadiusOfGyration = 4.5700E-10;      % m
    out.SolubilityParameter = 1.5050E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.08871;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.2330E09;        % m^2/kmol
    out.RefractiveIndex = 1.3926;           % -
    out.FlashPoint = 277.15;                % K
    out.FlammabilityLimits = [0.9 5.8];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 520.00;          % K
    
%% 6 
    out.solDensity.A = 8.1000E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [164.16,8.1000E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.2966E-01;
    out.liqDensity.B = 2.6042E-01;
    out.liqDensity.C = 5.5964E02;
    out.liqDensity.D = 2.7016E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [164.16, 6.9236E00; 559.64, 2.0339E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  9.0727E01;
    out.vapPressure.B = -7.4715E03;
    out.vapPressure.C = -1.0243E01;
    out.vapPressure.D =  6.8840E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [164.16, 1.0606E-03; 559.64, 2.4757E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.3089E07;
    out.HeatOfVap.B = 3.8112E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [164.16, 4.6509E07; 559.64, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -2.0180E04;
    out.solHeatCapacity.B = 1.8980E03;
    out.solHeatCapacity.C = -9.1850E00;
    out.solHeatCapacity.D = 2.3100E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [20.00, 1.4291E04; 164.16, 1.4606E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.0489E05;
    out.liqHeatCapacity.B = -1.2073E02;
    out.liqHeatCapacity.C = 9.2846E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [164.16, 2.1009E05; 500.00, 3.7664E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.2215E05;
    out.IG_HeatCapacity.B = 4.7270E05;
    out.IG_HeatCapacity.C = 1.5866E03;
    out.IG_HeatCapacity.D = 3.1650E05;
    out.IG_HeatCapacity.E = 6.9960E02;
    out.IG_HeatCapacity.range = [150.00, 1.2460E05; 1500.00, 5.0777E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.7362E-01;
    out.SecondViralCoef.B = -3.4177E02;
    out.SecondViralCoef.C = -5.4179E07;
    out.SecondViralCoef.D = -9.1933E19;
    out.SecondViralCoef.E = 1.3873E22;
    out.SecondViralCoef.range = [279.82, -4.5475E00; 1500.00, 2.9717E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1342E01;
    out.liqViscosity.B = 1.0740E03;
    out.liqViscosity.C = 1.3051E-02;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [164.16, 8.7998E-03; 390.80, 2.0026E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 4.4600E-07;
    out.vapViscosity.B = 5.5355E-01;
    out.vapViscosity.C = 2.2223E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [164.16, 3.1903E-06; 1000.00, 1.6704E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.9945E-01;
    out.liqThermalConductivity.B = -2.5820E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [164.16, 1.5706E-01; 448.00 8.3776E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 5.4211E-05;
    out.vapThermalConductivity.B = 1.0809E00;
    out.vapThermalConductivity.C = -1.3266E02;
    out.vapThermalConductivity.D = 1.7358E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [390.80, 1.9106E-02; 1000.0, 9.1069E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0322E-02;
    out.SurfaceTension.B = 1.2008E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [164.16, 3.3166E-02; 559.64, 0];
    out.SurfaceTension.eqn = 106;    
    
end


