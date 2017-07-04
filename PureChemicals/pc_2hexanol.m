function out = pc_2hexanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2hexanol','2-hexanol','2-Hexanol','2Hexanol'}; 
    out.CAS = '626-93-7';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 102.77;   % kg/kmol    
    out.Tcrit = 586.20;             % K    
    out.pcrit = 3.4000E06;          % Pa 
    out.vcrit = 0.38000;            % m^3/kmol 
    out.CritCompressFactor = 0.265; % -
    
%% 2 
    out.MeltingPoint = 223.00;              % K
    out.TripletPointT = 223.00;             % K
    out.TripletPointp = 9.4925E-02;         % Pa
    out.NormalBoilingPoint = 413.04;        % K
    out.liqMolarVolume = 0.126073;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.3350E08;            % J/kmol
    out.IG_GibbsOfFormation = -1.5020E08;           % J/kmol
    out.IG_AbsoluteEntropy = 4.3600E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = NaN;                % J/kmol
    out.StandardNetHeatOfCombustion = -3.6660E09;   % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5655;            % -
    out.RadiusOfGyration = 3.7870E-10;      % m
    out.SolubilityParameter = 2.1115E04;    % (J/m^3)^.5 
    out.DipoleMoment = NaN;                 % C^m
    out.VanDerWaalsVolume = 0.07285;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.0320E09;        % m^2/kmol
    out.RefractiveIndex = 1.4128;           % -
    out.FlashPoint = 314.15;                % K
    out.FlammabilityLimits = [1.2 8.3];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = NaN;             % K
    
%% 6 
    out.solDensity.A = 9.9600E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [223.00, 9.9600E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.9120E-01;
    out.liqDensity.B = 2.6275E-01;
    out.liqDensity.C = 5.8620E02;
    out.liqDensity.D = 2.6945E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [223.00, 8.5168E00; 586.20, 2.6306E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.2842E02;
    out.vapPressure.B = -1.0914E04;
    out.vapPressure.C = -1.5189E01;
    out.vapPressure.D =  5.9636E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [223.00, 9.4925E-02; 586.20, 3.4076E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 9.1260E07;
    out.HeatOfVap.B = 6.2140E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [223.00, 6.7778E07; 586.20, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 1.9600E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [223.00, 1.9600E05; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.6063E05;
    out.liqHeatCapacity.B = -1.9047E03;
    out.liqHeatCapacity.C = 5.2540E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [230.00, 2.0049E05; 330.00, 3.0424E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.2575E05;
    out.IG_HeatCapacity.B = 3.5890E05;
    out.IG_HeatCapacity.C = 1.8215E03;
    out.IG_HeatCapacity.D = 2.5557E05;
    out.IG_HeatCapacity.E = 8.1300E02;
    out.IG_HeatCapacity.range = [300.00, 1.5898E05; 1200.00, 3.7756E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.1711E-01;
    out.SecondViralCoef.B = -7.5358E01;
    out.SecondViralCoef.C = -1.0028E08;
    out.SecondViralCoef.D = 9.8894E20;
    out.SecondViralCoef.E = -4.2519E23;
    out.SecondViralCoef.range = [293.10, -1.2600E01; 2923.10, 8.7315E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.2444E02;
    out.liqViscosity.B = 8.3425E03;
    out.liqViscosity.C = 1.5972E01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [230.00, 2.6948E-01; 413.04, 3.2420E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.0200E-07;
    out.vapViscosity.B = 7.7710E-01;
    out.vapViscosity.C = 1.0050E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [223.00, 4.6978E-06; 1000.00, 1.9875E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1560E-01;
    out.liqThermalConductivity.B = -2.6510E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [223.00, 1.5648E-01; 463.00, 9.2859E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -7.9554E06;
    out.vapThermalConductivity.B = -1.6390E-01;
    out.vapThermalConductivity.C = 2.4484E09;
    out.vapThermalConductivity.D = -2.4038E13;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [413.04, 2.1961E-02; 1000.0, 1.1877E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1290E-02;
    out.SurfaceTension.B = 1.0435E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [223.00, 3.1124E-02; 586.20, 0];
    out.SurfaceTension.eqn = 106;    
    
end


