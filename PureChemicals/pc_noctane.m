function out = pc_noctane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'noctane','n-octane','octyl hydride'}; 
    out.CAS = '111-65-9';
    % components
    out.ChemicalFormular.C = 8;         
    out.ChemicalFormular.H = 18;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 114.231;   % kg/kmol    
    out.Tcrit = 568.70;              % K    
    out.pcrit = 2.4900E6;            % Pa 
    out.vcrit = 0.48600;             % m^3/kmol 
    out.CritCompressFactor = 0.256;  % -
    
%% 2 
    out.MeltingPoint = 216.38;             % K
    out.TripletPointT = 216.38;            % K
    out.TripletPointp = 2.1083E00;         % Pa
    out.NormalBoilingPoint = 398.83;       % K
    out.liqMolarVolume  = 0.163374;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.0875E08;            % J/kmol
    out.IG_GibbsOfFormation = 1.6000E07;            % J/kmol
    out.IG_AbsoluteEntropy = 4.6723E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 2.0740E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -5.0741E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3996;            % -
    out.RadiusOfGyration = 4.5460E-10;      % m
    out.SolubilityParameter = 1.5400E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.08872;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.3690E09;        % m^2/kmol
    out.RefractiveIndex = 1.3951;           % -
    out.FlashPoint = 286.00;                % K
    out.FlammabilityLimits = [0.8 6.5];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 479.00;          % K
    
%% 6 
    out.solDensity.A = 8.3409E00;
    out.solDensity.B = -3.1515E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [133.15, 7.9213E00; 216.38, 7.6590E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.3731E-01;
    out.liqDensity.B = 2.6115E-01;
    out.liqDensity.C = 5.6870E02;
    out.liqDensity.D = 2.8034E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [216.38, 6.6558E00; 568.70, 2.0575E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  9.6084E01;
    out.vapPressure.B = -7.9002E03;
    out.vapPressure.C = -1.1003E01;
    out.vapPressure.D =  7.1802E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [216.38, 2.1083E00; 568.70, 2.4673E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.5180E07;
    out.HeatOfVap.B = 3.8467E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [216.38, 4.5898E07; 568.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -2.4000E04;
    out.solHeatCapacity.B = 1.9472E03;
    out.solHeatCapacity.C = -8.5360E00;
    out.solHeatCapacity.D = 1.3400E-02;
    out.solHeatCapacity.E = 2.0940E-05;
    out.solHeatCapacity.range = [20.00, 1.1640E04; 216.37, 1.7933E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.2483E05;
    out.liqHeatCapacity.B = -1.8663E02;
    out.liqHeatCapacity.C = 9.5891E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [216.38, 2.2934E05; 460.00, 3.4189E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.3554E05;
    out.IG_HeatCapacity.B = 4.4310E05;
    out.IG_HeatCapacity.C = 1.6356E03;
    out.IG_HeatCapacity.D = 3.0540E05;
    out.IG_HeatCapacity.E = 7.4640E02;
    out.IG_HeatCapacity.range = [200.00, 1.4529E05; 1500.00, 4.9764E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.7390E-01;
    out.SecondViralCoef.B = -3.2144E02;
    out.SecondViralCoef.C = -6.6830E07;
    out.SecondViralCoef.D = -1.2710E20;
    out.SecondViralCoef.E = 1.6103E22;
    out.SecondViralCoef.range = [284.38, -5.4098E00; 1500.00, 3.9801E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.0463E01;
    out.liqViscosity.B = 1.4974E03;
    out.liqViscosity.C = 1.3790E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [216.38, 2.1813E-03; 398.83, 2.1383E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.1191E-08;
    out.vapViscosity.B = 9.2925E-01;
    out.vapViscosity.C = 5.5092E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [216.38, 3.6772E-06; 1000.00, 1.8134E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1560E-01;
    out.liqThermalConductivity.B = -2.9483E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [216.38, 1.5180E-01; 398.83, 9.8013E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -8.7580E03;
    out.vapThermalConductivity.B = 8.4480E-01;
    out.vapThermalConductivity.C = -2.7121E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [339.00, 1.5025E-02; 1000.00, 1.1053E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.2789E-02;
    out.SurfaceTension.B = 1.2323E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [216.38, 2.9261E-02; 568.70, 0];
    out.SurfaceTension.eqn = 106;    
    
end


