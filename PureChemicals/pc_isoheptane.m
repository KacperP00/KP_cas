function out = pc_isoheptane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'isoheptane','iso-heptane'}; 
    out.CAS = '591-76-4';
    % components
    out.ChemicalFormular.C = 7;         
    out.ChemicalFormular.H = 16;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 100.204;  % kg/kmol    
    out.Tcrit = 530.4;              % K    
    out.pcrit = 2.7400E6;           % Pa 
    out.vcrit = 0.421;              % m^3/kmol 
    out.CritCompressFactor = 0.262; % -
    
%% 2 
    out.MeltingPoint = 154.9;              % K
    out.TripletPointT = 154.9;             % K
    out.TripletPointp = 4.2989E-03;         % Pa
    out.NormalBoilingPoint = 363.199;      % K
    out.liqMolarVolume = 0.14869;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.9460E08;            % J/kmol
    out.IG_GibbsOfFormation = 3.4700E06;            % J/kmol
    out.IG_AbsoluteEntropy =  4.2041E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  9.1840E6;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.4596E09;   % J/kmol
    out.StandardGibbsOfFormation = -2.4762E06;      % J/kmol
    out.StandardHeatOfFormation = -2.2949E08;       % J/kmol
    out.StandardAbsoluteEntropy = 3.2334E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.328006;          % -
    out.RadiusOfGyration = 4.1670E-10;      % m
    out.SolubilityParameter = 1.4730E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.07848;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.0980E09;        % m^2/kmol
    out.RefractiveIndex = 1.38228;          % -
    out.FlashPoint = 250.0;                 % K
    out.FlammabilityLimits = [1.0 6.0];     % vol-%
    out.FlammabilityLimitTemps = [260 291]; % K
    out.AutoignitionTemp = 553.15;          % K
    
%% 6 
    out.solDensity.A = 9.1400E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [154.90, 9.1400E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.3993E-01;
    out.liqDensity.B = 2.6950E-01;
    out.liqDensity.C = 5.3040E02;
    out.liqDensity.D = 2.7930E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [154.90, 7.8100E00; 530.40, 2.3745E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.9531E01;
    out.vapPressure.B = -5.8264E03;
    out.vapPressure.C = -5.4269E00;
    out.vapPressure.D =  1.4542E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [154.90, 4.2989E-03; 530.40, 2.7401E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.7126E07;
    out.HeatOfVap.B = 3.6975E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [154.90, 4.1476E07; 530.40, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 3.7436E04;
    out.solHeatCapacity.B = 4.1437E02;
    out.solHeatCapacity.C = 1.3235E00;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [80.00, 7.9068E04; 154.90, 1.3342E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.7401E05;
    out.liqHeatCapacity.B = -1.0578E02;
    out.liqHeatCapacity.C = 9.0525E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [154.90, 1.7935E05; 363.20, 2.5501E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0600E05;
    out.IG_HeatCapacity.B = 4.2507E05;
    out.IG_HeatCapacity.C = 1.5770E03;
    out.IG_HeatCapacity.D = 2.7580E05;
    out.IG_HeatCapacity.E = 6.9080E02;
    out.IG_HeatCapacity.range = [200.00, 1.1915E05; 1500.00, 4.5148E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.4123E-01;
    out.SecondViralCoef.B = -2.9400E02;
    out.SecondViralCoef.C = -3.8010E07;
    out.SecondViralCoef.D = -4.3597E19;
    out.SecondViralCoef.E = 6.3133E21;
    out.SecondViralCoef.range = [265.20, -3.7141E00; 1500.00, 3.3966E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.2224E01;
    out.liqViscosity.B = 1.0213E03;
    out.liqViscosity.C = 1.5187E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [154.90, 7.7129E-03; 363.20, 2.0000E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.0132E-06;
    out.vapViscosity.B = 4.5610E-01;
    out.vapViscosity.C = 3.5980E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [154.90, 3.0414E-06; 1000.00, 1.7399E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.9690E-01;
    out.liqThermalConductivity.B = -2.6920E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [154.90, 1.5520E-01; 420.00, 8.3836E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 5.9926E-05;
    out.vapThermalConductivity.B = 1.0785E00;
    out.vapThermalConductivity.C = -8.5924E01;
    out.vapThermalConductivity.D = 1.5868E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [363.20, 1.7582E-02; 1000.0, 9.6076E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0288E-02;
    out.SurfaceTension.B = 1.1923E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [154.90, 3.3314E-02; 530.40, 0];
    out.SurfaceTension.eqn = 106;    
    
end


