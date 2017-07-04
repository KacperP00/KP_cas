function out = pc_nbutane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'nbutane','n-butane','methylethylmethane'}; 
    out.CAS = '106-97-8';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 58.123;    % kg/kmol    
    out.Tcrit = 425.12;              % K    
    out.pcrit = 3.7960E6;            % Pa 
    out.vcrit = 0.25500;             % m^3/kmol 
    out.CritCompressFactor =0.274;   % -
    
%% 2 
    out.MeltingPoint = 134.86;             % K
    out.TripletPointT = 134.86;            % K
    out.TripletPointp = 6.7358E-01;        % Pa
    out.NormalBoilingPoint = 272.65;       % K
    out.liqMolarVolume  = 0.096484;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.2579E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.6700E07;           % J/kmol
    out.IG_AbsoluteEntropy =  3.0991E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  4.6610E06;         % J/kmol
    out.StandardNetHeatOfCombustion = - 2.6573E9;   % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2002;            % -
    out.RadiusOfGyration = 2.8860E-10;      % m
    out.SolubilityParameter = 1.3700E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.04780;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 6.9400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3292;           % -
    out.FlashPoint = NaN;                   % K
    out.FlammabilityLimits = [1.5 9.0];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 561.00;          % K
    
%% 6 
    out.solDensity.A = 1.5194E01;
    out.solDensity.B = -6.8446E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [23.15, 1.5036E01; 134.86, 1.4271E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.0677E00;
    out.liqDensity.B = 2.7188E-01;
    out.liqDensity.C = 4.2512E02;
    out.liqDensity.D = 2.8688E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [134.86, 1.2620E01; 425.12, 3.9271E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  6.6343E01;
    out.vapPressure.B = -4.3632E03;
    out.vapPressure.C = -7.0460E00;
    out.vapPressure.D =  9.4509E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [134.86, 6.7441E-01; 425.12, 3.7699E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 3.6238E07;
    out.HeatOfVap.B = 8.3370E-01;
    out.HeatOfVap.C = -8.2274E-01;
    out.HeatOfVap.D = 3.9613E-01;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [134.86, 2.8684E07; 425.12, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.8858E04;
    out.solHeatCapacity.B = 1.5065E03;
    out.solHeatCapacity.C = -9.8874E00;
    out.solHeatCapacity.D = 3.2073E-02;
    out.solHeatCapacity.E = 1.6746E-05;
    out.solHeatCapacity.range = [22.35, 1.0235E04; 103.54, 6.8653E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 6.4730E01;
    out.liqHeatCapacity.B = 1.6184E05;
    out.liqHeatCapacity.C = 9.8341E02;
    out.liqHeatCapacity.D = -1.4315E03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [134.86, 1.1380E-05; 420.00, 5.0822E05];
    out.liqHeatCapacity.eqn = 114;    
    
    out.IG_HeatCapacity.A = 7.1340E04;
    out.IG_HeatCapacity.B = 2.4300E05;
    out.IG_HeatCapacity.C = 1.6300E03;
    out.IG_HeatCapacity.D = 1.5033E05;
    out.IG_HeatCapacity.E = 7.3042E02;
    out.IG_HeatCapacity.range = [200.00, 7.6733E04; 1500.00, 2.6602E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.4900E-01;
    out.SecondViralCoef.B = -1.4300E02;
    out.SecondViralCoef.C = -1.0010E07;
    out.SecondViralCoef.D = -3.0300E18;
    out.SecondViralCoef.E = 4.0300E20;
    out.SecondViralCoef.range = [212.58,-1.8377E00; 1500.00, 5.0701E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -7.2471E00;
    out.liqViscosity.B = 5.3482E02;
    out.liqViscosity.C = -5.7469E-01;
    out.liqViscosity.D = -4.6625E-27;
    out.liqViscosity.E = 1E01;
    out.liqViscosity.range = [134.86, 2.2434E-03; 420.00, 3.5664E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.2980E-07;
    out.vapViscosity.B = 6.9442E-01;
    out.vapViscosity.C = 2.2766E02;
    out.vapViscosity.D = -1.4610E04;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [134.86, 3.6738E-06; 1000.00, 2.2947E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.7349E-01;
    out.liqThermalConductivity.B = -7.1267E-04;
    out.liqThermalConductivity.C = 5.1555E-07;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [134.86, 1.8676E-01; 400.00, 7.0910E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 5.1094E-02;
    out.vapThermalConductivity.B = 4.5253E-01;
    out.vapThermalConductivity.C = 5.4555E03;
    out.vapThermalConductivity.D = 1.9798E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [272.65, 1.3570E-02; 1000.00, 1.3799E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1960E-02;
    out.SurfaceTension.B = 1.2181E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [134.86, 3.2644E-02; 425.12, 0];
    out.SurfaceTension.eqn = 106;    
    
end


