function out = pc_nhexane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'nhexane','n-hexane','hexyl hydride'}; 
    out.CAS = '110-54-3';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 86.177;    % kg/kmol    
    out.Tcrit = 507.60;              % K    
    out.pcrit = 3.0250E6;            % Pa 
    out.vcrit = 0.37100;             % m^3/kmol 
    out.CritCompressFactor = 0.266;  % -
    
%% 2 
    out.MeltingPoint = 177.83;             % K
    out.TripletPointT = 177.83;            % K
    out.TripletPointp = 9.0170E-01;        % Pa
    out.NormalBoilingPoint = 341.88;       % K
    out.liqMolarVolume  = 0.131362;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.6694E08;            % J/kmol
    out.IG_GibbsOfFormation = -6.6340E04;           % J/kmol
    out.IG_AbsoluteEntropy = 3.8874E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 1.3080E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.8551E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3013;            % -
    out.RadiusOfGyration = 3.7690E-10;      % m
    out.SolubilityParameter = 1.4900E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.06826;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 9.6400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3723;           % -
    out.FlashPoint = 251.50;                % K
    out.FlammabilityLimits = [1.1 7.7];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 498.00;          % K
    
%% 6 
    out.solDensity.A = 1.0657E01;
    out.solDensity.B = -4.6530E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [133.15, 1.0037E01; 177.83, 9.8296E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 7.0824E-01;
    out.liqDensity.B = 2.6411E-01;
    out.liqDensity.C = 5.0760E02;
    out.liqDensity.D = 2.7537E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [177.83, 8.7470E00; 507.60, 2.6816E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.0465E02;
    out.vapPressure.B = -6.9955E03;
    out.vapPressure.C = -1.2702E01;
    out.vapPressure.D =  1.2381E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [177.83, 9.0169E-01; 507.60, 3.0449E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.4544E07;
    out.HeatOfVap.B = 3.9002E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [177.83, 3.7647E07; 507.60, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.9370E04;
    out.solHeatCapacity.B = 1.6560E03;
    out.solHeatCapacity.C = -8.3720E00;
    out.solHeatCapacity.D = 2.0200E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [20.00, 1.0563E04; 177.86, 1.2398E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.7212E05;
    out.liqHeatCapacity.B = -1.8378E02;
    out.liqHeatCapacity.C = 8.8734E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [177.83, 1.6750E05; 460.00, 2.7534E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0440E05;
    out.IG_HeatCapacity.B = 3.5230E05;
    out.IG_HeatCapacity.C = 1.6946E03;
    out.IG_HeatCapacity.D = 2.3690E05;
    out.IG_HeatCapacity.E = 7.6160E02;
    out.IG_HeatCapacity.range = [200.00, 1.1117E05; 1500.00, 3.8620E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.9700E-01;
    out.SecondViralCoef.B = -2.0000E02;
    out.SecondViralCoef.C = -3.7100E07;
    out.SecondViralCoef.D = 1.6700E19;
    out.SecondViralCoef.E = -6.1800E21;
    out.SecondViralCoef.range = [253.70, -3.3097E00; 1500.00, 5.2675E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.0715E01;
    out.liqViscosity.B = 1.2075E03;
    out.liqViscosity.C = 1.4993E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [177.83, 2.1182E-03; 343.15, 2.1541E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.7514E-07;
    out.vapViscosity.B = 7.0737E-01;
    out.vapViscosity.C = 1.5714E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [177.83, 3.6305E-06; 1000.00, 2.0050E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2492E-01;
    out.liqThermalConductivity.B = -3.5300E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [177.83, 1.6209E-01; 370.00, 9.4199E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -6.5050E02;
    out.vapThermalConductivity.B = 8.0530E-01;
    out.vapThermalConductivity.C = -1.4121E09;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [339.09, 1.7036E-02; 1000.00, 1.2003E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.5003E-02;
    out.SurfaceTension.B = 1.2674E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [177.83, 3.1841E-02; 507.60, 0];
    out.SurfaceTension.eqn = 106;    
    
end


