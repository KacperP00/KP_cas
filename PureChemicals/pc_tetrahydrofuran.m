function out = pc_tetrahydrofuran

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'tetrahydrofuran','hydrofuran','tetramethylene oxide'}; 
    out.CAS = '109-99-9';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 8;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 72.107;    % kg/kmol    
    out.Tcrit = 540.15;              % K    
    out.pcrit = 5.1678E6;            % Pa 
    out.vcrit = 0.22393;             % m^3/kmol 
    out.CritCompressFactor = 0.259;  % -
    
%% 2 
    out.MeltingPoint = 164.65;             % K
    out.TripletPointT = 164.65;            % K
    out.TripletPointp = 1.2753E-01;        % Pa
    out.NormalBoilingPoint = 338.00;       % K
    out.liqMolarVolume  = 0.081942;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.8418E08;            % J/kmol
    out.IG_GibbsOfFormation = -7.9680E07;           % J/kmol
    out.IG_AbsoluteEntropy = 2.9729E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 8.5400E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.3250E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2264;            % -
    out.RadiusOfGyration = 2.6940E-10;      % m
    out.SolubilityParameter = 1.9129E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.4371E-30;          % C^m
    out.VanDerWaalsVolume = 0.04462;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 6.0000E08;        % m^2/kmol
    out.RefractiveIndex = 1.4050;           % -
    out.FlashPoint = 258.11;                % K
    out.FlammabilityLimits = [2.0 11.8];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 594.26;          % K
    
%% 6 
    out.solDensity.A = 1.6400E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [164.65, 1.6400E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.2543E00;
    out.liqDensity.B = 2.8084E-01;
    out.liqDensity.C = 5.4015E02;
    out.liqDensity.D = 2.9120E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [164.65, 1.3998E01; 540.15, 4.4662E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  8.5184E01;
    out.vapPressure.B = -6.3373E03;
    out.vapPressure.C = -9.5958E00;
    out.vapPressure.D =  8.1661E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [164.65, 1.2753E-01; 540.15, 5.1720E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.4430E07;
    out.HeatOfVap.B = 3.9080E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [164.65, 3.8545E07; 540.15, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.8200E04;
    out.solHeatCapacity.B = 1.4470E03;
    out.solHeatCapacity.C = -9.4700E00;
    out.solHeatCapacity.D = 2.6620E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [20.00, 7.1650E03; 164.65, 8.2141E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.7165E05;
    out.liqHeatCapacity.B = -7.9835E02;
    out.liqHeatCapacity.C = 2.8812E00;
    out.liqHeatCapacity.D = -2.4813E-03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [164.65, 1.0723E05; 320.00, 1.2991E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 4.0430E04;
    out.IG_HeatCapacity.B = 2.5410E05;
    out.IG_HeatCapacity.C = 1.3866E03;
    out.IG_HeatCapacity.D = 1.4200E05;
    out.IG_HeatCapacity.E = -6.3880E02;
    out.IG_HeatCapacity.range = [50.00, 4.0430E04; 1500.00, 2.5462E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.2580E-01;
    out.SecondViralCoef.B = -1.4320E02;
    out.SecondViralCoef.C = -2.0790E07;
    out.SecondViralCoef.D = 6.3140E19;
    out.SecondViralCoef.E = -1.9480E22;
    out.SecondViralCoef.range = [304.58, -1.0912E00; 2701.00, 7.1728E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.0335E01;
    out.liqViscosity.B = 8.8360E02;
    out.liqViscosity.C = -5.2650E-02;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [164.65, 5.3152E-03; 350.00, 2.9787E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.7780E-07;
    out.vapViscosity.B = 6.5330E-01;
    out.vapViscosity.C = 2.7101E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [338.00, 9.4123E-06; 998.00, 2.7056E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.9428E-01;
    out.liqThermalConductivity.B = -2.4900E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [164.65, 1.5328E-01; 338.00, 1.1012E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 9.5521E-06;
    out.vapThermalConductivity.B = 1.4561E00;
    out.vapThermalConductivity.C = 6.6222E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [338.00, 1.5534E-02; 998.00, 1.3369E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.7130E-02;
    out.SurfaceTension.B = 1.2310E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [164.65, 4.2908E-02; 540.15, 0];
    out.SurfaceTension.eqn = 106;    
    
end


