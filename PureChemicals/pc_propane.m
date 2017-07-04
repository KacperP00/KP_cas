function out = pc_propane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'propane','dimethylmethane','propyl hydride'}; 
    out.CAS = '74-98-6';
    % components
    out.ChemicalFormular.C = 3;         
    out.ChemicalFormular.H = 8;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 44.096;    % kg/kmol    
    out.Tcrit = 369.83;              % K    
    out.pcrit = 4.2480E6;            % Pa 
    out.vcrit = 0.20000;             % m^3/kmol 
    out.CritCompressFactor = 0.276;  % -
    
%% 2 
    out.MeltingPoint = 85.47;              % K
    out.TripletPointT = 85.47;             % K
    out.TripletPointp = 1.6850E-04;        % Pa
    out.NormalBoilingPoint = 231.11;       % K
    out.liqMolarVolume  = 0.075700;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.0468E08;            % J/kmol
    out.IG_GibbsOfFormation = -2.4390E07;           % J/kmol
    out.IG_AbsoluteEntropy = 2.7020E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 3.5240E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.0431E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.1523;            % -
    out.RadiusOfGyration = 2.4310E-10;      % m
    out.SolubilityParameter = 1.3100E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.03757;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 5.5900E08;        % m^2/kmol
    out.RefractiveIndex = 1.2861;           % -
    out.FlashPoint = NaN;                   % K
    out.FlammabilityLimits = [2.0 9.5];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 723.00;          % K
    
%% 6 
    out.solDensity.A = 1.8861E01;
    out.solDensity.B = -2.0332E-02;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [23.15, 1.8390E01; 85.47, 1.7123E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.3757E00;
    out.liqDensity.B = 2.7453E-01;
    out.liqDensity.C = 3.6983E02;
    out.liqDensity.D = 2.9359E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [85.47, 1.6583E01; 369.83, 5.0111E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.9078E01;
    out.vapPressure.B = -3.4926E03;
    out.vapPressure.C = -6.0669E00;
    out.vapPressure.D =  1.0919E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [85.47, 1.6788E-04; 369.83, 4.2135E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 2.9209E07;
    out.HeatOfVap.B = 7.8237E-01;
    out.HeatOfVap.C = -7.7319E-01;
    out.HeatOfVap.D = 3.9246E-01;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [85.47, 2.4787E07; 369.83, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.1230E04;
    out.solHeatCapacity.B = 1.0590E03;
    out.solHeatCapacity.C = -3.6000E00;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [30.00, 1.7300E04; 84.00, 5.2324E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 6.2983E01;
    out.liqHeatCapacity.B = 1.1363E05;
    out.liqHeatCapacity.C = 6.3321E02;
    out.liqHeatCapacity.D = -8.7346E02;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [85.47, 8.4879E04; 360.00, 2.6079E05];
    out.liqHeatCapacity.eqn = 114;    
    
    out.IG_HeatCapacity.A = 5.1920E04;
    out.IG_HeatCapacity.B = 1.9245E05;
    out.IG_HeatCapacity.C = 1.6265E03;
    out.IG_HeatCapacity.D = 1.1680E05;
    out.IG_HeatCapacity.E = 7.2360E02;
    out.IG_HeatCapacity.range = [200.00, 5.6322E04; 1500.00, 2.0556E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.1270E-01;
    out.SecondViralCoef.B = -9.9200E01;
    out.SecondViralCoef.C = -4.5100E06;
    out.SecondViralCoef.D = 3.0900E17;
    out.SecondViralCoef.E = -7.0500E19;
    out.SecondViralCoef.range = [184.91 ,-1.1900E00; 1500.00, 4.5230E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -6.9281E00;
    out.liqViscosity.B = 4.2076E02;
    out.liqViscosity.C = -6.3276E-01;
    out.liqViscosity.D = -1.7130E-26;
    out.liqViscosity.E = 1E01;
    out.liqViscosity.range = [85.47, 8.0682E-03; 360.00, 4.0666E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.4993E-07;
    out.vapViscosity.B = 6.8612E-01;
    out.vapViscosity.C = 1.7934E02;
    out.vapViscosity.D = -8.2546E03;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [85.47, 2.6865E-06; 1000.00, 2.4411E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.6755E-01;
    out.liqThermalConductivity.B = -6.6457E-04;
    out.liqThermalConductivity.C = 2.7740E-07;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [85.47, 2.1278E-01; 350.00, 6.8932E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -1.1200E00;
    out.vapThermalConductivity.B = 1.0972E-01;
    out.vapThermalConductivity.C = -9.8346E03;
    out.vapThermalConductivity.D = -7.5358E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [231.11, 1.1142E-02; 1000.00, 1.4599E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0920E-02;
    out.SurfaceTension.B = 1.2197E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [85.47, 3.6956E-02; 369.83, 0];
    out.SurfaceTension.eqn = 106;    
    
end


