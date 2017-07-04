function out = pc_isobutane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'isobutane','iso-butane','i-butane'}; 
    out.CAS = '75-28-5';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 58.123;   % kg/kmol    
    out.Tcrit = 407.80;             % K    
    out.pcrit = 3.6400E6;           % Pa 
    out.vcrit = 0.25900;            % m^3/kmol 
    out.CritCompressFactor = 0.278; % -
    
%% 2 
    out.MeltingPoint = 113.54;              % K
    out.TripletPointT = 113.54;             % K
    out.TripletPointp = 1.2062E-2;          % Pa
    out.NormalBoilingPoint = 261.43;        % K
    out.liqMolarVolume = 0.105350;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.3499E8;             % J/kmol
    out.IG_GibbsOfFormation = -2.1440E7;            % J/kmol
    out.IG_AbsoluteEntropy =  2.9550E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  4.5400E6;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.6490E9;    % J/kmol
    out.StandardGibbsOfFormation = -2.0760E07;      % J/kmol
    out.StandardHeatOfFormation = - 1.3418E08;      % J/kmol
    out.StandardAbsoluteEntropy = 2.9539E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.1835;            % -
    out.RadiusOfGyration = 2.9480E-10;      % m
    out.SolubilityParameter = 1.2560E4;     % (J/m^3)^.5 
    out.DipoleMoment = 4.4000E-31;          % C^m
    out.VanDerWaalsVolume = 0.04779;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 6.9300E08;        % m^2/kmol
    out.RefractiveIndex = 1.3175;           % -
    out.FlashPoint = NaN;                   % K
    out.FlammabilityLimits = [1.8 8.4];     % vol-%
    out.FlammabilityLimitTemps = [191 212]; % K
    out.AutoignitionTemp = 733.15;          % K
    
%% 6 
    out.solDensity.A = 1.4700E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [113.54, 1.4700E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.0631E00;
    out.liqDensity.B = 2.7506E-01;
    out.liqDensity.C = 4.0780E02;
    out.liqDensity.D = 2.7580E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [113.54, 1.2574E01; 407.80, 3.8650E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.0843E02;
    out.vapPressure.B = -5.0399E03;
    out.vapPressure.C = -1.5012E01;
    out.vapPressure.D =  2.2725E-02;
    out.vapPressure.E =  1.0000E00;
    out.vapPressure.range = [113.54, 1.2062E-02; 407.80, 3.6299E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 3.1880E07;
    out.HeatOfVap.B = 3.9006E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [113.54, 2.8070E07; 407.80, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 1.2546E05;
    out.solHeatCapacity.B = -2.1580E03;
    out.solHeatCapacity.C = 1.5616E01;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [77.00, 5.1881E04; 107.00, 7.3342E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.7237E05;
    out.liqHeatCapacity.B = -1.7839E03;
    out.liqHeatCapacity.C = 1.4759E01;
    out.liqHeatCapacity.D = -4.7909E-02;
    out.liqHeatCapacity.E = 5.8050E-05;
    out.liqHeatCapacity.range = [113.54, 9.9613E04; 380.00, 2.0725E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 6.5490E04;
    out.IG_HeatCapacity.B = 2.4776E05;
    out.IG_HeatCapacity.C = 1.5870E03;
    out.IG_HeatCapacity.D = 1.5750E05;
    out.IG_HeatCapacity.E = -7.0699E02;
    out.IG_HeatCapacity.range = [200.00, 7.2181E04; 1500.00, 2.6656E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.4990E-01;
    out.SecondViralCoef.B = -1.3056E02;
    out.SecondViralCoef.C = -8.7980E06;
    out.SecondViralCoef.D = -9.9500E17;
    out.SecondViralCoef.E = 5.2500E19;
    out.SecondViralCoef.range = [204.70, -1.7532E00; 1500.00, 6.0253E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.3912E01;
    out.liqViscosity.B = 7.9709E02;
    out.liqViscosity.C = 4.5508E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [110.00 , 1.0716E-02; 310.95, 1.5877E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.0871E-07;
    out.vapViscosity.B = 7.8135E-01;
    out.vapViscosity.C = 7.0639E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [150.00, 3.7065E-06; 1000.00, 2.2422E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0455E-01;
    out.liqThermalConductivity.B = -3.6589E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [113.54, 1.6301E-01; 400.00, 5.8194E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 8.9772E-02;
    out.vapThermalConductivity.B = 1.8501E-01;
    out.vapThermalConductivity.C = 6.3923E02;
    out.vapThermalConductivity.D = 1.1147E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [261.43, 1.2726E-02; 1000.0, 1.1701E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1359E-02;
    out.SurfaceTension.B = 1.2532E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [113.54 , 3.4121E-02; 407.80, 0];
    out.SurfaceTension.eqn = 106;    
    
end


