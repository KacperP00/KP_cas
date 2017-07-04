function out = pc_acetone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'acetone','Acetone','Allylic Alcohol','Dimethylketal','3-Hydroxypropene','Propanone','beta-Ketopropane'}; 
    out.CAS = '67-64-1';
    % components
    out.ChemicalFormular.C = 3;         
    out.ChemicalFormular.H = 6;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 58.080;   % kg/kmol    
    out.Tcrit = 508.20;             % K    
    out.pcrit = 4.7015E6;          % Pa 
    out.vcrit = 0.20900;            % m^3/kmol 
    out.CritCompressFactor = 0.233; % -
    
%% 2 
    out.MeltingPoint = 178.45;              % K
    out.TripletPointT = 178.45;             % K
    out.TripletPointp = 2.5938E00;          % Pa
    out.NormalBoilingPoint = 329.44;        % K
    out.liqMolarVolume = 0.073931;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.1715E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.5272E08;            % J/kmol
    out.IG_AbsoluteEntropy =  2.9535E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  5.6912E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.6592E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3064;            % -
    out.RadiusOfGyration = 2.7460E-10;      % m
    out.SolubilityParameter = 1.9774E04;    % (J/m^3)^.5 
    out.DipoleMoment = 9.6066E-30;          % C^m
    out.VanDerWaalsVolume = 0.03904;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 5.8400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3560;           % -
    out.FlashPoint = 255.37;                % K
    out.FlammabilityLimits = [2.6 12.8];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 810.93;          % K
    
%% 6 
    out.solDensity.A = 1.6677E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [174.15, 1.6677E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.2298E00;
    out.liqDensity.B = 2.5760E-01;
    out.liqDensity.C = 5.0820E02;
    out.liqDensity.D = 2.9903E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [178.45, 1.5721E01; 508.20, 4.7741E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.0720E01;
    out.vapPressure.B = -5.6850E03;
    out.vapPressure.C = -7.3510E00;
    out.vapPressure.D =  6.300E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [178.45, 2.5938E00; 508.20, 4.6733E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.9170E07;
    out.HeatOfVap.B = 1.0360E00;
    out.HeatOfVap.C = -1.2940E00;
    out.HeatOfVap.D = 6.7200E-01;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [178.45, 3.6887E07; 508.2, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.2000E04;
    out.solHeatCapacity.B = 1.2760E03;
    out.solHeatCapacity.C = -6.1700E00;
    out.solHeatCapacity.D = 1.3200E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [30.00, 2.1083E04; 178.45, 9.4233E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.3560E05;
    out.liqHeatCapacity.B = -1.7700E02;
    out.liqHeatCapacity.C = 2.8370E-01;
    out.liqHeatCapacity.D = 6.8900E-04;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [178.45,1.169E05; 329.44, 1.3271E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 4.8120E04;
    out.IG_HeatCapacity.B = 1.6440E05;
    out.IG_HeatCapacity.C = 1.2500E03;
    out.IG_HeatCapacity.D = 7.1700E04;
    out.IG_HeatCapacity.E = -5.244E02;
    out.IG_HeatCapacity.range = [100.00, 4.8340E04; 1500.00, 1.8699E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.1070E-01;
    out.SecondViralCoef.B = -9.9200E01;
    out.SecondViralCoef.C = -3.0540E07;
    out.SecondViralCoef.D = 9.8900E18;
    out.SecondViralCoef.E = -1.6454E22;
    out.SecondViralCoef.range = [254.10,-5.298E00; 2540.50, 6.9790E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.4918E01;
    out.liqViscosity.B = 1.0234E03;
    out.liqViscosity.C = 5.9610E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [223.15, 8.1834E-04; 329.44, 2.3506E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.1005E-08;
    out.vapViscosity.B = 9.7620E-01;
    out.vapViscosity.C = 2.3139E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [300.00, 7.5393E-06; 650.00, 1.6680E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.5020E-01;
    out.liqThermalConductivity.B = -2.9800E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [270.00, 1.6974E-01; 315.00, 1.5633E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -2.6800E01;
    out.vapThermalConductivity.B = 9.0980E-01;
    out.vapThermalConductivity.C = -1.2650E08;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.15, 9.5298E-03; 571.70, 3.9056E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.2200E-02;
    out.SurfaceTension.B = 1.1240E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [178.45, 3.8251E-02; 508.20, 0];
    out.SurfaceTension.eqn = 106;    
    
end


