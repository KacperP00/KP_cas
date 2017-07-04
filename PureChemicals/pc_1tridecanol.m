function out = pc_1tridecanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1tridecanol','n-tridecanol','n-Tridecanol','1-Tridecanol'}; 
    out.CAS = '112-70-9';
    % components
    out.ChemicalFormular.C = 13;         
    out.ChemicalFormular.H = 28;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 200.365;  % kg/kmol    
    out.Tcrit = 734.00;             % K    
    out.pcrit = 1.9350E06;          % Pa 
    out.vcrit = 0.74900;            % m^3/kmol 
    out.CritCompressFactor = 0.237; % -
    
%% 2 
    out.MeltingPoint = 303.75;              % K
    out.TripletPointT = 303.75;             % K
    out.TripletPointp = 9.4584E-02;          % Pa
    out.NormalBoilingPoint = 553.60;        % K
    out.liqMolarVolume = 0.241531;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -4.6040E08;            % J/kmol
    out.IG_GibbsOfFormation = -7.5660E07;            % J/kmol
    out.IG_AbsoluteEntropy =  7.1470E05;            % J/kmolK
    out.HeatFusionAtMeltPoint =  2.3300E07;         % J/kmol
    out.StandardNetHeatOfCombustion = -7.9014E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.7124;            % -
    out.RadiusOfGyration = 6.4170E-10;      % m
    out.SolubilityParameter = 1.9760E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.500E-30;           % C^m
    out.VanDerWaalsVolume = 0.14447;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.9780E09;        % m^2/kmol
    out.RefractiveIndex = 1.4433;           % -
    out.FlashPoint = 394.15;                % K
    out.FlammabilityLimits = [0.6 5.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 550.00;          % K
    
%% 6 
    out.solDensity.A = 4.7913E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [294.00, 4.7213E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.5093E-01;
    out.liqDensity.B = 2.6286E-01;
    out.liqDensity.C = 7.3400E02;
    out.liqDensity.D = 3.1074E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [303.75, 4.1402E00; 734.00, 1.3320E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  2.6599E02;
    out.vapPressure.B = -2.1388E04;
    out.vapPressure.C = -3.4855E01;
    out.vapPressure.D =  1.4125E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [303.75, 9.4584E-02; 734.00, 1.9240E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 2.3170E08;
    out.HeatOfVap.B = 2.3543E00;
    out.HeatOfVap.C = -1.7420E00;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [303.75, 9.6831E07; 734.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 3.8300E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [293.15, 3.8300E08; 293.15, 3.8300E08];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.0247E06;
    out.liqHeatCapacity.B = -3.2237E04;
    out.liqHeatCapacity.C = 1.3859E02;
    out.liqHeatCapacity.D = -2.4455E-01;
    out.liqHeatCapacity.E = 1.5407E-04;
    out.liqHeatCapacity.range = [303.75, 4.7757E05; 553.60, 6.3232E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.1780E05;
    out.IG_HeatCapacity.B = 6.8080E05;
    out.IG_HeatCapacity.C = 1.5702E03;
    out.IG_HeatCapacity.D = 5.0750E05;
    out.IG_HeatCapacity.E = 7.2440E02;
    out.IG_HeatCapacity.range = [150.00, 2.2082E05; 1500.00, 7.2076E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.4000E-02;
    out.SecondViralCoef.B = -5.5700E01;
    out.SecondViralCoef.C = -1.2900E07;
    out.SecondViralCoef.D = -6.4900E19;
    out.SecondViralCoef.E = -2.5800E22;
    out.SecondViralCoef.range = [367.00, -8.2103E00; 1500.00, 5.3664E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.6661E01;
    out.liqViscosity.B = 3.7541E03;
    out.liqViscosity.C = 1.6684E-02;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [303.75, 1.4899E-02; 553.60, 5.6891E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.2436E-08;
    out.vapViscosity.B = 8.8190E-01;
    out.vapViscosity.C = -4.0200E00;
    out.vapViscosity.D = 2.7480E04;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [303.75, 4.3038E-06; 1000.00, 1.4531E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1542E-01;
    out.liqThermalConductivity.B = -2.1522E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [303.75, 1.5005E-01; 553.60, 9.6274E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 3.4379E03;
    out.vapThermalConductivity.B = 9.5884E-01;
    out.vapThermalConductivity.C = 3.0868E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [553.60, 2.6319E-02; 1000.0, 8.3812E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.9791E-02;
    out.SurfaceTension.B = 1.2646E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [303.75, 3.0428E-02; 734.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


