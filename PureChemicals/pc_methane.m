function out = pc_methane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'methane','Methane'}; 
    out.CAS = '74-82-8';
    % components
    out.ChemicalFormular.C = 1;         
    out.ChemicalFormular.H = 4;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 16.043;    % kg/kmol    
    out.Tcrit = 190.56;              % K    
    out.pcrit = 4.5990E06;           % Pa 
    out.vcrit = 0.09860;             % m^3/kmol 
    out.CritCompressFactor = 0.286;  % -
    
%% 2 
    out.MeltingPoint = 90.69;              % K
    out.TripletPointT = 90.69;             % K
    out.TripletPointp = 1.1696E04;         % Pa
    out.NormalBoilingPoint = 111.66;       % K
    out.liqMolarVolume  = 0.037969;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -7.4520E07;            % J/kmol
    out.IG_GibbsOfFormation = -5.0490E07;           % J/kmol
    out.IG_AbsoluteEntropy =  1.8627E05;            % J/kmolK
    out.HeatFusionAtMeltPoint =  9.4140E05;         % J/kmol
    out.StandardNetHeatOfCombustion = -8.0262E08;   % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.0115;            % -
    out.RadiusOfGyration = 1.1180E-10;      % m
    out.SolubilityParameter = 1.1600E04;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.01705;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 2.8800E08;        % m^2/kmol
    out.RefractiveIndex = 1.0004;           % -
    out.FlashPoint = NaN;                   % K
    out.FlammabilityLimits = [5.0 15.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 810.00;          % K
    
%% 6 
    out.solDensity.A = 3.3022E01;
    out.solDensity.B = -1.5870E-02;
    out.solDensity.C = -1.5500E-04;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [23.15, 3.2572E01; 90.67, 3.0309E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 2.9214E00;
    out.liqDensity.B = 2.8976E-01;
    out.liqDensity.C = 1.9056E02;
    out.liqDensity.D = 2.8818E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [90.69, 2.8180E01; 190.56, 1.0082E01];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  3.9205E01;
    out.vapPressure.B = -1.3244E03;
    out.vapPressure.C = -3.4366E00;
    out.vapPressure.D =  3.1019E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [90.69, 1.1687E04; 190.56, 4.5897E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 1.0194E07;
    out.HeatOfVap.B = 2.6087E-01;
    out.HeatOfVap.C = -1.4694E-01;
    out.HeatOfVap.D = 2.2154E-01;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [90.69, 8.7235E06; 190.56, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -2.5500E03;
    out.solHeatCapacity.B = 1.2499E03;
    out.solHeatCapacity.C = -1.4170E01;
    out.solHeatCapacity.D = 6.6480E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [22.85, 1.9405E04; 90.67, 4.3841E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 6.5708E01;
    out.liqHeatCapacity.B = 3.8883E04;
    out.liqHeatCapacity.C = -2.5795E02;
    out.liqHeatCapacity.D = 6.1407E02;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [90.69, 5.3605E04; 190.00, 1.4978E06];
    out.liqHeatCapacity.eqn = 114;    
    
    out.IG_HeatCapacity.A = 3.3298E04;
    out.IG_HeatCapacity.B = 7.9933E04;
    out.IG_HeatCapacity.C = 2.0869E03;
    out.IG_HeatCapacity.D = 4.1602E04;
    out.IG_HeatCapacity.E = 9.9196E02;
    out.IG_HeatCapacity.range = [50.00, 3.3298E04; 1500.00, 8.8904E04];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 5.1075E-02;
    out.SecondViralCoef.B = -2.5181E01;
    out.SecondViralCoef.C = -2.5601E05;
    out.SecondViralCoef.D = 5.9777E15;
    out.SecondViralCoef.E =-5.7697E17;
    out.SecondViralCoef.range = [110.83, -3.3028E-01; 1500.00, 3.4212E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -6.1572E00;
    out.liqViscosity.B = 1.7815E02;
    out.liqViscosity.C = -9.5239E-01;
    out.liqViscosity.D = -9.0606E-24;
    out.liqViscosity.E = 1.0E01;
    out.liqViscosity.range = [90.69, 2.0634E-04; 188.00, 2.2625E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 5.2546E-07;
    out.vapViscosity.B = 5.9006E-01;
    out.vapViscosity.C = 1.0567E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [90.69, 3.4684E-06; 1000.00, 2.7996E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 4.1768E-01;
    out.liqThermalConductivity.B = -2.4528E-03;
    out.liqThermalConductivity.C = 3.5588E-06;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [90.69,  2.2451E-01; 180.00, 9.1481E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 6.3252E03;
    out.vapThermalConductivity.B = 4.3041E-01;
    out.vapThermalConductivity.C = 7.7040E08;
    out.vapThermalConductivity.D = -3.8725E10;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [90.69, 1.1625E-02; 1000.0, 1.6904E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 3.6557E-02;
    out.SurfaceTension.B = 1.1466E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [90.69, 1.7428E-02; 190.56, 0];
    out.SurfaceTension.eqn = 106;    
    
end


