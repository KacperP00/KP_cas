function out = pc_1octanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1octanol','1-octanol','1-Octanol','n-octanol','n-Octanol, octanol'}; 
    out.CAS = '111-87-5';
    % components
    out.ChemicalFormular.C = 8;         
    out.ChemicalFormular.H = 18;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 130.230;  % kg/kmol    
    out.Tcrit = 652.50;             % K    
    out.pcrit = 2.7770E06;          % Pa 
    out.vcrit = 0.49700;            % m^3/kmol 
    out.CritCompressFactor = 0.254; % -
    
%% 2 
    out.MeltingPoint = 257.65;              % K
    out.TripletPointT = 257.65;             % K
    out.TripletPointp = 8.5564E-02;         % Pa
    out.NormalBoilingPoint = 468.35;        % K
    out.liqMolarVolume = 0.158231;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.5550E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.1550E08;            % J/kmol
    out.IG_AbsoluteEntropy =  5.1850E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.2500E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.8983E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.3690E08;       % J/kmol
    out.StandardHeatOfFormation = -4.2650E08;        % J/kmol
    out.StandardAbsoluteEntropy = 3.5200E05;         % J/kmolK
    
%% 4 
    out.AcentricFactor = 0.5829;            % -
    out.RadiusOfGyration = 4.7870E-10;      % m
    out.SolubilityParameter = 2.1020E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.5000E-30;          % C^m
    out.VanDerWaalsVolume = 0.09332;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.3030E09;        % m^2/kmol
    out.RefractiveIndex = 1.4276;           % -
    out.FlashPoint = 354.00;                % K
    out.FlammabilityLimits = [0.9 6.4];     % vol-%
    out.FlammabilityLimitTemps = [355 392]; % K
    out.AutoignitionTemp = 555.00;          % K
    
%% 6 
    out.solDensity.A = 7.6600E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [257.65, 7.6600E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.1836E-01;
    out.liqDensity.B = 2.5771E-01;
    out.liqDensity.C = 6.5250E02;
    out.liqDensity.D = 2.7713E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [257.65, 6.5440E00; 652.50, 2.0114E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.5390E02;
    out.vapPressure.B = -1.4181E+04;
    out.vapPressure.C = -1.8251E+01;
    out.vapPressure.D =  1.2261E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [257.65, 8.5564E-02; 652.50, 2.7625E+06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 1.0730E+08;
    out.HeatOfVap.B = 6.4524E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [257.65, 7.7596E+07; 652.50, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 4.1740E+04;
    out.solHeatCapacity.B = 5.7900E+02;
    out.solHeatCapacity.C = 2.7800E-01;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [101.90, 1.0363E+05; 216.90 1.8040E+05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 5.7137E+05;
    out.liqHeatCapacity.B = -4.8490E+03;
    out.liqHeatCapacity.C = 1.9725E+01;
    out.liqHeatCapacity.D = -2.1532E-02;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [250.00, 2.5550E+05; 468.35, 4.1500E+05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.3805E+05;
    out.IG_HeatCapacity.B = 4.4590E+05;
    out.IG_HeatCapacity.C = 1.5751E+03;
    out.IG_HeatCapacity.D = 3.2016E+05;
    out.IG_HeatCapacity.E = 7.1880E+02;
    out.IG_HeatCapacity.range = [200.00, 1.5055E+05; 1500.00, 5.0965E+05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.2384E-01;
    out.SecondViralCoef.B = -3.5035E+02;
    out.SecondViralCoef.C = -1.2731E+08;
    out.SecondViralCoef.D = -3.0518E+20;
    out.SecondViralCoef.E = -1.1146E+22;
    out.SecondViralCoef.range = [326.25, -7.0600E+00; 1500.00, 5.2540E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -5.0674E+01;
    out.liqViscosity.B = 4.7251E+03;
    out.liqViscosity.C = 5.2499E+00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [280.00, 1.4749E-02; 425.00, 4.1665E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.7520E-07;
    out.vapViscosity.B = 6.9410E-01;
    out.vapViscosity.C = 2.0680E+02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [257.65, 4.5827E-06; 1000.00, 1.7547E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.3160E-01;
    out.liqThermalConductivity.B = -2.4070E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [257.65, 1.6958E-01; 570.15, 9.4365E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -3.0238E-03;
    out.vapThermalConductivity.B = 8.7450E-01;
    out.vapThermalConductivity.C = -1.3352E+04;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [468.35, 2.3796E-02; 1000.00, 1.0288E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1595E-02;
    out.SurfaceTension.B = 1.0483E+00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [257.65, 3.0355E-02; 652.50, 0.0000E+00];
    out.SurfaceTension.eqn = 106;    
    
end


