function out = pc_1undecanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1undecanol','1-undecanol','1-Undecanol','1Undecanol','n-Undecanol','n-undecanol'}; 
    out.CAS = '112-42-5';
    % components
    out.ChemicalFormular.C = 11;         
    out.ChemicalFormular.H = 24;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 172.311;   % kg/kmol    
    out.Tcrit = 703.60;             % K    
    out.pcrit = 2.1470E06;           % Pa 
    out.vcrit = 0.71800;            % m^3/kmol 
    out.CritCompressFactor = 0.264; % -
    
%% 2 
    out.MeltingPoint = 288.45;              % K
    out.TripletPointT = 288.45;             % K
    out.TripletPointp = 1.2889E-01;          % Pa
    out.NormalBoilingPoint = 518.15;        % K
    out.liqMolarVolume = 0.207672;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -4.1840E08;             % J/kmol
    out.IG_GibbsOfFormation = -9.1530E07;            % J/kmol
    out.IG_AbsoluteEntropy =  6.3620E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.7110E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -6.7261E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.2220E08;              % J/kmol
    out.StandardHeatOfFormation = -5.0480E08;               % J/kmol
    out.StandardAbsoluteEntropy = 4.4920E05;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.6247;            % -
    out.RadiusOfGyration = 5.8080E-10;      % m
    out.SolubilityParameter = 1.9920E04;     % (J/m^3)^.5 
    out.DipoleMoment = 5.5700E-30;          % C^m
    out.VanDerWaalsVolume = 0.12401;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.7080E09;        % m^2/kmol
    out.RefractiveIndex = 1.4386;           % -
    out.FlashPoint = 383.15;                % K
    out.FlammabilityLimits = [0.7 5.2];    % vol-%
    out.FlammabilityLimitTemps = [387 428]; % K
    out.AutoignitionTemp = 538.00;          % K
    
%% 6 
    out.solDensity.A = 5.6800E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [288.45, 5.6800E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.3013E-01;
    out.liqDensity.B = 2.3703E-01;
    out.liqDensity.C = 7.0360E02;
    out.liqDensity.D = 2.7001E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [288.45, 4.8537E00; 703.60, 1.3928E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  2.2521E02;
    out.vapPressure.B = -1.8532E04;
    out.vapPressure.C = -2.8939E01;
    out.vapPressure.D =  1.0981E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [288.45, 1.2889E-01; 703.60, 2.1419E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 1.2015E08;
    out.HeatOfVap.B = 6.3047E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [288.46, 8.6152E07; 703.60, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 3.3000E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [289.05, 3.3000E05; 289.05, 3.3000E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.2945E05;
    out.liqHeatCapacity.B = -3.0395E03;
    out.liqHeatCapacity.C = 2.7927E01;
    out.liqHeatCapacity.D = -6.1847E-02;
    out.liqHeatCapacity.E = 4.3042E-05;
    out.liqHeatCapacity.range = [289.05, 3.9103E05; 520.00, 5.5125E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.8590E05;
    out.IG_HeatCapacity.B = 5.8690E05;
    out.IG_HeatCapacity.C = 1.5718E03;
    out.IG_HeatCapacity.D = 4.3260E05;
    out.IG_HeatCapacity.E = 7.2270E02;
    out.IG_HeatCapacity.range = [200.00, 2.0232E05; 1500.00, 6.7834E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 5.5925E-01;
    out.SecondViralCoef.B = -7.1346E02;
    out.SecondViralCoef.C = -1.2861E08;
    out.SecondViralCoef.D = 5.6156E21;
    out.SecondViralCoef.E = -2.1970E24;
    out.SecondViralCoef.range = [351.80, -7.1054E00; 1500.00, 4.5665E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -7.3477E01;
    out.liqViscosity.B = 5.9928E03;
    out.liqViscosity.C = 8.6131E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [288.45, 2.0000E-02; 590.15, 2.3230E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 5.9537E-08;
    out.vapViscosity.B = 8.1842E-01;
    out.vapViscosity.C = 9.0245E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [288.45, 4.6767E-06; 1000.00, 1.5578E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1211E-01;
    out.liqThermalConductivity.B = -2.1815E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [288.45, 1.4918E-01; 561.20, 8.9684E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.4988E03;
    out.vapThermalConductivity.B = 9.5209E-01;
    out.vapThermalConductivity.C = 2.0167E10;
    out.vapThermalConductivity.D = -2.6800E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [518.15, 2.4658E-02; 1000.0, 8.8994E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.3432E-02;
    out.SurfaceTension.B = 1.1218E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [288.45, 2.9565E-02; 703.60, 0];
    out.SurfaceTension.eqn = 106;    
    
end


