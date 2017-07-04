function out = pc_1pentanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1pentanol','1-pentanol','1-Pentanol','1Pentanol','n-Pentanol','n-pentanol'}; 
    out.CAS = '71-41-0';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 88.150;   % kg/kmol    
    out.Tcrit = 588.10;             % K    
    out.pcrit = 3.8970E06;          % Pa 
    out.vcrit = 0.32600;            % m^3/kmol 
    out.CritCompressFactor = 0.260; % -
    
%% 2 
    out.MeltingPoint = 195.56;              % K
    out.TripletPointT = 195.56;             % K
    out.TripletPointp = 3.4525E-4;          % Pa
    out.NormalBoilingPoint = 410.95;        % K
    out.liqMolarVolume = 0.108538;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.9874E08;            % J/kmol
    out.IG_GibbsOfFormation = -1.4602E08;           % J/kmol
    out.IG_AbsoluteEntropy =  4.0250E05;            % J/kmolK
    out.HeatFusionAtMeltPoint =  9.7906E06;         % J/kmol
    out.StandardNetHeatOfCombustion = -3.0605E09;   % J/kmol
    out.StandardGibbsOfFormation = -1.5610E08;      % J/kmol
    out.StandardHeatOfFormation = -3.5160E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.5890E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5731;            % -
    out.RadiusOfGyration = 3.6790E-10;      % m
    out.SolubilityParameter = 2.2520E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.6706E-30;          % C^m
    out.VanDerWaalsVolume = 0.06263;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.9800E08;        % m^2/kmol
    out.RefractiveIndex = 1.4080;           % -
    out.FlashPoint = 305.93;                % K
    out.FlammabilityLimits = [1.2 10.0];    % vol-%
    out.FlammabilityLimitTemps = [316 353]; % K
    out.AutoignitionTemp = 573.15;          % K
    
%% 6 
    out.solDensity.A = 1.1800E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [195.56, 1.1800E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.1754E-01;
    out.liqDensity.B = 2.6732E-01;
    out.liqDensity.C = 5.8810E02;
    out.liqDensity.D = 2.5348E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [195.56, 1.0061E01; 588.10, 3.0583E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.6264E02;
    out.vapPressure.B = -1.2413E04;
    out.vapPressure.C = -2.0383E01;
    out.vapPressure.D =  1.0482E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [195.56, 3.4524E-04; 588.10, 3.9062E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 8.3972E07;
    out.HeatOfVap.B = 5.3507E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [195.56, 6.7638E07; 588.10, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.8050E04;
    out.solHeatCapacity.B = 1.6164E03;
    out.solHeatCapacity.C = -9.4040E00;
    out.solHeatCapacity.D = 2.4970E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [31.56, 2.4382E04; 178.06, 1.1258E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.0120E05;
    out.liqHeatCapacity.B = -6.5130E02;
    out.liqHeatCapacity.C = 2.2750E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [200.14, 1.6198E05; 389.15, 2.9227E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 9.0600E04;
    out.IG_HeatCapacity.B = 3.0620E05;
    out.IG_HeatCapacity.C = 1.6054E03;
    out.IG_HeatCapacity.D = 2.1150E05;
    out.IG_HeatCapacity.E = -7.1797E02;
    out.IG_HeatCapacity.range = [200.00, 9.8902E04; 1500.00, 3.4133E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.1875E-01;
    out.SecondViralCoef.B = -2.1461E02;
    out.SecondViralCoef.C = -5.9275E07;
    out.SecondViralCoef.D = -1.3203E20;
    out.SecondViralCoef.E = 9.8814E21;
    out.SecondViralCoef.range = [294.05, -4.6034E00; 1500.00, 5.8109E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.2758E01;
    out.liqViscosity.B = 2.9169E03;
    out.liqViscosity.C = 1.2839E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [253.15, 1.6075E-02; 420.00, 3.1659E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.8903E-07;
    out.vapViscosity.B = 7.0310E-01;
    out.vapViscosity.C = 1.7590E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [410.95, 9.1106E-06; 1000.00, 2.0676E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0060E-01;
    out.liqThermalConductivity.B = -1.6030E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [273.15, 1.5681E-01; 353.15, 1.4399E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.8960E03;
    out.vapThermalConductivity.B = 8.9850E-01;
    out.vapThermalConductivity.C = 1.2735E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [410.95, 2.0849E-02; 990.95, 1.1087E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0660E-02;
    out.SurfaceTension.B = -8.3920E-05;
    out.SurfaceTension.C = -3.8300E-09;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [283.15, 2.6591E-02; 373.15, 1.8812E-02];
    out.SurfaceTension.eqn = 100;    
    
end


