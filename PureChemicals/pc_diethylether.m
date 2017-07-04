function out = pc_diethylether

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'diethylether','aether','ether','ethyl ether'}; 
    out.CAS = '60-29-7';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 74.123;   % kg/kmol    
    out.Tcrit = 466.70;             % K    
    out.pcrit = 3.6400E6;           % Pa 
    out.vcrit = 0.28000;            % m^3/kmol 
    out.CritCompressFactor = 0.263; % -
    
%% 2 
    out.MeltingPoint = 156.85;              % K
    out.TripletPointT = 156.85;             % K
    out.TripletPointp = 3.9545E-01;         % Pa
    out.NormalBoilingPoint = 307.58;        % K
    out.liqMolarVolume = 0.104692;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.5210E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.2210E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.4230E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  7.1900E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.2032E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2811;            % -
    out.RadiusOfGyration = 3.1770E-10;      % m
    out.SolubilityParameter = 1.5420E04;    % (J/m^3)^.5 
    out.DipoleMoment = 3.8400E-30;          % C^m
    out.VanDerWaalsVolume = 0.05150;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.5400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3495;           % -
    out.FlashPoint = 228.15;                % K
    out.FlammabilityLimits = [1.9 48.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp =433.15;           % K
    
%% 6 
    out.solDensity.A = 1.3400E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [156.85, 1.3400E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.5540E-01;
    out.liqDensity.B = 2.6847E-01;
    out.liqDensity.C = 4.6670E02;
    out.liqDensity.D = 2.8140E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [156.85, 1.1487E01; 466.70, 3.5587E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A = 1.3690E02;
    out.vapPressure.B = -6.9543E03;
    out.vapPressure.C = -1.9254E01;
    out.vapPressure.D =  2.4508E-02;
    out.vapPressure.E =  1.0000E00;
    out.vapPressure.range = [156.85, 3.9545E-01; 466.70, 3.6412E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.0600E07;
    out.HeatOfVap.B = 3.8680E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [156.85, 3.4651E07; 466.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -2.3760E04;
    out.solHeatCapacity.B = 1.8611E03;
    out.solHeatCapacity.C = -1.2638E01;
    out.solHeatCapacity.D = 3.7650E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [30.00, 2.1715E04; 156.85 , 1.0252E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 4.4400E04;
    out.liqHeatCapacity.B = 1.3010E03;
    out.liqHeatCapacity.C = -5.5000E00;
    out.liqHeatCapacity.D = 8.7630E-03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [156.92, 1.4698E05; 460.00, 3.3202E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 8.6210E04;
    out.IG_HeatCapacity.B = 2.5510E05;
    out.IG_HeatCapacity.C = 1.5413E03;
    out.IG_HeatCapacity.D = 1.4370E05;
    out.IG_HeatCapacity.E = -6.8890E02;
    out.IG_HeatCapacity.range = [200.00, 9.3157E04; 1500.00, 2.9244E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.6330E-01;
    out.SecondViralCoef.B = -1.3680E02;
    out.SecondViralCoef.C = -2.0080E07;
    out.SecondViralCoef.D = -1.7500E19;
    out.SecondViralCoef.E = 2.6570E21;
    out.SecondViralCoef.range = [233.35,-2.6987E00; 1500.00, 6.6150E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = 1.0197E01;
    out.liqViscosity.B = -6.3800E01;
    out.liqViscosity.C = -3.2260E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [200.00,7.3593E-00; 373.15, 1.1412E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.9480E-06;
    out.vapViscosity.B = 4.1000E-01;
    out.vapViscosity.C = 4.9580E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [156.85, 3.7200E-06; 1000.00, 2.2116E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.4950E-01;
    out.liqThermalConductivity.B = -4.0700E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [156.85, 1.8566E-01; 433.15, 7.3208E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -4.4894E-03;
    out.vapThermalConductivity.B = 6.1550E-01;
    out.vapThermalConductivity.C = -3.2663E03;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [200.00, 7.6364E-03; 600.00, 5.1806E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.7928E-02;
    out.SurfaceTension.B = 1.2376E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [156.85, 3.4893E-02; 466.70, 0];
    out.SurfaceTension.eqn = 106;    
    
end


