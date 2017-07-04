function out = pc_2propanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2propanol','2-propanol','dimethyl carbinol','isohol'}; 
    out.CAS = '67-63-0';
    % components
    out.ChemicalFormular.C = 3;         
    out.ChemicalFormular.H = 8;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 60.096;   % kg/kmol    
    out.Tcrit = 508.31;             % K    
    out.pcrit = 4.7643E6;           % Pa 
    out.vcrit = 0.22013;            % m^3/kmol 
    out.CritCompressFactor = 0.248; % -
    
%% 2 
    out.MeltingPoint = 185.28;              % K
    out.TripletPointT = 185.28;             % K
    out.TripletPointp = 3.2000E-2;          % Pa
    out.NormalBoilingPoint = 355.41;        % K
    out.liqMolarVolume = 0.076784;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.7242E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.7339E8;            % J/kmol
    out.IG_AbsoluteEntropy =  3.0991E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  5.4099E6;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.8300E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.6689;            % -
    out.RadiusOfGyration = 2.8070E-10;      % m
    out.SolubilityParameter = 2.3575E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.5372E-30;          % C^m
    out.VanDerWaalsVolume = 0.04216;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 6.2700E08;        % m^2/kmol
    out.RefractiveIndex = 1.3752;           % -
    out.FlashPoint = 284.82;                % K
    out.FlammabilityLimits = [2 12.0];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 672.04;          % K
    
%% 6 
    out.solDensity.A = 1.71E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [185.28, 1.7100E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.18E00;
    out.liqDensity.B = 2.6475E-01;
    out.liqDensity.C = 5.0831E02;
    out.liqDensity.D = 2.43E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [185.28, 1.4656E01; 508.31, 4.4570E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  9.2935E01;
    out.vapPressure.B = -8.1771E03;
    out.vapPressure.C = -1.0031E01;
    out.vapPressure.D =  3.9988E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [185.28, 3.2000E-02; 508.31, 4.7673E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.6980E07;
    out.HeatOfVap.B = 8.7E-02;
    out.HeatOfVap.C = 3.0070E-01;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [185.28, 5.2121E07; 508.31, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -8.5E03;
    out.solHeatCapacity.B = 8.7600E02;
    out.solHeatCapacity.C = -8.100E-01;
    out.solHeatCapacity.D = -3.025E-02;
    out.solHeatCapacity.E = 1.4130E-04;
    out.solHeatCapacity.range = [12.00, 1.8460E03; 176.66, 9.1821E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 4.6640E05;
    out.liqHeatCapacity.B = -4.1086E03;
    out.liqHeatCapacity.C = 1.4506E01;
    out.liqHeatCapacity.D = -1.4126E-02;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [188.32, 1.1277E05; 463.15, 2.7174E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 4.7460E04;
    out.IG_HeatCapacity.B = 1.9350E05;
    out.IG_HeatCapacity.C = 1.1240E03;
    out.IG_HeatCapacity.D = 9.3800E04;
    out.IG_HeatCapacity.E = 4.6000E02;
    out.IG_HeatCapacity.range = [100.00, 4.8262E04; 1500.00, 2.1652E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.0200E-02;
    out.SecondViralCoef.B = -6.2900E00;
    out.SecondViralCoef.C = -4.3100E07;
    out.SecondViralCoef.D = 1.9860E20;
    out.SecondViralCoef.E = -8.6100E22;
    out.SecondViralCoef.range = [279.58, -4.8942E00; 2541.60, 2.5100E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -8.2300E00;
    out.liqViscosity.B = 2.2822E03;
    out.liqViscosity.C = -9.8495E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [187.35, 3.0037E-01; 354.45, 5.1389E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.9930E-07;
    out.vapViscosity.B = 7.2330E-01;
    out.vapViscosity.C = 1.7800E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [200.00, 4.8683E-06; 1000.00, 2.5019E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0290E-01;
    out.liqThermalConductivity.B = -2.2780E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [185.28, 1.6069E-01; 410.00, 1.0950E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -8.0642E01;
    out.vapThermalConductivity.B = -1.4549E00;
    out.vapThermalConductivity.C = -6.0442E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [355.41, 2.2388E-02; 450.00, 3.2427E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 3.818E-02;
    out.SurfaceTension.B = -3.8180E-05;
    out.SurfaceTension.C = -6.5100E-08;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [273.15, 2.2894E-02; 355.41, 1.6387E-02];
    out.SurfaceTension.eqn = 100;    
    
end


