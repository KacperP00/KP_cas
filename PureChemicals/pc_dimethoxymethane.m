function out = pc_dimethoxymethane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'OME1','Methylal','ome1','methylal'}; 
    out.CAS = '	109-87-5';
    % components
    out.ChemicalFormular.C = 3;         
    out.ChemicalFormular.H = 8;
    out.ChemicalFormular.O = 2;
%% 1  
    out.MolecularWeight = 76.095;   % kg/kmol    
    out.Tcrit = 480.6;             % K    
    out.pcrit = 3.9517E06;          % Pa 
    out.vcrit = 0.21300;            % m^3/kmol 
    out.CritCompressFactor = 0.211; % -
    
%% 2 
    out.MeltingPoint = 168.35;              % K
    out.TripletPointT = 168.35;             % K
    out.TripletPointp = 1.9634E00;         % Pa
    out.NormalBoilingPoint = 315.00;        % K
    out.liqMolarVolume = 0.089127;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.4820E08;             % J/kmol
    out.IG_GibbsOfFormation = -2.2631E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.3572E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  8.3303E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.7998E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2904;            % -
    out.RadiusOfGyration = 2.8910E-10;      % m
    out.SolubilityParameter = 1.7189E04;    % (J/m^3)^.5 
    out.DipoleMoment = 2.4684E-30;          % C^m
    out.VanDerWaalsVolume = 0.04497;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 6.7900E08;        % m^2/kmol
    out.RefractiveIndex = 1.3504;           % -
    out.FlashPoint = 255.37;                % K
    out.FlammabilityLimits = [1.6 17.6];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 510.35;          % K

%% 6 
    out.solDensity.A = 1.5600E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [168.35, 1.5600E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.4355E00;
    out.liqDensity.B = 3.0576E-01;
    out.liqDensity.C = 4.8060E02;
    out.liqDensity.D = 3.17755E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [168.35, 1.3194E01; 480.6, 4.6949E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.6361E01;
    out.vapPressure.B = -5.0483E03;
    out.vapPressure.C = -5.0136E00;
    out.vapPressure.D =  2.3659E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [168.35, 1.9634E00; 480.6, 3.9500E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.1300E07;
    out.HeatOfVap.B = 1.7170E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [168.35, 3.5183E07; 480.6, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.6050E04;
    out.solHeatCapacity.B = 1.3984E03;
    out.solHeatCapacity.C = -5.8400E00;
    out.solHeatCapacity.D = 8.5900E-03;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [15.00, 3.6410E03; 168.35, 9.4841E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.6050E05;
    out.liqHeatCapacity.B = -2.2930E02;
    out.liqHeatCapacity.C = 7.86600E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [168.35, 1.4417E05; 310.00, 1.6495E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 7.4976E04;
    out.IG_HeatCapacity.B = 1.6166E05;
    out.IG_HeatCapacity.C = 8.6287E02;
    out.IG_HeatCapacity.D = 7.8964E05;
    out.IG_HeatCapacity.E = 4.6718E03;
    out.IG_HeatCapacity.range = [298.15, 9.1671E04; 1000.00, 2.0789E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.6669E-01;
    out.SecondViralCoef.B = -1.6805E02;
    out.SecondViralCoef.C = -1.9195E07;
    out.SecondViralCoef.D = 2.9849E19;
    out.SecondViralCoef.E = -8.7864E21;
    out.SecondViralCoef.range = [240.30, -2.5200; 2400.00, 9.5290E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -6.0611E00;
    out.liqViscosity.B = 6.3620E02;
    out.liqViscosity.C = -7.1926E-1;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [168.35, 2.5567E-03; 378.50, 1.7514E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 4.3034E-07;
    out.vapViscosity.B = 6.1720E-01;
    out.vapViscosity.C = 2.2554E+02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [315, 8.7348E-06; 975.00, 2.4449E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.5239E-01;
    out.liqThermalConductivity.B = -3.6877E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [168.35, 1.9031E-01; 315, 1.3623E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 1.1366E-03;
    out.vapThermalConductivity.B = 6.9110E-01;
    out.vapThermalConductivity.C = 4.1879E+02;
    out.vapThermalConductivity.D = 1.8761E+05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [315.00, 1.4350E-02; 1000.00, 8.3763E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 7.0967E-02;
    out.SurfaceTension.B = 1.2676;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [168.35, 4.1083E-02; 480.60, 0.00E00];
    out.SurfaceTension.eqn = 106;    
    
end


