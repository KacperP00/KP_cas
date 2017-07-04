function out = pc_2pentanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2pentanol','2-pentanol','pentanol-2'}; 
    out.CAS = '6032-29-7';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 88.1497;  % kg/kmol    
    out.Tcrit = 560.3;              % K    
    out.pcrit = 3.6750E6;           % Pa 
    out.vcrit = 0.329;              % m^3/kmol 
    out.CritCompressFactor = 0.26;  % -
    
%% 2 
    out.MeltingPoint = 200.0;              % K
    out.TripletPointT = 200.0;             % K
    out.TripletPointp = 5.9659E-3;         % Pa
    out.NormalBoilingPoint = 392.15;       % K
    out.liqMolarVolume = 0.10957;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.1380E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.5930E8;            % J/kmol
    out.IG_AbsoluteEntropy =  3.9640E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  8.4800E06;         % J/kmol
    out.StandardNetHeatOfCombustion = -3.0515E09;   % J/kmol
    out.StandardGibbsOfFormation = -1.6850E08;      % J/kmol
    out.StandardHeatOfFormation = -3.6520E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.5500E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.558707;          % -
    out.RadiusOfGyration = 3.6190E-10;      % m
    out.SolubilityParameter = 2.1700E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.3900E-30;          % C^m
    out.VanDerWaalsVolume = 0.06262;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.9700E08;        % m^2/kmol
    out.RefractiveIndex = 1.4044;           % -
    out.FlashPoint = 307.0;                 % K
    out.FlammabilityLimits = [1.5 9.7];     % vol-%
    out.FlammabilityLimitTemps = [307 338]; % K
    out.AutoignitionTemp = 616.0;           % K
    
%% 6 
    out.solDensity.A = 1.1700E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [NaN, NaN];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.0126E-01;
    out.liqDensity.B = 2.6371E-01;
    out.liqDensity.C = 5.6030E02;
    out.liqDensity.D = 2.5302E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [200.00, 1.0008E01; 560.30, 3.0384E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.1227E02;
    out.vapPressure.B = -1.0221E04;
    out.vapPressure.C = -1.2511E01;
    out.vapPressure.D =  8.5670E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [200.00, 5.9659E-03; 560.30, 3.6767E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 8.1306E07;
    out.HeatOfVap.B = 5.3685E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [200.00, 6.4147E07; 560.30, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 0;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [NaN,NaN; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.7400E05;
    out.liqHeatCapacity.B = -6.5134E01;
    out.liqHeatCapacity.C = 9.5351E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [280.15,2.3059E05; 420.22, 3.1500E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0385E05;
    out.IG_HeatCapacity.B = 2.9140E05;
    out.IG_HeatCapacity.C = 1.6383E03;
    out.IG_HeatCapacity.D = 1.9435E05;
    out.IG_HeatCapacity.E = 7.5094E02;
    out.IG_HeatCapacity.range = [298.15, 1.3604E05; 1500.00, 3.4084E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.3739E-01;
    out.SecondViralCoef.B = -2.4456E02;
    out.SecondViralCoef.C = -3.2955E07;
    out.SecondViralCoef.D = 3.8515E20;
    out.SecondViralCoef.E = -1.6083E23;
    out.SecondViralCoef.range = [280.15, -7.1139E00; 1500.00, 6.4596E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.0809E02;
    out.liqViscosity.B = 7.5182E03;
    out.liqViscosity.C = 1.3552E01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [250.00, 4.1157E-02; 392.15, 3.3865E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.1030E-07;
    out.vapViscosity.B = 6.4060E-01;
    out.vapViscosity.C = 2.3282E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [200.00, 4.2711E-06; 1000.00, 2.1022E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2037E-01;
    out.liqThermalConductivity.B = -2.8791E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [200.00, 1.6279E-01; 435.60, 9.4956E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.8508E03;
    out.vapThermalConductivity.B = 9.0265E-01;
    out.vapThermalConductivity.C = 1.2060E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [392.15, 2.0326E-02; 1000.0, 1.2066E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.5779E-02;
    out.SurfaceTension.B = 1.1409E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [200.00, 3.3705E-02; 560.30, 0];
    out.SurfaceTension.eqn = 106;    
    
end


