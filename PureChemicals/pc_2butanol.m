function out = pc_2butanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2butanol','2-butanol','2-Butanol','2Butanol','Butane'}; 
    out.CAS = '78-92-2';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 74.123;   % kg/kmol    
    out.Tcrit = 536.05;             % K    
    out.pcrit = 4.1790E06;           % Pa 
    out.vcrit = 0.26900;            % m^3/kmol 
    out.CritCompressFactor = 0.252; % -
    
%% 2 
    out.MeltingPoint = 158.45;              % K
    out.TripletPointT = 158.45;             % K
    out.TripletPointp = 1.1323E-06;          % Pa
    out.NormalBoilingPoint = 372.70;        % K
    out.liqMolarVolume = 0.092017;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.9290E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.6960E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.6469E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  5.9710E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.4408E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5722;            % -
    out.RadiusOfGyration = 3.2030E-10;      % m
    out.SolubilityParameter = 2.2540E04;     % (J/m^3)^.5 
    out.DipoleMoment = 5.5400E-30;          % C^m
    out.VanDerWaalsVolume = 0.05239;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.6200E08;        % m^2/kmol
    out.RefractiveIndex = 1.3949;           % -
    out.FlashPoint = 297.00;                % K
    out.FlammabilityLimits = [1.7 9.8];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 679.00;          % K
    
%% 6 
    out.solDensity.A = 1.4700E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [158.45, 1.4700E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.6600E-01;
    out.liqDensity.B = 2.6064E-01;
    out.liqDensity.C = 5.3605E02;
    out.liqDensity.D = 2.7460E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [158.45, 1.2570E01; 536.05, 3.7063E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A = 1.5254E02;
    out.vapPressure.B = -1.1111E04;
    out.vapPressure.C = -1.9025E01;
    out.vapPressure.D = 1.0426E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [158.45, 1.1323E-06; 536.05, 4.2014E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 7.2560E07;
    out.HeatOfVap.B = 4.7740E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [158.45, 6.1383E07; 536.05, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.0530E04;
    out.solHeatCapacity.B = 1.2313E03;
    out.solHeatCapacity.C = -6.8730E00;
    out.solHeatCapacity.D = 1.9800E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [15.00, 6.4599E03; 150.00, 8.6347E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.0670E05;
    out.liqHeatCapacity.B = -1.0204E03;
    out.liqHeatCapacity.C = 3.2900E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [158.45, 1.2762E05; 372.70, 2.8340E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 8.2020E04;
    out.IG_HeatCapacity.B = 2.5220E05;
    out.IG_HeatCapacity.C = 1.6010E03;
    out.IG_HeatCapacity.D = 1.5864E05;
    out.IG_HeatCapacity.E = -7.0415E02;
    out.IG_HeatCapacity.range = [200.00, 8.8896E04; 1500.00, 2.8513E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.8777E-01;
    out.SecondViralCoef.B = -1.7264E02;
    out.SecondViralCoef.C = -3.6937E07;
    out.SecondViralCoef.D = -5.7699E19;
    out.SecondViralCoef.E = 4.7178E21;
    out.SecondViralCoef.range = [268.03, -3.8800E00; 1500.00, 6.1730E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.0638E02;
    out.liqViscosity.B = 7.4343E03;
    out.liqViscosity.C = 1.3285E01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [288.15, 4.7932E-03; 373.15, 4.1659E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.2534E-07;
    out.vapViscosity.B = 7.6560E-01;
    out.vapViscosity.C = 9.2847E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [158.45, 3.8197E-06; 1000.00, 2.2716E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0460E-01;
    out.liqThermalConductivity.B = -2.3542E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [200.00, 1.5752E-01; 398.15, 1.1087E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.0086E03;
    out.vapThermalConductivity.B = 9.7450E-01;
    out.vapThermalConductivity.C = 1.5008E10;
    out.vapThermalConductivity.D = -9.2508E11;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [372.70, 1.9153E-02; 1000.0, 1.1959E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 4.6450E-02;
    out.SurfaceTension.B = -7.4730E-05;
    out.SurfaceTension.C = -1.3000E-08;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [273.15, 2.5068E-02; 372.70, 1.6792E-02];
    out.SurfaceTension.eqn = 100;    
    
end


