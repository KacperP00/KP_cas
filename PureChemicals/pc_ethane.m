function out = pc_ethane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'ethane','bimethyl','ethyl hydride','methylmethane'}; 
    out.CAS = '74-64- 0';
    % components
    out.ChemicalFormular.C = 2;         
    out.ChemicalFormular.H = 6;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 30.070;   % kg/kmol    
    out.Tcrit = 305.32;             % K    
    out.pcrit = 4.8720E6;           % Pa 
    out.vcrit = 0.14550;            % m^3/kmol 
    out.CritCompressFactor = 0.279; % -
    
%% 2 
    out.MeltingPoint = 90.35;              % K
    out.TripletPointT = 90.35;             % K
    out.TripletPointp = 1.1300E00;          % Pa
    out.NormalBoilingPoint = 184.55;        % K
    out.liqMolarVolume = 0.055229;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -8.3820E07;             % J/kmol
    out.IG_GibbsOfFormation = -3.1920E07;            % J/kmol
    out.IG_AbsoluteEntropy =  2.2912E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.6590E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.4286E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.0995;            % -
    out.RadiusOfGyration = 1.8260E-10;      % m
    out.SolubilityParameter = 1.2400E04;    % (J/m^3)^.5 
    out.DipoleMoment = 0.0;                 % C^m
    out.VanDerWaalsVolume = 0.02734;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 4.2400E08;        % m^2/kmol
    out.RefractiveIndex = 1.1849;           % -
    out.FlashPoint = NaN;                   % K
    out.FlammabilityLimits = [2.9 13.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 745.00;          % K
    
%% 6 
    out.solDensity.A = 2.5263E01;
    out.solDensity.B = -1.0950E-02;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [23.15, 2.5010E01; 90.35, 2.4274E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.9122E00;
    out.liqDensity.B = 2.7937E-01;
    out.liqDensity.C = 3.0532E02;
    out.liqDensity.D = 2.9187E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [90.35, 2.1640E01; 305.32, 6.8447E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.1857E01;
    out.vapPressure.B = -2.5987E03;
    out.vapPressure.C = -5.1283E00;
    out.vapPressure.D =  1.4913E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [90.35, 1.1273E00; 305.32, 4.8522E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 2.1091E07;
    out.HeatOfVap.B = 6.0646E-01;
    out.HeatOfVap.C = -5.5492E-01;
    out.HeatOfVap.D = 3.2799E-01;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [90.35, 1.7879E07; 305.32, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 8.1584E03;
    out.solHeatCapacity.B = -1.4778E03;
    out.solHeatCapacity.C = 9.7183E01;
    out.solHeatCapacity.D = -1.6114E00;
    out.solHeatCapacity.E = 9.1233E-03;
    out.solHeatCapacity.range = [15.53, 3.1420E03; 63.03, 4.1591E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 4.4009E01;
    out.liqHeatCapacity.B = 8.9718E04;
    out.liqHeatCapacity.C = 9.1877E02;
    out.liqHeatCapacity.D = -1.8860E03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [92.00, 6.8554E04; 290.00, 1.2444E05];
    out.liqHeatCapacity.eqn = 114;    
    
    out.IG_HeatCapacity.A = 4.0326E04;
    out.IG_HeatCapacity.B = 1.3422E05;
    out.IG_HeatCapacity.C = 1.6555E03;
    out.IG_HeatCapacity.D = 7.3223E04;
    out.IG_HeatCapacity.E = 7.5287E02;
    out.IG_HeatCapacity.range = [200.00, 4.2557E04; 1500.00, 1.4562E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 8.0950E-02;
    out.SecondViralCoef.B = -6.1710E01;
    out.SecondViralCoef.C = -1.4350E06;
    out.SecondViralCoef.D = 6.7600E16;
    out.SecondViralCoef.E = -9.7400E18;
    out.SecondViralCoef.range = [152.71, -7.1318E-01; 1500.00, 3.9385E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -3.4134E00;
    out.liqViscosity.B = 1.9705E02;
    out.liqViscosity.C = -1.2193E00;
    out.liqViscosity.D = -9.2023E-26;
    out.liqViscosity.E = 1E01;
    out.liqViscosity.range = [90.35, 1.2020E-03; 300.00, 3.5196E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.1005E-08;
    out.vapViscosity.B = 9.7620E-01;
    out.vapViscosity.C = 2.3139E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range =[90.35, 2.6429E-06; 1000.00, 2.5827E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 3.5758E-01;
    out.liqThermalConductivity.B = -1.1458E-03;
    out.liqThermalConductivity.C = 6.1866E-07;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [90.35, 2.5911E-01; 300.00, 6.9519E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 7.3869E-05;
    out.vapThermalConductivity.B = 1.1689E00;
    out.vapThermalConductivity.C = 5.0073E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [164.55, 8.8628E-03; 1000.00, 1.5807E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 4.8643E-02;
    out.SurfaceTension.B = 1.1981E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [90.35, 3.1949E-02; 305.32, 0];
    out.SurfaceTension.eqn = 106;    
    
end


