function out = pc_dinbutylether

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'dinbutylether','di-n-butyl-ether','butyl oxide','1,1-oxybisbutane'}; 
    out.CAS = '142-96-1';
    % components
    out.ChemicalFormular.C = 8;         
    out.ChemicalFormular.H = 18;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 130.23;   % kg/kmol    
    out.Tcrit = 584.1;              % K    
    out.pcrit = 2.4600E6;           % Pa 
    out.vcrit = 0.487;              % m^3/kmol 
    out.CritCompressFactor = 0.247; % -
    
%% 2 
    out.MeltingPoint = 175.3;               % K
    out.TripletPointT = 175.3;              % K
    out.TripletPointp = 7.1413E-04;         % Pa
    out.NormalBoilingPoint = 414.15;        % K
    out.liqMolarVolume = 0.170405;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.3340E08;             % J/kmol
    out.IG_GibbsOfFormation = -8.8270E07;            % J/kmol
    out.IG_AbsoluteEntropy =  5.0140E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.6300E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.9469E09;    % J/kmol
    out.StandardGibbsOfFormation = -9.9490E07;       % J/kmol
    out.StandardHeatOfFormation = -3.7790E08;        % J/kmol
    out.StandardAbsoluteEntropy = 3.8980E05;         % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.447646;           % -
    out.RadiusOfGyration = 4.7970E-10;      % m
    out.SolubilityParameter = 1.5830E04;    % (J/m^3)^.5 
    out.DipoleMoment = 3.9000E-30;          % C^m
    out.VanDerWaalsVolume = 0.09192;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.2880E09;        % m^2/kmol
    out.RefractiveIndex = 1.3968;           % -
    out.FlashPoint = 298.15;                % K
    out.FlammabilityLimits = [1.5 7.6];     % vol-%
    out.FlammabilityLimitTemps = [308 341]; % K
    out.AutoignitionTemp =467.15;           % K
    
%% 6 
    out.solDensity.A = 7.7300E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [175.30, 7.730E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.5941E-01;
    out.liqDensity.B = 2.7243E-01;
    out.liqDensity.C = 5.8410E02;
    out.liqDensity.D = 2.9932E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [175.30, 6.6071E00; 584.10,2.0534E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A = 7.2227E01;
    out.vapPressure.B = -7.5376E03;
    out.vapPressure.C = -7.0596E00;
    out.vapPressure.D =  9.1442E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [175.30, 7.1414E-04; 584.10, 2.4592E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.9616E07;
    out.HeatOfVap.B = 3.8833E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [175.30, 5.1902E07; 584.10, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 0;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [NaN, NaN; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.7072E05;
    out.liqHeatCapacity.B = -2.5983E02;
    out.liqHeatCapacity.C = 9.5427E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [175.30, 2.5450E05; 450.00, 3.4704E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.6122E05;
    out.IG_HeatCapacity.B = 4.4777E05;
    out.IG_HeatCapacity.C = 1.6831E03;
    out.IG_HeatCapacity.D = 2.9180E05;
    out.IG_HeatCapacity.E = 7.8160E02;
    out.IG_HeatCapacity.range = [200.00, 1.6841E05; 1500.00, 5.2145E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.3911E-01;
    out.SecondViralCoef.B = -3.7344E02;
    out.SecondViralCoef.C = -7.7353E06;
    out.SecondViralCoef.D = -1.3080E20;
    out.SecondViralCoef.E = 4.4188E21;
    out.SecondViralCoef.range = [292.05, -6.2304E00; 1500.00, 6.7226E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = 1.0027E01;
    out.liqViscosity.B = 2.0600E02;
    out.liqViscosity.C = -3.1607E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [175.30, 5.9306E-03; 414.15, 1.9890E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 7.7147E-08;
    out.vapViscosity.B = 7.9906E-01;
    out.vapViscosity.C = 8.0765E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [175.30, 3.2784E-06; 1000.00, 1.7814E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.9418E-01;
    out.liqThermalConductivity.B = -2.2246E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [175.30, 1.5518E-01; 523.15, 7.7800E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 3.2694E-03;
    out.vapThermalConductivity.B = 5.8633E-01;
    out.vapThermalConductivity.C = 1.2599E03;
    out.vapThermalConductivity.D = 3.0089E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [323.15, 1.2440E-02; 1000.00, 7.3297E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1346E-02;
    out.SurfaceTension.B = 1.1604E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [175.30, 3.3937E-02; 584.10, 0];
    out.SurfaceTension.eqn = 106;    
    
end


