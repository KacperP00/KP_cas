function out = pc_1pentanal

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1pentanal','1-pentanal','1-Pentanal','1Pentanal','n-Pentanal','n-pentanal'}; 
    out.CAS = '110-62-3';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 86.134;   % kg/kmol    
    out.Tcrit = 566.10;             % K    
    out.pcrit = 3.9700E06;          % Pa 
    out.vcrit = 0.31300;            % m^3/kmol 
    out.CritCompressFactor = 0.264; % -
    
%% 2 
    out.MeltingPoint = 182.00;              % K
    out.TripletPointT = 182.00;             % K
    out.TripletPointp = 5.2282E-02;         % Pa
    out.NormalBoilingPoint = 376.15;        % K
    out.liqMolarVolume = 0.107035;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.2780E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.0710E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.8289E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.5000E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.9100E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3472;            % -
    out.RadiusOfGyration = 3.5160E-10;      % m
    out.SolubilityParameter = 1.8440E04;    % (J/m^3)^.5 
    out.DipoleMoment = 8.5700E-30;          % C^m
    out.VanDerWaalsVolume = 0.05950;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.5400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3917;           % -
    out.FlashPoint = 285.15;                % K
    out.FlammabilityLimits = [1.5 9.5];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 480.00;          % K
    
%% 6 
    out.solDensity.A = 1.2300E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [182.00, 1.2300E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.3871E-01;
    out.liqDensity.B = 2.6252E-01;
    out.liqDensity.C = 5.6610E02;
    out.liqDensity.D = 2.9444E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [182.00, 1.0534E01; 566.10, 3.1948E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.4958E02;
    out.vapPressure.B = -8.8900E03;
    out.vapPressure.C = -2.0697E01;
    out.vapPressure.D =  2.2101E-02;
    out.vapPressure.E =  1.0000E00;
    out.vapPressure.range = [182.00, 5.2282E-02; 566.10, 3.9685E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.1478E07;
    out.HeatOfVap.B = 3.7541E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [182.00, 4.4502E07; 566.10, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 9.3786E03;
    out.solHeatCapacity.B = 7.4571E02;
    out.solHeatCapacity.C = -8.5714E-01;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [80.00, 6.3550E04; 150.00, 1.0195E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.1205E05;
    out.liqHeatCapacity.B = 2.5778E02;
    out.liqHeatCapacity.C = 0;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [200.00, 1.6361E05; 376.15, 2.0901E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0743E05;
    out.IG_HeatCapacity.B = 2.8363E05;
    out.IG_HeatCapacity.C = 1.9549E03;
    out.IG_HeatCapacity.D = 2.0146E05;
    out.IG_HeatCapacity.E = 8.9044E02;
    out.IG_HeatCapacity.range = [200.00, 1.0960E05; 1500.00, 3.2404E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.9740E-01;
    out.SecondViralCoef.B = -2.2484E02;
    out.SecondViralCoef.C = -4.0107E07;
    out.SecondViralCoef.D = -1.9523E20;
    out.SecondViralCoef.E = 3.0063E22;
    out.SecondViralCoef.range = [283.05, -4.5262E00; 1500.00, 3.5616E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.0846E01;
    out.liqViscosity.B = 9.8001E02;
    out.liqViscosity.C = -5.4565E-03;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [182.00, 4.1290E-03; 376.15, 2.5532E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.2700E-07;
    out.vapViscosity.B = 6.7670E-01;
    out.vapViscosity.C = 1.9174E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [182.00, 3.7404E-06; 1000.00, 2.0415E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2697E-01;
    out.liqThermalConductivity.B = -3.3227E-04;
    out.liqThermalConductivity.C = 1.1770E-07;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [182.00, 1.7040E-01; 513.15, 8.7459E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -4.9187E06;
    out.vapThermalConductivity.B = -1.0297E-01;
    out.vapThermalConductivity.C = 2.6911E09;
    out.vapThermalConductivity.D = -2.3179E13;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [376.15, 1.7048E-02; 1000.0, 1.1788E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.4779E-02;
    out.SurfaceTension.B = 1.2556E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [182.00, 3.9804E-02; 566.10, 0.0000E00];
    out.SurfaceTension.eqn = 106;    
    
end


