function out = pc_2butanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2butanone','2-butanone','2-Butanone','2-butanone','butanone'}; 
    out.CAS = '78-93-3';
    % components
    out.ChemicalFormular.C = 2;         
    out.ChemicalFormular.H = 6;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 72.107;   % kg/kmol    
    out.Tcrit = 535.50;             % K    
    out.pcrit = 4.1543E06;           % Pa 
    out.vcrit = 0.26700;            % m^3/kmol 
    out.CritCompressFactor = 0.249; % -
    
%% 2 
    out.MeltingPoint = 186.48;              % K
    out.TripletPointT = 186.48;             % K
    out.TripletPointp = 1.3941E00;          % Pa
    out.NormalBoilingPoint = 352.79;        % K
    out.liqMolarVolume = 0.090204;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.3495E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.4606E8;            % J/kmol
    out.IG_AbsoluteEntropy =  3.3811E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  8.4391E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.2616E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3241;            % -
    out.RadiusOfGyration = 3.1350E-10;      % m
    out.SolubilityParameter = 1.8796E04;    % (J/m^3)^.5 
    out.DipoleMoment = 9.2064E-30;          % C^m
    out.VanDerWaalsVolume = 0.04927;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.1900E08;        % m^2/kmol
    out.RefractiveIndex = 1.3764;           % -
    out.FlashPoint = 267.04;                % K
    out.FlammabilityLimits = [1.8 10.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 788.71;          % K
    
%% 6 
    out.solDensity.A = 1.4700E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [186.48, 1.4700E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.3300E-01;
    out.liqDensity.B = 2.514E-01;
    out.liqDensity.C = 5.3550E02;
    out.liqDensity.D = 2.8570E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [186.48, 1.2593E01; 535.50, 3.7112E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.1474E02;
    out.vapPressure.B = -7.1300E03;
    out.vapPressure.C = -1.5184E01;
    out.vapPressure.D =  1.7234E-02;
    out.vapPressure.E =  1.0000E00;
    out.vapPressure.range = [186.48, 1.3941E00; 535.50, 4.1976E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.1400E07;
    out.HeatOfVap.B = 7.5450E-01;
    out.HeatOfVap.C = -5.9500E-01;
    out.HeatOfVap.D = 2.3300E-01;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [186.48, 4.0176E07; 535.50, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -9.712E03;
    out.solHeatCapacity.B = 1.3640E03;
    out.solHeatCapacity.C = -7.1303E00;
    out.solHeatCapacity.D = 1.6570E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [10.00, 3.2315E03; 186.48, 1.0415E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.3230E05;
    out.liqHeatCapacity.B = 2.0087E02;
    out.liqHeatCapacity.C = -9.5970E-01;
    out.liqHeatCapacity.D = 1.9533E-03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [186.48, 1.4905E05; 373.15, 1.7511E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 7.8400E04;
    out.IG_HeatCapacity.B = 2.1032E05;
    out.IG_HeatCapacity.C = 1.5488E03;
    out.IG_HeatCapacity.D = 1.1855E05;
    out.IG_HeatCapacity.E = 6.9300E02;
    out.IG_HeatCapacity.range = [200.00, 8.3967E04; 1500.00, 2.4816E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 9.0100E-02;
    out.SecondViralCoef.B = -8.6700E01;
    out.SecondViralCoef.C = -6.0500E07;
    out.SecondViralCoef.D = 3.4650E20;
    out.SecondViralCoef.E = -1.0330E23;
    out.SecondViralCoef.range = [277.75, -3.7632E00; 2677.50, 5.4567E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -8.1340E-01;
    out.liqViscosity.B = 5.0978E02;
    out.liqViscosity.C = -1.5324E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [188.71, 2.1504E-03; 533.15, 7.6450E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.9759E-08;
    out.vapViscosity.B = 1.0252E00;
    out.vapViscosity.C = -1.6480E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [273.15, 6.6159E-06; 573.15, 1.3684E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1920E-01;
    out.liqThermalConductivity.B = -2.4840E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [186.48, 1.7288E-01; 352.79, 1.3157E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -4.2027E06;
    out.vapThermalConductivity.B = -1.5240E-01;
    out.vapThermalConductivity.C = 2.0846E09;
    out.vapThermalConductivity.D = -1.4577E13;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [352.79, 1.5457E-02; 992.79, 1.1571E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.9441E-02;
    out.SurfaceTension.B = 1.1165E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [186.48, 2.4101E-02; 355.50, 0];
    out.SurfaceTension.eqn = 106;    
    
end


