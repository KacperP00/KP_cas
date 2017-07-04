function out = pc_1dodecanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1dodecanol','1-dodecanol','1-Dodecanol','n-dodecanol','Lauryl 24'}; 
    out.CAS = '112-53-8';
    % components
    out.ChemicalFormular.C = 12;         
    out.ChemicalFormular.H = 26;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 186.338;  % kg/kmol    
    out.Tcrit = 719.4;              % K    
    out.pcrit = 1.9940E06;          % Pa 
    out.vcrit = 0.716;              % m^3/kmol 
    out.CritCompressFactor = 0.239; % -
    
%% 2 
    out.MeltingPoint = 296.95;              % K
    out.TripletPointT = 296.95;             % K
    out.TripletPointp = 0.101105;           % Pa
    out.NormalBoilingPoint = 536.95;        % K
    out.liqMolarVolume = 0.224322;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -4.3670E08;             % J/kmol
    out.IG_GibbsOfFormation = -8.0820E07;            % J/kmol
    out.IG_AbsoluteEntropy =  6.7521E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = 3.1380E07;           % J/kmol
    out.StandardNetHeatOfCombustion = -7.3380E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.1500E08;       % J/kmol
    out.StandardHeatOfFormation = -5.2850E08;        % J/kmol
    out.StandardAbsoluteEntropy = 4.8300E05;         % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.666353;          % -
    out.RadiusOfGyration = 6.1190E-10;      % m
    out.SolubilityParameter = 1.9750E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.6400E-30;          % C^m
    out.VanDerWaalsVolume = 0.13424;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.8430E09;        % m^2/kmol
    out.RefractiveIndex = 1.4413;           % -
    out.FlashPoint = 400.0;                 % K
    out.FlammabilityLimits = [0.6 5.1];    % vol-%
    out.FlammabilityLimitTemps = [397 441]; % K
    out.AutoignitionTemp = 548.15;          % K
    
%% 6 
    out.solDensity.A = 5.5881E00;
    out.solDensity.B = -1.6655E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [78.15, 5.4579E00; 296.95, 5.0935E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.5818E-01;
    out.liqDensity.B = 2.5626E-01;
    out.liqDensity.C = 7.1940E02;
    out.liqDensity.D = 2.9963E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [296.95, 4.4623E00; 719.40, 1.3977E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  3.8678E02;
    out.vapPressure.B = -2.4382E04;
    out.vapPressure.C = -5.6170E01;
    out.vapPressure.D =  4.3253E-02;
    out.vapPressure.E = 1.0000E00;
    out.vapPressure.range = [296.95,1.011E-01; 719.40,1.9716E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 1.2925E08;
    out.HeatOfVap.B = 6.7609E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [296.95, 9.0182E07; 719.40, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 3.5700E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [293.15, 3.5700E05; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 5.1422E06;
    out.liqHeatCapacity.B = -5.4970E04;
    out.liqHeatCapacity.C = 2.2872E02;
    out.liqHeatCapacity.D = -4.0331E-01;
    out.liqHeatCapacity.E = 2.5844E-04;
    out.liqHeatCapacity.range = [296.95,4.3614E05; 486.00, 5.7106E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.0150E05;
    out.IG_HeatCapacity.B = 6.3270E05;
    out.IG_HeatCapacity.C = 1.5625E03;
    out.IG_HeatCapacity.D = 4.6733E05;
    out.IG_HeatCapacity.E = 7.2120E02;
    out.IG_HeatCapacity.range = [200.00, 2.1943E05; 1500.00, 7.3406E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 6.2391E-01;
    out.SecondViralCoef.B = -7.8675E02;
    out.SecondViralCoef.C = -1.6374E08;
    out.SecondViralCoef.D = 7.1921E21;
    out.SecondViralCoef.E = -2.8712E24;
    out.SecondViralCoef.range = [359.70,-7.9011E00; 1500.00, 5.1100E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -5.3664E01;
    out.liqViscosity.B = 5.2265E03;
    out.liqViscosity.C = 5.6140E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [293.15, 1.9373E-02; 536.95, 1.7667E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 5.5812E-08;
    out.vapViscosity.B = 8.2000E-01;
    out.vapViscosity.C = 5.4330E01;
    out.vapViscosity.D = 1.4110E04;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [296.95, 4.4285E-06; 1000.00, 1.5065E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0444E-01;
    out.liqThermalConductivity.B = -1.8378E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [296.95, 1.4987E-01; 536.95, 1.0576E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -2.8300E-03;
    out.vapThermalConductivity.B = 9.1261E-01;
    out.vapThermalConductivity.C = -1.9181E04;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [536.95, 2.5266E-02; 1000.00, 8.5114E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.3034E-02;
    out.SurfaceTension.B = 1.1008E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [293.15, 2.9808E-02; 719.40, 0];
    out.SurfaceTension.eqn = 106;    
    
end


