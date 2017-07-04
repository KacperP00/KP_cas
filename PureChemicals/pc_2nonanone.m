function out = pc_2nonanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2nonanone','2-nonanone','2-Nonanone'}; 
    out.CAS = '821-55-6';
    % components
    out.ChemicalFormular.C = 9;         
    out.ChemicalFormular.H = 18;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 142.241;  % kg/kmol    
    out.Tcrit = 652.50;             % K    
    out.pcrit = 2.4100E6;           % Pa 
    out.vcrit = 0.56000;            % m^3/kmol 
    out.CritCompressFactor = 0.249; % -
    
%% 2 
    out.MeltingPoint = 265.65;              % K
    out.TripletPointT = 265.65;             % K
    out.TripletPointp = 4.7226E00;          % Pa
    out.NormalBoilingPoint = 467.15;        % K
    out.liqMolarVolume = 0.173947;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.4070E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.0410E8;            % J/kmol
    out.IG_AbsoluteEntropy =  5.3560E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.5600E7;          % J/kmol
    out.StandardNetHeatOfCombustion = -5.3212E9;    % J/kmol
    out.StandardGibbsOfFormation = - 1.2270E8;              % J/kmol
    out.StandardHeatOfFormation = -3.9680E8;               % J/kmol
    out.StandardAbsoluteEntropy = 4.0830E5;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4979;            % -
    out.RadiusOfGyration = 5.3210E-10;      % m
    out.SolubilityParameter = 1.7560E4;     % (J/m^3)^.5 
    out.DipoleMoment = 9.1400E-30;          % C^m
    out.VanDerWaalsVolume = 0.10042;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.3940E09;        % m^2/kmol
    out.RefractiveIndex = 1.4187;           % -
    out.FlashPoint = 337.15;                % K
    out.FlammabilityLimits = [0.9 5.9];    % vol-%
    out.FlammabilityLimitTemps = [340 380]; % K
    out.AutoignitionTemp = 633.15;          % K
    
%% 6 
    out.solDensity.A = 6.9300E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [265.65, 6.9300E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 4.4319E-01;
    out.liqDensity.B = 2.4818E-01;
    out.liqDensity.C = 6.5250E02;
    out.liqDensity.D = 2.8763E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [265.65, 5.9233E00; 652.50, 1.7858E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.9160E01;
    out.vapPressure.B = -7.8739E03;
    out.vapPressure.C = -5.0102E00;
    out.vapPressure.D =  9.1592E-19;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [265.65, 4.7226E00; 652.50, 2.4054E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 7.1729E07;
    out.HeatOfVap.B = 4.0270E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [265.65, 5.8112E07; 652.50, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 2.5700E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [265.65, 2.5700E05; NaN,NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.5400E05;
    out.liqHeatCapacity.B = -5.5822E02;
    out.liqHeatCapacity.C = 1.3135E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [265.65, 2.9840E05; 510.00, 4.1095E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.5347E05;
    out.IG_HeatCapacity.B = 4.2140E05;
    out.IG_HeatCapacity.C = 1.3160E03;
    out.IG_HeatCapacity.D = 2.4060E05;
    out.IG_HeatCapacity.E = 6.4420E02;
    out.IG_HeatCapacity.range = [200.00, 1.6947E05; 1500.00, 5.1866E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.8301E-01;
    out.SecondViralCoef.B = -4.4922E02;
    out.SecondViralCoef.C = -1.3676E08;
    out.SecondViralCoef.D = -8.4430E20;
    out.SecondViralCoef.E = 1.1643E23;
    out.SecondViralCoef.range = [326.25, -8.7298E00; 1500.00, 4.2979E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1678E01;
    out.liqViscosity.B = 1.5184E03;
    out.liqViscosity.C = -2.9833E-02;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [265.65, 2.1793E-03; 467.15, 1.8209E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 7.3861E-08;
    out.vapViscosity.B = 7.9750E-01;
    out.vapViscosity.C = 9.7893E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [265.65, 4.6297E-06; 1000.00, 1.6609E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.9050E-01;
    out.liqThermalConductivity.B = -1.9230E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [265.65, 1.3942E-01; 593.15, 7.6437E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -2.2439E-03;
    out.vapThermalConductivity.B = 1.0196E00;
    out.vapThermalConductivity.C = -2.6430E04;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [467.15, 2.1276E-02; 1000.0, 1.0103E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0312E-02;
    out.SurfaceTension.B = 1.0822E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [265.65, 2.8574E-02; 652.50, 0];
    out.SurfaceTension.eqn = 106;    
    
end


