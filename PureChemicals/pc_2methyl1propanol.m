function out = pc_2methyl1propanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2methyl1propanol','2methyl-1propanol','Fermation Butyl Alcohol','Isobutyl Alcohol'}; 
    out.CAS = '78-83-1';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 74.123;   % kg/kmol    
    out.Tcrit = 547.78;             % K    
    out.pcrit = 4.3000E06;          % Pa 
    out.vcrit = 0.27300;            % m^3/kmol 
    out.CritCompressFactor = 0.258; % -
    
%% 2 
    out.MeltingPoint = 165.15;              % K
    out.TripletPointT = 165.15;             % K
    out.TripletPointp = 6.2545E-07;         % Pa
    out.NormalBoilingPoint = 380.81;        % K
    out.liqMolarVolume = 0.092913;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.8320E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.5490E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.4800E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  6.3220E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.4490E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5848;            % -
    out.RadiusOfGyration = 3.3320E-10;      % m
    out.SolubilityParameter = 2.2910E+04;     % (J/m^3)^.5 
    out.DipoleMoment = 5.4700E-30;          % C^m
    out.VanDerWaalsVolume = 0.05239;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.6200E08;        % m^2/kmol
    out.RefractiveIndex = 1.3938;           % -
    out.FlashPoint = 301.00;                % K
    out.FlammabilityLimits = [1.7 10.9];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 681.15;          % K
    
%% 6 
    out.solDensity.A = 1.4100E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [165.15, 1.4100E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.9680E-01;
    out.liqDensity.B = 2.7226E-01;
    out.liqDensity.C = 5.4778E02;
    out.liqDensity.D = 2.3890E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [165.15, 1.2084E01; 547.78, 3.6612E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.8779E02;
    out.vapPressure.B = -1.2955E04;
    out.vapPressure.C = -2.4285E01;
    out.vapPressure.D =  1.4262E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [165.15, 6.2544E-07; 547.78, 4.3312E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 7.6190E07;
    out.HeatOfVap.B = 5.0470E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [165.15, 6.3570E07; 547.78, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.0520E04;
    out.solHeatCapacity.B = 1.1474E03;
    out.solHeatCapacity.C = -5.7430E00;
    out.solHeatCapacity.D = 1.5260E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [15.00, 5.4503E+03; 160.00, 8.8548E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.8020E05;
    out.liqHeatCapacity.B = -7.3800E02;
    out.liqHeatCapacity.C = 2.4930E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [165.15, 1.2631E05; 380.81, 2.6069E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 8.7940E04;
    out.IG_HeatCapacity.B = 2.4160E05;
    out.IG_HeatCapacity.C = 1.7180E03;
    out.IG_HeatCapacity.D = 1.6540E05;
    out.IG_HeatCapacity.E = 7.9870E02;
    out.IG_HeatCapacity.range = [300.00, 1.1089E05; 1200.00, 2.6355E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.4800E-01;
    out.SecondViralCoef.B = -1.3400E02;
    out.SecondViralCoef.C = -2.8400E07;
    out.SecondViralCoef.D = 1.7000E18;
    out.SecondViralCoef.E = -7.7520E22;
    out.SecondViralCoef.range = [301.25, -5.1043E00; 1500.00, 5.0250E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -4.8035E01;
    out.liqViscosity.B = 4.3067E03;
    out.liqViscosity.C = 4.8948E00;
    out.liqViscosity.D = -3.5000E-28;
    out.liqViscosity.E = 1.0E01;
    out.liqViscosity.range = [210.93, 2.4105E-01; 544.26, 4.1665E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.1945E-07;
    out.vapViscosity.B = 6.4700E-01;
    out.vapViscosity.C = 2.2314E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [165.15, 3.6991E-06; 1000.00, 2.2800E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.8503E-01;
    out.liqThermalConductivity.B = -1.7855E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [165.15, 1.5554E-01; 380.81, 1.1704E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 1.8980E03;
    out.vapThermalConductivity.B = 8.6980E-01;
    out.vapThermalConductivity.C = 6.9300E09;
    out.vapThermalConductivity.D = 1.2574E07;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [380.81, 1.8322E-02; 1000.0, 1.1142E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 4.7506E-02;
    out.SurfaceTension.B = -8.6704E-05;
    out.SurfaceTension.C = 1.0150E-08;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [283.15, 2.3770E-02; 380.81, 1.5960E-02];
    out.SurfaceTension.eqn = 100;    
    
end


