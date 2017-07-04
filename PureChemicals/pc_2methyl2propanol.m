function out = pc_2methyl2propanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2methyl2propanol','2methyl-2propanol','t-Butanol'}; 
    out.CAS = '75-65-0';
    % components
    out.ChemicalFormular.C = 2;         
    out.ChemicalFormular.H = 6;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 74.123;   % kg/kmol    
    out.Tcrit = 506.21;             % K    
    out.pcrit = 3.9730E6;           % Pa 
    out.vcrit = 0.27500;            % m^3/kmol 
    out.CritCompressFactor = 0.260; % -
    
%% 2 
    out.MeltingPoint = 298.97;              % K
    out.TripletPointT = 298.97;             % K
    out.TripletPointp = 5.9356E03;          % Pa
    out.NormalBoilingPoint = 355.57;        % K
    out.liqMolarVolume = 0.094744;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.1240E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.7760E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.2630E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  6.7030E6;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.4239E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.6115;            % -
    out.RadiusOfGyration = 3.0670E-10;      % m
    out.SolubilityParameter = 2.1600E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.5700E-30;          % C^m
    out.VanDerWaalsVolume = 0.05238;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.6200E08;        % m^2/kmol
    out.RefractiveIndex = 1.3852;           % -
    out.FlashPoint = 284.26;                % K
    out.FlammabilityLimits = [2.4 8.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 751.00;          % K
    
%% 6 
    out.solDensity.A = 1.3000E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [298.97, 1.3000E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.2120E-01;
    out.liqDensity.B = 2.5440E-01;
    out.liqDensity.C = 5.0621E02;
    out.liqDensity.D = 2.7600E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [298.97, 1.0555E01; 506.21, 3.6211E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.7231E02;
    out.vapPressure.B = -1.1590E04;
    out.vapPressure.C = -2.2118E01;
    out.vapPressure.D =  1.3709E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [298.97, 5.9356E03; 506.21, 3.9910E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 7.7320E07;
    out.HeatOfVap.B = 5.6450E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [298.97, 4.6703E07; 506.21, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 4.6000E03;
    out.solHeatCapacity.B = 4.1200E02;
    out.solHeatCapacity.C = 4.1750E00;
    out.solHeatCapacity.D = -3.4270E-02;
    out.solHeatCapacity.E = 8.0080E-05;
    out.solHeatCapacity.range = [25.00, 1.7005E04; 286.13, 1.9826E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = -9.2546E05;
    out.liqHeatCapacity.B = 7.8949E03;
    out.liqHeatCapacity.C = -1.7661E01;
    out.liqHeatCapacity.D = 1.3617E-02;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [298.96, 2.2016E05; 460.00, 2.9455E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 7.7040E04;
    out.IG_HeatCapacity.B = 2.5390E05;
    out.IG_HeatCapacity.C = 1.5502E03;
    out.IG_HeatCapacity.D = 1.6690E05;
    out.IG_HeatCapacity.E = -6.7930E02;
    out.IG_HeatCapacity.range = [200.00, 8.5670E04; 1500.00, 2.8508E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.8698E-01;
    out.SecondViralCoef.B = -1.5686E02;
    out.SecondViralCoef.C = -3.2431E07;
    out.SecondViralCoef.D = -3.7866E19;
    out.SecondViralCoef.E =-2.8168E21;
    out.SecondViralCoef.range = [253.11, -4.0200E00; 1500.00, 7.2796E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.1640E02;
    out.liqViscosity.B = 1.3205E04;
    out.liqViscosity.C = 2.9254E01;
    out.liqViscosity.D = -2.4616E-27;
    out.liqViscosity.E = 1.0000E01;
    out.liqViscosity.range = [298.15, 4.3119E-03; 451.21, 1.0137E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 9.6050E-07;
    out.vapViscosity.B = 4.8560E-01;
    out.vapViscosity.C = 3.8100E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [298.97, 6.7267E-06; 600.00, 1.3124E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.5960E-01;
    out.liqThermalConductivity.B = -1.4700E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [293.15, 1.1651E-01; 479.15, 8.9165E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -1.4000E-02;
    out.vapThermalConductivity.B = 3.6250E-01;
    out.vapThermalConductivity.C = -2.3760E03;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [355.57, 2.0717E-02; 1000.0, 1.2445E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 4.0190E-02;
    out.SurfaceTension.B = -5.3100E-05;
    out.SurfaceTension.C = -5.1200E-08;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [298.15, 1.9807E-02; 355.57, 1.4836E-02];
    out.SurfaceTension.eqn = 100;    
    
end


