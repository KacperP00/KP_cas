function out = pc_methanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'methanol','Methanol','methyl alcohol'}; 
    out.CAS = '67-56-1';
    % components
    out.ChemicalFormular.C = 1;         
    out.ChemicalFormular.H = 4;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 32.042;    % kg/kmol    
    out.Tcrit = 512.64;              % K    
    out.pcrit = 8.0970E06;           % Pa 
    out.vcrit = 0.11800;             % m^3/kmol 
    out.CritCompressFactor = 0.224;  % -
    
%% 2 
    out.MeltingPoint = 175.47;             % K
    out.TripletPointT = 175.47;            % K
    out.TripletPointp = 1.1147E-01;        % Pa
    out.NormalBoilingPoint = 337.85;       % K
    out.liqMolarVolume  = 0.040581;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.0094E08;            % J/kmol
    out.IG_GibbsOfFormation = -1.6232E08;           % J/kmol
    out.IG_AbsoluteEntropy =  2.3988E05;            % J/kmolK
    out.HeatFusionAtMeltPoint =  3.2150E06;         % J/kmol
    out.StandardNetHeatOfCombustion = -6.38200E8;   % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5640;            % -
    out.RadiusOfGyration = 1.5520E-10;      % m
    out.SolubilityParameter = 2.9590E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.6700E-30;          % C^m
    out.VanDerWaalsVolume = 0.02171;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 3.5800E08;        % m^2/kmol
    out.RefractiveIndex = 1.3265;           % -
    out.FlashPoint = 284.00;                % K
    out.FlammabilityLimits = [7.3 36.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 737.00;          % K
    
%% 6 
    out.solDensity.A = 3.0585E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [163.15, 3.0585E01; 90.67, 3.0309E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 2.2880E00;
    out.liqDensity.B = 2.6850E-01;
    out.liqDensity.C = 5.1264E02;
    out.liqDensity.D = 2.4530E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [175.47, 2.7912E01; 512.64, 8.5214E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  8.1768E01;
    out.vapPressure.B = -6.8760E03;
    out.vapPressure.C = -8.7078E00;
    out.vapPressure.D =  7.1926E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [175.47, 1.1147E-01; 512.64, 8.1402E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.2390E07;
    out.HeatOfVap.B = 3.6820E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [175.47, 4.4900E07; 512.64, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.2706E04;
    out.solHeatCapacity.B = 9.5878E02;
    out.solHeatCapacity.C = -5.2332E00;
    out.solHeatCapacity.D = 1.3152E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [25.00, 8.1982E03; 150.00, 5.7752E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.0580E05;
    out.liqHeatCapacity.B = -3.6223E02;
    out.liqHeatCapacity.C = 9.3790E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [175.47, 7.1117E04; 400.00, 1.1097E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 3.9252E04;
    out.IG_HeatCapacity.B = 8.7900E04;
    out.IG_HeatCapacity.C = 1.9165E03;
    out.IG_HeatCapacity.D = 5.3654E03;
    out.IG_HeatCapacity.E = 8.9670E02;
    out.IG_HeatCapacity.range = [200.00, 3.9802E04; 1500.00, 1.0533E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.3000E-02;
    out.SecondViralCoef.B = -1.4400E01;
    out.SecondViralCoef.C = -2.0600E07;
    out.SecondViralCoef.D = 7.6000E19;
    out.SecondViralCoef.E =-4.6200E22;
    out.SecondViralCoef.range = [298.15, -2.0770E00; 1500.00, -2.7019E-03];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.5317E01;
    out.liqViscosity.B = 1.7892E03;
    out.liqViscosity.C = 2.0690E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [175.47, 1.1928E-02; 337.85, 3.4423E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 3.0663E-07;
    out.vapViscosity.B = 6.9655E-01;
    out.vapViscosity.C = 2.0500E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [240.00, 7.5232E-06; 1000.00, 3.1281E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.8370E-01;
    out.liqThermalConductivity.B = -2.8100E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [175.47, 2.3439E-01; 337.85, 1.8876E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -7.7630E00;
    out.vapThermalConductivity.B = 1.0279E00;
    out.vapThermalConductivity.C = -7.4360E07;
    out.vapThermalConductivity.D = 6.7700E09;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.00, 1.3652E-02; 1000.0, 1.3927E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 3.5130E-02;
    out.SurfaceTension.B = -7.0400E-06;
    out.SurfaceTension.C = -1.2160E-07;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [273.10, 2.4138E-02; 503.15, 8.0358E-04];
    out.SurfaceTension.eqn = 100;    
    
end


