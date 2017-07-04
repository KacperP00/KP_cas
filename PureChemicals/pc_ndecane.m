function out = pc_ndecane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'ndecane','n-decane','decyl hydride'}; 
    out.CAS = '124-18-5';
    % components
    out.ChemicalFormular.C = 10;         
    out.ChemicalFormular.H = 22;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 142.285;   % kg/kmol    
    out.Tcrit = 617.70;              % K    
    out.pcrit = 2.1100E6;            % Pa 
    out.vcrit = 0.60000;             % m^3/kmol 
    out.CritCompressFactor =0.247;   % -
    
%% 2 
    out.MeltingPoint = 243.51;             % K
    out.TripletPointT = 243.51;            % K
    out.TripletPointp = 1.3930;            % Pa
    out.NormalBoilingPoint = 447.30;       % K
    out.liqMolarVolume  = 0.195827;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.4946E8;             % J/kmol
    out.IG_GibbsOfFormation = 3.3180E07;            % J/kmol
    out.IG_AbsoluteEntropy = 5.4570E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 2.8710E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -6.2942E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4923;            % -
    out.RadiusOfGyration = 5.1480E-10;      % m
    out.SolubilityParameter = 1.5700E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.10920;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.5040E09;        % m^2/kmol
    out.RefractiveIndex = 1.4097;           % -
    out.FlashPoint = 319.00;                % K
    out.FlammabilityLimits = [0.7 5.4];     % vol-%
    out.FlammabilityLimitTemps = [284 313]; % K
    out.AutoignitionTemp = 474.00;          % K
    
%% 6 
    out.solDensity.A = 6.8978E00;
    out.solDensity.B = -2.4600E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [193.15, 6.4227; 243.15, 6.2988E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 4.2831E-01;
    out.liqDensity.B = 2.5745E-01;
    out.liqDensity.C = 6.1770E02;
    out.liqDensity.D = 2.8912E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [243.15, 5.3811E00; 617.70, 1.6637E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.1273E02;
    out.vapPressure.B = -9.7496E03;
    out.vapPressure.C = -1.3245E01;
    out.vapPressure.D =  7.1266E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [243.51, 1.3930E00; 617.70, 2.0908E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.6126E07;
    out.HeatOfVap.B = 3.9797E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [243.51, 5.4168E07; 617.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -3.4900E04;
    out.solHeatCapacity.B = 2.5280E03;
    out.solHeatCapacity.C = -1.2650E01;
    out.solHeatCapacity.D = 2.8480E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [19.93, 1.0684E04; 240.00, 2.3689E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.7862E05;
    out.liqHeatCapacity.B = -1.9791E02;
    out.liqHeatCapacity.C = 1.0737E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [243.51, 2.9409E05; 460.00, 4.1478E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.6720E05;
    out.IG_HeatCapacity.B = 5.3530E05;
    out.IG_HeatCapacity.C = 1.6141E03;
    out.IG_HeatCapacity.D = 3.7820E05;
    out.IG_HeatCapacity.E = 7.4200E02;
    out.IG_HeatCapacity.range = [200.00, 1.7967E05; 1500.00, 6.0932E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.800E-01;
    out.SecondViralCoef.B = -5.8850E02;
    out.SecondViralCoef.C = -7.6060E07;
    out.SecondViralCoef.D = -1.6130E21;
    out.SecondViralCoef.E = 3.9940E23;
    out.SecondViralCoef.range = [309.23, -7.8395E00; 1500.00, 6.5078E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.6468E01;
    out.liqViscosity.B = 1.5335E03;
    out.liqViscosity.C = 7.5110E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [243.51, 2.3742E-03; 448.15, 2.1165E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.6400E-08;
    out.vapViscosity.B = 9.4870E-01;
    out.vapViscosity.C = 7.1000E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [243.51, 3.7547E-06; 1000.00, 1.7295E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0630E-01;
    out.liqThermalConductivity.B = -2.5000E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [243.51, 1.4542E-01; 447.30, 9.4475E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -6.6840E02;
    out.vapThermalConductivity.B = 9.3230E-01;
    out.vapThermalConductivity.C = -4.0710E09;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [447.30, 2.1731E-02; 1000.00, 1.0286E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.5435E-02;
    out.SurfaceTension.B = 1.3095E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [243.51, 2.8756E-02; 617.70, 0];
    out.SurfaceTension.eqn = 106;    
    
end


