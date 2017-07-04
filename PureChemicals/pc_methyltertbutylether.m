function out = pc_nbutane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'methyltertbutylether','methyl tert-buty lether','propane'}; 
    out.CAS = '1634-04-4';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 88.150;    % kg/kmol    
    out.Tcrit = 497.10;              % K    
    out.pcrit = 3.4300E6;            % Pa 
    out.vcrit = 0.32900;             % m^3/kmol 
    out.CritCompressFactor =0.273;  % -
    
%% 2 
    out.MeltingPoint = 164.55;             % K
    out.TripletPointT = 164.55;            % K
    out.TripletPointp = 5.3566E-01;        % Pa
    out.NormalBoilingPoint = 328.35;       % K
    out.liqMolarVolume  = 0.119887;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.8350E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.1750E8;            % J/kmol
    out.IG_AbsoluteEntropy =  3.5780E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  7.6000E06;         % J/kmol
    out.StandardNetHeatOfCombustion = -3.1049E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2661;            % -
    out.RadiusOfGyration = 3.1790E-10;      % m
    out.SolubilityParameter = 1.5070E4;     % (J/m^3)^.5 
    out.DipoleMoment = 4.5400E-30;          % C^m
    out.VanDerWaalsVolume = 0.06171;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 9.0800E08;        % m^2/kmol
    out.RefractiveIndex = 1.3663;           % -
    out.FlashPoint = 245.00;                % K
    out.FlammabilityLimits = [2.0 15.1];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 733.15;          % K
    
%% 6 
    out.solDensity.A = 1.1400E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [164.55, 1.1400E01; 90.67, 3.0309E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.2157E-01;
    out.liqDensity.B = 2.7032E-01;
    out.liqDensity.C = 4.9710E02;
    out.liqDensity.D = 2.8290E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [164.55, 9.7682E00; 497.10, 3.0392E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.5875E01;
    out.vapPressure.B = -5.1316E03;
    out.vapPressure.C = -4.9604E00;
    out.vapPressure.D =  1.9123E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [164.55, 5.3566E-01; 497.10 , 3.4106E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.2024E07;
    out.HeatOfVap.B = 3.7826E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [164.55, 3.6096E07; 497.10, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.2579E04;
    out.solHeatCapacity.B = 1.4658E03;
    out.solHeatCapacity.C = -7.7409E00;
    out.solHeatCapacity.D = 9.6626E-03;
    out.solHeatCapacity.E = 9.9006E-05;
    out.solHeatCapacity.range = [11.50, 3.2704E03; 164.55, 1.3466E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.4012E05;
    out.liqHeatCapacity.B = -9.0000E00;
    out.liqHeatCapacity.C = 5.6300E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [164.55, 1.5388E05; 328.35, 1.9786E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 9.9328E04;
    out.IG_HeatCapacity.B = 3.0667E05;
    out.IG_HeatCapacity.C = 1.7426E03;
    out.IG_HeatCapacity.D = 2.0764E05;
    out.IG_HeatCapacity.E = 7.9559E02;
    out.IG_HeatCapacity.range = [200.00, 1.0394E05; 1500.00, 3.4321E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.8083E-01;
    out.SecondViralCoef.B = -1.8999E02;
    out.SecondViralCoef.C = -2.2566E07;
    out.SecondViralCoef.D = -9.3110E18;
    out.SecondViralCoef.E = 9.3054E18;
    out.SecondViralCoef.range = [248.60, -2.8896E00; 1500.00, 4.6565E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -8.4217E00;
    out.liqViscosity.B = 8.5730E02;
    out.liqViscosity.C = -4.3178E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [180.00, 2.7363E-03; 449.93, 1.0579E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.5439E-07;
    out.vapViscosity.B = 7.3601E-01;
    out.vapViscosity.C = 1.0824E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [164.55, 3.9839E-06; 1000.00, 2.2491E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2555E-01;
    out.liqThermalConductivity.B = -3.2780E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [164.55, 1.7161E-01; 328.35, 1.1792E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 3.1710E-04;
    out.vapThermalConductivity.B = 8.8303E-01;
    out.vapThermalConductivity.C = 4.6940E02;
    out.vapThermalConductivity.D = 8.6314E04;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.00, 1.1584E-02; 1000.00, 9.0857E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.9810E-02;
    out.SurfaceTension.B = 1.2352E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [164.55, 3.6402E-02; 497.10, 0];
    out.SurfaceTension.eqn = 106;    
    
end


