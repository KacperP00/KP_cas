function out = pc_2octanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2octanone','2-octanone','methyl hexyl ketone'}; 
    out.CAS = '111-13-7';
    % components
    out.ChemicalFormular.C = 8;         
    out.ChemicalFormular.H = 16;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 128.214;  % kg/kmol    
    out.Tcrit = 632.70;             % K    
    out.pcrit = 2.6400E06;          % Pa 
    out.vcrit = 0.49700;            % m^3/kmol 
    out.CritCompressFactor = 0.249; % -
    
%% 2 
    out.MeltingPoint = 252.85;              % K
    out.TripletPointT = 252.85;             % K
    out.TripletPointp = 4.6773E00;          % Pa
    out.NormalBoilingPoint = 446.15;        % K
    out.liqMolarVolume = 0.157287;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.2160E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.1380E08;            % J/kmol
    out.IG_AbsoluteEntropy =  4.9620E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.4419E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.6984E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.2917E08;              % J/kmol
    out.StandardHeatOfFormation = -3.7270E08;               % J/kmol
    out.StandardAbsoluteEntropy = 3.7384E05;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4549;            % -
    out.RadiusOfGyration = 4.8920E-10;      % m
    out.SolubilityParameter = 1.7580E4;    % (J/m^3)^.5 
    out.DipoleMoment = 8.2100E-30;                   % C^m
    out.VanDerWaalsVolume = 0.09019;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.2590E09;        % m^2/kmol
    out.RefractiveIndex = 1.4133;           % -
    out.FlashPoint = 324.15;                % K
    out.FlammabilityLimits = [1.0 6.1];     % vol-%
    out.FlammabilityLimitTemps = [326 363]; % K
    out.AutoignitionTemp = 653.00;          % K
    
%% 6 
    out.solDensity.A = 7.7800E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [252.85, 7.7800E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.0006E-01;
    out.liqDensity.B = 2.48511E-01;
    out.liqDensity.C = 6.3270E02;
    out.liqDensity.D = 2.9942E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [192.78, 6.1230E00; 582.80, 1.8484E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  6.3775E01;
    out.vapPressure.B = -7.7113E03;
    out.vapPressure.C = -5.7359E00;
    out.vapPressure.D =  3.0902E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [252.85, 4.6773E00; 632.70, 2.6466E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.5363E07;
    out.HeatOfVap.B = 3.8718E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [252.85, 5.3646E07; 632.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.9270E04;
    out.solHeatCapacity.B = 1.9170E03;
    out.solHeatCapacity.C = -8.2019E00;
    out.solHeatCapacity.D = 1.6844E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [13.00, 4.3019E03; 252.85, 2.1336E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.0040E05;
    out.liqHeatCapacity.B = -4.2620E02;
    out.liqHeatCapacity.C = 1.1172E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [252.86, 2.6406E05; 500.00, 3.6660E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.3901E05;
    out.IG_HeatCapacity.B = 3.8060E05;
    out.IG_HeatCapacity.C = 1.3717E03;
    out.IG_HeatCapacity.D = 2.2573E05;
    out.IG_HeatCapacity.E = 6.6096E02;
    out.IG_HeatCapacity.range = [150.00, 1.4162E05; 1500.00, 4.6547E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.3845E-01;
    out.SecondViralCoef.B = -3.9878E02;
    out.SecondViralCoef.C = -1.0424E08;
    out.SecondViralCoef.D = -5.3291E20;
    out.SecondViralCoef.E = 7.0168E22;
    out.SecondViralCoef.range = [316.35, -7.3161E00; 1500.00, 4.1692E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1736E01;
    out.liqViscosity.B = 1.4152E03;
    out.liqViscosity.C = 3.6180E-04;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [252.85, 2.1614E-03; 446.15, 1.9128E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 8.0901E-08;
    out.vapViscosity.B = 7.9062E-01;
    out.vapViscosity.C = 9.9338E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [252.85, 4.6110E-06; 1000.00, 1.7325E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1320E-01;
    out.liqThermalConductivity.B = -2.4940E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [252.85, 1.5014E-01; 499.00, 8.8749E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -2.0184E-03;
    out.vapThermalConductivity.B = 1.0027E00;
    out.vapThermalConductivity.C = -2.0406E04;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [446.15, 2.0463E-02; 1000.0, 1.0597E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.9062E-02;
    out.SurfaceTension.B = 1.2967E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [252.85, 3.0488E-02; 632.70, 0];
    out.SurfaceTension.eqn = 106;    
    
end


