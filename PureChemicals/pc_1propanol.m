function out = pc_1propanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1propanol','1-propanol','1Ppropanol','n-propanol','n-Propanol', 'propanol-1'}; 
    out.CAS = '71-23-8';
    % components
    out.ChemicalFormular.C = 3;         
    out.ChemicalFormular.H = 8;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 60.096;   % kg/kmol    
    out.Tcrit = 536.78;             % K    
    out.pcrit = 5.1750E+06;         % Pa 
    out.vcrit = 0.21900;            % m^3/kmol 
    out.CritCompressFactor = 0.254; % -
    
%% 2 
    out.MeltingPoint = 146.95;              % K
    out.TripletPointT = 146.95;             % K
    out.TripletPointp = 3.0828E-07;         % Pa
    out.NormalBoilingPoint = 370.35;        % K
    out.liqMolarVolume = 0.075089;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.5520E+08;            % J/kmol
    out.IG_GibbsOfFormation = -1.5990E+08;           % J/kmol
    out.IG_AbsoluteEntropy =  3.2247E+05;            % J/kmolK
    out.HeatFusionAtMeltPoint = 5.3720E+06;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.8438E+09;   % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK
    
%% 4 
    out.AcentricFactor = 0.6218;            % -
    out.RadiusOfGyration = 2.359E-10;       % m
    out.SolubilityParameter = 2.4450E+04;   % (J/m^3)^.5 
    out.DipoleMoment = 5.6000E-30;          % C^m
    out.VanDerWaalsVolume = 0.04217;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 6.2800E+08;       % m^2/kmol
    out.RefractiveIndex = 1.3837;           % -
    out.FlashPoint = 288.15;                % K
    out.FlammabilityLimits = [2.0 12.0];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 644.00;          % K
    
%% 6 
    out.solDensity.A = 1.7800E+01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [146.95, 1.7800E+01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.2350E+00;
    out.liqDensity.B = 2.7136E-01;
    out.liqDensity.C = 5.3678E+02;
    out.liqDensity.D = 2.4000E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [146.95, 1.5231E+01; 536.78, 4.5511E+00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  8.8134E+01;
    out.vapPressure.B = -8.4986E+03;
    out.vapPressure.C = -9.0766E+00;
    out.vapPressure.D =  8.3303E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [146.95, 3.0828E-07; 536.78, 5.1214E+06];
    out.vapPressure.eqn = 101;
    
    out.HeatOfVap.A = 6.3300E+07;
    out.HeatOfVap.B = 3.5750E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [146.95, 5.6460E+07; 536.78, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 8.7700E+03;
    out.solHeatCapacity.B = 8.1580E+02;
    out.solHeatCapacity.C = 2.6180E+00;
    out.solHeatCapacity.D = -7.6720E-02;
    out.solHeatCapacity.E = 3.3530E-04;
    out.solHeatCapacity.range = [15.00, 3.8141E+03; 140.00, 7.5044E+04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.5876E+05;
    out.liqHeatCapacity.B = -6.3500E+02;
    out.liqHeatCapacity.C = 1.9690E+00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [146.95, 1.0797E+05; 400.00 2.1980E+05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 6.1900E+04;
    out.IG_HeatCapacity.B = 2.0213E+05;
    out.IG_HeatCapacity.C = 1.6293E+03;
    out.IG_HeatCapacity.D = 1.2956E+05;
    out.IG_HeatCapacity.E = 7.2740E+02;
    out.IG_HeatCapacity.range = [200.00, 6.6651E+04; 1500.00, 2.2458E+05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.6200E-02;
    out.SecondViralCoef.B = -5.5600E+01;
    out.SecondViralCoef.C = -3.3370E+07;
    out.SecondViralCoef.D = -3.2000E+18;
    out.SecondViralCoef.E = -2.5610E+22;
    out.SecondViralCoef.range = [295.19, -2.9998E+00; 1500.00, -7.5487E-04];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.7528E02;
    out.liqViscosity.B = 6.1044E03;
    out.liqViscosity.C = 2.8761E01;
    out.liqViscosity.D = -5.1323E-02;
    out.liqViscosity.E = 1.0000E00;
    out.liqViscosity.range = [146.95, 9.3820E00; 370.35, 4.5317E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 7.9420E-07;
    out.vapViscosity.B = 5.4910E-01;
    out.vapViscosity.C = 4.1580E+02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [200.00, 4.7317E-06; 1000.00, 2.4902E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2030E-01;
    out.liqThermalConductivity.B = -2.1550E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [200.00, 1.7720E-01; 370.35, 1.4049E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -6.1384E+02;
    out.vapThermalConductivity.B = 7.9270E-01;
    out.vapThermalConductivity.C = -1.1574E+09;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [372.13, 2.1531E-02; 720.25, 7.0337E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 4.5330E-02;
    out.SurfaceTension.B = -6.8800E-05;
    out.SurfaceTension.C = -1.6000E-08;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [283.15, 2.4566E-02; 370.35, 1.7655E-02];
    out.SurfaceTension.eqn = 100;    
    
end


