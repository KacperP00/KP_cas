function out = pc_1decanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1decanol','1-decanol','1-Decanol','n-decanol','n-Decanol', 'decyl alcohol'}; 
    out.CAS = '112-30-1';
    % components
    out.ChemicalFormular.C = 10;         
    out.ChemicalFormular.H = 22;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 158.284;  % kg/kmol    
    out.Tcrit = 687.30;             % K    
    out.pcrit = 2.3150E06;          % Pa 
    out.vcrit = 0.64900;            % m^3/kmol 
    out.CritCompressFactor = 0.263; % -
    
%% 2 
    out.MeltingPoint = 280.05;              % K
    out.TripletPointT = 280.05;             % K
    out.TripletPointp = 1.2019E-01;         % Pa
    out.NormalBoilingPoint = 504.07;        % K
    out.liqMolarVolume = 0.191854;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.9670E08;             % J/kmol
    out.IG_GibbsOfFormation = -9.8790E07;            % J/kmol
    out.IG_AbsoluteEntropy =  5.9700E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = 3.7656E07;           % J/kmol
    out.StandardNetHeatOfCombustion = -6.1159E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.2650E08;       % J/kmol
    out.StandardHeatOfFormation = -4.7810E08;        % J/kmol
    out.StandardAbsoluteEntropy = 4.1690E05;         % J/kmolK
    
%% 4 
    out.AcentricFactor = 0.6219;            % -
    out.RadiusOfGyration = 5.4990E-10;      % m
    out.SolubilityParameter = 2.02500E04;   % (J/m^3)^.5 
    out.DipoleMoment = 5.4000E-30;          % C^m
    out.VanDerWaalsVolume = 0.11378;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.5730E09;        % m^2/kmol
    out.RefractiveIndex = 1.4350;           % -
    out.FlashPoint = 355.37;                % K
    out.FlammabilityLimits = [0.7 5.5];     % vol-%
    out.FlammabilityLimitTemps = [376 416]; % K
    out.AutoignitionTemp = 561.15;          % K
    
%% 6 
    out.solDensity.A = 6.1900E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [280.05, 6.1900E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.7384E-01;
    out.liqDensity.B = 2.4241E-01;
    out.liqDensity.C = 6.8730E02;
    out.liqDensity.D = 2.6646E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [280.05, 5.2902E00; 687.30, 1.5422E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  2.5059E02;
    out.vapPressure.B = -1.9169E04;
    out.vapPressure.C = -3.2903E01;
    out.vapPressure.D =  1.4627E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [280.05, 1.2019E-01; 687.30, 2.3302E06];
    out.vapPressure.eqn = 101;
    
    out.HeatOfVap.A = 1.1750E08;
    out.HeatOfVap.B = 6.5112E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [280.05, 8.3570E07; 687.30, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 4.11026E05;
    out.solHeatCapacity.B = -2.4283E03;
    out.solHeatCapacity.C = 7.3148E00;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [180.00, 2.8017E05; 270.00, 3.5787E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 4.9885E06;
    out.liqHeatCapacity.B = -5.2898E04;
    out.liqHeatCapacity.C = 2.1635E02;
    out.liqHeatCapacity.D = -3.7538E-01;
    out.liqHeatCapacity.E = 2.3674E-04;
    out.liqHeatCapacity.range = [280.05, 3.5373E05; 504.07, 5.0206E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.6984E05;
    out.IG_HeatCapacity.B = 5.5920E05;
    out.IG_HeatCapacity.C = 1.5680E03;
    out.IG_HeatCapacity.D = 3.9380E05;
    out.IG_HeatCapacity.E = 7.2050E02;
    out.IG_HeatCapacity.range = [200.00, 1.8502E05; 1500.00, 6.2186E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.8836E-01;
    out.SecondViralCoef.B = -5.8491E02;
    out.SecondViralCoef.C = -1.1672E08;
    out.SecondViralCoef.D = 3.8936E21;
    out.SecondViralCoef.E = -1.5303E24;
    out.SecondViralCoef.range = [343.65, -6.9662E00; 1500.00, 6.3948E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -8.0656E01;
    out.liqViscosity.B = 6.3255E03;
    out.liqViscosity.C = 9.6460E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [285.00, 1.9500E-02; 504.07, 3.0879E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 5.5065E-08;
    out.vapViscosity.B = 8.3410E-01;
    out.vapViscosity.C = 7.9560E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [280.05, 4.7154E-06; 1000.00, 1.6215E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2800E-01;
    out.liqThermalConductivity.B = -2.2300E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [280.05, 1.6555E-01; 503.35, 1.1575E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -3.0720E-01;
    out.vapThermalConductivity.B = 4.8900E-01;
    out.vapThermalConductivity.C = -6.7500E04;
    out.vapThermalConductivity.D = -2.9400E07;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [504.07, 2.5906E-02; 1000.00, 9.3886E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1263E-02;
    out.SurfaceTension.B = 1.0395E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [280.05, 2.9754E-02; 687.30, 0.0000E00];
    out.SurfaceTension.eqn = 106;    
    
end


