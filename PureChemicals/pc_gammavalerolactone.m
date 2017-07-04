function out = pc_gammavalerolactone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'gammavalerolactone','gamma-pentalactone','4-pentanolide'}; 
    out.CAS = '108-29-2';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 8;
    out.ChemicalFormular.O = 2;
%% 1  
    out.MolecularWeight = 100.117;  % kg/kmol    
    out.Tcrit = 727.00;             % K    
    out.pcrit = 4.8400E6;           % Pa 
    out.vcrit = 0.21900;            % m^3/kmol 
    out.CritCompressFactor = 0.223; % -
    
%% 2 
    out.MeltingPoint = 242.15;              % K
    out.TripletPointT = 242.15;             % K
    out.TripletPointp = 3.1332E-01;         % Pa
    out.NormalBoilingPoint = 480.65;        % K
    out.liqMolarVolume = 0.095430;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -4.1800E08;             % J/kmol
    out.IG_GibbsOfFormation = -2.9700E07;            % J/kmol
    out.IG_AbsoluteEntropy =  3.5000E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = NaN;                 % J/kmol
    out.StandardNetHeatOfCombustion = -2.4929E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4030;            % -
    out.RadiusOfGyration = 3.4319E-10;      % m
    out.SolubilityParameter = 2.287E04;     % (J/m^3)^.5 
    out.DipoleMoment = 1.4343E-29;          % C^m
    out.VanDerWaalsVolume = 0.05631;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.5900E08;        % m^2/kmol
    out.RefractiveIndex = 1.4330;           % -
    out.FlashPoint = 354.00;                % K
    out.FlammabilityLimits = [NaN 10.6];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = NaN;             % K
    
%% 6 
    out.solDensity.A = 1.2900E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [242.15, 1.2900E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.8680E-01;
    out.liqDensity.B = 2.7532E-01;
    out.liqDensity.C = 7.2700E02;
    out.liqDensity.D = 3.4900E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [242.15, 1.0983E01; 727.00, 3.5842E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  8.3640E01;
    out.vapPressure.B = -8.7843E03;
    out.vapPressure.C = -8.8878E00;
    out.vapPressure.D =  4.5338E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [242.15, 3.1332E-01; 727.00, 4.8391E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.1320E07;
    out.HeatOfVap.B = 2.9680E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [242.15, 5.4373E07; 727.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 1.4830E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [242.15, 1.4830E05; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 9.5800E04;
    out.liqHeatCapacity.B = 2.7750E02;
    out.liqHeatCapacity.C = 0;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [363.50, 1.9667E05; 543.50, 2.4662E05];
    out.liqHeatCapacity.eqn = 107;    
    
    out.IG_HeatCapacity.A = 8.9700E04;
    out.IG_HeatCapacity.B = 1.8340E05;
    out.IG_HeatCapacity.C = 6.7500E02;
    out.IG_HeatCapacity.D = 9.3000E04;
    out.IG_HeatCapacity.E = 2.2500E03;
    out.IG_HeatCapacity.range = [300.00 , 1.3190E05; 1500.00, 2.9902E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.1157E-01;
    out.SecondViralCoef.B = -2.9963E02;
    out.SecondViralCoef.C = -8.9927E07;
    out.SecondViralCoef.D = -2.9311E20;
    out.SecondViralCoef.E = -2.6964E22;
    out.SecondViralCoef.range = [363.50, -3.6900E00; 1500.00, -1.4841E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.8539E00;
    out.liqViscosity.B = 5.6334E02;
    out.liqViscosity.C = -1.1445E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [242.15, 1.1024E-03; 572.15, 1.0770E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 8.0120E-08;
    out.vapViscosity.B = 8.5930E-01;
    out.vapViscosity.C = 6.1350E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range =[242.15, 7.1501E-06; 1000, 2.8562E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2568E-01;
    out.liqThermalConductivity.B = -2.1770E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [242.15, 1.7296E-01; 480.65, 1.2104E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 3.8300E-04;
    out.vapThermalConductivity.B = 8.7240E-01;
    out.vapThermalConductivity.C = 7.0550E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [480.65, 3.3925E-02; 1000.00, 1.5807E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.6815E-02;
    out.SurfaceTension.B = 1.3153E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [242.15, 3.9217E-02; 727.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


