function out = pc_3petanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'3petanol','3-petanol','pentanol-3','diethyl carbinol'}; 
    out.CAS = '584-02-1';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 88.150;  % kg/kmol    
    out.Tcrit = 559.60;             % K    
    out.pcrit = 3.8800E6;           % Pa 
    out.vcrit = 0.32500;            % m^3/kmol 
    out.CritCompressFactor = 0.271; % -
    
%% 2 
    out.MeltingPoint = 204.15;              % K
    out.TripletPointT = 204.15;             % K
    out.TripletPointp = 1.4819E-02;         % Pa
    out.NormalBoilingPoint = 388.45;        % K
    out.liqMolarVolume = 0.108015;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.1673E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.5816E8;            % J/kmol
    out.IG_AbsoluteEntropy =  3.8200E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  9.0800E6;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.0483E9;    % J/kmol
    out.StandardGibbsOfFormation = -1.7050E08;      % J/kmol
    out.StandardHeatOfFormation = -3.6890E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.4910E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5383;            % -
    out.RadiusOfGyration = 3.2530E-10;      % m
    out.SolubilityParameter = 2.1780E4;     % (J/m^3)^.5 
    out.DipoleMoment = 5.4700E-30;          % C^m
    out.VanDerWaalsVolume = 0.06262;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.9700E08;        % m^2/kmol
    out.RefractiveIndex = 1.4079;           % -
    out.FlashPoint = 313.15;                % K
    out.FlammabilityLimits = [1.2 9.0];     % vol-%
    out.FlammabilityLimitTemps = [300 332]; % K
    out.AutoignitionTemp = 708.15;          % K
    
%% 6 
    out.solDensity.A = 1.1900E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [204.15, 1.1900E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.0198E-01;
    out.liqDensity.B = 2.6069E-01;
    out.liqDensity.C = 5.5960E02;
    out.liqDensity.D = 2.6161E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [204.15, 1.0152E01; 559.60, 3.0764E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.3775E02;
    out.vapPressure.B = -1.1249E04;
    out.vapPressure.C = -1.6331E01;
    out.vapPressure.D =  2.8016E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [204.15, 1.4819E-02; 559.60, 3.9083E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 8.4954E07;
    out.HeatOfVap.B = 6.0261E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [204.15, 6.4626E07; 559.60, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 0;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [NaN,NaN; NaN,NaN ];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = -1.6787E06;
    out.liqHeatCapacity.B = 1.1243E04;
    out.liqHeatCapacity.C = -1.6088E01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [303.02, 2.5094E05; 326.16, 2.7687E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0751E05;
    out.IG_HeatCapacity.B = 3.0130E05;
    out.IG_HeatCapacity.C = 1.8100E03;
    out.IG_HeatCapacity.D = 2.1755E05;
    out.IG_HeatCapacity.E = 8.1170E02;
    out.IG_HeatCapacity.range = [300, 1.3596E05; 1200.00, 3.2093E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.2585E-01;
    out.SecondViralCoef.B = -2.3629E02;
    out.SecondViralCoef.C = -3.0124E07;
    out.SecondViralCoef.D = 3.5131E20;
    out.SecondViralCoef.E = -1.5281E23;
    out.SecondViralCoef.range = [279.80, -7.1804E00; 1500.00, 5.9407E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.5290E02;
    out.liqViscosity.B = 1.4984E04;
    out.liqViscosity.C = 3.4620E01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [250.00, 1.6349E-01; 388.45, 3.6484E-04];
    out.liqViscosity.eqn = 101;
    
    %%Ab hier Scan abgeschnitten, Koeffizienten bereits richtig eingetragen aber Grenzen noch nicht.
    out.vapViscosity.A = 1.1891E-07;
    out.vapViscosity.B = 7.6401E-01;
    out.vapViscosity.C = 1.0337E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [NaN, NaN; NaN, NaN];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1707E-01;
    out.liqThermalConductivity.B = -2.7901E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [NaN, NaN; NaN, NaN];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 3.1596E03;
    out.vapThermalConductivity.B = 9.0260E-01;
    out.vapThermalConductivity.C = 1.2930E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [NaN, NaN; NaN, NaN];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.8608E-02;
    out.SurfaceTension.B = 1.1579E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [NaN, NaN; NaN, NaN];
    out.SurfaceTension.eqn = 106;    
    
end


