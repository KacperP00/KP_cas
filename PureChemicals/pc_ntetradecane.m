function out = pc_ntetradecane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'ntetradecane','n-tetradecane','tetradecyl hydride'}; 
    out.CAS = '629-59-4';
    % components
    out.ChemicalFormular.C = 14;         
    out.ChemicalFormular.H = 30;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 198.392;  % kg/kmol    
    out.Tcrit = 693.00;             % K    
    out.pcrit = 1.5700E06;          % Pa 
    out.vcrit = 0.83000;            % m^3/kmol 
    out.CritCompressFactor = 0.226; % -
    
%% 2 
    out.MeltingPoint = 279.01;              % K
    out.TripletPointT = 279.01;             % K
    out.TripletPointp = 2.5269E-01;         % Pa
    out.NormalBoilingPoint = 526.73;        % K
    out.liqMolarVolume = 0.261271;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.3244E08;            % J/kmol
    out.IG_GibbsOfFormation = 6.5990E07;            % J/kmol
    out.IG_AbsoluteEntropy = 7.0259E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = 4.5070E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -8.7328E09;   % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.6430;            % -
    out.RadiusOfGyration = 6.4270E-10;      % m
    out.SolubilityParameter = 1.6100E04;    % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.15010;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 2.0440E09;        % m^2/kmol
    out.RefractiveIndex = 1.4269;           % -
    out.FlashPoint = 373.15;                % K
    out.FlammabilityLimits = [0.5 4.7];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 473.00;          % K
    
%% 6 
    out.solDensity.A = 4.6071E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [0, 4.6071E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.0545E-01;
    out.liqDensity.B = 2.5350E-01;
    out.liqDensity.C = 6.9300E02;
    out.liqDensity.D = 3.0538E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [279.01, 3.8924E00; 693.00, 1.2049E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.4047E02;
    out.vapPressure.B = -1.3231E04;
    out.vapPressure.C = -1.6859E01;
    out.vapPressure.D =  6.5877E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [279.01, 2.5268E-01; 693.00, 1.5693E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 9.0539E07;
    out.HeatOfVap.B = 4.4467E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [279.01, 7.2002E07; 693.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -4.9782E04;
    out.solHeatCapacity.B = 3.3010E03;
    out.solHeatCapacity.C = -1.5807E01;
    out.solHeatCapacity.D = 3.6919E-02;
    out.solHeatCapacity.E = -1.3660E-05;
    out.solHeatCapacity.range = [29.98, 3.5958E04; 279.01, 3.5981E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.5314E05;
    out.liqHeatCapacity.B = 2.9130E01;
    out.liqHeatCapacity.C = 8.6116E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [279.01, 4.2831E05; 526.73, 6.0741E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.3082E05;
    out.IG_HeatCapacity.B = 7.8678E05;
    out.IG_HeatCapacity.C = 1.6823E03;
    out.IG_HeatCapacity.D = 5.4486E05;
    out.IG_HeatCapacity.E = 7.4310E02;
    out.IG_HeatCapacity.range = [200.00, 2.4864E05; 1500.00, 8.6225E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.8209E-01;
    out.SecondViralCoef.B = -5.0132E02;
    out.SecondViralCoef.C = -2.8594E08;
    out.SecondViralCoef.D = -2.5720E21;
    out.SecondViralCoef.E = 4.5523E23;
    out.SecondViralCoef.range = [346.50, -1.3893E01; 1500.00, 6.3065E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.0486E01;
    out.liqViscosity.B = 2.0884E03;
    out.liqViscosity.C = 1.2852E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [279.01, 3.1395E-03; 528.15, 2.0874E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 4.4565E-08;
    out.vapViscosity.B = 8.6840E-01;
    out.vapViscosity.C = 2.2816E02;
    out.vapViscosity.D = -4.3472E03;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [279.01, 3.3635E-06; 1000.00, 1.4672E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.02936E-01;
    out.liqThermalConductivity.B = -2.1798E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [279.01, 1.4211E-01; 526.73, 8.8113E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -1.6362E02;
    out.vapThermalConductivity.B = 9.1930E-01;
    out.vapThermalConductivity.C = -1.0876E09;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [526.73, 2.5172E-02; 1000.0, 8.6152E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.6276E-02;
    out.SurfaceTension.B = 1.3631E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [279.01, 2.7883E-02; 693.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


