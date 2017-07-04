function out = pc_nnonane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'nnonane','n-nonane','nonyl hydride'}; 
    out.CAS = '111-84-2';
    % components
    out.ChemicalFormular.C = 9;         
    out.ChemicalFormular.H = 20;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 128.258;   % kg/kmol    
    out.Tcrit = 594.60;              % K    
    out.pcrit = 2.2900E6;            % Pa 
    out.vcrit = 0.54400;             % m^3/kmol 
    out.CritCompressFactor = 0.252;  % -
    
%% 2 
    out.MeltingPoint = 219.66;             % K
    out.TripletPointT = 219.66;            % K
    out.TripletPointp = 4.3058E-01;        % Pa
    out.NormalBoilingPoint = 423.97;       % K
    out.liqMolarVolume  = 0.179559;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.2874E08;            % J/kmol
    out.IG_GibbsOfFormation = 2.4980E07;            % J/kmol
    out.IG_AbsoluteEntropy = 5.0640E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 1.5470E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -5.6846E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4435;            % -
    out.RadiusOfGyration = 4.8140E-10;      % m
    out.SolubilityParameter = 1.5600E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.09895;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.3690E09;        % m^2/kmol
    out.RefractiveIndex = 1.4031;           % -
    out.FlashPoint = 304.00;                % K
    out.FlammabilityLimits = [0.7 5.6];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 478.00;          % K
    
%% 6 
    out.solDensity.A = 7.4283E00;
    out.solDensity.B = -2.3699E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [193.15, 6.9706E00; 219.66, 6.9077E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 4.8387E-01;
    out.liqDensity.B = 2.6147E-01;
    out.liqDensity.C = 5.9460E02;
    out.liqDensity.D = 2.8281E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [219.66, 6.0070E00; 594.60, 1.8506E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.0935E02;
    out.vapPressure.B = -9.0304E03;
    out.vapPressure.C = -1.2882E01;
    out.vapPressure.D =  7.8544E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [219.66, 4.3058E-01; 594.60, 2.3054E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.0370E07;
    out.HeatOfVap.B = 3.8522E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [219.66, 5.0545E07; 594.60, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -3.1900E04;
    out.solHeatCapacity.B = 2.3720E03;
    out.solHeatCapacity.C = -1.2440E01;
    out.solHeatCapacity.D = 3.0180E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [20.00, 1.0805E04; 217.21, 2.0569E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 3.8308E05;
    out.liqHeatCapacity.B = -1.1398E03;
    out.liqHeatCapacity.C = 2.7101E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [219.66, 2.6346E05; 325.00, 2.9890E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.5175E05;
    out.IG_HeatCapacity.B = 4.9150E05;
    out.IG_HeatCapacity.C = 1.6448E03;
    out.IG_HeatCapacity.D = 3.4700E05;
    out.IG_HeatCapacity.E = 7.4960E02;
    out.IG_HeatCapacity.range = [200.00, 1.6257E05; 1500.00, 5.5407E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.9060E-01;
    out.SecondViralCoef.B = -4.6860E02;
    out.SecondViralCoef.C = -6.6800E07;
    out.SecondViralCoef.D = -7.8930E20;
    out.SecondViralCoef.E = 1.8200E23;
    out.SecondViralCoef.range = [297.83, -6.5900E00; 1500.00, 5.8381E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.1149E01;
    out.liqViscosity.B = 1.6580E03;
    out.liqViscosity.C = 1.4540E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [219.66, 3.1483E-03; 423.97, 2.1559E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.0344E-07;
    out.vapViscosity.B = 7.7301E-01;
    out.vapViscosity.C = 2.2047E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [219.66, 3.3347E-06; 1000.00, 1.7668E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0900E-01;
    out.liqThermalConductivity.B = -2.6400E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [219.66, 1.5101E-01; 423.97, 9.7072E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -6.5771E-02;
    out.vapThermalConductivity.B = 2.7198E-01;
    out.vapThermalConductivity.C = -3.4823E03;
    out.vapThermalConductivity.D = -1.5803E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [423.97, 2.1299E-02; 1000.00, 1.0597E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.4975E-02;
    out.SurfaceTension.B = 1.2897E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [219.66, 3.0331E-02; 594.60, 0];
    out.SurfaceTension.eqn = 106;    
    
end


