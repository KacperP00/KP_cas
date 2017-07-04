function out = pc_1butanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1butanol','1-butanol','1-Butanol','n-butanol','n-Butanol, Butan-1-ol'}; 
    out.CAS = '71-36-3';
    % components
    out.ChemicalFormular.C = 4;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 74.123;   % kg/kmol    
    out.Tcrit = 563.05;             % K    
    out.pcrit = 4.4230E06;          % Pa 
    out.vcrit = 0.27400;            % m^3/kmol 
    out.CritCompressFactor = 0.260; % -
    
%% 2 
    out.MeltingPoint = 183.85;              % K
    out.TripletPointT = 184.51;             % K
    out.TripletPointp = 5.7220E-04;         % Pa
    out.NormalBoilingPoint = 390.81;        % K
    out.liqMolarVolume = 0.091991;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.7460E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.5030E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.6148E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  9.3720E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.4560E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5935;            % -
    out.RadiusOfGyration = 3.2250E-10;      % m
    out.SolubilityParameter = 2.3350E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.5400E-30;          % C^m
    out.VanDerWaalsVolume = 0.05240;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.6200E08;        % m^2/kmol
    out.RefractiveIndex = 1.3971;           % -
    out.FlashPoint = 302.00;                % K
    out.FlammabilityLimits = [1.4 11.2];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 616.00;          % K

%% 6 
    out.solDensity.A = 1.4100E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [184.51, 1.4100E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.6500E-01;
    out.liqDensity.B = 2.6660E-01;
    out.liqDensity.C = 5.6305E02;
    out.liqDensity.D = 2.4419E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [184.51, 1.2016E01; 563.05, 3.6197E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  9.3173E01;
    out.vapPressure.B = -9.1859E03;
    out.vapPressure.C = -9.7464E00;
    out.vapPressure.D =  4.7796E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [184.51, 5.7220E-04; 563.05, 4.3392E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.7390E07;
    out.HeatOfVap.B = 1.7300E-01;
    out.HeatOfVap.C = 2.9150E-01;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [184.51, 6.0575E07; 563.05, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.1680E04;
    out.solHeatCapacity.B = 1.0262E03;
    out.solHeatCapacity.C = -1.7010E00;
    out.solHeatCapacity.D = -1.5820E-02;
    out.solHeatCapacity.E = 6.6500E-05;
    out.solHeatCapacity.range = [20.00, 8.0477E03; 184.51, 9.7456E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.9120E05;
    out.liqHeatCapacity.B = -7.3040E02;
    out.liqHeatCapacity.C = 2.2998E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [184.51, 1.3473E05; 390.81, 2.5701E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 7.4540E04;
    out.IG_HeatCapacity.B = 2.5907E05;
    out.IG_HeatCapacity.C = 1.6073E03;
    out.IG_HeatCapacity.D = 1.7320E05;
    out.IG_HeatCapacity.E = 7.1240E02;
    out.IG_HeatCapacity.range = [200.00, 8.1616E04; 1500.00, 2.8509E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.8900E-01;
    out.SecondViralCoef.B = -1.8200E02;
    out.SecondViralCoef.C = -4.0500E07;
    out.SecondViralCoef.D = -2.2700E20;
    out.SecondViralCoef.E = 4.5400E22;
    out.SecondViralCoef.range = [282.00, 3.9128E00; 1500.00, 5.5659E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -3.5426E01;
    out.liqViscosity.B = 3.1845E03;
    out.liqViscosity.C = 3.2965E00;
    out.liqViscosity.D = -3.0000E-27;
    out.liqViscosity.E = 1.0000E01;
    out.liqViscosity.range = [190.00, 2.5440E-01; 390.81, 3.8873E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.4031E-06;
    out.vapViscosity.B = 4.6110E-01;
    out.vapViscosity.C = 5.3700E+02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [184.51, 3.9786E-06; 1000.00, 2.2066E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1360E-01;
    out.liqThermalConductivity.B = -2.0340E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [184.51, 1.7607E-01; 390.81, 1.3411E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -4.4940E-02;
    out.vapThermalConductivity.B = 4.4600E-02;
    out.vapThermalConductivity.C = -1.3552E+03;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [370.00, 2.1971E-02; 800.00, 8.7247E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 4.9830E-02;
    out.SurfaceTension.B = -8.5400E-05;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [273.15, 2.6503E-02; 413.15, 1.4547E-02];
    out.SurfaceTension.eqn = 100;    
    
end


