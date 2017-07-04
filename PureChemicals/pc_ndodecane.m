function out = pc_ndodecane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'ndodecane','n-dodecane','duodecane','adakane 12'}; 
    out.CAS = '112-40-3';
    % components
    out.ChemicalFormular.C = 12;         
    out.ChemicalFormular.H = 26;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 170.338;   % kg/kmol    
    out.Tcrit = 658.00;              % K    
    out.pcrit = 1.8200E6;            % Pa 
    out.vcrit = 0.71600;             % m^3/kmol 
    out.CritCompressFactor = 0.238;  % -
    
%% 2 
    out.MeltingPoint = 263.57;             % K
    out.TripletPointT = 263.57;            % K
    out.TripletPointp = 6.1520E-01;        % Pa
    out.NormalBoilingPoint = 489.47;       % K
    out.liqMolarVolume  = 0.228605;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.9072E08;            % J/kmol
    out.IG_GibbsOfFormation = 3.3180E07;            % J/kmol
    out.IG_AbsoluteEntropy = 6.2415E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 3.6840E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -7.5137E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5764;            % -
    out.RadiusOfGyration = 5.9140E-10;      % m
    out.SolubilityParameter = 1.5900E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.12960;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.7740E09;        % m^2/kmol
    out.RefractiveIndex = 1.4151;           % -
    out.FlashPoint = 347.00;                % K
    out.FlammabilityLimits = [0.6 4.7];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 476.00;          % K
    
%% 6 
    out.solDensity.A = 5.2800E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [263.57, 5.2800E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 3.5541E-01;
    out.liqDensity.B = 2.5511E-01;
    out.liqDensity.C = 6.5800E02;
    out.liqDensity.D = 2.9368E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [263.57, 4.5132E00; 658.00, 1.3932E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.3747E02;
    out.vapPressure.B = -1.1976E04;
    out.vapPressure.C = -1.6698E01;
    out.vapPressure.D =  8.0906E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [263.57, 6.1534E-01; 658.00, 1.8223E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 7.7337E07;
    out.HeatOfVap.B = 4.0681E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [263.57, 6.2802E07; 658.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -4.0970E04;
    out.solHeatCapacity.B = 2.8390E03;
    out.solHeatCapacity.C = -1.3232E01;
    out.solHeatCapacity.D = 2.7972E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [29.97, 3.2983E04; 258.41, 2.9175E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 5.0821E05;
    out.liqHeatCapacity.B = -1.3687E03;
    out.liqHeatCapacity.C = 3.1015E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [263.57, 3.6292E05; 330.00, 3.9429E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.1295E05;
    out.IG_HeatCapacity.B = 6.6330E05;
    out.IG_HeatCapacity.C = 1.7155E03;
    out.IG_HeatCapacity.D = 4.5161E05;
    out.IG_HeatCapacity.E = 7.7750E02;
    out.IG_HeatCapacity.range = [200.00, 2.2442E05; 1500.00, 7.4325E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 8.800E-01;
    out.SecondViralCoef.B = -1.0910E03;
    out.SecondViralCoef.C = -5.0300E07;
    out.SecondViralCoef.D = -5.4871E21;
    out.SecondViralCoef.E = 1.4959E24;
    out.SecondViralCoef.range = [329.00, -1.0701E01; 1500.00, 1.3759E-01];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.0607E01;
    out.liqViscosity.B = 1.9430E03;
    out.liqViscosity.C = 1.3205E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [263.57, 2.8109E-03; 489.47, 2.1196E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 6.3440E-08;
    out.vapViscosity.B = 8.2870E-01;
    out.vapViscosity.C = 2.1950E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [263.57, 3.5111E-06; 1000.00, 1.5932E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0470E-01;
    out.liqThermalConductivity.B = -2.3260E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [263.57, 1.4339E-01; 489.47, 9.0849E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 5.7190E-06;
    out.vapThermalConductivity.B = 1.4699E00;
    out.vapThermalConductivity.C = 5.7940E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [489.47, 2.3537E-02; 1000.00, 9.3010E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.5493E-02;
    out.SurfaceTension.B = 1.3262E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [263.57, 2.8150E-02; 658.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


