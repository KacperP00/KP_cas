function out = pc_3pentanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'3pentanone','3-pentanone','metacetone','ethyl ketone'}; 
    out.CAS = '96-22-0';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 86.134;   % kg/kmol    
    out.Tcrit = 560.95;             % K    
    out.pcrit = 3.7400E6;           % Pa 
    out.vcrit = 0.33600;            % m^3/kmol 
    out.CritCompressFactor = 0.269; % -
    
%% 2 
    out.MeltingPoint = 234.18;              % K
    out.TripletPointT = 234.18;             % K
    out.TripletPointp = 7.3422E01;          % Pa
    out.NormalBoilingPoint = 375.14;        % K
    out.liqMolarVolume = 0.106410;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.5790E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.3440E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.7000E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.1590E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.8804E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.4200E08;       % J/kmol
    out.StandardHeatOfFormation = -2.9651E08;        % J/kmol
    out.StandardAbsoluteEntropy = 2.6600E05;         % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3448;            % -
    out.RadiusOfGyration = 3.5830E-10;      % m
    out.SolubilityParameter = 1.8410E04;    % (J/m^3)^.5 
    out.DipoleMoment = 9.4100E-30;          % C^m
    out.VanDerWaalsVolume = 0.05950;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.5400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3900;           % -
    out.FlashPoint = 286.00;                % K
    out.FlammabilityLimits = [1.5 8.0];     % vol-%
    out.FlammabilityLimitTemps = [277 308]; % K
    out.AutoignitionTemp = 725.00;          % K
    
%% 6 
    out.solDensity.A = 1.1800E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [234.18, 1.1800E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 7.1811E-01;
    out.liqDensity.B = 2.4129E-01;
    out.liqDensity.C = 5.6095E02;
    out.liqDensity.D = 2.7996E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [234.18, 1.0102E01; 560.95, 2.9761E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A = 4.4286E01;
    out.vapPressure.B = -5.4151E03;
    out.vapPressure.C = -3.0913E00;
    out.vapPressure.D =  1.8580E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [234.18, 7.3422E01; 560.95, 3.6993E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.2359E07;
    out.HeatOfVap.B = 4.0465E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [234.18, 4.2075E07; 560.95, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -9.2023E03;
    out.solHeatCapacity.B = 1.1108E03;
    out.solHeatCapacity.C = 5.6357E00;
    out.solHeatCapacity.D = -1.4318E-01;
    out.solHeatCapacity.E = 6.8742E-04;
    out.solHeatCapacity.range = [11.16, 3.7078E03; 110.39, 9.1569E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.9302E05;
    out.liqHeatCapacity.B = -1.7643E02;
    out.liqHeatCapacity.C = 5.6690E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [234.18, 1.8279E05; 375.14, 2.0661E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 9.6896E04;
    out.IG_HeatCapacity.B = 2.4907E05;
    out.IG_HeatCapacity.C = 1.4177E03;
    out.IG_HeatCapacity.D = 1.3010E05;
    out.IG_HeatCapacity.E = 6.4670E02;
    out.IG_HeatCapacity.range = [200.00, 1.0536E05; 1500.00, 3.0358E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.9680E-01;
    out.SecondViralCoef.B = -2.1873E02;
    out.SecondViralCoef.C = -4.5060E07;
    out.SecondViralCoef.D = -1.5690E20;
    out.SecondViralCoef.E = 1.9653E22;
    out.SecondViralCoef.range = [280.47,-4.8930E00; 1500.00, 3.7623E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.8695E00;
    out.liqViscosity.B = 5.9632E02;
    out.liqViscosity.C = -1.2025E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [234.18,1.0240E-03; 375.14, 2.2319E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.1640E-07;
    out.vapViscosity.B = 7.6150E-01;
    out.vapViscosity.C = 1.0794E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [234.18, 5.0785E-06; 1000.00, 2.0227E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1569E-01;
    out.liqThermalConductivity.B = -2.4081E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [234.18, 1.5930E-01; 375.14, 1.2535E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.2775E01;
    out.vapThermalConductivity.B = 1.0019E00;
    out.vapThermalConductivity.C = 1.9100E08;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.00, 8.9821E-03; 1000.00, 1.2082E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.0021E-02;
    out.SurfaceTension.B = 1.1704E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [234.18, 3.1888E-02; 560.95, 0];
    out.SurfaceTension.eqn = 106;    
    
end


