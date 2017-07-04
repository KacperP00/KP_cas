function out = pc_nbutylformate

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'nbutylformate','n-butyl formate','butyl methanoate'}; 
    out.CAS = '592-84-7';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 2;
%% 1  
    out.MolecularWeight = 102.133;   % kg/kmol    
    out.Tcrit = 559.0;               % K    
    out.pcrit = 3.5100E6;            % Pa 
    out.vcrit = 0.336;               % m^3/kmol 
    out.CritCompressFactor =0.254;   % -
    
%% 2 
    out.MeltingPoint = 181.25;             % K
    out.TripletPointT = 181.25;            % K
    out.TripletPointp = 0.043691;          % Pa
    out.NormalBoilingPoint = 379.25;       % K
    out.liqMolarVolume  = 0.115158;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -4.2710E8;             % J/kmol
    out.IG_GibbsOfFormation = -2.8440E08;           % J/kmol
    out.IG_AbsoluteEntropy = 4.0800E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 1.3900E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -2.7099E9;    % J/kmol
    out.StandardGibbsOfFormation = -2.9390E08;      % J/kmol
    out.StandardHeatOfFormation = -4.6920E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.9860E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.387793;          % -
    out.RadiusOfGyration = 3.7850E-10;      % m
    out.SolubilityParameter = 1.8230E4;     % (J/m^3)^.5 
    out.DipoleMoment = 6.7700E-30;          % C^m
    out.VanDerWaalsVolume = 0.0632;         % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 9.1400E08;        % m^2/kmol
    out.RefractiveIndex = 1.3874;           % -
    out.FlashPoint = 291.0;                 % K
    out.FlammabilityLimits = [1.7 8.0];     % vol-%
    out.FlammabilityLimitTemps = [284 313]; % K
    out.AutoignitionTemp = 595.0;           % K
    
%% 6 
    out.solDensity.A = 1.1400E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [181.25,1.1400E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 7.9185E-01;
    out.liqDensity.B = 2.6606E-01;
    out.liqDensity.C = 5.5900E02;
    out.liqDensity.D = 2.7850E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [181.25, 9.7549E00; 559.00, 2.9762E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.1554E02;
    out.vapPressure.B = -7.9070E03;
    out.vapPressure.C = -1.4945E01;
    out.vapPressure.D =  1.4711E-02;
    out.vapPressure.E =  1.0000E00;
    out.vapPressure.range = [181.25, 4.3690E-02; 559.00, 3.5210E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.5800E07;
    out.HeatOfVap.B = 4.1300E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [181.25, 4.7461E07; 559.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 0;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [NaN, NaN; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.0210E05;
    out.liqHeatCapacity.B = 3.2795E02;
    out.liqHeatCapacity.C = 0;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [290.15, 1.9725E05; 379.25, 2.2648E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0596E05;
    out.IG_HeatCapacity.B = 2.9480E05;
    out.IG_HeatCapacity.C = 1.8990E03;
    out.IG_HeatCapacity.D = 2.2542E05;
    out.IG_HeatCapacity.E = -8.2219E02;
    out.IG_HeatCapacity.range = [298.15, 1.3348E05; 1500.00, 3.3415E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.0700E-01;
    out.SecondViralCoef.B = -2.3200E02;
    out.SecondViralCoef.C = -4.0100E07;
    out.SecondViralCoef.D = -1.7314E20;
    out.SecondViralCoef.E = 2.9240E20;
    out.SecondViralCoef.range = [279.50, -4.2995E00; 1500.00, 4.0446E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.7801E01;
    out.liqViscosity.B = 1.4978E03;
    out.liqViscosity.C = 9.5330E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [250.00, 1.4356E-03; 379.25,2.7721E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 6.2067E-07;
    out.vapViscosity.B = 5.7800E-01;
    out.vapViscosity.C = 3.6965E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [181.25, 4.1243E-06; 1000.00, 2.4561E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1050E-01;
    out.liqThermalConductivity.B = -2.3200E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [181.25, 1.6845E-01; 503.15, 9.3769E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 1.3460E03;
    out.vapThermalConductivity.B = 9.9860E-01;
    out.vapThermalConductivity.C = 1.2465E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [379.25, 1.5403E-02; 999.25, 1.0678E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.7032E-02;
    out.SurfaceTension.B = 1.1073E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [181.25, 3.6953E-02; 559.00, 0];
    out.SurfaceTension.eqn = 106;    
    
end


