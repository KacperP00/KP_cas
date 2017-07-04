function out = pc_1nonanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1nonanol','1-nonanol','1Nonanol','1-Nonanol','Alcohol C-9','n-Nonyl Alcohol'}; 
    out.CAS = '143-08-8';
    % components
    out.ChemicalFormular.C = 9;         
    out.ChemicalFormular.H = 20;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 144.257;   % kg/kmol    
    out.Tcrit = 670.70;             % K    
    out.pcrit = 2.5280E06;           % Pa 
    out.vcrit = 0.57200;            % m^3/kmol 
    out.CritCompressFactor = 0.259; % -
    
%% 2 
    out.MeltingPoint = 268.15;              % K
    out.TripletPointT = 268.15;             % K
    out.TripletPointp = 7.1702E-02;         % Pa
    out.NormalBoilingPoint = 486.25;        % K
    out.liqMolarVolume = 0.174942;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.8120E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.1220E08;            % J/kmol
    out.IG_AbsoluteEntropy =  5.5770E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  2.8800E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -5.5007E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.3300E08;       % J/kmol
    out.StandardHeatOfFormation = -4.5360E08;        % J/kmol
    out.StandardAbsoluteEntropy = 3.8460E05;         % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5997;            % -
    out.RadiusOfGyration = 5.1080E-10;      % m
    out.SolubilityParameter = 2.0720E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.3700E-30;          % C^m
    out.VanDerWaalsVolume = 0.10355;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.4380E09;        % m^2/kmol
    out.RefractiveIndex = 1.4319;           % -
    out.FlashPoint = 347.00;                % K
    out.FlammabilityLimits = [0.8 6.1];    % vol-%
    out.FlammabilityLimitTemps = [366 405]; % K
    out.AutoignitionTemp = 550.00;          % K
    
%% 6 
    out.solDensity.A = 6.8500E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [268.15, 6.8500E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 4.4231E-01;
    out.liqDensity.B = 2.5300E-01;
    out.liqDensity.C = 6.7070E+02;
    out.liqDensity.D = 2.5261E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [268.15, 5.8515E00; 670.70, 1.7483E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  2.5234E02;
    out.vapPressure.B = -1.8775E04;
    out.vapPressure.C = -3.3278E01;
    out.vapPressure.D =  1.5535E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [268.15, 7.1702E-02; 670.70, 2.5489E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 1.1484E08;
    out.HeatOfVap.B = 6.6695E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [268.15, 8.1701E07; 670.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 3.8373E05;
    out.solHeatCapacity.B = -1.4621E03;
    out.solHeatCapacity.C = 4.4071E00;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [200.00, 2.6759E05; 260.00, 3.0150E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.0483E07;
    out.liqHeatCapacity.B = -1.1522E05;
    out.liqHeatCapacity.C = 4.7687E02;
    out.liqHeatCapacity.D = -8.5381E-01;
    out.liqHeatCapacity.E = 5.6246E-04;
    out.liqHeatCapacity.range = [310.00, 3.5059E05; 460.00, 4.6494E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.5400E05;
    out.IG_HeatCapacity.B = 4.9360E05;
    out.IG_HeatCapacity.C = 1.580E03;
    out.IG_HeatCapacity.D = 3.5880E05;
    out.IG_HeatCapacity.E = 7.2111E02;
    out.IG_HeatCapacity.range = [200.00, 1.6777E05; 1500.00, 5.6606E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 4.3008E-01;
    out.SecondViralCoef.B = -5.1453E02;
    out.SecondViralCoef.C = -9.3063E07;
    out.SecondViralCoef.D = 2.9028E+21;
    out.SecondViralCoef.E = - 1.1382E24;
    out.SecondViralCoef.range = [335.35, -6.6432E00; 1500.00, 5.9569E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -7.1348E00;
    out.liqViscosity.B = 2.7763E03;
    out.liqViscosity.C = -1.2064E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [280.00, 1.8002E-02; 486.25, 1.3792E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.2000E-07;
    out.vapViscosity.B = 7.4000E-01;
    out.vapViscosity.C = 1.8000E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [268.15, 4.4992E-06; 1000.00, 1.6877E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2920E-01;
    out.liqThermalConductivity.B = -2.3000E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [268.15, 1.6753E-01; 578.65, 9.6111E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -3.0715E01;
    out.vapThermalConductivity.B = -1.0750E-01;
    out.vapThermalConductivity.C = 8.1070E03;
    out.vapThermalConductivity.D = -1.5683E08;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [486.25, 2.4464E-02; 1000.0, 9.8947E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.1752E-02;
    out.SurfaceTension.B = 1.0562E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [268.15, 3.0183E-02; 670.70, 0.0000E00];
    out.SurfaceTension.eqn = 106;    
    
end


