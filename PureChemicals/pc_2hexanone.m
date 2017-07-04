function out = pc_2hexanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2hexanone','2-hexanone','2-Hexanone','2Hexanone'}; 
    out.CAS = '591-78-6';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 100.161;   % kg/kmol    
    out.Tcrit = 587.61;             % K    
    out.pcrit = 3.2870E06;          % Pa 
    out.vcrit = 0.37800;            % m^3/kmol 
    out.CritCompressFactor = 0.254; % -
    
%% 2 
    out.MeltingPoint = 217.35;              % K
    out.TripletPointT = 217.35;             % K
    out.TripletPointp = 1.4456E00;          % Pa
    out.NormalBoilingPoint = 400.70;        % K
    out.liqMolarVolume = 0.124128;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.7983E08;            % J/kmol
    out.IG_GibbsOfFormation = -1.3008E08;           % J/kmol
    out.IG_AbsoluteEntropy = 4.1786E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = 1.4900E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.4900E09;   % J/kmol
    out.StandardGibbsOfFormation = -1.3960E08;      % J/kmol
    out.StandardHeatOfFormation = -3.2201E08;       % J/kmol
    out.StandardAbsoluteEntropy = 3.0811E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3846;            % -
    out.RadiusOfGyration = 4.0930E-10;      % m
    out.SolubilityParameter = 1.8050E04;    % (J/m^3)^.5 
    out.DipoleMoment = 8.9400E-30;          % C^m
    out.VanDerWaalsVolume = 0.06973;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 9.8900E08;        % m^2/kmol
    out.RefractiveIndex = 1.3987;           % -
    out.FlashPoint = 296.15;                % K
    out.FlammabilityLimits = [1.2 8.0];     % vol-%
    out.FlammabilityLimitTemps = [294 330]; % K
    out.AutoignitionTemp = 696.15;          % K
    
%% 6 
    out.solDensity.A = 1.0200E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [217.35, 1.0200E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.7816E-01;
    out.liqDensity.B = 2.5634E-01;
    out.liqDensity.C = 5.8761E02;
    out.liqDensity.D = 2.8365E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [217.35, 8.7319E00; 587.61, 2.6455E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.0744E02;
    out.vapPressure.B = -8.5286E03;
    out.vapPressure.C = -1.2679E01;
    out.vapPressure.D =  8.4606E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [217.35, 1.4456E00; 587.61, 3.2859E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.6231E07;
    out.HeatOfVap.B = 3.8207E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [217.35, 4.7135E07; 587.61, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 5.4305E03;
    out.solHeatCapacity.B = 8.2129E02;
    out.solHeatCapacity.C = 5.0967E00;
    out.solHeatCapacity.D = -7.0876E-02;
    out.solHeatCapacity.E = 2.1969E-04;
    out.solHeatCapacity.range = [39.76, 4.2236E04; 217.69, 1.8794E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.0825E05;
    out.liqHeatCapacity.B = -1.0747E02;
    out.liqHeatCapacity.C = 2.0620E-01;
    out.liqHeatCapacity.D = 7.0293E-04;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [217.35, 2.0185E05; 460.00, 2.7087E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0940E05;
    out.IG_HeatCapacity.B = 1.8070E05;
    out.IG_HeatCapacity.C = 6.8900E02;
    out.IG_HeatCapacity.D = 1.4740E05;
    out.IG_HeatCapacity.E = 1.7720E02;
    out.IG_HeatCapacity.range = [200.00, 1.1815E05; 1200.00, 3.3207E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.4931E-01;
    out.SecondViralCoef.B = -2.8704E02;
    out.SecondViralCoef.C = -5.8288E07;
    out.SecondViralCoef.D = -2.8661E20;
    out.SecondViralCoef.E = 4.2633E22;
    out.SecondViralCoef.range = [293.80, -5.5749E00; 1500.00, 4.0669E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1445E01;
    out.liqViscosity.B = 1.1872E03;
    out.liqViscosity.C = 2.9076E-03;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [217.35, 2.5614E-03; 400.70, 2.1077E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 9.7820E-08;
    out.vapViscosity.B = 7.7720E-01;
    out.vapViscosity.C = 9.9530E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [217.35, 4.3968E-06; 1000.00, 1.9091E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1076E-01;
    out.liqThermalConductivity.B = -2.4000E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [217.35, 1.5860E-01; 400.85, 1.1456E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -1.2158E00;
    out.vapThermalConductivity.B = 2.6637E-02;
    out.vapThermalConductivity.C = -1.7116E03;
    out.vapThermalConductivity.D = -1.3176E07;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.00, 7.7543E-03; 1000.0, 1.0523E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.2446E-02;
    out.SurfaceTension.B = 1.2674E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [217.35, 3.4777E-02; 587.61, 0];
    out.SurfaceTension.eqn = 106;    
    
end


