function out = pc_2methyloctane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2methyloctane','2-methyloctane'}; 
    out.CAS = '3221-61 -2';
    % components
    out.ChemicalFormular.C = 9;         
    out.ChemicalFormular.H = 20;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 128.258;  % kg/kmol    
    out.Tcrit = 582.80;             % K    
    out.pcrit = 2.3100E06;          % Pa 
    out.vcrit = 0.54100;            % m^3/kmol 
    out.CritCompressFactor = 0.258; % -
    
%% 2 
    out.MeltingPoint = 192.78;              % K
    out.TripletPointT = 192.78;             % K
    out.TripletPointp = 1.8181E-02;         % Pa
    out.NormalBoilingPoint = 416.43;        % K
    out.liqMolarVolume = 0.180734;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.3590E08;             % J/kmol
    out.IG_GibbsOfFormation = 2.0080E7;            % J/kmol
    out.IG_AbsoluteEntropy =  4.9880E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.8000E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -5.6791E09;    % J/kmol
    out.StandardGibbsOfFormation = 8.2190E06;              % J/kmol
    out.StandardHeatOfFormation = -2.8082E08;               % J/kmol
    out.StandardAbsoluteEntropy = 3.8790E05;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4583;            % -
    out.RadiusOfGyration = 4.9730E-10;      % m
    out.SolubilityParameter = 1.53540E4;    % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.09894;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.3680E09;        % m^2/kmol
    out.RefractiveIndex = 1.4008;           % -
    out.FlashPoint = 296.00;                % K
    out.FlammabilityLimits = [0.9 5.4];     % vol-%
    out.FlammabilityLimitTemps = [298 333]; % K
    out.AutoignitionTemp = 493.15;          % K
    
%% 6 
    out.solDensity.A = 7.1600E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [192.78, 7.1600E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 4.8369E-01;
    out.liqDensity.B = 2.6168E-01;
    out.liqDensity.C = 5.8280E02;
    out.liqDensity.D = 2.8064E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [192.78, 6.1230E00; 582.80, 1.8484E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.8095E01;
    out.vapPressure.B = -7.7507E03;
    out.vapPressure.C = -7.9631E00;
    out.vapPressure.D =  1.4327E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [192.78, 1.8180E-02; 582.80, 2.3026E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.9830E07;
    out.HeatOfVap.B = 3.9749E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [192.78, 5.1002E07; 582.80, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 3.4830E04;
    out.solHeatCapacity.B = 7.6538E02;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [80.0, 9.6060E04; 160.00, 1.5729E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.0720E05;
    out.liqHeatCapacity.B = 2.0160E01;
    out.liqHeatCapacity.C = 7.5925E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [192.78, 2.3930E05; 500.00, 4.0709E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.3743E05;
    out.IG_HeatCapacity.B = 5.1693E05;
    out.IG_HeatCapacity.C = 1.5604E03;
    out.IG_HeatCapacity.D = 3.4726E05;
    out.IG_HeatCapacity.E = 6.9593E02;
    out.IG_HeatCapacity.range = [150.00, 1.4022E05; 1500.00, 5.6295E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 3.0640E-01;
    out.SecondViralCoef.B = -3.8576E02;
    out.SecondViralCoef.C = -7.2988E07;
    out.SecondViralCoef.D = -1.6901E20;
    out.SecondViralCoef.E = 2.6260E22;
    out.SecondViralCoef.range = [291.40, -5.4846E00; 1500.00, 2.7595E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1283E01;
    out.liqViscosity.B = 1.1519E03;
    out.liqViscosity.C = 1.7447E-03;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [192.78, 4.9985E-03; 416.43, 2.0218E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 5.3451E-07;
    out.vapViscosity.B = 5.2741E-01;
    out.vapViscosity.C = 2.7431E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [192.78, 3.5382E-06; 1000.00, 1.6029E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0102E-01;
    out.liqThermalConductivity.B = -2.4765E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [192.78, 1.5328E-01; 469.00, 8.4872E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 4.6511E-05;
    out.vapThermalConductivity.B = 1.0987E00;
    out.vapThermalConductivity.C = -1.1797E02;
    out.vapThermalConductivity.D = 1.6726E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [416.43, 2.0894E-02; 1000.0, 8.7652E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.3154E-02;
    out.SurfaceTension.B = 1.2669E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [192.78, 3.1956E-02; 582.80, 0];
    out.SurfaceTension.eqn = 106;    
    
end


