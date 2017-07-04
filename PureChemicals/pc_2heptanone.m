function out = pc_2heptanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2heptanone','2-heptanone','2Heptanone','2-Heptanone'}; 
    out.CAS = '110-43-0';
    % components
    out.ChemicalFormular.C = 7;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 114.188;   % kg/kmol    
    out.Tcrit = 611.40;             % K    
    out.pcrit = 2.9400E06;           % Pa 
    out.vcrit = 0.43400;            % m^3/kmol 
    out.CritCompressFactor = 0.251; % -
    
%% 2 
    out.MeltingPoint = 238.15;              % K
    out.TripletPointT = 238.15;             % K
    out.TripletPointp = 3.5448E00;          % Pa
    out.NormalBoilingPoint = 424.18;        % K
    out.liqMolarVolume = 0.140675;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.0045E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.1960E08;            % J/kmol
    out.IG_AbsoluteEntropy =  2.8064E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.8400E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.0995E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.3480E08;              % J/kmol
    out.StandardHeatOfFormation = -3.4820E08;               % J/kmol
    out.StandardAbsoluteEntropy = 3.4100E05;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4190;            % -
    out.RadiusOfGyration = 4.4980E-10;      % m
    out.SolubilityParameter = 1.7800E04;     % (J/m^3)^.5 
    out.DipoleMoment = 8.7100E-30;          % C^m
    out.VanDerWaalsVolume = 0.07996;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.1240E09;        % m^2/kmol
    out.RefractiveIndex = 1.4066;           % -
    out.FlashPoint = 312.15;                % K
    out.FlammabilityLimits = [1.1 7.9];     % vol-%
    out.FlammabilityLimitTemps = [311 350]; % K
    out.AutoignitionTemp = 666.15;          % K
    
%% 6 
    out.solDensity.A = 8.8400E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [238.15, 8.8400E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.8247E-01;
    out.liqDensity.B = 2.5279E-01;
    out.liqDensity.C = 6.1140E02;
    out.liqDensity.D = 2.9818E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [238.15, 7.5514E00; 611.40, 2.3042E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.5494E01;
    out.vapPressure.B = -7.8965E03;
    out.vapPressure.C = -7.5047E00;
    out.vapPressure.D =  8.9130E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [238.15, 3.5448E00; 611.40, 2.9456E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.1425E07;
    out.HeatOfVap.B = 3.9802E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [238.15, 5.0471E07; 611.40, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 2.0390E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [238.15, 2.0390E05; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.6504E05;
    out.liqHeatCapacity.B = -3.7568E02;
    out.liqHeatCapacity.C = 1.0024E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [238.15, 2.3242E05; 490.00, 3.2163E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.2507E05;
    out.IG_HeatCapacity.B = 2.1480E05;
    out.IG_HeatCapacity.C = 6.9120E02;
    out.IG_HeatCapacity.D = 1.6190E05;
    out.IG_HeatCapacity.E = 1.7593E03;
    out.IG_HeatCapacity.range = [150.00, 1.2688E05; 1200.00, 3.8446E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.9183E-01;
    out.SecondViralCoef.B = -3.4109E02;
    out.SecondViralCoef.C = -7.8913E07;
    out.SecondViralCoef.D = -4.0026E20;
    out.SecondViralCoef.E = 5.6762E22;
    out.SecondViralCoef.range = [305.70, -6.3996E00; 1500.00, 4.1041E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.3929E01;
    out.liqViscosity.B = 1.3219E03;
    out.liqViscosity.C = 4.0382E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [250.00, 1.6422E-03; 424.18, 2.3184E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 8.8629E-08;
    out.vapViscosity.B = 7.8376E-01;
    out.vapViscosity.C = 1.0018E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [238.15, 4.5495E-06; 1000.00, 1.8088E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1080E-01;
    out.liqThermalConductivity.B = -2.4600E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [236.15, 1.5222E-01; 424.05, 1.0648E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.0493E03;
    out.vapThermalConductivity.B = 1.0323E00;
    out.vapThermalConductivity.C = 2.2983E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [424.18, 1.9506E-02; 1000.0, 1.1145E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.1492E-02;
    out.SurfaceTension.B = 1.2802E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [236.15, 3.2692E-02; 611.40, 0];
    out.SurfaceTension.eqn = 106;    
    
end


