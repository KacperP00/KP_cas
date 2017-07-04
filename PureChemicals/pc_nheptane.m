function out = pc_nheptane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'nheptane','n-heptane','dipropymethane'}; 
    out.CAS = '142-82-5';
    % components
    out.ChemicalFormular.C = 7;         
    out.ChemicalFormular.H = 16;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 100.204;   % kg/kmol    
    out.Tcrit = 540.20;              % K    
    out.pcrit = 2.7400E6;            % Pa 
    out.vcrit = 0.42800;             % m^3/kmol 
    out.CritCompressFactor = 0.261;  % -
    
%% 2 
    out.MeltingPoint = 182.57;             % K
    out.TripletPointT = 182.57;            % K
    out.TripletPointp = 1.8269E-01;        % Pa
    out.NormalBoilingPoint = 371.58;       % K
    out.liqMolarVolume  = 0.147024;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.8765E08;            % J/kmol
    out.IG_GibbsOfFormation = 8.1650E07;            % J/kmol
    out.IG_AbsoluteEntropy = 4.2798E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 1.4050E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.4647E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3495;            % -
    out.RadiusOfGyration = 4.1730E-10;      % m
    out.SolubilityParameter = 1.5200E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.07849;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.0990E09;        % m^2/kmol
    out.RefractiveIndex = 1.3851;           % -
    out.FlashPoint = 269.00;                % K
    out.FlammabilityLimits = [1.0 7.0];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 477.00;          % K
    
%% 6 
    out.solDensity.A = 9.1435E00;
    out.solDensity.B = -2.8756E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [133.15, 8.7606E00; 182.57, 8.6185E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.1259E-01;
    out.liqDensity.B = 2.6211E-01;
    out.liqDensity.C = 5.4020E02;
    out.liqDensity.D = 2.8141E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [182.57, 7.6998E00; 540.20, 2.3371E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  8.7829E01;
    out.vapPressure.B = -6.9964E03;
    out.vapPressure.C = -9.8802E00;
    out.vapPressure.D =  7.2099E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [182.57, 1.8269E-01; 540.20, 2.7192E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.0014E07;
    out.HeatOfVap.B = 3.8795E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [182.57, 4.2619E07; 540.20, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -2.3900E04;
    out.solHeatCapacity.B = 1.9408E03;
    out.solHeatCapacity.C = -1.0480E01;
    out.solHeatCapacity.D = 2.7190E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [20.00, 1.0942E04; 182.58, 1.4658E+05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 6.1260E01;
    out.liqHeatCapacity.B = 3.1441E05;
    out.liqHeatCapacity.C = 1.8246E03;
    out.liqHeatCapacity.D = -2.5479E03;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [182.57, 1.9989E05; 520.00, 4.0657E05];
    out.liqHeatCapacity.eqn = 114;    
    
    out.IG_HeatCapacity.A = 1.2015E05;
    out.IG_HeatCapacity.B = 4.0010E05;
    out.IG_HeatCapacity.C = 1.6766E03;
    out.IG_HeatCapacity.D = 2.7400E05;
    out.IG_HeatCapacity.E = 7.5640E02;
    out.IG_HeatCapacity.range = [200.00, 1.2828E05; 1500.00, 4.4283E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.7460E-01;
    out.SecondViralCoef.B = -2.9100E02;
    out.SecondViralCoef.C = -4.4180E07;
    out.SecondViralCoef.D = -8.8000E19;
    out.SecondViralCoef.E = 1.2850E22;
    out.SecondViralCoef.range = [270.10, - 4.4720E00; 1500.00, 6.7507E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.4451E01;
    out.liqViscosity.B = 1.5331E03;
    out.liqViscosity.C = 2.0087E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [182.57, 3.7197E-03; 373.15, 2.1455E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 6.6720E-08;
    out.vapViscosity.B = 8.2837E-01;
    out.vapViscosity.C = 8.5752E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [182.57, 3.3910E-06; 1000.00, 1.8777E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1500E-01;
    out.liqThermalConductivity.B = -3.0300E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [182.57, 1.5968E-01; 371.58, 1.0241E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -7.0028E-02;
    out.vapThermalConductivity.B = 3.8068E-01;
    out.vapThermalConductivity.C = -7.0499E03;
    out.vapThermalConductivity.D = -2.4005E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [339.15, 1.5827E-02; 1000.00, 1.1493E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.4143E-02;
    out.SurfaceTension.B = 1.2512E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [182.57, 3.2317E-02; 540.20, 0];
    out.SurfaceTension.eqn = 106;    
    
end


