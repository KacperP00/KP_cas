function out = pc_npentane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'npentane','n-pentane','pental hydride'}; 
    out.CAS = '109-66-0';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 72.150;    % kg/kmol    
    out.Tcrit = 469.10;              % K    
    out.pcrit = 3.3700E6;            % Pa 
    out.vcrit = 0.31300;             % m^3/kmol 
    out.CritCompressFactor = 0.270;  % -
    
%% 2 
    out.MeltingPoint = 143.42;             % K
    out.TripletPointT = 143.42;            % K
    out.TripletPointp = 6.8642E-02;        % Pa
    out.NormalBoilingPoint = 309.22;       % K
    out.liqMolarVolume  = 0.116045;        % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.4676E08;            % J/kmol
    out.IG_GibbsOfFormation = 8.8130E06;            % J/kmol
    out.IG_AbsoluteEntropy = 3.4945E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = 8.4015E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.2449E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2515;            % -
    out.RadiusOfGyration = 3.3370E-10;      % m
    out.SolubilityParameter = 1.4400E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.05803;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.2900E08;        % m^2/kmol
    out.RefractiveIndex = 1.3547;           % -
    out.FlashPoint = 233.15;                % K
    out.FlammabilityLimits = [1.3 8.0];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 516.00;          % K
    
%% 6 
    out.solDensity.A = 1.3533E01;
    out.solDensity.B = -6.0608E-03;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [23.15, 1.3393E01; 143.42, 1.2664E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.4947E-01;
    out.liqDensity.B = 2.6726E-01;
    out.liqDensity.C = 4.6970E02;
    out.liqDensity.D = 2.7789E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [143.42, 1.0474E01; 469.70, 3.1784E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.8741E01;
    out.vapPressure.B = -5.4203E03;
    out.vapPressure.C = -8.8253E00;
    out.vapPressure.D =  9.6171E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [143.42, 6.8642E-02; 469.70, 3.3642E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 3.9109E07;
    out.HeatOfVap.B = 3.8681E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [143.42, 3.3968E07; 469.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.0050E04;
    out.solHeatCapacity.B = 1.0103E03;
    out.solHeatCapacity.C = 4.2700E00;
    out.solHeatCapacity.D = -1.0169E-01;
    out.solHeatCapacity.E = 4.2200E-04;
    out.solHeatCapacity.range = [12.29, 2.8324E03; 134.60, 9.3832E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.5908E05;
    out.liqHeatCapacity.B = -2.7050E02;
    out.liqHeatCapacity.C = 9.9537E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [143.42, 1.4076E05; 390.00, 2.0498E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 8.8050E04;
    out.IG_HeatCapacity.B = 3.0110E05;
    out.IG_HeatCapacity.C = 1.6502E03;
    out.IG_HeatCapacity.D = 1.8920E05;
    out.IG_HeatCapacity.E = 7.4760E02;
    out.IG_HeatCapacity.range = [200.00, 9.4039E04; 1500.00, 3.2927E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.9160E-01;
    out.SecondViralCoef.B = -1.9800E02;
    out.SecondViralCoef.C = -1.6400E07;
    out.SecondViralCoef.D = -1.7920E19;
    out.SecondViralCoef.E = 2.6850E21;
    out.SecondViralCoef.range = [234.80,-2.6205E00; 1500.00, 5.4740E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -2.0383E01;
    out.liqViscosity.B = 1.0504E03;
    out.liqViscosity.C = 1.4874E00;
    out.liqViscosity.D = -2.0167E-27;
    out.liqViscosity.E = 1E01;
    out.liqViscosity.range = [143.42, 3.4377E-03; 465.15, 4.7971E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 6.3412E-08;
    out.vapViscosity.B = 8.4758E-01;
    out.vapViscosity.C = 4.1718E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [143.42, 3.3051E-06; 1000.00, 2.1240E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.5370E-01;
    out.liqThermalConductivity.B = -5.7600E-04;
    out.liqThermalConductivity.C = 3.4400E-07;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [143.42, 1.7817E-01; 445.00, 6.5501E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -6.8440E02;
    out.vapThermalConductivity.B = 7.6400E-01;
    out.vapThermalConductivity.C = -1.0550E09;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.15, 1.2879E-02; 1000.00, 1.2707E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.2020E-02;
    out.SurfaceTension.B = 1.2041E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [143.42, 3.3546E-02; 469.70, 0];
    out.SurfaceTension.eqn = 106;    
    
end


