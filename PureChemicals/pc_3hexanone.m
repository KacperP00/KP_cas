function out = pc_3hexanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'3hexanone','3-hexanone','ethyl-n-propylketone'}; 
    out.CAS = '589-38-8';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 100.161;  % kg/kmol    
    out.Tcrit = 582.82;             % K    
    out.pcrit = 3.3200E6;           % Pa 
    out.vcrit = 0.37800;            % m^3/kmol 
    out.CritCompressFactor = 0.259; % -
    
%% 2 
    out.MeltingPoint = 217.15;              % K
    out.TripletPointT = 217.15;             % K
    out.TripletPointp = 2.2202E00;          % Pa
    out.NormalBoilingPoint = 396.65;        % K
    out.liqMolarVolume = 0.123630;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.7760E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.2600E8;            % J/kmol
    out.IG_AbsoluteEntropy =  4.0920E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.3490E7;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.4920E9;    % J/kmol
    out.StandardGibbsOfFormation = -1.3670E08;      % J/kmol
    out.StandardHeatOfFormation = -3.2020E08;       % J/kmol
    out.StandardAbsoluteEntropy = 3.0531E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.3801;            % -
    out.RadiusOfGyration = 4.0090E-10;      % m
    out.SolubilityParameter = 1.7930E4;     % (J/m^3)^.5 
    out.DipoleMoment = 9.5700E-30;          % C^m
    out.VanDerWaalsVolume = 0.06973;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 9.8900E08;        % m^2/kmol
    out.RefractiveIndex = 1.3980;           % -
    out.FlashPoint = 287.15;                % K
    out.FlammabilityLimits = [1.0 8.0];     % vol-%
    out.FlammabilityLimitTemps = [288 327]; % K
    out.AutoignitionTemp = 712.00;          % K
    
%% 6 
    out.solDensity.A = 1.0200E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [217.50, 1.0200E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.7666E-01;
    out.liqDensity.B = 2.5578E-01;
    out.liqDensity.C = 5.8282E02;
    out.liqDensity.D = 2.7746E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [217.50, 8.7631E00; 582.82, 2.6455E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.3155E01;
    out.vapPressure.B = -7.2429E03;
    out.vapPressure.C = -7.2569E00;
    out.vapPressure.D =  1.2741E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [217.50, 2.2202E00; 582.82, 3.3216E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.6232E07;
    out.HeatOfVap.B = 3.9972E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [217.50, 4.6655E07; 582.82, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -4.1000E03;
    out.solHeatCapacity.B = 4.2800E02;
    out.solHeatCapacity.C = 2.6860E01;
    out.solHeatCapacity.D = -3.5340E-01;
    out.solHeatCapacity.E = 1.3463E-03;
    out.solHeatCapacity.range = [12.52, 4.8084E03; 140.13, 1.2999E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.3596E05;
    out.liqHeatCapacity.B = -3.4594E02;
    out.liqHeatCapacity.C = 9.4278E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [217.50, 2.0532E05; 460.00, 2.7632E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.1237E05;
    out.IG_HeatCapacity.B = 2.9360E05;
    out.IG_HeatCapacity.C = 1.4010E03;
    out.IG_HeatCapacity.D = 1.6010E05;
    out.IG_HeatCapacity.E = 6.5050E02;
    out.IG_HeatCapacity.range = [150.00, 1.1443E05; 1500.00, 3.5874E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.4188E-01;
    out.SecondViralCoef.B = -2.7496E02;
    out.SecondViralCoef.C = -5.8174E07;
    out.SecondViralCoef.D = -2.1114E20;
    out.SecondViralCoef.E = 2.4535E22;
    out.SecondViralCoef.range = [297.00, -5.0275E00; 1500.00, 4.1329E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.3682E01;
    out.liqViscosity.B = 1.2834E03;
    out.liqViscosity.C = 3.3755E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [217.50, 2.5631E-03; 396.65, 2.1848E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 9.8882E-08;
    out.vapViscosity.B = 7.7550E-01;
    out.vapViscosity.C = 9.9825E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [217.50, 4.4032E-06; 1000.00, 1.9067E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.3493E-01;
    out.liqThermalConductivity.B = -2.9120E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [217.50, 1.7159E-01; 466.00, 9.9231E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -3.3262E-01;
    out.vapThermalConductivity.B = 1.2054E-01;
    out.vapThermalConductivity.C = -2.4726E03;
    out.vapThermalConductivity.D = -5.4934E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.00, 7.9990E-03; 1000.0, 1.0980E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.7740E-02;
    out.SurfaceTension.B = 1.1689E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [217.50, 3.3447E-02; 582.82, 0];
    out.SurfaceTension.eqn = 106;    
    
end


