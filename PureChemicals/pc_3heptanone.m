function out = pc_3heptanone

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'3heptanone','3-heptanone','butyl ethyl ketone','heptane-3-one'}; 
    out.CAS = '106-35-4';
    % components
    out.ChemicalFormular.C = 7;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 114.188;  % kg/kmol    
    out.Tcrit = 606.60;             % K    
    out.pcrit = 2.9200E6;           % Pa 
    out.vcrit = 0.43300;            % m^3/kmol 
    out.CritCompressFactor = 0.251; % -
    
%% 2 
    out.MeltingPoint = 234.15;              % K
    out.TripletPointT = 234.15;             % K
    out.TripletPointp = 2.2992E00;          % Pa
    out.NormalBoilingPoint = 420.55;        % K
    out.liqMolarVolume = 0.140257;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.0100E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.2250E8;            % J/kmol
    out.IG_AbsoluteEntropy =  4.5800E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.8100E7;          % J/kmol
    out.StandardNetHeatOfCombustion = -4.0980E9;    % J/kmol
    out.StandardGibbsOfFormation = -1.3580E08;      % J/kmol
    out.StandardHeatOfFormation = -3.4860E08;       % J/kmol
    out.StandardAbsoluteEntropy = 3.4140E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.4076;            % -
    out.RadiusOfGyration = 4.5390E-10;      % m
    out.SolubilityParameter = 1.7940E4;     % (J/m^3)^.5 
    out.DipoleMoment = 9.3700E-30;          % C^m
    out.VanDerWaalsVolume = 0.07996;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.1240E09;        % m^2/kmol
    out.RefractiveIndex = 1.4066;           % -
    out.FlashPoint = 310.15;                % K
    out.FlammabilityLimits = [1.1 6.8];     % vol-%
    out.FlammabilityLimitTemps = [309 344]; % K
    out.AutoignitionTemp = 683.00;          % K
    
%% 6 
    out.solDensity.A = 8.8600E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [234.15, 8.8600E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 5.9268E-01;
    out.liqDensity.B = 2.5663E-01;
    out.liqDensity.C = 6.0660E02;
    out.liqDensity.D = 2.7766E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [234.15, 7.5751E00; 606.60, 2.3095E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  7.8463E01;
    out.vapPressure.B = -8.0772E03;
    out.vapPressure.C = -7.9062E00;
    out.vapPressure.D =  8.0521E-18;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [234.15, 2.2992E-04; 606.60, 2.9190E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.3357E07;
    out.HeatOfVap.B = 4.2167E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [234.15, 5.1579E07; 606.60, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 2.0400E05;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [234.15, 2.0400E05; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 2.7073E05;
    out.liqHeatCapacity.B = -3.99893E02;
    out.liqHeatCapacity.C = 1.0601E00;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [234.15, 2.3522E05; 480.00, 3.2303E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.2768E05;
    out.IG_HeatCapacity.B = 3.3810E05;
    out.IG_HeatCapacity.C = 1.3831E03;
    out.IG_HeatCapacity.D = 1.8880E05;
    out.IG_HeatCapacity.E = 6.5030E02;
    out.IG_HeatCapacity.range = [200.00, 1.3968E05; 1500.00, 4.1386E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.9056E-01;
    out.SecondViralCoef.B = -3.3949E02;
    out.SecondViralCoef.C = -7.6495E07;
    out.SecondViralCoef.D = -4.1164E20;
    out.SecondViralCoef.E = 6.0973E22;
    out.SecondViralCoef.range = [303.3, -6.5114E00; 1500.00, 4.1554E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -9.3877E00;
    out.liqViscosity.B = 1.2049E03;
    out.liqViscosity.C = -3.2618E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [234.15, 2.4268E-03; 421.15, 2.0397E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 8.9656E-08;
    out.vapViscosity.B = 7.8236E-01;
    out.vapViscosity.C = 1.0014E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [234.15, 4.4848E-06; 1000.00, 1.8122E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.0260E-01;
    out.liqThermalConductivity.B = -2.2340E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [234.15, 1.5029E-01; 553.15, 7.9026E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 1.3486E03;
    out.vapThermalConductivity.B = 1.0313E00;
    out.vapThermalConductivity.C = 1.4832E10;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [420.55, 1.9429E-02; 1000.0, 1.1287E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.8308E-02;
    out.SurfaceTension.B = 1.2123E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [234.15, 3.2279E-02; 606.60, 8.3816E-04];
    out.SurfaceTension.eqn = 106;    
    
end


