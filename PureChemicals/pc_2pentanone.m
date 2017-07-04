function out = pc_2pentanone
%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2pentanone','2-pentanone','Ethylacetone'}; 
    out.CAS = '107-87-9';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 86.1338;  % kg/kmol    
    out.Tcrit = 561.08;             % K    
    out.pcrit = 3.6940E06;          % Pa 
    out.vcrit = 0.301;              % m^3/kmol 
    out.CritCompressFactor = 0.238; % -
    
%% 2 
    out.MeltingPoint = 196.29;              % K
    out.TripletPointT = 196.29;             % K
    out.TripletPointp = 0.752347;           % Pa
    out.NormalBoilingPoint = 375.46;        % K
    out.liqMolarVolume = 0.10725;           % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.5920E8;             % J/kmol
    out.IG_GibbsOfFormation = -1.3830E8;            % J/kmol
    out.IG_AbsoluteEntropy =  3.7860E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.0620E07;         % J/kmol
    out.StandardNetHeatOfCombustion = -2.8796E09;   % J/kmol
    out.StandardGibbsOfFormation = -1.4520E08;      % J/kmol
    out.StandardHeatOfFormation = -2.9729E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.7410E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.343288;          % -
    out.RadiusOfGyration = 3.6180E-10;      % m
    out.SolubilityParameter = 1.8290E04;    % (J/m^3)^.5 
    out.DipoleMoment = 9.2400E-30;          % C^m
    out.VanDerWaalsVolume = 0.0595;         % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.5400E08;        % m^2/kmol
    out.RefractiveIndex = 1.388;            % -
    out.FlashPoint = 280.0;                 % K
    out.FlammabilityLimits = [1.5 8.2];    % vol-%
    out.FlammabilityLimitTemps = [278 309]; % K
    out.AutoignitionTemp = 725.15;          % K
    
%% 6 
    out.solDensity.A = 1.2200E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [196.29, 1.2200E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.0411E-01;
    out.liqDensity.B = 2.7207E-01;
    out.liqDensity.C = 5.6108E02;
    out.liqDensity.D = 3.0669E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [196.29, 1.0398E01; 561.08, 3.3231E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  8.4635E01;
    out.vapPressure.B = -7.0784E03;
    out.vapPressure.C = -9.3000E00;
    out.vapPressure.D =  6.2702E-06;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [196.29, 7.5235E-01; 561.08, 3.7062E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 5.1740E07;
    out.HeatOfVap.B = 3.9422E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [196.29, 4.3663E07; 561.08, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -2.5097E04;
    out.solHeatCapacity.B = 2.0373E03;
    out.solHeatCapacity.C = -1.4289E01;
    out.solHeatCapacity.D = 4.4874E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [30.00, 2.4373E04; 110.00, 8.5836E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.9459E05;
    out.liqHeatCapacity.B = -2.6386E02;
    out.liqHeatCapacity.C = 7.6808E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [196.29, 1.7239E05; 375.46, 2.0380E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 9.0053E04;
    out.IG_HeatCapacity.B = 2.7085E05;
    out.IG_HeatCapacity.C = 1.6592E03;
    out.IG_HeatCapacity.D = 1.8012E05;
    out.IG_HeatCapacity.E = 7.4396E02;
    out.IG_HeatCapacity.range = [200.00, 9.5908E04; 1500.00, 3.0797E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.2325E-01;
    out.SecondViralCoef.B = -2.5789E02;
    out.SecondViralCoef.C = -3.5849E07;
    out.SecondViralCoef.D = -3.0892E20;
    out.SecondViralCoef.E = 5.9158E22;
    out.SecondViralCoef.range = [280.54, -4.8752E00; 1500.00, 4.0691E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.1055E01;
    out.liqViscosity.B = 1.0053E03;
    out.liqViscosity.C = 3.9301E-03;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [250.00, 9.0091E-04; 375.46, 2.3541E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.4630E-07;
    out.vapViscosity.B = 6.6530E-01;
    out.vapViscosity.C = 2.0870E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [196.29, 4.0030E-06; 1000.00, 2.0186E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1610E-01;
    out.liqThermalConductivity.B = -2.4866E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [196.29, 1.6729E-01; 375.46, 1.2274E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -1.7190E-02;
    out.vapThermalConductivity.B = 4.8320E-01;
    out.vapThermalConductivity.C = -3.7980E03;
    out.vapThermalConductivity.D = -1.2350E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [273.00, 8.7672E-03; 1000.0, 1.2002E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.3802E-02;
    out.SurfaceTension.B = 1.0733E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [196.29, 3.3893E-02; 561.08, 0];
    out.SurfaceTension.eqn = 106;    
    
end


