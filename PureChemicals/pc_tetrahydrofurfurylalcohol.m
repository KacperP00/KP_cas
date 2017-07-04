function out = pc_tetrahydrofurfurylalcohol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'tetrahydrofurfurylalcohol','tetrahydrofuryl carbinol','oxolan-2-methanol'}; 
    out.CAS = '97-99-4';
    % components
    out.ChemicalFormular.C = 5;         
    out.ChemicalFormular.H = 10;
    out.ChemicalFormular.O = 2;
%% 1  
    out.MolecularWeight = 102.133;   % kg/kmol    
    out.Tcrit = 639.00;              % K    
    out.pcrit = 4.6600E6;            % Pa 
    out.vcrit = 0.29000;             % m^3/kmol 
    out.CritCompressFactor = 0.254;  % -
    
%% 2 
    out.MeltingPoint = NaN;             % K
    out.TripletPointT = NaN;            % K
    out.TripletPointp = NaN;            % Pa
    out.NormalBoilingPoint = 450.80;    % K
    out.liqMolarVolume  = 0.097319;     % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.6920E08;            % J/kmol
    out.IG_GibbsOfFormation = -2.1600E08;           % J/kmol
    out.IG_AbsoluteEntropy = 3.7400E5;              % J/kmolK
    out.HeatFusionAtMeltPoint = NaN;                % J/kmol
    out.StandardNetHeatOfCombustion = -2.7412E9;    % J/kmol
    out.StandardGibbsOfFormation = NaN;             % J/kmol
    out.StandardHeatOfFormation = NaN;              % J/kmol
    out.StandardAbsoluteEntropy = NaN;              % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.7044;            % -
    out.RadiusOfGyration = 3.4810E-10;      % m
    out.SolubilityParameter = 2.2920E4;     % (J/m^3)^.5 
    out.DipoleMoment = 7.3400E-30;          % C^m
    out.VanDerWaalsVolume = 0.05830;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 7.4600E08;        % m^2/kmol
    out.RefractiveIndex = 1.4499;           % -
    out.FlashPoint = 343.15;                % K
    out.FlammabilityLimits = [1.5 9.7];     % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp = 553.15;          % K
    
%% 6 
    out.solDensity.A = 1.2900E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [0.0, 1.2900E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 9.7001E-01;
    out.liqDensity.B = 2.8130E-01;
    out.liqDensity.C = 6.3900E02;
    out.liqDensity.D = 2.3837E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [193.00, 1.1046E01; 639.00, 3.4483E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  1.5059E02;
    out.vapPressure.B = -1.1574E04;
    out.vapPressure.C = -1.9025E01;
    out.vapPressure.D =  1.4141E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [193.00, 1.2655E-04; 639.00, 4.6469E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 6.4109E07;
    out.HeatOfVap.B = 2.8538E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [193.00, 5.7857E07; 639.00, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 0;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [NaN, NaN; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 5.2700E04;
    out.liqHeatCapacity.B = 4.3580E02;
    out.liqHeatCapacity.C = 0;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [296.65, 1.8198E05; 469.50, 2.5731E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 2.6823E05;
    out.IG_HeatCapacity.B = 1.1670E05;
    out.IG_HeatCapacity.C = 1.6860E03;
    out.IG_HeatCapacity.D = -3.5200E05;
    out.IG_HeatCapacity.E = 2.4520E02;
    out.IG_HeatCapacity.range = [300.00, 1.3998E05; 1500.00, 3.3693E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.9766E-01;
    out.SecondViralCoef.B = -1.8908E02;
    out.SecondViralCoef.C = -8.1632E07;
    out.SecondViralCoef.D = -5.5857E20;
    out.SecondViralCoef.E = 8.4850E22;
    out.SecondViralCoef.range = [319.50, -5.5954E00; 1500.00, 4.7400E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -7.9742E00;
    out.liqViscosity.B = 2.7454E03;
    out.liqViscosity.C = -1.1468E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [220.00, 1.8628E-01; 503.00, 6.4431E-05];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.6196E-07;
    out.vapViscosity.B = 7.4453E-01;
    out.vapViscosity.C = 1.4579E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [193.00, 4.6420E-06; 1000.00, 2.4204E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1782E-01;
    out.liqThermalConductivity.B = -2.4050E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [193.00, 1.7140E-01; 450.80, 1.0940E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 2.8570E-04;
    out.vapThermalConductivity.B = 9.0229E-01;
    out.vapThermalConductivity.C = 7.0294E02;
    out.vapThermalConductivity.D = 0;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [450.80, 2.7698E-02; 1000.00, 8.5424E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.5073E-02;
    out.SurfaceTension.B = -8.8258E-05;
    out.SurfaceTension.C = -1.6138E-08;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [294.75, 3.7657E-02; 450.80, 2.2007E-02];
    out.SurfaceTension.eqn = 100;    
    
end


