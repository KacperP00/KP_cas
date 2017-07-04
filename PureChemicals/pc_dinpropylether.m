function out = pc_dinpropylether

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'dinpropylether','di-n-propyl-ether','n-propyl ether'}; 
    out.CAS = '1 11-43-3';
    % components
    out.ChemicalFormular.C = 8;         
    out.ChemicalFormular.H = 18;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 102.177;  % kg/kmol    
    out.Tcrit = 530.6;              % K    
    out.pcrit = 3.0280E6;           % Pa 
    out.vcrit = 0.382;              % m^3/kmol 
    out.CritCompressFactor = 0.262; % -
    
%% 2 
    out.MeltingPoint = 149.95;              % K
    out.TripletPointT = 149.95;             % K
    out.TripletPointp = 7.6265E-04;         % Pa
    out.NormalBoilingPoint = 363.23;        % K
    out.liqMolarVolume = 0.137649;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -2.9290E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.0550E08;            % J/kmol
    out.IG_AbsoluteEntropy =  4.2250E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  1.0770E07;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.7250E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.1204E07;       % J/kmol
    out.StandardHeatOfFormation = -3.2882E08;        % J/kmol
    out.StandardAbsoluteEntropy = 3.2388E05;         % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.368756;          % -
    out.RadiusOfGyration = 3.9140E-10;      % m
    out.SolubilityParameter = 1.5550E04;    % (J/m^3)^.5 
    out.DipoleMoment = 4.0400E-30;          % C^m
    out.VanDerWaalsVolume = 0.07196;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.0220E09;        % m^2/kmol
    out.RefractiveIndex = 1.378;            % -
    out.FlashPoint = 294.15;                % K
    out.FlammabilityLimits = [1.2 9.5];     % vol-%
    out.FlammabilityLimitTemps = [264 301]; % K
    out.AutoignitionTemp =488.15;           % K
    
%% 6 
    out.solDensity.A = 9.8800E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [149.95, 9.8800E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 7.2390E-01;
    out.liqDensity.B = 2.7653E-01;
    out.liqDensity.C = 5.3060E02;
    out.liqDensity.D = 2.7940E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [149.95, 8.4474E00; 530.60, 2.6178E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A = 6.5935E01;
    out.vapPressure.B = -6.1985E03;
    out.vapPressure.C = -6.3423E00;
    out.vapPressure.D =  2.0909E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [149.95, 7.6265E-04; 530.60, 3.0416E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.8430E07;
    out.HeatOfVap.B = 3.6670E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [149.95, 4.2877E07; 530.6, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.8140E04;
    out.solHeatCapacity.B = 1.6964E03;
    out.solHeatCapacity.C = -8.0340E00;
    out.solHeatCapacity.D = 1.8100E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [30.00, 2.6010E04; 149.95, 1.1662E05];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.9374E05;
    out.liqHeatCapacity.B = -1.2470E02;
    out.liqHeatCapacity.C = 7.3500E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [149.95, 1.9157E05; 363.23, 2.4542E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.2030E05;
    out.IG_HeatCapacity.B = 3.4123E05;
    out.IG_HeatCapacity.C = 1.5455E03;
    out.IG_HeatCapacity.D = 2.0480E05;
    out.IG_HeatCapacity.E = 7.2320E02;
    out.IG_HeatCapacity.range = [200.00, 1.2805E05; 1500.00, 4.0072E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.4620E-01;
    out.SecondViralCoef.B = -2.6016E02;
    out.SecondViralCoef.C = -3.9596E07;
    out.SecondViralCoef.D = 1.0688E19;
    out.SecondViralCoef.E = -1.1007E22;
    out.SecondViralCoef.range = [265.50, -4.0957E00; 1500.00, 6.1028E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -8.2119E00;
    out.liqViscosity.B = 8.9352E02;
    out.liqViscosity.C = -4.5940E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [149.95, 1.0516E-02; 363.23, 2.1174E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 4.0850E-08;
    out.vapViscosity.B = 8.9530E-01;
    out.vapViscosity.C = 0;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [149.95,3.6251E-06; 1000.00, 1.9819E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.1370E-01;
    out.liqThermalConductivity.B = -2.9200E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [149.95, 1.6991E-01; 420.00, 9.1060E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 7.9334E-02;
    out.vapThermalConductivity.B = 2.8330E-01;
    out.vapThermalConductivity.C = 4.5613E03;
    out.vapThermalConductivity.D = 1.3546E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [363.23, 1.7690E-02; 1000.00, 8.1192E-02];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.4950E-02;
    out.SurfaceTension.B = 1.2254E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [149.95, 3.6578E-02; 530.60, 0];
    out.SurfaceTension.eqn = 106;    
    
end


