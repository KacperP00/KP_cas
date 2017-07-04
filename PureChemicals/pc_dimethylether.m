function out = pc_dimethylether

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'dimethylether','methyl ether','ether','ethyl ether'}; 
    out.CAS = '115-10-6';
    % components
    out.ChemicalFormular.C = 2;         
    out.ChemicalFormular.H = 6;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 46.069;   % kg/kmol    
    out.Tcrit = 400.10;             % K    
    out.pcrit = 5.3700E6;           % Pa 
    out.vcrit = 0.17000;            % m^3/kmol 
    out.CritCompressFactor = 0.274; % -
    
%% 2 
    out.MeltingPoint = 131.66;              % K
    out.TripletPointT = 131.65;             % K
    out.TripletPointp = 3.0496E00;          % Pa
    out.NormalBoilingPoint = 248.31;        % K
    out.liqMolarVolume = 0.070231;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.8410E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.1280E08;            % J/kmol
    out.IG_AbsoluteEntropy =  2.6670E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  4.9370E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -1.3284E09;    % J/kmol
    out.StandardGibbsOfFormation = NaN;              % J/kmol
    out.StandardHeatOfFormation = NaN;               % J/kmol
    out.StandardAbsoluteEntropy = NaN;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.2002;            % -
    out.RadiusOfGyration = 2.1540E-10;      % m
    out.SolubilityParameter = 1.5120E04;    % (J/m^3)^.5 
    out.DipoleMoment = 4.3400E-30;          % C^m
    out.VanDerWaalsVolume = 0.03104;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 4.8400E08;        % m^2/kmol
    out.RefractiveIndex = 1.2984;           % -
    out.FlashPoint = 232.00;                % K
    out.FlammabilityLimits = [3.3 27.3];    % vol-%
    out.FlammabilityLimitTemps = [NaN NaN]; % K
    out.AutoignitionTemp =623.15;           % K
    
%% 6 
    out.solDensity.A = 2.2200E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [131.65, 2.2200E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 1.5693E00;
    out.liqDensity.B = 2.6790E-01;
    out.liqDensity.C = 4.0010E02;
    out.liqDensity.D = 2.8820E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [131.65, 1.8950E01; 400.10, 5.8578E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A = 4.4704E01;
    out.vapPressure.B = -3.5256E03;
    out.vapPressure.C = -3.4444E00;
    out.vapPressure.D =  5.4574E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [131.65, 3.0496E00; 400.10, 5.2735E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 2.9940E07;
    out.HeatOfVap.B = 3.5050E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [131.65, 2.6032E07; 400.10, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.3880E04;
    out.solHeatCapacity.B = 1.2580E03;
    out.solHeatCapacity.C = -8.2200E00;
    out.solHeatCapacity.D = 2.6752E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [24.26, 1.2183E04; 131.64 , 7.0305E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.1010E05;
    out.liqHeatCapacity.B = -1.5747E02;
    out.liqHeatCapacity.C = 5.1853E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [131.65, 9.8356E04; 250.00, 1.0314E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 5.1480E04;
    out.IG_HeatCapacity.B = 1.4420E05;
    out.IG_HeatCapacity.C = 1.6034E03;
    out.IG_HeatCapacity.D = 7.7470E04;
    out.IG_HeatCapacity.E = 7.2540E02;
    out.IG_HeatCapacity.range = [200.00, 5.4363E04; 1500.00, 1.6581E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 1.1570E-01;
    out.SecondViralCoef.B = -1.1010E02;
    out.SecondViralCoef.C = -5.2550E06;
    out.SecondViralCoef.D = -4.6140E18;
    out.SecondViralCoef.E = 7.6350E20;
    out.SecondViralCoef.range = [200.00, -1.4028E00; 1500.00, 4.0743E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.0620E01;
    out.liqViscosity.B = 4.4899E02;
    out.liqViscosity.C = 8.3967E-05;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [131.65, 7.3982E-04; 248.31, 1.4903E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.6800E-06;
    out.vapViscosity.B = 3.9750E-01;
    out.vapViscosity.C = 5.3400E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [131.65, 3.6879E-06; 1000.00, 2.7215E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 3.1174E-01;
    out.liqThermalConductivity.B = -5.6380E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [131.65, 2.3752E-01; 320.03, 1.3131E-01];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 5.9975E-02;
    out.vapThermalConductivity.B = 2.6670E-01;
    out.vapThermalConductivity.C = 1.0186E03;
    out.vapThermalConductivity.D = 1.0988E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [248.31, 1.1388E-02; 1500.00, 1.9458E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 6.0960E-02;
    out.SurfaceTension.B = 1.2286E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [131.65, 3.7336E-02; 400.10, 0];
    out.SurfaceTension.eqn = 106;    
    
end


