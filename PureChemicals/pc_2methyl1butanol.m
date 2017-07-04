function out = pc_2methyl1butanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'2methyl1butanol','2methyl-1butanol','amylalcohol'}; 
    out.CAS = '137-32-6';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 12;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 88.150;   % kg/kmol    
    out.Tcrit = 575.40;             % K    
    out.pcrit = 3.9400E06;           % Pa 
    out.vcrit = 0.32700;            % m^3/kmol 
    out.CritCompressFactor = 0.269; % -
    
%% 2 
    out.MeltingPoint = 195.00;              % K
    out.TripletPointT = 195.00;             % K
    out.TripletPointp = 1.2845E-03;         % Pa
    out.NormalBoilingPoint = 401.85;        % K
    out.liqMolarVolume = 0.108257;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -5.0209E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.4671E08;            % J/kmol
    out.IG_AbsoluteEntropy =  3.9351E05;             % J/kmolK
    out.HeatFusionAtMeltPoint =  8.2700E06;          % J/kmol
    out.StandardNetHeatOfCombustion = -3.0620E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.5940E08;              % J/kmol
    out.StandardHeatOfFormation = -3.5660E08;               % J/kmol
    out.StandardAbsoluteEntropy = 2.5330E05;               % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.5736;            % -
    out.RadiusOfGyration = 3.6120E-10;      % m
    out.SolubilityParameter = 2.2090E+04;     % (J/m^3)^.5 
    out.DipoleMoment = 6.2700E-30;          % C^m
    out.VanDerWaalsVolume = 0.06262;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 8.9700E08;        % m^2/kmol
    out.RefractiveIndex = 1.4086;           % -
    out.FlashPoint = 316.15;                % K
    out.FlammabilityLimits = [1.4 9.0];    % vol-%
    out.FlammabilityLimitTemps = [313 345]; % K
    out.AutoignitionTemp = 658.15;          % K
    
%% 6 
    out.solDensity.A = 1.1800E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [195.0, 1.1800E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 8.1662E-01;
    out.liqDensity.B = 2.6703E-01;
    out.liqDensity.C = 5.7540E02;
    out.liqDensity.D = 2.4336E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [195.00, 1.0092E01; 575.40, 3.0582E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  9.6196E01;
    out.vapPressure.B = -9.6589E03;
    out.vapPressure.C = -1.0112E01;
    out.vapPressure.D =  8.4773E-19;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [195.00, 1.2846E-03; 575.40, 3.9057E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 7.7920E07;
    out.HeatOfVap.B = 4.6991E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [195.00, 6.4149E07; 575.40, 0];
    out.HeatOfVap.eqn = 106;     

    %keine Daten in Stoffdatenblatt
    out.solHeatCapacity.A = 0;
    out.solHeatCapacity.B = 0;
    out.solHeatCapacity.C = 0;
    out.solHeatCapacity.D = 0;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [NaN,NaN; NaN, NaN];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 8.3025E04;
    out.liqHeatCapacity.B = 4.5971E02;
    out.liqHeatCapacity.C = 0;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [195.00, 1.7267E05; 401.85, 2.6776E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.0890E05;
    out.IG_HeatCapacity.B = 2.1850E05;
    out.IG_HeatCapacity.C = 8.5300E02;
    out.IG_HeatCapacity.D = 1.4000E05;
    out.IG_HeatCapacity.E = 2.9060E03;
    out.IG_HeatCapacity.range = [298.15, 1.3247E05; 1500.00, 3.4718E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.1048E-01;
    out.SecondViralCoef.B = -2.0521E02;
    out.SecondViralCoef.C = -5.4378E07;
    out.SecondViralCoef.D = -2.0676E20;
    out.SecondViralCoef.E = 2.8200E22;
    out.SecondViralCoef.range = [287.70, -5.1030E00; 1500.00, 5.7554E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -6.4568E01;
    out.liqViscosity.B = 5.2770E03;
    out.liqViscosity.C = 7.2659E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [273.15, 1.1253E-02; 403.15, 3.7514E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.0470E-07;
    out.vapViscosity.B = 7.8178E-01;
    out.vapViscosity.C = 9.2644E01;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [195.00, 4.3795E-06; 1000.00, 2.1223E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2553E-01;
    out.liqThermalConductivity.B = -2.8271E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [195.00, 1.7040E-01; 449.50, 9.8452E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 7.5670E00;
    out.vapThermalConductivity.B = -2.7146E-01;
    out.vapThermalConductivity.C = - 3.8945E03;
    out.vapThermalConductivity.D = 1.2574E07;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [401.85, 2.1482E-02; 1000.0, 1.1986E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.5638E-02;
    out.SurfaceTension.B = 1.0827E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [195.00, 3.5545E-02; 575.40, 0];
    out.SurfaceTension.eqn = 106;    
    
end


