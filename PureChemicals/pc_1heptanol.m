function out = pc_1heptanol

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'1heptanol','1-heptanol','1-Heptanol','n-heptanol','n-Heptanol', 'heptyl alcohol'}; 
    out.CAS = '111-70-6';
    % components
    out.ChemicalFormular.C = 7;         
    out.ChemicalFormular.H = 16;
    out.ChemicalFormular.O = 1;
%% 1  
    out.MolecularWeight = 116.203;  % kg/kmol    
    out.Tcrit = 632.60;             % K    
    out.pcrit = 3.0580E06;          % Pa 
    out.vcrit = 0.43500;            % m^3/kmol 
    out.CritCompressFactor = 0.253; % -
    
%% 2 
    out.MeltingPoint = 239.15;              % K
    out.TripletPointT = 239.15;             % K
    out.TripletPointp = 1.1773E-02;         % Pa
    out.NormalBoilingPoint = 449.45;        % K
    out.liqMolarVolume = 0.141801;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -3.3640E08;             % J/kmol
    out.IG_GibbsOfFormation = -1.2530E08;            % J/kmol
    out.IG_AbsoluteEntropy =  4.7919E05;             % J/kmolK
    out.HeatFusionAtMeltPoint = 1.8180E07;           % J/kmol
    out.StandardNetHeatOfCombustion = -4.2860E09;    % J/kmol
    out.StandardGibbsOfFormation = -1.4650E08;       % J/kmol
    out.StandardHeatOfFormation = -4.0330E08;        % J/kmol
    out.StandardAbsoluteEntropy = 3.2590E05;         % J/kmolK
    
%% 4 
    out.AcentricFactor = 0.5670;            % -
    out.RadiusOfGyration = 4.3800E-10;      % m
    out.SolubilityParameter = 2.1540E04;    % (J/m^3)^.5 
    out.DipoleMoment = 5.8000E-30;          % C^m
    out.VanDerWaalsVolume = 0.08309;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 1.1680E09;        % m^2/kmol
    out.RefractiveIndex = 1.4223;           % -
    out.FlashPoint = 350.00;                % K
    out.FlammabilityLimits = [1.0 7.2];     % vol-%
    out.FlammabilityLimitTemps = [344 379]; % K
    out.AutoignitionTemp = 555.00;          % K
    
%% 6 
    out.solDensity.A = 8.6800E00;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [239.15, 8.6800E00];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 6.0531E-01;
    out.liqDensity.B = 2.6329E-01;
    out.liqDensity.C = 6.3260E02;
    out.liqDensity.D = 2.7375E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [239.15, 7.4209E00; 632.60, 2.2990E+00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  2.4331E02;
    out.vapPressure.B = -1.7320E04;
    out.vapPressure.C = -3.2183E01;
    out.vapPressure.D =  1.6443E-05;
    out.vapPressure.E =  2.0000E00;
    out.vapPressure.range = [239.15, 1.1773E-02; 632.60, 3.0663E+06];
    out.vapPressure.eqn = 101;
    
    out.HeatOfVap.A = 1.0075E08;
    out.HeatOfVap.B = 6.1052E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [239.15, 7.5393E07; 632.60, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = 1.9300E04;
    out.solHeatCapacity.B = 9.8150E02;
    out.solHeatCapacity.C = -3.0750E00;
    out.solHeatCapacity.D = 6.9530E-03;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [80.00, 8.1700E04; 239.15, 1.7326E05];
    out.solHeatCapacity.eqn = 100;
    
      
    out.liqHeatCapacity.A = 2.4168E06;
    out.liqHeatCapacity.B = -2.6105E04;
    out.liqHeatCapacity.C = 1.1003E02;
    out.liqHeatCapacity.D = -1.9172E-01;
    out.liqHeatCapacity.E = 1.1968E-04;
    out.liqHeatCapacity.range = [239.15, 2.3590E05; 449.45, 3.8773E+05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 1.2215E05;
    out.IG_HeatCapacity.B = 3.9910E05;
    out.IG_HeatCapacity.C = 1.5800E03;
    out.IG_HeatCapacity.D = 2.8350E05;
    out.IG_HeatCapacity.E = 7.1770E02;
    out.IG_HeatCapacity.range = [200.00, 1.3330E05; 1500.00, 4.5346E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.8979E-01;
    out.SecondViralCoef.B = -3.0797E02;
    out.SecondViralCoef.C = -9.8592E07;
    out.SecondViralCoef.D = -2.7797E20;
    out.SecondViralCoef.E = 1.3614E22;
    out.SecondViralCoef.range = [316.30, -6.1382E00; 1500.00, 5.5313E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -7.5245E01;
    out.liqViscosity.B = 5.7567E03;
    out.liqViscosity.C = 8.9146E00;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [239.15, 9.5615E-02; 449.45, 3.4016E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 2.5720E-07;
    out.vapViscosity.B = 6.5020E-01;
    out.vapViscosity.C = 2.4860E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [239.15, 4.4397E-06; 1000.00, 1.8384E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 2.2390E-01;
    out.liqThermalConductivity.B = -2.2600E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [239.15, 1.6985E-01; 573.15, 9.4368E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = -6.1993E-02;
    out.vapThermalConductivity.B = 2.7920E-01;
    out.vapThermalConductivity.C = -3.3360E03;
    out.vapThermalConductivity.D = -1.6420E06;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [504.07, 2.3447E-02; 1000.00, 1.0722E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.3765E-02;
    out.SurfaceTension.B = 1.0974E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [239.15, 3.1928E-02; 632.60, 0.0000E00];
    out.SurfaceTension.eqn = 106;    
    
end


