function out = pc_isohexane

%% Source: Daubert & Danner, "Physical and thermodynamic properties of pure chemicals", Revision 2000

%% 0
    % all alternative names
    out.Names = {'isohexane','iso-hexane','2-methylpentane'}; 
    out.CAS = '107-83-5';
    % components
    out.ChemicalFormular.C = 6;         
    out.ChemicalFormular.H = 14;
    out.ChemicalFormular.O = 0;
%% 1  
    out.MolecularWeight = 86.1772;  % kg/kmol    
    out.Tcrit = 497.7;              % K    
    out.pcrit = 3.0400E6;           % Pa 
    out.vcrit = 0.368;              % m^3/kmol 
    out.CritCompressFactor = 0.27;  % -
    
%% 2 
    out.MeltingPoint = 119.55;              % K
    out.TripletPointT = 119.55;             % K
    out.TripletPointp = 2.0656E-5;          % Pa
    out.NormalBoilingPoint = 333.41;        % K
    out.liqMolarVolume = 0.132925;          % m^3/kmol
    
%% 3
    out.IG_HeatOfFormation = -1.7455E8;             % J/kmol
    out.IG_GibbsOfFormation = -5.3380E6;            % J/kmol
    out.IG_AbsoluteEntropy =  3.8089E5;             % J/kmolK
    out.HeatFusionAtMeltPoint =  6.2680E6;          % J/kmol
    out.StandardNetHeatOfCombustion =-3.8492E9;     % J/kmol
    out.StandardGibbsOfFormation = -8.4999E06;      % J/kmol
    out.StandardHeatOfFormation = -2.0464E08;       % J/kmol
    out.StandardAbsoluteEntropy = 2.9058E05;        % J/kmolK    
    
%% 4 
    out.AcentricFactor = 0.279149;          % -
    out.RadiusOfGyration = 3.7840E-10;      % m
    out.SolubilityParameter = 1.4400E4;     % (J/m^3)^.5 
    out.DipoleMoment = 0;                   % C^m
    out.VanDerWaalsVolume = 0.06825;        % m^3/kmol
    
%% 5    
    out.VanDerWaalsArea = 9.6300E08;        % m^2/kmol
    out.RefractiveIndex = 1.36873;          % -
    out.FlashPoint = 238.0;                 % K
    out.FlammabilityLimits = [1.2 7.0];     % vol-%
    out.FlammabilityLimitTemps = [239 269]; % K
    out.AutoignitionTemp = 579.0;           % K
    
%% 6 
    out.solDensity.A = 1.0800E01;
    out.solDensity.B = 0;
    out.solDensity.C = 0;
    out.solDensity.D = 0;
    out.solDensity.E = 0;
    out.solDensity.range = [119.55,1.0800E01];
    out.solDensity.eqn = 100;
    
    out.liqDensity.A = 7.2701E-01;
    out.liqDensity.B = 2.6754E-01;
    out.liqDensity.C = 4.9770E02;
    out.liqDensity.D = 2.8268E-01;
    out.liqDensity.E = 0;
    out.liqDensity.range = [119.55, 9.2041E00; 497.70, 2.7174E00];
    out.liqDensity.eqn = 105;
    
    out.vapPressure.A =  5.3579E01;
    out.vapPressure.B = -5.0412E03;
    out.vapPressure.C = -4.6404E00;
    out.vapPressure.D =  1.9443E-17;
    out.vapPressure.E =  6.0000E00;
    out.vapPressure.range = [119.55, 2.0656E-05; 497.70, 3.0439E06];
    out.vapPressure.eqn = 101;  
    
    out.HeatOfVap.A = 4.2522E07;
    out.HeatOfVap.B = 3.8070E-01;
    out.HeatOfVap.C = 0;
    out.HeatOfVap.D = 0;
    out.HeatOfVap.E = 0;
    out.HeatOfVap.range = [119.55, 3.8300E07; 497.70, 0];
    out.HeatOfVap.eqn = 106;     

    out.solHeatCapacity.A = -1.6127E04;
    out.solHeatCapacity.B = 1.6109E03;
    out.solHeatCapacity.C = -9.0041E00;
    out.solHeatCapacity.D = 2.6825E-02;
    out.solHeatCapacity.E = 0;
    out.solHeatCapacity.range = [8.00,1.0108E04; 119.53,9.3590E04];
    out.solHeatCapacity.eqn = 100;
    
    out.liqHeatCapacity.A = 1.4222E05;
    out.liqHeatCapacity.B = -4.7830E01;
    out.liqHeatCapacity.C = 7.3900E-01;
    out.liqHeatCapacity.D = 0;
    out.liqHeatCapacity.E = 0;
    out.liqHeatCapacity.range = [119.55, 1.4706E05; 333.41,2.0842E05];
    out.liqHeatCapacity.eqn = 100;    
    
    out.IG_HeatCapacity.A = 9.0300E04;
    out.IG_HeatCapacity.B = 3.8010E05;
    out.IG_HeatCapacity.C = 1.6020E03;
    out.IG_HeatCapacity.D = 2.4530E05;
    out.IG_HeatCapacity.E = -6.9160E02;
    out.IG_HeatCapacity.range = [200.00,1.0192E05; 1500.00, 3.9617E05];
    out.IG_HeatCapacity.eqn = 107;
    
    out.SecondViralCoef.A = 2.0487E-01;
    out.SecondViralCoef.B = -2.0373E02;
    out.SecondViralCoef.C = -3.1708E07;
    out.SecondViralCoef.D = -4.4141E18;
    out.SecondViralCoef.E = -3.7253E21;
    out.SecondViralCoef.range = [248.73,-3.3958E00; 1500.00, 5.9655E-02];
    out.SecondViralCoef.eqn = 104;  
    
    out.liqViscosity.A = -1.2860E01;
    out.liqViscosity.B = 9.4691E02;
    out.liqViscosity.C = 2.6191E-01;
    out.liqViscosity.D = 0;
    out.liqViscosity.E = 0;
    out.liqViscosity.range = [119.55, 2.5060E-02; 333.41, 2.0380E-04];
    out.liqViscosity.eqn = 101;
    
    out.vapViscosity.A = 1.1164E-06;
    out.vapViscosity.B = 4.5370E-01;
    out.vapViscosity.C = 3.7474E02;
    out.vapViscosity.D = 0;
    out.vapViscosity.E = 0;
    out.vapViscosity.range = [119.55, 2.3658E-06; 1000.00, 1.8651E-05];
    out.vapViscosity.eqn = 102;
    
    out.liqThermalConductivity.A = 1.9334E-01;
    out.liqThermalConductivity.B = -2.8038E-04;
    out.liqThermalConductivity.C = 0;
    out.liqThermalConductivity.D = 0;
    out.liqThermalConductivity.E = 0;
    out.liqThermalConductivity.range = [119.55, 1.5982E-01; 389.25, 8.4202E-02];
    out.liqThermalConductivity.eqn = 100;
       
    out.vapThermalConductivity.A = 6.1119E-05;
    out.vapThermalConductivity.B = 1.0861E00;
    out.vapThermalConductivity.C = -5.9592E01;
    out.vapThermalConductivity.D = 1.4126E05;
    out.vapThermalConductivity.E = 0;
    out.vapThermalConductivity.range = [333.41, 1.6063E-02; 1000.0, 1.0242E-01];
    out.vapThermalConductivity.eqn = 102; 
    
    out.SurfaceTension.A = 5.0048E-02;
    out.SurfaceTension.B = 1.1899E00;
    out.SurfaceTension.C = 0;
    out.SurfaceTension.D = 0;
    out.SurfaceTension.E = 0;
    out.SurfaceTension.range = [119.55, 3.6093E-02; 497.70, 0];
    out.SurfaceTension.eqn = 106;    
    
end


