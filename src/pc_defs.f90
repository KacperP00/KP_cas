! Pure Chemicals Database adapted from matlab code of J. Palmer (WSA)
! Written by Abhishek Deshmukh 06/2017
module pc_defs
  use precision

  implicit none

  ! Definition of pure chemical object
  type pc_t

     ! Temperature and pressure
     real(WP), pointer :: p, T

     ! Name
     character(len=128), pointer :: name, CAS

     ! Constants
     type(ChemicalFormula_t), pointer :: ChemicalFormula
     
     real(WP), pointer :: MolecularWeight, Tcrit, pcrit, vcrit, CritCompressFactor
     
     real(WP), pointer :: MeltingPoint, TripletPointT, TripletPointp, NormalBoilingPoint, liqMolarVolume

     real(WP), pointer :: IG_HeatOfFormation, IG_GibbsOfFormation, IG_AbsoluteEntropy, HeatFusionAtMeltPoint, &
                          StandardNetHeatOfCombustion, StandardGibbsOfFormation, StandardHeatOfFormation, StandardAbsoluteEntropy

     real(WP), pointer :: AcentricFactor, RadiusOfGyration, SolubilityParameter, DipoleMoment, VanDerWaalsVolume

     real(WP), pointer :: VanDerWaalsArea, RefractiveIndex, FlashPoint, AutoignitionTemp

     real(WP), dimension(:), pointer :: FlammabilityLimits, FlammabilityLimitTemps
     
     real(WP), pointer :: solDensity, liqDensity, IG_vapDensity, HeatOfVap, solHeatCapacity, liqHeatCapacity, &
                          IG_HeatCapacity, SecondViralCoef, DiffusionCoefficientWilkeLee, DiffusionCoefficientFuller, &
			  DiffusionCoefficientTheory, DiffusionCoefficientHighPres

     type(pc_prop_t), pointer :: solDensityMol, liqDensityMol, vapPressure, HeatOfVapMol, solHeatCapacityMol, &
                                 liqHeatCapacityMol, IG_HeatCapacityMol, SecondViralCoefMol, liqViscosity, vapViscosity, &
                                 liqThermalConductivity, vapThermalConductivity, SurfaceTension
     
  end type pc_t

  type ChemicalFormula_t
     
     integer, pointer :: C, H, O

  end type ChemicalFormula_t

  type pc_prop_t
     
     real(WP), pointer :: A, B, C, D, E
     real(WP), dimension(:,:), pointer :: range
     integer, pointer :: eqn
     real(WP), pointer :: val

  end type pc_prop_t

contains

  ! Allocate pc with default values
  subroutine nullify_pc(pc)
    implicit none

    ! ---------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------

    nullify(pc%p,pc%T)

    nullify(pc%name)
    nullify(pc%CAS)

    nullify(pc%ChemicalFormula%C)
    nullify(pc%ChemicalFormula%H)
    nullify(pc%ChemicalFormula%O)
    nullify(pc%ChemicalFormula)

    nullify(pc%MolecularWeight)
    nullify(pc%Tcrit)
    nullify(pc%pcrit)
    nullify(pc%vcrit)
    nullify(pc%CritCompressFactor)

    nullify(pc%MeltingPoint)
    nullify(pc%TripletPointT)
    nullify(pc%TripletPointp)
    nullify(pc%NormalBoilingPoint)
    nullify(pc%liqMolarVolume)

    nullify(pc%IG_HeatOfFormation)
    nullify(pc%IG_GibbsOfFormation)
    nullify(pc%IG_AbsoluteEntropy)
    nullify(pc%HeatFusionAtMeltPoint)
    nullify(pc%StandardNetHeatOfCombustion)
    nullify(pc%StandardGibbsOfFormation)
    nullify(pc%StandardHeatOfFormation)
    nullify(pc%StandardAbsoluteEntropy)

    nullify(pc%AcentricFactor)
    nullify(pc%RadiusOfGyration)
    nullify(pc%SolubilityParameter)
    nullify(pc%DipoleMoment)
    nullify(pc%VanDerWaalsVolume)

    nullify(pc%VanDerWaalsArea)
    nullify(pc%RefractiveIndex)
    nullify(pc%FlashPoint)
    nullify(pc%AutoignitionTemp)

    nullify(pc%FlammabilityLimits)
    nullify(pc%FlammabilityLimitTemps)

    nullify(pc%solDensityMol%A)
    nullify(pc%solDensityMol%B)
    nullify(pc%solDensityMol%C)
    nullify(pc%solDensityMol%D)
    nullify(pc%solDensityMol%E)
    nullify(pc%solDensityMol%range)
    nullify(pc%solDensityMol%eqn)
    nullify(pc%solDensityMol%val)
    nullify(pc%solDensityMol)
    nullify(pc%solDensity)

    nullify(pc%liqDensityMol%A)
    nullify(pc%liqDensityMol%B)
    nullify(pc%liqDensityMol%C)
    nullify(pc%liqDensityMol%D)
    nullify(pc%liqDensityMol%E)
    nullify(pc%liqDensityMol%range)
    nullify(pc%liqDensityMol%eqn)
    nullify(pc%liqDensityMol%val)
    nullify(pc%liqDensityMol)
    nullify(pc%liqDensity)

    nullify(pc%IG_vapDensity)

    nullify(pc%vapPressure%A)
    nullify(pc%vapPressure%B)
    nullify(pc%vapPressure%C)
    nullify(pc%vapPressure%D)
    nullify(pc%vapPressure%E)
    nullify(pc%vapPressure%range)
    nullify(pc%vapPressure%eqn)
    nullify(pc%vapPressure%val)
    nullify(pc%vapPressure)

    nullify(pc%HeatOfVapMol%A)
    nullify(pc%HeatOfVapMol%B)
    nullify(pc%HeatOfVapMol%C)
    nullify(pc%HeatOfVapMol%D)
    nullify(pc%HeatOfVapMol%E)
    nullify(pc%HeatOfVapMol%range)
    nullify(pc%HeatOfVapMol%eqn)
    nullify(pc%HeatOfVapMol%val)
    nullify(pc%HeatOfVapMol)

    nullify(pc%solHeatCapacityMol%A)
    nullify(pc%solHeatCapacityMol%B)
    nullify(pc%solHeatCapacityMol%C)
    nullify(pc%solHeatCapacityMol%D)
    nullify(pc%solHeatCapacityMol%E)
    nullify(pc%solHeatCapacityMol%range)
    nullify(pc%solHeatCapacityMol%eqn)
    nullify(pc%solHeatCapacityMol%val)
    nullify(pc%solHeatCapacityMol)

    nullify(pc%liqHeatCapacityMol%A)
    nullify(pc%liqHeatCapacityMol%B)
    nullify(pc%liqHeatCapacityMol%C)
    nullify(pc%liqHeatCapacityMol%D)
    nullify(pc%liqHeatCapacityMol%E)
    nullify(pc%liqHeatCapacityMol%range)
    nullify(pc%liqHeatCapacityMol%eqn)
    nullify(pc%liqHeatCapacityMol%val)
    nullify(pc%liqHeatCapacityMol)

    nullify(pc%IG_HeatCapacityMol%A)
    nullify(pc%IG_HeatCapacityMol%B)
    nullify(pc%IG_HeatCapacityMol%C)
    nullify(pc%IG_HeatCapacityMol%D)
    nullify(pc%IG_HeatCapacityMol%E)
    nullify(pc%IG_HeatCapacityMol%range)
    nullify(pc%IG_HeatCapacityMol%eqn)
    nullify(pc%IG_HeatCapacityMol%val)
    nullify(pc%IG_HeatCapacityMol)

    nullify(pc%SecondViralCoefMol%A)
    nullify(pc%SecondViralCoefMol%B)
    nullify(pc%SecondViralCoefMol%C)
    nullify(pc%SecondViralCoefMol%D)
    nullify(pc%SecondViralCoefMol%E)
    nullify(pc%SecondViralCoefMol%range)
    nullify(pc%SecondViralCoefMol%eqn)
    nullify(pc%SecondViralCoefMol%val)
    nullify(pc%SecondViralCoefMol)

    nullify(pc%liqViscosity%A)
    nullify(pc%liqViscosity%B)
    nullify(pc%liqViscosity%C)
    nullify(pc%liqViscosity%D)
    nullify(pc%liqViscosity%E)
    nullify(pc%liqViscosity%range)
    nullify(pc%liqViscosity%eqn)
    nullify(pc%liqViscosity%val)
    nullify(pc%liqViscosity)

    nullify(pc%vapViscosity%A)
    nullify(pc%vapViscosity%B)
    nullify(pc%vapViscosity%C)
    nullify(pc%vapViscosity%D)
    nullify(pc%vapViscosity%E)
    nullify(pc%vapViscosity%range)
    nullify(pc%vapViscosity%eqn)
    nullify(pc%vapViscosity%val)
    nullify(pc%vapViscosity)

    nullify(pc%liqThermalConductivity%A)
    nullify(pc%liqThermalConductivity%B)
    nullify(pc%liqThermalConductivity%C)
    nullify(pc%liqThermalConductivity%D)
    nullify(pc%liqThermalConductivity%E)
    nullify(pc%liqThermalConductivity%range)
    nullify(pc%liqThermalConductivity%eqn)
    nullify(pc%liqThermalConductivity%val)
    nullify(pc%liqThermalConductivity)

    nullify(pc%vapThermalConductivity%A)
    nullify(pc%vapThermalConductivity%B)
    nullify(pc%vapThermalConductivity%C)
    nullify(pc%vapThermalConductivity%D)
    nullify(pc%vapThermalConductivity%E)
    nullify(pc%vapThermalConductivity%range)
    nullify(pc%vapThermalConductivity%eqn)
    nullify(pc%vapThermalConductivity%val)
    nullify(pc%vapThermalConductivity)

    nullify(pc%SurfaceTension%A)
    nullify(pc%SurfaceTension%B)
    nullify(pc%SurfaceTension%C)
    nullify(pc%SurfaceTension%D)
    nullify(pc%SurfaceTension%E)
    nullify(pc%SurfaceTension%range)
    nullify(pc%SurfaceTension%eqn)
    nullify(pc%SurfaceTension%val)  
    nullify(pc%SurfaceTension)

    nullify(pc%DiffusionCoefficientWilkeLee)
    nullify(pc%DiffusionCoefficientFuller)
    nullify(pc%DiffusionCoefficientTheory)
    nullify(pc%DiffusionCoefficientHighPres)

    nullify(pc)

  end subroutine nullify_pc

  ! Allocate pc with default values
  subroutine allocate_pc(pc)
    implicit none

    ! ---------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------

    !call nullify_pc(pc)

    allocate(pc)
    allocate(pc%p); pc%p = -9999.0_WP
    allocate(pc%T); pc%T = -9999.0_WP
    allocate(pc%name); pc%name = 'pc'
    allocate(pc%CAS); pc%CAS = 'CAS'

    allocate(pc%ChemicalFormula); 
    allocate(pc%ChemicalFormula%C); pc%ChemicalFormula%C = -9999.0_WP
    allocate(pc%ChemicalFormula%H); pc%ChemicalFormula%H = -9999.0_WP
    allocate(pc%ChemicalFormula%O); pc%ChemicalFormula%O = -9999.0_WP

    allocate(pc%MolecularWeight); pc%MolecularWeight = -9999.0_WP
    allocate(pc%Tcrit); pc%Tcrit = -9999.0_WP
    allocate(pc%pcrit); pc%pcrit = -9999.0_WP
    allocate(pc%vcrit); pc%vcrit = -9999.0_WP
    allocate(pc%CritCompressFactor); pc%CritCompressFactor = -9999.0_WP

    allocate(pc%MeltingPoint); pc%MeltingPoint = -9999.0_WP
    allocate(pc%TripletPointT); pc%TripletPointT = -9999.0_WP
    allocate(pc%TripletPointp); pc%TripletPointp = -9999.0_WP
    allocate(pc%NormalBoilingPoint); pc%NormalBoilingPoint = -9999.0_WP
    allocate(pc%liqMolarVolume); pc%liqMolarVolume = -9999.0_WP

    allocate(pc%IG_HeatOfFormation); pc%IG_HeatOfFormation = -9999.0_WP
    allocate(pc%IG_GibbsOfFormation); pc%IG_GibbsOfFormation = -9999.0_WP
    allocate(pc%IG_AbsoluteEntropy); pc%IG_AbsoluteEntropy = -9999.0_WP
    allocate(pc%HeatFusionAtMeltPoint); pc%HeatFusionAtMeltPoint = -9999.0_WP
    allocate(pc%StandardNetHeatOfCombustion); pc%StandardNetHeatOfCombustion = -9999.0_WP
    allocate(pc%StandardGibbsOfFormation); pc%StandardGibbsOfFormation = -9999.0_WP
    allocate(pc%StandardHeatOfFormation); pc%StandardHeatOfFormation = -9999.0_WP
    allocate(pc%StandardAbsoluteEntropy); pc%StandardAbsoluteEntropy = -9999.0_WP

    allocate(pc%AcentricFactor); pc%AcentricFactor = -9999.0_WP
    allocate(pc%RadiusOfGyration); pc%RadiusOfGyration = -9999.0_WP
    allocate(pc%SolubilityParameter); pc%SolubilityParameter = -9999.0_WP
    allocate(pc%DipoleMoment); pc%DipoleMoment = -9999.0_WP
    allocate(pc%VanDerWaalsVolume); pc%VanDerWaalsVolume = -9999.0_WP

    allocate(pc%VanDerWaalsArea); pc%VanDerWaalsArea = -9999.0_WP
    allocate(pc%RefractiveIndex); pc%RefractiveIndex = -9999.0_WP
    allocate(pc%FlashPoint); pc%FlashPoint = -9999.0_WP
    allocate(pc%AutoignitionTemp); pc%AutoignitionTemp = -9999.0_WP

    allocate(pc%FlammabilityLimits(2)); pc%FlammabilityLimits(2) = -9999.0_WP
    allocate(pc%FlammabilityLimitTemps(2)); pc%FlammabilityLimitTemps(2) = -9999.0_WP

    allocate(pc%solDensity); pc%solDensity = -9999.0_WP
    allocate(pc%solDensityMol)
    allocate(pc%solDensityMol%A); pc%solDensityMol%A = -9999.0_WP
    allocate(pc%solDensityMol%B); pc%solDensityMol%B = -9999.0_WP
    allocate(pc%solDensityMol%C); pc%solDensityMol%C = -9999.0_WP
    allocate(pc%solDensityMol%D); pc%solDensityMol%D = -9999.0_WP
    allocate(pc%solDensityMol%E); pc%solDensityMol%E = -9999.0_WP
    allocate(pc%solDensityMol%range(2,2)); pc%solDensityMol%range = -9999.0_WP
    allocate(pc%solDensityMol%eqn); pc%solDensityMol%eqn = -9999
    allocate(pc%solDensityMol%val); pc%solDensityMol%val = -9999.0_WP

    allocate(pc%liqDensity); pc%liqDensity = -9999.0_WP
    allocate(pc%liqDensityMol)
    allocate(pc%liqDensityMol%A); pc%liqDensityMol%A = -9999.0_WP
    allocate(pc%liqDensityMol%B); pc%liqDensityMol%B = -9999.0_WP
    allocate(pc%liqDensityMol%C); pc%liqDensityMol%C = -9999.0_WP
    allocate(pc%liqDensityMol%D); pc%liqDensityMol%D = -9999.0_WP
    allocate(pc%liqDensityMol%E); pc%liqDensityMol%E = -9999.0_WP
    allocate(pc%liqDensityMol%range(2,2)); pc%liqDensityMol%range = -9999.0_WP
    allocate(pc%liqDensityMol%eqn); pc%liqDensityMol%eqn = -9999
    allocate(pc%liqDensityMol%val); pc%liqDensityMol%val = -9999.0_WP

    allocate(pc%IG_vapDensity); pc%IG_vapDensity = -9999.0_WP

    allocate(pc%vapPressure)
    allocate(pc%vapPressure%A); pc%vapPressure%A = -9999.0_WP
    allocate(pc%vapPressure%B); pc%vapPressure%B = -9999.0_WP
    allocate(pc%vapPressure%C); pc%vapPressure%C = -9999.0_WP
    allocate(pc%vapPressure%D); pc%vapPressure%D = -9999.0_WP
    allocate(pc%vapPressure%E); pc%vapPressure%E = -9999.0_WP
    allocate(pc%vapPressure%range(2,2)); pc%vapPressure%range = -9999.0_WP
    allocate(pc%vapPressure%eqn); pc%vapPressure%eqn = -9999
    allocate(pc%vapPressure%val); pc%vapPressure%val = -9999.0_WP

    allocate(pc%HeatOfVap); pc%HeatOfVap =-9999.0_WP
    allocate(pc%HeatOfVapMol)
    allocate(pc%HeatOfVapMol%A); pc%HeatOfVapMol%A = -9999.0_WP
    allocate(pc%HeatOfVapMol%B); pc%HeatOfVapMol%B = -9999.0_WP
    allocate(pc%HeatOfVapMol%C); pc%HeatOfVapMol%C = -9999.0_WP
    allocate(pc%HeatOfVapMol%D); pc%HeatOfVapMol%D = -9999.0_WP
    allocate(pc%HeatOfVapMol%E); pc%HeatOfVapMol%E = -9999.0_WP
    allocate(pc%HeatOfVapMol%range(2,2)); pc%HeatOfVapMol%range = -9999.0_WP
    allocate(pc%HeatOfVapMol%eqn); pc%HeatOfVapMol%eqn = -9999
    allocate(pc%HeatOfVapMol%val); pc%HeatOfVapMol%val = -9999.0_WP

    allocate(pc%solHeatCapacity); pc%solHeatCapacity =-9999.0_WP
    allocate(pc%solHeatCapacityMol)
    allocate(pc%solHeatCapacityMol%A); pc%solHeatCapacityMol%A = -9999.0_WP
    allocate(pc%solHeatCapacityMol%B); pc%solHeatCapacityMol%B = -9999.0_WP
    allocate(pc%solHeatCapacityMol%C); pc%solHeatCapacityMol%C = -9999.0_WP
    allocate(pc%solHeatCapacityMol%D); pc%solHeatCapacityMol%D = -9999.0_WP
    allocate(pc%solHeatCapacityMol%E); pc%solHeatCapacityMol%E = -9999.0_WP
    allocate(pc%solHeatCapacityMol%range(2,2)); pc%solHeatCapacityMol%range = -9999.0_WP
    allocate(pc%solHeatCapacityMol%eqn); pc%solHeatCapacityMol%eqn = -9999
    allocate(pc%solHeatCapacityMol%val); pc%solHeatCapacityMol%val = -9999.0_WP

    allocate(pc%liqHeatCapacity); pc%liqHeatCapacity =-9999.0_WP
    allocate(pc%liqHeatCapacityMol)
    allocate(pc%liqHeatCapacityMol%A); pc%liqHeatCapacityMol%A = -9999.0_WP
    allocate(pc%liqHeatCapacityMol%B); pc%liqHeatCapacityMol%B = -9999.0_WP
    allocate(pc%liqHeatCapacityMol%C); pc%liqHeatCapacityMol%C = -9999.0_WP
    allocate(pc%liqHeatCapacityMol%D); pc%liqHeatCapacityMol%D = -9999.0_WP
    allocate(pc%liqHeatCapacityMol%E); pc%liqHeatCapacityMol%E = -9999.0_WP
    allocate(pc%liqHeatCapacityMol%range(2,2)); pc%liqHeatCapacityMol%range = -9999.0_WP
    allocate(pc%liqHeatCapacityMol%eqn); pc%liqHeatCapacityMol%eqn = -9999
    allocate(pc%liqHeatCapacityMol%val); pc%liqHeatCapacityMol%val = -9999.0_WP

    allocate(pc%IG_HeatCapacity); pc%IG_HeatCapacity =-9999.0_WP
    allocate(pc%IG_HeatCapacityMol)
    allocate(pc%IG_HeatCapacityMol%A); pc%IG_HeatCapacityMol%A = -9999.0_WP
    allocate(pc%IG_HeatCapacityMol%B); pc%IG_HeatCapacityMol%B = -9999.0_WP
    allocate(pc%IG_HeatCapacityMol%C); pc%IG_HeatCapacityMol%C = -9999.0_WP
    allocate(pc%IG_HeatCapacityMol%D); pc%IG_HeatCapacityMol%D = -9999.0_WP
    allocate(pc%IG_HeatCapacityMol%E); pc%IG_HeatCapacityMol%E = -9999.0_WP
    allocate(pc%IG_HeatCapacityMol%range(2,2)); pc%IG_HeatCapacityMol%range = -9999.0_WP
    allocate(pc%IG_HeatCapacityMol%eqn); pc%IG_HeatCapacityMol%eqn = -9999
    allocate(pc%IG_HeatCapacityMol%val); pc%IG_HeatCapacityMol%val = -9999.0_WP

    allocate(pc%SecondViralCoef); pc%SecondViralCoef =-9999.0_WP
    allocate(pc%SecondViralCoefMol)
    allocate(pc%SecondViralCoefMol%A); pc%SecondViralCoefMol%A = -9999.0_WP
    allocate(pc%SecondViralCoefMol%B); pc%SecondViralCoefMol%B = -9999.0_WP
    allocate(pc%SecondViralCoefMol%C); pc%SecondViralCoefMol%C = -9999.0_WP
    allocate(pc%SecondViralCoefMol%D); pc%SecondViralCoefMol%D = -9999.0_WP
    allocate(pc%SecondViralCoefMol%E); pc%SecondViralCoefMol%E = -9999.0_WP
    allocate(pc%SecondViralCoefMol%range(2,2)); pc%SecondViralCoefMol%range = -9999.0_WP
    allocate(pc%SecondViralCoefMol%eqn); pc%SecondViralCoefMol%eqn = -9999
    allocate(pc%SecondViralCoefMol%val); pc%SecondViralCoefMol%val = -9999.0_WP

    allocate(pc%liqViscosity)
    allocate(pc%liqViscosity%A); pc%liqViscosity%A = -9999.0_WP
    allocate(pc%liqViscosity%B); pc%liqViscosity%B = -9999.0_WP
    allocate(pc%liqViscosity%C); pc%liqViscosity%C = -9999.0_WP
    allocate(pc%liqViscosity%D); pc%liqViscosity%D = -9999.0_WP
    allocate(pc%liqViscosity%E); pc%liqViscosity%E = -9999.0_WP
    allocate(pc%liqViscosity%range(2,2)); pc%liqViscosity%range = -9999.0_WP
    allocate(pc%liqViscosity%eqn); pc%liqViscosity%eqn = -9999
    allocate(pc%liqViscosity%val); pc%liqViscosity%val = -9999.0_WP

    allocate(pc%vapViscosity)
    allocate(pc%vapViscosity%A); pc%vapViscosity%A = -9999.0_WP
    allocate(pc%vapViscosity%B); pc%vapViscosity%B = -9999.0_WP
    allocate(pc%vapViscosity%C); pc%vapViscosity%C = -9999.0_WP
    allocate(pc%vapViscosity%D); pc%vapViscosity%D = -9999.0_WP
    allocate(pc%vapViscosity%E); pc%vapViscosity%E = -9999.0_WP
    allocate(pc%vapViscosity%range(2,2)); pc%vapViscosity%range = -9999.0_WP
    allocate(pc%vapViscosity%eqn); pc%vapViscosity%eqn = -9999
    allocate(pc%vapViscosity%val); pc%vapViscosity%val = -9999.0_WP

    allocate(pc%liqThermalConductivity)
    allocate(pc%liqThermalConductivity%A); pc%liqThermalConductivity%A = -9999.0_WP
    allocate(pc%liqThermalConductivity%B); pc%liqThermalConductivity%B = -9999.0_WP
    allocate(pc%liqThermalConductivity%C); pc%liqThermalConductivity%C = -9999.0_WP
    allocate(pc%liqThermalConductivity%D); pc%liqThermalConductivity%D = -9999.0_WP
    allocate(pc%liqThermalConductivity%E); pc%liqThermalConductivity%E = -9999.0_WP
    allocate(pc%liqThermalConductivity%range(2,2)); pc%liqThermalConductivity%range = -9999.0_WP
    allocate(pc%liqThermalConductivity%eqn); pc%liqThermalConductivity%eqn = -9999
    allocate(pc%liqThermalConductivity%val); pc%liqThermalConductivity%val = -9999.0_WP

    allocate(pc%vapThermalConductivity)
    allocate(pc%vapThermalConductivity%A); pc%vapThermalConductivity%A = -9999.0_WP
    allocate(pc%vapThermalConductivity%B); pc%vapThermalConductivity%B = -9999.0_WP
    allocate(pc%vapThermalConductivity%C); pc%vapThermalConductivity%C = -9999.0_WP
    allocate(pc%vapThermalConductivity%D); pc%vapThermalConductivity%D = -9999.0_WP
    allocate(pc%vapThermalConductivity%E); pc%vapThermalConductivity%E = -9999.0_WP
    allocate(pc%vapThermalConductivity%range(2,2)); pc%vapThermalConductivity%range = -9999.0_WP
    allocate(pc%vapThermalConductivity%eqn); pc%vapThermalConductivity%eqn = -9999
    allocate(pc%vapThermalConductivity%val); pc%vapThermalConductivity%val = -9999.0_WP

    allocate(pc%SurfaceTension)
    allocate(pc%SurfaceTension%A); pc%SurfaceTension%A = -9999.0_WP
    allocate(pc%SurfaceTension%B); pc%SurfaceTension%B = -9999.0_WP
    allocate(pc%SurfaceTension%C); pc%SurfaceTension%C = -9999.0_WP
    allocate(pc%SurfaceTension%D); pc%SurfaceTension%D = -9999.0_WP
    allocate(pc%SurfaceTension%E); pc%SurfaceTension%E = -9999.0_WP
    allocate(pc%SurfaceTension%range(2,2)); pc%SurfaceTension%range = -9999.0_WP
    allocate(pc%SurfaceTension%eqn); pc%SurfaceTension%eqn = -9999
    allocate(pc%SurfaceTension%val); pc%SurfaceTension%val = -9999.0_WP

    allocate(pc%DiffusionCoefficientWilkeLee); pc%DiffusionCoefficientWilkeLee =-9999.0_WP
    allocate(pc%DiffusionCoefficientFuller); pc%DiffusionCoefficientFuller =-9999.0_WP
    allocate(pc%DiffusionCoefficientTheory); pc%DiffusionCoefficientTheory =-9999.0_WP
    allocate(pc%DiffusionCoefficientHighPres); pc%DiffusionCoefficientHighPres =-9999.0_WP

  end subroutine allocate_pc

  subroutine deallocate_pc(pc)
    implicit none
    ! ---------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    ! ---------------------------------

    deallocate(pc%p,pc%T)

    deallocate(pc%name)
    deallocate(pc%CAS)

    deallocate(pc%ChemicalFormula%C)
    deallocate(pc%ChemicalFormula%H)
    deallocate(pc%ChemicalFormula%O)
    deallocate(pc%ChemicalFormula)

    deallocate(pc%MolecularWeight)
    deallocate(pc%Tcrit)
    deallocate(pc%pcrit)
    deallocate(pc%vcrit)
    deallocate(pc%CritCompressFactor)

    deallocate(pc%MeltingPoint)
    deallocate(pc%TripletPointT)
    deallocate(pc%TripletPointp)
    deallocate(pc%NormalBoilingPoint)
    deallocate(pc%liqMolarVolume)

    deallocate(pc%IG_HeatOfFormation)
    deallocate(pc%IG_GibbsOfFormation)
    deallocate(pc%IG_AbsoluteEntropy)
    deallocate(pc%HeatFusionAtMeltPoint)
    deallocate(pc%StandardNetHeatOfCombustion)
    deallocate(pc%StandardGibbsOfFormation)
    deallocate(pc%StandardHeatOfFormation)
    deallocate(pc%StandardAbsoluteEntropy)

    deallocate(pc%AcentricFactor)
    deallocate(pc%RadiusOfGyration)
    deallocate(pc%SolubilityParameter)
    deallocate(pc%DipoleMoment)
    deallocate(pc%VanDerWaalsVolume)

    deallocate(pc%VanDerWaalsArea)
    deallocate(pc%RefractiveIndex)
    deallocate(pc%FlashPoint)
    deallocate(pc%AutoignitionTemp)

    deallocate(pc%FlammabilityLimits)
    deallocate(pc%FlammabilityLimitTemps)

    deallocate(pc%solDensityMol%A)
    deallocate(pc%solDensityMol%B)
    deallocate(pc%solDensityMol%C)
    deallocate(pc%solDensityMol%D)
    deallocate(pc%solDensityMol%E)
    deallocate(pc%solDensityMol%range)
    deallocate(pc%solDensityMol%eqn)
    deallocate(pc%solDensityMol%val)
    deallocate(pc%solDensityMol)
    deallocate(pc%solDensity)

    deallocate(pc%liqDensityMol%A)
    deallocate(pc%liqDensityMol%B)
    deallocate(pc%liqDensityMol%C)
    deallocate(pc%liqDensityMol%D)
    deallocate(pc%liqDensityMol%E)
    deallocate(pc%liqDensityMol%range)
    deallocate(pc%liqDensityMol%eqn)
    deallocate(pc%liqDensityMol%val)
    deallocate(pc%liqDensityMol)
    deallocate(pc%liqDensity)

    deallocate(pc%IG_vapDensity)

    deallocate(pc%vapPressure%A)
    deallocate(pc%vapPressure%B)
    deallocate(pc%vapPressure%C)
    deallocate(pc%vapPressure%D)
    deallocate(pc%vapPressure%E)
    deallocate(pc%vapPressure%range)
    deallocate(pc%vapPressure%eqn)
    deallocate(pc%vapPressure%val)
    deallocate(pc%vapPressure)

    deallocate(pc%HeatOfVapMol%A)
    deallocate(pc%HeatOfVapMol%B)
    deallocate(pc%HeatOfVapMol%C)
    deallocate(pc%HeatOfVapMol%D)
    deallocate(pc%HeatOfVapMol%E)
    deallocate(pc%HeatOfVapMol%range)
    deallocate(pc%HeatOfVapMol%eqn)
    deallocate(pc%HeatOfVapMol%val)
    deallocate(pc%HeatOfVapMol)

    deallocate(pc%solHeatCapacityMol%A)
    deallocate(pc%solHeatCapacityMol%B)
    deallocate(pc%solHeatCapacityMol%C)
    deallocate(pc%solHeatCapacityMol%D)
    deallocate(pc%solHeatCapacityMol%E)
    deallocate(pc%solHeatCapacityMol%range)
    deallocate(pc%solHeatCapacityMol%eqn)
    deallocate(pc%solHeatCapacityMol%val)
    deallocate(pc%solHeatCapacityMol)

    deallocate(pc%liqHeatCapacityMol%A)
    deallocate(pc%liqHeatCapacityMol%B)
    deallocate(pc%liqHeatCapacityMol%C)
    deallocate(pc%liqHeatCapacityMol%D)
    deallocate(pc%liqHeatCapacityMol%E)
    deallocate(pc%liqHeatCapacityMol%range)
    deallocate(pc%liqHeatCapacityMol%eqn)
    deallocate(pc%liqHeatCapacityMol%val)
    deallocate(pc%liqHeatCapacityMol)

    deallocate(pc%IG_HeatCapacityMol%A)
    deallocate(pc%IG_HeatCapacityMol%B)
    deallocate(pc%IG_HeatCapacityMol%C)
    deallocate(pc%IG_HeatCapacityMol%D)
    deallocate(pc%IG_HeatCapacityMol%E)
    deallocate(pc%IG_HeatCapacityMol%range)
    deallocate(pc%IG_HeatCapacityMol%eqn)
    deallocate(pc%IG_HeatCapacityMol%val)
    deallocate(pc%IG_HeatCapacityMol)

    deallocate(pc%SecondViralCoefMol%A)
    deallocate(pc%SecondViralCoefMol%B)
    deallocate(pc%SecondViralCoefMol%C)
    deallocate(pc%SecondViralCoefMol%D)
    deallocate(pc%SecondViralCoefMol%E)
    deallocate(pc%SecondViralCoefMol%range)
    deallocate(pc%SecondViralCoefMol%eqn)
    deallocate(pc%SecondViralCoefMol%val)
    deallocate(pc%SecondViralCoefMol)

    deallocate(pc%liqViscosity%A)
    deallocate(pc%liqViscosity%B)
    deallocate(pc%liqViscosity%C)
    deallocate(pc%liqViscosity%D)
    deallocate(pc%liqViscosity%E)
    deallocate(pc%liqViscosity%range)
    deallocate(pc%liqViscosity%eqn)
    deallocate(pc%liqViscosity%val)
    deallocate(pc%liqViscosity)

    deallocate(pc%vapViscosity%A)
    deallocate(pc%vapViscosity%B)
    deallocate(pc%vapViscosity%C)
    deallocate(pc%vapViscosity%D)
    deallocate(pc%vapViscosity%E)
    deallocate(pc%vapViscosity%range)
    deallocate(pc%vapViscosity%eqn)
    deallocate(pc%vapViscosity%val)
    deallocate(pc%vapViscosity)

    deallocate(pc%liqThermalConductivity%A)
    deallocate(pc%liqThermalConductivity%B)
    deallocate(pc%liqThermalConductivity%C)
    deallocate(pc%liqThermalConductivity%D)
    deallocate(pc%liqThermalConductivity%E)
    deallocate(pc%liqThermalConductivity%range)
    deallocate(pc%liqThermalConductivity%eqn)
    deallocate(pc%liqThermalConductivity%val)
    deallocate(pc%liqThermalConductivity)

    deallocate(pc%vapThermalConductivity%A)
    deallocate(pc%vapThermalConductivity%B)
    deallocate(pc%vapThermalConductivity%C)
    deallocate(pc%vapThermalConductivity%D)
    deallocate(pc%vapThermalConductivity%E)
    deallocate(pc%vapThermalConductivity%range)
    deallocate(pc%vapThermalConductivity%eqn)
    deallocate(pc%vapThermalConductivity%val)
    deallocate(pc%vapThermalConductivity)

    deallocate(pc%SurfaceTension%A)
    deallocate(pc%SurfaceTension%B)
    deallocate(pc%SurfaceTension%C)
    deallocate(pc%SurfaceTension%D)
    deallocate(pc%SurfaceTension%E)
    deallocate(pc%SurfaceTension%range)
    deallocate(pc%SurfaceTension%eqn)
    deallocate(pc%SurfaceTension%val)  
    deallocate(pc%SurfaceTension)

    deallocate(pc%DiffusionCoefficientWilkeLee)
    deallocate(pc%DiffusionCoefficientFuller)
    deallocate(pc%DiffusionCoefficientTheory)
    deallocate(pc%DiffusionCoefficientHighPres)

    deallocate(pc)

  end subroutine deallocate_pc

end module pc_defs
