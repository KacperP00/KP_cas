! Pure Chemicals Database adapted from matlab code of J. Palmer (WSA)
module pc_func
  use pc_defs
  use pc_database

  public :: choosePC, computeSolDensity, computeLiqDensity, computeIG_vapDensity, &
            computeVapPressure, computeHeatOfVap, computeSolHeatCapacity, computeLiqHeatCapacity, &
            computeIG_HeatCapacity, computeSeconfViralCoef, computeLiqViscosity, &
            computeVapViscosity, computeLiqThermalConductivity, computeSurfaceTension, &
            computeDiffusionCoeffcientWilkeLee, computeDiffusionCoeffcientFuller

contains

  ! Select the fuel from database
  subroutine choosePC(pc,fuel)
    implicit none

    ! ---------------------------------------------
    type(pc_t), pointer, intent(inout) :: pc
    character(len=128), intent(in) :: fuel
    ! ---------------------------------------------

    select case(trim(fuel))
    case('diesel')

       call pc_diesel(pc)
       write(*,*) 'Warning! Not all physical properties of diesel are available!'
       write(*,*) 'Some properties are taken from ndodecane. Please check pc_database.f90 for details!'

    case('ndodecane','n-dodecane','duodecane','adakane 12')

       call pc_ndodecane(pc)

    case('1butanol','1-butanol','1-Butanol','n-butanol','n-Butanol, Butan-1-ol')

       call pc_1butanol(pc)

    case('1decanol','1-decanol','1-Decanol','n-decanol','n-Decanol', 'decyl alcohol')

       call pc_1decanol(pc)

    case('1dodecanol','1-dodecanol','1-Dodecanol','n-dodecanol','Lauryl 24')

       call pc_1dodecanol(pc)

    case('1heptanol','1-heptanol','1-Heptanol','n-heptanol','n-Heptanol', 'heptyl alcohol')

       call pc_1heptanol(pc)

    case('1hexanol','1-hexanol','1-Hexanol','1Hexanol','n-Hexanol','n-hexanol')

       call pc_1hexanol(pc)

    case('1nonanol','1-nonanol','1Nonanol','1-Nonanol','Alcohol C-9','n-Nonyl Alcohol')

       call pc_1nonanol(pc)

    case('1octanol','1-octanol','1-Octanol','n-octanol','n-Octanol, octanol')

       call pc_1octanol(pc)

    case('1pentanal','1-pentanal','1-Pentanal','1Pentanal','n-Pentanal','n-pentanal')

       call pc_1pentanal(pc)

    case('1pentanol','1-pentanol','1-Pentanol','1Pentanol','n-Pentanol','n-pentanol')

       call pc_1pentanol(pc)

    case('1propanol','1-propanol','1Ppropanol','n-propanol','n-Propanol', 'propanol-1')

       call pc_1propanol(pc)

    case('1tetradecanol','n-tetradecanol','n-Tetradecanol','1Tetradecanol','1-tetradecanol','1-Tetradecanol')

       call pc_1tetradecanol(pc)

    case('1tridecanol','n-tridecanol','n-Tridecanol','1-Tridecanol')

       call pc_1tridecanol(pc)

    case('1undecanol','1-undecanol','1-Undecanol','1Undecanol','n-Undecanol','n-undecanol')

       call pc_1undecanol(pc)

    case('2butanol','2-butanol','2-Butanol','2Butanol','Butane')

       call pc_2butanol(pc)

    case('2butanone','2-butanone','2-Butanone','butanone')

       call pc_2butanone(pc)

    case('2heptanone','2-heptanone','2Heptanone','2-Heptanone')

       call pc_2heptanone(pc)

    case('2hexanol','2-hexanol','2-Hexanol','2Hexanol')

       call pc_2hexanol(pc)

    case('2hexanone','2-hexanone','2-Hexanone','2Hexanone')

       call pc_2hexanone(pc)

    case('2methyl1butanol','2methyl-1butanol','amylalcohol')

       call pc_2methyl1butanol(pc)

    case('2methyl1propanol','2methyl-1propanol','Fermation Butyl Alcohol','Isobutyl Alcohol')

       call pc_2methyl1propanol(pc)

    case('2methyl2propanol','2methyl-2propanol','t-Butanol')

       call pc_2methyl2propanol(pc)

    case('2methylnonane','2-methylnonane')

       call pc_2methylnonane(pc)

    case('2methyloctane','2-methyloctane')

       call pc_2methyloctane(pc)

    case('2nonanone','2-nonanone','2-Nonanone')

       call pc_2nonanone(pc)

    case('2octanone','2-octanone','methyl hexyl ketone')

       call pc_2octanone(pc)

    case('2pentanol','2-pentanol','pentanol-2')

       call pc_2pentanol(pc)

    case('2pentanone','2-pentanone','Ethylacetone')

       call pc_2pentanone(pc)

    case('2propanol','2-propanol','dimethyl carbinol','isohol')

       call pc_2propanol(pc)

    case('3heptanone','3-heptanone','butyl ethyl ketone','heptane-3-one')

       call pc_3heptanone(pc)

    case('3hexanone','3-hexanone','ethyl-n-propylketone')

       call pc_3hexanone(pc)

    case('3pentanone','3-pentanone','metacetone','ethyl ketone')

       call pc_3pentanone(pc)

    case('3petanol','3-petanol','pentanol-3','diethyl carbinol')

       call pc_3petanol(pc)

    case('acetone','Acetone','Allylic Alcohol','Dimethylketal','3-Hydroxypropene','Propanone','beta-Ketopropane')

       call pc_acetone(pc)

    case('diethylether','aether','ether','ethyl ether')

       call pc_diethylether(pc)

    case('OME1','Methylal','ome1','methylal')

       call pc_dimethoxymethane(pc)

    case('dimethylether','methyl ether')

       call pc_dimethylether(pc)

    case('dinbutylether','di-n-butyl-ether','butyl oxide','1,1-oxybisbutane')

       call pc_dinbutylether(pc)

    case('dinpropylether','di-n-propyl-ether','n-propyl ether')

       call pc_dinpropylether(pc)

    case('ethane','bimethyl','ethyl hydride','methylmethane')

       call pc_ethane(pc)

    case('ethanol','Ethanol')

       call pc_ethanol(pc)

    case('gammavalerolactone','gamma-pentalactone','4-pentanolide')

       call pc_gammavalerolactone(pc)

    case('isobutane','iso-butane','i-butane')

       call pc_isobutane(pc)

    case('isoheptane','iso-heptane')

       call pc_isoheptane(pc)

    case('isohexane','iso-hexane','2-methylpentane')

       call pc_isohexane(pc)

    case('isooctane','iso-octane')

       call pc_isooctane(pc)

    case('isopentane','iso-pentane','butane','2-methylbutane')

       call pc_isopentane(pc)

    case('methane','Methane')

       call pc_methane(pc)

    case('methanol','Methanol','methyl alcohol')

       call pc_methanol(pc)

    case('methyltertbutylether','methyl tert-buty lether')

       call pc_methyltertbutylether(pc)

    case('nbutane','n-butane','methylethylmethane')

       call pc_nbutane(pc)

    case('nbutylformate','n-butyl formate','butyl methanoate')

       call pc_nbutylformate(pc)

    case('ndecane','n-decane','decyl hydride')

       call pc_ndecane(pc)

    case('nheptane','n-heptane','dipropymethane')

       call pc_nheptane(pc)

    case('nhexane','n-hexane','hexyl hydride')

       call pc_nhexane(pc)

    case('nnonane','n-nonane','nonyl hydride')

       call pc_nnonane(pc)

    case('noctane','n-octane','octyl hydride')

       call pc_noctane(pc)

    case('npentane','n-pentane','pental hydride')

       call pc_npentane(pc)

    case('ntetradecane','n-tetradecane','tetradecyl hydride')

       call pc_ntetradecane(pc)

    case('ntridecane','n-tridecane','tridecyl hydride')

       call pc_ntridecane(pc)

    case('propane','dimethylmethane','propyl hydride')

       call pc_propane(pc)

    case('tetrahydrofuran','hydrofuran','tetramethylene oxide')

       call pc_tetrahydrofuran(pc)

    case('tetrahydrofurfurylalcohol','tetrahydrofuryl carbinol','oxolan-2-methanol')

       call pc_tetrahydrofurfurylalcohol(pc)

    case default

       write(*,*) 'Error:: Fuel not available in the database'

    end select

  end subroutine choosePC

! Solid Density - kg/m^3 (kmol/m^3)
subroutine computeSolDensity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  pc%T = pc%solDensityMol%range(1,1)

  pc%solDensityMol%val = computeVal(pc,pc%solDensityMol%eqn,pc%T,pc%solDensityMol%A,pc%solDensityMol%B,pc%solDensityMol%C,pc%solDensityMol%D,pc%solDensityMol%E)
  pc%solDensity = pc%solDensityMol%val*pc%MolecularWeight
  
end subroutine computeSolDensity

! Liquid Density - kg/m^3 (kmol/m^3)
subroutine computeLiqDensity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%liqDensityMol%range(1,1) .and. pc%T < pc%liqDensityMol%range(2,1) ) then
     pc%liqDensityMol%val = computeVal(pc,pc%liqDensityMol%eqn,pc%T,pc%liqDensityMol%A,pc%liqDensityMol%B,pc%liqDensityMol%C,pc%liqDensityMol%D,pc%liqDensityMol%E)
  else if (pc%T < pc%liqDensityMol%range(1,1)) then
     !Write(*,*) 'Warning! liqDensity out of range. Setting it for minimum of the range...'
     !pc%T = pc%liqDensityMol%range(1,1)
     pc%liqDensityMol%val = pc%liqDensityMol%range(1,2)
  else if (pc%T > pc%liqDensityMol%range(2,1)) then
     !Write(*,*) 'Warning! liqDensity out of range. Setting it for maximum of the range...'
     !pc%T = pc%liqDensityMol%range(2,1)
     pc%liqDensityMol%val = pc%liqDensityMol%range(2,2)
  end if

  pc%liqDensity = pc%liqDensityMol%val*pc%MolecularWeight

end subroutine computeLiqDensity

! Vapor Density - kg/m^3 (kmol/m^3) Ideal Gas Law
subroutine computeIG_vapDensity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------
  real(WP) :: R_gas

  R_gas = 8.3144598E03_WP   ! J/K/kmol

  pc%IG_vapDensity = pc%p*pc%MolecularWeight/(R_gas*pc%T)

end subroutine computeIG_vapDensity

! Vapor Pressure - Pa
subroutine computeVapPressure(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%vapPressure%range(1,1) .and. pc%T < pc%vapPressure%range(2,1) ) then
     pc%vapPressure%val = computeVal(pc,pc%vapPressure%eqn,pc%T,pc%vapPressure%A,pc%vapPressure%B,pc%vapPressure%C,pc%vapPressure%D,pc%vapPressure%E)
  else if (pc%T < pc%vapPressure%range(1,1)) then
     !Write(*,*) 'Warning! vapPressure out of range. Setting it for minimum of the range...'
     !pc%T = pc%vapPressure%range(1,1)
     pc%vapPressure%val = pc%vapPressure%range(1,2)
  else if (pc%T > pc%vapPressure%range(2,1)) then
     !Write(*,*) 'Warning! vapPressure out of range. Setting it for maximum of the range...'
     !pc%T = pc%vapPressure%range(2,1)
     pc%vapPressure%val = pc%vapPressure%range(2,2)
  end if

end subroutine computeVapPressure

! Heat of Vaporization - J/kg (J/kmol)
subroutine computeHeatOfVap(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%HeatOfVapMol%range(1,1) .and. pc%T < pc%HeatOfVapMol%range(2,1) ) then
     pc%HeatOfVapMol%val = computeVal(pc,pc%HeatOfVapMol%eqn,pc%T,pc%HeatOfVapMol%A,pc%HeatOfVapMol%B,pc%HeatOfVapMol%C,pc%HeatOfVapMol%D,pc%HeatOfVapMol%E)
  else if (pc%T < pc%HeatOfVapMol%range(1,1)) then
     !Write(*,*) 'Warning! HeatOfVap out of range. Setting it for minimum of the range...'
     !pc%T = pc%HeatOfVapMol%range(1,1)
     pc%HeatOfVapMol%val = pc%HeatOfVapMol%range(1,2)
  else if (pc%T > pc%HeatOfVapMol%range(2,1)) then
     !Write(*,*) 'Warning! HeatOfVap out of range. Setting it for maximum of the range...'
     !pc%T = pc%HeatOfVapMol%range(2,1)
     pc%HeatOfVapMol%val = pc%HeatOfVapMol%range(2,2)
  end if

  pc%HeatOfVap = pc%HeatOfVapMol%val/pc%MolecularWeight

end subroutine computeHeatOfVap

! Solid Heat Capacity - J/kg K (J/kmol K)
subroutine computeSolHeatCapacity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%solHeatCapacityMol%range(1,1) .and. pc%T < pc%solHeatCapacityMol%range(2,1) ) then
     pc%solHeatCapacityMol%val = computeVal(pc,pc%solHeatCapacityMol%eqn,pc%T,pc%solHeatCapacityMol%A,pc%solHeatCapacityMol%B,pc%solHeatCapacityMol%C,pc%solHeatCapacityMol%D,pc%solHeatCapacityMol%E)
  else if (pc%T < pc%solHeatCapacityMol%range(1,1)) then
     !Write(*,*) 'Warning! HeatCapacity out of range. Setting it for minimum of the range...'
     !pc%T = pc%solHeatCapacityMol%range(1,1)
     pc%solHeatCapacityMol%val = pc%solHeatCapacityMol%range(1,2)
  else if (pc%T > pc%solHeatCapacityMol%range(2,1)) then
     !Write(*,*) 'Warning! HeatCapacity out of range. Setting it for maximum of the range...'
     !pc%T = pc%solHeatCapacityMol%range(2,1)
     pc%solHeatCapacityMol%val = pc%solHeatCapacityMol%range(2,2)
  end if

  pc%solHeatCapacity = pc%solHeatCapacityMol%val*pc%MolecularWeight

end subroutine computeSolHeatCapacity

! Liquid Heat Capacity - J/kg K (J/kmol K)
subroutine computeLiqHeatCapacity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%liqHeatCapacityMol%range(1,1) .and. pc%T < pc%liqHeatCapacityMol%range(2,1) ) then
     pc%liqHeatCapacityMol%val = computeVal(pc,pc%liqHeatCapacityMol%eqn,pc%T,pc%liqHeatCapacityMol%A,pc%liqHeatCapacityMol%B,pc%liqHeatCapacityMol%C,pc%liqHeatCapacityMol%D,pc%liqHeatCapacityMol%E)
  else if (pc%T < pc%liqHeatCapacityMol%range(1,1)) then
     !Write(*,*) 'Warning! liqHeatCapacity out of range. Setting it for minimum of the range...'
     !pc%T = pc%liqHeatCapacityMol%range(1,1)
     pc%liqHeatCapacityMol%val = pc%liqHeatCapacityMol%range(1,2)
  else if (pc%T > pc%liqHeatCapacityMol%range(2,1)) then
     !Write(*,*) 'Warning! liqHeatCapacity out of range. Setting it for maximum of the range...'
     !pc%T = pc%liqHeatCapacityMol%range(2,1)
     pc%liqHeatCapacityMol%val = pc%liqHeatCapacityMol%range(2,2)
  end if

  pc%liqHeatCapacity = pc%liqHeatCapacityMol%val/pc%MolecularWeight

end subroutine computeLiqHeatCapacity

! Ideal Gas heat Capacity - J/kg K (J/kmol K)
subroutine computeIG_HeatCapacity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%IG_HeatCapacityMol%range(1,1) .and. pc%T < pc%IG_HeatCapacityMol%range(2,1) ) then
     pc%IG_HeatCapacityMol%val = computeVal(pc,pc%IG_HeatCapacityMol%eqn,pc%T,pc%IG_HeatCapacityMol%A,pc%IG_HeatCapacityMol%B,pc%IG_HeatCapacityMol%C,pc%IG_HeatCapacityMol%D,pc%IG_HeatCapacityMol%E)
  else if (pc%T < pc%IG_HeatCapacityMol%range(1,1)) then
     !Write(*,*) 'Warning! IG_HeatCapacity out of range. Setting it for minimum of the range...'
     !pc%T = pc%IG_HeatCapacityMol%range(1,1)
     pc%IG_HeatCapacityMol%val = pc%IG_HeatCapacityMol%range(1,2)
  else if (pc%T > pc%IG_HeatCapacityMol%range(2,1)) then
     !Write(*,*) 'Warning! IG_HeatCapacity out of range. Setting it for maximum of the range...'
     !pc%T = pc%IG_HeatCapacityMol%range(2,1)
     pc%IG_HeatCapacityMol%val = pc%IG_HeatCapacityMol%range(2,2)
  end if

  pc%IG_HeatCapacity = pc%IG_HeatCapacityMol%val/pc%MolecularWeight

end subroutine computeIG_HeatCapacity

! Second Viral Coefficient - m^3/kg (m^3/kmol)
subroutine computeSeconfViralCoef(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%SecondViralCoefMol%range(1,1) .and. pc%T < pc%SecondViralCoefMol%range(2,1) ) then
     pc%SecondViralCoefMol%val = computeVal(pc,pc%SecondViralCoefMol%eqn,pc%T,pc%SecondViralCoefMol%A,pc%SecondViralCoefMol%B,pc%SecondViralCoefMol%C,pc%SecondViralCoefMol%D,pc%SecondViralCoefMol%E)
  else if (pc%T < pc%SecondViralCoefMol%range(1,1)) then
     !Write(*,*) 'Warning! SecondViralCoef out of range. Setting it for minimum of the range...'
     !pc%T = pc%SecondViralCoefMol%range(1,1)
     pc%SecondViralCoefMol%val = pc%SecondViralCoefMol%range(1,2)
  else if (pc%T > pc%SecondViralCoefMol%range(2,1)) then
     !Write(*,*) 'Warning! SecondViralCoef out of range. Setting it for maximum of the range...'
     !pc%T = pc%SecondViralCoefMol%range(2,1)
     pc%SecondViralCoefMol%val = pc%SecondViralCoefMol%range(2,2)
  end if

  pc%SecondViralCoef = pc%SecondViralCoefMol%val/pc%MolecularWeight

end subroutine computeSeconfViralCoef

! Liquid Viscosity - Pa s
subroutine computeLiqViscosity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%liqViscosity%range(1,1) .and. pc%T < pc%liqViscosity%range(2,1) ) then
     pc%liqViscosity%val = computeVal(pc,pc%liqViscosity%eqn,pc%T,pc%liqViscosity%A,pc%liqViscosity%B,pc%liqViscosity%C,pc%liqViscosity%D,pc%liqViscosity%E)
  else if (pc%T < pc%liqViscosity%range(1,1)) then
     !Write(*,*) 'Warning! liqViscosity out of range. Setting it for minimum of the range...'
     !pc%T = pc%liqViscosity%range(1,1)
     pc%liqViscosity%val = pc%liqViscosity%range(1,2)
  else if (pc%T > pc%liqViscosity%range(2,1)) then
     !Write(*,*) 'Warning! liqViscosity out of range. Setting it for maximum of the range...'
     !pc%T = pc%liqViscosity%range(2,1)
     pc%liqViscosity%val = pc%liqViscosity%range(2,2)
  end if

end subroutine computeLiqViscosity

! Vapor Viscosity - Pa s
subroutine computeVapViscosity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%vapViscosity%range(1,1) .and. pc%T < pc%vapViscosity%range(2,1) ) then
     pc%vapViscosity%val = computeVal(pc,pc%vapViscosity%eqn,pc%T,pc%vapViscosity%A,pc%vapViscosity%B,pc%vapViscosity%C,pc%vapViscosity%D,pc%vapViscosity%E)
  else if (pc%T < pc%vapViscosity%range(1,1)) then
     !Write(*,*) 'Warning! vapViscosity out of range. Setting it for minimum of the range...'
     !pc%T = pc%vapViscosity%range(1,1)
     pc%vapViscosity%val = pc%vapViscosity%range(1,2)
  else if (pc%T > pc%vapViscosity%range(2,1)) then
     !Write(*,*) 'Warning! vapViscosity out of range. Setting it for maximum of the range...'
     !pc%T = pc%vapViscosity%range(2,1)
     pc%vapViscosity%val = pc%vapViscosity%range(2,2)
  end if

end subroutine computeVapViscosity

! Liquid Thermal Conductivity - W/mK
subroutine computeLiqThermalConductivity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%liqThermalConductivity%range(1,1) .and. pc%T < pc%liqThermalConductivity%range(2,1) ) then
     pc%liqThermalConductivity%val = computeVal(pc,pc%liqThermalConductivity%eqn,pc%T,pc%liqThermalConductivity%A,pc%liqThermalConductivity%B,pc%liqThermalConductivity%C,pc%liqThermalConductivity%D,pc%liqThermalConductivity%E)
  else if (pc%T < pc%liqThermalConductivity%range(1,1)) then
     !Write(*,*) 'Warning! liqThermalConductivity out of range. Setting it for minimum of the range...'
     !pc%T = pc%liqThermalConductivity%range(1,1)
     pc%liqThermalConductivity%val = pc%liqThermalConductivity%range(1,2)
  else if (pc%T > pc%liqThermalConductivity%range(2,1)) then
     !Write(*,*) 'Warning! liqThermalConductivity out of range. Setting it for maximum of the range...'
     !pc%T = pc%liqThermalConductivity%range(2,1)
     pc%liqThermalConductivity%val = pc%liqThermalConductivity%range(2,2)
  end if

end subroutine computeLiqThermalConductivity

! Vapor Thermal Conductivity - W/mK
subroutine computeVapThermalConductivity(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%vapThermalConductivity%range(1,1) .and. pc%T < pc%vapThermalConductivity%range(2,1) ) then
     pc%vapThermalConductivity%val = computeVal(pc,pc%vapThermalConductivity%eqn,pc%T,pc%vapThermalConductivity%A,pc%vapThermalConductivity%B,pc%vapThermalConductivity%C,pc%vapThermalConductivity%D,pc%vapThermalConductivity%E)
  else if (pc%T < pc%vapThermalConductivity%range(1,1)) then
     !Write(*,*) 'Warning! vapThermalConductivity out of range. Setting it for minimum of the range...'
     !pc%T = pc%vapThermalConductivity%range(1,1)
     pc%vapThermalConductivity%val = pc%vapThermalConductivity%range(1,2)
  else if (pc%T > pc%vapThermalConductivity%range(2,1)) then
     !Write(*,*) 'Warning! vapThermalConductivity out of range. Setting it for maximum of the range...'
     !pc%T = pc%vapThermalConductivity%range(2,1)
     pc%vapThermalConductivity%val = pc%vapThermalConductivity%range(2,2)
  end if

end subroutine computeVapThermalConductivity

! Surface Tension - N/m
subroutine computeSurfaceTension(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------

  if (pc%T > pc%SurfaceTension%range(1,1) .and. pc%T < pc%SurfaceTension%range(2,1) ) then
     pc%SurfaceTension%val = computeVal(pc,pc%SurfaceTension%eqn,pc%T,pc%SurfaceTension%A,pc%SurfaceTension%B,pc%SurfaceTension%C,pc%SurfaceTension%D,pc%SurfaceTension%E)
  else if (pc%T < pc%SurfaceTension%range(1,1)) then
     !Write(*,*) 'Warning! SurfaceTension out of range. Setting it for minimum of the range...'
     !pc%T = pc%SurfaceTension%range(1,1)
     pc%SurfaceTension%val = pc%SurfaceTension%range(1,2)
  else if (pc%T > pc%SurfaceTension%range(2,1)) then
     !Write(*,*) 'Warning! SurfaceTension out of range. Setting it for maximum of the range...'
     !pc%T = pc%SurfaceTension%range(2,1)
     pc%SurfaceTension%val = pc%SurfaceTension%range(2,2)
  end if

end subroutine computeSurfaceTension

! Diffusion Coefficients - m^2/s
! According to J. Poling, The Properties of Gases and Liquids, 5th edition

subroutine computeDiffusionCoeffcientWilkeLee(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------
  real(WP) :: M_fuel, M_air, M_AB, V_b_fuel, sigma_air, sigma_fuel, sigma_AB, &
              epsilon_air, epsilon_fuel, epsilon_AB, T_star, A, B, C, D, E, F,&
              G, H, Omega_D

  M_fuel = pc%MolecularWeight;
  M_air = 28.96_WP;
  M_AB = 2.0_WP* (((1/M_air) + (1/M_fuel))**(-1));
  
  ! V_b = liquid molar volume @ normal boiling point 
  ! tabled values according to Schroeder (Chaper 4):
  V_b_fuel =  7.0_WP * ( pc%ChemicalFormula%C + &
                         pc%ChemicalFormula%H + &
                         pc%ChemicalFormula%O )

  sigma_air = 3.617_WP ! air /Angström
  sigma_fuel = 1.18_WP*V_b_fuel**(1/3)  ! fuel /Angström
  sigma_AB = 0.5_WP* (sigma_air + sigma_fuel)
  
  epsilon_air = 97.0_WP ! air /K
  epsilon_fuel = 1.15_WP*pc%NormalBoilingPoint; ! normal boiling point /K
  epsilon_AB = (epsilon_air * epsilon_fuel)**(0.5_WP); ! /K

  T_star = pc%T / epsilon_AB

  A = 1.06036_WP; B = 0.15610_WP; C = 0.19300_WP; D = 0.47635_WP;
  E = 1.03587_WP; F = 1.52996_WP; G = 1.76474_WP; H = 3.89411_WP;

  Omega_D = (A/(T_star**B)) + (C/(exp(D * T_star))) + &
            (E/(exp(F * T_star))) + (G/(exp(H * T_star)))

  pc%DiffusionCoefficientWilkeLee = ((3.03_WP - 0.98_WP/(M_AB**(0.5_WP)))*(1E-3_WP) &
                                  * pc%T**(1.5_WP))/((pc%p/(100000.0_WP))*(M_AB**(0.5_WP)) &
                                  * (sigma_AB**2)*Omega_D)/(10000.0_WP)

end subroutine computeDiffusionCoeffcientWilkeLee

subroutine computeDiffusionCoeffcientFuller(pc)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  ! ---------------------------------------------
  real(WP) :: M_fuel, M_air, M_AB, Sigma_A, Sigma_B

  M_fuel = pc%MolecularWeight
  M_air = 28.96_WP
  M_AB = 2.0_WP* (((1/M_air) + (1/M_fuel))**(-1))
  Sigma_A = 19.7_WP; ! air, Table 11.1
  Sigma_B = 15.9_WP * pc%ChemicalFormula%C + &
            2.31_WP * pc%ChemicalFormula%H + &
            6.11_WP * pc%ChemicalFormula%O
 
  ! aromatic rings should be considered differently!!! with -18.3; heterocyclinc rings with also -18.3
              
  pc%DiffusionCoefficientFuller = (0.00143_WP*pc%T**1.75_WP) / ((pc%p/(100000.0_WP)) &
                                * (M_AB**0.5_WP) * (Sigma_A**(1.0_WP/3.0_WP) &
                                + Sigma_B**(1.0_WP/3.0_WP))**2)/(10000.0_WP)
 
end subroutine computeDiffusionCoeffcientFuller

real(WP) function computeVal(pc,eqn,T,A,B,C,D,E)
 implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  integer, intent(in) :: eqn
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------

  select case (eqn)
  case (100)
     computeVal = eqn100(pc,T,A,B,C,D,E)
  case (101)
     computeVal = eqn101(pc,T,A,B,C,D,E)
  case (102)
     computeVal = eqn102(pc,T,A,B,C,D,E)
  case (104)
     computeVal = eqn104(pc,T,A,B,C,D,E)
  case (105)
     computeVal = eqn105(pc,T,A,B,C,D,E)
  case (106)
     computeVal = eqn106(pc,T,A,B,C,D,E)
  case (107)
     computeVal = eqn107(pc,T,A,B,C,D,E)
  case (114)
     computeVal = eqn114(pc,T,A,B,C,D,E)
  case (115)
     computeVal = eqn115(pc,T,A,B,C,D,E)
  case default 
     write(*,*) 'No equation available to compute property.'
  end select

end function computeVal

real(WP) function eqn100(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------

  eqn100 = A + B*T + C*T**2 + D*T**3 + E*T**4

end function eqn100

real(WP) function eqn101(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------

  eqn101 = exp(A + B/T + C*log(T) + D*(T**E))

end function eqn101

real(WP) function eqn102(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------

  eqn102 = (A * T**B) / (1 + C/T + D/(T**2))

end function eqn102

real(WP) function eqn104(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------

  eqn104 = A + (B/T) + (C/((T)**3)) + (D/((T)**8)) + (E/((T)**9))

end function eqn104

real(WP) function eqn105(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------

  eqn105 = A/(B**(1+(1-T/C)**D))

end function eqn105

real(WP) function eqn106(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------
  real(WP) :: Tred
  

  Tred = T / pc%Tcrit;
  eqn106 = A * (1 - Tred)**(B + C*Tred + D*Tred**2 + E*Tred**3)

end function eqn106

real(WP) function eqn107(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------
  real(WP) :: n, m

  n = C/T
  m = E/T
  eqn107 = A + B*(n/sinh(n))**2 + D*(m/cosh(m))**2;

end function eqn107

real(WP) function eqn114(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------
  real(WP) :: n, m, Tred, t1

  Tred = T / pc%Tcrit
  t1 = (1.0_WP-Tred)
  eqn114 = (A**2)/t1+B-2.0_WP*A*C*t1-2.0_WP*A*D*(t1**2)-(C**2)*(t1**3)/3.0_WP-C*D*(t1**4)/2.0_WP-(D**2)*(t1**5)/5.0_WP

end function eqn114

real(WP) function eqn115(pc,T,A,B,C,D,E)
  implicit none
  ! --------------------------------------------
  type(pc_t), pointer, intent(in) :: pc
  real(WP), intent(in) :: T, A, B, C, D, E
  ! --------------------------------------------
  real(WP) :: t1

  t1 = T/1000.0_WP

  eqn115 = A + B*t1 + C*t1**2 + D*t1**3 + E*t1**(-2)

end function eqn115

end module pc_func
