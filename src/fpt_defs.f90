!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com                                           !
!--------------------------------------------------------------------------------!

! Data Structure for Fuel Property Table
! Written by Abhishek Deshmukh 02/2020
module fpt_defs
  use precision

  implicit none

  ! Definition of fuel property table object
  type fpt_t

     ! Name
     character(len=128), pointer :: name, SMILES

     ! Constants
     type(ChemicalFormula_t), pointer :: ChemicalFormula
     
     real(WP), pointer :: MolecularWeight, Tcrit, pcrit, vcrit
     
     real(WP), pointer :: MeltingPoint, NormalBoilingPoint
     real(WP), pointer :: DiffusionCoefficientWilkeLee, DiffusionCoefficientFuller, DiffusionCoefficientTheory

     real(WP), dimension(:), pointer :: lTemp, vTemp

     real(WP), dimension(:), pointer :: liqViscosity, SurfaceTension, HeatOfVap, vapPressure, &
                          liqThermalConductivity, liqDensity, liqHeatCapacity
     real(WP), dimension(:), pointer :: vapViscosity, vapThermalConductivity, vapDensity, vapHeatCapacity, &
                          DiffusionCoefficient

     character(len=256), pointer :: file
     integer, pointer :: nlines, lrows, vrows
        
  end type fpt_t

  type ChemicalFormula_t
     
     integer, pointer :: C, H, O

  end type ChemicalFormula_t

contains

  ! Allocate fpt with default values
  subroutine nullify_fpt(fpt)
    implicit none

    ! ---------------------------------
    type(fpt_t), pointer, intent(inout) :: fpt
    ! ---------------------------------

    nullify(fpt%name)
    nullify(fpt%SMILES)

    nullify(fpt%ChemicalFormula%C)
    nullify(fpt%ChemicalFormula%H)
    nullify(fpt%ChemicalFormula%O)
    nullify(fpt%ChemicalFormula)

    nullify(fpt%MolecularWeight)
    nullify(fpt%Tcrit)
    nullify(fpt%pcrit)
    nullify(fpt%vcrit)

    nullify(fpt%MeltingPoint)
    nullify(fpt%NormalBoilingPoint)
    nullify(fpt%DiffusionCoefficientWilkeLee)
    nullify(fpt%DiffusionCoefficientFuller)
    nullify(fpt%DiffusionCoefficientTheory)

    nullify(fpt%lTemp,fpt%vTemp)

    nullify(fpt%liqViscosity)
    nullify(fpt%SurfaceTension)
    nullify(fpt%HeatOfVap)
    nullify(fpt%vapPressure)
    nullify(fpt%liqThermalConductivity)
    nullify(fpt%liqDensity)
    nullify(fpt%liqHeatCapacity)

    nullify(fpt%vapViscosity)
    nullify(fpt%vapThermalConductivity)
    nullify(fpt%vapDensity)
    nullify(fpt%vapHeatCapacity)
    nullify(fpt%DiffusionCoefficient)
    
    nullify(fpt%file)
    nullify(fpt%nlines)
    nullify(fpt%lrows)
    nullify(fpt%vrows)

    nullify(fpt)

  end subroutine nullify_fpt

  ! Allocate fpt with default values
  subroutine allocate_fpt(fpt)
    implicit none

    ! ---------------------------------
    type(fpt_t), pointer, intent(inout) :: fpt
    ! ---------------------------------

    !call nullify_fpt(fpt)

    allocate(fpt)
    allocate(fpt%name); fpt%name = 'fpt'
    allocate(fpt%SMILES); fpt%SMILES = 'SMILES'

    allocate(fpt%ChemicalFormula); 
    allocate(fpt%ChemicalFormula%C); fpt%ChemicalFormula%C = -9999.0_WP
    allocate(fpt%ChemicalFormula%H); fpt%ChemicalFormula%H = -9999.0_WP
    allocate(fpt%ChemicalFormula%O); fpt%ChemicalFormula%O = -9999.0_WP

    allocate(fpt%MolecularWeight); fpt%MolecularWeight = -9999.0_WP
    allocate(fpt%Tcrit); fpt%Tcrit = -9999.0_WP
    allocate(fpt%pcrit); fpt%pcrit = -9999.0_WP
    allocate(fpt%vcrit); fpt%vcrit = -9999.0_WP

    allocate(fpt%MeltingPoint); fpt%MeltingPoint = -9999.0_WP
    allocate(fpt%NormalBoilingPoint); fpt%NormalBoilingPoint = -9999.0_WP

!!$    allocate(fpt%DiffusionCoefficientWilkeLee); fpt%DiffusionCoefficientWilkeLee =-9999.0_WP
!!$    allocate(fpt%DiffusionCoefficientFuller); fpt%DiffusionCoefficientFuller =-9999.0_WP
!!$
!!$    allocate(fpt%lTemp); fpt%lTemp = -9999.0_WP
!!$    allocate(fpt%vTemp); fpt%vTemp = -9999.0_WP
!!$
!!$    allocate(fpt%liqViscosity); fpt%liqViscosity = -9999.0_WP
!!$    allocate(fpt%SurfaceTension); fpt%SurfaceTension = -9999.0_WP
!!$    allocate(fpt%HeatOfVap); fpt%HeatOfVap =-9999.0_WP
!!$    allocate(fpt%vapPressure); fpt%vapPressure = -9999.0_WP
!!$    allocate(fpt%liqThermalConductivity); fpt%liqThermalConductivity = -9999.0_WP
!!$    allocate(fpt%liqDensity); fpt%liqDensity = -9999.0_WP
!!$    allocate(fpt%liqHeatCapacity); fpt%liqHeatCapacity =-9999.0_WP
!!$
!!$    allocate(fpt%vapViscosity); fpt%vapViscosity = -9999.0_WP
!!$    allocate(fpt%vapThermalConductivity); fpt%vapThermalConductivity = -9999.0_WP
!!$    allocate(fpt%vapDensity); fpt%vapDensity = -9999.0_WP
!!$    allocate(fpt%vapHeatCapacity); fpt%vapHeatCapacity =-9999.0_WP
!!$    allocate(fpt%DiffusionCoefficient); fpt%DiffusionCoefficient =-9999.0_WP

    allocate(fpt%file); fpt%file = 'noname'
    allocate(fpt%nlines); fpt%nlines = -9999
    allocate(fpt%lrows); fpt%lrows = -9999
    allocate(fpt%vrows); fpt%vrows = -9999

  end subroutine allocate_fpt

  subroutine deallocate_fpt(fpt)
    implicit none
    ! ---------------------------------
    type(fpt_t), pointer, intent(inout) :: fpt
    ! ---------------------------------

    deallocate(fpt%name)
    deallocate(fpt%SMILES)

    deallocate(fpt%ChemicalFormula%C)
    deallocate(fpt%ChemicalFormula%H)
    deallocate(fpt%ChemicalFormula%O)
    deallocate(fpt%ChemicalFormula)

    deallocate(fpt%MolecularWeight)
    deallocate(fpt%Tcrit)
    deallocate(fpt%pcrit)
    deallocate(fpt%vcrit)

    deallocate(fpt%MeltingPoint)
    deallocate(fpt%NormalBoilingPoint)

    deallocate(fpt%DiffusionCoefficientWilkeLee)
    deallocate(fpt%DiffusionCoefficientFuller)
    deallocate(fpt%DiffusionCoefficientTheory)

    deallocate(fpt%lTemp,fpt%vTemp)

    deallocate(fpt%liqViscosity)
    deallocate(fpt%SurfaceTension)
    deallocate(fpt%HeatOfVap)
    deallocate(fpt%vapPressure)
    deallocate(fpt%liqThermalConductivity)
    deallocate(fpt%liqDensity)
    deallocate(fpt%liqHeatCapacity)

    deallocate(fpt%vapViscosity)
    deallocate(fpt%vapThermalConductivity)
    deallocate(fpt%vapDensity)
    deallocate(fpt%vapHeatCapacity)
    deallocate(fpt%DiffusionCoefficient)
    
    deallocate(fpt%file)
    deallocate(fpt%nlines)
    deallocate(fpt%lrows)
    deallocate(fpt%vrows)

    deallocate(fpt)

  end subroutine deallocate_fpt

end module fpt_defs
