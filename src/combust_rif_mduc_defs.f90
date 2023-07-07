!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com or a.saha at itv.rwth-aachen.de           !
!--------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------!
! Adapted from CIAO code at Institute for Combustion Technology, RWTH Aachen     !
! University. Due credit goes to Marco Davidovic and other writers of the CIAO   !
! code                                                                           !
!--------------------------------------------------------------------------------!
module combust_mduc_data
   use precision
   use string
   use spray_defs

  ! ==================================================================================== !
  ! define mduc object in order to call model from different modules
  ! e.g. spark, combustion, ...
  ! ==================================================================================== !


   implicit none

   integer, PARAMETER :: max_mduc_dim = 2
   integer, PARAMETER :: TYPE_MDUC_NO_SOOT = 0, &
        TYPE_MDUC_1D_SOOT = 1, &
        TYPE_MDUC_2D_SOOT = 2
 
   integer, PARAMETER :: TYPE_2DMDUC_NO_TRANSF = 0, &
        TYPE_2DMDUC_DORAN_TRANSF = 1, &
        TYPE_2DMDUC_HASSE_TRANSF = 2
 
   ! flamelet scalar dissipation rate functional form
   integer, PARAMETER :: TYPE_MDUC_CHI_F_LOG = 1, &
        TYPE_MDUC_CHI_F_ERFC = 2, &
        TYPE_MDUC_CHI_F_QUAD = 3
 
 
!!$   type mduc_t
!!$   ! ==================================================================================== !
!!$   !!! this derived datatype has a corresponding MPI datatype                           !!!
!!$   !!! any changes must be consistent with create_mpi_.. !!!
!!$   ! ==================================================================================== !
!!$      SEQUENCE
!!$ 
!!$      integer(8) :: mem          ! memory location
!!$      !integer :: index          ! mainly to distinguish variables
!!$      integer :: nx1,nx2         ! number of grid points in each direction
!!$      integer :: hdim            ! enthalpy dimension in mduc model
!!$      integer :: transf2d        ! coordinate transformation type in 2D
!!$      ! activate(1)/deactivate(0) source terms : 
!!$      ! Nucleation,Condensation,Coagulation,SurfaceGrowth,Oxidation,Fragmentation
!!$      integer, dimension(:), pointer :: sootActive => null()
!!$   ! ==================================================================================== !
!!$   !!! this derived datatype has a corresponding MPI datatype                           !!!
!!$   !!! any changes must be consistent with create_mpi_.. !!!
!!$   ! ==================================================================================== !
!!$   end type mduc_t
 
   type mduc_a_
 
      type(mduc_t), pointer :: d => null()
 
   end type mduc_a_
 
   integer, PARAMETER :: TYPE_GLOBAL_FLAMELET = 1, &
        TYPE_LOCAL_FLAMELET = 2

   integer, parameter, public :: TYPE_PDF_BETA                   = 0, &
                                 TYPE_PDF_SML                    = 1, &
                                 TYPE_PDF_GAUSS                  = 2, &
                                 TYPE_PDF_LOG                    = 3, &
                                 TYPE_PDF_DIRICHLET              = 4, &
                                 TYPE_PDF_JBETA                  = 5, &
                                 TYPE_PDF_INDEPENDENT_BETA       = 6, &
                                 TYPE_PDF_INDEPENDENT_SML        = 7

contains

  subroutine initialize_mduc_struct(mduc_in)
    implicit none
    type(mduc_t), intent(inout) :: mduc_in
    ! ---------------------------------
    ! ---------------------------------

    call set_default_mduc(mduc_in)
    
  end subroutine initialize_mduc_struct

  subroutine initialize_mduc_array(mduc_a,n_mduc)
    implicit none

    type(mduc_a_), dimension(:), pointer, intent(inout) :: mduc_a
    integer, intent(in) :: n_mduc
    ! ---------------------------------
    integer :: i
    ! ---------------------------------
    
    if (associated(mduc_a)) &
         call destroy_mduc_array(mduc_a)

    allocate(mduc_a(n_mduc))

    do i=1,n_mduc
       allocate(mduc_a(i)%d)
       call set_default_mduc(mduc_a(i)%d)
    end do

  end subroutine initialize_mduc_array

  subroutine set_default_mduc(mduc_in)
    implicit none
    type(mduc_t), intent(inout) :: mduc_in
    ! ---------------------------------
    ! ---------------------------------

    mduc_in%mem = 0
    mduc_in%nx1 = 0
    mduc_in%nx2 = 0
    mduc_in%hdim = -1
    mduc_in%transf2d = TYPE_2DMDUC_DORAN_TRANSF
    nullify(mduc_in%sootActive)

  end subroutine set_default_mduc

  subroutine destroy_mduc(mduc_in)
    implicit none
    type(mduc_t), intent(inout), pointer :: mduc_in
    ! ---------------------------------
    ! ---------------------------------

    if (associated(mduc_in%sootActive)) &
         deallocate(mduc_in%sootActive)

    deallocate(mduc_in)
    nullify(mduc_in)

  end subroutine destroy_mduc


  subroutine destroy_mduc_array(mduc_a)
    implicit none
    
    type(mduc_a_), dimension(:), pointer, intent(inout) :: mduc_a
    ! ---------------------------------
    integer :: i
    ! ---------------------------------
    
    do i=1,ubound(mduc_a(:),1)
       if (associated(mduc_a(i)%d)) &
            call destroy_mduc(mduc_a(i)%d)
    end do

    deallocate(mduc_a)
    nullify(mduc_a)

  end subroutine destroy_mduc_array

end module combust_mduc_data
 
module combust_rif_defs
  use precision
  use combust_mduc_data

  implicit none

  !!!-------------------------------------------------!!!
  !!!--- Definition of RIF model related variables ---!!!
  !!!-------------------------------------------------!!!

  ! ---------------------------------!
  !     Parameter definitions        !
  ! ---------------------------------!
  ! rif solution state (needs to be consistent w/ rif_check_injection)
  integer, PARAMETER :: TYPE_RIF_I_1D = 10, &
       TYPE_RIF_R_1D  = 11, &
       TYPE_RIF_II_2D = 20, &
       TYPE_RIF_RI_2D = 21, &
       TYPE_RIF_RR_2D = 22

  ! rif scalar dissipation rate calculation method
  integer, PARAMETER :: TYPE_RIF_CHI_DISTR = 1, &
       TYPE_RIF_CHI_BIN = 2, &
       TYPE_RIF_CHI_PITSCH = 3

  ! rif scalar dissipation rate scaling method
  integer, PARAMETER :: TYPE_RIF_NO_CHI_NORM = 0, &
       TYPE_RIF_G_CHI_NORM = 1, &
       TYPE_RIF_L_CHI_NORM = 2, &
       TYPE_RIF_GZ_CHI_NORM = 3, &
       TYPE_RIF_LZ_CHI_NORM = 4

  ! number of species in hmom soot surface and oxidation reactions
  integer, PARAMETER :: nSpecSoot = 7, &
                        nProcSoot = 6 ! physical/chemical processes

  ! soot source term calculation method
  integer, PARAMETER :: TYPE_NO_SOOT_SRC = 0, &
                        TYPE_LOCAL_SOOT_SRC = 1, &
                        TYPE_CONDZ_SOOT_SRC = 2

  ! ---------------------------------!
  !            Data types            !
  ! ---------------------------------!
  type rif_flow_t
     character(len=str_medium) :: sc_name
     logical :: update_bounds
     integer :: type ! flamelet representative for local cell or entire domain
     real(WP), dimension(:,:,:), pointer :: SC
     real(WP), dimension(:,:,:), pointer :: SC_VAR
     real(WP), dimension(:,:,:), pointer :: SC_CHI
     real(WP) :: min_val,max_val ! flow scalar bounds
     ! inverse of max_val - min_val
     real(WP) :: deltaI
     ! monitor variables
     real(WP) :: chi_min_all,chi_max_all
  end type rif_flow_t

  type pdf_t
     integer :: nmean
     integer :: nvar
     real(WP) :: meanBounds(2)
     real(WP) :: varBounds(2)
     real(WP),dimension(:),pointer :: mean_grid
     real(WP),dimension(:),pointer :: dmeanI
     real(WP),dimension(:),pointer :: var_grid
     real(WP),dimension(:),pointer :: dvarI
     real(WP),dimension(:,:,:), pointer :: pdf
     integer(WP),dimension(:,:,:), pointer :: pdfB
  end type pdf_t

  type rif_t
     ! flamelet number
     integer :: i
     ! is this flamelet currently solved?
     logical :: active
     ! mduc data
     type(mduc_t) :: mduc
     ! flow data
     type(rif_flow_t), dimension(:),allocatable :: f
     ! major species fields
     real(WP), dimension(:), pointer :: Yi_SC => null() 
     ! flamelet marker (tracking scalar for MRIF)
     real(WP), dimension(:), pointer :: mark => null()
     ! flamelet weight (MRIF)
     real(WP), dimension(:), pointer :: weight => null()
     ! mask for locally disabling rif
     integer, dimension(:,:,:), pointer :: mask => null()
     ! progress variable field
     integer, dimension(:,:,:), pointer :: prog => null()
     ! flamelet dimension
     integer :: ndim
     ! number of species in chemical mechanism
     integer :: Nspec
     ! number of non-steady-state species 
     integer :: NspecS 
     ! number of major species
     integer :: NspecT
     ! soot method (local / conditional)
     integer :: soot_method 
     ! number of soot moments
     integer :: nSootEq
     ! soot species names
     character(len=str_medium), dimension(:), pointer :: soot_names => null()
     ! sott species index
     integer, dimension(:), pointer :: soot_ispec => null()
     ! phi threshold
     real(WP) :: phimin_thresh(2) 
     ! weight threshold (MRIF)
     real(WP) :: wmin_thresh
     ! pdf type (recently changed to integer values)
     !character(len=str_medium) :: pdf_type(2)
     integer :: pdf_type(2)
     ! injection status
     integer :: inj
     ! start of chemistry calculation
     real(WP), dimension(:), pointer :: tinit => null()
     ! scalar dissipation rate mode (distribution or global-zero)
     integer, dimension(:), pointer :: chi_mode => null()
     ! scalar dissipation rate functional form
     integer :: chi_form  
     !
     integer, dimension(:), pointer :: chi_norm => null()
     ! calculate reference scalar dissipation to scale functional form?
     logical :: calc_chi_ref
     ! variables for visc computation
     real(WP),dimension(:),pointer :: kOverEps, mucoeff
     real(WP),dimension(:,:),pointer :: phi_W
     ! write frequency in steps
     integer :: dump_freq
     ! root process which solves 1d flamelet
     integer :: root, nproc, irank
     ! flamelet solution arrays
     ! temperature 1d
     real(WP), dimension(:), pointer :: fT_1d => null()
     ! all species 1d
     real(WP), dimension(:,:), pointer :: fY_1d => null()
     ! major species 1d
     real(WP), dimension(:,:), pointer :: fYT_1d => null()
     ! soot source terms 1d
     real(WP), dimension(:,:), pointer :: fS_1d => null()
     ! temperature 2d
     real(WP), dimension(:,:), pointer :: fT_2d => null()
     ! major species 2d
     real(WP), dimension(:,:,:), pointer :: fY_2d => null()
     ! next flamelet (MRIF)
     type(rif_t), pointer :: next => null()
  end type rif_t

  ! global variables
  type(rif_t), pointer, public :: first_rif => null()

  ! pdf table
  logical :: precomp_pdf
  type(pdf_t), pointer :: prePDF

  real(WP)   , dimension(:), pointer :: switch_times
  integer(WP), dimension(:), pointer :: flamelet_number
  integer                            :: n_flamelets
  integer                            :: n_warning_

  real(WP), dimension(:,:), pointer :: fYtmp
  real(WP), dimension(:,:), pointer :: fYOx

  !!!-----------------------------------------------------!!!
  !!!--- End Definition of RIF model related variables ---!!!
  !!!-----------------------------------------------------!!!

  contains

!TODO should this be here?
!  subroutine set_default_rif_t(var)
!
!    implicit none
!
!    type(rif_t), intent(inout) :: var
!    ! ---------------------------------
!
!    ! Pointer
!    nullify( var%mduc, var%f, var%Yi_SC, &
!             var%mask, &
!             var%soot_names, var%soot_ispec, var%tinit, &
!             var%chi_mode, var%chi_norm, var%kOverEps, &
!             var%mucoeff, var%phi_W, var%fT_1d, &
!             var%fY_1d, var%fS_1d, var%fT_2d, &
!             var%fY_2d, var%next )
!     
!    ! initialize
!    call initialize_mduc_struct(var%mduc)
!
!    ! Integer
!    var%ndim = 0
!    var%Nspec = 0
!    var%NspecS = 0
!    var%NspecT = 0
!    var%soot_method = TYPE_NO_SOOT_SRC
!    var%nSootEq = 0
!    var%inj = 0
!    var%chi_form = 0
!    var%dump_freq = 0
!    var%root = 0
!
!    ! Reals
!    var%phimin_thresh(:) = 1.0e-6_WP
!
!    ! Logicals
!    var%calc_chi_ref = .false.
!
!    ! Characters
!    var%pdf_type = 'N/A'
!
!    return
!  end subroutine set_default_rif_t
!
!  subroutine cleanup_rif_t(var)
!
!    implicit none
!
!    type(rif_t), intent(inout) :: var
!    ! ---------------------------------
!
!    !!TODO does there need to be any check if any pointer from %mduc is still
!    ! associated?
!    call deallocate_ptr(var%mduc)
!    call deallocate_ptr(var%f)
!    call deallocate_ptr(var%Yi_SC)
!    call deallocate_ptr(var%rif_mask)
!    call deallocate_ptr(var%soot_names)
!    call deallocate_ptr(var%soot_ispec)
!    call deallocate_ptr(var%tinit)
!    call deallocate_ptr(var%chi_mode)
!    call deallocate_ptr(var%chi_norm)
!    call deallocate_ptr(var%kOverEps)
!    call deallocate_ptr(var%mucoeff)
!    call deallocate_ptr(var%phi_W)
!    call deallocate_ptr(var%fT_1d)
!    call deallocate_ptr(var%fY_1d)
!    call deallocate_ptr(var%fS_1d)
!    call deallocate_ptr(var%fT_2d)
!    call deallocate_ptr(var%fY_2d)
!    call deallocate_ptr(var%next)
!     
!    return
!  end subroutine cleanup_rif_t
!
!  subroutine destroy_rif_t(var)
!
!    implicit none
!
!    type(rif_t), pointer :: var
!    ! ---------------------------------
!
!    if (associated(var%mduc)) call cleanup_mduc(var%mduc)
!    deallocate(var%mduc)
!    if (associated(var)) call cleanup_rif(var)
!    deallocate(var)
!
!    return
!  end subroutine destroy_rif_t
!
!  subroutine set_default_rif_flow_t(var)
!
!    implicit none
!
!    type(rif_flow_t), intent(inout) :: var
!    ! ---------------------------------
!
!    ! Pointer
!    nullify(var%SC, var%SC_VAR, var%SC_CHI)
!
!    ! Integer
!    var%type = 0
!
!    ! Reals
!    var%min_val = 0.0_WP
!    var%max_val  = 0.0_WP
!    var%chi_min_all = 0.0_WP
!    var%chi_max_all = 0.0_WP
!
!    ! Logicals
!    var%update_bounds = .false.
!
!    ! Characters
!    var%sc_name = 'N/A'
!
!   end subroutine set_default_rif_flow_t
!
!   subroutine cleanup_rif_flow_t(var)
! 
!     implicit none
! 
!     type(rif_flow_t), intent(inout) :: var
!     ! ---------------------------------
! 
!     !!TODO does there need to be any check if any pointer from %mduc is still
!     ! associated?
!     call deallocate_ptr(var%SC)
!     call deallocate_ptr(var%SC_VAR)
!     call deallocate_ptr(var%SC_CHI)
!      
!     return
!   end subroutine cleanup_rif_flow_t
! 
!   subroutine destroy_rif_flow_t(var)
! 
!     implicit none
! 
!     type(rif_flow_t), pointer :: var
!     ! ---------------------------------
! 
!     if (associated(var)) call cleanup_rif_flow_t(var)
!     deallocate(var)
! 
!     return
!   end subroutine destroy_rif_flow_t

end module combust_rif_defs
