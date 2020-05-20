!> @file combustion_mduc.f90

! ==================================================================================== !
!> @brief actual interface functions to mduc code
! ==================================================================================== !
module combust_mduc_interface
  use combust_mduc_data
  implicit none

  private  ! default property

  ! public members
  ! =================================
  public :: mduc_init,mduc_init_flow,mduc_set_data_dir,mduc_set_pressure

  public :: mduc_get_time,mduc_update_time,mduc_get_grid_size,mduc_get_grid_1D,mduc_get_grid_2D, &
       mduc_grid_locate,mduc_get_point_H,mduc_get_point_T,mduc_get_point_Y, mduc_is_parallel

  public :: mduc_get_n_species,mduc_get_nnons_species,mduc_get_species_index, &     
       mduc_get_species_names

  public :: mduc_calc_pdf_1d,mduc_calc_pdf_2d

  public :: mduc_get_var,mduc_calc_mean_var_1d,mduc_calc_mean_var_2d,mduc_calc_Y_mean_1d,mduc_calc_Y_mean_2d

  public :: mduc_update_wmol_all,mduc_comp_thermo_data,mduc_comp_transp_data

  public :: mduc_set_sol_1d,mduc_set_sol_2d

  public :: mduc_advance_flamelet_1d,mduc_solve_2d

  public :: mduc_get_solution_1d,mduc_get_solution_2d

  public :: mduc_calc_chi_form,mduc_calc_chi_log_form,mduc_calc_chi_erfc_form

  public :: mduc_write_flow_data_1d,mduc_write_flow_data_2d,mduc_write_scalars,mduc_write_flamelet, &
       mduc_write_flow_data_mem

  public mduc_soot_init,mduc_get_mean_soot_src,mduc_set_soot_moments,mduc_get_soot_solution,mduc_soot_hmom_src

  public mduc_update_timestep

  public mduc_finalize
  ! ---------------------------------

  ! MDUC Functions
  !------------------------------------------------------------
#ifdef MDUC
  integer, external :: MDUCgetGridSize, MDUCgetNdim, MDUCgetHdim
  integer, external :: MDUCgetNSpecs, MDUCgetNSpecies,MDUCgetSpeciesIndex
  integer, external :: MDUCcalcPDF_1D,MDUCcalcPDF_2D,MDUCcalcMeanY_1D,MDUCcalcMeanY_2D
  integer, external :: MDUCisParallel
  real(WP), external :: MDUCcalcTfromH,MDUCgetTime,MDUCcalcChiLogForm,MDUCcalcChiErfcForm, &
                        MDUCgetPointH,MDUCgetPointT
#endif

contains
  ! ========================================================== !
  !> @brief initialize mduc
  ! ========================================================== !
  subroutine mduc_init(mduc,filename,root, nprocs,com)

    implicit none

    character(len=*), intent(in) :: filename
    type(mduc_t),pointer, intent(inout) :: mduc
    integer, intent(in) :: nprocs,com, root
    real, pointer, dimension(:) :: nullvec => null()
    ! ---------------------------------
    character(len=str_medium) :: cfilename
    integer ::  is_root
    ! ---------------------------------
#ifdef MDUC

    call add_end_char(filename,cfilename)

    ! Temporarily set print processor 
    is_root = 1
    call MDUCtempprintproc(is_root)

    call MDUCinit(cfilename,mduc%mem)
    
    if (mduc%mem.eq.0) then
       write (*,*) 'MDUC not initialized correctly'
       call abort()
    end if

    ! get flamelet grid information
    mduc%nx1 = mduc_get_grid_size(mduc,1)
    mduc%nx2 = mduc_get_grid_size(mduc,2)

    ! get enthalpy direction
    mduc%hdim = MDUCgetHdim(mduc%mem)

    ! Setup parallel infrastructure if needed
    if (MDUCisParallel()>0) then
       call MDUCsetParallelmRIF(com, mduc%mem, nprocs, root)
    end if

#endif
    return
  end subroutine mduc_init

  ! --------------------------------------------- !
  !> @brief define path for flamelet data files   !
  ! --------------------------------------------- !
  subroutine mduc_set_data_dir(mduc,dir)
    implicit none

    type(mduc_t),pointer, intent(inout) :: mduc
    character(len=*), intent(in) :: dir
    ! ---------------------------------
    character(len=str_medium) :: cdir
    ! ---------------------------------

#ifdef MDUC  

    call add_end_char(dir,cdir)

    call MDUCsetDataDir(cdir, mduc%mem)

#endif

    return
  end subroutine mduc_set_data_dir

  ! ========================================================== !
  !> @brief initialize mduc flow quantities
  ! ========================================================== !
  subroutine mduc_init_flow(mduc,start_time,pressure,var1_min,var1_max,var2_min_opt,var2_max_opt)

    implicit none

    type(mduc_t),pointer, intent(in) :: mduc
    real(WP), intent(in) :: var1_min,var1_max,pressure,start_time
    real(WP), intent(in), optional :: var2_min_opt,var2_max_opt
    ! ---------------------------------
    real(WP) :: var2_min, var2_max
    ! ---------------------------------
#ifdef MDUC

    if (PRESENT(var2_min_opt).and.PRESENT(var2_max_opt)) then
       var2_min = var2_min_opt
       var2_max = var2_max_opt
    else
       var2_min = var1_min
       var2_max = var1_max
       if(mduc%nx2.gt.1) then
          write(*,*) 'WARNING: only bounds for 1st flamelet variable specified'
       end if
    end if

    ! initialize flow field data in mduc
    call MDUCinitFlow(start_time, pressure, var1_min, var1_max, var2_min, var2_max, mduc%mem)

#endif
    return
  end subroutine mduc_init_flow

  ! --------------------------------------------- !
  !> @brief reset pressure in flamelet domain     !
  ! --------------------------------------------- !
  subroutine mduc_set_pressure(mduc,rel_step,p_in)
    implicit none

    type(mduc_t),pointer, intent(in) :: mduc
    integer, intent(in) :: rel_step
    real(WP), intent(in) :: p_in
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCsetPressure(rel_step,p_in,mduc%mem)
#endif

    return
  end subroutine mduc_set_pressure

  ! ------------------------------------------------------------ !
  !> @brief reset species and temperature in flamelet domain     !
  ! ------------------------------------------------------------ !
  subroutine mduc_set_sol_1d(mduc,time_in, p_in,fY,fT)
    implicit none

    type(mduc_t),pointer, intent(in) :: mduc
    real(WP), intent(in) :: time_in, p_in
    real(WP), dimension(:,:), intent(in) :: fY
    real(WP), dimension(:), intent(in) :: fT
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCsetsolution1d(time_in, p_in,fY,fT,mduc%mem)
    ! MDUCresetODE is internally called
#endif

    return
  end subroutine mduc_set_sol_1d

  ! ------------------------------------------------------------- !
  !> @brief reset species and temperature in flamelet domain     !
  ! ------------------------------------------------------------- !
  subroutine mduc_set_sol_2d(mduc,time_in, p_in,fY,fT)
    implicit none

    type(mduc_t),pointer, intent(in) :: mduc
    real(WP), intent(in) :: time_in, p_in
    real(WP), dimension(:,:,:), intent(in) :: fY
    real(WP), dimension(:,:), intent(in) :: fT
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCsetsolution2d(time_in, p_in,fY,fT,mduc%mem)
    ! MDUCresetODE is internally called
#endif

    return
  end subroutine mduc_set_sol_2d


  ! ========================================================== !
  !> @brief solver routines
  !                                                            !
  ! ========================================================== !

  ! -------------------------------------------------- !
  !> @brief advance solution in 1d (inert or reactive)
  ! -------------------------------------------------- !
  subroutine mduc_advance_flamelet_1d(mduc,tstep,pavg,var_min,var_max,chi,inert,reset,fY,fT)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), intent(in) :: tstep,pavg,var_min,var_max
    real(WP), dimension(:), intent(in) :: chi
    integer, intent(in) :: inert
    logical, intent(in) :: reset
    real(WP), dimension(:,:), intent(inout) :: fY
    real(WP), dimension(:), intent(in) :: fT
    ! ---------------------------------

#ifdef MDUC
    if (reset) &
         call MDUCresetODE(mduc%mem)

    call MDUCflamelet1D(tstep,pavg,var_min,var_max,chi,inert,fY,fT,mduc%mem)

#endif

    return
  end subroutine mduc_advance_flamelet_1d

  ! -------------------------------------------------- !
  !> @brief advance solution in 1d (inert or reactive)
  ! -------------------------------------------------- !
!!$  subroutine mduc_advance_flamelet_1d(mduc,tstep,pavg,var_min,var_max,chi,inert,reset,fY,fT,fS)
!!$    implicit none
!!$
!!$    type(mduc_t), intent(in) :: mduc
!!$    real(WP), intent(in) :: tstep,pavg,var_min,var_max
!!$    real(WP), dimension(:), intent(in) :: chi
!!$    integer, intent(in) :: inert
!!$    logical, intent(in) :: reset
!!$    real(WP), dimension(:,:), intent(inout) :: fY, fS
!!$    real(WP), dimension(:), intent(in) :: fT
!!$    ! ---------------------------------
!!$
!!$    integer :: nS = 4
!!$    ! ---------------------------------
!!$
!!$#ifdef MDUC
!!$    if (reset) &
!!$         call MDUCresetODE(mduc%mem)
!!$
!!$    call MDUCflamelet1D(tstep,pavg,var_min,var_max,chi,inert,fY,fT,mduc%mem)
!!$
!!$    call MDUCgetsootsolution1D(nS, fS, mduc%mem)
!!$#endif
!!$
!!$    return
!!$  end subroutine mduc_advance_flamelet_1d

  ! ---------------------------------------- !
  !> @brief advance solution in 2d (reactive)
  ! ---------------------------------------- !
  subroutine mduc_solve_2d(mduc,tstep,pavg,var1_min,var1_max,var2_min,var2_max,chi1,chi12,chi2,fY,fT)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), intent(in) :: tstep,pavg,var1_min,var1_max,var2_min,var2_max
    real(WP), dimension(:,:), intent(in) :: chi1,chi12,chi2
    real(WP), dimension(:,:,:), intent(inout) :: fY
    real(WP), dimension(:,:), intent(in) :: fT
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    !call MDUCresetODE(mduc%mem) ! Eric does not reinitialize in 2d..
    call MDUCsolve2D(tstep,pavg,var1_min,var1_max,var2_min,var2_max, &
         chi1,chi12,chi2,fY,fT,mduc%mem)
#endif

    return
  end subroutine mduc_solve_2d

  ! ================================================================== !
  !> @brief some utility interface functions to obtain RIF quantities  !
  !                                                                    !
  ! ================================================================== !
  ! -------------------------------------------------------- !
  !> @brief get time at current(n=0) or prev(n=-1) time step !
  ! -------------------------------------------------------- !
  function mduc_get_time(mduc,nt) result(time_mduc)
    implicit none

    integer, intent(in) :: nt
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------
    real(WP) :: time_mduc
    ! ---------------------------------

#ifdef MDUC
    time_mduc = MDUCgetTime(nt,mduc%mem)
#endif

  end function mduc_get_time
  ! -------------------------------------------------------- !
  !> @brief checks whether MDUC was compiled for MPI usage   !
  ! -------------------------------------------------------- !
  function mduc_is_parallel() result(parallel)
    implicit none

    ! ---------------------------------
    integer :: parallel
    ! ---------------------------------

#ifdef MDUC
    parallel = MDUCisParallel()
#endif

  end function mduc_is_parallel
  ! ------------------------------------- !
  !> @brief update time information       !
  ! ------------------------------------- !
  subroutine mduc_update_time(time,timestep,iteration,mduc)
    implicit none
    
    real(WP), intent(in)     :: time,timestep
    integer, intent(in)      :: iteration
    type(mduc_t), intent(in) :: mduc
    
#ifdef MDUC
    call mducupdatetime(time, timestep, iteration, mduc%mem)
#endif
  
  end subroutine mduc_update_time

  ! ------------------------------------- !
  !> @brief get number of grid points     !
  ! ------------------------------------- !
  function mduc_get_grid_size(mduc,dir) result(npts)
    implicit none

    integer, intent(in) :: dir
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------
    integer :: npts
    ! ---------------------------------

#ifdef MDUC
    npts = MDUCgetGridSize(dir,mduc%mem)
#endif

  end function mduc_get_grid_size

  ! ------------------------------------- !
  !> @brief get flamlet coordinates       !
  ! ------------------------------------- !
  subroutine mduc_get_grid_1D(mduc,phi1)
    implicit none

    real(WP), dimension(:), intent(inout) :: phi1
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCgetGrid(phi1, mduc%mem)
#endif

    return
  end subroutine mduc_get_grid_1D

  subroutine mduc_get_grid_2D(mduc,phi1,phi2)
    implicit none

    real(WP), dimension(:), intent(inout) :: phi1,phi2
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCgetGrid2D(phi1, phi2, mduc%mem)
#endif

    return
  end subroutine mduc_get_grid_2D

  ! ------------------------------------------------------------------ !
  !> @brief get species mass fractions at certain flamelet point       !
  ! ------------------------------------------------------------------ !
  subroutine mduc_get_point_Y(mduc,phi1,phi2,Ypt)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), intent(in) :: phi1,phi2
    real(WP), dimension(:), intent(inout) :: Ypt
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call mducgetpointy(phi1,phi2,Ypt,mduc%mem)
#endif

  end subroutine mduc_get_point_Y

  ! ------------------------------------------------ !
  !> @brief get temperature at certain grid point    !
  ! ------------------------------------------------ !
  function mduc_get_point_T(mduc,phi1,phi2) result(temp_out)
    implicit none

    real(WP), intent(in) :: phi1,phi2
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------
    real(WP) :: temp_out
    ! ---------------------------------

#ifdef MDUC
    temp_out =  MDUCgetPointT(phi1,phi2,mduc%mem)
#endif

  end function mduc_get_point_T

  ! ------------------------------------------------ !
  !> @brief get enthalpy at certain grid point       !
  ! ------------------------------------------------ !
  function mduc_get_point_H(mduc,phi1,phi2) result(enth_out)
    implicit none

    real(WP), intent(in) :: phi1,phi2
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------
    real(WP) :: enth_out
    ! ---------------------------------

#ifdef MDUC
    enth_out =  MDUCgetPointH(phi1,phi2,mduc%mem)
#endif

  end function mduc_get_point_H

  ! --------------------------------------------- !
  !> @brief obtain number of species in mechanism !
  ! --------------------------------------------- !
  function mduc_get_n_species(mduc) result(nspec)
    implicit none

    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------
    integer :: nspec
    ! ---------------------------------

#ifdef MDUC
    nspec = MDUCgetNspecies(mduc%mem)
#endif

  end function mduc_get_n_species

  ! --------------------------------------------- !
  !> @brief obtain number of species in mechanism !
  ! --------------------------------------------- !
  function mduc_get_nnons_species(mduc) result(nspec)
    implicit none

    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------
    integer :: nspec
    ! ---------------------------------

#ifdef MDUC
    nspec = MDUCgetNspecs(mduc%mem)
#endif

  end function mduc_get_nnons_species

  ! --------------------------------------------- !
  !> @brief get mechanism index for species name  !
  ! --------------------------------------------- !
  function mduc_get_species_index(mduc,spec_name) result(indspec)
    implicit none

    type(mduc_t), intent(in) :: mduc
    character(len=*), intent(in) :: spec_name
    ! ---------------------------------
    integer :: indspec
    character(len=str_medium) :: cspec_name
    ! ---------------------------------

    call add_end_char(spec_name, cspec_name)

#ifdef MDUC
    indspec = mducgetspeciesindex(cspec_name,mduc%mem)
#endif

  end function mduc_get_species_index

  ! ------------------------------------------------------------------ !
  !> @brief update molecular weight of all (also steady-state) species !
  ! ------------------------------------------------------------------ !
  subroutine mduc_update_wmol_all(mduc,wmol)
    implicit none

    real(WP), dimension(:), intent(inout) :: wmol
    type(mduc_t), intent(in) :: mduc
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCgetMolWeightAll(wmol,mduc%mem)
#endif

    return
  end subroutine mduc_update_wmol_all


  ! ------------------------------------------------------------------ !
  !> @brief update molecular weight of all (also steady-state) species !
  ! ------------------------------------------------------------------ !
  subroutine mduc_get_species_names(mduc,names,length)
    implicit none

    type(mduc_t), intent(in) :: mduc    
    character(len=str_medium), dimension(:), intent(inout) :: names
    integer, intent(in) :: length
    ! ---------------------------------
    integer :: i
    ! ---------------------------------

#ifdef MDUC
    call MDUCgetSpeciesNames(%ref(names),length,mduc%mem)
#endif

    ! remove the end characters
    do i = 1,size(names)
       call remove_end_char(names(i))
    end do

    return
  end subroutine mduc_get_species_names

  ! ---------------------------------------------------------- !
  !> @brief get flamelet solution vectors for Yi, T
  ! ---------------------------------------------------------- !
  subroutine mduc_get_solution_1d(mduc,fY,fT)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), dimension(:,:), intent(inout) :: fY
    real(WP), dimension(:), intent(inout) :: fT
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    fY = 0.0_WP
    fT = 0.0_WP
    call mducgetsolution1d(fY, fT, mduc%mem)
#endif

    return
  end subroutine mduc_get_solution_1d

  ! ---------------------------------------------------------- !
  !> @brief get flamelet solution arrays for Yi, T
  ! ---------------------------------------------------------- !
  subroutine mduc_get_solution_2d(mduc,fY,fT)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), dimension(:,:,:), intent(inout) :: fY
    real(WP), dimension(:,:), intent(inout) :: fT
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    fY = 0.0_WP
    fT = 0.0_WP
    call mducgetsolution2d(fY, fT, mduc%mem)
#endif

    return
  end subroutine mduc_get_solution_2d

  ! ------------------------------------------------------------------ !
  !> @brief compute species enthalpy and specific heat at given temp   !
  ! ------------------------------------------------------------------ !
  subroutine  mduc_comp_thermo_data(mduc,hi,cpi,tcalc)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), intent(in) :: tcalc
    real(WP), dimension(:), intent(inout) :: hi,cpi
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call mduccomputethermodata(hi,cpi,tcalc,mduc%mem)
#endif

    return
  end subroutine mduc_comp_thermo_data

  ! ------------------------------------------------------------------ !
  !> @brief compute species enthalpy and specific heat at given temp   !
  ! ------------------------------------------------------------------ !
  subroutine  mduc_comp_transp_data(mduc,Nspec,tcalc,Yi,mu,rhoD)
    implicit none

    type(mduc_t), intent(in) :: mduc
    integer, intent(in) :: Nspec
    real(WP), intent(in) :: tcalc
    real(WP), dimension(:), intent(in) :: Yi
    real(WP), intent(in) :: mu,rhoD
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC    
    ! TODO: results of this function need to be verified!!
    !call MDUCcomputeTranspMix(Nspec,tcalc,Yi,mu,rhoD,mduc%mem)
#endif

    return
  end subroutine mduc_comp_transp_data
  ! ---------------------------------------------------------- !
  !> @brief convolute species flamelet solution with pdf  
  !! ,phi1 grid needs to be supplied                           !
  ! ---------------------------------------------------------- !
  subroutine mduc_calc_Y_mean_1d(mduc,NspecS,pdf_method,phi1,phi1_mean,phi1_var,fY,Ymean,pdf_in,pdfB)
    implicit none

    type(mduc_t), intent(in) :: mduc
    integer, intent(in) :: NspecS
    integer, intent(in) :: pdf_method
    real(WP), intent(in) :: phi1_mean,phi1_var
    real(WP), dimension(mduc%nx1), intent(in) :: phi1
    real(WP), dimension(:,:), intent(in) :: fY
    real(WP), dimension(:), intent(inout) :: Ymean
    real(WP), dimension(mduc%nx1), intent(inout), optional :: pdf_in
    integer, intent(inout), optional :: pdfB(2)
    ! ---------------------------------
    integer :: pdf_flag
    ! ---------------------------------

#ifdef MDUC


    Ymean = 0.0_WP
    pdf_flag = MDUCcalcMeanY_1D(pdf_method,mduc%nx1,phi1,NspecS, &
         phi1_mean,phi1_var,fY,Ymean, &
         pdf_in,pdfB)

#endif

    return
  end subroutine mduc_calc_Y_mean_1d

  ! ---------------------------------------------------------- !
  !> @brief convolute species flamelet solution with pdf
  !   ,1 phi1 & phi2 grids need to be supplied                 !
  ! ---------------------------------------------------------- !
  subroutine mduc_calc_Y_mean_2d(mduc,NspecS,pdf_method,phi1,phi2, &
       phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,fY,Ymean)
    implicit none

    type(mduc_t), intent(in) :: mduc
    integer, intent(in) :: NspecS
    integer, intent(in) :: pdf_method
    real(WP), intent(in) :: phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var
    real(WP), dimension(:), intent(in) :: phi1,phi2
    real(WP), dimension(:,:,:), intent(in) :: fY
    real(WP), dimension(:), intent(inout) :: Ymean
    ! ---------------------------------
    integer :: pdf_flag
    real(WP) :: pdf(mduc%nx1,mduc%nx2)
    ! ---------------------------------

#ifdef MDUC


    Ymean = 0.0_WP
    pdf_flag = MDUCcalcMeanY_2D(pdf_method,mduc%nx1,phi1,mduc%nx2,phi2,NspecS, &
         phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,fY,Ymean, &
         pdf,mduc%transf2d)

#endif

    return
  end subroutine mduc_calc_Y_mean_2d


  ! ---------------------------------------------------------- !
  !> @brief get variable values on flamelet grid   
  !! ,supply either 1d or 2d solution array                    !
  ! ---------------------------------------------------------- !
  subroutine mduc_get_var(mduc,varname,var_1d,var_2d)
    implicit none

    character(len=str_medium), intent(in) :: varname
    type(mduc_t), intent(in) :: mduc
    real(WP), dimension(:), intent(inout), optional :: var_1d
    real(WP), dimension(:,:), intent(inout), optional :: var_2d
    ! ---------------------------------
    character(len=str_medium) :: cvarname
    ! ---------------------------------

#ifdef MDUC

    call add_end_char(varname,cvarname)
    if (mduc%nx2.eq.1) then
       var_1d = 0.0_WP
       call MDUCgetVar(cvarname, var_1d, mduc%mem)
    else
       var_2d = 0.0_WP
       call MDUCgetVar(cvarname, var_2d, mduc%mem)
    end if

#endif

    return
  end subroutine mduc_get_var

  ! ---------------------------------------------------------- !
  !> @brief convolute supplied flamelet variable with pdf 
  !! ,phi1 grid and variable vector need to be supplied        !
  ! ---------------------------------------------------------- !
  function mduc_calc_mean_var_1d(nx,var_val,pdf_method,phi1,phi1_mean,phi1_var,pdf) result(VARmean)
    implicit none

    integer,intent(in) :: nx
    integer, intent(in) :: pdf_method
    real(WP), intent(in) :: phi1_mean,phi1_var
    real(WP), dimension(nx), intent(in) :: phi1
    real(WP), dimension(nx), intent(in) :: var_val
    real(WP), dimension(nx), intent(inout) :: pdf
    ! ---------------------------------
    real(WP) :: VARmean
    ! ---------------------------------

#ifdef MDUC

    ! compute the local pdf
    call mduc_calc_pdf_1d(nx,pdf_method,phi1,phi1_mean,phi1_var,pdf)

    ! obtain enthalpy from flamelet
    !call mduc_get_var(mduc,varname,var_1d=var_buf)

    !MDUCcalcLocalMeanVar

    ! convolute local variable with pdf*dx
    VARmean = sum(pdf(:)*var_val(:))

#endif

  end function mduc_calc_mean_var_1d

  ! ---------------------------------------------------------- !
  !> @brief convolute supplied flamelet variable with pdf 
  !! ,phi1 & phi2 grids and var field need to be supplied      !
  ! ---------------------------------------------------------- !
  function mduc_calc_mean_var_2d(mduc,var_val,pdf_method,phi1,phi2, &
       phi1_mean,phi1_var,phi2_mean,phi2_var, &
       phi0_var,pdf) result(VARmean)
    implicit none

    type(mduc_t), intent(in) :: mduc    
    integer, intent(in) :: pdf_method
    real(WP), dimension(mduc%nx1), intent(in) :: phi1
    real(WP), dimension(mduc%nx2), intent(in) :: phi2
    real(WP), intent(in) :: phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var
    real(WP), dimension(mduc%nx1,mduc%nx2), intent(in) :: var_val
    real(WP), dimension(mduc%nx1,mduc%nx2), intent(inout) :: pdf
    ! ---------------------------------
    real(WP) :: VARmean
    integer :: i,j
    ! ---------------------------------

#ifdef MDUC

    ! compute the local pdf
    call mduc_calc_pdf_2d(mduc%nx1,mduc%nx2,pdf_method,phi1,phi2, &
         phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,pdf,mduc%transf2d)

    ! obtain enthalpy from flamelet
    !    call mduc_get_var(mduc,varname,var_2d=var_buf)

    ! convolute local variable with pdf*dx
    VARmean = 0.0_WP
    do j = 1, mduc%nx2
       do i = 1, mduc%nx1
          VARmean = VARmean + var_val(i,j)*pdf(i,j)
       end do
    end do

#endif

  end function mduc_calc_mean_var_2d

  ! ---------------------------------------------------------- !
  !> @brief compute Pdx on flamelet grid  
  !! ,phi1 grid needs to be supplied                           !
  ! ---------------------------------------------------------- !

  subroutine mduc_calc_pdf_1d(nx,pdf_method,phi1,phi1_mean,phi1_var,pdf)
    implicit none

    integer, intent(in) :: nx
    integer :: pdf_method
    real(WP), dimension(nx), intent(in) :: phi1
    real(WP), intent(in) :: phi1_mean,phi1_var
    real(WP), dimension(nx), intent(inout) :: pdf

    ! ---------------------------------
    integer :: pdf_flag,pdf_bound(2)
    ! ---------------------------------

#ifdef MDUC


    pdf_flag = mduccalcpdf_1d(pdf_method,nx,phi1,phi1_mean,phi1_var,pdf,pdf_bound)

#endif

    return
  end subroutine mduc_calc_pdf_1d

  ! ---------------------------------------------------------- !
  !> @brief compute Pdx on flamelet grid 
  !! ,phi1 & phi2 grids need to be supplied                    !
  ! ---------------------------------------------------------- !
  subroutine mduc_calc_pdf_2d(nx1,nx2,pdf_method,phi1,phi2,phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,pdf,transf2d)
    implicit none

!    character(len=str_medium) :: pdf_method
    integer :: pdf_method
    real(WP), dimension(nx1), intent(in) :: phi1
    real(WP), dimension(nx2), intent(in) :: phi2
    real(WP), intent(in) :: phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var
    real(WP), dimension(nx1,nx2), intent(inout) :: pdf
    integer, intent(in) :: nx1,nx2,transf2d
    ! ---------------------------------
    integer :: pdf_flag
    ! ---------------------------------

#ifdef MDUC


    pdf_flag = mduccalcpdf_2d(pdf_method,nx1,phi1,nx2,phi2,phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,pdf,transf2d)

#endif

    return
  end subroutine mduc_calc_pdf_2d
  ! --------------------------------------------- !
  !> @brief convolute soot_src(Z) with local pdf  !
  ! --------------------------------------------- !
  subroutine mduc_soot_init(mduc,nprocess)
    implicit none

    type(mduc_t), intent(inout) :: mduc
    integer, intent(in) :: nprocess
    ! ---------------------------------

    ! ---------------------------------

    allocate(mduc%sootActive(nprocess))

    ! TODO: take user input
    mduc%sootActive = 1
       !param => get_param('RIF Major Species')
       !if (associated(param)) then
       !   call in_param_getsize(rif_in%NspecT,param)
       !else
          
       !end if


  end subroutine mduc_soot_init

  ! --------------------------------------------- !
  !> @brief convolute soot_src(Z) with local pdf  !
  ! --------------------------------------------- !
  subroutine mduc_get_mean_soot_src(mduc,phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,soot_src)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP), intent(in) :: phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var
    real(WP), dimension(:), intent(inout) :: soot_src
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call mducmeansootsrc(phi1_mean,phi1_var,phi2_mean,phi2_var,phi0_var,soot_src,mduc%mem)
#endif

    return
  end subroutine mduc_get_mean_soot_src

  ! ---------------------------------------------------------- !
  !> @brief get flamelet soot source term solution array
  ! ---------------------------------------------------------- !
  subroutine mduc_get_soot_solution(mduc,nS,fS)
    implicit none

    type(mduc_t), intent(in) :: mduc
    integer, intent(in) :: nS
    real(WP), dimension(:,:), intent(inout) :: fS
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    fS = 0.0_WP
    call mducgetsootsolution_(nS,fS,mduc%mem)
#endif

    return
  end subroutine mduc_get_soot_solution
  ! ---------------------------------------------------------- !
  !> @brief compute hmom soot model source terms directly
  ! 
  ! ---------------------------------------------------------- !
  subroutine mduc_soot_hmom_src(mduc,nMom,temp_cv,visc_cv,rho_cv,wmix_cv,c_spec,mom_in,mom_src, &
                                rxnrates,drho_dt,mom_srcP)
    implicit none

    type(mduc_t), intent(in) :: mduc
    integer, intent(in) :: nMom
    ! temperature [K], viscosity [Pa*s], density [kg/m^3], molar mass [kg/m^3]
    real(WP), intent(in) :: temp_cv,visc_cv,rho_cv,wmix_cv
    ! species concentrations in order A4 (pyrene), OH, H2O, H, H2, C2H2, O2 [kmol/m^3]
    real(WP), intent(in) :: c_spec(7)
    ! vector of all moments solved
    real(WP), intent(in) :: mom_in(nMom)
    ! vector of all moments source terms to be computed
    real(WP), intent(inout) :: mom_src(nMom)
    ! output species reaction rates in order A4 (pyrene), OH, H2O, H, H2, C2H2, O2 [kmol/m^3/s]
    real(WP), intent(inout) :: rxnrates(7)
    ! mass change due to solid particle formation/oxidation
    real(WP), intent(inout) :: drho_dt
    ! contribution of each of the physical/chemical processes to every moment source term
    real(WP), dimension(nMom,6), intent(inout) :: mom_srcP
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call mducsoothmom(nMom,temp_cv,visc_cv,rho_cv,wmix_cv,c_spec, mom_in, mom_src, rxnrates, &
                       drho_dt, mom_srcP, mduc%sootActive)
#endif

    return
  end subroutine mduc_soot_hmom_src

  ! ------------------------------------------------ !
  !> @brief pass conditional soot moments to solver  !
  ! ------------------------------------------------ !
  subroutine mduc_set_soot_moments(mduc,nx1,nSootEq,sootMxy)
    implicit none

    type(mduc_t), intent(in) :: mduc
    integer, intent(in) :: nx1,nSootEq
    real(WP), dimension(nx1,nSootEq), intent(in) :: sootMxy
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call mducsetsootmom(nx1,nSootEq,sootMxy,mduc%mem)
#endif

    return
  end subroutine mduc_set_soot_moments

  ! ============================================================= !
  !> @brief procedures for calculation of scalar dissipation rate
  ! ============================================================= !
  function mduc_calc_chi_form(type,phi,phi_ref,phi_max) result(f_phi)
    implicit none

    integer, intent(in) :: type
    real(WP), intent(in) :: phi,phi_ref,phi_max
    ! ---------------------------------
    real(WP) :: f_phi
    ! ---------------------------------

#ifdef MDUC
    select case (type)
    case (TYPE_MDUC_CHI_F_LOG)
       f_phi = MDUCcalcChiLogForm(phi,phi_ref,phi_max)
    case (TYPE_MDUC_CHI_F_ERFC)
       f_phi = MDUCcalcChiErfcForm(phi,phi_ref,phi_max)
    case (TYPE_MDUC_CHI_F_QUAD)
       write(*,*) 'Error: quad form interface not yet implemented'
       call abort()
    end select
#endif

  end function mduc_calc_chi_form


  function mduc_calc_chi_log_form(phi,phi_ref,phi_max) result(f_phi)
    implicit none

    real(WP), intent(in) :: phi,phi_ref,phi_max
    ! ---------------------------------
    real(WP) :: f_phi
    ! ---------------------------------

#ifdef MDUC
    f_phi = MDUCcalcChiLogForm(phi,phi_ref,phi_max)
#endif

  end function mduc_calc_chi_log_form


  function mduc_calc_chi_erfc_form(phi,phi_ref,phi_max) result(f_phi)
    implicit none

    real(WP), intent(in) :: phi,phi_ref,phi_max
    ! ---------------------------------
    real(WP) :: f_phi
    ! ---------------------------------

#ifdef MDUC
    f_phi = MDUCcalcChiErfcForm(phi,phi_ref,phi_max)
#endif

  end function mduc_calc_chi_erfc_form

  ! ========================================================== !
  !> @brief write functions
  ! ========================================================== !
  subroutine  mduc_write_flow_data_1d(mduc,chi,pdfG)
    implicit none

    type(mduc_t), intent(in) :: mduc    
    real(WP), dimension(mduc%nx1), intent(in) :: chi,pdfG
    ! ---------------------------------
    logical, parameter :: trueval = .true., &
         falseval = .false.
    ! ---------------------------------

#ifdef MDUC
    call MDUCwriteFlowData(trueval,chi,falseval,chi,trueval,pdfG, &
         falseval,chi,mduc%mem)
#endif

    return
  end subroutine mduc_write_flow_data_1d

  subroutine  mduc_write_flow_data_2d(mduc,chi1,chi2,pdfG)
    implicit none

    type(mduc_t), intent(in) :: mduc    
    real(WP), dimension(mduc%nx1,mduc%nx2), intent(in) :: chi1,chi2,pdfG
    ! ---------------------------------
    logical, parameter :: trueval = .true., &
         falseval = .false.
    ! ---------------------------------

#ifdef MDUC
    call MDUCwriteFlowData(trueval,chi1,trueval,chi2,trueval,pdfG, &
         falseval,chi1,mduc%mem)
#endif

    return
  end subroutine mduc_write_flow_data_2d

  subroutine  mduc_write_scalars(mduc,t,pdf_1d,pdf_2d)
    implicit none

    type(mduc_t), intent(in) :: mduc
    real(WP),intent(in) :: t
    real(WP), dimension(mduc%nx1), intent(in), optional :: pdf_1d
    real(WP), dimension(mduc%nx1,mduc%nx2), intent(in), optional :: pdf_2d
    ! ---------------------------------
    logical, parameter :: trueval = .true., &
         falseval = .false.
    ! ---------------------------------

#ifdef MDUC
    if (PRESENT(pdf_1d)) then
       call MDUCwriteScalars(t,pdf_1d,mduc%mem)
    else if (PRESENT(pdf_2d)) then
       call MDUCwriteScalars(t,pdf_2d,mduc%mem)
    end if

#endif

    return
  end subroutine mduc_write_scalars

  subroutine mduc_write_flamelet(mduc,filename)
    implicit none

    type(mduc_t), intent(in) :: mduc
    character(len=str_medium), intent(in) :: filename
    ! ---------------------------------
    character(len=str_medium) :: cfilename
    ! ---------------------------------

#ifdef MDUC
    call add_end_char(filename,cfilename)

    call MDUCwriteFlamelet(cfilename,mduc%mem)
#endif

    return
  end subroutine mduc_write_flamelet

  ! ---------------------------------------------------------- !
  !> @brief write flow data that is in mduc memory for restart
  ! ---------------------------------------------------------- !
  subroutine mduc_write_flow_data_mem(mduc,filename)
    implicit none

    type(mduc_t), intent(in) :: mduc
    character(len=str_medium), intent(in) :: filename
    ! ---------------------------------
    character(len=str_medium) :: cfilename
    ! ---------------------------------

#ifdef MDUC
    call add_end_char(filename,cfilename)

    call MDUCwriteFlowDataMem(cfilename,mduc%mem)
#endif

    return
  end subroutine mduc_write_flow_data_mem

  ! ========================================================== !
  !> @brief  update the time step,                             !
  !! (Avoid overwriting data in case of restart)               !
  ! ========================================================== !

  subroutine mduc_update_timestep(ns,mduc)
    implicit none

    type(mduc_t),pointer, intent(inout) :: mduc
    integer, intent(inout) :: ns
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCupdatetimestep(ns, mduc%mem)
#endif
    
    return
  end subroutine mduc_update_timestep

  ! ========================================================== !
  !> @brief cleanup mduc memory                                !
  !                                                            !
  ! ========================================================== !

  subroutine mduc_finalize(mduc)
    implicit none

    type(mduc_t),pointer, intent(inout) :: mduc
    ! ---------------------------------

    ! ---------------------------------

#ifdef MDUC
    call MDUCdestroy(mduc%mem)
#endif
    
    call destroy_mduc(mduc)

    return
  end subroutine mduc_finalize

  ! ========================================================== !
  !> @brief utility functions that do not actually call mduc
  ! ========================================================== !

  ! --------------------------------------------------------------- !
  !> @brief find grid point to the left of some (normalized) value
  ! --------------------------------------------------------------- !
  function mduc_grid_locate(nx,x,phi) result(ind)
    implicit none

    integer, intent(in) :: nx
    real(WP), dimension(nx), intent(in) :: x
    real(WP), intent(in) :: phi
    ! ---------------------------------
    integer :: jL,jR,jM,ind
    real(WP) :: asc
    ! ---------------------------------

    ! use bisection method to find nearest grid point for phi
    jL = 0; jR = nx - 1
    ! account for non-ascending grids
    asc = sign(1.0_WP,x(nx) - x(1))

    do while(jR-jL.gt.1)
       ! calc midpoint
       jM = (jL + jR)/2
       if (asc*(phi-x(jM)).ge.0.0_WP) then
          jL = jM
       else
          jR = jM
       end if
    end do

    ind = max(1,min(nx-1,jL))
    return
  end function  mduc_grid_locate

  ! ================================================================== !
  !> @brief some utility functions to ensure C-fortran compatibility   !
  !                                                                    !
  ! ================================================================== !

  ! ------------------------------------------------------------------------ !
  !> @brief modify characters before passing to C !  don't use trim on the
  !! output string afterwards !
  ! ------------------------------------------------------------------------ !
  subroutine add_end_char(string_in, string_out)
    implicit none
    character(len=*), intent(in) :: string_in
    character(len=*), intent(inout) :: string_out
    ! ---------------------------------

    ! ---------------------------------

    if (len(string_out) .le. len_trim(string_in)) then
       write(*,*) 'Error: output string too short to append C end character'
       call abort()
    end if

    ! trim and add end character
    string_out = trim(string_in)//char(0)

    return
  end subroutine add_end_char

  ! ---------------------------------------------------------- !
  !> @brief modify characters after receiving from C           !
  ! ---------------------------------------------------------- !
  subroutine remove_end_char(string)
    implicit none

    character(len=*), intent(inout) :: string
    ! ---------------------------------
    integer :: tmp
    ! ---------------------------------

    ! remove end character
    tmp = scan(adjustl(string),char(0)) - 1
    string = string(1:tmp)

    return
  end subroutine remove_end_char
  

end module combust_mduc_interface
