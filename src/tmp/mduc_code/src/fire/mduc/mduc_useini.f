c=======================================================================
      SUBROUTINE mduc_useini(mat,mph)
c=======================================================================
c
c     USEINI IS A SPECIAL PURPOSE ROUTINE FOR INITIALISAION 
c     BASED ON USER CODING
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c-----------------------------------------------------------------------
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : nsp,nep,ncell,nsb,neb,nbfac,itst,i_useini,
     x                     i_out_scal,num_actscl,num_passcl,time,dt,
     x                     iter,pref,nreg,iener,tref_h,ispray,
     x                     resta,iters,timo,irezon,iampro,numpro
      USE comm1,    ONLY : hs,hsb,t,den,tb,denb,actphi,actphb,pasphi,
     x                     pasphb,isol_actscl,iout_actscl,str_actscl,
     x                     isol_passcl,iout_passcl,str_passcl,actphr,
     x                     actpho,actpoo,pasphr,paspho,paspoo,iter,
     x                     p,vol,u,ub,ner,nsr,cpcof,cpcofb,ibc,te,ed
      USE cthmod,   ONLY : ymfrac,ymfrab,index_fuel,index_o2,index_n2,
     x                     index_co2,index_h2o,index_co,numscl,ymfstr,
     x                     iscsol,wmolar,enthfo,ymolfr,ymolfb
      USE usedefpar
      USE mduc
      USE mduc_sdm
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mat,mph
c-----
      INTEGER    :: nc,nb,istat,isc,tmp_test,i1,w_sum
      INTEGER    :: nZ1,nZ2,nSpecies,allSpecies,jZ,iY
      INTEGER    :: ir,ib,i
      REAL(prec) :: val0,val1,Havg_ref
      REAL(prec) :: hs_chem,wmix,usqr,cpmix,cpref
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: hi,cpi,Yox,Yfuel
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: Z1,Z2
      CHARACTER(len=256) :: species_id
      CHARACTER(len=256) :: stri1, stri2, text1, text2
c-----
c----- local variables for MDUC
      REAL(prec) :: Pavg, Zmin, Zmax
      REAL(prec) :: mass, totVolume, totMass
c-----------------------------------------------------------------------
c-----
c----- activation of useini (=1) in the GUI 
      IF(I_USEINI == 1 .AND. itst < 2) THEN ! 
c-----------------------------------------------------------------------
c----- set up the user defined scalars
c-----------------------------------------------------------------------
        CALL FIO_UUMESS('USEINI','I','Initialising Scalars for MDUC','')
        IF(num_actscl > mduc_actoff) THEN ! number of active scalar > 0
c-----solve active scalar transport equation
          isol_actscl(:)= 1  ! 0= solving no transport equation; 1= solve it
c-----
c----- active scalar output strings       ! string names for post-processing
          i = 1
          str_actscl(mduc_actoff+i)= 'Z1_mean';     i=i+1
          str_actscl(mduc_actoff+i)= 'Z1_variance'; i=i+1
          str_actscl(mduc_actoff+i)= 'H_mean';      i=i+1
          IF (mduc_nInj .EQ. 2) THEN
            str_actscl(mduc_actoff+i)= 'Z2_mean';     i=i+1
            str_actscl(mduc_actoff+i)= 'Z2_variance'; i=i+1
            str_actscl(mduc_actoff+i)= 'Z0_variance'; i=i+1
          END IF
          IF (mduc_soot .GT. 0) THEN
            str_actscl(mduc_actoff+i) = 'soot_n_density';        i=i+1
            str_actscl(mduc_actoff+i) = 'soot_volume_frac';      i=i+1
            str_actscl(mduc_actoff+i) = 'soot_surface_area';     i=i+1
            IF (mduc_soot .GT. 1) THEN
              str_actscl(mduc_actoff+i) = 'soot_volume_var';     i=i+1
              str_actscl(mduc_actoff+i) = 'soot_vol_surf_covar'; i=i+1
              str_actscl(mduc_actoff+i) = 'soot_surface_var';    i=i+1
            END IF
            str_actscl(mduc_actoff+i) = 'soot_delta_func';       i=i+1
          END IF
c----- active scalar indices
          CALL mduc_set_actphi_index
        END IF
c-----
        IF(num_passcl > mduc_pasoff) THEN       ! number of passive scalar > 0
c----- solve passive scalar transport equation
          isol_passcl(:)= 0  ! 0= solving no transport equation; 1= solve it
c-----      
c----- passive scalar output strings      ! string names for post-processing
          i = 1
c         fuel properties
          str_passcl(mduc_pasoff+i) = 'liq_vol_1';   i=i+1
          str_passcl(mduc_pasoff+i) = 'mass_evap_1'; i=i+1
c         enthalpy terms
          str_passcl(mduc_pasoff+i) = 'h_chem';      i=i+1
          str_passcl(mduc_pasoff+i) = 'HR_rate';     i=i+1
          str_passcl(mduc_pasoff+i) = 'hs';          i=i+1
c         reference scalar dissipation rates for each mixture
          str_passcl(mduc_pasoff+i) = 'chi1_ref';    i=i+1
c         variance source terms
          str_passcl(mduc_pasoff+i) = 'Z1_var_src';  i=i+1
          str_passcl(mduc_pasoff+i) = 'Z1_var_eps';  i=i+1
c         turb. diffusivity coefficients
          str_passcl(mduc_pasoff+i) = 'fdt_1';       i=i+1
c         exponent in SDM evap.
          str_passcl(mduc_pasoff+i) = 'xi_1';        i=i+1
c         coefficient in SDM evap.
          str_passcl(mduc_pasoff+i) = 'alpha_1';     i=i+1
c         Spalding number for evaporation
          str_passcl(mduc_pasoff+i) = 'spalding_no'; i=i+1
c         variance source terms (production, evaporation, dissipation)
          IF (mduc_outSpraySrc .EQ. 1) THEN
            str_passcl(mduc_pasoff+i) = 'Z1V_prod';    i=i+1
            str_passcl(mduc_pasoff+i) = 'Z1V_evap' ;   i=i+1    
            str_passcl(mduc_pasoff+i) = 'Z1V_diss';    i=i+1   
          END IF
          IF (mduc_nInj .EQ. 2) THEN
c           fuel properties
            str_passcl(mduc_pasoff+i) = 'liq_vol_2';   i=i+1
            str_passcl(mduc_pasoff+i) = 'mass_evap_2'; i=i+1
c           reference scalar dissipation rates for each mixture
            str_passcl(mduc_pasoff+i) = 'chi0_ref';    i=i+1
            str_passcl(mduc_pasoff+i) = 'chi2_ref';    i=i+1
c           variance source terms
            str_passcl(mduc_pasoff+i) = 'Z0_var_src';  i=i+1
            str_passcl(mduc_pasoff+i) = 'Z2_var_src';  i=i+1
c           variance dissipation terms
            str_passcl(mduc_pasoff+i) = 'Z0_var_eps';  i=i+1
            str_passcl(mduc_pasoff+i) = 'Z2_var_eps';  i=i+1
c           turb. diffusivity coefficients
            str_passcl(mduc_pasoff+i) = 'fdt_0';       i=i+1
            str_passcl(mduc_pasoff+i) = 'fdt_2';       i=i+1
c           exponent in SDM evap.
            str_passcl(mduc_pasoff+i) = 'xi_0';        i=i+1
            str_passcl(mduc_pasoff+i) = 'xi_2';        i=i+1
c           coefficient in SDM evap.
            str_passcl(mduc_pasoff+i) = 'alpha_0';     i=i+1
            str_passcl(mduc_pasoff+i) = 'alpha_2';     i=i+1
c           variance source terms (production, evaporation, dissipation)
            IF (mduc_outSpraySrc .EQ. 1) THEN
              str_passcl(mduc_pasoff+i) = 'Z0V_prod';    i=i+1
              str_passcl(mduc_pasoff+i) = 'Z2V_prod';    i=i+1
              str_passcl(mduc_pasoff+i) = 'Z0V_evap';    i=i+1
              str_passcl(mduc_pasoff+i) = 'Z2V_evap';    i=i+1
              str_passcl(mduc_pasoff+i) = 'Z0V_diss';    i=i+1
              str_passcl(mduc_pasoff+i) = 'Z2V_diss';    i=i+1
            END IF
          END IF
          IF (mduc_soot .GT. 0) THEN
            str_passcl(mduc_pasoff+i) = 'soot_n_particle'; i=i+1
            str_passcl(mduc_pasoff+i) = 'soot_volume';     i=i+1
            str_passcl(mduc_pasoff+i) = 'soot_area';       i=i+1
          END IF
          IF (mduc_outSootSrc .EQ. 1) THEN
            str_passcl(mduc_pasoff+i) = 'soot_src_n_density';   i=i+1
            str_passcl(mduc_pasoff+i) = 'soot_src_v_frac';      i=i+1
            str_passcl(mduc_pasoff+i) = 'soot_src_s_area';      i=i+1
            IF (mduc_soot .GT. 1) THEN
              str_passcl(mduc_pasoff+i) = 'soot_src_v_var';     i=i+1
              str_passcl(mduc_pasoff+i) = 'soot_src_v_s_covar'; i=i+1
              str_passcl(mduc_pasoff+i) = 'soot_src_s_var';     i=i+1
            END IF
            str_passcl(mduc_pasoff+i) = 'soot_src_delta_func';  i=i+1
          END IF
          IF (mduc_outSootTP .EQ. 1) THEN
            str_passcl(mduc_pasoff+i) = 'soot_src_n_density_tp';   i=i+1
            str_passcl(mduc_pasoff+i) = 'soot_src_v_frac_tp';      i=i+1
            str_passcl(mduc_pasoff+i) = 'soot_src_s_area_tp';      i=i+1
            IF (mduc_soot .GT. 1) THEN
              str_passcl(mduc_pasoff+i) = 'soot_src_v_var_tp';     i=i+1
              str_passcl(mduc_pasoff+i) = 'soot_src_v_s_covar_tp'; i=i+1
              str_passcl(mduc_pasoff+i) = 'soot_src_s_var_tp';     i=i+1
            END IF
            str_passcl(mduc_pasoff+i) = 'soot_src_delta_func_tp';  i=i+1
          END IF
          IF (mduc_all_species .EQ. 0) THEN
            str_passcl(mduc_pasoff+i) = 'Cpmix';     i=i+1
            str_passcl(mduc_pasoff+i) = 'dCpmix_dT'; i=i+1
            str_passcl(mduc_pasoff+i) = 'T_Cp';      i=i+1
          END IF

c----- passive scalar indices
          CALL mduc_set_pasphi_index
        END IF
c-----
        IF(num_actscl > 0) THEN           ! number of active scalar > 0
          DO nc=nsp(mat),nep(mat)         ! loop over internal cell
            DO isc=1,num_actscl
              actphi(nc+(isc-1)*ncell) = 0.0
            END DO
          END DO
          DO nb=nsb(mat),neb(mat)         ! loop over boundary faces
            DO isc=1,num_actscl
              actphb(nb+(isc-1)*nbfac) = 0.0
            END DO
          END DO
        END IF
c-----
        IF(num_passcl > 0) THEN           ! number of passive scalar > 0
          DO nc=nsp(mat),nep(mat)         ! loop over internal cell
            DO isc=1,num_passcl
               pasphi(nc+(isc-1)*ncell)= 0.0 ! set all passive scalars to zero
            END DO
          END DO
          DO nb=nsb(mat),neb(mat)       ! loop over boundary faces
            DO isc=1,num_passcl
               pasphb(nb+(isc-1)*nbfac)= 0.0 ! set all passive scalars to zero
            END DO
          END DO
        END IF
c-----
c----- make sure that all processes have an updated species list
        IF (numpro .gt. 1) THEN
          DO iY=1,numscl
            CALL cbcast(0, ymfstr(iY), len(ymfstr(iY)))
          END DO
        END IF
c-----------------------------------------------------------------------
c----- call initialisation of MDUC
c-----------------------------------------------------------------------
c----- get parameters from the GUI
c      input file
        i1 = UDP_GET_INDEX('MDUC_INPUT_FILE')
        IF (i1 > 0 .AND. i1 <= n_udp) THEN
          mduc_infile = udp(i1)%c_value(1:udp(i1)%i_lstr_value)
        END IF
c      check the enthalpy formulation
        i1 = UDP_GET_INDEX('Enthalpy_formulation')
        IF (i1 > 0 .AND. i1 <= n_udp) THEN
          H_method_name = udp(i1)%c_value(1:udp(i1)%i_lstr_value)
        END IF
        IF (H_method_name .eq. "total") THEN
          H_method = 0
        ELSE IF (H_method_name .eq. "sensible") THEN
          H_method = 1
        ELSE IF (H_method_name .eq. "scalar") THEN
          H_method = 2
        END IF
c      check for explicit initialization times
        i1 = UDP_GET_INDEX('t_init_pilot')
        IF (i1 > 0 .AND. i1 <= n_udp) THEN
          mduc_tpilot = UDP_DOUBLE_VALUE(i1)
        END IF
        i1 = UDP_GET_INDEX('t_init_main')
        IF (i1 > 0 .AND. i1 <= n_udp) THEN
          mduc_tmain = UDP_DOUBLE_VALUE(i1)
        END IF
c----- call the initialisation function
        text1 = trim(mduc_infile)//char(0) ! trim and add end character
        CALL MDUCinit(text1,mduc_mem)
c----- call parallel initialisation if necessary
        IF (numpro > 1) THEN
          MPI_COMM_MDUC = 0
          rootpro = 1
          CALL MDUCsetParallel(MPI_COMM_MDUC, mduc_mem)
        ELSE
          rootpro = 0
        END IF
c----- set the time from the solver to the flamelet
c     TODO: check if update time necessary here
        CALL MDUCupdateTime(time,dt,itst,mduc_mem)
c----- index the species in the flow domain to the flamelet
        ALLOCATE(indexY(numscl)); indexY    = -1
c       turn off solution of all scalars and find names
        index_fuel = 1
        DO iY=1,numscl
          iscsol(iY) = 0
          species_id = trim(ymfstr(iY))//char(0)
c          write(*,*) 'useini: species: ', species_id
          indexY(iY) = MDUCgetSpeciesIndex(species_id,mduc_mem)
c         set the molecular weights to be consistent with the mechanism
          IF (indexY(iY) .GE. 0) THEN
            wmolar(iY) = MDUCgetMolWeight(species_id,mduc_mem)
          END IF
        END DO
c----- store the indices for the necessary species for the soot model
        IF (mduc_soot .GT. 0) THEN
c         number of soot moments set in usedef
c         define which frequency the soot source terms should be updated
          IF (sootSrcUpdateDef .EQ. 'timestep') THEN
            sootSrcUpdate = 0
          ELSE IF (sootSrcUpdateDef .EQ. 'iteration') THEN
            sootSrcUpdate = 1
          END IF
c         define the species necessary for the soot model
          nSootSpecies = 7
c         set the species names (these must match the order in hmom.h)
          sootSpecies(1) = "H"
          sootSpecies(2) = "H2"
          sootSpecies(3) = "OH"
          sootSpecies(4) = "H2O"
          sootSpecies(5) = "O2"
          sootSpecies(6) = "C2H2"
          sootSpecies(7) = "A4-C16H10"
c         initialize indices to -1
          DO i=1,nSootSpecies
            sootIndx(i) = -1
          END DO
c         search for species
          DO iY=1,numscl
            DO i=1,nSootSpecies
              IF (ymfstr(iY) .EQ. sootSpecies(i)) sootIndx(i) = iY
            END DO
          END DO
c         error check that all necessary species are present
          DO i=1,nSootSpecies
            IF (sootIndx(i) .EQ. -1) THEN
              WRITE (*,*) 'Error: species ', sootSpecies(i),
     x             ' necessary for soot model'
            END IF
          END DO
          CALL mduc_soot_init_moments(mat)
        END IF
c----- compute the average flow field quantities required
        Pavg      = 0.0
        totVolume = 0.0; totMass   = 0.0
        DO nc=nsp(mat),nep(mat)
          Pavg      = Pavg      + (pref(mat) + p(nc))*vol(nc)
          totVolume = totVolume + vol(nc)
          totMass   = totMass   + den(nc)*vol(nc)
        END DO
        CALL dglsum(Pavg)
        CALL dglsum(totVolume)
        Pavg = Pavg/totVolume
c       set the global average pressure 
        Pavg_old = Pavg
        Zmin = 0.0
        Zmax = 1.0
        dPdt = 0.0
        CALL MDUCinitFlow(time,Pavg,Zmin,Zmax,Zmin,Zmax,mduc_mem)
c----- write some info to the file
c----- get the grid size and number of species
        nZ1      = MDUCgetGridSize(1, mduc_mem)
        nZ2      = MDUCgetGridSize(2, mduc_mem)
        nSpecies = MDUCgetNspecies(mduc_mem)
c        mduc_nspecies = nSpecies
c        CALL FIO_UUMESS2('USEINI','Grid points','R',nZ1,'-')
c        CALL FIO_UUMESS2('USEINI','Species','R',nSpecies,'-')
c----- allocate and retrieve the grid
        ALLOCATE(Z1(nZ1), STAT=istat); Z1 = 0.0
        ALLOCATE(Z2(nZ2), STAT=istat); Z2 = 0.0
c      set the zero point for the grid
        IF (nZ2 .eq. 1) THEN
          CALL MDUCgetGrid(Z1,mduc_mem)
          Z1_zero = 0.5*(Z1(2)-Z1(1))
          Z2_zero = 1.0e-6
        ELSE
          CALL MDUCgetGrid2D(Z1, Z2, mduc_mem)
          Z1_zero = 0.5*(Z1(2)-Z1(1))
          Z2_zero = 0.5*(Z2(2)-Z2(1))
        END IF
c
        IF(ALLOCATED(Z1)) DEALLOCATE(Z1)
        IF(ALLOCATED(Z2)) DEALLOCATE(Z2)
c----- compute the reference enthalpy of the fuel vapour
        dPdt = 0.0
        ALLOCATE(hi(nSpecies),     STAT=istat); hi     = 0.0
        ALLOCATE(cpi(nSpecies),    STAT=istat); cpi    = 0.0
        ALLOCATE(Yox(nSpecies),    STAT=istat); Yox    = 0.0
        ALLOCATE(Yfuel(nSpecies),  STAT=istat); Yfuel  = 0.0
        val0 = 0.0
        val1 = 1.0
        CALL MDUCgetPointY(val0,val0,Yox,mduc_mem)
        CALL MDUCgetPointY(val1,val0,Yfuel,mduc_mem)
        T_drop_ref = 298.15
        H_vap_ref  = 0.0
        CALL MDUCcomputeThermoData(hi,cpi,T_drop_ref,mduc_mem)
        DO iY=1,nSpecies
          H_vap_ref = H_vap_ref + hi(iY)*Yfuel(iY)
        END DO
c-----
c     set the injection to 0 (not solve)
        mduc_inj = 0
c-----------------------------------------------------------------------
c----- initialise the Single Droplet Model
c-----------------------------------------------------------------------
        IF (ispray == 1) THEN
          CALL sdmInit()
        END IF
c-----------------------------------------------------------------------
c----- initialise total enthalpy (scalar) everywhere in domain
c-----------------------------------------------------------------------
        Havg_ref = 0.0
        totMass  = 0.0
        DO nc=nsp(mat),nep(mat)
          actphi(nc+iHTOT) = 0.0
          CALL MDUCcomputeThermoData(hi,cpi,t(nc),mduc_mem)
          DO i=1,numscl
            mass = den(nc)*vol(nc)
            IF (indexY(i) >= 0) THEN
              actphi(nc+iHTOT) = actphi(nc+iHTOT)+
     x             hi(indexY(i))*ymfrac(nc+(i-1)*ncell)
            END IF
          END DO
          Havg_ref  = Havg_ref + actphi(nc+iHTOT)*mass
          totMass   = totMass  + mass
        END DO
        Havg_ref = Havg_ref/totMass
c----- update the boundary faces
        DO ir=0,nreg
          IF (ibc(2,ir) == mat) THEN
            DO ib=nsr(ir),ner(ir)
              actphb(ib+bHTOT) = 0.0
              CALL MDUCcomputeThermoData(hi,cpi,tb(ib),mduc_mem)
              DO i=1,numscl
                IF (indexY(i) >= 0) THEN
                  actphb(ib+bHTOT) = actphb(ib+bHTOT)+
     x                 hi(indexY(i))*ymfrab(ib+(i-1)*nbfac)
                END IF
              END DO
            END DO
          END IF
        END DO
c-----------------------------------------------------------------------
c----- set the initial enthalpy the enthalpy
c-----------------------------------------------------------------------
c$$$        open(unit=17,file="enthalpy_compare.dat")
c$$$        write(17,*) '# time[s] itst iter cell h_tot actphi perc_diff hs
c$$$     x       usqr pasphi '
c$$$        close(unit=17)
c        tref_h = 298.15
        IF (H_method .EQ. 0) THEN
          IF (iampro < 2) THEN
            WRITE (*,*) 'MDUC_init start: updating initial enthalpy'
          END IF
          CALL update_enthalpy_total(mat)
        ELSE IF (H_method .EQ. 1 .OR. H_method .EQ. 2) THEN
          CALL update_enthalpy_sensible(mat)
        END IF
c----- store the FIRE internal enthalpy in a passive scalar to output it
        DO nc=nsp(mat),nep(mat)
          pasphi(nc+iHS) = hs(nc)
        END DO
        DO ir=0,nreg
          IF (ibc(2,ir) == mat) THEN
            DO ib=nsr(ir),ner(ir)
              pasphb(ib+bHS) = hsb(ib)
            END DO
          END IF
        END DO
c-----
        IF(ALLOCATED(hi))     DEALLOCATE(hi);
        IF(ALLOCATED(cpi))    DEALLOCATE(cpi);
        IF(ALLOCATED(Yox))    DEALLOCATE(Yox);
        IF(ALLOCATED(Yfuel))  DEALLOCATE(Yfuel);
c-----------------------------------------------------------------------
c-----  write a summary
c-----------------------------------------------------------------------
        IF (iampro < 2) THEN
          w_sum = 80
          WRITE(*,*) ''
          CALL mduc_out_header(w_sum, "MDUC Initialization Summary")
          WRITE(text2, '(A)') REPEAT(' ', w_sum-2)
          CALL mduc_out(w_sum, text2)

          WRITE(text1, '(A30,A)') '  Input File : ',trim(mduc_infile)
          CALL mduc_out(w_sum, text1)
          CALL mduc_out(w_sum, text2)
c          WRITE (*,*) '* MDUC grid points: ', nZ1, nZ2
          WRITE(text1,'(A30,A)') '  Enthalpy form. : ',
     x         trim(H_method_name)
          CALL mduc_out(w_sum, text1)
          WRITE(text1,'(A30,I)') '  Enthalpy method : ', H_method
          CALL mduc_out(w_sum, text1)
c
          WRITE(text1,'(A30,ES12.5,A)') '  Initial Pressure : ',
     x         Pavg, ' [Pa]'
          CALL mduc_out(w_sum, text1)
c
          CALL mduc_out(w_sum, text2)
          WRITE(text1,'(A30,F8.4,A)') '  Pilot time : ',
     x         mduc_tpilot,' [s]'
          CALL mduc_out(w_sum, text1)
          WRITE(text1,'(A30,F8.4,A)') '  Main time : ',
     x         mduc_tmain,' [s]'
          CALL mduc_out(w_sum, text1)
          WRITE(text1,'(A30,F8.2,A)') '  Ref. Vap. Enthalpy : ',
     x         H_vap_ref/1e3, ' [kJ/kg]'
          CALL mduc_out(w_sum, text1)
c
          CALL mduc_out(w_sum, text2)
          WRITE(text1, '(A)') REPEAT('*', w_sum-2)
          CALL mduc_out(w_sum, text1)
c
          WRITE(*,*) ''
          WRITE(*,*) ''
c
c         write the chemistry info
c
          w_sum = 40
          CALL mduc_out_header(w_sum, "Chemical Species Mapping")
          CALL mduc_out(w_sum, text2)
          DO iY=1,numscl
            WRITE(text1,'(A25,A,I5)') trim(ymfstr(iY)),' : ',indexY(iY)
            CALL mduc_out(w_sum, text1)
          END DO
          CALL mduc_out(w_sum, text2)
          WRITE(text1, '(A)') REPEAT('*', w_sum-2)
          CALL mduc_out(w_sum, text1)
c
c         write soot info
c
          IF (mduc_soot .GT. 0) THEN
            w_sum = 60
            WRITE(*,*) ''
            WRITE(*,*) ''
            CALL mduc_out_header(w_sum, "Soot Model")
            CALL mduc_out(w_sum, text2)
            WRITE(text1,'(A28,I)') '  Soot moments : ', nSootMom
            CALL mduc_out(w_sum, text1)
            CALL mduc_out(w_sum, text2)
            WRITE(text1,'(A28,I)') '  Soot species : ', nSootSpecies
            CALL mduc_out(w_sum, text1)
            DO iY=1,nSootSpecies
              WRITE(text1,'(A25,A,I5)') trim(sootSpecies(iY)), ' : ',
     x             sootIndx(iY)
              CALL mduc_out(w_sum, text1)
            END DO
            
            WRITE(text1, '(A)') REPEAT('*', w_sum-2)
            CALL mduc_out(w_sum, text1)
            CALL mduc_out(w_sum, text2)
          END IF
c
        END IF

c        DO nc=nsp(mat),nep(mat)
c          write(*,*) 'tke: ', te(nc), ' eps: ', ed(nc)
c        END DO
c-----------------------------------------------------------------------
c----- mduc_useini call after restart file has been read
c-----------------------------------------------------------------------
      ELSE IF ((I_USEINI == 1 .OR. I_USEINI == 99) .AND. resta == 1 
     x       .AND. iter == iters) THEN
        IF (iampro < 2) THEN
          WRITE(*,*) 'MDUC_init_restart: updating to restart conditions'
        END IF
c-----make sure that the scalar indices are up-to-date
        CALL mduc_set_scalar_index
c-----make sure that time step number is updated for restart
        isc = itst-1
        CALL MDUCupdateTime(timo,dt,isc,mduc_mem)
c-----set Pavg_old to the right value for restart
        Pavg_old  = 0.0
        totVolume    = 0.0
        DO nc=nsp(mat),nep(mat)
          Pavg_old = Pavg_old + (pref(mat)+p(nc))*vol(nc)
          totVolume   = totVolume + vol(nc)
        END DO
c     reduce across processors
        CALL dglsum(Pavg_old)
        CALL dglsum(totVolume)
c     finish normalization
        Pavg_old = Pavg_old/totVolume
        Zmin = 0.0
        Zmax = 1.0
        CALL MDUCinitFlow(timo,Pavg_old,Zmin,Zmax,
     x       Zmin,Zmax,mduc_mem)
c
c-----------------------------------------------------------------------
c----- initialise total enthalpy (scalar) everywhere in domain
c-----------------------------------------------------------------------
        nSpecies = MDUCgetNspecies(mduc_mem)
        ALLOCATE(hi(nSpecies),     STAT=istat); hi     = 0.0
        ALLOCATE(cpi(nSpecies),     STAT=istat); cpi     = 0.0
        Havg_ref = 0.0
        totMass  = 0.0
        DO nc=nsp(mat),nep(mat)
          actphi(nc+iHTOT) = 0.0
          CALL MDUCcomputeThermoData(hi,cpi,t(nc),mduc_mem)
          DO i=1,numscl
            mass = den(nc)*vol(nc)
            IF (indexY(i) >= 0) THEN
              actphi(nc+iHTOT) = actphi(nc+iHTOT)+
     x             hi(indexY(i))*ymfrac(nc+(i-1)*ncell)
            END IF
          END DO
          Havg_ref  = Havg_ref + actphi(nc+iHTOT)*mass
          totMass   = totMass  + mass
        END DO
        Havg_ref = Havg_ref/totMass
c----- update the boundary faces
        DO ir=0,nreg
          IF (ibc(2,ir) == mat) THEN
            DO ib=nsr(ir),ner(ir)
              actphb(ib+bHTOT) = 0.0
              CALL MDUCcomputeThermoData(hi,cpi,tb(ib),mduc_mem)
              DO i=1,numscl
                IF (indexY(i) >= 0) THEN
                  actphb(ib+bHTOT) = actphb(ib+bHTOT)+
     x                 hi(indexY(i))*ymfrab(ib+(i-1)*nbfac)
                END IF
              END DO
            END DO
          END IF
        END DO
c-----------------------------------------------------------------------
c----- set the initial enthalpy
c-----------------------------------------------------------------------
        IF (H_method .EQ. 0) THEN
          IF (iampro < 2) THEN
            WRITE (*,*) 'MDUC_init_restart: updating initial enthalpy'
          END IF
          CALL update_enthalpy_total(mat)
        ELSE IF (H_method .EQ. 1 .OR. H_method .EQ. 2) THEN
          CALL update_enthalpy_sensible(mat)
        END IF
c----- store the FIRE internal enthalpy in a passive scalar to output it
        DO nc=nsp(mat),nep(mat)
          pasphi(nc+iHS) = hs(nc)
        END DO
        DO ir=0,nreg
          IF (ibc(2,ir) == mat) THEN
            DO ib=nsr(ir),ner(ir)
              pasphb(ib+bHS) = hsb(ib)
            END DO
          END IF
        END DO
c-----
        IF(ALLOCATED(hi))     DEALLOCATE(hi);
        IF(ALLOCATED(cpi))    DEALLOCATE(cpi);
c
      ELSE IF (irezon .eq. 1) THEN ! update after rezoning
        WRITE (*,*) 'MDUC: updating after rezone'
c-----make sure that the scalar indices are up-to-date
        CALL mduc_set_scalar_index
c-----update the enthalpy (this should be pointer function)        
        IF (H_method .eq. 0) THEN
          CALL update_enthalpy_total(mat)
        ELSE IF (H_method .eq. 1 .or. H_method .eq. 2) THEN
          CALL update_enthalpy_sensible(mat)
        END IF
c     ELSE
c         WRITE (*,*) 'MDUC check: I_USEINI =', I_USEINI 
      END IF ! the following code is called before each time step
c-----  reset mean mixture fraction sources for next spray loop
      DO nc=nsp(mat),nep(mat)
        pasphi(nc+iMEVAP1) = 0.0
        IF (mduc_nInj .EQ. 2) pasphi(nc+iMEVAP2) = 0.0
      END DO
c----- set up default values for scalar dissipation rate model
      chiNorm(1) = 3     
c     TODO: need to set up user options      
      RETURN
      END SUBROUTINE mduc_useini
c-----------------------------------------------------------------------
c-----
c=======================================================================
      SUBROUTINE update_enthalpy_total(mat)
c=======================================================================
c
c     Routine to compute the enthalpy from the current temperature and
c     species mass fractions.  This formulation updates the variable
c     hs to be the total enthalpy, defined by:
c
c        hs = cp*(T-T_ref) + 0.5*u^2 + h_chem
c
c-----
c-----------------------------------------------------------------------
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : nsp,nep,ncell,nsb,neb,nbfac,nreg,iener,
     x                     tref_h,msin,i_usecps,msin
      USE comm1,    ONLY : hs,hsb,t,tb,actphi,actphb,pasphi,pasphb,
     x                     p,vol,u,ub,ner,nsr,cpcof,cpcofb,ibc
      USE cthmod,   ONLY : numscl,ymfrac,ymfrab,wmolar,enthfo,
     x                     ymolfr,ymolfb
      USE mduc
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mat
      INTEGER    :: istat,nc,i,ib,ir,nSpecies
      REAL(prec) :: usqr,h_chem,wmix,cpmix,cpP,cpM,delT,tmp
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: hi,cpi,hr,cpr,cpiP,cpiM
c-----
c----- allocate memory for species enthalpy/cp
      nSpecies = MDUCgetNspecies(mduc_mem)
      ALLOCATE(hi(nSpecies),   STAT=istat); hi  = 0.0
      ALLOCATE(cpi(nSpecies),  STAT=istat); cpi = 0.0
      ALLOCATE(hr(nSpecies),   STAT=istat); hr  = 0.0
      ALLOCATE(cpr(nSpecies),  STAT=istat); cpr = 0.0
      ALLOCATE(cpiP(nSpecies), STAT=istat); cpiP = 0.0
      ALLOCATE(cpiM(nSpecies), STAT=istat); cpiM = 0.0
c----- calc formation enthalpies
      H_Tref = 298.15
      CALL cthc_h_ck(1,numscl,H_Tref,enthfo) ! should already be defined
c----- compute reference thermo data from MDUC mechanism
      CALL MDUCcomputeThermoData(hr,cpr,H_Tref,mduc_mem)
      DO i=1,numscl
        IF (indexY(i) .GE. 0) THEN
          enthfo(i) = hr(indexY(i))*wmolar(i)/1.d3
        END IF
      END DO
c-----  update internal cells
      DO nc=nsp(mat),nep(mat)
        usqr    = DOT_PRODUCT(u(:,nc), u(:,nc))
        usqr    = usqr*iener
c-----  compute the chemical enthalpy
        h_chem = 0.0
        wmix   = 0.0
        CALL MDUCcomputeThermoData(hi,cpi,t(nc),mduc_mem)
        DO i=1,numscl
          IF (indexY(i) .GE. 0) THEN
            h_chem = h_chem + ymolfr(nc+(i-1)*ncell)*enthfo(i) ! kJ/kmol
            wmix   = wmix   + ymolfr(nc+(i-1)*ncell)*wmolar(i) ! kg/kmol
          END IF
        END DO
        h_chem = (h_chem/wmix)*1.d3 ! J/kg
c-----  store chemical enthalpy for temperature update
        pasphi(nc+iHCHEM) = h_chem
c-----  store the specific heat if necessary
        IF (mduc_all_species .EQ. 0 ) THEN
          delT   = 10
          cpmix  = 0.0
          cpP    = 0.0
          cpM    = 0.0
          CALL MDUCcomputeThermoData(hi,cpiP,t(nc)+delT,mduc_mem)
          CALL MDUCcomputeThermoData(hi,cpiM,t(nc)-delT,mduc_mem)
          DO i=1,numscl
            IF (indexY(i) .GE. 0) THEN
              cpmix  = cpmix  + ymfrac(nc+(i-1)*ncell)*cpi(indexY(i))
              cpP    = cpP    + ymfrac(nc+(i-1)*ncell)*cpiP(indexY(i))
              cpM    = cpM    + ymfrac(nc+(i-1)*ncell)*cpiM(indexY(i))
            END IF
          END DO
          pasphi(nc+iTCP)   = t(nc)
          pasphi(nc+iCPMIX) = cpmix
          pasphi(nc+iDCPDT) = (cpP - cpM)/(2*delT)
          cpcof(nc)         = cpmix
        END IF
c-----  compute the total enthalpy              
        hs(nc)=(t(nc)-H_Tref)*cpcof(nc) + 0.5*usqr + h_chem
      END DO
      CALL exchng(hs,1,1)
c-----
c----- update the boundary faces
      DO ir=0,nreg
        IF (ibc(2,ir) == mat) THEN
          DO ib=nsr(ir),ner(ir)
            usqr    = DOT_PRODUCT(ub(:,ib), ub(:,ib))
            usqr    = usqr*iener
c-----      compute the chemical enthalpy
            h_chem = 0.0
            wmix   = 0.0
            cpmix  = 0.0
            CALL MDUCcomputeThermoData(hi,cpi,tb(ib),mduc_mem)
            DO i=1,numscl
              IF (indexY(i) .GE. 0) THEN
                h_chem = h_chem + ymolfb(ib+(i-1)*nbfac)*enthfo(i) ! kJ/kmol
                wmix   = wmix   + ymolfb(ib+(i-1)*nbfac)*wmolar(i) ! kg/kmol
                cpmix  = cpmix  + ymfrab(ib+(i-1)*nbfac)*cpi(indexY(i))
              END IF
            END DO
            h_chem = (h_chem/wmix)*1.d3 ! J/kg
c-----      store chemical enthalpy for temperature update
            pasphb(ib+bHCHEM) = h_chem
c            cpcofb(ib)        = cpmix
c-----      compute the total enthalpy              
            hsb(ib)=(tb(ib)-H_Tref)*cpcofb(ib) + 0.5*usqr + h_chem
          END DO
        END IF
      END DO
      IF(ALLOCATED(hi))   DEALLOCATE(hi);
      IF(ALLOCATED(cpi))  DEALLOCATE(cpi);
      IF(ALLOCATED(hr))   DEALLOCATE(hr);
      IF(ALLOCATED(cpr))  DEALLOCATE(cpr);
      IF(ALLOCATED(cpiP)) DEALLOCATE(cpiP);
      IF(ALLOCATED(cpiM)) DEALLOCATE(cpiM);
c-----
      END SUBROUTINE update_enthalpy_total
c-----
c=======================================================================
      SUBROUTINE update_enthalpy_sensible(mat)
c=======================================================================
c
c     Routine to compute the enthalpy from the current temperature and
c     species mass fractions. This formulation implements the variable
c     hs as the sensible enthalpy directly, defined as:
c
c       hs = sum_i( Y_i*(h_i - h_i_ref) )
c
c-----
c-----------------------------------------------------------------------
      USE prec_mod, ONLY : prec
      USE comm0,  ONLY : nsp,nep,ncell,nsb,neb,nbfac,nreg,iener,tref_h,
     x                   i_usecps,msin
      USE comm1,  ONLY : hs,hsb,t,tb,actphi,actphb,pasphi,pasphb,
     x                   p,vol,u,ub,ner,nsr,cpcof,cpcofb,ibc
      USE cthmod, ONLY : numscl,ymfrac,ymfrab,wmolar,enthfo,
     x                   ymolfr,ymolfb
      USE mduc
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mat
      INTEGER    :: istat,nc,i,ib,ir,nSpecies
      REAL(prec) :: usqr,h_chem,Ym
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: hi,cpi,hr,cpr
c-----
c----- allocate memory for species enthalpy/cp
      nSpecies = MDUCgetNspecies(mduc_mem)
      ALLOCATE(hi(nSpecies),  STAT=istat); hi  = 0.0
      ALLOCATE(cpi(nSpecies), STAT=istat); cpi = 0.0
      ALLOCATE(hr(nSpecies),  STAT=istat); hr  = 0.0
      ALLOCATE(cpr(nSpecies), STAT=istat); cpr = 0.0
c----- calc formation enthalpies
      H_Tref = 298.15
      CALL cthc_h_ck(1,numscl,H_Tref,enthfo) ! should already be defined
c----- update specific heats
      IF (I_USECPS >= 1) THEN
        CALL usecps(mat,msin)
      ELSE
        CALL cthucp(1,ir)
      END IF
c----- compute reference thermo data from MDUC mechanism
      CALL MDUCcomputeThermoData(hr,cpr,H_Tref,mduc_mem)
c----- update the internal cells
      DO nc=nsp(mat),nep(mat)
        usqr    = DOT_PRODUCT(u(:,nc), u(:,nc))
        usqr    = usqr*iener
c-----  compute the sensible enthalpy
        hs(nc) = 0.0
        h_chem = 0.0
        CALL MDUCcomputeThermoData(hi,cpi,t(nc),mduc_mem)
        DO i=1,numscl
          IF (indexY(i) .GE. 0) THEN
            Ym     = ymfrac(nc+(i-1)*ncell)
            h_chem = h_chem + Ym*hr(indexY(i))
            hs(nc) = hs(nc) + Ym*(hi(indexY(i))-hr(indexY(i)))
          END IF
        END DO
c-----  add the kinetic energy (for now?)
        hs(nc) = hs(nc) + 0.5*usqr
c-----  compute a specific heat that is consistent with the enthalpy
c        cpcof(nc) = hs(nc)/(t(nc) - H_Tref)
c-----  store the chemical enthalpy for use in the temperature calc.
        pasphi(nc+iHCHEM) = h_chem
      END DO
c-----
c----- update the boundary faces
      DO ir=0,nreg
        IF (ibc(2,ir) == mat) THEN
          DO ib=nsr(ir),ner(ir)
            usqr    = DOT_PRODUCT(ub(:,ib), ub(:,ib))
            usqr    = usqr*iener
c-----      compute the sensible enthalpy
            hsb(ib) = 0.0
            h_chem  = 0.0
            CALL MDUCcomputeThermoData(hi,cpi,tb(ib),mduc_mem)
            DO i=1,numscl
              IF (indexY(i) .GE. 0) THEN
                Ym      = ymfrab(ib+(i-1)*nbfac)
                h_chem  = h_chem  + Ym*hr(indexY(i))
                hsb(ib) = hsb(ib) + Ym*(hi(indexY(i))-hr(indexY(i)))
              END IF
            END DO
c-----      add the kinetic energy (for now?)
            hsb(ib) = hsb(ib) + 0.5*usqr
c-----      compute a specific heat that is consistent with the enthalpy
c            cpcofb(ib) = hsb(ib)/(tb(ib) - H_Tref)
c-----      store the chemical enthalpy for use in the temperature calc.
            pasphb(ib+bHCHEM) = h_chem
          END DO
        END IF
      END DO
c-----
      IF(ALLOCATED(hi))  DEALLOCATE(hi);
      IF(ALLOCATED(cpi)) DEALLOCATE(cpi);
      IF(ALLOCATED(hr))  DEALLOCATE(hr);
      IF(ALLOCATED(cpr)) DEALLOCATE(cpr);
c-----
      END SUBROUTINE update_enthalpy_sensible
c-----
c=======================================================================
      SUBROUTINE mduc_set_scalar_index
c=======================================================================
c
c     Function to initialize indices corresponding to the offsets of
c     each scalar
c-----
c-----------------------------------------------------------------------
      USE comm0,  ONLY : ncell,nbfac
      USE mduc
      IMPLICIT NONE
c-----      
c     call the functions to set both active and passive scalar indices
      CALL mduc_set_actphi_index
      CALL mduc_set_pasphi_index
c-----
      END SUBROUTINE mduc_set_scalar_index
c-----------------------------------------------------------------------
c-----
c=======================================================================
      SUBROUTINE mduc_set_actphi_index
c=======================================================================
c
c     Function to initialize indices corresponding to the offsets of
c     active scalars
c-----
c-----------------------------------------------------------------------
      USE comm0,  ONLY : ncell,nbfac
      USE mduc
      IMPLICIT NONE
c-----
      INTEGER i,j
c-----------------------------------------------------------------------
      i = 0
      j = mduc_actoff
      iZ1MEAN = (i+j)*ncell; bZ1MEAN = (i+j)*nbfac; i=i+1
      iZ1VAR  = (i+j)*ncell; bZ1VAR  = (i+j)*nbfac; i=i+1
      iHTOT   = (i+j)*ncell; bHTOT   = (i+j)*nbfac; i=i+1
      IF (mduc_nInj .EQ. 2) THEN
        iZ2MEAN = (i+j)*ncell; bZ2MEAN = (i+j)*nbfac; i=i+1
        iZ2VAR  = (i+j)*ncell; bZ2VAR  = (i+j)*nbfac; i=i+1
        iZ0VAR  = (i+j)*ncell; bZ0VAR  = (i+j)*nbfac; i=i+1
      END IF
      IF (mduc_soot .GT. 0) THEN
        iSOOT_MOM = (i+j)*ncell; bSOOT_MOM = (i+j)*nbfac;
        iSOOT_M00 = (i+j)*ncell; bSOOT_M00 = (i+j)*nbfac; i=i+1
        iSOOT_M10 = (i+j)*ncell; bSOOT_M10 = (i+j)*nbfac; i=i+1
        iSOOT_M01 = (i+j)*ncell; bSOOT_M01 = (i+j)*nbfac; i=i+1
        IF (mduc_soot .GT. 1) THEN
          iSOOT_M20 = (i+j)*ncell; bSOOT_M20 = (i+j)*nbfac; i=i+1
          iSOOT_M02 = (i+j)*ncell; bSOOT_M02 = (i+j)*nbfac; i=i+1
          iSOOT_M11 = (i+j)*ncell; bSOOT_M11 = (i+j)*nbfac; i=i+1
        END IF
        iSOOT_N0 = (i+j)*ncell; bSOOT_N0 = (i+j)*nbfac; i=i+1
      END IF
c-----
      END SUBROUTINE mduc_set_actphi_index
c-----------------------------------------------------------------------
c-----
c=======================================================================
      SUBROUTINE mduc_set_pasphi_index
c=======================================================================
c
c     Function to initialize indices corresponding to the offsets of
c     passive scalars
c-----
c-----------------------------------------------------------------------
      USE comm0,  ONLY : ncell,nbfac
      USE mduc
      IMPLICIT NONE
c-----
      INTEGER i,j
c-----------------------------------------------------------------------
c     passive scalar offsets
      i = 0
      j = mduc_pasoff
      iLIQV1  = (i+j)*ncell;                          i=i+1
      iMEVAP1 = (i+j)*ncell;                          i=i+1      
      iHCHEM  = (i+j)*ncell; bHCHEM  = (i+j)*nbfac;   i=i+1
      iHRATE  = (i+j)*ncell; bHRATE  = (i+j)*nbfac;   i=i+1
      iHS     = (i+j)*ncell; bHS     = (i+j)*nbfac;   i=i+1
      iCHI1   = (i+j)*ncell; bCHI1   = (i+j)*nbfac;   i=i+1
      iZ1VSRC = (i+j)*ncell; bZ1VSRC = (i+j)*nbfac;   i=i+1
      iZ1VEPS = (i+j)*ncell; bZ1VEPS = (i+j)*nbfac;   i=i+1
      iFDT1   = (i+j)*ncell; bFDT1   = (i+j)*nbfac;   i=i+1
      iXI1    = (i+j)*ncell; bXI1    = (i+j)*nbfac;   i=i+1
      iALPH1  = (i+j)*ncell; bALPH1  = (i+j)*nbfac;   i=i+1
      iSPALD  = (i+j)*ncell; bSPALD  = (i+j)*nbfac;   i=i+1
      IF (mduc_outSpraySrc .EQ. 1) THEN
        iZ1VPR  = (i+j)*ncell; bZ1VPR  = (i+j)*nbfac;   i=i+1
        iZ1VEV  = (i+j)*ncell; bZ1VEV  = (i+j)*nbfac;   i=i+1
        iZ1VDI  = (i+j)*ncell; bZ1VDI  = (i+j)*nbfac;   i=i+1
      END IF
      IF (mduc_nInj .EQ. 2) THEN
        iLIQV2  = (i+j)*ncell;                          i=i+1
        iMEVAP2 = (i+j)*ncell;                          i=i+1
        iCHI0   = (i+j)*ncell; bCHI0   = (i+j)*nbfac;   i=i+1
        iCHI2   = (i+j)*ncell; bCHI2   = (i+j)*nbfac;   i=i+1
        iZ0VSRC = (i+j)*ncell; bZ0VSRC = (i+j)*nbfac;   i=i+1
        iZ2VSRC = (i+j)*ncell; bZ2VSRC = (i+j)*nbfac;   i=i+1
        iZ0VEPS = (i+j)*ncell; bZ0VEPS = (i+j)*nbfac;   i=i+1
        iZ2VEPS = (i+j)*ncell; bZ2VEPS = (i+j)*nbfac;   i=i+1
        iFDT0   = (i+j)*ncell; bFDT0   = (i+j)*nbfac;   i=i+1
        iFDT2   = (i+j)*ncell; bFDT2   = (i+j)*nbfac;   i=i+1
        iXI0    = (i+j)*ncell; bXI0    = (i+j)*nbfac;   i=i+1
        iXI2    = (i+j)*ncell; bXI2    = (i+j)*nbfac;   i=i+1
        iALPH0  = (i+j)*ncell; bALPH0  = (i+j)*nbfac;   i=i+1
        iALPH2  = (i+j)*ncell; bALPH2  = (i+j)*nbfac;   i=i+1
        IF (mduc_outSpraySrc .EQ. 1) THEN
          iZ0VPR  = (i+j)*ncell; bZ0VPR  = (i+j)*nbfac;   i=i+1
          iZ2VPR  = (i+j)*ncell; bZ2VPR  = (i+j)*nbfac;   i=i+1
          iZ0VEV  = (i+j)*ncell; bZ0VEV  = (i+j)*nbfac;   i=i+1
          iZ2VEV  = (i+j)*ncell; bZ2VEV  = (i+j)*nbfac;   i=i+1
          iZ0VDI  = (i+j)*ncell; bZ0VDI  = (i+j)*nbfac;   i=i+1
          iZ2VDI  = (i+j)*ncell; bZ2VDI  = (i+j)*nbfac;   i=i+1
        END IF
      END IF
      IF (mduc_soot .GT. 0) THEN
        iSOOT_N    = (i+j)*ncell; i=i+1
        iSOOT_VOL  = (i+j)*ncell; i=i+1
        iSOOT_AREA = (i+j)*ncell; i=i+1
      END IF
      IF (mduc_outSootSrc .GT. 0) THEN
        iSOOT_SRC = (i+j)*ncell;  i=i+nSootMom
      END IF
      IF (mduc_outSootTP .GT. 0) THEN
        iSOOT_TP  = (i+j)*ncell;  i=i+nSootMom
      END IF
      IF (mduc_all_species .EQ. 0) THEN
        iCPMIX  = (i+j)*ncell; bCPMIX  = (i+j)*nbfac; i=i+1
        iDCPDT  = (i+j)*ncell; bDCPDT  = (i+j)*nbfac; i=i+1
        iTCP    = (i+j)*ncell; bTCP    = (i+j)*nbfac; i=i+1
      END IF
c-----
      END SUBROUTINE mduc_set_pasphi_index
c-----------------------------------------------------------------------
