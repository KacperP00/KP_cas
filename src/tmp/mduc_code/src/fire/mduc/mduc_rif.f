c=======================================================================
      SUBROUTINE mduc_check_injection(mat)
c=======================================================================
c     
c     Purpose:
c       function to check whether to start solving reactive flamelet 
c    
c     Variables:
c       input  - mat            material domain
c       update - mduc_inj       injection flag
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,time
      USE comm1,    ONLY : actphi
      USE mduc,     ONLY : mduc_inj,mduc_mix,mduc_tmain,mduc_tpilot,
     x                     iZ1MEAN,iZ2MEAN,Z1_zero,Z2_zero,mducgetndim,
     x                     mduc_nInj,mduc_mem
c-----
      IMPLICIT NONE
c-----
c----- arguments
      INTEGER, INTENT(IN) :: mat
c----- local variables
      INTEGER :: nc, mixture, solution, prev_inj
c-----------------------------------------------------------------------
c     save the current injection
      prev_inj = mduc_inj
c     set default to inert 1D
      mduc_inj = 10
c----- check for existence of each mixture fraction
c     check if mixture fraction 1 is present
      mduc_mix = 0
      checkZ1: DO nc=nsp(mat),nep(mat)
      IF (actphi(nc+iZ1MEAN) > Z1_zero) THEN
        mduc_mix = 1
        EXIT checkZ1
      END IF
      END DO checkZ1
c     check if mixture fraction 2 is present
      IF (mduc_nInj .EQ. 2 .AND. mduc_mix .EQ. 1) THEN
        checkZ2: DO nc=nsp(mat),nep(mat)
        IF (actphi(nc+iZ2MEAN) > Z2_zero) THEN
          mduc_mix = 2
          EXIT checkZ2
        END IF
      END DO checkZ2
      END IF
c----- check status for solving flamelet
      solution = 0
c     check whether to solve pilot injection
      IF (mduc_tpilot >= 0.0) THEN
        IF (time > mduc_tpilot) THEN
          solution = 1
        END IF
      ELSE IF (mduc_tpilot >= -1.0) THEN
        IF (mduc_mix .ge. 1) THEN
          solution = 1
        END IF
      END IF
c     for restart cases with mduc_tmain = -2.0
      IF (solution > 0) THEN
        IF (mduc_tmain < -1.0 .AND. MDUCgetNdim(mduc_mem) > 1) THEN
          mduc_tmain = -1.0
        END IF
      END IF
c     check whether to solve main injection
      IF (mduc_tmain >= 0.0) THEN
        IF (time > mduc_tmain) THEN
          solution = 2
        END IF
      ELSE IF (mduc_tmain >= -1.0) THEN
        IF (mduc_mix .EQ. 2) THEN
          solution = 2
        END IF
      END IF
c----- set the solution flag:
c         AB -> A, number of mixtures, B, solution dimension
      SELECT CASE (solution)
      CASE (0)                  ! inert
        IF (mduc_mix .LT. 2) THEN
          mduc_inj = 10
        ELSE IF (mduc_mix .EQ. 2) THEN
          mduc_inj = 20
        END IF
      CASE (1)                  ! reactive 1D
        IF (mduc_mix .le. 1) THEN
          mduc_inj = 11
        ELSE IF (mduc_mix .EQ. 2) THEN
          mduc_inj = 21
        END IF
      CASE (2)                  ! reactive 2D
        mduc_inj = 22
      END SELECT
c----- reduce the injections across processors
      CALL iglmax(mduc_inj)
      CALL iglmax(mduc_mix)
c-----
      mduc_inj = max(mduc_inj, prev_inj)
c-----
      END SUBROUTINE mduc_check_injection
c-----
c=======================================================================
      SUBROUTINE mduc_inert_1d(mat,Pavg)
c=======================================================================
c     
c     Purpose:
c       update the flamelet based on the current boundary conditions
c    
c     Variables:
c       input  - mat            material domain
c              - Pavg           volume average pressure in domain
c       update - ymfrac         current species mass fraction
c              - ymfrao         previous species mass fraction
c              - ymfrab         boundary species mass fraction
c              - pasphi(HCHEM)  chemical enthalpy
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,time,dt,iampro,numpro
      USE comm1,    ONLY : den,vol,t,hs,actphi,pasphi
      USE cthmod,   ONLY : numscl,ymfrac
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN) :: mat
      REAL(prec), INTENT(IN) :: Pavg
c----- declare local variables
      INTEGER    :: nc,i,j,n1,nSpecies
      REAL(prec) :: Tox,Tfuel,mass,totMass,chiRefMax
      REAL(prec) :: val_0 = 0.0, val_1 = 1.0
      REAL(prec), DIMENSION(:),   ALLOCATABLE :: x1,chi,fT,pdfG
      REAL(prec), DIMENSION(:,:), ALLOCATABLE :: fY
      INTEGER, PARAMETER :: outID = 3, outW = 80
      CHARACTER(len=256) :: text1,text2
c-----------------------------------------------------------------------
c----- allocate/initialise the local memory
      n1       = MDUCgetGridSize(1, mduc_mem)
      nSpecies = MDUCgetNspecs(mduc_mem)
      ALLOCATE(x1(n1),chi(n1),fT(n1),pdfG(n1),fY(n1,nSpecies))
c----- get the flamelet grid
      CALL MDUCgetGrid(x1,mduc_mem)
c----- scalar dissipation rate is zero
      chi = 0.0
      CALL mduc_scalardiss_1d(mat,n1,x1,chi)
c----- check to see whether to start solving reactive flamelet
      IF (mduc_tpilot .LT. -1.0) THEN
c       compute the max chiref
        chiRefMax = 0.0
        DO nc=nsp(mat),nep(mat)
          IF (pasphi(nc+iCHI1) > chiRefMax) THEN
            chiRefMax = pasphi(nc+iCHI1)
          END IF
        END DO
c-----
        IF (iampro < 2) THEN
          WRITE (*,*) 'Max ChiRef is :',  chiRefMax
        END IF
        IF (chiRefMax .gt. 10.0) THEN
          mduc_inj = 1
        END IF
      END IF
c ----- compute the mass averaged temperature
      Tox     = 0.0
      totMass = 0.0
      fY      = 0.0             ! set all species to zero
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)
        DO i=1,numscl
          IF (indexY(i) >= 1) THEN
            fY(1,indexY(i))  = ymfrac(nc+(i-1)*ncell)
          END IF
        END DO
c       compute for averaging
        Tox     = Tox     + mass*t(nc)
        totMass = totMass + mass
      END DO
c----- reduce variables over all processors
      CALL dglsum(Tox)
      CALL dglsum(totMass)
c----- finish normalization
      Tox   = Tox/totMass
      Tfuel = MDUCgetPointT(val_1,val_0,mduc_mem)
c-----
      IF (iampro .EQ. rootpro) THEN
        CALL MDUCupdateInert1D(dt,PAvg,val_0,val_1,Tox,Tfuel,chi,
     x       fY,fT,mduc_mem)
      END IF
c----- make sure that all the processors have the current solution
      IF (numpro > 1) THEN
        CALL dbcast(rootpro-1, fT, n1)
        CALL dbcast(rootpro-1, fY, n1*nSpecies)
      END IF
c----- update the cfd flow field based on the flamelet solution
      CALL mduc_updatecfd_1d(mat,n1,nSpecies,x1,fY,pdfG)
c----- ensure all species at boundaries are updated to cell values
      CALL mduc_update_boundY(mat)
c----- update all the mole fractions based on the current mass fractions
      CALL cthumf(mat)
c----- store the chemical enthalpy term for the temperature update
      CALL mduc_update_h_chem(mat)
c ----- write flow data to file IF specified
      IF (iampro .EQ. rootpro) THEN
        CALL MDUCwriteFlowData(.TRUE.,chi,.FALSE.,chi,.TRUE.,pdfG,
     x       .FALSE.,chi,mduc_mem)
        CALL MDUCwriteScalars(time,pdfG,mduc_mem)
      END IF
c-----
      DEALLOCATE(x1,chi,fT,pdfG,fY)
c----- write a summary
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_header(outW, "Inert Flamelet (1D) Output")
        CALL mduc_out_real(outW,outID,'Oxidizer temperature','K',Tox,
     x       val_0)
        CALL mduc_out_real(outW,outID,'Fuel temperature','K',Tfuel,
     x       val_0)
c
        WRITE(text1, '(A)') repeat ('*', outW-2)
        CALL mduc_out(outW, text1)
      END IF
      END SUBROUTINE mduc_inert_1d
c-----
c
c=======================================================================
      SUBROUTINE mduc_rif_1d(mat,Pavg,Zmin,Zmax,inert)
c=======================================================================
c     
c     Purpose:
c       solve a 1D flamelet and use in flow field
c    
c     Variables:
c       input  - mat            material domain
c              - Pavg           volume average pressure in domain
c              - Zmin           min mixture fraction in domain
c              - Zmax           max mixture fraction in domain
c       update - ymfrac         current species mass fraction
c              - ymfrao         previous species mass fraction (in case)
c              - pasphi(HCHEM)  chemical enthalpy 
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nbfac,nsp,nep,nsb,neb,time,dt,iampro,
     x                     numpro
      USE comm1,    ONLY : den,vol,t,actphi,pasphi
      USE cthmod,   ONLY : numscl,ymfrac
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN) :: mat,inert
      REAL(prec), INTENT(IN) :: Pavg,Zmin,Zmax
c----- declare local variables
      INTEGER    :: nc,nb,i,j
      INTEGER    :: n1,nSpecies
      INTEGER    :: FALSEVAL=0,TRUEVAL=1
      REAL(prec) :: Tox,Tfuel,mass,totMass,chiRefMax
      REAL(prec) :: val_0 = 0.0, val_1 = 1.0
      REAL(prec), DIMENSION(:),   ALLOCATABLE :: x1,chi,fT,pdfG
      REAL(prec), DIMENSION(:,:), ALLOCATABLE :: fY,sootM_Z
      INTEGER, PARAMETER :: outID = 3, outW = 80
      CHARACTER(len=256) :: text1,text2
c-----------------------------------------------------------------------
c----- allocate/initialise the local memory
      n1       = MDUCgetGridSize(1, mduc_mem)
      nSpecies = MDUCgetNspecs(mduc_mem)
      ALLOCATE(x1(n1),chi(n1),fT(n1),pdfG(n1))
      ALLOCATE(fY(n1,nSpecies))
c----- get the flamelet grid
      CALL MDUCgetGrid(x1,mduc_mem)
c----- compute scalar dissipation rate
      chi = 0.0
      CALL mduc_scalardiss_1d(mat,n1,x1,chi)
c---- check oxidizer temperature and max chi      
      chiRefMax = 0.0
      Tox     = 0.0
      totMass = 0.0
      DO nc=nsp(mat),nep(mat)
        mass    = den(nc)*vol(nc)
        Tox     = Tox     + mass*t(nc)
        totMass = totMass + mass
c----- check to see whether to start solving reactive flamelet
        IF (pasphi(nc+iCHI1) > chiRefMax) THEN
          chiRefMax = pasphi(nc+iCHI1)
        END IF
      END DO
c---- reduce variables over all processors
      CALL dglsum(Tox)
      CALL dglsum(totMass)
      CALL dglmax(chiRefMax)
c---- finish normalization
      Tox   = Tox/totMass
      Tfuel = MDUCgetPointT(val_1,val_0,mduc_mem)
c----- if soot model is active, compute and pass conditional soot moments
      IF (mduc_soot .GT. 0) THEN
        ALLOCATE(sootM_Z(n1,nSootMom))
        CALL mduc_soot_cond_1d(mat,n1,x1,sootM_Z)
c       pass conditional moments to flamelet solver
        CALL MDUCsetSootMom(n1,nSootMom,sootM_Z,mduc_mem)
      END IF
c----- CALL the main flamelet solver
      fT = 0.0
      fY = 0.0
      IF (iampro .EQ. rootpro) THEN
        IF (inert .EQ. 1) THEN
          CALL mduc_out_info("updating inert 1D flamelet")
        ELSE
          CALL mduc_out_info("solving 1D flamelet")
        END IF
      END IF
      IF (inert .EQ. 1) THEN
        CALL MDUCsetTempBounds(Tox,Tfuel,val_0,mduc_mem)
      END IF
      CALL MDUCresetODE(mduc_mem)
      CALL MDUCflamelet1D(dt,Pavg,Zmin,Zmax,chi,inert,fY,fT,mduc_mem)
c----- make sure that all the processors have the current solution
      IF (numpro > 1) THEN
        CALL dbcast(rootpro-1, fT, n1)
        CALL dbcast(rootpro-1, fY, n1*nSpecies)
      END IF
c----- update the CFD solution
      CALL mduc_updatecfd_1d(mat,n1,nSpecies,x1,fY,pdfG)
c----- ensure all species at boundaries are updated to cell values
      CALL mduc_update_boundY(mat)
c----- update all the mole fractions based on the current mass fractions
      CALL cthumf(mat)
c----- store the chemical enthalpy term for the temperature update
      CALL mduc_update_h_chem(mat)
c----- write flow data to file if specified
      IF (iampro .EQ. rootpro) THEN
        CALL MDUCwriteFlowData(TRUEVAL,chi,FALSEVAL,chi,TRUEVAL,pdfG,
     x       FALSEVAL,chi,mduc_mem)
        CALL MDUCwriteScalars(time,pdfG,mduc_mem)
      END IF
c-----      
      DEALLOCATE(x1,chi,fT,pdfG,fY)
      IF (ALLOCATED(sootM_Z)) DEALLOCATE(sootM_Z)
      IF (mduc_soot .GT. 0) THEN
c     TODO: compute and store the moments here for once per timestep
      END IF
c----- write a summary
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_header(outW, "MDUC RIF (1D) Summary")
        WRITE(text1, '(A)') repeat (' ', outW-2)
        CALL mduc_out(outW, text1)
c     TODO: finish summary
        WRITE(text2,'(A20)') 'TODO: finish summary '
        CALL mduc_out(outW, text1)
        WRITE(text1, '(A)') repeat ('*', outW-2)
        CALL mduc_out(outW, text1)
      END IF
      END SUBROUTINE mduc_rif_1d
c-----
c=======================================================================
      SUBROUTINE mduc_rif_21d(mat,Pavg,Zmin,Zmax)
c=======================================================================
c     
c     Purpose:
c       solve a 1D flamelet and use in flow field
c        update chi2 in flamelet model after 2D initialization
c    
c     Variables:
c       input  - mat            material domain
c              - Pavg           volume average pressure in domain
c              - Zmin           min mixture fraction in domain
c              - Zmax           max mixture fraction in domain
c       update - ymfrac         current species mass fraction
c              - ymfrao         previous species mass fraction (in case)
c              - pasphi(HCHEM)  chemical enthalpy 
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nbfac,nsp,nep,nsb,neb,time,dt,iampro,
     x                     numpro,itst,dtnw
      USE comm1,    ONLY : den,vol,actphi,pasphi
      USE cthmod,   ONLY : numscl,ymfrac
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN) :: mat
      REAL(prec), INTENT(IN) :: Pavg,Zmin,Zmax
c----- declare local variables
      INTEGER    :: nc,nb,i,j
      INTEGER    :: n1,n2,nSpecies,init2D,nZ1,nZ2,istat
      INTEGER    :: w_sum,FALSEVAL=0,TRUEVAL=1
      REAL(prec) :: mass,totMass
      REAL(prec) :: chi2max
      REAL(prec), DIMENSION(:),   ALLOCATABLE :: x1,x2,chi,fT,pdfG
      REAL(prec), DIMENSION(:,:),   ALLOCATABLE :: chi1,chi12,chi2
      REAL(prec), DIMENSION(:,:), ALLOCATABLE :: fY
      REAL(prec), DIMENSION(:), ALLOCATABLE :: Z1,Z2
      CHARACTER(len=256) :: text1,text2
      CHARACTER(len=25), DIMENSION(2) :: fuelName2
      REAL(prec) :: tstart,tend
c-----------------------------------------------------------------------
c----- allocate/initialise the local memory
      n1       = MDUCgetGridSize(1, mduc_mem)
c     TODO: assume same grid spacing as in x2 and x1 direction
      n2 = n1
      Z2_zero = Z1_zero
      nSpecies = MDUCgetNspecs(mduc_mem)
      ALLOCATE(x1(n1),x2(n2),chi(n1),fT(n1),pdfG(n1))
      ALLOCATE(fY(n1,nSpecies))
      ALLOCATE(chi1(n1,n2),chi12(n1,n2),chi2(n1,n2))
c----- get the flamelet grid
      CALL MDUCgetGrid(x1,mduc_mem)
      DO j=1,n1
        x2(j) = x1(j)
      END DO
      Z2_zero = Z1_zero
c----- compute scalar dissipation rate
      chi = 0.0
      chi1 = 0.0; chi12 = 0.0; chi2 = 0.0
      CALL mduc_scalardiss_2d(mat,n1,n2,x1,x2,chi1,chi12,chi2)
      DO j=1,n1
        chi(j) = chi1(j,1)
      END DO
c----- CALL the main flamelet solver
      fT = 0.0
      fY = 0.0
      CALL MDUCresetODE(mduc_mem)
      CALL MDUCflamelet1D(dt,Pavg,Zmin,Zmax,chi,fY,fT,mduc_mem)
c----- make sure that all the processors have the current solution
      IF (numpro > 1) THEN
        CALL dbcast(rootpro-1, fT, n1)
        CALL dbcast(rootpro-1, fY, n1*nSpecies)
      END IF
c----- update the CFD solution
      CALL mduc_updatecfd_1d(mat,n1,nSpecies,x1,fY,pdfG)
c----- ensure all species at boundaries are updated to cell values
      CALL mduc_update_boundY(mat)
c----- update all the mole fractions based on the current mass fractions
      CALL cthumf(mat)
c----- store the chemical enthalpy term for the temperature update
      CALL mduc_update_h_chem(mat)
c----- write flow data to file if specified
      IF (iampro .EQ. rootpro) THEN
        CALL MDUCwriteFlowData(TRUEVAL,chi,FALSEVAL,chi,TRUEVAL,pdfG,
     x       FALSEVAL,chi,mduc_mem)
        CALL MDUCwriteScalars(time,pdfG,mduc_mem)
      END IF
c-----      
c----- check to see whether to start solving 2D flamelet
      init2D = 0
      IF (mduc_tmain < -1.0) THEN
        chi2max = 0.0
        DO j=1,n2
          DO i=1,n1
            IF (chi2(i,j) .GT. chi2max) THEN
              chi2max = chi2(i,j)
            END IF
          END DO
        END DO
        IF (chi2max > 0.0) THEN
          IF (iampro .LE. 1) THEN
            WRITE(text1,*)'chi2 = ',chi2max,' exceeds specified minimum'
            CALL mduc_out_info(text1)
          END IF
          init2D = 1
        END IF
      ELSE IF (mduc_tmain >= 0.0 .AND. time+dtnw > mduc_tmain) THEN
        init2D=1
      END IF
c
      IF (init2D .EQ. 1 .AND. MDUCgetNdim(mduc_mem) .EQ. 1) THEN
        IF (iampro .LE. 1) THEN
          text1 = 'Initializing a 2D flamelet from a 1D flamelet'
          CALL mduc_out_info(text1)
          WRITE(text1,*) 'f_init nt: ', itst
          CALL mduc_out_info(text1)
        END IF
c     TODO: make 2D initialization more general
c     for now, must be specified in MDUC input file
        nZ2 = -1
        CALL MDUCinit2Dfrom_BC_and_1D(nZ2,-1,2,fuelName2,%VAL(0),
     x       -1.0,mduc_mem)
c     reinitialize the parallel settings
        CALL MDUCsetParallel(MPI_COMM_MDUC, mduc_mem)
c     update the zero point for the grids
        nZ1 = MDUCgetGridSize(1, mduc_mem)
        nZ2 = MDUCgetGridSize(2, mduc_mem)
        ALLOCATE(Z1(nZ1), STAT=istat); Z1 = 0.0
        ALLOCATE(Z2(nZ2), STAT=istat); Z2 = 0.0
        CALL MDUCgetGrid2D(Z1, Z2, mduc_mem)
        Z1_zero = 0.5*(Z1(2)-Z1(1))
        Z2_zero = 0.5*(Z2(2)-Z2(1))
        DEALLOCATE(Z1); DEALLOCATE(Z2)
c     write the initial solution
c     TODO: recompute total enthalpy before doing so
        IF (iampro .LE. 1) THEN
          CALL mduc_out_info('Writing 2D initialisation')
          CALL MDUCwriteFlamelet('mduc_initial_2D.rif',mduc_mem)
          CALL mduc_out_info('Finished flamlet write')
        END IF
        n2 = MDUCgetGridSize(2, mduc_mem)
c     update chi_m1
        CALL MDUCsetScalarDiss2D(n1,n2,chi1,chi12,chi2,mduc_mem)
      END IF
c
      DEALLOCATE(x1,x2,chi,chi1,chi12,chi2,fT,pdfG,fY)
c----- write a summary
      IF (iampro .EQ. rootpro) THEN
        w_sum = 80
        CALL mduc_out_header(w_sum, "MDUC RIF (2-1D) Summary")
        WRITE(text1, '(A)') repeat (' ', w_sum-2)
        CALL mduc_out(w_sum, text1)
c     TODO: finish summary
        WRITE(text2,'(A20)') 'TODO: finish summary '
        CALL mduc_out(w_sum, text1)
        WRITE(text1, '(A)') repeat ('*', w_sum-2)
        CALL mduc_out(w_sum, text1)
      END IF
      END SUBROUTINE mduc_rif_21d
c-----
c=======================================================================
      SUBROUTINE mduc_inert_2d(mat,Pavg)
c=======================================================================
c     
c     Purpose:
c       update the flamelet based on the current boundary conditions
c    
c     Variables:
c       input  - mat            material domain
c              - Pavg           volume average pressure in domain
c       update - ymfrac         current species mass fraction
c              - ymfrao         previous species mass fraction
c              - ymfrab         boundary species mass fraction
c              - pasphi(HCHEM)  chemical enthalpy
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,time,dt,iampro,numpro
      USE comm1,    ONLY : den,vol,t,hs,actphi,pasphi
      USE cthmod,   ONLY : numscl,ymfrac
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN) :: mat
      REAL(prec), INTENT(IN) :: Pavg
c----- declare local variables
      INTEGER    :: nc,i,j,n1,n2,nSpecies,w_sum
      CHARACTER(len=256) :: text1,text2
      REAL(prec) :: Tox,Tfuel1,Tfuel2,mass,totMass,chiRefMax
      REAL(prec) :: val_0 = 0.0, val_1 = 1.0
      REAL(prec), DIMENSION(:),     ALLOCATABLE :: x1,x2
      REAL(prec), DIMENSION(:,:),   ALLOCATABLE :: chi1,chi12,chi2
      REAL(prec), DIMENSION(:,:),   ALLOCATABLE :: pdfG,fT
      REAL(prec), DIMENSION(:,:,:), ALLOCATABLE :: fY
c-----------------------------------------------------------------------
c----- allocate/initialise the local memory
      n1       = MDUCgetGridSize(1, mduc_mem)
      n2       = MDUCgetGridSize(2, mduc_mem)
      nSpecies = MDUCgetNspecs(mduc_mem)
      ALLOCATE(x1(n1),x2(n2))
      ALLOCATE(chi1(n1,n2),chi12(n1,n2),chi2(n1,n2))
      ALLOCATE(fT(n1,n2),pdfG(n1,n2),fY(n1,n2,nSpecies))
c----- get the flamelet grid
      CALL MDUCgetGrid2D(x1,x2,mduc_mem)
c----- compute scalar dissipation rate
      chi1 = 0.0; chi12 = 0.0; chi2 = 0.0
      CALL mduc_scalardiss_2d(mat,n1,n2,x1,x2,chi1,chi12,chi2)
c----- DEBUGGING: write scalar dissipation rate from fortran
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Writing 2D domain scalar dissipation rate')
        OPEN(unit=1,file='chiF.dat')
        WRITE(1,*) '# x1 x2 chi1 chi12 chi2', n1, n2
        DO j=1,n2
          DO i=1,n1
            WRITE(1,'(F9.6,1X,F9.6,ES15.6,ES15.6,ES15.6)')
     x           x1(i), x2(j), chi1(i,j),chi12(i,j),chi2(i,j)
          END DO
          WRITE(1,*) ''
        END DO
        CLOSE(unit=1)
      END IF
c----- compute the mass averaged temperature
      Tox     = 0.0
      totMass = 0.0
      fY      = 0.0
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)
        DO i=1,numscl
          IF (indexY(i) >= 1) THEN
            fY(1,1,indexY(i)) = ymfrac(nc+(i-1)*ncell)
          END IF
        END DO
c----- compute for averaging
        Tox     = Tox     + mass*t(nc)
        totMass = totMass + mass
      END DO
c----- reduce variables over all processors
      CALL dglsum(Tox)
      CALL dglsum(totMass)
c----- finish normalization
      Tox    = Tox/totMass
      Tfuel1 = MDUCgetPointT(val_1,val_0,mduc_mem)
      Tfuel2 = MDUCgetPointT(val_0,val_1,mduc_mem)
c----- update the flamelet solver
      fT = 0.0
      fY = 0.0
      IF (iampro .EQ. rootpro) THEN
        CALL MDUCupdateInert2D(dt,PAvg,val_0,val_1,val_0,val_1,Tox,
     x       Tfuel1,Tfuel2,chi1,chi12,chi2,fY,fT,mduc_mem)
c----- DEBUGGING: write species data from fortran to check consistency
        CALL mduc_out_info('Writing fortran species')
        OPEN(unit=1,file='f_species.dat')
        DO j=1,nSpecies
          WRITE(1,'(ES15.6)') fY(10,5,j)
        END DO
        WRITE(1,*) ''
        CLOSE(unit=1)
      END IF
c----- make sure that all the processors have the current solution
      IF (numpro > 1) THEN
        CALL dbcast(rootpro-1, fT, n1*n2)
        CALL dbcast(rootpro-1, fY, n1*n2*nSpecies)
      END IF
c----- update the cfd flow field based on the flamelet solution
      CALL mduc_updatecfd_2d(mat,n1,n2,nSpecies,x1,x2,fY,fT,pdfG)
c----- DEBUGGING: write fortran pdf for consistency check
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Writing 2D domain PDF')
        OPEN(unit=1,file='pdfF.dat')
        WRITE(1,*) '# x1 x2 pdf ', n1, n2
        DO j=1,n2
          DO i=1,n1
            WRITE(1,'(F9.6,1X,F9.6,ES15.6)') x1(i), x2(j), pdfG(i,j)
          END DO
          WRITE(1,*) ''
        END DO
        CLOSE(unit=1)
      END IF
c----- ensure all species at boundaries are updated to cell values
      CALL mduc_update_boundY(mat)
c----- update all the mole fractions based on the current mass fractions
      CALL cthumf(mat)
c----- store the chemical enthalpy term for the temperature update
      CALL mduc_update_h_chem(mat)
c----- write flow data to file IF specified
      if (iampro .EQ. rootpro) THEN
        CALL MDUCwriteFlowData(.TRUE.,chi1,.FALSE.,chi2,.TRUE.,pdfG,
     x       .FALSE.,chi1,mduc_mem)
        CALL MDUCwriteScalars(time,pdfG,mduc_mem)
      END IF
c-----
      DEALLOCATE(x1,x2)
      DEALLOCATE(chi1,chi12,chi2,fT,pdfG)
      DEALLOCATE(fY)
c----- write a summary (TODO)
      IF (iampro .EQ. rootpro) THEN
        w_sum = 80
        CALL mduc_out_header(w_sum, "MDUC RIF (1D) Summary")
        WRITE(text1, '(A)') repeat (' ', w_sum-2)
        CALL mduc_out(w_sum, text1)
c     TODO: finish summary
        WRITE(text2,'(A20)') 'TODO: finish summary '
        CALL mduc_out(w_sum, text1)
        WRITE(text1, '(A)') repeat ('*', w_sum-2)
        CALL mduc_out(w_sum, text1)
      END IF
      END SUBROUTINE mduc_inert_2d
c-----
c=======================================================================
      SUBROUTINE mduc_rif_2d(mat,Pavg,Z1min,Z1max,Z2min,Z2max)
c=======================================================================
c     
c     Purpose:
c       solve a 2D flamelet and use in flow field
c    
c     Variables:
c       input  - mat            material domain
c              - Pavg           volume average pressure in domain
c              - Z1min,Z2min    min mixture fraction in domain
c              - Z1max,Z2max    max mixture fraction in domain
c       update - ymfrac         current species mass fraction
c              - ymfrao         previous species mass fraction (in case)
c              - pasphi(HCHEM)  chemical enthalpy 
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nbfac,nsp,nep,nsb,neb,time,dt,
     x                     iampro,numpro,itst
      USE comm1,    ONLY : den,vol,actphi,pasphi
      USE cthmod,   ONLY : numscl,ymfrac
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN) :: mat
      REAL(prec), INTENT(IN) :: Pavg,Z1min,Z1max,Z2min,Z2max
c----- declare local variables
      INTEGER    :: nc,nb,i,j,ir,ib
      INTEGER    :: n1,n2,nSpecies,nSolve,pdfFlag
      INTEGER    :: nZ1, counter
      INTEGER    :: w_sum,FALSEVAL=0, TRUEVAL=1
      INTEGER, DIMENSION(2) :: pdfBound
      CHARACTER(len=256) :: text1,text2
      REAL(prec) :: Tox,mass,totMass,meanPDF,wmix,h_chem,totMol
      REAL(prec) :: chiRefMax
      REAL(prec), DIMENSION(:),     ALLOCATABLE :: x1,x2
      REAL(prec), DIMENSION(:,:),   ALLOCATABLE :: chi1,chi12,chi2
      REAL(prec), DIMENSION(:,:),   ALLOCATABLE :: fT,pdfG
      REAL(prec), DIMENSION(:,:,:), ALLOCATABLE :: fY
c-----------------------------------------------------------------------
      IF (iampro .le. 1) THEN
        CALL mduc_out_info('Starting 2D Flamelet coupling')
      END IF
c----- allocate/initialise the local memory
      n1       = MDUCgetGridSize(1, mduc_mem)
      n2       = MDUCgetGridSize(2, mduc_mem)
      nSpecies = MDUCgetNspecs(mduc_mem)
      ALLOCATE(x1(n1),x2(n2))
      ALLOCATE(chi1(n1,n2),chi12(n1,n2),chi2(n1,n2))
      ALLOCATE(fT(n1,n2),pdfG(n1,n2),fY(n1,n2,nSpecies))
c----- get the flamelet grid
      CALL MDUCgetGrid2D(x1,x2,mduc_mem)
c----- compute scalar dissipation rate
      chi1 = 0.0; chi12 = 0.0; chi2 = 0.0
      CALL mduc_scalardiss_2d(mat,n1,n2,x1,x2,chi1,chi12,chi2)
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Writing 2D domain scalar dissipation rate')
        OPEN(unit=1,file='chiF.dat')
        WRITE(1,*) '# x1 x2 chi1 chi12 chi2'
        DO j=1,n2
          DO i=1,n1
            WRITE(1,'(F9.6,1X,F9.6,ES15.6,ES15.6,ES15.6)')
     x           x1(i), x2(j), chi1(i,j),chi12(i,j),chi2(i,j)
          END DO
          WRITE(1,*) ''
        END DO
        CLOSE(unit=1)
      END IF
c----- call the main flamelet solver
      fT = 0.0
      fY = 0.0
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Starting 2D flamelet solve')
      END IF
      CALL MDUCsolve2D(dt,Pavg,Z1min,Z1max,Z2min,Z2max,chi1,chi12,chi2,
     x     fY,fT,mduc_mem)
c----- make sure that all the processors have the current solution
      IF (numpro > 1) THEN
        CALL dbcast(rootpro-1, fT, n1*n2)
        CALL dbcast(rootpro-1, fY, n1*n2*nSpecies)
      END IF
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Finished 2D flamelet solve')
        CALL mduc_out_info('Starting 2D cell update')
      END IF
c----- update the cfd solution
      CALL mduc_updatecfd_2d(mat,n1,n2,nSpecies,x1,x2,fY,fT,pdfG)
c----- DEBUGGING: write fortran pdf for consistency check
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Finished 2D cell update')
        CALL mduc_out_info('Writing 2D domain PDF')
        OPEN(unit=1,file='pdfF.dat')
        WRITE(1,*) '# x1 x2 pdf ', n1, n2
        DO j=1,n2
          DO i=1,n1
            WRITE(1,'(F9.6,1X,F9.6,ES15.6)') x1(i), x2(j), pdfG(i,j)
          END DO
          WRITE(1,*) ''
        END DO
        CLOSE(unit=1)
      END IF
c----- ensure all species at boundaries are updated to cell values
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Updating boundary species')
      END IF
      CALL mduc_update_boundY(mat)
c----- update all the mole fractions based on the current mass fractions
      CALL cthumf(mat)
c----- store the chemical enthalpy term for the temperature update
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Updating chemical enthalpy')
      END IF
      CALL mduc_update_h_chem(mat)
c----- write flow data to file if specified
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Writing flow statistics')
        CALL MDUCwriteFlowData(TRUEVAL,chi1,FALSEVAL,chi2,TRUEVAL,pdfG,
     x       FALSEVAL,chi1,mduc_mem)
        CALL MDUCwriteScalars(time,pdfG,mduc_mem)
      END IF
c-----
      DEALLOCATE(x1,x2)
      DEALLOCATE(chi1,chi12,chi2)
      DEALLOCATE(fT,fY,pdfG)
c----- write a summary (TODO)
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_info('Finished 2D Flamelet coupling')
      END IF
      END SUBROUTINE mduc_rif_2d
c-----
c=======================================================================
      SUBROUTINE mduc_updatecfd_1d(mat,n1,nSpecies,x1,fY,pdfG)
c=======================================================================
c     
c     Purpose:
c       update the flow field using the 1D flamelet results
c    
c     Variables:
c       input  - mat            material domain (from FIRE)
c              - fY             species mass fractions condition on
c                               conditioned on mixture fraction
c              - n1             size of mixture fraction grid
c              - x1             mixture fraction grid
c              - nSpecies       number of species present
c       update - ymfrac         current species mass fraction
c                ymfrao         previous species mass fraction (in case)
c       output - pdfG           global mass averaged pdf
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,iampro,dt
      USE comm1,    ONLY : t,den,vol,actphi,pasphi,paspho,cpcof
      USE cthmod,   ONLY : numscl,ymfrac,ymfrao
      USE mduc,     ONLY : mduc_mem,iZ1MEAN,iZ1VAR,MDUCcalcMeanY_1D,
     x                     indexY,Z1_zero,rootpro,iCPMIX,iTCP,iDCPDT,
     x                     mduc_all_species,iHCHEM
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN)                          :: mat,n1,nSpecies
      REAL(prec), INTENT(IN),  DIMENSION(n1)          :: x1
      REAL(prec), INTENT(IN),  DIMENSION(n1,nSpecies) :: fY
      REAL(prec), INTENT(OUT), DIMENSION(n1)          :: pdfG
c----- declare local variables
      INTEGER    :: i,j,nc,totCell
      INTEGER    :: nZ1
      REAL(prec) :: mass,totMass,cpmix,cpP,cpM,dCpdT,delT
      REAL(prec) :: Zmean,Zmin,Zmax,Zvar
      REAL(prec) :: Y_O2,Y_CO,Y_CO2,Y_H2O,Y_NO,Y_NO2,Y_NOX
      REAL(prec), DIMENSION(nSpecies) :: Ymean,h,cp,cpiP,cpiM
      INTEGER, PARAMETER :: outID = 3, outW = 80
      CHARACTER(len=256) :: text1,text2
c----- pdf variables
      INTEGER                   :: pdfFlag
      INTEGER,    DIMENSION(2)  :: pdfBound
      REAL(prec), DIMENSION(n1) :: pdf
c-----------------------------------------------------------------------
      nZ1     = 0
      totMass = 0.0
      pdfG    = 0.0
      Zmean   = 0.0
      Zvar    = 0.0
      Zmin    =  1e30
      Zmax    = -1e30
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)
c----- compute the mean mass fraction of all species
        Ymean = 0.0
        pdf   = 0.0
        IF (actphi(nc+iZ1MEAN) > Z1_zero) THEN
c         compute the mean species mass fractions assuming distribution
          pdfFlag = MDUCcalcMeanY_1D('beta',n1,x1,nSpecies,
     x         actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),fY,Ymean,
     x         pdf,pdfBound)
c         add contribution to global pdf
          DO j=pdfBound(1),pdfBound(2)
            pdfG(j) = pdfG(j) + mass*pdf(j)
          END DO
          nZ1 = nZ1 + 1  ! count cells with Z1
        ELSE             ! no mixture fraction, delta function oxidizer
          DO i=1,nSpecies
            Ymean(i) = fY(1,i)
          END DO
          pdfG(1) = pdfG(1) + mass
        END IF
c----- update the cell properties
        DO i=1,numscl
c         update mean species
          IF (indexY(i) >= 0) THEN
            ymfrac(nc+(i-1)*ncell) = Ymean(indexY(i))
            ymfrao(nc+(i-1)*ncell) = ymfrac(nc+(i-1)*ncell)
          END IF
        END DO
        IF (mduc_all_species .EQ. 0) THEN
c          update the mixture specific heat
          delT = 10
          CALL MDUCcomputeThermoData(h,cp,  t(nc),     mduc_mem)
          CALL MDUCcomputeThermoData(h,cpiP,t(nc)+delT,mduc_mem)
          CALL MDUCcomputeThermoData(h,cpiM,t(nc)-delT,mduc_mem)
          cpmix = 0.0
          cpP   = 0.0
          cpM   = 0.0
          DO i=1,nSpecies
            cpmix = cpmix + Ymean(i)*cp(i)
            cpP   = cpP   + Ymean(i)*cpiP(i)
            cpM   = cpM   + Ymean(i)*cpiM(i)
          END DO
          pasphi(nc+iTCP)   = t(nc)
          pasphi(nc+iCPMIX) = cpmix
          pasphi(nc+iDCPDT) = (cpP - cpM)/(2*delT)
          cpcof(nc) = cpmix
        END IF
c-----  compute domain averaged quantities     
        Zmean   = Zmean + mass*actphi(nc+iZ1MEAN)
        Zvar    = Zvar  + mass*actphi(nc+iZ1VAR)
        IF (actphi(nc+iZ1MEAN) .lt. Zmin) Zmin = actphi(nc+iZ1MEAN)
        IF (actphi(nc+iZ1MEAN) .gt. Zmax) Zmax = actphi(nc+iZ1MEAN)
        totMass = totMass + mass
      END DO
c----- get the total cells for this processor
      totCell = ncell
c----- reduce variables across processors
      CALL iglsum(totCell)
      CALL iglsum(nZ1)
      CALL dglsum(Zmean)
      CALL dglsum(Zvar)
      CALL dglmin(Zmin)
      CALL dglmax(Zmax)
      CALL dglsum(totMass)
      CALL dglsumvec(pdfG,n1)
      Zmean = Zmean/totMass
      Zvar  = Zvar/totMass
c----- normalise the global pdf
      pdfG = pdfG/totMass
c----- write summary
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_header(outW,'Flamelet (1D) Combustion Output')
        CALL mduc_out_int(outW,outID,'Mixture cells updated',nZ1)
        CALL mduc_out_real(outW,outID,'Mixture percent domain','%',
     x       real(nZ1/totCell*100.0),0)
        CALL mduc_out_real(outW,outID,'Mixture Fraction Mean','-',
     x       Zmean,1)
        CALL mduc_out_real(outW,outID,'Mixture Fraction Variance','-',
     x       Zvar,1)
        CALL mduc_out_real(outW,outID,'Minimum Mixture Fraction','-',
     x       Zmin,1)
        CALL mduc_out_real(outW,outID,'Maximum Mixture Fraction','-',
     x       Zmax,1)
        WRITE(text1, '(A)') repeat ('*', outW-2)
        CALL mduc_out(outW, text1)
      END IF
c-----
      END SUBROUTINE mduc_updatecfd_1d
c-----      
c=======================================================================
      SUBROUTINE mduc_updatecfd_2d(mat,n1,n2,nSpecies,x1,x2,fY,fT,pdfG)
c=======================================================================
c     
c     Purpose:
c       update the flow field using the 2D flamelet results
c    
c     Variables:
c       input  - mat            material domain (from FIRE)
c              - fY             species mass fractions conditioned on
c                               mixture fractions
c              - n1,n2          size of mixture fraction grid
c              - x1,x2          mixture fraction grid
c              - nSpecies       number of species present
c       update - ymfrac         current species mass fraction
c                ymfrao         previous species mass fraction (in case)
c       output - pdfG           global mass averaged pdf
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nbfac,nsp,nep,iampro,dt
      USE comm1,    ONLY : t,den,vol,actphi,pasphi,paspho,cpcof
      USE cthmod,   ONLY : numscl,ymfrac,ymfrao
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN)                             :: mat,n1,n2
      INTEGER,    INTENT(IN)                             :: nSpecies
      REAL(prec), INTENT(IN),  DIMENSION(n1)             :: x1
      REAL(prec), INTENT(IN),  DIMENSION(n2)             :: x2
      REAL(prec), INTENT(IN),  DIMENSION(n1,n2)          :: fT
      REAL(prec), intent(out), DIMENSION(n1,n2)          :: pdfG
      REAL(prec), INTENT(IN),  DIMENSION(n1,n2,nSpecies) :: fY
c----- declare local variables
      INTEGER :: i,j,k,nc
      INTEGER :: nZ1,nZ2,nZ,counter,nCount,lineCount    ! counters
      REAL(prec), DIMENSION(3) :: Zmean,Zmin,Zmax,Zvar
      INTEGER :: totCell
      REAL(prec) :: mass,totMass,cpmix,tmp
      REAL(prec), DIMENSION(nSpecies) :: Ymean,h,cp
      INTEGER, PARAMETER :: outID = 3, outW = 80
      CHARACTER(len=256) :: text1,text2
c----- pdf variables
      INTEGER                            :: pdfFlag
      INTEGER,    DIMENSION(2)           :: pdfBound1D
      REAL(prec), DIMENSION(n1,n2)       :: pdf
      REAL(prec), DIMENSION(n1)          :: pdfX1
      REAL(prec), DIMENSION(n2)          :: pdfX2
      REAL(prec), DIMENSION(n1,nSpecies) :: fY1
      REAL(prec), DIMENSION(n2,nSpecies) :: fY2
c-----------------------------------------------------------------------
      nZ1 = 0; nZ2 = 0; nZ = 0; counter = 0
      pdfG    = 0.0
      totMass = 0.0
      Zmean   = 0.0
      Zvar    = 0.0
      Zmin    = 0.0
      Zmax    = 0.0
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)  ! local cell mass
c----- compute the mean mass fraction of all species
        Ymean = 0.0
        pdf   = 0.0
        IF (actphi(nc+iZ1MEAN) > Z1_zero) THEN
c         check IF both mixture fractions are present
          IF (actphi(nc+iZ2MEAN) > Z2_zero) THEN
c           compute the mean species (use dirichlet for now)
            pdfFlag = MDUCcalcMeanY_2D('dirichlet',n1,x1,n2,x2,nSpecies,
     x           actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),                  
     x           actphi(nc+iZ2MEAN),actphi(nc+iZ2VAR),actphi(nc+iZ0VAR),
     x           fY,Ymean,pdf,1)
c     TODO: could check here for convergence, pdf errors, etc
c           add contribution to global pdf
            DO j=1,n2
              DO i=1,n1
                pdfG(i,j) = pdfG(i,j) + mass*pdf(i,j)
              END DO
            END DO
            nZ = nZ + 1
          ELSE                  ! cell only has Z1
c     compute the mean species
c     NOTE: not sure if need a local fY yet
            fY1 = fY(:,1,:)
            pdfFlag = MDUCcalcMeanY_1D('beta',n1,x1,nSpecies,
     x           actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),fY1,Ymean,
     x           pdfX1,pdfBound1D)
c           store for the global pdf
            DO i=pdfBound1D(1),pdfBound1D(2)
              pdfG(i,1) = pdfG(i,1) + mass*pdfX1(i)
            END DO
            nZ1 = nZ1 + 1       ! count cells with Z1
          END IF
        ELSE IF (actphi(nc+iZ2MEAN) > Z2_zero) THEN ! cell only has Z2
c----- compute the mean species
          fY2 = fY(1,:,:)
          pdfFlag = MDUCcalcMeanY_1D('beta',n2,x2,nSpecies,
     x         actphi(nc+iZ2MEAN),actphi(nc+iZ2VAR),fY2,Ymean,
     x         pdfX2,pdfBound1D)
c         store for the global pdf
          DO j=pdfBound1D(1),pdfBound1D(2)
            pdfG(1,j) = pdfG(1,j) + mass*pdfX2(j)
          END DO
          nZ2 = nZ2 + 1 ! count cells with Z1
        ELSE            ! no mixture fraction, delta function oxidizer
          DO i=1,nSpecies
            Ymean(i) = fY(1,1,i)
          END DO
          pdfG(1,1) = pdfG(1,1) + mass
        END IF
c----- update the cell properties
        DO i=1,numscl
          IF (indexY(i) >= 0) THEN
            ymfrac(nc+(i-1)*ncell) = Ymean(indexY(i))
            ymfrao(nc+(i-1)*ncell) = ymfrac(nc+(i-1)*ncell)
          END IF
        END DO
c       update the mixture specific heat
        IF (mduc_all_species .EQ. 0) THEN
          CALL MDUCcomputeThermoData(h,cp,t(nc),mduc_mem)
          cpmix = 0.0
          DO i=1,nSpecies
            cpmix = cpmix + Ymean(i)*cp(i)
          END DO
          cpcof(nc) = cpmix
        END IF
c-----  compute domain average quantities
        Zmean(1) = Zmean(1) + mass*actphi(nc+iZ1MEAN)
        Zvar(1)  = Zvar(1)  + mass*actphi(nc+iZ1VAR)
        Zmean(2) = Zmean(2) + mass*actphi(nc+iZ2MEAN)
        Zvar(2)  = Zvar(2)  + mass*actphi(nc+iZ2VAR)
        Zmean(3) = Zmean(3) +
     x       mass*(1-actphi(nc+iZ1MEAN)-actphi(nc+iZ2MEAN))
        IF (actphi(nc+iZ1MEAN).LT.Zmin(1)) Zmin(1) = actphi(nc+iZ1MEAN)
        IF (actphi(nc+iZ1MEAN).GT.Zmax(1)) Zmax(1) = actphi(nc+iZ1MEAN)
        IF (actphi(nc+iZ2MEAN).LT.Zmin(2)) Zmin(2) = actphi(nc+iZ2MEAN)
        IF (actphi(nc+iZ2MEAN).GT.Zmax(2)) Zmax(2) = actphi(nc+iZ2MEAN)
        IF (1.0-actphi(nc+iZ1MEAN)-actphi(nc+iZ2MEAN).LT.Zmin(3)) THEN
          Zmin(3) = 1.0-actphi(nc+iZ1MEAN)-actphi(nc+iZ2MEAN)
        END IF
        IF (1.0-actphi(nc+iZ1MEAN)-actphi(nc+iZ2MEAN).GT.Zmax(2)) THEN 
          Zmax(3) = 1.0-actphi(nc+iZ1MEAN)-actphi(nc+iZ2MEAN)
        END IF
        totMass = totMass + mass
      END DO
      totCell = ncell
c----- reduce variables across processors
      CALL dglsumvec(pdfG,n1*n2)
      CALL dglsum(totMass)
      CALL iglsum(nZ1)
      CALL iglsum(nZ2)
      CALL iglsum(nZ)
      CALL dglsumvec(Zmean,3)
      CALL dglsumvec(Zvar,3)
      CALL dglminvec(Zmin,3)
      CALL dglmaxvec(Zmax,3)
      CALL iglsum(totCell)
      DO i=1,3
        Zmean(i) = Zmean(i)/totMass
        Zvar(i)  = Zvar(i)/totMass
      END DO
c----- normalise the global pdf
      pdfG = pdfG/totMass
c---- write a summary
      tmp = real(100.0/totCell)
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_header(outW,'Flamelet (2D) Combustion Output')
        CALL mduc_out_int(outW,outID,'Total Mixture cells updated',
     x       nZ+nZ1+nZ2)
        CALL mduc_out_int(outW,outID,'Mixture 1 cells updated',nZ1)
        CALL mduc_out_int(outW,outID,'Mixture 1 cells updated',nZ2)
        CALL mduc_out_int(outW,outID,'Mixture 1+2 cells updated',
     x       nZ1+nZ2)
        CALL mduc_out_real(outW,outID,'Mixture percent domain','%',
     x       (nZ+nZ1+nZ2)*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture 1 percent domain','%',
     x       nZ1*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture 2 percent domain','%',
     x       nZ2*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture 1+2 percent domain','%',
     x       (nZ1+nZ2)*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture Fraction 1 Mean','-',
     x       Zmean(1),1)
        CALL mduc_out_real(outW,outID,'Mixture Fraction 1 Variance','-',
     x       Zvar(1),1)
        CALL mduc_out_real(outW,outID,'Minimum Mixture 1 Fraction','-',
     x       Zmin(1),1)
        CALL mduc_out_real(outW,outID,'Maximum Mixture 1 Fraction','-',
     x       Zmax(1),1)
        CALL mduc_out_real(outW,outID,'Mixture Fraction 2 Mean','-',
     x       Zmean(2),1)
        CALL mduc_out_real(outW,outID,'Mixture Fraction 2 Variance','-',
     x       Zvar(2),1)
        CALL mduc_out_real(outW,outID,'Minimum Mixture 2 Fraction','-',
     x       Zmin(2),1)
        CALL mduc_out_real(outW,outID,'Maximum Mixture 2 Fraction','-',
     x       Zmax(2),1)
        WRITE(text1, '(A)') repeat ('*', outW-2)
        CALL mduc_out(outW, text1)
      END IF
      END SUBROUTINE mduc_updatecfd_2d
c-----
c=======================================================================
      SUBROUTINE mduc_update_h_chem(mat)
c=======================================================================
c     
c     Purpose:
c       stores the chemical enthalpy for the current composition
c    
c     Variables:
c       input  - mat            material domain
c       udpate - pasphi(iHCHEM) current chemical enthalpy
c              - paspho(iHCHEM) previous chemical enthalpy
c       output - N/A
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nbfac,nsp,nep,nsb,neb,nreg,dt,iampro
      USE comm1,    ONLY : nsr,ner,lb,ibc,pasphi,pasphb,paspho,vol,den
      USE cthmod,   ONLY : numscl,ymfrac,ymfrab,ymolfr,ymolfb,enthfo,
     x                     wmolar
      USE mduc,     ONLY : mduc_mem,iHCHEM,bHCHEM,H_Tref,H_method,
     x                     indexY,MDUCgetNspecies,rootpro,iHRATE
c-----
      IMPLICIT NONE
c-----
c----- function variables
      INTEGER, INTENT(IN) :: mat
c----- local variables
      INTEGER    :: i,nc,nb,ir,ib,nSpecies
      REAL(prec) :: h_chem,wmix,Qc,qL
      REAL(prec), DIMENSION(:), ALLOCATABLE :: hr, cpr
      INTEGER, PARAMETER :: outID = 4, outW = 80
      CHARACTER(len=256) :: text1
c----- allocate working memory
      nSpecies = MDUCgetNspecies(mduc_mem)
      ALLOCATE(hr(nSpecies), cpr(nSpecies))
c----- store chemical enthalpy term for the temperature update
      Qc = 0.0
      IF (H_method .EQ. 0) THEN
        DO nc=nsp(mat),nep(mat) ! cell update
          h_chem = 0.0
          wmix   = 0.0
          DO i=1,numscl
            h_chem = h_chem + ymolfr(nc+(i-1)*ncell)*enthfo(i) ! kJ/kmol
            wmix   = wmix   + ymolfr(nc+(i-1)*ncell)*wmolar(i) ! kg/kmol
          END DO
          h_chem = (h_chem/wmix)*1.d3 ! J/kg
          qL = den(nc)*vol(nc)*(pasphi(nc+iHCHEM) - h_chem)
          Qc = Qc + qL
          pasphi(nc+iHCHEM) = h_chem
          pasphi(nc+iHRATE) = qL/dt
        END DO
        DO ir=0,nreg ! boundary update
          IF (ibc(2,ir) == mat) THEN
            DO ib=nsr(ir),ner(ir)
              h_chem = 0.0
              wmix    = 0.0
              DO i=1,numscl
                h_chem = h_chem + ymolfb(ib+(i-1)*nbfac)*enthfo(i)
                wmix   = wmix   + ymolfb(ib+(i-1)*nbfac)*wmolar(i)
              END DO
              pasphb(ib+bHCHEM) = (h_chem/wmix)*1.d3 ! J/kg
            END DO
          END IF
        END DO
      ELSE IF (H_method .GE. 1) THEN
c----- compute the enthalpies of formation (should store this somewhere)
        CALL MDUCcomputeThermoData(hr,cpr,H_Tref,mduc_mem)
c----- store the chemical enthalpy term for the temperature update
        DO nc=nsp(mat),nep(mat) ! cell update
          h_chem = 0.0
          DO i=1,numscl
            h_chem = h_chem + ymfrac(nc+(i-1)*ncell)*hr(indexY(i))
          END DO
          qL = den(nc)*vol(nc)*(pasphi(nc+iHCHEM) - h_chem)
          Qc = Qc + qL
          pasphi(nc+iHCHEM) = h_chem
          pasphi(nc+iHRATE) = qL/dt
        END DO
        DO ir=0,nreg ! boundary update
          IF (ibc(2,ir) == mat) THEN
            DO ib=nsr(ir),ner(ir)
              h_chem = 0.0
              DO i=1,numscl
                h_chem = h_chem + ymfrab(ib+(i-1)*nbfac)*hr(indexY(i))
              END DO
              pasphb(ib+bHCHEM) = h_chem
            END DO
          END IF
        END DO
      END IF
c-----      
      CALL dglsum(Qc)
      DEALLOCATE(hr,cpr)
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_header(outW,'Chemical Energy Update')
        CALL mduc_out_real(outW,outID,'Heat Release','J',Qc,1)
        CALL mduc_out_real(outW,outID,'Heat Release Rate','J/s',Qc/dt,1)
        WRITE(text1, '(A)') repeat ('*', outW-2)
        CALL mduc_out(outW,text1)
      END IF
c-----
      END SUBROUTINE mduc_update_h_chem
c-----
c=======================================================================
      SUBROUTINE mduc_update_boundY(mat)
c=======================================================================
c 
c     Purpose:
c       ensure that boundary mass fractions have cell values
c    
c     Variables:
c       input  - mat            material domain
c       udpate - ymfrab         species mass fraction at boundary
c       output - N/A
c                 
c=======================================================================
      USE comm0,    ONLY : ncell,nbfac,nreg
      USE comm1,    ONLY : nsr,ner,lb,ibc,cpcof,cpcofb
      USE cthmod,   ONLY : numscl,ymfrac,ymfrab
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER, INTENT(IN) :: mat
c----- declare local variables
      INTEGER i,nc,ib,ir
c-----------------------------------------------------------------------
c----- loop over all boundary cells
      DO ir=0,nreg
        IF (ibc(2,ir) == mat) THEN
          DO ib=nsr(ir),ner(ir)
            nc = lb(ib)
            DO i=1,numscl
              ymfrab(ib+(i-1)*nbfac) = ymfrac(nc+(i-1)*ncell)
            END DO
            cpcofb(ib) = cpcof(nc)
          END DO
        END IF
      END DO
      END SUBROUTINE mduc_update_boundY
c-----
c=======================================================================
      SUBROUTINE mduc_out(w,text)
c=======================================================================
c 
c     Purpose:
c       find the grid point adjacent to a specified value
c
c     Variables:
c       input  - w              width of output
c              - npad           number of padding characters
c              - cpad           character to use as padding
c              - text           formatted string to be output
c       output - stderr         writes text to stderr
c             
c=======================================================================
      IMPLICIT NONE
c-----
c----- declare function variables
      INTEGER,          INTENT(IN) :: w
      CHARACTER(len=*), INTENT(IN) :: text
c----- declare local variables
      INTEGER :: len_text,npre,npost
c-----------------------------------------------------------------------
c----- get the length of the text
      len_text = len(trim(text))
c----- calculate the padding
      npre  = (132-w)/2
      npost = w - (len_text+2)
      WRITE(*,'(A,A,A,A,A)') repeat(' ',npre),'*', trim(text),
     x     repeat(' ',npost),'*'
c-----
      END SUBROUTINE mduc_out
c-----
c=======================================================================
      SUBROUTINE mduc_out_int(w,id,name,val)
c=======================================================================
c 
c     Purpose:
c       find the grid point adjacent to a specified value
c
c     Variables:
c       input  - w              width of output
c              - id             code section ident
c              - name           variable name
c              - val            variable value
c       output - stderr         writes text to stderr
c             
c=======================================================================
      IMPLICIT NONE
c-----
c----- declare function variables
      INTEGER,          INTENT(IN) :: w,id,val
      CHARACTER(len=*), INTENT(IN) :: name
c----- declare local variables
      INTEGER :: len_text,npre,nmid
      CHARACTER(len=64) :: forma
c-----------------------------------------------------------------------
c----- get the length of the text
      len_text = len(trim(name))+24
c----- calculate the padding
      npre  = (132-w)/2
      nmid  = w - (len_text+2)
      forma = '(2A,I0.2,2X,3A,I12,1X,A)'
      WRITE(*,forma) repeat(' ',npre),'*  %M',id,trim(name),
     x     repeat(' ',nmid),' = ',val,'*'
c-----
      END SUBROUTINE mduc_out_int
c-----
c=======================================================================
      SUBROUTINE mduc_out_real(w,id,name,unit,val,ex)
c=======================================================================
c 
c     Purpose:
c       find the grid point adjacent to a specified value
c
c     Variables:
c       input  - w              width of output
c              - id             code section ident
c              - name           variable name
c              - unit           variable units
c              - val            variable value
c              - ex             exponential format?
c       output - stderr         writes text to stderr
c             
c=======================================================================
      USE prec_mod, ONLY : prec
      IMPLICIT NONE
c-----
c----- declare function variables
      INTEGER,          INTENT(IN) :: w,id,ex
      CHARACTER(len=*), INTENT(IN) :: name,unit
      REAL(prec),       INTENT(IN) :: val
c----- declare local variables
      INTEGER :: len_text,npre,nmid
      CHARACTER(len=64) :: forma,text
c-----------------------------------------------------------------------
c----- get the length of the text
      text = trim(name)//' ['//trim(unit)//']'
      len_text = len(trim(text))+24
c----- calculate the padding
      npre  = (132-w)/2
      nmid  = w - (len_text+2)
      IF (ex .eq. 1) THEN
        forma = '(2A,I0.2,2X,3A,ES12.4,1X,A)'
      ELSE
        forma = '(2A,I0.2,2X,3A,F12.2,1X,A)'
      END IF
      WRITE(*,forma) repeat(' ',npre),'*  %M',id,trim(text),
     x     repeat(' ',nmid),' = ',val,'*'
c-----
      END SUBROUTINE mduc_out_real
c-----
c=======================================================================
      SUBROUTINE mduc_out_header(w,text)
c=======================================================================
c 
c     Purpose:
c       find the grid point adjacent to a specified value
c
c     Variables:
c       input  - w              width of output
c              - text           formatted string to be output
c       output - stderr         writes text to stderr
c             
c=======================================================================
      IMPLICIT NONE
c-----
c----- declare function variables
      INTEGER,          INTENT(IN) :: w
      CHARACTER(len=*), INTENT(IN) :: text
c----- declare local variables
      INTEGER :: len_text, npre, npad, pad1, pad2
c-----
c----- get the length of the text
      len_text = len(trim(text))
c----- calculate the padding
      npre = (132-w)/2
      npad = w - (len_text+2)
      pad1 = npad/2+1
      pad2 = npad - pad1
c-----      
      WRITE(*,*) ''
      WRITE(*,'(2A)') repeat(' ',npre), repeat('*', w)
      WRITE(*,'(6A)') repeat(' ',npre),'*',
     x     repeat(' ',pad1),trim(text),repeat(' ',pad2),'*'
      WRITE(*,'(2A)') repeat(' ',npre), repeat('*', w)
      END SUBROUTINE mduc_out_header
c-----
c=======================================================================
      SUBROUTINE mduc_out_info(text)
c=======================================================================
c 
c     Purpose:
c       find the grid point adjacent to a specified value
c
c     Variables:
c       input  - text           formatted string to be output
c       output - stderr         writes text to stderr
c             
c=======================================================================
      IMPLICIT NONE
c-----
c----- declare function variables
      CHARACTER(len=*), INTENT(IN) :: text
c----- declare local variables
      INTEGER :: w,len_text,npre,npad
c-----
c----- get the length of the text
      len_text = len(trim(text))
      IF (len_text+12+14 .LT. 80) THEN
        w = 80
      ELSE
        w = len_text + 12 + 14
      END IF
      npad = w - (len_text+12+14)
c----- calculate the padding
      npre = (132-w)/2
c-----      
      WRITE(*,'(2A,2X,3A,2X,A)') repeat(' ',npre),repeat('=',5),
     x     'MDUC call : ',trim(text),repeat(' ',npad),repeat('=',5)
      END SUBROUTINE mduc_out_info
c-----
c=======================================================================
      function mduc_grid_locate(n,x,val) result (indx)
c=======================================================================
c 
c     Purpose:
c       find the grid point adjacent to a specified value
c
c     Variables:
c       input  - n              number of grid points
c              - x              grid points
c             - val            value to search for in grid
c       output - indx           grid point index next to value
c             
c=======================================================================
      USE prec_mod, ONLY : prec
      IMPLICIT NONE
c-----
c----- function variables
      INTEGER,    INTENT(IN)               :: n   ! size of grid   
      REAL(prec), INTENT(IN)               :: val ! value to locate
      REAL(prec), INTENT(IN), DIMENSION(n) :: x   ! grid  
c----- local variables
      INTEGER :: indx           ! return variable
      INTEGER :: jL, jR, jM     ! left, right, and mid
      LOGICAL :: ascnd
c-----------------------------------------------------------------------
      jL = 1
      jR = n
      IF (x(n) .ge. x(1)) THEN
        ascnd = .TRUE.
      ELSE
        ascnd = .FALSE.
      END IF
c----- use bisection to find nearest grid point
      DO WHILE (jR - jL > 1)
c     compute a midpoint
        jM = (jL+jR)/2;
        IF ((val .ge. x(jM)) .AND. (ascnd) ) THEN
          jL = jM;
        ELSE
          jR = jM;
        END IF
      END DO
      indx = max(1, min(n-1, jL))
      RETURN
c-----
      END FUNCTION mduc_grid_locate
c-----
c=======================================================================
