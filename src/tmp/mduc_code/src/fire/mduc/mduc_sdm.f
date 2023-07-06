c=======================================================================
      MODULE mduc_sdm
c=======================================================================
c
c     This module computes a source term from the spray evaporation
c     for the scalar variance equation.  It is based on the Single
c     Droplet Model (SDM) developed in:
c
c       REVEILLON, J. and VERVISCH, L., Spray Vaporization in
c       Nonpremixed Turbulent Combustion Modeling: A Single Droplet
c       Model, COMBUST. FLAME 121:75-90 (2000)
c
c     Required modules:
c       usemduc.f, useroot.f
c     Dependent modules:
c       useini.f, usepsc.f
c
c.....contact emdoran@stanford.edu
c-----
c=======================================================================
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : iampro,numpro
      USE usedefpar
      USE mduc
      IMPLICIT NONE
c-----
      INTEGER    :: nNIST       ! number of NIST data points
      INTEGER    :: nZsdm       ! local grid points
      REAL(prec) :: Yf_inf      ! ambient fuel mass fraction            [-]
      REAL(prec) :: T_inf       ! ambient temperature                   [K]
      REAL(prec) :: P_inf       ! ambient pressure                      [Pa]
      REAL(prec) :: Psat        ! saturation pressure                   [Pa]
      REAL(prec) :: Zsat        ! saturation mixture fraction           [-]
      REAL(prec) :: hfg         ! enthalpy of vaporisation              [J/kg]
      REAL(prec) :: cpL, cpG    ! specific heat of liquid/gas           [J/kg]
      REAL(prec) :: sdmW        ! total fuel source term in cell        [kg/s]
      REAL(prec) :: alpha       ! constant in ersatz function           [?]
      REAL(prec) :: Zdrop       ! mixture fraction of liquid drop       [-]
      REAL(prec) :: T_liq_0     ! initial liquid temperature            [K]
      REAL(prec) :: Mfuel, Mair ! molecular weights of fuel/air         [g/mol]
      REAL(prec) :: massSDM     ! cell mass                             [kg]
      REAL(prec) :: tau         ! turbulent eddy turnover time          [s]
      REAL(prec) :: tau_tauZ    ! turblent/scalar time scale ratio      [-]
      REAL(prec) :: Re_t        ! cell turbulent Reynolds number        [-]
      REAL(prec) :: Sc_t = 0.7  ! turbulent Schmidt number              [-]
      REAL(prec) :: C_x  = 2.0  ! timescale proportionality constant    [-]
c----- storage for local arrays
      REAL(prec), DIMENSION(:), allocatable :: eta, pdf
c----- storage to read in NIST DATA
      REAL(prec),DIMENSION(:),allocatable :: NIST_T, NIST_P
      REAL(prec),DIMENSION(:),allocatable :: NIST_CPL,NIST_CPG,NIST_HFG
c-----some constants to save time
      REAL(prec) :: M_PI        = 3.141592653589793
      REAL(prec) :: one_third   = 0.333333333333333
      REAL(prec) :: six_over_pi = 1.909859317102744
      INTEGER    :: out_w
      CHARACTER(len=256) :: sdm_txt1, sdm_txt2
c-----
      CONTAINS
c-----
c-----------------------------------------------------------------------
c     intialise all the single droplet model parameters
      SUBROUTINE sdmInit()
c-----
      INTEGER    :: j,i1
      REAL(prec) :: dZ
c----- get the location of the NIST file (otherwise specified in csf)
      i1 = UDP_GET_INDEX('SaturationData')
      IF (i1 > 0 .AND. i1 <= n_udp) THEN
        mduc_satData = udp(i1)%c_value(1:udp(i1)%i_lstr_value)
      END IF
c----- read in the NIST thermo file 
      CALL readNISTthermo(mduc_satData)
c----- set the fuel liquid properties
c     TODO: get these properties from the solver
      T_liq_0 = 320.0
      Mfuel   = 100.204
      Mair    = 28.9
c----- set the mixture fraction grid
      nZsdm = 201
      allocate(eta(nZsdm), pdf(nZsdm))
      dZ = 1.0/(nZsdm-1)
      eta(1) = 0.0
      DO j=2,nZsdm
        eta(j) = eta(j-1) + dZ
      END DO
c----- write a summary
      out_w = 80
      IF (iampro < 2) THEN
        WRITE (*,*) ''
          CALL mduc_out_header(out_w, "Single Droplet Model Summary")

          WRITE(sdm_txt1, '(A)') REPEAT(' ', out_w-2)
          CALL mduc_out(out_w, sdm_txt1)

c          WRITE(sdm_txt2,'(A20, A)') '  Saturation Data : ', NIST_File
c          CALL mduc_out(out_w, sdm_txt2)
          WRITE(sdm_txt2,'(A20, F12.2, A)') '  Liquid temperature : ',
     x         T_liq_0, ' [K]'
          CALL mduc_out(out_w, sdm_txt2)
          WRITE(sdm_txt2,'(A20, F12.2, A)') '  Fuel molar mass : ',
     x         Mfuel, ' [g/mol]'
          CALL mduc_out(out_w, sdm_txt2)
          CALL mduc_out(out_w, sdm_txt1)

          WRITE(sdm_txt1, '(A)') REPEAT('*', out_w-2)
          CALL mduc_out(out_w, sdm_txt1)
          WRITE (*,*) ''
      END IF
c-----
      END SUBROUTINE sdmInit
c-----
c-----------------------------------------------------------------------
c----- rhs function for Spalding number calculation
c-----------------------------------------------------------------------
      FUNCTION sdmCalcSpaldingFunc(Ts) result (res)
c-----
      REAL(prec), INTENT(in) :: Ts
      REAL(prec) :: Xfs, Yfs    ! fuel mole/mass fraction at droplet surface
      REAL(prec) :: res, BT, Bf, Q
c-----
c----- get the saturation properties for the given surface temperature
      CALL interpNISTthermo(hfg, cpG, cpL, Psat, Ts)
c----- compute mass fraction at the surface based on saturation pressure
      Xfs = Psat/P_inf
      Yfs = Mfuel*Xfs/(Mfuel*Xfs+(1-Xfs)*Mair)
c----- compute the heat transfer
      Q = hfg + cpL*(Ts - T_liq_0)
c----- compute both species and temperature transfer coefficients
      BT = cpG*(T_inf - Ts)/Q
      Bf = (Yfs - Yf_inf)/(1.0 - Yfs)
c----- return the residual
      res = BT - Bf
      RETURN
      END FUNCTION sdmCalcSpaldingFunc
c-----
c-----------------------------------------------------------------------
c----- rhs function for the conditional source term
c        input  : xi  ..... current estimate of exponent [-]
c        output : res ..... residual of current source term estimate
c                            and target source, sdmW
c        module : pdf, eta, alpha, Zdrop, Zsat, sdmW
c                 
c-----------------------------------------------------------------------
      FUNCTION sdmCalcMeanSrc(xi) result (res)
c-----
      REAL(prec), INTENT(in) :: xi
      INTEGER    :: j
      REAL(prec) :: W, res
c-----
      W = 0.0
      DO j=2,nZsdm
        W = W + pdf(j)*alpha*(eta(j)*Zsat/Zdrop)**xi
      END DO
c----- return the residual
      res = W - sdmW
      RETURN
      END FUNCTION sdmCalcMeanSrc
c-----
c-----------------------------------------------------------------------
c----- read in the NIST thermo table
c        input  : fileName  ..... data file to read from
c        output : NIST_T    ..... temperature array              [K]
c                 NIST_P    ..... saturation pressure array      [Pa]
c                 NIST_CPG  ..... specific heat of vapour array  [J/kg]
c                 NIST_CPL  ..... specific heat of liquid array  [J/kg]
c                 NIST_HFG  ..... enthalpy of vaporisation array [J/kg]
c-----------------------------------------------------------------------
      SUBROUTINE readNISTthermo(fileName)
c-----
      CHARACTER(len=512), INTENT(in) :: fileName
      CHARACTER(len=800) :: buffer
      INTEGER    :: ifile, iostat, EOF
      INTEGER    :: i,j,nLines,nVars
      INTEGER    :: sT, sPsat,sCpL,sCpG,sHG,sHL
      REAL(prec) :: HL, HG, dT
      REAL(prec), DIMENSION(:), allocatable :: rowVars
c----- standard NIST file dump
      nVars = 25;
      sT = 1; sPsat = 2; sHL = 6; sCpL = 9; sHG = 18; sCpG = 21;
      ALLOCATE(rowVars(nVars))
c----- get the number of data lines in the file
      ifile = 10
      open(ifile,file=trim(fileName),form='formatted',action='read')
c-----
      nLines = 0
      EOF = 0
      DO while(EOF.eq.0)
        read(ifile, '(A800)', iostat=EOF) buffer 
c        print*, buffer
        IF (len(buffer) > 0) THEN
          nLines = nLines + 1
        END IF
      END DO
      nLines = nLines-1
      nNIST = nLines-1
c----- rewind the file
      REWIND(ifile)
c----- allocate the memory
      ALLOCATE(NIST_T(nNIST), NIST_P(nNIST))
      ALLOCATE(NIST_CPL(nNIST),NIST_CPG(nNIST),NIST_HFG(nNIST))
c----- read in the header
      READ(ifile, '(A800)', iostat=EOF) buffer
c----- read in the rest of the data
      DO i=1,nNIST
        READ(ifile,*,iostat=EOF) (rowVars(j), j=1,nVars)
        NIST_T(i) = rowVars(sT)
        NIST_P(i) = rowVars(sPsat)*1.0E6
        NIST_CPG(i) = rowVars(sCpG)
        NIST_CPL(i) = rowVars(sCpL)
        NIST_HFG(i) = (rowVars(sHG)-rowVars(sHL))*1.0e3
      END DO
c-----
      close(ifile)
c-----
      DEALLOCATE(rowVars)
      END SUBROUTINE readNISTthermo
c-----
c-----------------------------------------------------------------------
c----- interpolate data from the NIST table
c        input  : T    ..... droplet temperature      [K]
c        output : hfg  ..... enthalpy of vaporisation [J/kg]
c                 cpG  ..... specific heat of vapour  [J/kg]
c                 cpL  ..... specific heat of liquid  [J/kg]
c                 Psat ..... saturation pressure      [Pa]
c-----------------------------------------------------------------------
      SUBROUTINE interpNISTthermo(hfg, cpG, cpL, Psat, T)
c-----
      REAL(prec), INTENT(out) :: hfg, cpG, cpL, Psat
      REAL(prec), INTENT(in)  :: T
      INTEGER    :: jL, jR
      REAL(prec) :: C
c-----
      IF ( T < NIST_T(1) .OR. T > NIST_T(nNIST) ) THEN
c     should be an error here
      END IF
c----- find the data point
      jR = 2
      DO WHILE (NIST_T(jR) < T .AND. jR < nNIST)
        jR = jR + 1
      END DO
      jL = MAX(1, jR - 1)
c----- compute the weighting coefficient and interpolate values
      C    = (T - NIST_T(jL))/(NIST_T(jR) - NIST_T(jL))
      hfg  = NIST_HFG(jL) + C*(NIST_HFG(jR)  - NIST_HFG(jL))
      cpG  = NIST_CPG(jL) + C*(NIST_CPG(jR)  - NIST_CPG(jL))
      cpL  = NIST_CPL(jL) + C*(NIST_CPL(jR)  - NIST_CPL(jL))
      Psat = NIST_P(jL)   + C*(NIST_P(jR)    - NIST_P(jL))
      END SUBROUTINE interpNISTthermo
c-----
c-----------------------------------------------------------------------
c----- function to compute the Spalding number
c        input  : Yf    ..... ambient mass fraction of fuel   [-]
c                 Yox   ..... ambient mass fraction of oxygen [-]
c                 T_amb ..... ambient temperature             [K]
c                 P_amb ..... ambient pressure                [Pa]
c        output : B     ..... Spalding number                 [-]
c-----------------------------------------------------------------------
      FUNCTION sdm_calcSpaldingNo(Yf,Yox,T_amb,P_amb) result (B)
c-----
      USE rootSolve, ONLY : rootBrent
      IMPLICIT NONE
c-----
      REAL(prec), INTENT(in) :: Yf, Yox, T_amb, P_amb
      REAL(prec) :: B           ! Spalding number 
      REAL(prec) :: Ts          ! droplet surface temperature
      REAL(prec) :: Yfs         ! droplet surface fuel mass fraction
      REAL(prec) :: Q           ! heat transfer
      REAL(prec) :: HI, LO, brent_tol=1.0e-4
      REAL(prec) :: nu_f = 3.6723703644565089 ! should be generalized
c----- set the module variables
      Yf_inf = Yf
      T_inf  = T_amb
      P_inf  = P_amb
c----- check whether a root for this temperature exists
      LO = sdmCalcSpaldingFunc(NIST_T(1))
      HI = sdmCalcSpaldingFunc(NIST_T(nNIST))
c----- compute the droplet surface temperature
      IF (LO*HI < 0.0) THEN
        Ts = rootBrent(sdmCalcSpaldingFunc, NIST_T(1), NIST_T(nNIST),
     x       brent_tol)
      else
        Ts = 447.1
      END IF
c----- get the thermo data from NIST
      CALL interpNISTthermo(hfg, cpG, cpL, Psat, Ts)
c----- compute the heat transfer
      Q = hfg + cpL*(Ts - T_liq_0)
c----- compute the Spalding number
      B = cpG*(T_inf - Ts)/Q
c----- set the saturation mixture fraction in the module
      Yfs = (B + Yf_inf)/(1.0 + B)
      Zsat = (nu_f*Yfs - Yox + 0.2153)/(nu_f + 0.2153)
      RETURN
      END FUNCTION sdm_calcSpaldingNo
c-----
c-----------------------------------------------------------------------
c----- function to compute variance sources given the cell conditions
c        input  : f_vol ..... volume of liquid fuel in cell    [m^3]
c                 f_src ..... fuel evaporation source term     [kg]
c                 dt    ..... current time step                [s]
c                 Zmean ..... mean mixture fraction            [-]
c                 Zvar  ..... mixture fraction variance        [-]
c                 mu_L  ..... laminar viscosity                [Pa-s]
c                 B     ..... Spalding number                  [-]
c        output : abOut ..... coefficient of ersatz function   [1/s]
c                 xiOut ..... exponent of ersatz function      [-]
c                 vSrc  ..... variance source from evaporation [1/s]
c                 vEps  ..... variance dissipation             [1/s]
c                 fdt   ..... diffusion factor                 [-]
c-----------------------------------------------------------------------
      SUBROUTINE sdm_calcVarSrc(f_vol,f_src,dt,Zmean,Zvar,mu_L,B,
     x     abOut,xiOut,vSrc,vEps,fdt)
c-----
      USE rootSolve, ONLY : rootBrac, rootBisect
      IMPLICIT NONE
c-----
      REAL(prec), INTENT(in)  :: f_vol, f_src, dt, Zmean, Zvar, mu_L, B
      REAL(prec), INTENT(out) :: vSrc, vEps, fdt, abOut, xiOut
      INTEGER :: j, pdfFlag
      INTEGER, DIMENSION(2) :: pdfBound
      LOGICAL :: bracket
      REAL(prec) :: root_tol=1.0e-4
      REAL(prec) :: mean, var
      REAL(prec) :: LO,HI
      REAL(prec) :: xi,ab,dropD,C
      REAL(prec) :: srcP,srcM,tmp
      REAL(prec) :: xiMin = 1.0 ! maximum exponent for ersatz function  [-]
      REAL(prec) :: xiMax = 4.0 ! maximum exponent for ersatz function  [-]
      REAL(prec), DIMENSION(nZsdm) :: Wv
c----- compute the effective diameter of a single droplet for the cell
      dropD = (six_over_pi*f_vol)**one_third
c----- get the change in mass fraction for this time set
      Zdrop = f_src/massSDM   ! [-]
c----- compute the constant part in alpha
      alpha = 2.0*M_PI*mu_L/(massSDM*Re_t*Sc_t)*LOG(1.0+B)*dropD
c----- set the source term constraint
      sdmW = f_src/(massSDM*dt) ! [-/s]
c----- compute the scaled pdf to be used for integration
      mean = Zmean/Zsat
      var  = Zvar/(Zsat**2)
      pdfFlag=MDUCcalcPDF_1D('beta',nZsdm,eta,mean,var,pdf,pdfBound)
c----- solve for the exponent
      LO = 1.0; HI = 10.0
      CALL rootBrac(sdmCalcMeanSrc, LO, HI, bracket)
      xi = 0.0
      IF (bracket) THEN
        xi = rootBisect(sdmCalcMeanSrc,LO,HI,root_tol)
      END IF
c----- compute the constant for ersatz functional form
      ab = 0.0
      IF (xi >= xiMin .AND. xi < xiMax) THEN
        ab = alpha*Zdrop**(-xi)
      END IF
c----- compute the dynamic timescale ratio
      tau_tauZ = C_x
      C = 0.0
      IF (Zmean > 1.0e-6 .AND. ab > 0.0) THEN
        C = ab*( xi*Zmean**(xi-1.0) - (xi+1.0)*Zmean**xi )
        IF (C*tau > 0.0) THEN
          tau_tauZ = 3.0*C_x*C*tau/(1.0-EXP(-3.0*C*tau))
          fdt      = tau_tauZ/(tau_tauZ - 2.0*C*tau)
        END IF
      END IF
c----- compute the conditional mean source of fuel
      Wv = 0.0
      IF (Zdrop > 0.0 .and. ab > 0.0) THEN
        DO j=2,nZsdm
          Wv(j) = alpha*(eta(j)*Zsat/Zdrop)**xi
        END DO
      END IF
c----- compute the SDM variance source terms
      srcP = 0.0; srcM = 0.0
      DO j=pdfBound(1),pdfBound(2)
        tmp  = (eta(j)*Zsat - Zmean)*Wv(j)*pdf(j)
        srcP = srcP + tmp
        srcM = srcM + (eta(j)*Zsat - Zmean)*tmp
      END DO
      srcP = srcP*Zsat
      srcM = srcM*Zsat
c----- set the source terms appropriately
      vSrc = 2*(1.0 - Zmean)*srcP - srcM
      vEps = MAX(tau_tauZ/tau, C_x/tau) ! make sure source is larger 
      fdt  = MIN(1.5, MAX(1.0, fdt)) ! ensure diffusion coefficient is bounded
      abOut = ab
      xiOut = xi
      END SUBROUTINE sdm_calcVarSrc
c-----
c-----------------------------------------------------------------------
c----- computes the actual source
c-----------------------------------------------------------------------
      SUBROUTINE SDMcalcSrc(mat)
c-----
      USE comm0,    ONLY : ncell,nbfac,nsp,nep,nsb,neb,dt,pref
      USE comm1,    ONLY : den,vol,actphi,pasphi,p,t,vim,te,ed
      USE cthmod,   ONLY : ymfrac,index_fuel,index_o2
      USE cydrnoz,  ONLY : voispr,nnozzl,ninjec,ndrops,sdropn,
     x                     cdropn,pi,ncelln,nozzln
      USE mduc
      IMPLICIT NONE
c-----
      INTEGER, intent(in) :: mat
c----- local variables
      INTEGER    :: nc,nd
      REAL(prec) :: B             ! Spalding number
      REAL(prec) :: TINY = 1.0E-6 ! mixture fraction threshhold
      REAL(prec) :: Pc            ! cell pressure
      REAL(prec) :: Z0mean, Z0var
      REAL(prec) :: pi_over_6, drop_vol
      REAL(prec), DIMENSION(3) :: f_vol, f_src
c-----
c----- set all the variables to the appropriate values
      DO nc=nsp(mat),nep(mat)
        pasphi(nc+iLIQV1)  = 0.0;
        pasphi(nc+iZ1VEPS) = 0.0; pasphi(nc+iZ1VSRC) = 0.0
        pasphi(nc+iFDT1)   = 1.0; pasphi(nc+iXI1)    = 0.0
        pasphi(nc+iALPH1)  = 0.0; pasphi(nc+iSPALD)  = 0.0
      END DO
      IF (mduc_nInj .EQ. 2) THEN
        DO nc=nsp(mat),nep(mat)
          pasphi(nc+iLIQV2)  = 0.0
          pasphi(nc+iZ0VEPS) = 0.0; pasphi(nc+iZ2VEPS) = 0.0
          pasphi(nc+iZ0VSRC) = 0.0; pasphi(nc+iZ2VSRC) = 0.0
          pasphi(nc+iFDT0)   = 1.0; pasphi(nc+iFDT2)   = 1.0
          pasphi(nc+iXI0)    = 0.0; pasphi(nc+iXI2)    = 0.0
          pasphi(nc+iALPH0)  = 0.0; pasphi(nc+iALPH2)  = 0.0
        END DO
      END IF
c----- compute the liquid volume in all cells
      pi_over_6 = pi/6.0
      DO nd=1,ndrops
        nc = ncelln(nd)
        drop_vol = pi_over_6*(sdropn(nd)**3)*cdropn(nd)
        IF (nozzln(nd) .EQ. 1) THEN
          pasphi(nc+iLIQV1) = pasphi(nc+iLIQV1) + drop_vol
        ELSE IF (nozzln(nd) .EQ. 2) THEN
          pasphi(nc+iLIQV2) = pasphi(nc+iLIQV2) + drop_vol
        END IF
      END DO
c----- loop over all cells to compute variance source terms
      DO nc=nsp(mat),nep(mat)
        massSDM = den(nc)*vol(nc)            ! local cell mass     [kg]
        tau     = te(nc)/(ed(nc)+TINY)       ! eddy turnover time  [s]
        Re_t    = den(nc)*te(nc)*tau/vim(nc) ! turb. Reynolds no.  [-]
        Pc      = p(nc)+pref(mat)            ! ambient pressure    [Pa]
c----- set the dissipation terms to the default
        pasphi(nc+iZ1VEPS) = C_x/(tau+1e-6)
c----- get the fuel volume and fuel source for each injection
        f_vol = 0.0
        f_src = 0.0
        f_vol(1) = pasphi(nc+iLIQV1)
        f_src(1) = pasphi(nc+iMEVAP1)
        IF (mduc_nInj .EQ. 2) THEN
          pasphi(nc+iZ0VEPS) = pasphi(nc+iZ1VEPS)
          pasphi(nc+iZ2VEPS) = pasphi(nc+iZ1VEPS)
          f_vol(2) = pasphi(nc+iLIQV2)
          f_src(2) = pasphi(nc+iMEVAP2)
          f_vol(3) = f_vol(1)+f_vol(2)
          f_src(3) = f_src(1)+f_src(2)
          Z0mean   = actphi(nc+iZ1MEAN) + actphi(nc+iZ2MEAN)
        END IF
c----- compute the Spalding number
        IF ( f_vol(1) > 0.0 .OR. f_vol(2) > 0.0 ) THEN
          B = sdm_calcSpaldingNo(ymfrac(nc+(index_fuel-1)*ncell),
     x         ymfrac(nc+(index_o2-1)*ncell), t(nc), Pc)
          pasphi(nc+iSPALD) = B
        END IF
c----- compute the source terms
c     Z1 variance
        IF (f_vol(1) > 0.0 .AND. f_src(1) > 0.0) THEN
          call sdm_calcVarSrc(f_vol(1),f_src(1),dt,
     x         actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),vim(nc),B,
     x         pasphi(nc+iALPH1),pasphi(nc+iXI1),
     x         pasphi(nc+iZ1VSRC),pasphi(nc+iZ1VEPS),pasphi(nc+iFDT1))
        END IF
c     Z2 variance
        IF (f_vol(2) > 0.0 .AND. f_src(2) > 0.0) THEN
          call sdm_calcVarSrc(f_vol(2),f_src(2),dt,
     x         actphi(nc+iZ2MEAN),actphi(nc+iZ2VAR),vim(nc),B,
     x         pasphi(nc+iALPH2),pasphi(nc+iXI2),
     x         pasphi(nc+iZ2VSRC),pasphi(nc+iZ2VEPS),pasphi(nc+iFDT2))
        END IF
c     Z0 variance
        IF (f_vol(3) > 0.0 .AND. f_src(3) > 0.0) THEN
          call sdm_calcVarSrc(f_vol(3),f_src(3),dt,
     x         Z0mean,actphi(nc+iZ0VAR),vim(nc),B,
     x         pasphi(nc+iALPH0),pasphi(nc+iXI0),
     x         pasphi(nc+iZ0VSRC),pasphi(nc+iZ0VEPS),pasphi(nc+iFDT0))
        END IF
      END DO
      END SUBROUTINE SDMcalcSrc
c-----------------------------------------------------------------------
c-----
      END MODULE mduc_sdm

