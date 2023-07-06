c=======================================================================
      SUBROUTINE mduc_scalardiss_1d(mat,n1,x1,chi)
c=======================================================================
c 
c     Purpose:
c       compute an average scalar dissipation rate profile
c
c     Variables:
c       input  - mat            material domain
c              - n1             size of mixture fraction grid
c              - x1             mixture fraction grid
c       update - N/A
c       output - chi            scalar dissipation rate
c             
c=======================================================================
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,numpro,iampro
      USE comm1,    ONLY : den,te,ed,vol,actphi,pasphi
      USE mduc
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    intent(in)                 :: mat,n1
      REAL(prec), intent(in),  DIMENSION(n1) :: x1
      REAL(prec), intent(out), DIMENSION(n1) :: chi
c----- declare local variables
      INTEGER    :: nc,j,jZmin,jZmax,pdfFlag,jChiMax,nZ,totalCells
      INTEGER, DIMENSION(2) :: pdfBound
      REAL(prec) :: taui                       ! turb. eddy turnover time
      REAL(prec) :: sigma                      ! standard deviation
      REAL(prec) :: c_chi = 2.0                ! dissipation coefficient
      REAL(prec) :: chiMean,chiLocal,f1,denom  ! storage variables
      REAL(prec) :: mass
      REAL(prec) :: tmp
      REAL(prec) :: Zmin, Zmax, DZ, xi, xiRef, xiMax, chiMax
      REAL(prec) :: Zplus, Zminus
      REAL(prec), DIMENSION(n1) :: pdf, chiDen
      INTEGER, PARAMETER :: outID = 2, outW = 80
      CHARACTER(len=256) :: text1, text2
c----- make sure that the array is zeroed out
      chi    = 0.0
      xiMax  = 1.0
      pdf    = 0.0
      chiDen = 0.0
c----- find the global min/max of the mixture fraction in the domain
      IF ((chiNorm(1) .eq. 1) .or. (chiNorm(1) .eq. 3)) THEN
        Zmin = 1.0
        Zmax = 0.0
        DO nc=nsp(mat),nep(mat)
          IF (actphi(nc+iZ1MEAN) > 0.0) THEN
c           compute the mixture fraction r.m.s.
            sigma = sqrt(max(actphi(nc+iZ1VAR),0.0))
c           assume min/max is +/- 2*sigma
            Zplus  = actphi(nc+iZ1MEAN) + 2*sigma
            Zminus = actphi(nc+iZ1MEAN) - 2*sigma
            IF (Zminus < Zmin) THEN
              Zmin = max(Zminus, 0.0)
            END IF
            IF (Zplus > Zmax) THEN
              Zmax = min(Zplus, 1.0)
            END IF
          END IF
        END DO
c       reduce across processors
        CALL dglmin(Zmin)
        CALL dglmax(Zmax)
c
        Zmin = max(Zmin, 0.0)
        Zmax = min(Zmax, 1.0)
        IF (chiNorm(1) == 3) Zmin = 0.0
        jZmin = 1
        DO WHILE (x1(jZmin) < Zmin)
          jZmin = jZmin + 1
        END DO
        jZmin = jZmin - 1
        jZmin = max(min(jZmin,n1),1)
c        
        jZmax = 1
        DO WHILE (x1(jZmax) < Zmax)
          jZmax = jZmax + 1
        END DO
        jZmax = max(min(jZmax,n1),1)
c
        Zmin = x1(jZmin)
        Zmax = x1(jZmax)
        DZ   = Zmax - Zmin
      END IF
c-----
c----- start loop over all cells
      nZ = 0
      DO nc=nsp(mat),nep(mat)
        pasphi(nc+iCHI1) = 0.0
        IF (actphi(nc+iZ1MEAN) > Z1_zero) THEN
          mass = den(nc)*vol(nc)
          taui = ed(nc)/(te(nc)+1.0e-6)
c         compute the conditional mean scalar dissipation rate
          chiMean = pasphi(nc+iZ1VEPS)*actphi(nc+iZ1VAR)
c         compute the PDF based on beta function
c         TODO: generalize pdf type
          pdfFlag = MDUCcalcPDF_1D('beta',n1,x1,                        
     x         actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),pdf,pdfBound)
c
          IF (chiNorm(1) == 2 .or. chiNorm(1) == 4) THEN
c           set the bounds to consider for the local pdf
c           this is taken for now as 2*sigma in each direction, as the
c           cutoff for zero given by the PDF calculation is arbitrarily
c           based on the cutoff threshold
            sigma = sqrt(max(actphi(nc+iZ1VAR),0.0))
c           set local cell min/max mixture fraction
            Zmin = 0.0
            IF (chiNorm(1) == 2) THEN
              Zmin = max(actphi(nc+iZ1MEAN)-2*sigma, 0.0)
            END IF
            Zmax = min(actphi(nc+iZ1MEAN)+2*sigma, 1.0)
c            
            jZmin = 1
            DO WHILE (x1(jZmin) < Zmin)
              jZmin = jZmin + 1
            END DO
            jZmin = jZmin - 1
            jZmin = max(min(jZmin,n1),1)
c            
            jZmax = 1
            DO WHILE (x1(jZmax) < Zmax)
              jZmax = jZmax + 1
            END DO
            jZmax = max(min(jZmax,n1),1)      
            Zmin = x1(jZmin)
            Zmax = x1(jZmax)
            DZ   = Zmax - Zmin
          END IF
c-----    set bounds based on min/max locations
          IF (pdfBound(1) < jZmin) pdfBound(1) = jZmin
          IF (pdfBound(2) > jZmax) pdfBound(2) = jZmax
c        
          denom = 0.0
          DO j=1,n1
            f1 = MDUCcalcChiLogForm(x1(j),actphi(nc+iZ1MEAN),Zmax)
            denom = denom + f1*pdf(j)
          END DO
c          
          IF (denom > 1.0E-9) pasphi(nc+iCHI1) = chiMean/denom
c     TODO: put bounds on pasphi directly
          DO j=1,n1
            f1 = MDUCcalcChiLogForm(x1(j),actphi(nc+iZ1MEAN),Zmax)
            chiLocal  = pasphi(nc+iCHI1)*f1
            chi(j)    = chi(j)    + mass*chiLocal
            chiDen(j) = chiDen(j) + mass
          END DO
          nZ = nZ+1
        END IF
      END DO 
c----- reduce variables across processors
      CALL iglsum(nZ)
      CALL dglsumvec(chi,n1)
      CALL dglsumvec(chiDen,n1)
c----- get the total cells across processors
      totalCells = ncell
      CALL iglsum(totalCells)
c----- finish normalization
      jChimax = 1
      chiMax = 0.0
      DO j=1,n1
        IF (chiDen(j) > 0.0) chi(j) = chi(j)/chiDen(j)
c       keep track of maximum scalar dissipation rate
        IF (chi(j) .GT. chiMax) THEN
          chiMax = chi(j)
          jChiMax = j
        END IF
      END DO
c----- write summary
      tmp = real(100.0/totalCells)
      IF (iampro .eq. rootpro) THEN
        CALL mduc_out_header(outW,'Scalar Dissipation Rate (1D) Output')
        CALL mduc_out_int(outW,outID,'Mixture cells',nZ)
        CALL mduc_out_real(outW,outID,'Mixture percent domain','%',
     x       real(nZ/totalCells*100.0),0)
        CALL mduc_out_real(outW,outID,'Maximum mixture fraction (Z1)',
     x       '-',Zmax,1)
        CALL mduc_out_real(outW,outID,'Maximum scalar dissipation rate',
     x       '1/s',chiMax,1)
        CALL mduc_out_real(outW,outID,'Mixture fraction at max chi',
     x       '-',x1(jChiMax),1)
        WRITE(text1,'(A)') repeat('*',outW-2)
        CALL mduc_out(outW,text1)
      END IF
      
      END SUBROUTINE mduc_scalardiss_1d

c=======================================================================
      SUBROUTINE mduc_scalardiss_2d(mat,n1,n2,x1,x2,chi1,chi12,chi2)
c=======================================================================
c 
c    Purpose:
c      compute an average scalar dissipation rate profile
c
c    Variables:
c      input  - mat            material domain
c             - n1,n2          size of mixture fraction grid
c             - x1,x2          mixture fraction grid
c      update - N/A
c      output - chi1           scalar dissipation rate of Z1
c             - chi12          cross scalar dissipation rate
c             - chi2           scalar dissipation rate of Z2
c             
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,iampro,numpro
      USE comm1,    ONLY : den,te,ed,vol,actphi,pasphi
      USE mduc
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    intent(in)                    :: n1,n2,mat
      REAL(prec), intent(in),  DIMENSION(n1)    :: x1
      REAL(prec), intent(in),  DIMENSION(n2)    :: x2
      REAL(prec), intent(out), DIMENSION(n1,n2) :: chi1, chi12, chi2
c-----  declare local variables
      INTEGER    :: nc,i,j,totalCells
      INTEGER    :: nZ, nZ1, nZ2                    ! counters
      INTEGER    :: nZ1clip, nZ2clip
      INTEGER    :: mduc_grid_locate
      REAL(prec) :: taui                            ! inv. turb. timescale
      REAL(prec) :: sigma1,sigma2                   ! standard deviations
      REAL(prec) :: c_chi = 2.0                     ! dissipation coefficient
      REAL(prec) :: mass
      REAL(prec) :: tmp
      REAL(prec) :: TINY_DEN = 1.0E-12              ! denominator threshhold
      REAL(prec), DIMENSION(3)     :: denom   
      INTEGER, PARAMETER :: outID = 2, outW = 80
      CHARACTER(len=256) :: text1, text2
c-----  mixture fraction variables
      REAL(prec)                   :: Zplus,Zminus
      REAL(prec), DIMENSION(3)     :: ZminG,ZmaxG
      REAL(prec), DIMENSION(n2)    :: Z1min,Z1max
      REAL(prec), DIMENSION(n1)    :: Z2min,Z2max
      REAL(prec), DIMENSION(n1,n2) :: Z1,Z2          ! untra
c-----  temporary chi variables
      REAL(prec) :: chi0Mean,chi1Mean,chi2Mean
      REAL(prec) :: chi0Local,chi1Local,chi2Local
      REAL(prec) :: chi0Max, chi1Max, chi2Max, chi12Max
      INTEGER,    DIMENSION(4)     :: iChiMax, jChiMax
      REAL(prec), DIMENSION(n1,n2) :: f0, f1, f2
      REAL(prec), DIMENSION(n1,n2) :: chi0,chi0D,chi1D,chi2D
c----- pdf variables
      INTEGER                      :: pdfFlag
      INTEGER,    DIMENSION(2)     :: pdfBound1D
      REAL(prec), DIMENSION(n1)    :: pdfX1
      REAL(prec), DIMENSION(n2)    :: pdfX2
      REAL(prec), DIMENSION(n1,n2) :: pdf
c
c-----  make sure that the arrays are zeroed
      chi1 = 0.0; chi12 = 0.0;  chi2 = 0.0
      chi0 = 0.0; chi0D = 0.0; chi1D = 0.0; chi2D = 0.0
      pdf  = 0.0; pdfX1 = 0.0; pdfX2 = 0.0
c     TODO: set switch for functional form of chi
c     fill the non-transformed value grid using appropriate transformation
      DO j=1,n2
        DO i=1,n1
c         standard transformation for now
c         TODO: make general transformation
          Z1(i,j) = x1(i)
          Z2(i,j) = x2(j)*(1.0 - x1(i))
        END DO
      END DO
c     find the global min/max of each mixture fraction in the domain
      Z1min = 1.0; Z1max = 0.0
      Z2min = 1.0; Z2max = 0.0
      IF (chiNorm(1) == 1 .or. chiNorm(1) == 3) THEN
        Z1min = 1.0
        Z1max = 0.0
        DO nc=nsp(mat),nep(mat)
c-----    check Z1
          IF (actphi(nc+iZ1MEAN) > 0.0) THEN
c           get r.m.s. value
            sigma1 = sqrt(max(actphi(nc+iZ1VAR),0.0))
c           find where to condition
            j = mduc_grid_locate(n2, x2, actphi(nc+iZ2MEAN))
c           assume min/max is +/- 2*sigma
            Zplus  = actphi(nc+iZ1MEAN) + 2*sigma1
            Zminus = actphi(nc+iZ1MEAN) - 2*sigma1
            IF (Zminus < Z1min(j)) THEN
              Z1min(j) = max(Zminus, 0.0)
            END IF
            IF (Zplus > Z1max(j)) THEN
              Z1max(j) = min(Zplus, 1.0)
            END IF
          END IF
c-----    check Z2
          IF (actphi(nc+iZ2MEAN) > 0.0) THEN
c           get r.m.s. value
            sigma2 = sqrt(max(actphi(nc+iZ2VAR),0.0))
c           find where to condition
            i = mduc_grid_locate(n1, x1, actphi(nc+iZ1MEAN))
c           assume min/max is +/- 2*sigma
            Zplus  = actphi(nc+iZ2MEAN) + 2*sigma2
            Zminus = actphi(nc+iZ2MEAN) - 2*sigma2
            IF (Zminus < Z2min(i)) THEN
              Z2min(i) = max(Zminus, 0.0)
            END IF
            IF (Zplus > Z2max(i)) THEN
              Z2max(i) = min(Zplus, 1.0)
            END IF
          END IF
        END DO
c----- get min/max across processors
        CALL dglminvec(Z1min, n2)
        CALL dglmaxvec(Z1max, n2)
        CALL dglminvec(Z2min, n1)
        CALL dglmaxvec(Z2max, n1)
c       get the global min/max values for the domain
        ZminG = 1.0; ZmaxG = 0.0
        DO j=1,n2
          IF (Z1min(j) < ZminG(1)) ZminG(1) = Z1min(j)
          IF (Z1max(j) > ZmaxG(1)) ZmaxG(1) = Z1max(j)
        END DO
        ZmaxG(1) = max(ZmaxG(1), 0.0)
        DO i=1,n1
          IF (Z2min(i) < ZminG(2)) ZminG(2) = Z2min(i)
          IF (Z2max(i) > ZmaxG(2)) ZmaxG(2) = Z2max(i)
        END DO
        ZmaxG(2) = max(ZmaxG(2), 0.0)    
        ZmaxG(3) = max(1.0 - ZmaxG(1) - ZmaxG(2), 0.0)
c       if specified, force zero as global min
        IF (chiNorm(1) == 3) ZminG = 0.0
      END IF
c     start loop over all cells
      nZ = 0; nZ1 = 0; nZ2 = 0
      nZ1clip = 0
      nZ2clip = 0
      DO nc=nsp(mat),nep(mat)
c     ensure all storage memory is reset to zero
        pasphi(nc+iCHI0) = 0.0
        pasphi(nc+iCHI1) = 0.0
        pasphi(nc+iCHI2) = 0.0
c     
        f0 = 0.0; f1 = 0.0; f2 = 0.0
c     
        IF (actphi(nc+iZ1MEAN) > Z1_zero) THEN
          mass = den(nc)*vol(nc)
          taui = ed(nc)/te(nc)
c         check if the cell contains both mixture fractions
          IF (actphi(nc+iZ2MEAN) > Z2_zero) THEN
c           compute the conditional mean scalar dissipation rate
            chi0Mean = pasphi(nc+iZ0VEPS)*actphi(nc+iZ0VAR)
            chi1Mean = pasphi(nc+iZ1VEPS)*actphi(nc+iZ1VAR)
            chi2Mean = pasphi(nc+iZ2VEPS)*actphi(nc+iZ2VAR)
c           compute the PDF based on beta function
            pdfFlag = MDUCcalcPDF_2D('dirichlet',n1, x1, n2, x2,        
     x           actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),                 
     x           actphi(nc+iZ2MEAN),actphi(nc+iZ2VAR),actphi(nc+iZ0VAR),
     x           pdf,1)
           
            denom = 0.0
c           compute the integral of the functional form to get mean 
c           scalar dissipation
            DO j=1,n2
              DO i=1,n1
                f0(i,j) = MDUCcalcChiLogForm(
     x               Z1(i,j)/ZmaxG(1) + Z2(i,j)/ZmaxG(2),
     x               actphi(nc+iZ1MEAN)/ZmaxG(1) +
     x               actphi(nc+iZ2MEAN)/ZmaxG(2), 1.0)
                f1(i,j) = MDUCcalcChiLogForm(Z1(i,j),
     x               actphi(nc+iZ1MEAN),ZmaxG(1))
                f2(i,j) = MDUCcalcChiLogForm(Z2(i,j),
     x               actphi(nc+iZ2MEAN),ZmaxG(2))
c                
                denom(1) = denom(1) + f1(i,j)*pdf(i,j)
                denom(2) = denom(2) + f2(i,j)*pdf(i,j)
                denom(3) = denom(3) + f0(i,j)*pdf(i,j)
              END DO
            END DO
c           compute the reference scalar dissipation rate
            IF (denom(1) > TINY_DEN) pasphi(nc+iCHI1)=chi1Mean/denom(1)
            IF (denom(2) > TINY_DEN) pasphi(nc+iCHI2)=chi2Mean/denom(2)
            IF (denom(3) > TINY_DEN) pasphi(nc+iCHI0)=chi0Mean/denom(3)
c
            IF (pasphi(nc+iCHI1) .GT. 100.0) nZ1clip = nZ1clip + 1
            IF (pasphi(nc+iCHI2) .GT. 100.0) nZ2clip = nZ2clip + 1
c           TODO: only use these limiters temporarily
c           pasphi(nc+iCHI0) = min(pasphi(nc+iCHI0), 100.0)
c           pasphi(nc+iCHI1) = min(pasphi(nc+iCHI1), 100.0)
c           pasphi(nc+iCHI2) = min(pasphi(nc+iCHI2), 100.0)
c-----      scale the distribution and add contribution to mass avg.
            chi0Local = 0
            chi1Local = 0
            chi2Local = 0
            DO j=1,n2
              DO i=1,n1
                IF (Z1(i,j) <= ZmaxG(1) .AND. Z2(i,j) <= ZmaxG(2)
     x               .AND. Z1(i,j)+Z2(i,j) <= 1.0) THEN
                  chi0Local = pasphi(nc+iCHI0)*f0(i,j)
                  chi1Local = pasphi(nc+iCHI1)*f1(i,j)
                  chi2Local = pasphi(nc+iCHI2)*f2(i,j)
c                  
                  chi0(i,j)  = chi0(i,j)  + mass*chi0Local
                  chi0D(i,j) = chi0D(i,j) + mass
                  chi1(i,j)  = chi1(i,j)  + mass*chi1Local
                  chi1D(i,j) = chi1D(i,j) + mass
                  chi2(i,j)  = chi2(i,j)  + mass*chi2Local
                  chi2D(i,j) = chi2D(i,j) + mass
                END IF
              END DO
            END DO
            nZ = nZ + 1
          else                  ! this cell only has Z1
c-----    compute the conditional mean scalar dissipation rate
            chi1Mean = pasphi(nc+iZ1VEPS)*actphi(nc+iZ1VAR)
c          get the 1D pdf
            pdfFlag = MDUCcalcPDF_1D('beta',n1,x1,
     x           actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),pdfX1,pdfBound1D)
c           compute f1(x) based on assumed functional form
            denom(1) = 0.0
            DO i=1,n1
c             compute functional form
              f1(i,1) = MDUCcalcChiLogForm(x1(i),actphi(nc+iZ1MEAN),
     x             ZmaxG(1))
c              compute contribution to integral
              denom(1) = denom(1) + f1(i,1)*pdfX1(i)
            END DO
c           compute the reference scalar dissipation rate
            IF (denom(1) > TINY_DEN) pasphi(nc+iCHI1)=chi1Mean/denom(1)
c             TODO: remove
c             pasphi(nc+iCHI1) = min(pasphi(nc+iCHI1), 100.0)
c             scale the distribution and add contribution to mass avg.
            DO i=1,n1
              chi1Local = pasphi(nc+iCHI1)*f1(i,1)
              chi0(i,1)  = chi0(i,1)  + mass*chi1Local
              chi0D(i,1) = chi0D(i,1) + mass
              chi1(i,1)  = chi1(i,1)  + mass*chi1Local
              chi1D(i,1) = chi1D(i,1) + mass
            END DO
            nZ1 = nZ1 + 1
          END IF
        else IF (actphi(nc+iZ2MEAN) > Z2_zero) THEN ! this cell only has Z2
          mass     = den(nc)*vol(nc)
          taui     = ed(nc)/te(nc)
          chi2Mean = pasphi(nc+iZ2VEPS)*actphi(nc+iZ2VAR)
c         get the 1D pdf
          pdfFlag = MDUCcalcPDF_1D('beta',n2,x2,actphi(nc+iZ2MEAN),
     x         actphi(nc+iZ2VAR),pdfX2,pdfBound1D)
c         compute f1(x) based on assumed functional form
          denom(2) = 0.0
          DO j=1,n2
c           compute functional form
            f1(1,j) = MDUCcalcChiLogForm(x2(j),actphi(nc+iZ2MEAN),
     x           ZmaxG(2))
c           compute contribution to integral
            denom(2) = denom(2) + f1(1,j)*pdfX2(j)
          END DO
c-----    compute the reference scalar dissipation rate
          IF (denom(2) > TINY_DEN) pasphi(nc+iCHI2) = chi2Mean/denom(2)
c         TODO: remove
c         pasphi(nc+iCHI2) = min(pasphi(nc+iCHI2), 100.0)
c-----    scale the distribution and add contribution to mass avg.
          DO j=1,n2
            chi2Local = pasphi(nc+iCHI2)*f1(1,j)
            chi0(1,j)  = chi0(1,j)  + mass*chi2Local
            chi0D(1,j) = chi0D(1,j) + mass
            chi2(1,j)  = chi2(1,j)  + mass*chi2Local
            chi2D(1,j) = chi2D(1,j) + mass
          END DO
          nZ2 = nZ2 + 1
        END IF
      END DO
c----- reduce variables across processors
      CALL iglsum(nZ)
      CALL iglsum(nZ1)
      CALL iglsum(nZ2)
      CALL iglsum(nZ1clip)
      CALL iglsum(nZ2clip)
      CALL dglsumvec(chi0,n1*n2)
      CALL dglsumvec(chi0D,n1*n2)
      CALL dglsumvec(chi1,n1*n2)
      CALL dglsumvec(chi1D,n1*n2)
      CALL dglsumvec(chi2,n1*n2)
      CALL dglsumvec(chi2D,n1*n2)
c-----  get the total cells across processors
      totalCells = ncell
      CALL iglsum(totalCells)
c-----  normalize the distribution
      chi0Max = 0.0; chi1Max = 0.0; chi2Max = 0.0; chi12Max = 0.0
      iChiMax = 1; jChiMax = 1
      DO j=1,n2
        DO i=1,n1
          IF (chi0D(i,j) > TINY_DEN) chi0(i,j) = chi0(i,j)/chi0D(i,j)
          IF (chi1D(i,j) > TINY_DEN) chi1(i,j) = chi1(i,j)/chi1D(i,j)
          IF (chi2D(i,j) > TINY_DEN) chi2(i,j) = chi2(i,j)/chi2D(i,j)
c         compute a consistent cross-dissipation rate
          chi12(i,j) = 0.5*(chi0(i,j) - chi1(i,j) - chi2(i,j))
c         keep track of maximum scalar dissipation rates
          IF (chi0(i,j) .GT. chi0Max) THEN
            chi0Max = chi0(i,j)
            iChiMax(3) = i
            jChiMax(3) = j
          END IF
          IF (chi1(i,j) .GT. chi1Max) THEN
            chi1Max = chi1(i,j)
            iChiMax(1) = i
            jChiMax(1) = j
          END IF
          IF (chi2(i,j) .GT. chi2Max) THEN 
            chi2Max = chi2(i,j)
            iChiMax(2) = i
            jChiMax(2) = j
          END IF
          IF (abs(chi12(i,j)) .GT. abs(chi12Max)) THEN
            chi12Max = chi12(i,j)
            iChiMax(4) = i
            jChiMax(4) = j
          END IF
        END DO
      END DO
c----- enforce zero boundary conditions as appropriate
      DO j=1,n2
        chi1(1,j)  = 0.0; chi1(n1,j) = 0.0
        chi2(n1,j) = 0.0
        chi12(1,j) = 0.0; chi1(n1,j) = 0.0
      END DO
c-----  write summary
      tmp = real(100.0/totalCells) ;
      IF (iampro .EQ. rootpro) THEN
        CALL mduc_out_header(outW,'Scalar Dissipation Rate (2D) Output')
        CALL mduc_out_int(outW,outID,'Total mixture cells',nZ+nZ1+nZ2)
        CALL mduc_out_int(outW,outID,'Mixture 1 cells (Z1)',nZ1)
        CALL mduc_out_int(outW,outID,'Mixture 2 cells (Z2)',nZ2)
        CALL mduc_out_int(outW,outID,'Mixture 1+2 cells (Z1+Z2)',
     x       nZ1+nZ2)
        CALL mduc_out_real(outW,outID,'Mixture percent domain','%',
     x       (nZ+nZ1+nZ2)*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture 1 percent domain','%',
     x       nZ1*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture 2 percent domain','%',
     x       nZ2*tmp,0)
        CALL mduc_out_real(outW,outID,'Mixture 1+2 percent domain','%',
     x       (nZ1+nZ2)*tmp,0)

        CALL mduc_out_real(outW,outID,'Maximum mixture 1 fraction (Z1)',
     x       '-',ZmaxG(1),1)
        CALL mduc_out_real(outW,outID,'Maximum mixture 2 fraction (Z2)',
     x       '-',ZmaxG(2),1)
        CALL mduc_out_real(outW,outID,'Maximum mixture 0 fraction (Z0)',
     x       '-',ZmaxG(3),1)

        WRITE(text1,'(A)') repeat('*',outW-2)
        CALL mduc_out(outW,text1)
c---- table for maximum scalar dissipation rates
        text1 = '  %M02  Max Scalar Dissipation Rates'
        CALL mduc_out(outW,text1)
        WRITE(text1,'(20X,A10,A10,A12)') 'Z1','Z2','Max'
        CALL mduc_out(outW,text1)
        WRITE(text1,'(A20,F10.2,F10.2,ES12.2)') 'chi0',
     x       Z1(iChiMax(3),jChiMax(3)),Z2(iChiMax(3),jChiMax(3)),
     x       chi0Max
        CALL mduc_out(outW,text1)
        WRITE(text1,'(A20,F10.4,F10.4,ES12.2)') 'chi1',
     x       Z1(iChiMax(1),jChiMax(1)),Z2(iChiMax(1),jChiMax(1)),
     x       chi1Max
        CALL mduc_out(outW,text1)
        WRITE(text1,'(A20,F10.4,F10.4,ES12.2)') 'chi12',
     x       Z1(iChiMax(4),jChiMax(4)),Z2(iChiMax(4),jChiMax(4)),
     x       chi12Max
        CALL mduc_out(outW,text1)
        WRITE(text1,'(A20,F10.4,F10.4,ES12.2)') 'chi2',
     x       Z1(iChiMax(2),jChiMax(2)),Z2(iChiMax(2),jChiMax(2)),
     x       chi2Max
        CALL mduc_out(outW,text1)
        WRITE(text1,'(A)') repeat('*',outW-2)
        CALL mduc_out(outW,text1)
      END IF
c-----      
      END SUBROUTINE mduc_scalardiss_2d
c-----
