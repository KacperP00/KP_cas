c=======================================================================
      SUBROUTINE mduc_soot_psc(mat)
c=======================================================================
c     
c     Purpose:
c       compute the source terms for the soot moment equations and add
c       them to the appropriate active scalars
c    
c     Variables:
c       input  - mat            
c       update - actphr(iSoot:iSoot+nSootMom)
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,num_actscl,num_passcl,dtin,
     x                     pref,pref_old,dt,iter,time,timo,nsb,neb,nreg,
     x                     nbfac,iampro
      USE comm1,    ONLY : actphi,actphb,actphr,pasphi,pasphr,
     x                     vol,volo,den,denb,t,tb,ed,edb,p,po,te,teb,
     x                     vit,vitb,nsr,ner,lb,vim
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER, INTENT(IN) :: mat
c----- define local soot variables
      INTEGER    :: istat,nc,iM
      REAL(prec) :: mass,sootRhoSrc,mom,tmp
      REAL(prec), DIMENSION(nSootMom) :: sootSrc,sootTP
      REAL(prec), ALLOCATABLE, DIMENSION(:,:)     :: grad_T,tpb_d
      REAL(prec), ALLOCATABLE, DIMENSION(:,:,:)   :: tp_in,tp_d
      REAL(prec), ALLOCATABLE, DIMENSION(:,:,:)   :: tp_g1,tp_g2,tp_g3
      REAL(prec), ALLOCATABLE, DIMENSION(:,:,:,:) :: tp_out
c-----------------------------------------------------------------------
c-----
c      compute the gradient of temperature field
      ALLOCATE(grad_T(3,ncell), STAT=istat); grad_T = 0.0
      CALL gradfi(t,tb,grad_T,1,1,mat)
c----- compute the gradients necessary for the thermophoresis
      ALLOCATE(tp_g1(nSootMom,3,ncell), STAT=istat); tp_g1   = 0.0
      ALLOCATE(tp_g2(nSootMom,3,ncell), STAT=istat); tp_g2   = 0.0
      ALLOCATE(tp_g3(nSootMom,3,ncell), STAT=istat); tp_g3   = 0.0
      ALLOCATE(tp_d(nSootMom,3,ncell),  STAT=istat); tp_d    = 0.0
      ALLOCATE(tpb_d(3,nbfac),          STAT=istat); tpb_d   = 0.0
c     get the vector of 0.55*mu/(rho*T)*M_(x,y)*dT/dx_j
      DO nc=nsp(mat),nep(mat)
        tmp = 0.55*vim(nc)/(den(nc)*t(nc))
        DO iM=1,nSootMom
          mom = actphi(nc+iSOOT_MOM+(iM-1)*ncell)
          tp_d(iM,1,nc) = tmp*mom*grad_T(1,nc)
          tp_d(iM,2,nc) = tmp*mom*grad_T(2,nc)
          tp_d(iM,3,nc) = tmp*mom*grad_T(3,nc)
        END DO
      END DO
c     compute the necessary gradients of the bracketed term
      DO iM=1,nSootMom
        CALL gradfi(tp_d(iM,:,:),tpb_d,tp_g1(iM,:,:),3,1,mat)
        CALL gradfi(tp_d(iM,:,:),tpb_d,tp_g2(iM,:,:),3,2,mat)
        CALL gradfi(tp_d(iM,:,:),tpb_d,tp_g3(iM,:,:),3,3,mat)
      END DO
c----- compute the soot source terms for Hybrid Method of Moments Eqs.
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)
c-----  call the external soot library to get the source terms
        IF (sootSrcUpdate .EQ. 1) THEN
          tmp = 0.0
          CALL MDUCmeanSootSrc(actphi(nc+iZ1MEAN),actphi(nc+iZ1VAR),
     x         tmp,tmp,tmp,sootSrc,mduc_mem)
        END IF
c       add the thermophoresis term
c         d/dx(C*dT/dx) + d/dy(C*dT/dy) + d/dz(C*dT/dz)            
        sootTP = 0.0
        DO iM=1,nSootMom
          sootTP(iM) = sootTP(iM) + 
     x         tp_g1(iM,1,nc) + tp_g2(iM,2,nc) + tp_g3(iM,3,nc)
c          sootSrc(iM) = sootSrc(iM) + sootTP(iM);
c         store the 
          IF (mduc_outSootSrc .EQ. 1) THEN
            pasphi(nc+iSOOT_SRC+(iM-1)*ncell) = den(nc)*sootSrc(iM)
          END IF
          IF (mduc_outSootTP .EQ. 1) THEN
            pasphi(nc+iSOOT_TP+(iM-1)*ncell) = den(nc)*sootTP(iM)
          END IF
        END DO
c-----  set the source terms for the first order moments
c       M_00: number density 1/(m^3-s)*kg/m^3
        actphr(nc+iSOOT_M00) = den(nc)*sootSrc(1)
c       M_10: volume fraction m^3/(m^3-s)*kg/m^3
        actphr(nc+iSOOT_M10) = den(nc)*sootSrc(2)
c       M_01: surface area density m^2/(m^3-s)*kg/m^3
        actphr(nc+iSOOT_M01) = den(nc)*sootSrc(3)
c-----  compute the source terms for the second order moments
        IF (mduc_soot .GT. 1) THEN
c         M_20: volume variance m^6/(m^3-s)*kg/m^3
          actphr(nc+iSOOT_M20) = den(nc)*sootSrc(4)
c         M_11: volume-surface covariance m^5/(m^3-s)*kg/m^3
          actphr(nc+iSOOT_M11) = den(nc)*sootSrc(5)
c         M_02: surface variance m^4/(m^3-s)*kg/m^3
          actphr(nc+iSOOT_M02) = den(nc)*sootSrc(6)
        END IF
c       N_0: weight of the delta function 1/(m^3-s)*kg/m^3
        actphr(nc+iSOOT_N0) = den(nc)*sootSrc(nSootMom)
c----- add the source term to the continuity equation kg/(m^3-s)
      END DO
c     
      IF(ALLOCATED(grad_T))  DEALLOCATE(grad_T)
      IF(ALLOCATED(tp_g1))   DEALLOCATE(tp_g1)
      IF(ALLOCATED(tp_g2))   DEALLOCATE(tp_g2)
      IF(ALLOCATED(tp_g3))   DEALLOCATE(tp_g3)
      IF(ALLOCATED(tp_d))    DEALLOCATE(tp_d)
      IF(ALLOCATED(tpb_d))   DEALLOCATE(tpb_d)
c-----      
      END SUBROUTINE mduc_soot_psc
c-----
c=======================================================================
      SUBROUTINE mduc_soot_init_moments(mat)
c=======================================================================
c     
c     Purpose:
c       initialize the soot moments to a small number density, with other
c       moments determined by M_xy = M_00*40^(x+2/3*y)
c    
c     Variables:
c       input  - mat            
c       update - actphi(iSoot:iSoot+nSootMom)
c                 
c=======================================================================
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,iampro
      USE comm1,    ONLY : actphi,actphb
      USE mduc,     ONLY : nSootMom,iSOOT_MOM
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER, INTENT(IN)      :: mat
c----- define local soot variables
      INTEGER                  :: nc,iM
      REAL(prec), DIMENSION(nSootMom) :: sootMom0
c-----------------------------------------------------------------------
c-----
c     initialize the soot moments to a small number density, with other
c     moments determined by M_xy = M_00*40^(x+2/3*y)
c----- first order moments
c     M_00: number density (1/m^3)
      sootMom0(1) = 1.0e-20
c     M_10: volume fraction (m^3/m^3)
      sootMom0(2) = 1.0e-20*40**(1.0+2.0/3.0*0.0)
c     M_01: surface area density (m^2/m^3)
      sootMom0(3) = 1.0e-20*40**(0.0+2.0/3.0*1.0)
c----- second order moments
      IF (nSootMom .GT. 4) THEN
c       M_20: volume variance (m^6/m^3)
        sootMom0(4) = 1.0e-20*40**(2.0+2.0/3.0*0.0)
c       M_11: volume-surface covariance (m^5/m^3)
        sootMom0(5) = 1.0e-20*40**(1.0+2.0/3.0*1.0)
c       M_02: surface variance (m^4/m^3)
        sootMom0(6) = 1.0e-20*40**(0.0+2.0/3.0*2.0)
      END IF
c     N_0: weight of the delta function (1/m^3)
      sootMom0(nSootMom) = 0.5e-20;
c----- assign the field variables to the initial values
      DO nc=nsp(mat),nep(mat)
        DO iM=1,nSootMom
          actphi(nc+iSOOT_MOM+(iM-1)*ncell) = sootMom0(iM)
        END DO
      END DO
c-----      
      END SUBROUTINE mduc_soot_init_moments
c-----
c=======================================================================
      SUBROUTINE mduc_soot_cond_1d(mat,n1,x1,mom_z)
c=======================================================================
c     
c     Purpose:
c       condition the soot moments on mixture fraction
c    
c     Variables:
c       input  - mat            material domain (from FIRE)
c              - n1             size of mixture fraction grid
c              - x1             mixture fraction grid
c       output - mom_z          global mass averaged pdf
c                 
c=======================================================================
      USE prec_mod,  ONLY : prec
      USE comm0,     ONLY : ncell,nsp,nep,iampro
      USE comm1,     ONLY : den,vol,actphi
      USE mduc,      ONLY : rootpro,iZ1MEAN,iZ1VAR,nSootMom,iSOOT_MOM
      USE rootsolve, ONLY : mduc_grid_locate
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN)                          :: mat,n1
      REAL(prec), INTENT(IN),  DIMENSION(n1)          :: x1
      REAL(prec), INTENT(OUT), DIMENSION(n1,nSootMom) :: mom_z
c----- declare local variables
      INTEGER    :: i,j,nc,totalCells,iM
      INTEGER    :: nZ1
      REAL(prec) :: mass,dx,mom,w1,w2,perc
      REAL(prec), DIMENSION(nSootMom)    :: momMax,momZmax
      REAL(prec), DIMENSION(n1,nSootMom) :: mom_z_den
      INTEGER, PARAMETER :: outID = 5, outW = 80
      CHARACTER(len=256) :: text1,text2
c-----------------------------------------------------------------------
      nZ1       = 0
      mom_z     = 0.0
      mom_z_den = 0.0
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)
c----- compute the mean mass fraction of all species
        IF (actphi(nc+iZ1MEAN) > 1.0E-6) THEN
c         search for the bin that this value fits in
          j  = mduc_grid_locate(n1,x1,actphi(nc+iZ1MEAN))
          dx = x1(j+1)-x1(j)
c         compute the weights
          w1 = (x1(j+1)            - actphi(nc+iZ1MEAN))/dx
          w2 = (actphi(nc+iZ1MEAN) - x1(j))/dx
c         add to each moment distribution
          DO iM=1,nSootMom
            mom = actphi(nc+iSOOT_MOM+(iM-1)*ncell)
            mom_z(j,iM)       = mom_z(j,iM)       + mass*w1*mom
            mom_z_den(j,iM)   = mom_z_den(j,iM)   + mass*w1
            mom_z(j+1,iM)     = mom_z(j+1,iM)     + mass*w2*mom
            mom_z_den(j+1,iM) = mom_z_den(j+1,iM) + mass*w2
          END DO
          nZ1 = nZ1 + 1
        END IF
      END DO
c----- get the total cells for this processor
      totalCells = ncell
c----- reduce variables across processors
      CALL iglsum(totalCells)
      CALL iglsum(nZ1)
      CALL dglsumvec(mom_z,n1*nSootMom)
      CALL dglsumvec(mom_z_den,n1*nSootMom)
c----- normalise the global pdf
      DO j=1,n1
        DO iM=1,nSootMom
          IF (mom_z_den(j,iM) .GT. 0.0) THEN
            mom_z(j,iM) = mom_z(j,iM)/mom_z_den(j,iM)
c          ELSE
c            mom_z(j,iM) = 0.0
          END IF
        END DO
      END DO
c      mom_z(1,:)  = 0.0
c      mom_z(n1,:) = 0.0
c----- write summary
      IF (iampro .EQ. rootpro) THEN
        perc = real(nZ1/totalCells*100.0)
        momMax = -1e20
        DO j=1,n1
          DO iM=1,nSootMOM
            IF (mom_z(j,iM) > momMax(iM)) THEN
              momMax(iM)  = mom_z(j,iM)
              momZmax(iM) = x1(j)
            END IF
          END DO
        END DO
        CALL mduc_out_header(outW, "Conditioned Soot Moments")
        CALL mduc_out_int(outW,outID,'Mixture cells updated',nZ1)
        CALL mduc_out_real(outW,outID,'Mixture percent domain','%',
     x      perc,0)
        CALL mduc_out_real(outW,outID,'Max soot number density','1/m^3',
     x       momMax(1),1)
        CALL mduc_out_real(outW,outID,'Max soot number density Z','-',
     x       momZmax(1),1)
        CALL mduc_out_real(outW,outID,'Max soot volume fraction',
     x       'm^3/m^3',momMax(2),1)
        CALL mduc_out_real(outW,outID,'Max soot volume fraction Z','-',
     x       momZmax(2),1)
        CALL mduc_out_real(outW,outID,'Max soot area','m^2/m^3',
     x       momMax(3),1)
        CALL mduc_out_real(outW,outID,'Max soot area Z','-',
     x       momZmax(3),1)

c        WRITE(text2, '(A20,I8)') 'nCell ',totalCells
c        CALL mduc_out(outW, text2)

        WRITE(text1, '(A)') REPEAT('*', outW-2)
        CALL mduc_out(outW, text1)
        WRITE (*,*) ''
c       write to file
        OPEN(5001,file='sootMom.dat',status='unknown',position='append')
        IF (nSootMom > 4) THEN
          WRITE(5001,'(8A15)') '# Z1<1>', 'n_M00<2>',
     x         'vol_M10<3>','area_M01<4>','vol_v_M20<5>',
     x         'area_v_M02<6>','vol_area_v_M11<7>','delta_N0<8>'
        ELSE
          WRITE(5001,'(5A15)') '# Z1<1>', 'n_M00<2>',
     x         'vol_M10<3>','area_M01<4>','delta_N0<5>'
        END IF
        DO j=1,n1
          IF (nSootMOM > 4) THEN
            WRITE(5001,'(8ES15.6)') x1(j),mom_z(j,1),mom_z(j,2),
     x           mom_z(j,3),mom_z(j,4),mom_z(j,5),mom_z(j,6),mom_z(j,7)
          ELSE
            WRITE(5001,'(5ES15.6)') x1(j),mom_z(j,1),mom_z(j,2),
     x           mom_z(j,3),mom_z(j,4)
          END IF
        END DO
        CLOSE(5001)
      END IF
c-----
      END SUBROUTINE mduc_soot_cond_1d
c-----    
c=======================================================================
      SUBROUTINE mduc_soot_write_integral(mat)
c=======================================================================
c     
c     Purpose:
c       write the integral soot values from the domain to a file
c    
c     Variables:
c       input  - mat            material domain (from FIRE)
c                 
c=======================================================================
      USE prec_mod,  ONLY : prec
      USE comm0,     ONLY : ncell,nsp,nep,iampro,time,itst,crank
      USE comm1,     ONLY : den,vol,actphi,pasphi
      USE mduc
c-----
      IMPLICIT NONE
c-----
c----- declare function arguments
      INTEGER,    INTENT(IN) :: mat
c----- declare local variables
      INTEGER    :: i,j,nc
      LOGICAL    :: fExist
      REAL(prec) :: mass,sootN,sootVol,sootArea,totVol,totDens
      REAL(prec) :: sootMass,sootDens
      INTEGER, PARAMETER :: outID = 6, outW = 80
      CHARACTER(len=256) :: text1,text2,sootFile
c-----------------------------------------------------------------------
      sootDens = 1800.0
      sootN    = 0.0
      sootVol  = 0.0
      sootArea = 0.0
      totVol   = 0.0
      DO nc=nsp(mat),nep(mat)
        mass = den(nc)*vol(nc)
c       update the soot integral values
        pasphi(nc+iSOOT_N)    = vol(nc)*actphi(nc+iSOOT_M00)
        pasphi(nc+iSOOT_VOL)  = vol(nc)*actphi(nc+iSOOT_M10)
        pasphi(nc+iSOOT_AREA) = vol(nc)*actphi(nc+iSOOT_M01)
c       add to integral
        sootN    = sootN    + pasphi(nc+iSOOT_N)
        sootVol  = sootVol  + pasphi(nc+iSOOT_VOL)
        sootArea = sootArea + pasphi(nc+iSOOT_AREA)
        totVol   = totVol   + vol(nc)
      END DO
c---- reduce variables across processors
      CALL dglsum(sootN)
      CALL dglsum(sootVol)
      CALL dglsum(sootArea)
      CALL dglsum(totVol)
c---- 
      sootMass = sootDens*sootVol
      totDens  = sootMass/totVol
c---- append to soot file
      sootFile = 'soot_integral.dat'
      IF (iampro .EQ. rootpro) THEN
        INQUIRE(file=sootFile,exist=fExist)
c---- write header if file doesn't exist
        IF (.NOT.fExist) THEN
          OPEN(5001,file=sootFile,status='unknown',position='append')
          WRITE(5001,'(3A10,6A20)') '# nt<1>', 'time[s]<2>',
     x         'CA[deg]<3>','cyl_vol[m^3]<4>',
     x         'number[-]<5>', 'volume[m^3]<6>', 'area[m^2]<7>',
     x         'mass[kg]<8>', 'density[kg/m^3]<9>'
          CLOSE(5001)
        END IF
        OPEN(5001,file=sootFile,status='unknown',position='append')
        WRITE(5001,'(I5,8E15.6)') itst,time,crank,totVol,
     x       sootN,sootVol,sootArea,
     x       sootMass,totDens
        CLOSE(5001)
c---- write summary to .fla file
        CALL mduc_out_header(outW,"Integral Soot Values")
        CALL mduc_out_real(outW,outID,'Soot Number','-',sootN,1)
        CALL mduc_out_real(outW,outID,'Soot Volume','m^3',sootVol,1)
        CALL mduc_out_real(outW,outID,'Soot Area','m^2',sootArea,1)
        CALL mduc_out_real(outW,outID,'Soot Mass','mg',sootMass*1e6,1)
        WRITE(text1, '(A)') REPEAT('*', outW-2)
        CALL mduc_out(outW, text1)
        WRITE (*,*) ''
      END IF
c-----
      END SUBROUTINE mduc_soot_write_integral
c-----
c-----------------------------------------------------------------------
