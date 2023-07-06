c=======================================================================
      SUBROUTINE mduc_usepsc(iconv,mph,kpsc)
c=======================================================================
c
c     USEPSC IS A SPECIAL PURPOSE ROUTINE TO MODIFY SOURCE TERM IN
c     SCALAR EQUATION
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AFTER EACH ITERATION
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c-----------------------------------------------------------------------
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : ncell,nsp,nep,num_actscl,num_passcl,dtin,
     x                     pref,pref_old,dt,iter,ispray,time,timo,
     x                     nsb,neb,nreg,nbfac,iampro
      USE comm1,    ONLY : actphi,actphb,actphr,pasphi,pasphr,su1,sp1,
     x                     isol_actscl,isol_passcl,vol,volo,den,denb,
     x                     ed,edb,p,po,te,teb,vit,vitb,nsr,ner,lb
      USE cydrnoz,  ONLY : ssmass,ssmasc,ieulag,senerg,nnozzl,nprops,
     x                     nspsca,isosca
      USE cthmod,   ONLY : numscl,ymfrac,ymfrao,ymfrar,index_fuel
      USE mduc
      USE mduc_sdm
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(IN) :: iconv,mph,kpsc
      INTEGER    :: mat,nc,istat,prop,ii,ir,ib,isc,gZ2i_0,gZ2i_1,
     x              gZ2b_0,gZ2b_1
      REAL(prec) :: C_prod   ! production coefficient for variance
      REAL(prec) :: C_eps    ! dissipation coefficient for variance
      REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: grad_Z0,grad_Z1,grad_Z2
      REAL(prec) :: mag2_gZ0, mag2_gZ1, mag2_gZ2
c-----------------------------------------------------------------------
c-----
      mat  = 1       ! number of materials
      Sc_t = 0.7     ! defined in mduc_sdm
c-----
      
c----- if solving with a spray, the variance source terms from the 
c      SDM model must be determined and stored before the source terms
c      since they are stored in passive scalars
      IF (num_passcl > mduc_pasoff) THEN
        IF (ispray == 1) THEN
c          CALL SDMcalcSrc(mat)
        END IF
      END IF
c      write(*,*) iampro, ': in mduc_usepsc'
c-----
c-----------------------------------------------------------------------
c----- calculation of scalar source terms of the mixture fraction
c-----------------------------------------------------------------------
      IF (num_actscl > mduc_actoff) THEN   ! number of active scalar > 0
c----- compute the gradients of the mixture fraction for use in variance
        ALLOCATE(grad_Z1(3,ncell), STAT=istat); grad_Z1 = 0.0
        CALL gradfi(actphi(1:ncell),actphb(1:nbfac),grad_Z1,1,1,mat)
        IF (mduc_nInj .EQ. 2) THEN
          ALLOCATE(grad_Z0(3,ncell), STAT=istat); grad_Z0 = 0.0
          ALLOCATE(grad_Z2(3,ncell), STAT=istat); grad_Z2 = 0.0
          gZ2i_0 = 1+iZ2MEAN
          gZ2i_1 = gZ2i_0+ncell-1
          gZ2b_0 = 1+bZ2MEAN
          gZ2b_1 = gZ2b_0+nbfac-1
          CALL gradfi(actphi(gZ2i_0:gZ2i_1),actphb(gZ2b_0:gZ2b_1),
     x         grad_Z2,1,1,mat)
          DO nc=nsp(mat),nep(mat)
            DO ii=1,3
              grad_Z0(ii,nc) = -( grad_Z1(ii,nc) + grad_Z2(ii,nc) )
            END DO
          END DO
        END IF
c-----
        IF (ispray == 0) THEN                           ! no spray
          DO nc=nsp(mat),nep(mat)
c----- coefficients
            C_prod= 2.0*vit(nc)/Sc_t                   ! production 
            C_eps = 2.0*den(nc)*ed(nc)/(te(nc)+1.0e-6) ! dissipation
c----- mixture fraction mean sources
            actphr(nc+iZ1MEAN)= 0.0
c----- mixture fraction variance sources
c      compute the magnitudes of dZ/dx
            mag2_gZ1 = DOT_PRODUCT(grad_Z1(:,nc),grad_Z1(:,nc))
c      variance sources
            actphr(nc+iZ1VAR)= C_prod*mag2_gZ1-C_eps*actphi(nc+iZ1VAR)
c      total enthalpy
            actphr(nc+iHTOT) = dPdt
            IF (mduc_nInj .EQ. 2) THEN
              actphr(nc+iZ2MEAN)= 0.0
c             compute the magnitudes of dZ_i/dx
              mag2_gZ0 = DOT_PRODUCT(grad_Z0(:,nc),grad_Z0(:,nc))
              mag2_gZ2 = DOT_PRODUCT(grad_Z2(:,nc),grad_Z2(:,nc))
c             variance sources
              actphr(nc+iZ0VAR)= C_prod*mag2_gZ0-C_eps*actphi(nc+iZ0VAR)
              actphr(nc+iZ2VAR)= C_prod*mag2_gZ2-C_eps*actphi(nc+iZ2VAR)
            END IF
          END DO
        ELSE IF (ispray == 1) THEN ! spray activated
c----- compute the variance source terms using Single Droplet Model (SDM)
c      of Reveillon and Vervisch
          CALL SDMcalcSrc(mat)
c-----
          DO nc=nsp(mat),nep(mat)
c----- coefficients
            C_prod= 2.0*vit(nc)/Sc_t                   ! production coeff
c----- mixture fraction mean sources
            IF ( kpsc > (numscl+200) .AND.
     x           kpsc <= ((numscl+200)+num_actscl)) THEN
              isc = kpsc-(numscl+200) ! index of active scalar
              IF (isc .EQ. mduc_actoff+1) THEN
                sp1(nc) = sp1(nc) + ssmass(nc)*dtin*ieulag(4)
                su1(nc) = su1(nc) + pasphi(nc+iMEVAP1)*dtin*ieulag(4)
              ELSE IF (mduc_nInj .EQ. 2 .AND.
     x                 isc .EQ. mduc_actoff+4) THEN
                sp1(nc) = sp1(nc) + ssmass(nc)*dtin*ieulag(4)
                su1(nc) = su1(nc) + pasphi(nc+iMEVAP2)*dtin*ieulag(4)
              END IF
            END IF
c----- compute the magnitudes of dZ/dx
            mag2_gZ1 = DOT_PRODUCT(grad_Z1(:,nc),grad_Z1(:,nc))
c----- mixture fraction variance sources
            actphr(nc+iZ1VAR) = C_prod*mag2_gZ1                 ! production
     x           + den(nc)*pasphi(nc+iZ1VSRC)                   ! evaporation
     x           - den(nc)*pasphi(nc+iZ1VEPS)*actphi(nc+iZ1VAR) ! dissipation
c-----
            actphr(nc+iHTOT) = dPdt                             ! pressure
     x           + senerg(nc)*dtin*den(nc)                      ! evaporation
     x           + H_vap_ref*ssmass(nc)*dtin*ieulag(4)/vol(nc)  ! species
            IF (mduc_outSpraySrc) THEN
              pasphi(nc+iZ1VPR) =  C_prod*mag2_gZ1
              pasphi(nc+iZ1VEV) =  den(nc)*pasphi(nc+iZ1VSRC)
              pasphi(nc+iZ1VDI) = -den(nc)*pasphi(nc+iZ1VEPS)*
     x                                     actphi(nc+iZ1VAR)
            END IF
            IF (mduc_nInj .EQ. 2) THEN
c----- compute the magnitudes of dZ_i/dx
              mag2_gZ0 = DOT_PRODUCT(grad_Z0(:,nc),grad_Z0(:,nc))
              mag2_gZ2 = DOT_PRODUCT(grad_Z2(:,nc),grad_Z2(:,nc))
c----- mixture fraction variance sources
              actphr(nc+iZ0VAR) = C_prod*mag2_gZ0                 ! production
     x             + den(nc)*pasphi(nc+iZ0VSRC)                   ! evaporation
     x             - den(nc)*pasphi(nc+iZ0VEPS)*actphi(nc+iZ0VAR) ! dissipation
c-----
              actphr(nc+iZ2VAR) = C_prod*mag2_gZ2                 ! production
     x             + den(nc)*pasphi(nc+iZ2VSRC)                   ! evaporation
     x             - den(nc)*pasphi(nc+iZ2VEPS)*actphi(nc+iZ2VAR) ! dissipation
              IF (mduc_outSpraySrc) THEN
c----- these are for debugging the SDM module
                pasphi(nc+iZ0VPR) =  C_prod*mag2_gZ0
                pasphi(nc+iZ2VPR) =  C_prod*mag2_gZ2
                pasphi(nc+iZ0VEV) =  den(nc)*pasphi(nc+iZ0VSRC)
                pasphi(nc+iZ2VEV) =  den(nc)*pasphi(nc+iZ2VSRC)
                pasphi(nc+iZ0VDI) = -den(nc)*pasphi(nc+iZ0VEPS)*
     x                                       actphi(nc+iZ0VAR)
                pasphi(nc+iZ2VDI) = -den(nc)*pasphi(nc+iZ2VEPS)*
     x                                       actphi(nc+iZ2VAR)
              END IF
            END IF
          END DO
c----- set the boundary faces to the cell centred values
c$$$          DO ir=0,nreg
c$$$            DO ib=nsr(ir),ner(ir)
c$$$              nc = lb(ib)
c$$$              actphb(ib+bZ1MEAN) = actphi(nc+iZ2MEAN)
c$$$              actphb(ib+bZ2MEAN) = actphi(nc+iZ2MEAN)
c$$$              actphb(ib+bZ0VAR)  = actphi(nc+iZ0VAR)
c$$$              actphb(ib+bZ1VAR)  = actphi(nc+iZ1VAR)
c$$$              actphb(ib+bZ2VAR)  = actphi(nc+iZ2VAR)
c$$$              actphb(ib+bHTOT)   = actphi(nc+iHTOT)
c$$$            END DO
c$$$          END DO
c-----
c-----  Update the boundary cells 
c----- 
c$$$          DO nb=nsp(mat),nep(mat)
c$$$c----- coefficients
c$$$            C_prod= 2.0*vitb(nb)/Sc_t                   ! production
c$$$            C_eps = 2.0*denb(nb)*edb(nb)/(teb(nb)+1.0e-6) ! dissipation
c$$$c----- mixture fraction mean sources
c$$$            actphb(nb+bZ1MEAN)= ssmasc(nc+0*ncell)*dtin/vol(nc)
c$$$            actphb(nb+bZ2MEAN)= ssmasc(nc+1*ncell)*dtin/vol(nc)
c$$$c----- Z1 mean gradient squared
c$$$c----- mixture fraction variance sources
c$$$c      compute the magnitudes of dZ/dx
c$$$            mag2_gZ1=grad_Z1(1,nc)**2+grad_Z1(2,nc)**2+grad_Z1(3,nc)**2
c$$$            mag2_gZ2=grad_Z2(1,nc)**2+grad_Z2(2,nc)**2+grad_Z2(3,nc)**2
c$$$            mag2_gZ0 = (grad_Z1(1,nc)+grad_Z2(1,nc))**2 +
c$$$     x                 (grad_Z1(2,nc)+grad_Z2(2,nc))**2 +
c$$$     x                 (grad_Z1(3,nc)+grad_Z2(3,nc))**2
c$$$c----- variance sources
c$$$            actphb(nb+bZ0VAR)= C_prod*mag2_gZ0                  ! production
c$$$     x           + denb(nb)*pasphb(nb+bZ0VSRC)                   ! dissipation
c$$$     x           - denb(nb)*actphb(nb+bZ0VAR)*pasphb(nb+bZ0VEPS) ! evaporation
c$$$c-----
c$$$            actphb(nb+bZ1VAR)= 
c$$$     x           C_prod*mag2_gZ1 - C_eps*actphb(nb+bZ1VAR)
c$$$c-----
c$$$            actphr(nb+bZ2VAR)= C_prod*mag2_gZ2                  ! production
c$$$     x           + denb(nb)*pasphb(nb+bZ2VSRC)                   ! dissipation
c$$$     x           - denb(nb)*actphb(nb+bZ2VAR)*pasphb(nb+bZ2VEPS) ! evaporation
c$$$          END DO
        END IF
        IF(ALLOCATED(grad_Z0))  DEALLOCATE(grad_Z0)
        IF(ALLOCATED(grad_Z1))  DEALLOCATE(grad_Z1)
        IF(ALLOCATED(grad_Z2))  DEALLOCATE(grad_Z2)
c-----
c-----  compute the equations necessary for the soot model
        IF (mduc_soot .GT. 0) THEN
          CALL mduc_soot_psc(mat)
        END IF
c-----  end soot
c-----
      END IF
c-----
      IF(num_passcl > mduc_pasoff) THEN  ! number of passive scalar > 0
      END IF
c-----
      RETURN
c-----
      END SUBROUTINE mduc_usepsc
