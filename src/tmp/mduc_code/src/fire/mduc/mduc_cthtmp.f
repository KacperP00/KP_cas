c=======================================================================
      SUBROUTINE mduc_cthtmp
c=======================================================================
c
c     description: user function for calculating the temperature 
c                  of gas mixture
c
c-----------------------------------------------------------------------
c
c     author:            Peter Priesching
c     date:              10-10-01
c     project:           Fire v2010
c     project manager:   Peter Priesching
c
c     literature reference:
c
c-----------------------------------------------------------------------
c
c     modifications:
c
c     author:            Eric Doran
c     date:              continuing
c     purpose:
c     identification:
c
c-----------------------------------------------------------------------
c
c     argument list:
c
c-----------------------------------------------------------------------
c
c     subroutine is called by: cthdis
c
c-----------------------------------------------------------------------
c
c     subroutine calls:
c
c     name: cthucp            purpose: determine specific heat
c
c=======================================================================
c-----
      USE mduc, ONLY : H_method
      USE comm0, ONLY : iampro
c----- call the correct function based on total or sensible enthalpy
c      write(*,*) iampro, ': in function mduc_cthtmp'
      IF (H_method .EQ. 0) THEN
        CALL mduc_cthtmp_total
      ELSE IF (H_method .EQ. 1 .OR. H_method .EQ. 2) THEN
        CALL mduc_cthtmp_sensible
      END IF
      
      END SUBROUTINE mduc_cthtmp
c-----
c=======================================================================
      SUBROUTINE mduc_cthtmp_total
c=======================================================================
c
c     description: iterate for the temperature based on the total
c                  enthalpy of AVL, i.e.
c                     h_tot = cp*(T-Tref) + 0.5*u^2 + h_chem
c
c=======================================================================
c----- modules
      USE prec_mod, ONLY : prec
      USE comm0
      USE comm1,    ONLY : u,hs,t,ub,hsb,cpcofb,tb,ibc,nsr,ner,lb,
     x                     cpcof,pasphi,pasphb
      USE rot_mod,  ONLY : i_rot,r_rot,omega,abs_vel,id_sel_rot
      USE mduc
c-----
      IMPLICIT NONE
c----- local variables
      INTEGER    :: it,nc,mat,it_max,ir,ib,ir1,i
      REAL(prec) :: terr,usqr,temper,cm8
      REAL(prec) :: hsum, cpsum, deltaT, Ym, c_mass
      REAL(prec) :: cpmix, cpmix_ref, Href
      REAL(prec) :: hs_chem,weightmix,hs_comp
c-----------------------------------------------------------------------
c-----
      mat = 1
c-----
      cm8=1.0E-08
c-----
      it_max = 100
c-----
c-----------------------------------------------------------------------
c----- internal cells
c-----------------------------------------------------------------------
c-----
c      write(*,*) 'in cthtmp'
      lit1: DO it=1,it_max
c-----
        IF (I_USECPS >= 1) THEN
          CALL usecps(mat,msin)
        ELSE IF (mduc_all_species .EQ. 0) THEN
          DO nc=nsp(mat),nep(mat)
            cpcof(nc) = pasphi(nc+iCPMIX) + 
     x           (t(nc)-pasphi(nc+iTCP))*pasphi(nc+iDCPDT)
          END DO
        ELSE IF (mduc_all_species .EQ. 1) THEN
c         now updated in functions updatecfd_1d, updatecfd_2d
          CALL cthucp(1,ir)
        END IF
        terr=0.
c-----
        lib1: DO nc=nsp(mat),nep(mat)
          usqr   = DOT_PRODUCT(u(:,nc), u(:,nc))
          usqr   = usqr*iener
          temper = (hs(nc)-0.5*usqr-pasphi(nc+iHCHEM))/cpcof(nc)
     x             + H_Tref
          terr   = MAX(terr, ABS(t(nc)-temper)/t(nc))
          t(nc)  = temper
        END DO lib1
c-----
        CALL exchng(t,1,1)
        CALL dglmax(terr)
        IF(terr<cm8) EXIT lit1
c-----
      END DO lit1
c-----------------------------------------------------------------------
c----- boundary faces
c-----------------------------------------------------------------------
c-----
      lit2: DO it=1,it_max
c-----
        IF (I_USECPS >= 1) THEN
          CALL usecps(mat,msin)
        ELSE IF (mduc_all_species .EQ. 0) THEN
c       TODO: compute boundary cp here
        ELSE IF (mduc_all_species .EQ. 1) THEN
c         now updated in function mduc_update_boundY
          CALL cthucp(2,ir)
        END IF
        terr=0.
c-----
        lir2: DO ir=0,nreg
c-----limit to current material
          IF(ibc(2,ir)/=mat) cycle lir2
c-----exclude fixed temperature bc
          IF(.NOT.(ibc(1,ir)==4 .AND. ibc(5,ir)==3)) THEN
c-----
            lib2: DO ib=nsr(ir),ner(ir)
              usqr   = DOT_PRODUCT(ub(:,ib), ub(:,ib))
              usqr   = usqr*iener
              temper = (hsb(ib)-0.5*usqr-pasphb(ib+bHCHEM))/cpcofb(ib) 
     x                 + H_Tref
              terr   = MAX(terr, ABS(tb(ib)-temper)/tb(ib))
              tb(ib) = temper
            END DO lib2
c-----
          END IF
c-----
        END DO lir2
c-----
        CALL dglmax(terr)
        IF(terr<cm8) EXIT lit2
c-----
      END DO lit2
c-----
      END SUBROUTINE mduc_cthtmp_total
c-----
c=======================================================================
      SUBROUTINE mduc_cthtmp_sensible
c=======================================================================
c
c     description: iterate for temperature based on definition of total
c                  enthalpy as:
c                      h_tot = sum(Yi*hi)
c                  First, the total enthalpy must be reconstructed from 
c                  the sensible enthalpy variable in FIRE, i.e. hs
c
c=======================================================================
c----- modules
      USE prec_mod, ONLY : prec
      USE comm0
      USE comm1,    ONLY : u,ub,hs,hsb,cpcof,cpcofb,t,tb,ibc,nsr,ner,lb,
     x                     pasphi,pasphb,actphi,actphb
      USE rot_mod,  ONLY : i_rot,r_rot,omega,abs_vel,id_sel_rot
      USE cthmod,   ONLY : ymolfr,ymolfb,enthfo,numscl,wmolar,ymfrac,
     x                     ymfrab
      USE mduc
c-----
      IMPLICIT NONE
c----- local variables
      INTEGER    :: it,nc,mat,it_max,ir,ib,ir1,i
      REAL(prec) :: temper,usqr,cm8
      REAL(prec) :: h_tot,hsum,cpsum,deltaT,Ym
      REAL(prec), DIMENSION(:), ALLOCATABLE :: hi, cpi
c-----------------------------------------------------------------------
c-----
      mat = 1
c----- set convergence tolerance
      cm8 = 1.0E-08
c----- set limit on maximum iterations
      it_max = 1000
c-----
      IF(.NOT.ALLOCATED(hi))  ALLOCATE(hi(numscl))
      IF(.NOT.ALLOCATED(cpi)) ALLOCATE(cpi(numscl))

c-----------------------------------------------------------------------
c----- internal cells
c-----------------------------------------------------------------------
      lib1: DO nc=nsp(mat),nep(mat)
c-----  kinetic energy
        usqr=DOT_PRODUCT(u(:,nc), u(:,nc))
        usqr= usqr*iener
c-----  reconstruct total enthalpy 
        IF (H_method .eq. 1) THEN
          h_tot = hs(nc) - 0.5*usqr + pasphi(nc+iHCHEM)
        ELSE IF (H_method .eq. 2) THEN
          h_tot = actphi(nc+iHTOT)
        END IF
c-----  use current temperature as initial guess        
        temper = t(nc)
c-----  iterate to find temperature
        tfromh: DO it=1,it_max
          CALL MDUCcomputeThermoData(hi, cpi, temper)
          hsum  = 0.0
          cpsum = 0.0
          DO i=1,numscl
            IF (indexY(i) .GE. 0) THEN
              Ym    = ymfrac(nc+(i-1)*ncell)
              hsum  = hsum  +  hi(indexY(i))*Ym
              cpsum = cpsum + cpi(indexY(i))*Ym
            END IF
          END DO
          deltaT = -(hsum - h_tot)/cpsum
          temper = temper + deltaT
          IF (ABS(deltaT/temper) < cm8) EXIT tfromh
        END DO tfromh
        IF (it .GE. it_max) THEN
          print*, 'Error: temperature not converged'
        END IF
c-----  update the temperature in this cell       
        t(nc) = temper
c-----  update the specific heat consistent with the sensible enthalpy
c        cpcof(nc) = (hs(nc) - 0.5*usqr)/(t(nc)-H_Tref)
        END DO lib1
        CALL exchng(t,1,1)
c-----------------------------------------------------------------------
c----- boundary faces
c-----------------------------------------------------------------------
        lir2: DO ir=0,nreg
c-----    limit to current material
          IF(ibc(2,ir)/=mat) cycle lir2
c-----    exclude fixed temperature bc
          IF(.NOT.(ibc(1,ir)==4 .AND. ibc(5,ir)==3)) THEN
c-----
            lib2: DO ib=nsr(ir),ner(ir)
c-----        kinetic energy
              usqr=DOT_PRODUCT(ub(:,ib), ub(:,ib))
              usqr= usqr*iener
c-----        reconstruct total enthalpy 
              IF (H_method .eq. 1) THEN
                h_tot = hsb(ib) - 0.5*usqr + pasphb(ib+bHCHEM)
              ELSE IF (H_method .eq. 2) THEN
                h_tot = actphb(ib+bHTOT)
              END IF
c-----        use current temperature as initial guess        
              temper = tb(ib)
c-----        iterate to find temperature
              tfromhb: DO it=1,it_max
                CALL MDUCcomputeThermoData(hi, cpi, temper)
                hsum  = 0.0
                cpsum = 0.0
                DO i=1,numscl
                  IF (indexY(i) .GE. 0) THEN
                    Ym    = ymfrab(ib+(i-1)*nbfac)
                    hsum  = hsum  +  hi(indexY(i))*Ym
                    cpsum = cpsum + cpi(indexY(i))*Ym
                  END IF
                END DO
                deltaT = -(hsum - h_tot)/cpsum
                temper = temper + deltaT
                IF (ABS(deltaT/temper) < cm8) EXIT tfromhb
              END DO tfromhb
              IF (it .GE. it_max) THEN
                print*, 'Error: temperature not converged'
              END IF
c-----        update the temperature in this cell       
              tb(ib) = temper  
c-----  update the specific heat consistent with the sensible enthalpy
c              cpcofb(ib) = (hsb(ib) - 0.5*usqr)/(tb(ib)-H_Tref)
            END DO lib2
c-----
          END IF
c-----
        END DO lir2
        CALL exchng(tb,1,1)
c-----  update the specific heats based on the new temperature
        IF (I_USECPS >= 1) THEN
          CALL usecps(mat,msin)
        ELSE
          CALL cthucp(1,ir)
        END IF
c-----
c-----------------------------------------------------------------------
c-----
      IF(ALLOCATED(hi))  DEALLOCATE(hi)
      IF(ALLOCATED(cpi)) DEALLOCATE(cpi)
c-----------------------------------------------------------------------
c      print*,'exit mduc_cthtmp'
      RETURN
      END SUBROUTINE mduc_cthtmp_sensible
c-----
