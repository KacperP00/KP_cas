c=======================================================================
      SUBROUTINE mduc_useout(mph)
c=======================================================================
c
c     USEOUT IS A SPECIAL PURPOSE ROUTINE FOR GENERATING USER
c     DEFINED OUTPUT OR FOR CHANGING THE RUN PARAMETERS
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AT THE END OF EVERY TIME-STEP
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c=======================================================================
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : i_useout,ncell,nsp,nep,pref,pref_old,
     x                     time,timo,dt,itst,iter,nreg,numpro,iampro
      USE comm1,    ONLY : p,po,vol,volo,pasphi,pasphb,hs,hsb,
     x                     ibc,ner,nsr
      USE mduc
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mph
c-----
      INTEGER :: mat,nc,istat,isc,ir,ib
      INTEGER :: nZ1,nZ2,out_w,outID,inert
c-----
      REAL(prec) :: Pavg,totVol
      REAL(prec) :: Z1min,Z1max,Z2min,Z2max
      REAL(prec), DIMENSION(:), ALLOCATABLE :: Z1,Z2
c-----
      CHARACTER(len=256) :: text1, text2
      CHARACTER(len=25), DIMENSION(2) :: fuelName2
c-----
c-----------------------------------------------------------------------
c-----
      mat = 1
      IF(I_USEOUT == 1) THEN          ! activation of useout in the GUI (= 1)
c-----------------------------------------------------------------------
c----- calculate the average thermodynamic pressure
c-----------------------------------------------------------------------
c     TODO: test whether we can get the old pressure here to avoid
c           storage of Pavg_old
        Pavg   = 0.0
        totVol = 0.0
        DO nc=nsp(mat),nep(mat)
          Pavg   = Pavg   + (pref(mat)+p(nc))*vol(nc)
          totVol = totVol + vol(nc)
        END DO
c     reduce the variables across all processors
        CALL dglsum(Pavg)
        CALL dglsum(totVol)
c     finish the normalization
        Pavg = Pavg/totVol
c     compute the average pressure rate
        dPdt = (Pavg-Pavg_old)/dt
c----- output the pressure data
c----- set min/max extents in domain
c     TODO: for now assume always between 0 and 1
        Z1min = 0.0
        Z1max = 1.0
c----- store the FIRE internal enthalpy in a passive scalar to output it
c     TODO: this should only be temporary
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
c-----------------------------------------------------------------------
c-----  check for the injection status
c-----------------------------------------------------------------------
        CALL mduc_check_injection(mat)
c-----
        outID = 1
        out_w = 80
        IF (iampro .LE. 1) THEN
          CALL mduc_out_header(out_w,'Flamelet Input Summary')
          CALL mduc_out_int(out_w,outID,'Mixture detected',mduc_mix)
          CALL mduc_out_int(out_w,outID,'Current injection',mduc_inj)
          CALL mduc_out_real(out_w,outID,'Pressure','Pa',Pavg,1)
          CALL mduc_out_real(out_w,outID,'Pressure Rate','Pa/s',dPdt,1)
          WRITE(text1,'(A)') repeat('*',out_w-2)
          CALL mduc_out(out_w,text1)
        END IF
c-----------------------------------------------------------------------
c-----  call  the appropriate solver based on the current status
c-----------------------------------------------------------------------
        IF (mduc_inj == 10) THEN     ! no mixture fraction present
          inert = 1
          CALL mduc_rif_1d(mat,Pavg,Z1min,Z1max,inert)
        ELSEIF (mduc_inj == 11) THEN ! pilot injection
          inert = 0
          CALL mduc_rif_1d(mat,Pavg,Z1min,Z1max,inert)
        ELSEIF (mduc_inj >= 20) THEN ! main injection
c         initialize 2D flamelet if necessary
          IF (MDUCgetNdim(mduc_mem) .EQ. 1 .AND. mduc_inj .NE. 21) THEN
            IF (iampro .LE. 1) THEN
              CALL mduc_out_info
     x             ('Initializing a 2D flamelet from a 1D flamelet')
              WRITE(*,*) 'f_init nt: ', itst
            END IF
c     TODO: make 2D initialization more general
c           for now, must be specified in MDUC input file
            nZ2 = -1
            CALL MDUCinit2Dfrom_BC_and_1D(nZ2,-1,2,fuelName2,%VAL(0),
     x           -1.0,mduc_mem)
c           reinitialize the parallel settings
            CALL MDUCsetParallel(MPI_COMM_MDUC, mduc_mem)
c           update the zero point for the grids
            nZ1 = MDUCgetGridSize(1, mduc_mem)
            nZ2 = MDUCgetGridSize(2, mduc_mem)
            ALLOCATE(Z1(nZ1), STAT=istat); Z1 = 0.0
            ALLOCATE(Z2(nZ2), STAT=istat); Z2 = 0.0
            CALL MDUCgetGrid2D(Z1, Z2, mduc_mem)
            Z1_zero = 0.5*(Z1(2)-Z1(1))
            Z2_zero = 0.5*(Z2(2)-Z2(1))
            DEALLOCATE(Z1); DEALLOCATE(Z2)
c           write the initial solution
c           TODO: recompute total enthalpy before doing so
            IF (iampro .LE. 1) THEN
              CALL MDUCwriteFlamelet('mduc_initial_2D.rif',mduc_mem)
            END IF
          END IF
          IF (mduc_inj == 20) THEN
            CALL mduc_inert_2d(mat,Pavg)
          ELSE IF (mduc_inj == 21) THEN
c           TODO: make more general, assumes that x2 grid is same as x1
            CALL mduc_rif_21d(mat,Pavg,Z1min,Z1max)
          ELSE IF (mduc_inj == 22) THEN
            CALL mduc_rif_2d(mat,Pavg,Z1min,Z1max,Z2min,Z2max)
          END IF
        END IF
c-----
c----- store the average pressure for use in the dP/dt term
        Pavg_old = Pavg
c----- update the soot values if necessary
        IF (mduc_soot .GT. 0) THEN
          CALL mduc_soot_write_integral(mat)
        END IF
c----- make sure that the temperature is up-to-date with the current enthalpy
        CALL use_cthtmp
      END IF
c-----
      END SUBROUTINE mduc_useout
