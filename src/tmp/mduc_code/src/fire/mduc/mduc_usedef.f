c=======================================================================
      SUBROUTINE mduc_usedef(mat,mph)
c=======================================================================
c
c     USEDEF IS A SPECIAL PURPOSE ROUTINE TO MODIFY DEFINITION
c     OF SOME VARIABLES AND CONSTANTS
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AT PROGRAM START/RESTART 
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c-----------------------------------------------------------------------
c-----
      USE comm0,   ONLY : num_actscl,num_passcl,ispray,iampro
      USE mduc,    ONLY : mduc_actoff,mduc_pasoff,mduc_soot,mduc_nInj,
     x                    mduc_all_species,mduc_outSpraySrc,
     x                    mduc_outSootSrc,mduc_outSootTP,nSootMom
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(IN) :: mat,mph
c-----------------------------------------------------------------------
c-----
c----- read the parameters from the csf file
      CALL mduc_csf_read('csf')
c----- store the number of active scalars from other libraries
      mduc_actoff = num_actscl
      mduc_pasoff = num_passcl
c----- add to number of active/passive scalar arrays
c      see mduc_useini.f for definitions
      num_actscl = num_actscl + 3
      num_passcl = num_passcl + 12
      IF (mduc_outSpraySrc .EQ. 1) THEN
        IF (ispray .EQ. 0) THEN
          WRITE(*,*) 'Warning: no spray'
          WRITE(*,*) '         will not output spray source terms'
        ELSE
          num_passcl = num_passcl + 3
        END IF
      END IF
c----- add split-injection passive scalars
      IF (mduc_nInj .EQ. 2) THEN
        num_actscl = num_actscl + 3
        num_passcl = num_passcl + 14
        IF (mduc_outSpraySrc .EQ. 1 .AND. ispray .GT. 0) THEN
          num_passcl = num_passcl + 6
        END IF
      END IF
c----- add soot equations      
      IF (mduc_soot .GT. 0) THEN
c       add first order moments
        nSootMom = 4
c         add second order moments
        IF (mduc_soot .GT. 1) THEN
          nSootMom = nSootMom + 3
        END IF
        num_actscl = num_actscl + nSootMom
        num_passcl = num_passcl + 3
      ELSE
        mduc_outSootSrc = 0
        mduc_outSootTP  = 0
      END IF
      IF (mduc_outSootSrc .GT. 0) THEN
        num_passcl = num_passcl + nSootMom
      END IF
      IF (mduc_outSootTP .GT. 0) THEN
        num_passcl = num_passcl + nSootMom
      END IF
c----- if using a subset of species, additional storage for Cp
      IF (mduc_all_species .EQ. 0) THEN
        num_passcl = num_passcl + 3
      END IF
c----- write summary
      IF (iampro < 2) THEN
        WRITE(*,*) '************* MDUC usedef summary **************'
        WRITE(*,*) '            number of injections: ', mduc_nInj
        IF (mduc_soot .GT. 0) THEN
          WRITE(*,*) '               soot moment order: ', mduc_soot
        END IF
        WRITE(*,*) ''
        WRITE(*,*) '   active scalars - total number: ', num_actscl
        WRITE(*,*) '                     mduc offset: ', mduc_actoff
        WRITE(*,*) '  passive scalars - total number: ', num_passcl
        WRITE(*,*) '                     mduc offset: ', mduc_pasoff
        WRITE(*,*) '************************************************'
      END IF
c-----
      RETURN
c-----
      END SUBROUTINE mduc_usedef
c-----
c=======================================================================
      SUBROUTINE mduc_csf_read(csf_filename)
c=======================================================================
c
c     Function to initialize indices corresponding to the offsets of
c     passive scalars
c-----
c-----------------------------------------------------------------------
      USE comm0, ONLY : iampro
      USE mduc
      IMPLICIT NONE
c-----
      INTEGER :: csf_numpara
      PARAMETER(csf_numpara=13)
c-----
      CHARACTER(len=20), INTENT(in) :: csf_filename ! filename of CSF
      CHARACTER(len=20) :: csf_segment='[MDUC]'     ! name of segment
      CHARACTER(len=20), DIMENSION(csf_numpara) ::
     x     csf_value, csf_content ! values and their description
      LOGICAL :: csf_status                     ! status
      INTEGER, DIMENSION(csf_numpara) :: icsf_f
c-----
      INTEGER            :: i, w_sum
      CHARACTER(len=256) :: text1,text2
c-----
c-----------------------------------------------------------------------
c-----
c----- define the keywords for the csf file
      csf_content(1)='inputFile'
      csf_content(2)='tPilot'
      csf_content(3)='tMain'
      csf_content(4)='enthalpy'
      csf_content(5)='satData'
      csf_content(6)='soot'
      csf_content(7)='sootSrcUpdate'
      csf_content(8)='species'
      csf_content(9)='nInjection'
      csf_content(10)='writeSpraySource'
      csf_content(11)='writeSootSource'
      csf_content(12)='writeSootTP'
      csf_content(13)='outputDebug'
c----- set the default values
      mduc_infile      = 'mduc.inp'
      H_method_name    = 'total'
      H_method         = 0
      mduc_satData     = 'saturation.dat'
      mduc_tpilot      = -1.0
      mduc_tmain       = -1.0
      mduc_nInj       = 1
      mduc_soot        = 0
      sootSrcUpdateDef = 'iteration'
      mduc_all_species = 0
      mduc_outSpraySrc = 0
      mduc_outSootSrc  = 0
      mduc_outSootTP   = 0
c----- read in the segment of the csf file
      CALL csf_read (csf_filename, csf_segment, csf_numpara,
     x     csf_content, csf_value, icsf_f, csf_status)
      IF (.NOT.csf_status) THEN
         WRITE(*,*) 'Error 1: mduc_csf_read -> csf not found!'
         STOP
      END IF
c----- assign the values to the appropriate flags
      IF (icsf_f(1).EQ.1)  READ(csf_value(1),  '(a)') mduc_infile
      IF (icsf_f(2).EQ.1)  READ(csf_value(2),  '(f)') mduc_tpilot 
      IF (icsf_f(3).EQ.1)  READ(csf_value(3),  '(f)') mduc_tmain
      IF (icsf_f(4).EQ.1)  READ(csf_value(4),  '(a)') H_method_name
      IF (icsf_f(5).EQ.1)  READ(csf_value(5),  '(a)') mduc_satData
      IF (icsf_f(6).EQ.1)  READ(csf_value(6),  '(i)') mduc_soot
      IF (icsf_f(7).EQ.1)  READ(csf_value(7),  '(a)') sootSrcUpdateDef
      IF (icsf_f(8).EQ.1)  READ(csf_value(8),  '(i)') mduc_all_species
      IF (icsf_f(9).EQ.1)  READ(csf_value(9),  '(i)') mduc_nInj
      IF (icsf_f(10).EQ.1) READ(csf_value(10), '(i)') mduc_outSpraySrc
      IF (icsf_f(11).EQ.1) READ(csf_value(11), '(i)') mduc_outSootSrc
      IF (icsf_f(12).EQ.1) READ(csf_value(12), '(i)') mduc_outSootTP
      IF (icsf_f(13).EQ.1) READ(csf_value(13), '(i)') mduc_outDebug
c----- write summary of input read
      IF (iampro < 2) THEN
        w_sum = 80
        WRITE(*,*) ''
        CALL mduc_out_header(w_sum, "MDUC csf read")
        WRITE(text2, '(A)') REPEAT(' ', w_sum-2)
        CALL mduc_out(w_sum, text2)
c        
        DO i=1,csf_numpara
          IF (icsf_f(i).EQ.1) THEN
            WRITE(text1, '(A20,A,A20)') trim(csf_content(i)), ' : ',
     x           trim(csf_value(i))
            CALL mduc_out(w_sum, text1)
          END IF
          
        END DO
c     
        CALL mduc_out(w_sum, text2)
        WRITE(text1, '(A)') REPEAT('*', w_sum-2)
        CALL mduc_out(w_sum, text1)
      END IF
c-----
      END SUBROUTINE mduc_csf_read
c-----------------------------------------------------------------------
