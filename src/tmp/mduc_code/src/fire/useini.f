c=======================================================================
      SUBROUTINE useini(mat,mph)
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
      USE comm0, ONLY : i_useini
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mat,mph
c-----------------------------------------------------------------------
c-----  activation of useini (=1) in the GUI 
      IF (I_USEINI .GT. 0) THEN
        CALL mduc_useini(mat,mph)
      END IF
c-----
      END SUBROUTINE useini
c-----------------------------------------------------------------------
c-----
c=======================================================================
      SUBROUTINE csf_read(csf_filename, csf_segment, csf_numpara,
     x     csf_content, csf_value, i_csf_found, csf_status)
c=======================================================================
c
c     Common Steering File - read-out subroutine
c
c       input  : csf_filename ..... string containing file to be read
c                csf_segment  ..... identifier of segment within file
c                csf_numpara  ..... number of possible parameters
c       output : csf_content  ..... description of variables
c                csf_value    ..... value assigned to variable
c                i_csf_found  ..... number of variables found
c                csf_status   ..... 
c
c-----------------------------------------------------------------------
      USE comm0, ONLY : iampro
      IMPLICIT NONE
c-----      
      INTEGER,           INTENT(in) :: csf_numpara
      CHARACTER(len=20), INTENT(in) :: csf_filename, csf_segment
c-----
      CHARACTER(len=20), DIMENSION (csf_numpara) :: csf_content
      CHARACTER(len=20), DIMENSION (csf_numpara) :: csf_value
      INTEGER          , DIMENSION (csf_numpara) :: i_csf_found !out
      LOGICAL :: csf_status 
c-----
      INTEGER :: i,j
      INTEGER :: slength        ! length of string
      CHARACTER(len=8) :: intro='# csf : '
      CHARACTER(len=20) :: string1, string2
      CHARACTER(len=40) :: stringtemp
      LOGICAL :: fexist, csf_segment_found=.FALSE.
      LOGICAL, DIMENSION (csf_numpara) :: csf_value_found
c-----
c-----------------------------------------------------------------------
c-----
      DO i=1,csf_numpara
        csf_value_found (i)=.FALSE.
        csf_value (i)='0'
        i_csf_found = 0
      END DO
c-----
      INQUIRE (file=csf_filename, exist=fexist)
      IF (.NOT.fexist) THEN
        IF (iampro < 2) THEN
          WRITE (*,'(4a)') intro, 'common steering file "',
     x         csf_filename, '" does not exist'
        END IF
        csf_status=.FALSE.
        RETURN
      END IF
c----- open the file
      OPEN(728, file=csf_filename, status='old',form='formatted')
c----- search until segment name is found
      DO
        READ (728,'(a20)',END=5411) stringtemp
        slength=len_trim(csf_segment)
        IF (stringtemp(1:slength) .EQ. csf_segment(1:slength)) THEN
          csf_segment_found=.TRUE.
          EXIT
        END IF
      END DO
 5411 CONTINUE
c-----
      csf_status=csf_segment_found
      IF (.NOT.csf_segment_found) THEN
        IF (iampro < 2) THEN
          WRITE (*,'(4a)') intro, 'CSF segment ',
     x         csf_segment,' not found'
        END IF
        RETURN
      END IF
      IF (iampro < 2) THEN
        WRITE (*,'(3a)') intro, 'found segment ', csf_segment
        WRITE (*,'(2a)') intro, 'reading values ...'
      END IF
c-----
c----- read contents of segment until EOF or next segments starts
      DO
        READ (728,'(a40)',END=5413) stringtemp
        IF (stringtemp (1:1).EQ.'[') EXIT
c----- read value (until blank or tab) into string1
        i=1
        string1=''
        DO
          IF ( (i.LE.20) .AND.
     x         (stringtemp (i:i).NE.char(9)) .AND.
     x         (stringtemp (i:i).NE.char(32)) ) THEN
            string1 (i:i)=stringtemp (i:i)
            i=i+1
            CYCLE
          ELSE
            EXIT
          END IF
        END DO
c----- throw away blanks and tabs until next letter and write the
c      result to string2
        DO
          IF ( (i.LE.20) .AND. (
     x         (stringtemp (i:i).EQ.char(9)) .OR.
     x         (stringtemp (i:i).EQ.char(32)) )) THEN
            i=i+1
            CYCLE
          ELSE
            EXIT
          END IF
        END DO
        string2=stringtemp (i:i+19)
c-----  blank out string2 after first tab or space
        i=1
        DO
          IF ( (i.LE.20) .AND.
     x         (string2 (i:i).NE.char(9)) .AND.
     x         (string2 (i:i).NE.char(32)) ) THEN
            i=i+1
            CYCLE
          ELSE
            EXIT
          END IF
        END DO
        slength=i
        DO j=i,20
          string2 (j:j)=' '
        END DO
c----- check for equal description within set of given descriptions
        DO i=1,csf_numpara
          IF (string2.EQ.csf_content (i)) THEN
            csf_value (i)=string1
            IF (iampro < 2) THEN
              WRITE (*,'(3a)') intro, string1, string2
            END IF
            csf_value_found (i)=.TRUE.
            i_csf_found(i) = 1  !RSk
            EXIT
          END IF
        END DO
      END DO                    ! infinite loop reading lines
 5413 CONTINUE
c----- check if all listed variables are found
      DO i=1,csf_numpara
        IF (.NOT.csf_value_found (i)) THEN
          slength=len_trim (csf_content (i))
          IF (iampro < 2) THEN
            WRITE (*,'(5a)') intro, 'CSF variable "',
     x           csf_content (i)(1:slength),
     x           '" not found in segment ', csf_segment
          END IF
        END IF
      END DO
c-----
      CLOSE(728)
c-----
      END SUBROUTINE csf_read
c-----
