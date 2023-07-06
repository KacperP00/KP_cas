c=======================================================================
      MODULE rootsolve
c=======================================================================
c
c     This module defines several root finding algorithms that are 
c     required to use the MDUC libraries in FIRE.
c     These functions are based on those found in NR
c
c     NOTE: this module must be compiled before the following user
c     routines that use it:
c        useini.f, usepsc.f, usetim.f
c     Put useroot.f before these in the list of routines compiled
c
c.....contact emdoran@stanford.edu
c-----
c=======================================================================
c-----
      USE prec_mod, ONLY : prec
      IMPLICIT NONE
c-----
      CONTAINS
c-----
      SUBROUTINE rootBrac(func,x1,x2,succes)
      implicit none
      REAL(prec), intent(inout) :: x1,x2
      LOGICAL, intent(out) :: succes
c-----
      INTERFACE
        FUNCTION func(x)
        use prec_mod, ONLY : prec
        implicit none
        REAL(prec), intent(in) :: x
        REAL(prec) :: func
        END FUNCTION func
      END INTERFACE
c-----
      INTEGER(SELECTED_INT_KIND(9)), PARAMETER :: NTRY=50
      REAL(prec), PARAMETER :: FACTOR=1.6_prec
      INTEGER(SELECTED_INT_KIND(9)) :: j
      REAL(prec) :: f1,f2
      IF (x1 == x2) THEN
        PRINT*, 'rootBrac: must guess an initial range'
        CALL EXIT(1)
      END IF
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      DO j=1,NTRY
        IF ((f1 > 0.0 .AND. f2 < 0.0) .OR.
     x       (f1 < 0.0 .AND. f2 > 0.0)) RETURN
        IF (ABS(f1) < ABS(f2)) THEN
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        ELSE
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        END IF
      END DO
      succes=.false.
      END subroutine rootBrac
c----------------------------------------------------------------------
c-----
c----- use Bisection method
      FUNCTION rootBisect(func,x1,x2,xacc) RESULT (rtbis)
      implicit none
      REAL(prec), intent(in) :: x1,x2,xacc
      REAL(prec) :: rtbis
c-----
      INTERFACE
        FUNCTION func(x)
        use prec_mod, ONLY : prec
        implicit none
        REAL(prec), intent(in) :: x
        REAL(prec) :: func
        END FUNCTION func
      END INTERFACE
c-----
      INTEGER(SELECTED_INT_KIND(9)), PARAMETER :: MAXIT=100
      INTEGER(SELECTED_INT_KIND(9)) :: j
      REAL(prec) :: dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      IF (f*fmid >= 0.0) THEN
        PRINT*, 'Error: rootBisect: root must be bracketed'
        CALL EXIT(1)
      END IF
      IF (f < 0.0) THEN
        rtbis=x1
        dx=x2-x1
      ELSE
        rtbis=x2
        dx=x1-x2
      END IF
      DO j=1,MAXIT
        dx=dx*0.5_prec
        xmid=rtbis+dx
        fmid=func(xmid)
        IF (fmid <= 0.0) rtbis=xmid
        IF (ABS(dx) < xacc .OR. fmid == 0.0) RETURN
      END DO
      PRINT*, 'rootBisect: too many bisections'
      CALL EXIT(1)
      END FUNCTION rootBisect
c-----------------------------------------------------------------------
c-----
c----- use Brents method for line solutions
      FUNCTION rootBrent(func,x1,x2,tol) RESULT (zbrent)
      implicit none
      REAL(prec), intent(in) :: x1,x2,tol
      REAL(prec) :: zbrent
c-----
      INTERFACE
        FUNCTION func(x)
        use prec_mod, ONLY : prec
        implicit none
        REAL(prec), intent(in) :: x
        REAL(prec) :: func
        END FUNCTION func
      END INTERFACE
c-----
      INTEGER(SELECTED_INT_KIND(9)), PARAMETER :: ITMAX=100
      REAL(prec), PARAMETER :: EPS=epsilon(x1)
      INTEGER(SELECTED_INT_KIND(9)) :: iter
      REAL(prec) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      IF ((fa > 0.0 .AND. fb > 0.0) .OR. (fa < 0.0 .AND. fb < 0.0)) THEN
        PRINT*, 'Error: root must be bracketed for Brents method'
        CALL EXIT(1)
      END IF

      c=b
      fc=fb
      DO iter=1,ITMAX
        IF ((fb > 0.0 .AND. fc > 0.0) .OR. 
     x       (fb < 0.0 .AND. fc < 0.0)) THEN
          c=a
          fc=fa
          d=b-a
          e=d
        END IF
        IF (ABS(fc) < ABS(fb)) THEN
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        END IF
        tol1=2.0_prec*EPS*ABS(b)+0.5_prec*tol
        xm=0.5_prec*(c-b)
        IF (ABS(xm) <= tol1 .OR. fb == 0.0) THEN
          zbrent=b
          RETURN
        END IF
        IF (ABS(xm) >= tol1 .AND. ABS(fa) > ABS(fb)) THEN
          s = fb/fa
          IF (a == c) THEN
            p=2.0_prec*xm*s
            q=1.0_prec-s
          ELSE
            q=fa/fc
            r=fb/fc
            p=s*(2.0_prec*xm*q*(q-r)-(b-a)*(r-1.0_prec))
            q=(q-1.0_prec)*(r-1.0_prec)*(s-1.0_prec)
          END IF
          IF (p > 0.0) q=-q
          p=ABS(p)
          IF (2.0_prec*p < min(3.0_prec*xm*q-ABS(tol1*q),ABS(e*q))) THEN
            e=d
            d=p/q
          ELSE
            d=xm
            e=d
          END IF
        ELSE
          d=xm
          e=d
        END IF
        a=b
        fa=fb
        b=b+merge(d,SIGN(tol1,xm), ABS(d) > tol1)
        fb=func(b)
      END DO
      PRINT*, 'rootBrent: exceeded maximum iterations'
      CALL EXIT(1)
      END FUNCTION rootBrent
c-----
c-----------------------------------------------------------------------
c----- bisection function to find nearest point in grid (to left)
c-----------------------------------------------------------------------
      FUNCTION mduc_grid_locate(n1,x1,x) result (loc)
      USE prec_mod, ONLY : prec
c-----
      INTEGER,    INTENT(in) :: n1
      REAL(prec), INTENT(in) :: x
      REAL(prec), INTENT(in), DIMENSION(n1) :: x1
c-----
      INTEGER    :: jL,jR,jM,ascnd,loc
      REAL(prec) :: res
c-----
      jL = 1
      jR = n1
      IF (x1(n1) .GE. x1(1)) THEN
        ascnd = 1
      ELSE
        ascnd = 0
      END IF
c-----
      DO WHILE (jR - jL > 1)
        jM = (jL+jR)/2
        IF (x .GE. x1(jM)) THEN
          IF (ascnd .EQ. 1) THEN
            jL = jM;
          ELSE
            jR = jM;
          END IF
        ELSE
          IF (ascnd .EQ. 1) THEN
            jR = jM;
          ELSE
            jL = jM;
          END IF
        END IF
      END DO
c-----
      loc = MAX(1,MIN(n1-1,jL))
      RETURN
      END FUNCTION mduc_grid_locate
c-----
c-----------------------------------------------------------------------
      END MODULE rootsolve
