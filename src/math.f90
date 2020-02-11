module math
  use precision

  ! Trigonometric parameters
  real(WP), parameter :: Pi    = 3.1415926535897932385_WP
  real(WP), parameter :: twoPi = 6.2831853071795864770_WP
  
  ! Bessel first zero
  real(WP), parameter :: bessj1_zero = 3.8317059702075123115_WP

  ! Blasius data points
  real(WP), dimension(0:9), private :: by0 = (/ &
       0.0_WP, 0.165571818583440_WP, 0.650024518764203_WP, 1.39680822972500_WP, &
       2.30574664618049_WP, 3.28327391871370_WP, 4.27962110517696_WP, &
       5.27923901129384_WP, 6.27921363832835_WP, 7.27921257797747_WP /)
  real(WP), dimension(0:9), private :: by1 = (/ &
       0.0_WP, 0.329780063306651_WP, 0.629765721178679_WP, 0.84604458266019_WP, &
       0.95551827831671_WP, 0.99154183259084_WP, 0.99897290050990_WP, &
        0.9999216098795_WP, 0.99999627301467_WP, 0.99999989265063_WP /)
  real(WP), dimension(0:9), private :: by2 = (/ &
       0.332057384255589_WP, 0.323007152241930_WP, 0.266751564401387_WP, 0.161360240845588_WP, &
       0.06423404047594_WP, 0.01590689966410_WP, 0.00240199722109_WP, &
       0.00022016340923_WP, 0.00001224984692_WP, 0.00000041090325_WP /)

contains
  
  ! Returns the beta distribution
  real(WP) function betaDist(xx, alpha, beta)
    implicit none

    real(WP), intent(in) :: xx, alpha, beta
    real(WP) :: tmp, fgam

    fgam = gammaln(alpha+beta) - gammaln(alpha) - gammaln(beta);
    tmp = (alpha-1.0_WP)*log(xx) + (beta-1.0_WP)*log(1.0_WP-xx) + fgam;
    
    betaDist = exp(tmp)

    return

  end function betaDist

  ! Returns the gamma function
  function gamma(xx)
    implicit none
    
    real(WP) :: gamma
    real(WP), intent(in) :: xx
    
    gamma=exp(gammaln(xx))
    
    return
  end function gamma
  
  
  ! Returns the log of the gamma function
  function gammaln(xx)
    implicit none
  
    real(WP) :: gammaln
    real(WP), intent(in) :: xx
    
    real(WP), parameter :: stp = 2.5066282746310005_WP
    real(WP), dimension(6), parameter :: cof = (/ 76.18009172947146_WP, &
         -86.50532032941677_WP, 24.01409824083091_WP,-1.231739572450155_WP, &
         .1208650973866179E-2_WP, -.5395239384953E-5_WP /)
    
    real(WP) :: ser,tmp,x,y
    integer :: j

    x = xx
    y = x
    tmp = x + 5.5_WP
    tmp = (x+0.5_WP)*log(tmp)-tmp
    ser = 1.000000000190015_WP
    do j=1,6
       y = y + 1.0_WP
       ser = ser+cof(j)/y
    end do
    gammaln = tmp + log(stp*ser/x)
    
    return
  end function gammaln

  ! Error function
  function erf(xx)
    implicit none
    real(WP), intent(in) :: xx
    real(WP) :: erf

    if (xx.ge.0.0_WP) then
       erf = 1.0_WP - erfccheb(xx);
    else
       erf = erfccheb(-xx) - 1.0_WP;
    end if
    return
  end function erf

  ! Complementary error function
  function erfc(xx)
    implicit none
    real(WP), intent(in) :: xx
    real(WP) :: erfc

    if (xx.ge.0.0_WP) then
       erfc = erfccheb(xx);
    else
       erfc = 2.0_WP - erfccheb(-xx);
    end if
    return
  end function erfc

  ! inverse complementary error function
  function inverfc(p)
    implicit none

    real(WP), intent(in) :: p
    real(WP) :: x, err, t, pp
    real(WP) :: inverfc
    integer :: j

    ! arbitrary large numbers for boundaries
    if (p.ge.2.0_WP) then
       inverfc = -100.0_WP
       return
    end if
    if (p.le.0.0_WP) then
       inverfc = 100.0_WP
       return
    end if

    if (p.lt.1.0_WP) then
       pp = p
    else
       pp = 2.0_WP - p
    end if

    t = sqrt(-2.0_WP*log(0.5_WP*pp)) ! initial guess
    x = -0.70711_WP*((2.30753_WP+t*0.27061_WP)/(1.0_WP+t*(0.99229_WP+t*0.04481_WP)) - t)
    do j=1,2
       err = erfc(x) - pp
       x = x + err/(1.12837916709551257_WP*exp(-x**2)-x*err)
    end do

    if (p.lt.1.0_WP) then 
       inverfc = x
    else
       inverfc = -x
    end if
    return
  end function inverfc

  ! inverse error function
  function inverf(p)
    implicit none
    real(WP), intent(in) :: p
    real(WP) :: inverf
    inverf = inverfc(1.0_WP-p)
    return
  end function inverf

  ! approximate error function using Chebyshev coefficients
  function erfccheb(z)
    implicit none
    real(WP) :: erfccheb
    integer :: j, ncof = 28
    real(WP), intent(in) :: z
    real(WP) :: t, ty, tmp, d, dd
    real(WP), dimension(28), parameter :: cof = (/ -1.3026537197817094_WP, &
         6.4196979235649026e-1_WP, 1.9476473204185836e-2_WP, &
         -9.561514786808631e-3_WP, -9.46595344482036e-4_WP, &
         3.66839497852761e-4_WP, 4.2523324806907e-5_WP, -2.0278578112534e-5_WP,&
         -1.624290004647e-6_WP, 1.303655835580e-6_WP, 1.5626441722e-8_WP, &
         -8.5238095915e-8_WP, 6.529054439e-9_WP, 5.059343495e-9_WP, &
         -9.91364156e-10_WP, -2.27365122e-10_WP, 9.6467911e-11_WP, &
         2.394038e-12_WP, -6.886027e-12_WP, 8.94487e-13_WP, 3.13092e-13_WP, &
         -1.12708e-13_WP, 3.81e-16_WP, 7.106e-15_WP, -1.523e-15_WP, &
         -9.4e-17_WP, 1.21e-16_WP, -2.8e-17_WP /)

    if (z.lt.0.0) then
       write (*,*) 'Error: erfccheb requires non-negative argument'
       call exit(1)
    end if
    d  = 0.0_WP
    dd = 0.0_WP
    t  = 2.0_WP/(2.0_WP + z);
    ty = 4.0_WP*t - 2.0_WP;
    do j=ncof,2,-1
       tmp = d;
       d   = ty*d - dd + cof(j);
       dd  = tmp;
    end do
    erfccheb = t*exp(-z*z + 0.5_WP*(cof(1) + ty*d) - dd);
    return
  end function erfccheb

  ! Returns the Bessel function J0(x) for any real x.
  ! [Numerical Recipes in Fortran]
  function bessj0(x)
    implicit none

    real(WP) :: bessj0
    real(WP), intent(in) :: x
    
    real(WP), dimension(6), parameter :: r = (/57568490574.0_WP,-13362590354.0_WP,651619640.7_WP, &
         -11214424.18_WP,77392.33017_WP,-184.9052456_WP/)
    real(WP), dimension(6), parameter :: s = (/57568490411.0_WP,1029532985.0_WP,9494680.718_WP, &
         59272.64853_WP,267.8532712_WP,1.0_WP/)
    real(WP), dimension(5), parameter :: p = (/1.0_WP,-.1098628627E-2_WP,.2734510407E-4_WP, &
         -.2073370639E-5_WP,.2093887211E-6_WP/)
    real(WP), dimension(5), parameter :: q = (/-.1562499995E-1_WP,.1430488765E-3_WP,-.6911147651E-5_WP, &
         .7621095161E-6_WP,-.934945152E-7_WP/)

    real(WP) :: ax,xx,z,y
    
    if(abs(x).lt.8.0_WP)then
       y=x**2
       bessj0=(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) / &
            (s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*s(6))))))
    else
       ax=abs(x)
       z=8.0_WP/ax
       y=z**2
       xx=ax-.785398164_WP
       bessj0=sqrt(.636619772_WP/ax)*(cos(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5))))) - &
            z*sin(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))))
    endif
    return
  end function bessj0

  ! Returns the Bessel function J1(x) for any real x.
  ! [Numerical Recipes in Fortran]
  function bessj1(x)
    implicit none

    real(WP) :: bessj1
    real(WP), intent(in) :: x
    
    real(WP), dimension(6), parameter :: r = (/72362614232.0_WP,-7895059235.0_WP,242396853.1_WP, &
         -2972611.439_WP,15704.48260_WP,-30.16036606_WP/)
    real(WP), dimension(6), parameter :: s = (/144725228442.0_WP,2300535178.0_WP,18583304.74_WP, &
         99447.43394_WP,376.9991397_WP,1.0_WP/)
    real(WP), dimension(5), parameter :: p = (/1.0_WP,.183105E-2_WP,-.3516396496E-4_WP,.2457520174E-5_WP, &
         -.240337019E-6_WP/)
    real(WP), dimension(5), parameter :: q = (/.04687499995_WP,-.2002690873E-3_WP,.8449199096E-5_WP, &
         -.88228987E-6_WP,.105787412E-6_WP/)

    real(WP) :: ax,xx,z,y
    
    if(abs(x).lt.8.0_WP)then
       y=x**2
       bessj1=x*(r(1)+y*(r(2)+y*(r(3)+y*(r(4)+y*(r(5)+y*r(6)))))) / &
            (s(1)+y*(s(2)+y*(s(3)+y*(s(4)+y*(s(5)+y*s(6))))))
    else
       ax=abs(x)
       z=8.0_WP/ax
       y=z**2
       xx=ax-2.356194491_WP
       bessj1=sqrt(.636619772_WP/ax)*(cos(xx)*(p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*p(5))))) - &
            z*sin(xx)*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*q(5))))))*sign(1.0_WP,x)
    endif
    return
  end function bessj1
  

  ! Return the Blasius function
  function blasius0(x)
    implicit none

    real(WP) :: blasius0
    real(WP), intent(in) :: x
    integer  :: i
    real(WP) :: f2l,f2r,f0l,f0r

    if (x.le.0.0_WP) then
       blasius0 = 0.0_WP
    else if (x.ge.9.0_WP) then
       blasius0 = by0(9) + (x-9.0_WP)
    else
       i = int(x)
       f0l = by0(i)
       f0r = by0(i+1)
       f2l = by2(i)
       f2r = by2(i+1)
       blasius0 = &
            1.0_WP/6.0_WP*f2l*(real(i+1,WP)-x)**3 + &
            1.0_WP/6.0_WP*f2r*(x-real(i,WP))**3 + &
            (f0l-1.0_WP/6.0_WP*f2l)*(real(i+1,WP)-x) + &
            (f0r-1.0_WP/6.0_WP*f2r)*(x-real(i,WP))
    end if
    
    return
  end function blasius0
  
  ! Return the first derivative of the Blasius function
  function blasius1(x)
    implicit none

    real(WP) :: blasius1
    real(WP), intent(in) :: x
    integer  :: i
    real(WP) :: f2l,f2r,f0l,f0r

    if (x.le.0.0_WP) then
       blasius1 = 0.0_WP
    else if (x.ge.9.0_WP) then
       blasius1 = 1.0_WP
    else
       i = int(x)
       f0l = by1(i)
       f0r = by1(i+1)
       f2l = -0.5_WP*by0(i)*by2(i)
       f2r = -0.5_WP*by0(i+1)*by2(i+1)
       blasius1 = &
            1.0_WP/6.0_WP*f2l*(real(i+1,WP)-x)**3 + &
            1.0_WP/6.0_WP*f2r*(x-real(i,WP))**3 + &
            (f0l-1.0_WP/6.0_WP*f2l)*(real(i+1,WP)-x) + &
            (f0r-1.0_WP/6.0_WP*f2r)*(x-real(i,WP))
    end if
    
    return
  end function blasius1
  
  ! Return the second derivative of the Blasius function
  function blasius2(x)
    implicit none

    real(WP) :: blasius2
    real(WP), intent(in) :: x
    integer  :: i
    real(WP) :: f2l,f2r,f0l,f0r

    if (x.le.0.0_WP) then
       blasius2 = by2(0)
    else if (x.ge.9.0_WP) then
       blasius2 = 0.0_WP
    else
       i = int(x)
       f0l = by2(i)
       f0r = by2(i+1)
       f2l = 0.25_WP*by0(i)**2*by2(i) - 0.5_WP*by1(i)*by2(i)
       f2r = 0.25_WP*by0(i+1)**2*by2(i+1) - 0.5_WP*by1(i+1)*by2(i+1)
       blasius2 = &
            1.0_WP/6.0_WP*f2l*(real(i+1,WP)-x)**3 + &
            1.0_WP/6.0_WP*f2r*(x-real(i,WP))**3 + &
            (f0l-1.0_WP/6.0_WP*f2l)*(real(i+1,WP)-x) + &
            (f0r-1.0_WP/6.0_WP*f2r)*(x-real(i,WP))
    end if
    
    return
  end function blasius2
  
  ! Return the third derivative of the Blasius function
  function blasius3(x)
    implicit none

    real(WP) :: blasius3
    real(WP), intent(in) :: x

    blasius3 = -0.5_WP*blasius0(x)*blasius2(x)
    
    return
  end function blasius3
  
  function arctan(dx,dy)
    implicit none
    
    real(WP), intent(in) :: dx,dy
    real(WP) :: arctan
    
    ! Evaluate atan
    if (abs(dx)+abs(dy).lt.1.0e-9_WP) then
       arctan = 0.0_WP
    else
       arctan = atan(dy/dx)
    end if
    
    ! Quadrant correction
    if (dx.le.0.0_WP) then
       arctan = Pi+arctan
    elseif (dy.le.0.0_WP .and. dx.gt.0.0_WP) then
       arctan = twoPi+arctan
    end if
    
    return
  end function arctan
  
  ! Inverse the linear system : AX=B
  ! c added for usage with IBM compiler
  ! otherwise matrices do overlap
  subroutine solve_linear_system(A,B,X,n)
    implicit none

    integer, intent(in) :: n
    real(WP) :: c
    real(WP), dimension(n,n) :: A
    real(WP), dimension(n) :: B
    real(WP), dimension(n) :: X
    integer :: i,j,l
    
    ! Forward elimination
    do i=1,n
       ! Pivoting
       j = i-1+maxloc(abs(A(i:n,i)),1)
       X = A(i,:)
       A(i,:) = A(j,:)
       A(j,:) = X
       X(1) = B(i)
       B(i) = B(j)
       B(j) = X(1)
       ! Elimination
       c = A(i,i)
       B(i) = B(i) / c
       A(i,:) = A(i,:) / c
       do l=i+1,n
          c = A(l,i)
          B(l) = B(l) - c*B(i)
          A(l,:) = A(l,:) - c*A(i,:)
       end do
    end do

    ! Backward substitution
    do i=n,1,-1
       X(i) = B(i)
       do l=i+1,n
          X(i) = X(i) - A(i,l)*X(l)
       end do
    end do
    
    return
  end subroutine solve_linear_system
  
  ! Inverse the linear system : AX=B
  ! c added for usage with IBM compiler
  ! otherwise matrices do overlap
  subroutine solve_linear_system_nopivot(A,B,X,n)
    implicit none

    integer, intent(in) :: n
    real(WP) :: c
    real(WP), dimension(n,n) :: A
    real(WP), dimension(n) :: B
    real(WP), dimension(n) :: X
    integer :: i,l

    ! Forward elimination
    do i=1,n
       c = A(i,i)
       B(i) = B(i) / c
       A(i,:) = A(i,:) / c
       do l=i+1,n
          c = A(l,i)
          B(l) = B(l) - c*B(i)
          A(l,:) = A(l,:) - c*A(i,:)
       end do
    end do

    ! Backward substitution
    do i=n,1,-1
       X(i) = B(i)
       do l=i+1,n
          X(i) = X(i) - A(i,l)*X(l)
       end do
    end do
    
    return
  end subroutine solve_linear_system_nopivot

  
  ! Inverse the matrix
  ! c added for usage with IBM compiler
  ! otherwise matrices do overlap
  subroutine inverse_matrix(A,B,n)
    implicit none

    integer, intent(in) :: n
    real(WP) :: c
    real(WP), dimension(n,n) :: A,B
    integer :: i,l

    B=0.0_WP

    ! Forward elimination
    do i=1,n
       B(i,i) = 1.0_WP
       c = A(i,i)
       B(i,:) = B(i,:) / c
       A(i,:) = A(i,:) / c
       do l=i+1,n
          c = A(l,i)
          B(l,:) = B(l,:) - c*B(i,:)
          A(l,:) = A(l,:) - c*A(i,:)
       end do
    end do

    ! Backward substitution
    do i=n,1,-1
       do l=i+1,n
          B(i,:) = B(i,:) - A(i,l)*B(l,:)
       end do
    end do
    
    return
  end subroutine inverse_matrix
  
  ! Third order finite volume interpolation
  subroutine poly3_coeff(xx,xp,coeff)
    implicit none
    
    real(WP), intent(in),  dimension(4) :: xx
    real(WP), intent(out), dimension(3) :: coeff
    real(WP), intent(in)     :: xp
    real(WP), dimension(3,3) :: A,B
    
    ! Form the matrix
    A(1,1) = (xx(2)**3-xx(1)**3)/(3.0_WP*(xx(2)-xx(1)))
    A(2,1) = (xx(3)**3-xx(2)**3)/(3.0_WP*(xx(3)-xx(2)))
    A(3,1) = (xx(4)**3-xx(3)**3)/(3.0_WP*(xx(4)-xx(3)))
    
    A(1,2) = (xx(2)**2-xx(1)**2)/(2.0_WP*(xx(2)-xx(1)))
    A(2,2) = (xx(3)**2-xx(2)**2)/(2.0_WP*(xx(3)-xx(2)))
    A(3,2) = (xx(4)**2-xx(3)**2)/(2.0_WP*(xx(4)-xx(3)))
    
    A(1,3) = (xx(2)**1-xx(1)**1)/(1.0_WP*(xx(2)-xx(1)))
    A(2,3) = (xx(3)**1-xx(2)**1)/(1.0_WP*(xx(3)-xx(2)))
    A(3,3) = (xx(4)**1-xx(3)**1)/(1.0_WP*(xx(4)-xx(3)))
    
    ! Invert it
    call inverse_matrix(A,B,3)
    
    ! Compute metrics
    coeff(1) = B(1,1)*xp**2+B(2,1)*xp**1+B(3,1)*xp**0
    coeff(2) = B(1,2)*xp**2+B(2,2)*xp**1+B(3,2)*xp**0
    coeff(3) = B(1,3)*xp**2+B(2,3)*xp**1+B(3,3)*xp**0
    
    return
  end subroutine poly3_coeff
  
  ! HO finite volume interpolation
  subroutine hofvi(n,xx,xp,coeff)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n+1) :: xx
    real(WP), intent(out), dimension(n)   :: coeff
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i+1)**(n+1-j)-xx(i)**(n+1-j))/(real(n+1-j,WP)*(xx(i+1)-xx(i)))
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=1,n
          coeff(i) = coeff(i) + B(j,i)*xp**(n-j)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofvi
  
  ! HO finite difference interpolation
  subroutine hofdi(n,xx,xp,coeff)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xp)**(j-1)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = B(1,:)
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdi
  
  ! HO finite difference interpolation
  ! with Dirichlet condition
  subroutine hofdi_d(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=1,n
          coeff(i) = coeff(i) + B(j,i)*(xp-xbc)**(j)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdi_d
  
  ! HO finite difference interpolation
  ! with Neumann condition
  subroutine hofdi_n(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       A(i,1) = 1.0_WP
       do j=2,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       coeff(i) = B(1,i)
       do j=2,n
          coeff(i) = coeff(i) + B(j,i)*(xp-xbc)**(j)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdi_n
  
  ! HO finite difference differenciation
  subroutine hofdd(n,xx,xp,coeff)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xp)**(j-1)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = B(2,:)
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdd
  
  ! HO finite difference differenciation
  ! with Dirichlet condition
  subroutine hofdd_d(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       do j=1,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=1,n
          coeff(i) = coeff(i) + real(j,WP)*B(j,i)*(xp-xbc)**(j-1)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdd_d
  
  ! HO finite difference differenciation
  ! with Neumann condition
  subroutine hofdd_n(n,xx,xp,coeff,xbc)
    implicit none
    
    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    
    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B
    
    integer :: i,j
    
    ! Allocate the matrices
    allocate(A(n,n))
    allocate(B(n,n))
    
    ! Form the matrix
    do i=1,n
       A(i,1) = 1.0_WP
       do j=2,n
          A(i,j) = (xx(i)-xbc)**(j)
       end do
    end do
    
    ! Invert it
    call inverse_matrix(A,B,n)
    
    ! Compute metrics
    coeff = 0.0_WP
    do i=1,n
       do j=2,n
          coeff(i) = coeff(i) + real(j,WP)*B(j,i)*(xp-xbc)**(j-1)
       end do
    end do
    
    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)
    
    return
  end subroutine hofdd_n

  ! HO finite difference differenciation
  ! with Neumann condition of a non-zero value of the gradient
  subroutine hofdd_n_bc(n,xx,xp,coeff,xbc,coeff_bc)
    implicit none

    integer,  intent(in) :: n
    real(WP), intent(in) :: xp
    real(WP), intent(in),  dimension(n) :: xx
    real(WP), intent(out), dimension(n) :: coeff
    real(WP), intent(in) :: xbc
    real(WP), intent(out) :: coeff_bc

    real(WP), dimension(:,:), pointer :: A
    real(WP), dimension(:,:), pointer :: B

    integer :: i,j

    ! Allocate the matrices
    allocate(A(n+1,n+1))
    allocate(B(n+1,n+1))

    ! Form the matrix
    do i=1,n
       A(i,1) = 1.0_WP
       A(i,2) = xx(i)
       do j=2,n
          A(i,j+1) = (xx(i)-xbc)**(j)
       end do
    end do
    A(n+1,:) = 0.0_WP
    A(n+1,2) = 1.0_WP

    ! Invert it
    call inverse_matrix(A,B,n+1)

    ! Compute metrics
    do i=1,n
       coeff(i) = B(2,i)
       do j=2,n
          coeff(i) = coeff(i) + real(j,WP)*B(j+1,i)*(xp-xbc)**(j-1)
       end do
    end do

    coeff_bc = B(2,n+1)
    do j=2,n
       coeff_bc = coeff_bc + real(j,WP)*B(j+1,n+1)*(xp-xbc)**(j-1)
    end do

    ! Deallocate the matrices
    deallocate(A); nullify(A)
    deallocate(B); nullify(B)

    return
  end subroutine hofdd_n_bc
  
  ! Eigenvalue/vector solver for a real symmetric NxN matrix A
  ! A is destroyed, d contains the eigenvalues, v contains the eigenvectors in columns
  subroutine eigensolver(A,n,d,v,nrot)
    use precision
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: nrot
    real(WP), dimension(:), intent(out) :: d
    real(WP), dimension(:,:), intent(inout) :: A
    real(WP), dimension(:,:), intent(out) :: v
    integer :: i,ip,iq
    real(WP) :: c,g,h,s,sm,t,tau,theta,tresh
    real(WP), dimension(size(d)) :: b,z
    
    ! Initialize v to be the identity matrix
    v = 0.0_WP
    do ip=1,n
       do iq=1,n
          v(ip,iq) = 0.0_WP
       end do
       v(ip,ip) = 1.0_WP
    end do
    
    ! Initialize b and d to the diagonal of A and z to 0
    do ip=1,n
       b(ip) = A(ip,ip)
       d(ip) = b(ip)
       z(ip) = 0.0_WP
    end do
    
    ! Eigenvalue solver
    nrot = 0
    do i=1,50
       sm=0.0_WP
       do ip=1,n-1
          do iq=ip+1,n
             sm=sm+abs(A(ip,iq))
          end do
       end do
       if (sm.eq.0.0_WP) return
       if (i.lt.4) then
          tresh=0.2_WP*sm/n**2
       else
          tresh=0.0_WP
       end if
       do ip=1,n-1
          do iq=ip+1,n
             g=100.0_WP*abs(A(ip,iq))
             if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
                A(ip,iq)=0.0_WP
             else if (abs(A(ip,iq)).gt.tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g.eq.abs(h)) then
                   t=A(ip,iq)/h
                else
                   theta=0.5_WP*h/A(ip,iq)
                   t=1.0_WP/(abs(theta)+sqrt(1.0_WP+theta**2))
                   if (theta.lt.0.0_WP) t=-t
                end if
                c=1.0_WP/sqrt(1.0_WP+t**2)
                s=t*c
                tau=s/(1.0_WP+c)
                h=t*A(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                A(ip,iq)=0.0_WP
                call jrotate(A(1:ip-1,ip),   A(1:ip-1,iq)   )
                call jrotate(A(ip,ip+1:iq-1),A(ip+1:iq-1,iq))
                call jrotate(A(ip,iq+1:n),   A(iq,iq+1:n)   )
                call jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             end if
          end do
       end do
       b(:)=b(:)+z(:)
       d(:)=b(:)
       z(:)=0.0_WP
    end do
    print*,'Too many iterations in eigensolver'
    
  contains
    
    subroutine jrotate(a1,a2)
      implicit none
      
      real(WP), dimension(:), intent(inout) :: a1,a2
      real(WP), dimension(size(a1)) :: wk1
      
      wk1(:)=a1(:)
      a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
      a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
      
      return
    end subroutine jrotate
    
  end subroutine eigensolver
  
  ! Eigenvalue solver for tri-diag matrices
  ! d contains the diagonal on input, the eigenvalues on output
  ! e contains the lower diag, e(1) is not used. e is destroyed
  subroutine tqli(d,n,e,z)
    use precision
    implicit none
    
    integer, intent(in) :: n
    real(WP), dimension(n), intent(inout) :: d,e
    real(WP), dimension(n,n), intent(inout) :: z
    integer :: i,iter,l,m
    real(WP) :: b,c,dd,f,g,p,r,s
    real(WP), dimension(n) :: ff
    
    ! Renumber e
    do i=2,n
       e(i-1)=e(i)
    end do
    
    ! Do it
    do l=1,n
       iter=0
       iterate: do
          do m=l,n-1
             dd=abs(d(m))+abs(d(m+1))
             if (abs(e(m))+dd.eq.dd) exit
          end do
          if (m.eq.l) exit iterate
          if (iter.eq.30) stop 'tqli: too many iterations.'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.0_WP*e(l))
          r=pythag(g,1.0_WP)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0_WP
          c=1.0_WP
          p=0.0_WP
          do i=m-1,l,-1
             f=s*e(i)
             b=c*e(i)
             r=pythag(f,g)
             e(i+1)=r
             if (r.eq.0.0_WP) then
                d(i+1)=d(i+1)-p
                e(m)=0.0_WP
                cycle iterate
             end if
             s=f/r
             c=g/r
             g=d(i+1)-p
             r=(d(i)-g)*s+2.0_WP*c*b
             p=s*r
             d(i+1)=g+p
             g=c*r-b
             ! Eigenvector
             ff(1:n)=z(1:n,i+1)
             z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
             z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
          end do
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0_WP
       end do iterate
    end do
    
  contains
    
    function pythag(a,b)
      implicit none
      real(WP) :: a,b
      real(WP) :: pythag
      real(WP) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if (absa.gt.absb) then
         pythag=absa*sqrt(1.0_WP+(absb/absa)**2)
      else
         if (absb.eq.0.0_WP) then
            pythag=0.0_WP
         else
            pythag=absb*sqrt(1.0_WP+(absa/absb)**2)
         end if
      end if
      return
    end function pythag
    
  end subroutine tqli

  ! Adapted from MDUC
  subroutine betaPDF(mean,var,nzz,zz,bpdf,pdfBound)
    implicit none

    ! ---------------------------------
    real(WP), intent(inout) :: mean, var
    real(WP), dimension(:), intent(out) :: zz, bpdf
    integer, intent(inout) :: nzz
    integer, dimension(:), intent(inout) :: pdfBound
    ! ---------------------------------
    real(WP) :: tmp, alpha, beta, fgam, sumPDF, sumPDFI, truncate = 1.0E-06_WP
    real(WP) :: dzz, dzz0, dzzN, dzzI
    integer :: i, j, j0, j1, jMean, pdf_flag

    bpdf = 0.0_WP

    dzz0 = 0.5*(zz(2)-zz(1));     ! for boundary approximation
    dzzN = 0.5*(zz(nzz)-zz(nzz-1));

    ! Ensure that the mean and variance are properly bounded
    call clipMoments(mean,var)

    pdf_flag = 0
    ! Check for special cases of pdf, e.g. delta function
    call checkPDF(nzz, zz, mean, var, bpdf, pdfbound, pdf_flag)

    ! If pdfFlag is greater than zero, special case was found, so exit
    if (pdf_flag > 0) then
       zz(1) = 1.0E-10_WP
       return
    end if

    ! Compute the shape parameters for the beta distribution
    tmp = mean*(1.0_WP-mean)/var - 1.0_WP
    alpha = mean*tmp
    beta  = (1.0_WP-mean)*tmp
    if (abs(alpha - 1.0_WP) < 1.0E-09_WP) alpha = 1.0_WP
    if (abs(beta - 1.0_WP) < 1.0E-09_WP) beta = 1.0_WP

    ! ----- compute a beta distribution based on alpha and beta ----- 
    sumPDF  = 0.0_WP
    if (alpha == 1.0_WP .and. beta == 1.0_WP) then ! uniform distribution (not likely)
       j0 = 1
       j1 = nzz
       dzzI = 1.0_WP / (zz(j1)-zz(j0))
       bpdf(j0) = 0.5_WP*(zz(2)-zz(1))*dzzI
       bpdf(j1) = 0.5_WP*(zz(j1+1)-zz(j1))*dzzI
       sumPDF = bpdf(j0) + bpdf(j1)
       do j=j0+1,j1-1

          dzz = 0.5_WP*(zz(j+1)-zz(j-1))
          bpdf(j) = dzz*dzzI
          sumPDF = sumPDF + bpdf(j)

       end do

       ! U shaped distribution
    else if ( (alpha < 1.0_WP .and. beta < 1.0_WP) ) then

       ! Need to calculate over whole domain, ignore the ends
       j0 = 1
       j1 = nzz

       ! Compute the PDF over all the interior points
       do j=j0+1,j1-1

          dzz = 0.5_WP*(zz(j+1)-zz(j-1));
          bpdf(j) = betaDist(zz(j), alpha, beta)*dzz;
          sumPDF = sumPDF + bpdf(j)

       end do

       ! Decreasing distribution
    else if ( (alpha < 1.0_WP .and. beta >= 1.0_WP) &
         .or.(  alpha == 1.0_WP .and. beta > 1.0_WP) ) then

       ! Start at the left boundary and march till less than truncate
       j0 = 1
       do j = j0+1, nzz-1

          dzz = 0.5_WP*(zz(j+1)-zz(j-1))
          bpdf(j) = betaDist(zz(j), alpha, beta)*dzz

          if (bpdf(j) < truncate) then
             bpdf(j) = 0.0_WP
             exit
          end if

          j1 = j
          sumPDF = sumPDF + bpdf(j)

       end do

       ! Increasing distribution
    else if ( (alpha > 1.0_WP .and. beta <= 1.0_WP) &
         .or. (  alpha == 1.0_WP .and. beta < 1.0_WP) ) then

       ! Start at the right boundary and march till less than truncate
       j1 = nzz
       j = j1-1
       do while (j>1)
          dzz = 0.5_WP*(zz(j+1)-zz(j-1))
          bpdf(j) = betaDist(zz(j), alpha, beta)*dzz

          if (bpdf(j) < truncate) then

             bpdf(j) = 0.0_WP
             exit
          end if

          j0 = j
          sumPDF = sumPDF + bpdf(j)
          j = j -1
       end do

       ! Unimodal distribution
    else

       j0 = 1
       j1 = nzz

       ! Find the location of the mesh point before the mean
       jMean = 1
       do while (zz(jMean) < mean) 
          jMean = jMean + 1
       end do
       !jMean = jMean -1

       ! Start from the mean and decrease zz till the pdf is negligible
       j = jMean
       do while(j>1)

          dzz = 0.5_WP*(zz(j+1)-zz(j-1))
          bpdf(j) = betaDist(zz(j), alpha, beta)*dzz

          if (bpdf(j) < truncate) then

             bpdf(j) = 0.0
             exit
          end if

          j0 = j
          sumPDF = sumPDF + bpdf(j)
          j = j - 1

       end do

       ! Start from the mean and increase zz till the pdf is negligible
       j = jMean+1
       do j=jMean+1,nzz-1

          dzz = 0.5_WP*(zz(j+1)-zz(j-1))
          bpdf(j) = betaDist(zz(j), alpha, beta)*dzz

          if (bpdf(j) < truncate) then

             bpdf(j) = 0.0_WP
             exit
          end if

          j1 = j
          sumPDF = sumPDF + bpdf(j)
       end do

    end if ! End different cases

    ! Compute the boundary points explicitly if necessary using approximate expansion 
    if (.not.(alpha == 1.0_WP .and. beta == 1.0_WP)) then

       fgam = gammaln(alpha+beta) - gammaln(alpha) - gammaln(beta)
       if (j0 == 1 .or. j0 == 2) then
          j0 = 1
          bpdf(j0) = exp(alpha*log(dzz0) + fgam)/alpha
          sumPDF = sumPDF + bpdf(j0)
       end if

       if (j1 == nzz .or. j1 == nzz-1) then
          j1 = nzz
          bpdf(j1) = exp(beta*log(dzzN) + fgam)/beta
          sumPDF = sumPDF + bpdf(j1)
       end if
    end if

    sumPDFI = 1.0_WP / sumPDF
    ! Normalize the PDF
    bpdf = bpdf*sumPDFI


    ! Store the limits of the non-zero pdf
    pdfBound(1) = j0
    pdfBound(2) = j1

    zz(1) = 1.0E-10_WP

  contains

    ! Subroutine based on MDUC function clipMoments1D
    subroutine clipMoments(mean, var)
      implicit none
      real(WP), intent(inout) :: mean, var
      ! ---------------------------------

      mean = min(1.0_WP, max(0.0_WP, mean))
      var  = max(0.0_WP, min(var, mean*(1.0_WP-mean)))

    end subroutine clipMoments

    ! Subroutine based on MDUC function checkPDF_1D
    subroutine checkPDF(n, x, mean, var, pdf, pdfbound, pdf_flag)
      implicit none
      integer, intent(inout) :: n, pdf_flag
      real(WP), intent(inout) :: mean, var
      real(WP), dimension(:), intent(inout) :: x, pdf
      integer, dimension(:), intent(inout) :: pdfbound
      ! ---------------------------------

      integer :: j0, j1, jjMean
      real(WP) :: dx, sigma, sumPDF

      ! delta at x = 0
      if (mean <= 1.0E-10_WP) then

         pdf(1)      = 1.0_WP
         pdfBound(1) = 1
         pdfBound(2) = 1
         pdf_flag    = 1

      end if

      ! delta at x = 1
      if (mean >= (1.0_WP - 1.0E-10_WP) ) then

         pdf(n)      = 1.0_WP
         pdfBound(1) = n
         pdfBound(2) = n
         pdf_flag    = 1

      end if

      ! double delta function
      if (var >= mean*(1.0_WP-mean)) then

         pdf(1)      = 1.0_WP - mean
         pdf(n)      = mean
         pdfBound(1) = 1
         pdfBound(2) = n
         pdf_flag    = 1

      end if

      ! find the location of the mesh point before the mean
      jjMean = 1
      do while (x(jjMean) < mean)
         jjMean = jjMean + 1
      end do
      !jjMean = jjMean - 1

      sigma = sqrt(var)    ! compute the standard deviation

      !if the mesh size about the mean is less than 1.5 standard deviations
      !compute at two points about the mean such that the mean is conserved
      if (x(jjMean+1) - x(jjMean) > 1.5_WP*sigma) then

         j0 = jjMean
         j1 = jjMean + 1

         dx = (x(j1)-x(j0))
         bpdf(j0) = abs(x(j1) - mean)/dx
         bpdf(j1) = abs(mean - x(j0))/dx

         ! normalise the PDF
         sumPDF  = bpdf(j0)+bpdf(j1)
         bpdf(j0) = bpdf(j0)/sumPDF
         bpdf(j1) = bpdf(j1)/sumPDF

         ! store bounds of non-zero PDF values
         pdfBound(1) = j0
         pdfBound(2) = j1
         pdf_flag    = 1

      end if

    end subroutine checkPDF

  end subroutine betaPDF
  
  subroutine RosinRammlerPDF(xm,x2,xmax,nx,x,f)
    implicit none

    ! ---------------------------------
    real(WP), intent(in) :: xm, x2, xmax
    integer, intent(in) :: nx
    real(WP), dimension(:), intent(out) :: x, f
    ! ---------------------------------
    
    real(WP), dimension(nx+1) :: x_
    real(WP), dimension(nx) :: dx_, dx, dx_adapt
    real(WP), dimension(nx-1) :: invdiff
    real(WP) :: dxx, xbar, r_min, r_max, q, ratio, norm
    integer :: i
    
    xbar = 0.0_WP

    dxx = 2.0_WP*xmax/(2.0_WP*nx-1.0_WP)

    x_ = (/ (dxx*real(i,WP),i=0,nx,1) /)
    f = (/ (real(0.0,WP),i=1,nx,1) /)
    x = (x_(:nx)+0.5_WP*dxx)

    r_min = 1.05_WP; r_max = 1.9_WP
    ratio = x2/xm**2

    if (xm > 0.0_WP .and. x2 > 0.0_WP) then

       if (ratio > r_min ) then
          ratio = min(r_max,ratio)
          call bisection(q,ratio)
       else
          q = 5.1334_WP
       end if

       xbar = xm/gamma(1.0_WP+1.0_WP/q)

       f = (q/xbar**q)*x**(q-1.0_WP)*exp(-(x/xbar)**q)

       norm = sum(f)
       if (norm > 0.0_WP) then
          f = f/norm;
       end if

       ! Adapt 
       invdiff = 1.0_WP/abs(f(:nx-1)-f(2:nx))

       dx_adapt(2:nx) = invdiff/maxval(invdiff)
       dx_adapt(1) = dx_adapt(2)
       
       dx_ = dx_adapt/sum(dx_adapt)
       
       x_ = 0.0_WP
       do i=2,nx+1
          x_(i) = x_(i-1) + dx_(i-1) 
       end do
      
       x_ = xmax*x_

       dx = x_(2:nx+1)-x_(1:nx)

       x = x_(:nx)+0.5_WP*dx

       f = (q/xbar**q)*x**(q-1.0_WP)*exp(-(x/xbar)**q)

       norm = sum(f*dx)
       if (norm > 0.0_WP) then
          f = f*dx/norm;
       end if

    end if
  contains

    subroutine bisection(q_final,ratio)
      implicit none

      real(WP), intent(inout) :: q_final
      real(WP), intent(in) :: ratio
      ! -------------------------------
      real(WP), parameter :: ql0 = 1.0, qh0 = 10.0_WP, tol = 1E-06
      real(WP) :: fs, fl, fh, qs, ql, qh, err
      integer :: iter
      ! -------------------------------

      iter = 0
      err = 1.0_WP
      ql = ql0; qh = qh0; qs = 0.5_WP*(ql+qh)
      do while ( err > tol)
         iter = iter + 1
         fs = gamma(1.0_WP+2.0_WP/qs)/(gamma(1.0_WP+1.0_WP/qs))**2 - ratio
         fl = gamma(1.0_WP+2.0_WP/ql)/(gamma(1.0_WP+1.0_WP/ql))**2 - ratio
         fh = gamma(1.0_WP+2.0_WP/qh)/(gamma(1.0_WP+1.0_WP/qh))**2 - ratio

         if (fl*fs < 0.0_WP) then
            qh = qs; qs = 0.5_WP*(ql+qh)
         elseif (fs*fh < 0.0_WP) then
            ql = qs; qs = 0.5_WP*(ql+qh)
         end if

         err = abs(fs)

      end do
      !write(*,*) 'Bisection iters:', iter
      q_final = qs
    end subroutine bisection

  end subroutine RosinRammlerPDF

  logical function isnan(a)
    real(WP) :: a
    if (a.ne.a) then
       isnan = .true.
    else
       isnan = .false.
    end if
    return
  end function isnan

  logical function isinf(a)
    real(WP) ::  a

    if ((a-1.0_WP).eq.a) then
       isinf = .true.
    else
       isinf = .false.
    end if
    return
  end function isinf

end module math
