MODULE MODH_Spline1D

!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part
! ----------------------------------------------------------------------------------------------------------------------------------
! Public Part
! ----------------------------------------------------------------------------------------------------------------------------------
!INTERFACE GetSplineVandermonde
!   MODULE PROCEDURE GetSplineVandermonde 
!END INTERFACE

INTERFACE Spline
   MODULE PROCEDURE Spline
END INTERFACE

INTERFACE EvalSpline
   MODULE PROCEDURE EvalSpline
END INTERFACE

PUBLIC::GetSplineVandermonde
PUBLIC::Spline
PUBLIC::EvalSpline


!===================================================================================================================================


CONTAINS


SUBROUTINE GetSplineVandermonde(nIn,nOut,SplVdm,xiOut_opt)
!============================================================================================================================
! a regular knot span spline to interpolate nIn points is interpolated
! at a regular reference space point distribution at nOut points
! only the "Lagrange" spline functions are evaluated, and the 
! interpolated values can be found by:
! yOut(1:nOut)=MATMUL(SplVdm(1:nOut,1:nIn),yIn(1:nIn))
!----------------------------------------------------------------------
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nIn !number of interpolation points
INTEGER,INTENT(IN) :: nOut !number of regular points for spline evaluation
REAL,INTENT(IN),OPTIONAL  :: xiOut_opt(nOut) ![-1,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)  :: SplVdm(nOut,nIn)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: i,j   
!REAL,DIMENSION(nIn):: xiIn,b,c,d
REAL,DIMENSION(nIn):: yi,si
REAL               :: xiOut(nOut),step
!============================================================================================================================
!DO j=1,nIn
!  xiIn(j)=REAL(j-1)
!END DO
IF(PRESENT(xiOut_opt))THEN
  !normalize to spline interval [0;Nin-1]
  xiOut=REAL(nIn-1)*0.5*(xiOut_opt+1.)
ELSE
  ! regular point distribution
  step=REAL(nIn-1)/REAL(nOut-1)
  DO i=1,nOut
    xiOut(i)=REAL(i-1)*step  
  END DO
END IF
DO j=1,nIn
  yi=0.
  yi(j)=1.
  !CALL spline (xiIn, yi, b, c, d,nIn) 
  CALL spline_regknot(nIn,yi,(/3,3/),si)
  DO i=1,nOut
    !SplVdm(i,j)= evalSpline(xiOut(i), xiIn, yi, b, c, d, nIn)
    SplVdm(i,j)=  evalSpline_regknot(xiOut(i),nIn,yi,si)
  END DO
END DO
  
END SUBROUTINE GetSplineVandermonde



SUBROUTINE spline (x, y, b, c, d, n)
!============================================================================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!============================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: n
REAL,INTENT(IN)    :: x(n), y(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: b(n), c(n), d(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i, j, gap
REAL    :: h
!=======================================================================
!============================================================================================================================

gap = n-1
! check input
if ( n < 2 ) RETURN
if ( n < 3 ) THEN
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  RETURN
END IF
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
DO i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
END DO
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
IF(n /= 3) THEN
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
END IF
!
! step 3: forward elimination 
!
DO i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
END DO
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
DO j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
END DO
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
DO i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
END DO
c(n) = 3.0*c(n)
d(n) = d(n-1)
END SUBROUTINE spline

FUNCTION evalSpline(u, x, y, b, c, d, n)
!============================================================================================================================
! function evalSpline evaluates the cubic spline interpolation at point z
! evalSpline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! evalSpline = interpolated value at point u
!============================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: u
REAL,INTENT(IN)    :: x(n), y(n)
REAL,INTENT(IN)    :: b(n), c(n), d(n)
INTEGER,INTENT(IN) :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL    :: evalSpline
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i, j, k
REAL    :: dx
!============================================================================================================================

! if u is ouside the x() interval take a boundary value (left or right)
IF(u <= x(1)) THEN
  evalSpline = y(1)
  RETURN
END IF
IF(u >= x(n)) THEN
  evalSpline = y(n)
  RETURN
END IF
!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
DO WHILE (j > i+1)
  k = (i+j)/2
  IF(u < x(k)) THEN
    j=k
  ELSE
    i=k
  END IF
END DO
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
evalSpline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
END FUNCTION evalSpline



SUBROUTINE spline_regknot(nIn,yi,whichBC,si)
!============================================================================================================================
! interpolate a spline trhough nIn data points, with a regular knot sequence 0,...nIn-1
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nIn  !number of interpolation points
INTEGER,INTENT(IN) :: whichBC(2) !boundary conditions left(1) and right(2)
                                 !BC 1: S''=q0/qn given (=0 for natural BC),
                                 !   2: parabolic,
                                 !   3:not-a-knot (free): S'''(1)=S'''(2),
                                 !            S'''(n)=S'''(n-1), continuous 3rd derivativeo
                                 !   4:clamped(t1,t2) !not used
                            !
REAL,INTENT(IN)    :: yi(nIn)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: si(nIn)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL, DIMENSION(nIn) :: a,b,c,r
REAL                 :: q0,qn,sbeta
INTEGER              :: i
!============================================================================================================================
a=0.
b=0.
c=0.

c=1. ! righ diag, since dt=const=1 easy!

a=c !left from diag 
b=4. !diagonal
b(1)=1.
b(nIn)=1.

r=0. !rhs, continuity condition
DO i=2,nIn-1 
  r(i)=3.*(yi(i+1)-yi(i-1))
END DO


!left BC  
SELECT CASE(whichBC(1))
  CASE(1)
    q0=0.
    ! S''=q0, natural BC in q0=0 
    r(1)=1.5*(yi(2)-yi(1))-0.25*q0
    c(1)=0.5
  CASE(2)
    ! parabolic: S'''=0
    r(1)=2.*(yi(2)-yi(1))
    c(1)=1
  CASE(3)
    ! not-a-knot
    r(1)=2.5*(yi(2)-yi(1))+0.5*(yi(3)-yi(2))
    c(1)=2.
  !CASE(4)
  !  !clamped
  !  r(1)=tangvec(1)
  !  c(1)=0.
END SELECT
!right BC
SELECT CASE( whichBC(2))
  CASE(1)
    qn=0.
    ! S''=qn, natural BC if qn=0 
    r(nIn)=1.5*(yi(nIn)-yi(nIn-1))+0.25*qn
    a(nIn)=0.5
  CASE(2)
    ! parabolic: S'''=0
    r(nIn)=2.*(yi(nIn)-yi(nIn-1))
    a(nIn)=1
  CASE(3)
    ! not-a-knot
    r(nIn)=0.5*(yi(nIn-1)-yi(nIn-2))+2.5*(yi(nIn)-yi(nIn-1))
    a(nIn)=2.
  !CASE(4)
  !  !clamped
  !  r(nIn)=tangvec(2)
  !  a(nIn)=0.
END SELECT

!tridiagonal solve (thomas algo) (c and r are changed)
si=0.

c(1)=c(1)/b(1)
r(1)=r(1)/b(1)
DO i=2,nIn
  sbeta=1./(b(i)-c(i-1)*a(i))
  c(i)=c(i)*sbeta
  r(i)=(r(i)-r(i-1)*a(i))*sbeta
END DO
si(nIn)=r(nIn)
DO i=nIn-1,1,-1
  si(i)=r(i)-c(i)*si(i+1)
END DO

END SUBROUTINE spline_regknot


FUNCTION evalSpline_regknot(u, nIn,yi,si)
!============================================================================================================================
! function evalSpline evaluates the cubic spline interpolation at point
! t=u-x(i)
! evalSpline = y(i)(1-t)^3+(y(i)+1/3s(i))*t*(1-t)^2+(y(i+1)-1/3s(i+1))*t^2(1-t)+y(i+1)*t**3
! where  x(i) <= u <= x(i+1), and a regular knot span with x(i)=i-1 
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! yi    = the array of given data points
! si    = derivatives at the given data points
! nIn       = the number of data points
! output:
! evalSpline = interpolated value at point u
!============================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nIn
REAL,INTENT(IN)    :: u,  yi(nIn), si(nIn)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL    :: evalSpline_regknot
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: k
REAL    :: t,t1
!============================================================================================================================

! if u is ouside the x() interval take a boundary value (left or right)
IF(u .LE. 0.) THEN
  evalSpline_regknot = yi(1)
  RETURN
END IF
IF(u .GE. REAL(nIn-1)) THEN
  evalSpline_regknot = yi(nIn)
  RETURN
END IF
 k=MAX(MIN(nIn-1,FLOOR(1.+u)),1)
 t=u-REAL(k-1) !t
 t1=1.-t   
 evalSpline_regknot=t1*t1*(t1*yi(k)+t*(3*yi(k)+si(k)))+t*t*(t1*(3*yi(k+1)-si(k+1))+t*yi(k+1))
END FUNCTION evalSpline_regknot

END MODULE MODH_spline1D
