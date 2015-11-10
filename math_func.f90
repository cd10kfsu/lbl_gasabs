module math_func
!$$$  program documentation block
!         .           .            .
!  program name: math_func
!    programmer: da,cheng        org: umd      date: 2015-Nov-08
!
!  purpose:
!
!  revision history:
!    2015-Nov-08     da    - creator
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block
  use type_define, only : i4, r4, r8
  implicit none

  private
  public :: linspace
  public :: linterp1d
  public :: trpint

Contains
subroutine linspace(v1,v2,n,varr)
! purpose:
!  Generate n value in the interval [v1, v2] with equal interval
!
  implicit none
  integer(i4),intent(in) :: n
  real(r8),intent(in)    :: v1, v2
  real(r8),allocatable   :: varr(:)
  !local 
  real(r8)    :: dv
  integer(i4) :: i

  allocate( varr(n) )
  dv = (v2-v1)*1./(n-1)
  varr(1) = v1; varr(n) = v2
  do i = 2, n-1 
     varr(i) = varr(i-1) + dv
  enddo
endsubroutine


function linterp1d( n, x, y, xint ) result( yint )
!$$$ subprogram documentation block
!             .                   .                    . 
!  Program name:
!    programmer: Da,Cheng       org: umd            date: 2015-MON-DD
!
!  Purpose:
!    
!
!  Revision history:
!    2015-MON-DD     Da,Cheng     -Creator
!
!  Input args:
!
!  Output args:
!
!$$$
  implicit none
! Passed args
  integer(i4),intent(in) :: n
  real(r8),   intent(in) :: x(n)
  real(r8),   intent(in) :: y(n)
  real(r8),   intent(in) :: xint
  real(r8)               :: yint
! Local vars
  integer(i4) :: i

  If ( xint < x(1) ) Then
     yint = y(1)
  Elseif ( xint >= x(n) ) Then
     yint = y(n)
  Else
    Do i = 1, n-1
       If ( xint >= x(i) .and. xint < x(i+1) ) Exit
    Enddo
    If ( ABS(x(i+1)-xint) <= ABS(xint-x(i)) ) Then
       yint = (y(i+1)-y(i))*(xint-x(i))/(x(i+1)-x(i)) + y(i)
    Else
       yint = (y(i+1)-y(i))*(xint-x(i+1))/(x(i+1)-x(i)) + y(i+1)
    Endif
    !write(6,*) "x0,xint,x1=", x(i), xint, x(i+1)
    !write(6,*) "y0,yint,y1=", y(i), yint, y(i+1)
  Endif

endfunction


function trpint( n, x, y ) result (integral)
  implicit none

  integer(i4),intent(in) :: n
  real(r8),   intent(in) :: x(n)
  real(r8),   intent(in) :: y(n)
  real(r8)               :: integral

  integer(i4) :: i

  integral = 0._r8
  Do i = 1, n-1
     integral = integral + (y(i+1)+y(i))*(x(i+1)-x(i))/2
  Enddo
  
endfunction

endmodule

