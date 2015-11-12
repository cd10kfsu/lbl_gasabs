module atmos
!$$$  program documentation block
!         .           .            .
!  program name: atmos
!    programmer: da,cheng        org: umd      date: 2015-Nov-09
!
!  purpose:
!
!  revision history:
!    2015-Nov-09     da    - creator
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block
  use type_define, only : i4, r4, r8, t_atmos
  implicit none

  private
  public :: create_atmos
  public :: destroy_atmos

Contains

subroutine create_atmos( atmos, nlevel, ngas )
  implicit none
  type(t_atmos), intent(inout) :: atmos
  integer(i4),   intent(in) :: nlevel
  integer(i4),   intent(in) :: ngas

  atmos%nlevel = nlevel
  atmos%nlayer = atmos%nlevel - 1
  atmos%ngas   = ngas
  Allocate( atmos%pres(atmos%nlevel), atmos%z(atmos%nlevel), &
            atmos%temp(atmos%nlevel), atmos%gas(atmos%nlevel,atmos%ngas) )
  Allocate( atmos%presl(atmos%nlayer), atmos%dz(atmos%nlayer), &
            atmos%templ(atmos%nlayer), atmos%gasl(atmos%nlayer,atmos%ngas) )
endsubroutine

subroutine destroy_atmos( atmos )
  implicit none

  type(t_atmos), intent(inout) :: atmos
  
  atmos%nlayer = 0
  atmos%nlevel = 0
  atmos%ngas   = 0
  If (Allocated(atmos%pres))  Deallocate( atmos%pres )
  If (Allocated(atmos%z))     Deallocate( atmos%z )
  If (Allocated(atmos%temp))  Deallocate( atmos%temp ) 
  If (Allocated(atmos%gas))   Deallocate( atmos%gas )
  If (Allocated(atmos%dz))    Deallocate( atmos%dz )
  If (Allocated(atmos%presl)) Deallocate( atmos%presl )
  If (Allocated(atmos%templ)) Deallocate( atmos%templ )
  If (Allocated(atmos%gasl))  Deallocate( atmos%gasl )

endsubroutine

endmodule

