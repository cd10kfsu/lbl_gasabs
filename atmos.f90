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
  public :: destroy_atmos

Contains

subroutine destroy_atmos( atms )
  implicit none

  type(t_atmos), intent(inout) :: atms
  
  atms%nlayer = 0
  atms%nlevel = 0
  atms%ngas   = 0
  If (Allocated(atms%pres))  Deallocate( atms%pres )
  If (Allocated(atms%z))     Deallocate( atms%z )
  If (Allocated(atms%temp))  Deallocate( atms%temp ) 
  If (Allocated(atms%gas))   Deallocate( atms%gas )
  If (Allocated(atms%dz))    Deallocate( atms%dz )
  If (Allocated(atms%presl)) Deallocate( atms%presl )
  If (Allocated(atms%templ)) Deallocate( atms%templ )
  If (Allocated(atms%gasl))  Deallocate( atms%gasl )

endsubroutine

endmodule

