module type_define
!$$$  program documentation block
!         .           .            .
!  program name: type_define
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
  implicit none

  public 

  integer,parameter :: i4 = 4
  integer,parameter :: i8 = 8
  
  integer,parameter :: r4 = 4
  integer,parameter :: r8 = 8


  ! TOA                    SFC
  !  o----------o-----------o
  ! lev1       lev2        lev3
  !      lay1        lay2
  type t_atmos
    integer(i4) :: nlayer
    integer(i4) :: nlevel
    integer(i4) :: ngas
    real(r8),allocatable :: pres(:)
    real(r8),allocatable :: z(:)
    real(r8),allocatable :: temp(:)
    real(r8),allocatable :: gas(:,:)

    real(r8),allocatable :: dz(:)
    real(r8),allocatable :: presl(:)
    real(r8),allocatable :: templ(:)
    real(r8),allocatable :: gasl(:,:)
  endtype

  type t_absline
    real(r8) :: frac   ! isotope fratction (0.0 ~ 1.0)
    ! total internal partition function table
    integer(i4) :: nq
    real(r8),allocatable :: qtable_t(:) ! temperature (K)
    real(r8),allocatable :: qtable_q(:) ! total internal partition function

    ! absorption lines from HITRAN
    integer(i4) :: nline
    real(r8),allocatable :: wv(:)   ! vacumn wavenumber (cm^-1)
    real(r8),allocatable :: s(:)    ! intensity at 296K ( cm^-1/(molecule cm^-2) )
    real(r8),allocatable :: gamma_air(:) ! air-broadened half-width at 296K (cm^-1 atm^-1)
    real(r8),allocatable :: gamma_self(:) ! self-broadened half-width at 296K (cm^-1 atm^-1)
    real(r8),allocatable :: wvshift(:) ! air pressure-induced line shift at 296K (cm^-1 atm^-1)
    real(r8),allocatable :: nair(:) ! temperature dependence exponent for gamma_air (N/A)
    real(r8),allocatable :: en(:) ! lower-state energy (cm^-1)
  endtype

  type t_srf
    integer(i4) :: npt
    real(r8),allocatable :: wv(:)   ! wavenumber (cm^-1)
    real(r8),allocatable :: weight(:)  ! weight at the wavenumber
  endtype
    

endmodule

