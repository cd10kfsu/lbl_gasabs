module gas_absorption
!$$$  program documentation block
!         .           .            .
!  program name: gas_absorption
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
  use type_define, only : i4, r8, t_absline, t_atmos
  use math_func, only : linterp1d

  private
  public :: co2, o3, h2o
  public :: calculate_layer_od
  public :: create_abslines
  public :: destroy_abslines
  public :: read_abslines
  public :: get_num_abslines

  type(t_absline), save :: co2 ! CO2/626
  type(t_absline), save :: o3 ! O3/666
  type(t_absline), save :: h2o ! H2O/161

  real(r8), parameter :: tref = 296._r8 ! reference temperature (K)
  real(r8), parameter :: hpa2atm = 1.0_r8/1013.25_r8
  real(r8), parameter :: c2 = 1.4388_r8 ! (cm K)
  real(r8), parameter :: kB = 1.3806488D-23
  real(r8), parameter :: pi = ACOS( -1.0_r8 )


Contains
  
function calculate_layer_od( wv, presl, templ, dz, ngas, gasvmr, lusegas ) &
                    result ( odl )
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
!
!
!$$$
  Implicit None
! Passed args
  real(r8),   intent(in) :: wv
  real(r8),   intent(in) :: presl  ! (hPa)
  real(r8),   intent(in) :: templ
  real(r8),   intent(in) :: dz
  integer(i4),intent(in) :: ngas
  real(r8),   intent(in) :: gasvmr(ngas)  ! 1-H2O, 2-O3, 3-CO2
  logical,    intent(in) :: lusegas(ngas) ! 
  real(r8)               :: odl
! Local vars
  real(r8) :: gamma0, q0, wv0, f0, s0, kabs0, n0
  real(r8) :: qref
  integer :: i

  odl = 0.0_r8
!---------------------------------------------------------------------
! 1. H2O absorption
!---------------------------------------------------------------------
  If ( abs(gasvmr(1)) > 1.d-10 .and. lusegas(1) ) Then
!---1.1 H2O_161
  ! find Q at layer temperature & reference tempearture
  q0 = linterp1d( h2o%nq, h2o%qtable_t, h2o%qtable_q, templ )
  qref = linterp1d( h2o%nq, h2o%qtable_t, h2o%qtable_q, tref )
  ! accumlate intensity over all lines
  Do i = 1, h2o%nline
     gamma0 = (tref/templ)**h2o%nair(i) * &
           ( h2o%gamma_air(i) *(1-gasvmr(1))*presl*hpa2atm + &
             h2o%gamma_self(i)*gasvmr(1)    *presl*hpa2atm )
     wv0 = h2o%wv(i) + h2o%wvshift(i)*presl*hpa2atm
     ! assuming Lorentz line
     ! line shape
     f0 = gamma0/( gamma0*gamma0+(wv-wv0)*(wv-wv0) )/pi
     ! line strength
     s0 = h2o%s(i) * &
         (qref*exp(-c2*h2o%en(i)/templ)*(1-exp(-c2*h2o%wv(i)/templ)) )/ &
         (q0  *exp(-c2*h2o%en(i)/tref )*(1-exp(-c2*h2o%wv(i)/tref )) )
     ! monochromatic absorption coeff ( 1/(molecule cm^-2) )
     kabs0 = s0 * f0
     ! number density of isotop
     n0 = h2o%frac*gasvmr(1)*presl*100/kB/templ
     ! accumulate optical depth
     odl = odl + n0*kabs0*1.D-4*dz
  Enddo

  Endif

!---------------------------------------------------------------------
! 2. O3 absorption
!---------------------------------------------------------------------
  If ( abs(gasvmr(2)) > 1.d-10 .and. lusegas(2) ) Then
!---2.1 O3_666
  ! find Q at layer temperature & reference tempearture
  q0 = linterp1d( o3%nq, o3%qtable_t, o3%qtable_q, templ )
  qref = linterp1d( o3%nq, o3%qtable_t, o3%qtable_q, tref )
  ! accumlate intensity over all lines
  Do i = 1, o3%nline
     gamma0 = (tref/templ)**o3%nair(i) * &
           ( o3%gamma_air(i) *(1-gasvmr(2))*presl*hpa2atm + &
             o3%gamma_self(i)*gasvmr(2)    *presl*hpa2atm )
     wv0 = o3%wv(i) + o3%wvshift(i)*presl*hpa2atm
     ! assuming Lorentz line
     ! line shape
     f0 = gamma0/( gamma0*gamma0+(wv-wv0)*(wv-wv0) )/pi
     ! line strength
     s0 = o3%s(i) * &
         (qref*exp(-c2*o3%en(i)/templ)*(1-exp(-c2*o3%wv(i)/templ)) )/ &
         (q0  *exp(-c2*o3%en(i)/tref )*(1-exp(-c2*o3%wv(i)/tref )) )
     ! monochromatic absorption coeff ( 1/(molecule cm^-2) )
     kabs0 = s0 * f0
     ! number density of isotop
     n0 = o3%frac*gasvmr(2)*presl*100/kB/templ
     ! accumulate optical depth
     odl = odl + n0*kabs0*1.D-4*dz
  Enddo

  Endif

!---------------------------------------------------------------------
! 3. CO2 absorption
!---------------------------------------------------------------------
  If ( abs(gasvmr(3)) > 1.0d-10 .and. lusegas(3) ) Then
!---3.1 CO2_626
  ! find Q at layer temperature & reference tempearture
  q0 = linterp1d( co2%nq, co2%qtable_t, co2%qtable_q, templ )
  qref = linterp1d( co2%nq, co2%qtable_t, co2%qtable_q, tref )
  ! accumlate intensity over all lines
  Do i = 1, co2%nline
     gamma0 = (tref/templ)**co2%nair(i) * &
           ( co2%gamma_air(i) *(1-gasvmr(3))*presl*hpa2atm + &
             co2%gamma_self(i)*gasvmr(3)    *presl*hpa2atm )
     wv0 = co2%wv(i) + co2%wvshift(i)*presl*hpa2atm
     ! assuming Lorentz line
     ! line shape
     f0 = gamma0/( gamma0*gamma0+(wv-wv0)*(wv-wv0) )/pi
     ! line strength
     s0 = co2%s(i) * &
         (qref*exp(-c2*co2%en(i)/templ)*(1-exp(-c2*co2%wv(i)/templ)) )/ &
         (q0  *exp(-c2*co2%en(i)/tref )*(1-exp(-c2*co2%wv(i)/tref )) )
     ! monochromatic absorption coeff ( 1/(molecule cm^-2) )
     kabs0 = s0 * f0
     ! number density of isotop
     n0 = co2%frac*gasvmr(3)*presl*100/kB/templ
     ! accumulate optical depth
     odl = odl + n0*kabs0*1.D-4*dz
  Enddo

  Endif


endfunction


subroutine destroy_abslines( gas )
  implicit none

  type(t_absline), intent(inout) :: gas

  gas%frac = 0.0_r8
  gas%nq   = 0
  If (Allocated(gas%qtable_t)) Deallocate(gas%qtable_t)
  If (Allocated(gas%qtable_q)) Deallocate(gas%qtable_q)
  gas%nline = 0
  If (Allocated(gas%wv))        Deallocate(gas%wv)
  If (Allocated(gas%s))         Deallocate(gas%s)
  If (Allocated(gas%gamma_air)) Deallocate(gas%gamma_air)
  If (Allocated(gas%gamma_self))Deallocate(gas%gamma_self)
  If (Allocated(gas%wvshift))   Deallocate(gas%wvshift)
  If (Allocated(gas%nair))      Deallocate(gas%nair)
  If (Allocated(gas%en))        Deallocate(gas%en)

endsubroutine


subroutine get_num_abslines(fn, uid, nline)
  implicit none

  character(*),intent(in) :: fn
  integer,     intent(in) :: uid
  integer,     intent(out) :: nline

  integer :: i4tmp, ios
  
  nline = 0
  Open(uid,file=trim(fn),action="read")
  Do while (.true.)
     read(uid,*,iostat=ios) i4tmp
     if (ios /= 0) exit
     nline = nline + 1
  Enddo
  Close(uid)

end subroutine


subroutine create_abslines(gas, nline)
  implicit none

  type(t_absline),intent(inout) :: gas
  integer,intent(in) :: nline

  gas%nline = nline
  Allocate( gas%wv(gas%nline), gas%s(gas%nline), &
            gas%gamma_air(gas%nline), &
            gas%gamma_self(gas%nline), gas%wvshift(gas%nline), &
            gas%nair(gas%nline), gas%en(gas%nline) )

end subroutine

subroutine read_abslines(fn, uid, gas)
  implicit none
  character(*),   intent(in) :: fn
  integer,        intent(in) :: uid
  type(t_absline),intent(inout) :: gas

  integer :: i, i4tmp
  real(r8) :: r8tmp

  Open(uid,file=trim(fn),action="read")
  Do i = 1, gas%nline
     Read(uid,*) i4tmp, gas%wv(i), gas%s(i), r8tmp, gas%gamma_air(i), &
                gas%gamma_self(i), gas%en(i), gas%nair(i), gas%wvshift(i)
     ! check point: gas absorption line
     !Write(300,*) i4tmp, gas%wv(i), gas%s(i), r8tmp, gas%gamma_air(i), &
     !        gas%gamma_self(i), gas%en(i), gas%nair(i), gas%wvshift(i)
  Enddo
  Close(uid)

endsubroutine

endmodule

