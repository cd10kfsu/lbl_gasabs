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
  public :: co2, h2o
  public :: calculate_layer_od
  public :: destroy_absline

  type(t_absline), save :: co2 ! CO2/626
  type(t_absline), save :: h2o ! H2O/161

  real(r8), parameter :: tref = 296._r8 ! reference temperature (K)
  real(r8), parameter :: hpa2atm = 1.0_r8/1013.25_r8
  real(r8), parameter :: c2 = 1.4388_r8 ! (cm K)
  real(r8), parameter :: kB = 1.3806488D-23
  real(r8), parameter :: pi = ACOS( -1.0_r8 )


Contains
  
function calculate_layer_od( wv, presl, templ, dz, ngas, gasvmr ) &
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
  real(r8),   intent(in) :: gasvmr(ngas)  ! 1-H2O, 2-CO2
  real(r8)               :: odl
! Local vars
  real(r8) :: gamma0, q0, wv0, f0, s0, kabs0, n0
  real(r8) :: qref
  integer :: i

  odl = 0.0_r8
!---------------------------------------------------------------------
! 1. CO2 absorption
!---------------------------------------------------------------------
!---1.1 CO2_626
  ! find Q at layer temperature & reference tempearture
  q0 = linterp1d( co2%nq, co2%qtable_t, co2%qtable_q, templ )
  qref = linterp1d( co2%nq, co2%qtable_t, co2%qtable_q, tref )
  ! accumlate intensity over all lines
  Do i = 1, co2%nline
     gamma0 = (tref/templ)**co2%nair(i) * &
           ( co2%gamma_air(i) *(1-gasvmr(2))*presl*hpa2atm + &
             co2%gamma_self(i)*gasvmr(2)    *presl*hpa2atm )
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
     n0 = co2%frac*gasvmr(2)*presl*100/kB/templ
     ! accumulate optical depth
     odl = odl + n0*kabs0*1.D-4*dz
  Enddo

!---------------------------------------------------------------------
! 2. H2O absorption
!---------------------------------------------------------------------
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


endfunction


subroutine destroy_absline( gas )
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



endmodule

