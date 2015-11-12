program lbl
!$$$  program documentation block
!         .           .            .
!  program name: lbl
!    programmer: da,cheng        org: umd      date: 2015-Nov-08
!
!  purpose:
!    a simplified line-by-line model for gas absorption. 
!    Using Lorentz line shape, and only considering CO2_626, H2O_161
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
  use type_define, only : i4, r8, &
                           t_atmos, t_absline, t_srf
  use math_func, only : linspace, trpint, linterp1d
  use atmos, only : destroy_atmos
  use gas_absorption, only : co2, h2o, calculate_layer_od, destroy_absline

  implicit none

  type(t_atmos) :: atms
  type(t_srf) :: chan, tmpchan

! vars related to SRF
  real(r8) :: wv_center
  real(r8) :: wv_half
  real(r8) :: wv0
  real(r8) :: wv1

  integer(i4) :: i, k

! vars related to optics
  real(r8),allocatable :: layer_od(:,:)      ! nlayer * nline
  real(r8),allocatable :: level_trans(:,:)   ! nlevel * nline
  real(r8),allocatable :: inst_trans(:) ! nlayer
  real(r8),allocatable :: inst_wf(:) ! nlayer

! temp vars
  integer(i4) :: i4tmp
  real(r8) :: r8tmp
  real(r8),allocatable :: r8tmp2d(:,:)



!---------------------------------------------------------------------
! 1. create & fill in atmosphere structure
!---------------------------------------------------------------------

!---1.1 create atmosphere structure
  atms%nlevel = 43
  atms%nlayer = atms%nlevel - 1
  atms%ngas   = 2

  Allocate( atms%pres(atms%nlevel), atms%z(atms%nlevel), &
            atms%temp(atms%nlevel), atms%gas(atms%nlevel,atms%ngas) )
  Allocate( atms%presl(atms%nlayer), atms%dz(atms%nlayer), &
            atms%templ(atms%nlayer), atms%gasl(atms%nlayer,atms%ngas) )

!---1.2 read in level vars
  Allocate( r8tmp2d(atms%nlevel,6) )

  Open(10,file="../sounding.dat",action="read")
  Read(10,*)
  Do i = 1, atms%nlevel
     Read(10,*) ( r8tmp2d(i,k), k = 1, 6 ) 
  Enddo
  Close(10)

  atms%pres(:)  = r8tmp2d(atms%nlevel:1:-1,1)/100  ! level pres (hPa)
  atms%temp(:)  = r8tmp2d(atms%nlevel:1:-1,2)  ! level temp (K)
  atms%z(:)     = r8tmp2d(atms%nlevel:1:-1,3)  ! level height (m)
  atms%gas(:,1) = r8tmp2d(atms%nlevel:1:-1,4)  ! H2O VMR
  atms%gas(:,2) = r8tmp2d(atms%nlevel:1:-1,6)  ! CO2 VMR
  
  Deallocate( r8tmp2d )
  ! check point: level values
  !Do i = 1, atms%nlevel
  !    Write(100,"(3(F11.2,1X),2(E10.3,1X))") &
  !    atms%pres(i), atms%temp(i), atms%z(i), atms%gas(i,1), atms%gas(i,2)
  !Enddo


!---1.3 intepolate level vars to layers
! TOA                    SFC
!  o----------o-----------o
! lev1       lev2        lev3
!      lay1        lay2
  Do i = 1, atms%nlayer
     atms%presl(i) = ( atms%pres(i) + atms%pres(i+1) )/2
     atms%templ(i)  = ( atms%temp(i) + atms%temp(i+1) )/2
     atms%dz(i)    = atms%z(i) - atms%z(i+1)
     atms%gasl(i,:) = ( atms%gas(i,:) + atms%gas(i+1,:) )/2
     ! check point: interpolated layer values
     !Write(200,"(3(F11.2,1X),2(E10.3,1X))") &
     !atms%presl(i), atms%templ(i), atms%dz(i), atms%gasl(i,1), atms%gasl(i,2)
  Enddo

  Write(6,*) "Finish loading atmospheric profile"

!---------------------------------------------------------------------
! 2. create & fill in absorption line structure
!---------------------------------------------------------------------

!---2.1 create & fill in absorption line stucture for CO2
!----2.1.1 create & fill in line strcture for CO2_626
  !co2%nline = 1212 ! 1X chan1
  !co2%nline = 2623 ! 1X chan2
  co2%nline = 2217 ! 5X chan1
  co2%frac  = .984204_r8
  Allocate( co2%wv(co2%nline), co2%s(co2%nline), &
            co2%gamma_air(co2%nline), &
            co2%gamma_self(co2%nline), co2%wvshift(co2%nline), &
            co2%nair(co2%nline), co2%en(co2%nline) )
  Open(20,file="../co2/CO2_chan1.out",action="read")
  !Open(20,file="../co2/CO2_chan1.txt",action="read")
  Do i = 1, co2%nline
     Read(20,*) i4tmp, co2%wv(i), co2%s(i), r8tmp, co2%gamma_air(i), &
                co2%gamma_self(i), co2%en(i), co2%nair(i), co2%wvshift(i)
     ! check point: co2 absorption line
     !Write(300,*) i4tmp, co2%wv(i), co2%s(i), r8tmp, co2%gamma_air(i), &
     !        co2%gamma_self(i), co2%en(i), co2%nair(i), co2%wvshift(i)
  Enddo
  Close(20)

!----2.1.2 read in total internal partition function table of CO2
!          only considering CO2_626
  co2%nq = 271  ! 70K ~ 340K
  Allocate( co2%qtable_t(co2%nq), co2%qtable_q(co2%nq) )
  Open(30,file="../parsum.dat",action="read")
  Read(30,*)
  Do i = 1, co2%nq
     Read(30,*) co2%qtable_t(i), (r8tmp,k=1,6), co2%qtable_q(i)
     ! check point: co2 Q table
     !Write(301,*) co2%qtable_t(i), co2%qtable_q(i)
  Enddo
  Close(30)

  Write(6,*) "Finish loading absoprtion lines for CO2_626"
         
!---2.2 create & fill in absorption line stucture for H2O
!----2.2.1 create & fill in line strcture for H2O_161
  !h2o%nline = 35 ! 1X chan1
  !h2o%nline = 141 ! 1X chan2
  h2o%nline = 97 ! 5X chan1
  h2o%frac  = .997317_r8
  Allocate( h2o%wv(h2o%nline), h2o%s(h2o%nline), &
            h2o%gamma_air(h2o%nline), &
            h2o%gamma_self(h2o%nline), h2o%wvshift(h2o%nline), &
            h2o%nair(h2o%nline), h2o%en(h2o%nline) )
  Open(40,file="../h2o/H2O_chan1.out",action="read")
  !Open(40,file="../h2o/H2O_chan1.txt",action="read")
  Do i = 1, h2o%nline
     Read(40,*) i4tmp, h2o%wv(i), h2o%s(i), r8tmp, h2o%gamma_air(i), &
                h2o%gamma_self(i), h2o%en(i), h2o%nair(i), h2o%wvshift(i)
     ! check point: h2o absorption line
     !Write(400,*) i4tmp, h2o%wv(i), h2o%s(i), r8tmp, h2o%gamma_air(i), &
     !        h2o%gamma_self(i), h2o%en(i), h2o%nair(i), h2o%wvshift(i)
  Enddo
  Close(40)

!----2.1.2 read in total internal partition function table of H2O
!          only considering H2O_161
  h2o%nq = 271  ! 70K ~ 340K
  Allocate( h2o%qtable_t(h2o%nq), h2o%qtable_q(h2o%nq) )
  Open(50,file="../parsum.dat",action="read")
  Read(50,*)
  Do i = 1, h2o%nq
     Read(50,*) h2o%qtable_t(i), h2o%qtable_q(i)
     ! check point: h2o Q table
     !Write(401,*) h2o%qtable_t(i), h2o%qtable_q(i)
  Enddo
  Close(50)

  Write(6,*) "Finish loading absoprtion lines for H2O_161"

!---------------------------------------------------------------------
! 3. load instrument Spectral Response Function
!---------------------------------------------------------------------
!---Option 1: artificial SRF
  !wv_center = 668.18_r8
  !wv_half = 1.5_r8
  !wv0 = wv_center - wv_half
  !wv1 = wv_center + wv_half
  !chan%npt = 301
  !Call linspace( wv0, wv1, chan%npt, chan%wv )
  !Allocate( chan%weight(chan%npt) )
  !chan%weight = 1
  ! check point: check Spectral Response Function (SRF)
  !Do i = 1, chan%npt
  !   Write(500,*) i, chan%wv(i), chan%weight(i)
  !Enddo

!---Option 2: operational SRF
  !Open(60,file="../SRF/rtcoef_noaa_18_hirs_srf_ch02.txt",action="read")
  !Read(60,*); Read(60,*)  !skip two rows
  !Read(60,*) chan%npt
  !chan%npt = 2*chan%npt -1 
  !chan%npt = chan%npt
  !Write(6,*) "SRF points=", chan%npt
  !Allocate( chan%wv(chan%npt), chan%weight(chan%npt) )
  !Read(60,*) !skip one row
  !Do i = 1, chan%npt, 2
  !   Read(60,*) chan%wv(i), chan%weight(i)
  !   Write(6,*) "wv, response=", chan%wv(i), chan%weight(i)
  !Enddo
  !Do i = 2, chan%npt - 1, 2
  !   chan%wv(i) = ( chan%wv(i-1) + chan%wv(i+1) )/2
  !   chan%weight(i) = ( chan%weight(i-1) + chan%weight(i+1) ) /2
  !Enddo
  !Close(60)

!---Option 3: high-density SRF interpolated from operational SRF
  Open(60,file="../SRF/rtcoef_noaa_18_hirs_srf_ch01.txt",action="read")
  Read(60,*); Read(60,*)  !skip two rows
  Read(60,*) tmpchan%npt
  Write(6,*) "original SRF points=", tmpchan%npt
  Allocate( tmpchan%wv(tmpchan%npt), tmpchan%weight(tmpchan%npt) )
  Read(60,*) !skip one row
  Do i = 1, tmpchan%npt
     Read(60,*) tmpchan%wv(i), tmpchan%weight(i)
     !Write(6,*) "wv, response=", tmpchan%wv(i), tmpchan%weight(i)
  Enddo
  Close(60)
  wv0 = tmpchan%wv(1)
  wv1 = tmpchan%wv(tmpchan%npt)
  chan%npt = 20 * tmpchan%npt 
  Write(6,*) "interpolated SRF points=", chan%npt
  Call linspace( wv0, wv1, chan%npt, chan%wv )
  Allocate( chan%weight(chan%npt) )
  Do i = 1, chan%npt
     chan%weight(i) = linterp1d( tmpchan%npt, tmpchan%wv, tmpchan%weight, &
                      chan%wv(i) )
  Enddo
  Deallocate( tmpchan%weight, tmpchan%wv )

  Write(6,*) "Finish SRF structure"

!---------------------------------------------------------------------
! 4. calculate total optical depth & transmittence of each layer
!---------------------------------------------------------------------

!---4.1 calcualte optical depth
  Allocate( layer_od(atms%nlayer,chan%npt) )
  layer_od = 0.0_r8
  Do i = 1, chan%npt
     if (MOD(i,100)==0) Write(6,*) "compute line ", i, " of total ", chan%npt, " lines"
     Do k = 1, atms%nlayer
        ! TOA                    SFC
        !  o----------o-----------o
        ! lev1       lev2        lev3
        !      lay1        lay2
        layer_od(k,i) = calculate_layer_od( chan%wv(i), atms%presl(k), &
             atms%templ(k), atms%dz(k), atms%ngas, atms%gasl(k,:) )
        !print*, "wv=", chan%wv(i)
        !print*, "presl=", atms%presl(k)
        !print*, "templ,dz, ngas,gas=", atms%templ(k), atms%dz(k), atms%ngas, atms%gasl(k,:)
        !pause 1
     Enddo
  Enddo
  
  !check point: layer OD
  !Do k = 1, atms%nlayer
  !   write(600,"(20(E15.4,1X))") (layer_od(k,i), i=1, 10)
  !Enddo

!---4.2 calculate transmission function assuming nadir looking (mu=1)
  Allocate( level_trans(atms%nlevel,chan%npt) )
  Do i = 1, chan%npt
     level_trans(1,i) = 1._r8
     Do k = 2, atms%nlevel
        level_trans(k,i) = exp( -SUM(layer_od(1:k-1,i) ) )
     Enddo
  Enddo

  !check point: layer transmission function
  !Do k = 1, atms%nlevel
  !   write(700,"(20(E15.4,1X))") (level_trans(k,i), i=1, 10)
  !Enddo

  Write(6,*) "Finish line OD & trans calculation"
     
!---------------------------------------------------------------------
! 5. convolve transmission function with channel SRF
!---------------------------------------------------------------------
  Allocate( inst_trans(atms%nlevel) )
  inst_trans(1) = 1.0_r8
  Do k = 2, atms%nlevel
     inst_trans(k) = &
       trpint(chan%npt, chan%wv(:), chan%weight(:)*level_trans(k,:))/ &
       trpint(chan%npt, chan%wv(:), chan%weight(:))
  Enddo

  Write(6,*) "Finish convolving SFR"

!---------------------------------------------------------------------
! 6. Calculate weighting function
!---------------------------------------------------------------------
  Allocate( inst_wf(atms%nlayer) )
  Do k = 1, atms%nlayer
     ! TOA                    SFC
     !  o----------o-----------o
     ! lev1       lev2        lev3
     !       lay1        lay2
     ! trans1     tans2       trans3
     !        od1         od2
     inst_wf(k) = ABS(inst_trans(k+1)-inst_trans(k)) / &
                 ( LOG(atms%pres(k+1)/atms%pres(k) ) )
     write(1000,*) atms%presl(k), inst_wf(k)
  Enddo

  Write(6,*) "Finish WF calculation"

!---------------------------------------------------------------------
! 7. destory all temporary vars
!---------------------------------------------------------------------
   Call destroy_atmos( atms )
   Call destroy_absline( co2 )
   Call destroy_absline( h2o )
   Deallocate( chan%wv, chan%weight )
   Deallocate( layer_od, level_trans )
   Deallocate( inst_trans, inst_wf )

endprogram
