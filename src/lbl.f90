program lbl
!$$$  program documentation block
!         .           .            .
!  program name: lbl
!    programmer: da,cheng        org: umd      date: 2015-Nov-08
!
!  purpose:
!    a simplified line-by-line model for gas absorption. 
!    (1) use Lorentz line shape
!    (2) only consider CO2_626, H2O_161, O3_666
!  
!    in order to run the program, you need to prepare 
!    (1) files containing absorption lines for the above gases in 
!        directory ../co2 ../h2o ../o3
!    (2) sounding data in ..
!    (3) instrument spectral response function in ../SRF (optional)
!
!    the calculated weighting function is written in fort.1000
!    column 1: layer pressure (hPa)
!    column 2: dTrans/d(lnp)
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
  use type_define, only : i4, r8, slen, &
                           t_atmos, t_absline, t_srf
  use math_func, only : linspace, trpint, linterp1d
  use atmos, only : create_atmos, destroy_atmos
  use gas_absorption, only : co2, h2o, o3, &
                             calculate_layer_od, destroy_abslines, &
                             get_num_abslines, read_abslines, create_abslines
#ifdef OMP
  use omp_lib
#endif
  implicit none

  type(t_atmos) :: atms
  type(t_srf) :: chan, tmpchan

! vars related to SRF
  real(r8) :: wv_center
  real(r8) :: wv_half
  real(r8) :: wv0
  real(r8) :: wv1
  integer  :: fac_oversample

  integer(i4) :: i, k

! vars related to optics
  real(r8),allocatable :: layer_od(:,:)      ! nlayer * nline
  real(r8),allocatable :: level_trans(:,:)   ! nlevel * nline
  real(r8),allocatable :: inst_trans(:) ! nlayer
  real(r8),allocatable :: inst_wf(:) ! nlayer
  logical :: use_o3, use_h2o, use_co2

! vars to record run time
  real(r8) :: tic, toc

! vars related to input files
  character(len=slen) :: fn_sounding, fn_qtable, fn_srf
  character(len=slen) :: fn_gas_h2o_161, fn_gas_co2_626, fn_gas_o3_666

! temp vars
  real(r8) :: r8tmp
  real(r8),allocatable :: r8tmp2d(:,:)

  namelist /config/ fn_sounding, fn_qtable, fn_srf, &
                    fn_gas_h2o_161, fn_gas_co2_626, fn_gas_o3_666, &
                    use_h2o, use_o3, use_co2, &
                    fac_oversample

!---------------------------------------------------------------------
! 0. config
!---------------------------------------------------------------------
  fn_sounding    = "../sounding.dat"
  fn_qtable      = "../parsum.dat"
  fn_srf         = "../SRF/rtcoef_noaa_18_hirs_srf_ch01.txt"
  fn_gas_h2o_161 = "../h2o/H2O_chan1.txt"
  fn_gas_co2_626 = "../co2/CO2_chan1.txt"
  fn_gas_o3_666  = "../o3/O3_chan1.out"

  use_h2o = .true.; use_co2 = .true.; use_o3 = .false.

  fac_oversample = 30
  write(*,nml=config)

! set timer 
#ifdef OMP
  print*, "enable OpenMP: nthreads = ", omp_get_max_threads()
#else
  print*, "use serial version"
#endif
  tic = omp_get_wtime()


! read in atmosphere/absorption line/SRF files in parallel if openmp enabled
!$omp parallel 
!$omp sections

!---------------------------------------------------------------------
! 1. create & fill in atmosphere structure
!---------------------------------------------------------------------
!$omp section

!---1.1 create atmosphere structure
  ! create atmos structure with 43 levels (42 layers), 
  ! and 3 kinds of gases
  ! gas 1: H2O, gas 2: O3 gas, 3: CO2
  Call create_atmos( atms, 43, 3 )

!---1.2 read in level vars
  Allocate( r8tmp2d(atms%nlevel,6) )
  ! sounding data is from surface to TOA
  Open(10,file=trim(fn_sounding),action="read")
  Read(10,*)
  Do i = 1, atms%nlevel
     Read(10,*) ( r8tmp2d(i,k), k = 1, 6 ) 
  Enddo
  Close(10)

  ! atmos structure requires input from TOA to sfc
  atms%pres(:)  = r8tmp2d(atms%nlevel:1:-1,1)/100  ! level pres (hPa)
  atms%temp(:)  = r8tmp2d(atms%nlevel:1:-1,2)  ! level temp (K)
  atms%z(:)     = r8tmp2d(atms%nlevel:1:-1,3)  ! level height (m)
  atms%gas(:,1) = r8tmp2d(atms%nlevel:1:-1,4)  ! H2O VMR
  atms%gas(:,2) = r8tmp2d(atms%nlevel:1:-1,5)  ! O3 VMR
  atms%gas(:,3) = r8tmp2d(atms%nlevel:1:-1,6)  ! CO2 VMR

  
  
  Deallocate( r8tmp2d )
  ! check point: level values
  !Do i = 1, atms%nlevel
  !    Write(100,"(3(F11.2,1X),2(E10.3,1X))") &
  !    atms%pres(i), atms%temp(i), atms%z(i), atms%gas(i,1), atms%gas(i,2)
  !Enddo


!---1.3 average level vars to layers
! TOA                    SFC
!  o----------o-----------o
! lev1       lev2        lev3
!      lay1        lay2
  Do i = 1, atms%nlayer
     atms%presl(i) = ( atms%pres(i) + atms%pres(i+1) )/2
     atms%templ(i)  = ( atms%temp(i) + atms%temp(i+1) )/2
     atms%dz(i)    = atms%z(i) - atms%z(i+1)
     atms%gasl(i,:) = ( atms%gas(i,:) + atms%gas(i+1,:) )/2
     ! check point: layer values
     Write(200,"(3(F11.2,1X),4(E10.3,1X))") &
     atms%presl(i), atms%templ(i), atms%dz(i), (atms%gasl(i,k),k=1,atms%ngas)
  Enddo

  Write(6,*) "Finish loading atmospheric profile"

!---------------------------------------------------------------------
! 2. create & fill in absorption line structure
!---------------------------------------------------------------------
!$omp section
!---2.1 create & fill in absorption line stucture for CO2
!----2.1.1 create & fill in line strcture for CO2_626
  If (use_co2) Then
      co2%frac  = .984204_r8
      call get_num_abslines(trim(fn_gas_co2_626), 20, co2%nline)
      call create_abslines(co2,co2%nline)
      call read_abslines(trim(fn_gas_co2_626), 20, co2)

    !----2.1.2 read in total internal partition function table of CO2
    !          only considering CO2_626
      co2%nq = 271  ! 70K ~ 340K
      Allocate( co2%qtable_t(co2%nq), co2%qtable_q(co2%nq) )
      Open(30,file=trim(fn_qtable),action="read")
      Read(30,*)
      Do i = 1, co2%nq
         Read(30,*) co2%qtable_t(i), (r8tmp,k=1,6), co2%qtable_q(i)
         ! check point: co2 Q table
         !Write(301,*) co2%qtable_t(i), co2%qtable_q(i)
      Enddo
      Close(30)

      Write(6,*) "Finish loading absoprtion lines for CO2_626 (", co2%nline, " lines)."
  Else
      Write(6,*) "Skip CO2_626"
  Endif

!$omp section
!---2.2 create & fill in absorption line stucture for H2O
!----2.2.1 create & fill in line strcture for H2O_161
  If (use_h2o) Then
      h2o%frac  = .997317_r8
      call get_num_abslines(trim(fn_gas_h2o_161), 40, h2o%nline)
      call create_abslines(h2o,h2o%nline)
      call read_abslines(trim(fn_gas_h2o_161), 40, h2o)

    !----2.2.2 read in total internal partition function table of H2O
    !          only considering H2O_161
      h2o%nq = 271  ! 70K ~ 340K
      Allocate( h2o%qtable_t(h2o%nq), h2o%qtable_q(h2o%nq) )
      Open(50,file=trim(fn_qtable),action="read")
      Read(50,*)
      Do i = 1, h2o%nq
         Read(50,*) h2o%qtable_t(i), h2o%qtable_q(i)
         ! check point: h2o Q table
         !Write(401,*) h2o%qtable_t(i), h2o%qtable_q(i)
      Enddo
      Close(50)

      Write(6,*) "Finish loading absoprtion lines for H2O_161 (", h2o%nline, " lines)."
  Else
      Write(6,*) "Skip H2O_161."
  Endif

!$omp section
!---2.3 create & fill in absorption line stucture for O3
!----2.3.1 create & fill in line strcture for O3_666
  If (use_o3) Then
      o3%frac  = .992901_r8

      call get_num_abslines(trim(fn_gas_o3_666), 41, o3%nline)
      call create_abslines(o3,o3%nline)
      call read_abslines(trim(fn_gas_o3_666), 41, o3)

    !----2.3.2 read in total internal partition function table of O3
    !          only considering O3_666
      o3%nq = 271  ! 70K ~ 340K
      Allocate( o3%qtable_t(o3%nq), o3%qtable_q(o3%nq) )
      Open(51,file=trim(fn_qtable),action="read")
      Read(51,*)
      Do i = 1, o3%nq
         Read(51,*) o3%qtable_t(i), (r8tmp,k=1,14), o3%qtable_q(i)
         ! check point: o3 Q table
         Write(510,*) o3%qtable_t(i), o3%qtable_q(i)
      Enddo
      Close(51)

      Write(6,*) "Finish loading absoprtion lines for O3_666 (", o3%nline, " lines)."
  Else
      Write(6,*) "Skip O3_666."
  Endif

!---------------------------------------------------------------------
! 3. load instrument Spectral Response Function
!---------------------------------------------------------------------
!$omp section
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
  !Open(60,file=trim(fn_srf),action="read")
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
  Open(60,file=trim(fn_srf),action="read")
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
  chan%npt = fac_oversample * tmpchan%npt 
  Write(6,*) "interpolated SRF points=", chan%npt
  Call linspace( wv0, wv1, chan%npt, chan%wv )
  Allocate( chan%weight(chan%npt) )
  Do i = 1, chan%npt
     chan%weight(i) = linterp1d( tmpchan%npt, tmpchan%wv, tmpchan%weight, &
                      chan%wv(i) )
  Enddo
  Deallocate( tmpchan%weight, tmpchan%wv )

  Write(6,*) "Finish SRF structure"

!$omp end sections
!$omp end parallel
!---------------------------------------------------------------------
! 4. calculate total optical depth & transmittence of each layer
!---------------------------------------------------------------------

!---4.1 calcualte optical depth
  Allocate( layer_od(atms%nlayer,chan%npt) )
  layer_od = 0.0_r8
  !Do i = 1, chan%npt
  !   if (MOD(i,100)==0) Write(6,*) "compute line ", i, " of total ", chan%npt, " lines"
  !   Do k = 1, atms%nlayer


!$omp parallel private(i,k) &
!$omp shared(chan, layer_od, atms)
!$omp do
  Do i = 1, chan%npt
#ifndef OMP
     if (mod(i,100)==0) print*, "finish ", i, " of ", chan%npt, " lines"
#endif
     Do k = 1, atms%nlayer
        !print*, "k=", k

        ! TOA                    SFC
        !  o----------o-----------o
        ! lev1       lev2        lev3
        !      lay1        lay2
        layer_od(k,i) = calculate_layer_od( chan%wv(i), atms%presl(k), &
             atms%templ(k), atms%dz(k), atms%ngas, atms%gasl(k,:), (/use_h2o,use_o3,use_co2/) )
        !print*, "wv=", chan%wv(i)
        !print*, "presl=", atms%presl(k)
        !print*, "templ,dz, ngas,gas=", atms%templ(k), atms%dz(k), atms%ngas, atms%gasl(k,:)
        !pause 1
     Enddo
  Enddo
!$omp end do
!$omp end parallel

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
     !print*, "k, tau:",k, inst_trans(k)
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
   Call destroy_abslines( co2 )
   Call destroy_abslines( h2o )
   Call destroy_abslines( o3 )
   Deallocate( chan%wv, chan%weight )
   Deallocate( layer_od, level_trans )
   Deallocate( inst_trans, inst_wf )

   toc = omp_get_wtime()

   Write(6,*) "run time=", toc-tic, "seconds."

#ifndef OMP
 contains
    function omp_get_wtime() result(t)
      real(r8) :: t
      Call cpu_time(t)
    end function
#endif

endprogram
