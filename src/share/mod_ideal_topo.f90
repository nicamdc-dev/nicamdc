!-------------------------------------------------------------------------------
module mod_ideal_topo
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for the Dyn Core Test Initialization.
  !
  !
  !++ Current Corresponding Author : R.Yoshida
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      12-10-19   Imported from mod_restart.f90 of NICAM
  !      0.01      13-06-12   Test cases in DCMIP2012 were imported
  !      0.02      14-01-28   Test case in Tomita and Satoh 2004 was imported
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IDEAL_topo

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: IDEAL_topo_JW
  private :: IDEAL_topo_Schar_Moderate
  private :: IDEAL_topo_Schar_Steep

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine IDEAL_topo( &
       lat, &
       lon, &
       Zsfc )
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_KNONE,     &
       ADM_gall,      &
       ADM_lall
    implicit none

    real(8), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    character(len=ADM_NSYS) :: topo_type = ''

    namelist / IDEALTOPOPARAM / &
       topo_type

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[ideal topo]/Category[nhm ideal]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=IDEALTOPOPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** IDEALTOPOPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist IDEALTOPOPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist IDEALTOPOPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=IDEALTOPOPARAM)

    if    ( topo_type == 'Schar_Moderate' ) then

       call IDEAL_topo_Schar_Moderate( lat (:,:,:), & !--- [IN]
                                       lon (:,:,:), & !--- [IN]
                                       Zsfc(:,:,:)  ) !--- [OUT]

    elseif( topo_type == 'Schar_Steep' ) then

       call IDEAL_topo_Schar_Steep( lat (:,:,:), & !--- [IN]
                                    lon (:,:,:), & !--- [IN]
                                    Zsfc(:,:,:)  ) !--- [OUT]

    elseif( topo_type == 'JW' ) then

       call IDEAL_topo_JW( lat (:,:,:), & !--- [IN]
                           lon (:,:,:), & !--- [IN]
                           Zsfc(:,:,:)  ) !--- [OUT]

    else
       write(*,          *) 'xxx Not appropriate topo_type. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate topo_type. STOP.'
       call ADM_proc_stop
    endif

    return
  end subroutine IDEAL_topo

  !-----------------------------------------------------------------------------
  !> Moderately-steep Schar-like circular mountain (Ref.: DCMIP 2012 eq.(48))
  subroutine IDEAL_topo_Schar_Moderate( &
       lat, &
       lon, &
       Zsfc )
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_KNONE,     &
       ADM_lall,      &
       ADM_gall
    use mod_cnst, only: &
       PI     => CNST_PI,      &
       D2R    => CNST_D2R,     &
       RADIUS => CNST_ERADIUS
    implicit none

    real(8), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: center_lon =  270.D0 ! Longitude of Schar-type mountain center point [deg]
    real(8) :: center_lat =    0.D0 ! Latitude  of Schar-type mountain center point [deg]
    real(8) :: H0         = 2000.D0 ! Maximum Schar-type mountain height [m]
    real(8) :: Rm_deg     =  135.D0 ! Schar-type mountain radius     [deg]
    real(8) :: QSIm_deg   = 11.25D0 ! Schar-type mountain wavelength [deg]

    namelist / IDEALTOPOPARAM_Schar_Moderate / &
       center_lon, &
       center_lat, &
       H0,         &
       Rm_deg,     &
       QSIm_deg

    real(8) :: LAMBDA,  PHI
    real(8) :: LAMBDAm, PHIm, Rm, QSIm
    real(8) :: sinPHIm, cosPHIm
    real(8) :: distance, mask

    integer :: g, l, K0
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[topo Schar Moderate]/Category[nhm ideal]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=IDEALTOPOPARAM_Schar_Moderate,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** IDEALTOPOPARAM_Schar_Moderate is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist IDEALTOPOPARAM_Schar_Moderate. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist IDEALTOPOPARAM_Schar_Moderate. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=IDEALTOPOPARAM_Schar_Moderate)

    K0 = ADM_KNONE

    LAMBDAm = center_lon * D2R ! [deg]->[rad]
    PHIm    = center_lat * D2R ! [deg]->[rad]
    Rm      = Rm_deg     * D2R ! [deg]->[rad]
    QSIm    = QSIm_deg   * D2R ! [deg]->[rad]
    sinPHIm = sin(PHIm)
    cosPHIm = cos(PHIm)

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       LAMBDA = lon(g,K0,l)
       PHI    = lat(g,K0,l)

       distance = acos( sinPHIm * sin(PHI)                       &
                      + cosPHIm * cos(PHI) * cos(LAMBDA-LAMBDAm) )

       mask = 0.5D0 - sign(0.5D0,distance-Rm) ! if distance > Rm, mask = 0

       Zsfc(g,ADM_KNONE,l) = H0/2.D0                              &
                           * ( 1.D0 + cos( PI * distance / Rm ) ) &
                           * cos( PI * distance / QSIm )**2       &
                           * mask
    enddo
    enddo

    return
  end subroutine IDEAL_topo_Schar_Moderate

  !-----------------------------------------------------------------------------
  !> Steep Schar-type mountain (Ref.: DCMIP 2012 eq.(76))
  subroutine IDEAL_topo_Schar_Steep( &
       lat, &
       lon, &
       Zsfc )
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_KNONE,     &
       ADM_lall,      &
       ADM_gall
    use mod_cnst, only: &
       PI     => CNST_PI,      &
       D2R    => CNST_D2R,     &
       RADIUS => CNST_ERADIUS
    implicit none

    real(8), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: center_lon =   45.D0 ! Longitude of Schar-type mountain center point [deg]
    real(8) :: center_lat =    0.D0 ! Latitude  of Schar-type mountain center point [deg]
    real(8) :: H0         =  250.D0 ! Maximum Schar-type mountain height [m]
    real(8) :: d          = 5000.D0 ! Schar-type mountain half-width [m]
    real(8) :: QSI        = 4000.D0 ! Schar-type mountain wavelength [m]

    namelist / IDEALTOPOPARAM_Schar_Steep / &
       center_lon, &
       center_lat, &
       H0,         &
       d,          &
       QSI

    real(8) :: LAMBDA,  PHI
    real(8) :: LAMBDAc, PHIc
    real(8) :: sinPHIc, cosPHIc
    real(8) :: distance

    integer :: g, l, K0
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[topo Schar Steep]/Category[nhm ideal]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=IDEALTOPOPARAM_Schar_Steep,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** IDEALTOPOPARAM_Schar_Steep is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist IDEALTOPOPARAM_Schar_Steep. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist IDEALTOPOPARAM_Schar_Steep. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=IDEALTOPOPARAM_Schar_Steep)

    K0 = ADM_KNONE

    LAMBDAc = center_lon * D2R ! [deg]->[rad]
    PHIc    = center_lat * D2R ! [deg]->[rad]
    sinPHIc = sin(PHIc)
    cosPHIc = cos(PHIc)

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       LAMBDA = lon(g,K0,l)
       PHI    = lat(g,K0,l)

       distance = RADIUS * acos( sinPHIc * sin(PHI)                       &
                               + cosPHIc * cos(PHI) * cos(LAMBDA-LAMBDAc) )

       Zsfc(g,ADM_KNONE,l) = H0                                  &
                           * exp( -(distance*distance) / (d*d) ) &
                           * cos( PI * distance / QSI )**2
    enddo
    enddo

    return
  end subroutine IDEAL_topo_Schar_Steep

  !-----------------------------------------------------------------------------
  !> mountain for JW06 testcase
  subroutine IDEAL_topo_JW( &
       lat, &
       lon, &
       Zsfc )
    use mod_adm, only: &
       ADM_KNONE, &
       ADM_lall,  &
       ADM_gall
    use mod_cnst, only: &
       PI     => CNST_PI,      &
       RADIUS => CNST_ERADIUS, &
       OHM    => CNST_EOHM,    &
       GRAV   => CNST_EGRAV
    implicit none

    real(8), intent(in)  :: lat (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(in)  :: lon (ADM_gall,ADM_KNONE,ADM_lall)
    real(8), intent(out) :: Zsfc(ADM_gall,ADM_KNONE,ADM_lall)

    real(8), parameter :: ETA0 = 0.252D0 ! Value of eta at a reference level (position of the jet)
    real(8), parameter :: ETAs =    1.D0 ! Value of eta at the surface
    real(8), parameter :: u0   =   35.D0 ! Maximum amplitude of the zonal wind

    real(8) :: PHI, ETAv, u0cos32ETAv
    real(8) :: f1, f2

    integer :: g, l, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    ETAv        = ( ETAs - ETA0 ) * PI/2.D0
    u0cos32ETAv = u0 * cos(ETAv)**(3.D0/2.D0)

    ! for globe
    do l = 1, ADM_lall
    do g = 1, ADM_gall
       PHI = lat(g,K0,l)

       f1 = -2.D0 * sin(PHI)**6 * ( cos(PHI)**2 + 1.D0/3.D0 ) + 10.D0/63.D0
       f2 = 8.D0/5.D0 * cos(PHI)**3 * ( sin(PHI)**2 + 2.D0/3.D0 ) - PI/4.D0

       Zsfc(g,k0,l) = u0cos32ETAv * ( u0cos32ETAv*f1 + RADIUS*OHM*f2 ) / GRAV
    enddo
    enddo

    return
  end subroutine IDEAL_topo_JW

end module mod_ideal_topo
