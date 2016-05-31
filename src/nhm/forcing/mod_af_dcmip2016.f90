!-------------------------------------------------------------------------------
!>
!! Module DCMIP2016 Physics forcing
!!
!! @par Description
!!         This module contains subroutines for physics for DCMIP2016.
!!
!! @author R.Yoshida
!!
!! @par History
!! @li      2016-04-28 (R.Yoshida)  [NEW]
!!
!<
module mod_af_dcmip2016
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: af_dcmip2016_init
  public :: af_dcmip2016

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  logical, private :: USE_SimpleMicrophys = .false.
  logical, private :: SM_Latdepend_SST    = .false.
  logical, private :: SM_LargeScaleCond   = .false.
  logical, private :: SM_PBL_Bryan        = .false.
  logical, private :: USE_Kessler         = .false.
  logical, private :: USE_ToyChemistry    = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_dcmip2016_init
    use mod_adm, only: &
       ADM_CTL_FID,  &
       ADM_proc_stop
    use mod_runconf, only: &
       CHEM_TYPE, &
       NCHEM_MAX, &
       NCHEM_STR, &
       NCHEM_END
    implicit none

    logical :: SET_RJ2012       = .false.
    logical :: SET_DCMIP2016_11 = .false.
    logical :: SET_DCMIP2016_12 = .false.
    logical :: SET_DCMIP2016_13 = .false.
    logical :: SET_DCMIP2016_DRY = .false.

    namelist /FORCING_DCMIP_PARAM/ &
       SET_RJ2012,          &
       SET_DCMIP2016_11,    &
       SET_DCMIP2016_12,    &
       SET_DCMIP2016_13,    &
       SET_DCMIP2016_DRY,   &
       USE_SimpleMicrophys, &
       SM_Latdepend_SST,    &
       SM_LargeScaleCond,   &
       SM_PBL_Bryan,        &
       USE_Kessler,         &
       USE_ToyChemistry

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[af_dcmip2016]/Category[nhm forcing]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=FORCING_DCMIP_PARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** FORCING_DCMIP_PARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist FORCING_DCMIP_PARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist FORCING_DCMIP_PARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=FORCING_DCMIP_PARAM)

    ! overwrite setting
    if ( SET_RJ2012 ) then
       write(ADM_LOG_FID,*) '*** Force setting of Reed and Jablonowski (2012)'
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .true.
       SM_PBL_Bryan        = .false.
       USE_Kessler         = .false.
       USE_ToyChemistry    = .false.
    elseif( SET_DCMIP2016_11 ) then
       write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016 Case 1-1 (Moist baroclinic wave with terminator chemistry)'
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       if ( SET_DCMIP2016_DRY ) then
          write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016: DRY condition'
          USE_Kessler      = .false.
       else
          USE_Kessler      = .true.
       endif
       USE_ToyChemistry    = .true.
    elseif( SET_DCMIP2016_12 ) then
       write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016 Case 1-2 (Idealized tropical cyclone)'
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       if ( SET_DCMIP2016_DRY ) then
          write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016: DRY condition'
          USE_Kessler      = .false.
       else
          USE_Kessler      = .true.
       endif
       USE_ToyChemistry    = .false.
    elseif( SET_DCMIP2016_13 ) then
       write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016 Case 1-3 (Mesoscale storm)'
       USE_SimpleMicrophys = .false.
       SM_Latdepend_SST    = .false.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       if ( SET_DCMIP2016_DRY ) then
          write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016: DRY condition'
          USE_Kessler      = .false.
       else
          USE_Kessler      = .true.
       endif
       USE_ToyChemistry    = .false.
    endif

    ! initial value of the tracer is set in mod_prgvar - mod_ideal_init

    if ( USE_ToyChemistry ) then
       if ( CHEM_TYPE == 'PASSIVE' ) then
          if ( NCHEM_MAX /= 2 ) then
             write(*,          *) 'xxx Not appropriate number of passive tracer. STOP.', NCHEM_MAX
             write(ADM_LOG_FID,*) 'xxx Not appropriate number of passive tracer. STOP.', NCHEM_MAX
             call ADM_proc_stop
          endif
       else
          write(*,          *) 'xxx CHEM_TYPE must be set to PASSIVE. STOP.', trim(CHEM_TYPE)
          write(ADM_LOG_FID,*) 'xxx CHEM_TYPE must be set to PASSIVE. STOP.', trim(CHEM_TYPE)
          call ADM_proc_stop
       endif
    endif

    return
  end subroutine af_dcmip2016_init

  !-----------------------------------------------------------------------------
  subroutine af_dcmip2016( &
       ijdim,   &
       lat,     &
       lon,     &
       alt,     &
       rho,     &
       pre,     &
       tem,     &
       vx,      &
       vy,      &
       vz,      &
       q,       &
       ein,     &
       pre_sfc, &
       fvx,     &
       fvy,     &
       fvz,     &
       fe,      &
       fq,      &
       precip,  &
       ix,      &
       iy,      &
       iz,      &
       jx,      &
       jy,      &
       jz,      &
       dt       )
    use mod_adm, only: &
       vlayer => ADM_vlayer, &
       kdim   => ADM_kall,   &
       kmin   => ADM_kmin,   &
       kmax   => ADM_kmax
    use mod_cnst, only: &
       d2r   => CNST_D2R,  &
       Rdry  => CNST_RAIR, &
       CPdry => CNST_CP,   &
       CVdry => CNST_CV,   &
       PRE00 => CNST_PRE00
    use mod_runconf, only: &
       TRC_VMAX,  &
       RAIN_TYPE, &
       I_QV,      &
       I_QC,      &
       I_QR,      &
       NCHEM_STR, &
       NCHEM_END, &
       CVW
    use Terminator, only: &
       tendency_Terminator
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: lat    (ijdim)
    real(RP), intent(in)  :: lon    (ijdim)
    real(RP), intent(in)  :: alt    (ijdim,kdim)
    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: vx     (ijdim,kdim)
    real(RP), intent(in)  :: vy     (ijdim,kdim)
    real(RP), intent(in)  :: vz     (ijdim,kdim)
    real(RP), intent(in)  :: q      (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: ein    (ijdim,kdim)
    real(RP), intent(in)  :: pre_sfc(ijdim)
    real(RP), intent(out) :: fvx    (ijdim,kdim)
    real(RP), intent(out) :: fvy    (ijdim,kdim)
    real(RP), intent(out) :: fvz    (ijdim,kdim)
    real(RP), intent(out) :: fe     (ijdim,kdim)
    real(RP), intent(out) :: fq     (ijdim,kdim,TRC_VMAX)
    real(RP), intent(out) :: precip (ijdim)
    real(RP), intent(in)  :: ix     (ijdim)
    real(RP), intent(in)  :: iy     (ijdim)
    real(RP), intent(in)  :: iz     (ijdim)
    real(RP), intent(in)  :: jx     (ijdim)
    real(RP), intent(in)  :: jy     (ijdim)
    real(RP), intent(in)  :: jz     (ijdim)
    real(RP), intent(in)  :: dt

    ! for kessler
    real(RP) :: theta(vlayer) ! potential temperature (K)
    real(RP) :: qv   (vlayer) ! water vapor mixing ratio (gm/gm)
    real(RP) :: qc   (vlayer) ! cloud water mixing ratio (gm/gm)
    real(RP) :: qr   (vlayer) ! rain  water mixing ratio (gm/gm)
    real(RP) :: rhod (vlayer) ! dry air density (not mean state as in KW) (kg/m^3)
    real(RP) :: pk   (vlayer) ! Exner function (p/p0)**(R/cp)
    real(RP) :: z    (vlayer) ! heights of thermo. levels in the grid column (m)
    real(RP) :: qd   (vlayer)
    real(RP) :: cv   (vlayer)

    ! for simple physics
    real(RP) :: t    (ijdim,vlayer)   ! Temperature at full-model level (K)
    real(RP) :: qvv  (ijdim,vlayer)   ! Specific Humidity at full-model level (kg/kg)
    real(RP) :: u    (ijdim,vlayer)   ! Zonal wind at full-model level (m/s)
    real(RP) :: v    (ijdim,vlayer)   ! Meridional wind at full-model level (m/s)
    real(RP) :: pmid (ijdim,vlayer)   ! Pressure is full-model level (Pa)
    real(RP) :: pint (ijdim,vlayer+1) ! Pressure at model interfaces (Pa)
    real(RP) :: pdel (ijdim,vlayer)   ! Layer thickness (Pa)
    real(RP) :: rpdel(ijdim,vlayer)   ! Reciprocal of layer thickness (1/Pa)
    real(RP) :: ps   (ijdim)          ! surface pressure output [dummy]
    integer  :: test

    ! for toy-chemistory
    real(RP) :: lat_deg, lon_deg
    real(RP) :: cl, cl2
    real(RP) :: cl_f, cl2_f

    integer :: ij, k, kk
    !---------------------------------------------------------------------------

    fvx(:,:)   = 0.0_RP
    fvy(:,:)   = 0.0_RP
    fvz(:,:)   = 0.0_RP
    fe (:,:)   = 0.0_RP
    fq (:,:,:) = 0.0_RP

    precip(:) = 0.0_RP

    if ( USE_Kessler ) then
       do ij = 1, ijdim
          qd   (:) = 1.0_RP               &
                   - q(ij,kmin:kmax,I_QV) &
                   - q(ij,kmin:kmax,I_QC) &
                   - q(ij,kmin:kmax,I_QR)

          qv   (:) = q  (ij,kmin:kmax,I_QV) / qd(:)
          qc   (:) = q  (ij,kmin:kmax,I_QC) / qd(:)
          qr   (:) = q  (ij,kmin:kmax,I_QR) / qd(:)
          rhod (:) = rho(ij,kmin:kmax)      * qd(:)

          pk   (:) = ( pre(ij,kmin:kmax) / PRE00 )**( Rdry / CPdry )
          theta(:) = tem(ij,kmin:kmax) / pk(:)
          z    (:) = alt(ij,kmin:kmax)

          call kessler( theta(:),  & ! [INOUT]
                        qv   (:),  & ! [INOUT]
                        qc   (:),  & ! [INOUT]
                        qr   (:),  & ! [INOUT]
                        rhod (:),  & ! [INOUT] but not changed
                        pk   (:),  & ! [IN]
                        dt,        & ! [IN]
                        z    (:),  & ! [IN]
                        vlayer,    & ! [IN]
                        precip(ij) ) ! [INOUT]

          qd(:) = 1.0_RP &
                / ( 1.0_RP + qv(:) + qc(:) + qr(:) )
          qv(:) = qv(:) * qd(:)
          qc(:) = qc(:) * qd(:)
          qr(:) = qr(:) * qd(:)

          cv(:) = qd(:) * CVdry     &
                + qv(:) * CVW(I_QV) &
                + qc(:) * CVW(I_QC) &
                + qr(:) * CVW(I_QR)

          fq(ij,kmin:kmax,I_QV) = ( qv(:) - q(ij,kmin:kmax,I_QV) ) / dt
          fq(ij,kmin:kmax,I_QC) = ( qc(:) - q(ij,kmin:kmax,I_QC) ) / dt
          fq(ij,kmin:kmax,I_QR) = ( qr(:) - q(ij,kmin:kmax,I_QR) ) / dt
          fe(ij,kmin:kmax)      = ( cv(:) * theta(:) * pk(:) - ein(ij,kmin:kmax) ) / dt
       enddo
    endif

    if ( USE_SimpleMicrophys ) then
       if ( SM_Latdepend_SST ) then
          test = 1
       else
          test = 0
       endif

       do k = 1, vlayer
          kk = kmax - k + 1 ! reverse order

          t   (:,k) = tem(:,kk)
          qvv (:,k) = q  (:,kk,I_QV)
          u   (:,k) = vx (:,kk) * ix(:) &
                    + vy (:,kk) * iy(:) &
                    + vz (:,kk) * iz(:)
          v   (:,k) = vx (:,kk) * jx(:) &
                    + vy (:,kk) * jy(:) &
                    + vz (:,kk) * jz(:)
          pmid(:,k) = pre(:,kk)
       enddo

       pint(:,1) = 0.0_RP
       do k = 2, vlayer
          pint(:,k) = 0.5_RP * ( pmid(:,k-1) + pmid(:,k) )
       enddo
       pint(:,vlayer+1) = pre_sfc(:)

       do k = 1, vlayer
          pdel (:,k) = pint(:,k) - pint(:,k+1)
          rpdel(:,k) = 1.0_RP / pdel(:,k)
       enddo

       ps(:) = pre_sfc(:)

       call simple_physics( ijdim,             & ! [IN]
                            vlayer,            & ! [IN]
                            dt,                & ! [IN]
                            lat   (:),         & ! [IN]
                            t     (:,:),       & ! [INOUT]
                            qvv   (:,:),       & ! [INOUT]
                            u     (:,:),       & ! [INOUT]
                            v     (:,:),       & ! [INOUT]
                            pmid  (:,:),       & ! [INOUT] but not changed
                            pint  (:,:),       & ! [INOUT] but not changed
                            pdel  (:,:),       & ! [INOUT] but not changed
                            rpdel (:,:),       & ! [INOUT] but not changed
                            ps    (:),         & ! [INOUT] but not changed
                            precip(:),         & ! [INOUT]
                            test,              & ! [IN]
                            SM_LargeScaleCond, & ! [IN]
                            SM_PBL_Bryan       ) ! [IN]

       do k = 1, vlayer
          kk = kmax - k + 1 ! reverse order

          fvx(:,kk) = ( u(:,k) * ix(:) + v(:,k) * jx(:) - vx(:,kk) ) / dt
          fvy(:,kk) = ( u(:,k) * iy(:) + v(:,k) * jy(:) - vy(:,kk) ) / dt
          fvz(:,kk) = ( u(:,k) * iz(:) + v(:,k) * jz(:) - vz(:,kk) ) / dt
       enddo

       do k = 1, vlayer
          kk = kmax - k + 1 ! reverse order

          do ij = 1, ijdim
             if    ( RAIN_TYPE == 'DRY' ) then
                qv(k) = qvv(ij,k)

                qd(k) = 1.0_RP &
                      - qv(k)

                cv(k) = qd(k) * CVdry     &
                      + qv(k) * CVW(I_QV)
             elseif( RAIN_TYPE == 'WARM' ) then
                qv(k) = qvv(ij,k)
                qc(k) = q  (ij,kk,I_QC)
                qr(k) = q  (ij,kk,I_QR)

                qd(k) = 1.0_RP &
                      - qv(k)  &
                      - qc(k)  &
                      - qr(k)

                cv(k) = qd(k) * CVdry     &
                      + qv(k) * CVW(I_QV) &
                      + qc(k) * CVW(I_QC) &
                      + qr(k) * CVW(I_QR)
             endif

             fq(ij,kk,I_QV) = ( qv(k) - q(ij,kk,I_QV) ) / dt
             fe(ij,kk)      = ( cv(k) * t(ij,k) - ein(ij,kk) ) / dt
          enddo
       enddo

    endif

    if ( USE_ToyChemistry ) then
       do k  = kmin, kmax
       do ij = 1,    ijdim
          lat_deg = lat(ij) / d2r
          lon_deg = lon(ij) / d2r

          cl  = q(ij,k,NCHEM_STR)
          cl2 = q(ij,k,NCHEM_END)

          call tendency_Terminator( lat_deg, lon_deg, cl, cl2, dt, cl_f, cl2_f )

          fq(ij,k,NCHEM_STR) = cl_f
          fq(ij,k,NCHEM_END) = cl2_f
       enddo
       enddo
    endif

  end subroutine af_dcmip2016

end module mod_af_dcmip2016
!-------------------------------------------------------------------------------
