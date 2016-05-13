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

    namelist /FORCING_DCMIP_PARAM/ &
       SET_RJ2012,          &
       SET_DCMIP2016_11,    &
       SET_DCMIP2016_12,    &
       SET_DCMIP2016_13,    &
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
       USE_Kessler         = .true.
       USE_ToyChemistry    = .true.
    elseif( SET_DCMIP2016_12 ) then
       write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016 Case 1-2 (Idealized tropical cyclone)'
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       USE_Kessler         = .true.
       USE_ToyChemistry    = .false.
    elseif( SET_DCMIP2016_13 ) then
       write(ADM_LOG_FID,*) '*** Force setting of DCMIP2016 Case 1-3 (Mesoscale storm)'
       USE_SimpleMicrophys = .false.
       SM_Latdepend_SST    = .false.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       USE_Kessler         = .true.
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
       ijdim, &
       lat,   &
       lon,   &
       pre,   &
       tem,   &
       vx,    &
       vy,    &
       vz,    &
       q,     &
       fvx,   &
       fvy,   &
       fvz,   &
       fe,    &
       fq,    &
       dt     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       d2r => CNST_D2R
    use mod_runconf, only: &
       TRC_VMAX,  &
       NCHEM_STR, &
       NCHEM_END
    use Terminator, only: &
       tendency_Terminator
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: lat(ijdim)
    real(RP), intent(in)  :: lon(ijdim)
    real(RP), intent(in)  :: pre(ijdim,kdim)
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(in)  :: vx (ijdim,kdim)
    real(RP), intent(in)  :: vy (ijdim,kdim)
    real(RP), intent(in)  :: vz (ijdim,kdim)
    real(RP), intent(in)  :: q  (ijdim,kdim,TRC_VMAX)
    real(RP), intent(out) :: fvx(ijdim,kdim)
    real(RP), intent(out) :: fvy(ijdim,kdim)
    real(RP), intent(out) :: fvz(ijdim,kdim)
    real(RP), intent(out) :: fe (ijdim,kdim)
    real(RP), intent(out) :: fq (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: dt

    real(RP) :: preh (ijdim,kdim+1)
    real(RP) :: dpre (ijdim,kdim)
    real(RP) :: rdpre(ijdim,kdim)

    real(RP) :: lat_deg, lon_deg
    real(RP) :: cl, cl2
    real(RP) :: cl_f, cl2_f

    integer :: n
    !---------------------------------------------------------------------------

    fvx(:,:)   = 0.0_RP
    fvy(:,:)   = 0.0_RP
    fvz(:,:)   = 0.0_RP
    fe (:,:)   = 0.0_RP
    fq (:,:,:) = 0.0_RP

    if ( USE_Kessler ) then

    endif

    if ( USE_SimpleMicrophys ) then

    endif

    if ( USE_ToyChemistry ) then
       do n = 1, ijdim
          lat_deg = lat(n) / d2r
          lon_deg = lon(n) / d2r
          cl  = q(n,kmin,NCHEM_STR)
          cl2 = q(n,kmin,NCHEM_END)

          call tendency_Terminator( lat_deg, lon_deg, cl, cl2, dt, cl_f, cl2_f )

          fq(n,:,NCHEM_STR) = cl_f
          fq(n,:,NCHEM_END) = cl2_f
       enddo
    endif

  end subroutine af_dcmip2016

end module mod_af_dcmip2016
!-------------------------------------------------------------------------------
