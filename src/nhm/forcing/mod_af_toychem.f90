!-------------------------------------------------------------------------------
!>
!! Module TOY-CHEM forcing
!!
!! @par Description
!!         This module contains subroutines for initial condition and forcing term of DCMIP2016.
!!         DCMIP2016 Test Case 1-1: Toy chemistory
!!
!! @author NICAM developers
!!
!<
module mod_af_toychem
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
  public :: af_toychem_init
  public :: af_toychem

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_toychem_init
    use mod_adm, only: &
       ADM_proc_stop
    use mod_runconf, only: &
       CHEM_TYPE, &
       NCHEM_MAX, &
       NCHEM_STR, &
       NCHEM_END
    implicit none
    !---------------------------------------------------------------------------

    ! initial value of the tracer is set in mod_prgvar - mod_ideal_init
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

    return
  end subroutine af_toychem_init

  !-----------------------------------------------------------------------------
  subroutine af_toychem( &
       ijdim, &
       lat,   &
       lon,   &
       q,     &
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
    real(RP), intent(in)  :: q  (ijdim,kdim,TRC_VMAX)
    real(RP), intent(out) :: fq (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: dt

    real(RP) :: lat_deg, lon_deg
    real(RP) :: cl, cl2
    real(RP) :: cl_f, cl2_f

    integer :: n
    !---------------------------------------------------------------------------

    fq(:,:,:) = 0.0_RP

    do n = 1, ijdim
       lat_deg = lat(n) / d2r
       lon_deg = lon(n) / d2r
       cl  = q(n,kmin,NCHEM_STR)
       cl2 = q(n,kmin,NCHEM_END)

       call tendency_Terminator( lat_deg, lon_deg, cl, cl2, dt, cl_f, cl2_f )

       fq(n,:,NCHEM_STR) = cl_f
       fq(n,:,NCHEM_END) = cl2_f
    enddo

  end subroutine af_toychem

end module mod_af_toychem
!-------------------------------------------------------------------------------
