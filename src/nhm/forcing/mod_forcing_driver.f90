!-------------------------------------------------------------------------------
!>
!! Module forcing driver
!!
!! @par Description
!!         This module is for the artificial forcing
!!
!! @author R.Yoshida
!!
!! @par History
!! @li      2012-10-11 (R.Yoshida) [NEW] extract from phystep
!! @li      2013-03-07 (H.Yashiro) marge, refactoring
!!
!<
module mod_forcing_driver
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: forcing_init
  public :: forcing

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
  subroutine forcing_init
    use mod_runconf, only: &
       AF_TYPE
    use mod_af_heldsuarez, only: &
       af_heldsuarez_init
    implicit none

    !---------------------------------------------------------------------------

    select case(AF_TYPE)
    case('NONE')
       !--- do nothing
    case('HELD-SUAREZ')
       call af_heldsuarez_init
    case default
       write(*,*) 'Msg : Sub[af_init]/Mod[af_driver]'
       write(*,*) ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
    end select

    return
  end subroutine forcing_init

  !-----------------------------------------------------------------------------
  subroutine forcing
    use mod_adm, only: &
       ADM_gall_in, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_time, only: &
       TIME_DTL
    use mod_grd, only: &
       GRD_vz, &
       GRD_Z
    use mod_gmtr, only: &
       GMTR_lat
    use mod_vmtr, only: &
       VMTR_GSGAM2,  &
       VMTR_GSGAM2H, &
       VMTR_PHI
    use mod_runconf, only: &
       AF_TYPE, &
       TRC_VMAX
    use mod_prgvar, only: &
       prgvar_get_in_withdiag, &
       prgvar_set_in
    use mod_gtl, only: &
       GTL_clip_region, &
       GTL_clip_region_1layer
    use mod_bndcnd, only: &
       bndcnd_thermo
    use mod_af_heldsuarez, only: &
       af_HeldSuarez
    use mod_history, only: &
       history_in
    implicit none

    real(8) :: rhog  (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogw (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhoge (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rhogq (ADM_gall_in,ADM_kall,ADM_lall,TRC_vmax)
    real(8) :: rho   (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: pre   (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: tem   (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vx    (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vy    (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: vz    (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: w     (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: q     (ADM_gall_in,ADM_kall,ADM_lall,TRC_vmax)

    ! forcing tendency 
    real(8) :: fvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: fvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: fvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: fw (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: fe (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: fq (ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)

    ! geometry, coordinate
    Real(8) :: gsgam2 (ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: gsgam2h(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: phi    (ADM_gall_in,ADM_kall,ADM_lall)
    Real(8) :: z      (ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: lat    (ADM_gall_in,ADM_lall)

    real(8) :: frhogq(ADM_gall_in,ADM_kall,ADM_lall)

    integer :: l, nq
    !---------------------------------------------------------------------------

    call GTL_clip_region(VMTR_GSGAM2 (:,:,:),gsgam2, 1,ADM_kall)
    call GTL_clip_region(VMTR_GSGAM2H(:,:,:),gsgam2h,1,ADM_kall)
    call GTL_clip_region(VMTR_PHI    (:,:,:),phi,    1,ADM_kall)
    call GTL_clip_region(GRD_vz(:,:,:,GRD_Z),z,      1,ADM_kall)

    call GTL_clip_region_1layer(GMTR_lat(:,:),lat)

    !--- get the prognostic and diagnostic variables
    call prgvar_get_in_withdiag( rhog,   & ! [IN]
                                 rhogvx, & ! [IN]
                                 rhogvy, & ! [IN]
                                 rhogvz, & ! [IN]
                                 rhogw,  & ! [IN]
                                 rhoge,  & ! [IN]
                                 rhogq,  & ! [IN]
                                 rho,    & ! [IN]
                                 pre,    & ! [IN]
                                 tem,    & ! [IN]
                                 vx,     & ! [IN]
                                 vy,     & ! [IN]
                                 vz,     & ! [IN]
                                 w,      & ! [IN]
                                 q       ) ! [IN]

    !--- boundary condition
    do l = 1, ADM_lall
       call bndcnd_thermo( ADM_gall_in, & ! [IN]
                           tem(:,:,l),  & ! [INOUT]
                           rho(:,:,l),  & ! [INOUT]
                           pre(:,:,l),  & ! [INOUT]
                           phi(:,:,l)   ) ! [IN]

       vx(:,ADM_kmax+1,l) = vx(:,ADM_kmax,l)
       vy(:,ADM_kmax+1,l) = vy(:,ADM_kmax,l)
       vz(:,ADM_kmax+1,l) = vz(:,ADM_kmax,l)
       vx(:,ADM_kmin-1,l) = vx(:,ADM_kmin,l)
       vy(:,ADM_kmin-1,l) = vy(:,ADM_kmin,l)
       vz(:,ADM_kmin-1,l) = vz(:,ADM_kmin,l)

       q(:,ADM_kmax+1,l,:) = 0.D0
       q(:,ADM_kmin-1,l,:) = 0.D0
    enddo

    ! forcing
    select case(trim(AF_TYPE))
    case('HELD-SUAREZ')

       do l = 1, ADM_lall
          call af_HeldSuarez( ADM_gall_in, & ! [IN]
                              lat(:,l),    & ! [IN]
                              pre(:,:,l),  & ! [IN]
                              tem(:,:,l),  & ! [IN]
                              vx (:,:,l),  & ! [IN]
                              vy (:,:,l),  & ! [IN]
                              vz (:,:,l),  & ! [IN]
                              fvx(:,:,l),  & ! [OUT]
                              fvy(:,:,l),  & ! [OUT]
                              fvz(:,:,l),  & ! [OUT]
                              fw (:,:,l),  & ! [OUT]
                              fe (:,:,l)   ) ! [OUT]

          call history_in( 'ml_af_fvx', fvx(:,:,l) )
          call history_in( 'ml_af_fvy', fvy(:,:,l) )
          call history_in( 'ml_af_fvz', fvz(:,:,l) )
          call history_in( 'ml_af_fw',  fw (:,:,l) )
          call history_in( 'ml_af_fe',  fe (:,:,l) )
       enddo
       fq(:,:,:,:) = 0.D0

    case default

       fvx(:,:,:) = 0.D0
       fvy(:,:,:) = 0.D0
       fvz(:,:,:) = 0.D0
       fw (:,:,:) = 0.D0
       fe (:,:,:) = 0.D0

       fq (:,:,:,:) = 0.D0

    end select

    rhogvx(:,:,:) = rhogvx(:,:,:) + TIME_DTL * fvx(:,:,:) * rho(:,:,:) * GSGAM2 (:,:,:)
    rhogvy(:,:,:) = rhogvy(:,:,:) + TIME_DTL * fvy(:,:,:) * rho(:,:,:) * GSGAM2 (:,:,:)
    rhogvz(:,:,:) = rhogvz(:,:,:) + TIME_DTL * fvz(:,:,:) * rho(:,:,:) * GSGAM2 (:,:,:)
    rhogw (:,:,:) = rhogw (:,:,:) + TIME_DTL * fw (:,:,:) * rho(:,:,:) * GSGAM2H(:,:,:)
    rhoge (:,:,:) = rhoge (:,:,:) + TIME_DTL * fe (:,:,:) * rho(:,:,:) * GSGAM2 (:,:,:)

    do nq = 1, TRC_VMAX
       frhogq(:,:,:) = fq(:,:,:,nq) * rho(:,:,:) * GSGAM2(:,:,:)

       rhog (:,:,:)    = rhog (:,:,:)    + TIME_DTL * frhogq(:,:,:)
       rhogq(:,:,:,nq) = rhogq(:,:,:,nq) + TIME_DTL * frhogq(:,:,:)
    enddo

    !--- set the prognostic variables
    call prgvar_set_in( rhog,   & ! [IN]
                        rhogvx, & ! [IN]
                        rhogvy, & ! [IN]
                        rhogvz, & ! [IN]
                        rhogw,  & ! [IN]
                        rhoge,  & ! [IN]
                        rhogq   ) ! [IN]

    return
  end subroutine forcing

end module mod_forcing_driver
!-------------------------------------------------------------------------------
