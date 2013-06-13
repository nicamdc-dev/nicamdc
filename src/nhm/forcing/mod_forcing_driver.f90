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
  public :: updating
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
  integer, private, parameter :: I_RHOG     = 1 ! Density x G^{1/2} x gamma^2
  integer, private, parameter :: I_RHOGVX   = 2 ! Density x G^{1/2} x gamma^2 x Horizontal velocity (X-direction)
  integer, private, parameter :: I_RHOGVY   = 3 ! Density x G^{1/2} x gamma^2 x Horizontal velocity (Y-direction)
  integer, private, parameter :: I_RHOGVZ   = 4 ! Density x G^{1/2} x gamma^2 x Horizontal velocity (Z-direction)
  integer, private, parameter :: I_RHOGW    = 5 ! Density x G^{1/2} x gamma^2 x Vertical   velocity
  integer, private, parameter :: I_RHOGE    = 6 ! Density x G^{1/2} x gamma^2 x Internal Energy
  integer, private, parameter :: I_RHOGETOT = 7 ! Density x G^{1/2} x gamma^2 x Total Energy
  !
  integer, private, parameter :: nmax_TEND     = 7
  integer, private, parameter :: nmax_PROG     = 6
  integer, private, parameter :: nmax_v_mean_c = 5
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

  ! [add; original by H.Miura] 20130613 R.Yoshida
  !-----------------------------------------------------------------------------
  subroutine updating( &
       PROG0, PROG0_pl,  &  !--- IN : prognostic variables for save
       PROG,  PROG_pl    &  !--- INOUT : prognostic variables for update
       )
       !
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_kmax,    &
       ADM_kmin,    &
       ADM_log_fid, &  ! R.Yoshida 13/06/12 [add]
       ADM_proc_stop   ! R.Yoshida 13/06/12 [add]
    use mod_time, only:  &
       TIME_DTL
    use mod_grd, only: &
       GRD_x,    &
       GRD_x_pl, &
       GRD_vz,   &
       GRD_vz_pl
    use mod_gmtr, only: &
       GMTR_lon,    &
       GMTR_lon_pl, &
       GMTR_lat,    &
       GMTR_lat_pl
    use mod_runconf, only: &
       RUN_TYPE,       & ! R.Yoshida 13/06/13 [add]
       TRC_VMAX,       &
       TRC_ADV_TYPE
    use mod_af_trcadv, only: & ![add] 20130612 R.Yoshida
       test11_velocity,  &
       test12_velocity
    implicit none
    !--- prognostic variables (save)
    real(8), intent(in) :: PROG0     (ADM_gall,   ADM_kall,ADM_lall,   nmax_PROG)
    real(8), intent(in) :: PROG0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_PROG)
    !--- prognostic variables
    real(8), intent(inout) :: PROG      (ADM_gall,   ADM_kall,ADM_lall,   nmax_PROG)
    real(8), intent(inout) :: PROG_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,nmax_PROG)

    !--- horizontal velocity_x  ( physical )
    real(8) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_y  ( physical )
    real(8) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_z  ( physical )
    real(8) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- vertical velocity ( physical )
    real(8) :: w   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: w_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- density deviation from the base state ( G^{1/2} X gamma2 )
    real(8) :: rhogd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: ij, k ,l

    ! for tracer advection test  [add; original by H.Miura] 20130612 R.Yoshida
    real(8), save :: time=0.d0


    !--- reset density
    rhogd(:,:,:)    = PROG0(:,:,:,I_rhog)
    rhogd_pl(:,:,:) = PROG0_pl(:,:,:,I_rhog)

    !--- update velocity
    time=time+TIME_DTL
    vx=0.d0; vx_pl=0.d0
    vy=0.d0; vy_pl=0.d0
    vz=0.d0; vz_pl=0.d0
    w=0.d0

    select case (RUN_TYPE)
    !------------------------------------------------------------------------
    case ("TRCADV-1") 
       do l=1,ADM_lall
       do k=ADM_kmin-1,ADM_kmax+1
          ! full (1): u,v
          ! half (2): w
          do ij=1,ADM_gall
             call test11_velocity (time,GMTR_lon(ij,l),GMTR_lat(ij,l),GRD_vz(ij,k,l,1),GRD_vz(ij,k,l,2), &
                                   vx(ij,k,l),vy(ij,k,l),vz(ij,k,l),w(ij,k,l))
          end do
       end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
          do k=ADM_kmin-1,ADM_kmax+1
          do ij=1,ADM_GALL_PL
             call test11_velocity (time,GMTR_lon_pl(ij,l),GMTR_lat_pl(ij,l),GRD_vz_pl(ij,k,l,1),GRD_vz_pl(ij,k,l,2), &
                                   vx_pl(ij,k,l),vy_pl(ij,k,l),vz_pl(ij,k,l),w_pl(ij,k,l)) 
          end do
          end do
          end do
       end if

    !------------------------------------------------------------------------
    case ("TRCADV-2")
       do l=1,ADM_lall
       do k=ADM_kmin-1,ADM_kmax+1
          ! full (1): u,v
          ! half (2): w
          do ij=1,ADM_gall
             call test12_velocity (time,GMTR_lon(ij,l),GMTR_lat(ij,l),GRD_vz(ij,k,l,1),GRD_vz(ij,k,l,2), &
                                   vx(ij,k,l),vy(ij,k,l),vz(ij,k,l),w(ij,k,l))
          end do
       end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
          do k=ADM_kmin-1,ADM_kmax+1
          do ij=1,ADM_GALL_PL
             call test12_velocity (time,GMTR_lon_pl(ij,l),GMTR_lat_pl(ij,l),GRD_vz_pl(ij,k,l,1),GRD_vz_pl(ij,k,l,2), &
                                   vx_pl(ij,k,l),vy_pl(ij,k,l),vz_pl(ij,k,l),w_pl(ij,k,l)) 
          end do
          end do
          end do
       end if

    !------------------------------------------------------------------------
    end select
    !
    PROG(:,:,:,I_RHOGVX)=vx(:,:,:)*rhogd(:,:,:); PROG_pl(:,:,:,I_RHOGVX)=vx_pl(:,:,:)*rhogd_pl(:,:,:)
    PROG(:,:,:,I_RHOGVY)=vy(:,:,:)*rhogd(:,:,:); PROG_pl(:,:,:,I_RHOGVY)=vy_pl(:,:,:)*rhogd_pl(:,:,:)
    PROG(:,:,:,I_RHOGVZ)=vz(:,:,:)*rhogd(:,:,:); PROG_pl(:,:,:,I_RHOGVZ)=vz_pl(:,:,:)*rhogd_pl(:,:,:)
    PROG(:,:,:,I_RHOGW) =w(:,:,:) *rhogd(:,:,:); PROG_pl(:,:,:,I_RHOGW) =w_pl(:,:,:) *rhogd_pl(:,:,:)
    !
    return
  end subroutine updating
  !
end module mod_forcing_driver
!-------------------------------------------------------------------------------
