!-------------------------------------------------------------------------------
!>
!! Module Held-Suarez forcing
!!
!! @par Description
!!         This module contains subroutines for forcing term of GCSS CASE1.
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  [NEW]
!!
!<
module mod_af_heldsuarez
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
  public :: af_heldsuarez_init
  public :: af_heldsuarez

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
  subroutine af_heldsuarez_init
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine af_heldsuarez_init

  !-----------------------------------------------------------------------------
  subroutine af_heldsuarez( &
       ijdim, &
       lat,   &
       pre,   &
       tem,   &
       vx,    &
       vy,    &
       vz,    &
       fvx,   &
       fvy,   &
       fvz,   &
       fw,    &
       fe     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CV    => CNST_CV,    &
       KAPPA => CNST_KAPPA, &
       PRE00 => CNST_PRE00
    implicit none

    integer, intent(in)  :: ijdim
    REAL(RP), intent(in)  :: lat(ijdim)
    REAL(RP), intent(in)  :: pre(ijdim,kdim)
    REAL(RP), intent(in)  :: tem(ijdim,kdim)
    REAL(RP), intent(in)  :: vx (ijdim,kdim)
    REAL(RP), intent(in)  :: vy (ijdim,kdim)
    REAL(RP), intent(in)  :: vz (ijdim,kdim)
    REAL(RP), intent(out) :: fvx(ijdim,kdim)
    REAL(RP), intent(out) :: fvy(ijdim,kdim)
    REAL(RP), intent(out) :: fvz(ijdim,kdim)
    REAL(RP), intent(out) :: fw (ijdim,kdim)
    REAL(RP), intent(out) :: fe (ijdim,kdim)

    REAL(RP) :: T_eq, acl, asl, ap0
    REAL(RP), parameter :: T_eq2 = 200.0_RP
    REAL(RP), parameter :: DT_y  =  60.0_RP ! [K]
    REAL(RP), parameter :: Dth_z =  10.0_RP ! [K]

    REAL(RP) :: k_t, k_v
    REAL(RP) :: sigma, fact_sig
    REAL(RP), parameter :: sigma_b = 0.7_RP
    REAL(RP), parameter :: k_f     = 1.0_RP / ( 1.0_RP * 86400.0_RP )
    REAL(RP), parameter :: k_a     = 1.0_RP / (40.0_RP * 86400.0_RP )
    REAL(RP), parameter :: k_s     = 1.0_RP / ( 4.0_RP * 86400.0_RP )

    integer :: n, k
    !---------------------------------------------------------------------------

    fvx(:,:) = 0.0_RP
    fvy(:,:) = 0.0_RP
    fvz(:,:) = 0.0_RP
    fw (:,:) = 0.0_RP
    fe (:,:) = 0.0_RP

    do k = kmin, kmax
    do n = 1,    ijdim
       asl = abs( sin(lat(n)) )
       acl = abs( cos(lat(n)) )
       ap0 = abs( pre(n,k) / PRE00 )

       T_eq = (315.0_RP - DT_y*asl*asl - Dth_z*log(ap0)*acl*acl ) * ap0**KAPPA
       T_eq = max(T_eq,T_eq2)

       sigma    = pre(n,k) / ( 0.5_RP * ( pre(n,kmin) + pre(n,kmin-1) ) )
       fact_sig = max( 0.0_RP, (sigma-sigma_b) / (1.0_RP-sigma_b) )

       k_v = k_f * fact_sig

       fvx(n,k) = -k_v * vx(n,k)
       fvy(n,k) = -k_v * vy(n,k)
       fvz(n,k) = -k_v * vz(n,k)

       k_t = k_a + ( k_s-k_a ) * fact_sig * acl**4

       fe(n,k) = -k_t * CV * ( tem(n,k)-T_eq )
    enddo
    enddo

  end subroutine af_heldsuarez

end module mod_af_heldsuarez
!-------------------------------------------------------------------------------
