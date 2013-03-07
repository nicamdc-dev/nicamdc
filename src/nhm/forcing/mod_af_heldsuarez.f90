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
    real(8), intent(in)  :: lat(ijdim)
    real(8), intent(in)  :: pre(ijdim,kdim)
    real(8), intent(in)  :: tem(ijdim,kdim)
    real(8), intent(in)  :: vx (ijdim,kdim)
    real(8), intent(in)  :: vy (ijdim,kdim)
    real(8), intent(in)  :: vz (ijdim,kdim)
    real(8), intent(out) :: fvx(ijdim,kdim)
    real(8), intent(out) :: fvy(ijdim,kdim)
    real(8), intent(out) :: fvz(ijdim,kdim)
    real(8), intent(out) :: fw (ijdim,kdim)
    real(8), intent(out) :: fe (ijdim,kdim)

    real(8) :: T_eq, acl, asl, ap0
    real(8), parameter :: T_eq2 = 200.D0
    real(8), parameter :: DT_y  =  60.D0 ! [K]
    real(8), parameter :: Dth_z =  10.D0 ! [K]

    real(8) :: k_t, k_v
    real(8) :: sigma, fact_sig
    real(8), parameter :: sigma_b = 0.7D0
    real(8), parameter :: k_f     = 1.D0 / ( 1.D0 * 86400.D0 )
    real(8), parameter :: k_a     = 1.D0 / (40.D0 * 86400.D0 )
    real(8), parameter :: k_s     = 1.D0 / ( 4.D0 * 86400.D0 )

    integer :: n, k
    !---------------------------------------------------------------------------

    fvx(:,:) = 0.D0
    fvy(:,:) = 0.D0
    fvz(:,:) = 0.D0
    fw (:,:) = 0.D0
    fe (:,:) = 0.D0

    do k = kmin, kmax
    do n = 1,    ijdim
       asl = abs( sin(lat(n)) )
       acl = abs( cos(lat(n)) )
       ap0 = abs( pre(n,k) / PRE00 )

       T_eq = (315.D0 - DT_y*asl*asl - Dth_z*log(ap0)*acl*acl ) * ap0**KAPPA
       T_eq = max(T_eq,T_eq2)

       sigma    = pre(n,k) / ( 0.5D0 * ( pre(n,kmin) + pre(n,kmin-1) ) )
       fact_sig = max( 0.D0, (sigma-sigma_b) / (1.D0-sigma_b) )

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
