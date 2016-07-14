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
  use mod_stdio
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
  real(RP), private, parameter :: T_eq2 = 200.0_RP
  real(RP), private, parameter :: DT_y  =  60.0_RP ! [K]
  real(RP), private, parameter :: Dth_z =  10.0_RP ! [K]

  real(RP), private, parameter :: sigma_b = 0.7_RP
  real(RP), private, parameter :: k_f     = 1.0_RP / ( 1.0_RP * 86400.0_RP )
  real(RP), private, parameter :: k_a     = 1.0_RP / (40.0_RP * 86400.0_RP )
  real(RP), private, parameter :: k_s     = 1.0_RP / ( 4.0_RP * 86400.0_RP )

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_heldsuarez_init
    implicit none
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[af_heldsuarez]/Category[nhm forcing]'

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
       fe     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_const, only: &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: lat(ijdim)
    real(RP), intent(in)  :: pre(ijdim,kdim)
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(in)  :: vx (ijdim,kdim)
    real(RP), intent(in)  :: vy (ijdim,kdim)
    real(RP), intent(in)  :: vz (ijdim,kdim)
    real(RP), intent(out) :: fvx(ijdim,kdim)
    real(RP), intent(out) :: fvy(ijdim,kdim)
    real(RP), intent(out) :: fvz(ijdim,kdim)
    real(RP), intent(out) :: fe (ijdim,kdim)

    real(RP) :: T_eq, coslat, sinlat, ap0
    real(RP) :: sigma, factor

    integer  :: ij, k
    !---------------------------------------------------------------------------

    fvx(:,:) = 0.0_RP
    fvy(:,:) = 0.0_RP
    fvz(:,:) = 0.0_RP
    fe (:,:) = 0.0_RP

    do k  = kmin, kmax
    do ij = 1,    ijdim
       sigma  = pre(ij,k) / ( 0.5_RP * ( pre(ij,kmin) + pre(ij,kmin-1) ) )
       factor = max( (sigma-sigma_b) / (1.0_RP-sigma_b), 0.0_RP )

       fvx(ij,k) = -k_f * factor * vx(ij,k)
       fvy(ij,k) = -k_f * factor * vy(ij,k)
       fvz(ij,k) = -k_f * factor * vz(ij,k)

       sinlat   = abs( sin(lat(ij)) )
       coslat   = abs( cos(lat(ij)) )
       ap0      = abs( pre(ij,k) / PRE00 )

       T_eq     = ( 315.0_RP - DT_y*sinlat*sinlat - Dth_z*log(ap0)*coslat*coslat ) * ap0**(Rdry/CPdry)
       T_eq     = max( T_eq, T_eq2 )

       fe(ij,k) = -( k_a + (k_s-k_a) * factor * coslat**4 ) * (tem(ij,k)-T_eq) * CVdry
    enddo
    enddo

  end subroutine af_heldsuarez

end module mod_af_heldsuarez
!-------------------------------------------------------------------------------
