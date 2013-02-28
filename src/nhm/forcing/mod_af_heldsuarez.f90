!-------------------------------------------------------------------------------
!
!+  HeldSuarez forcing module
!
!-------------------------------------------------------------------------------
module mod_af_heldsuarez
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for forcing term
  !       of GCSS CASE1.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Add this module.
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: af_heldsuarez_init
  public :: af_heldsuarez
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_heldsuarez_init
    !
    return
    !
  end subroutine af_heldsuarez_init
  !-----------------------------------------------------------------------------
  subroutine af_heldsuarez(     &
       ijdim,                   &
       lat,                     &  !--- IN : latitude  
       rho,                     &  !--- IN : density
       vx,                      &  !--- IN : Vx
       vy,                      &  !--- IN : Vy
       vz,                      &  !--- IN : Vz
       w,                       &  !--- IN : w
       pre,                     &  !--- IN : pressure
       tem,                     &  !--- IN : temperature
       fvx,                     &  !--- INOUT : tend. of rhovx
       fvy,                     &  !--- INOUT : tend. of rhovy
       fvz,                     &  !--- INOUT : tend. of rhovz
       fw,                      &  !--- INOUT : tend. of rhow
       fe                       &  !--- INOUT : tend. of rhoe
       )
    use mod_adm, only :  &
         kdim => ADM_kall,&
         kmin => ADM_kmin,&
         kmax => ADM_kmax
    use mod_cnst, only : &
         CNST_PRE00,     &
         CNST_KAPPA,     &
         CNST_CP
    implicit none
    !
    integer, intent(in) :: ijdim 
    real(8), intent(in) :: lat(1:ijdim)
    real(8), intent(in) :: rho(1:ijdim,1:kdim)
    real(8), intent(in) :: vx(1:ijdim,1:kdim)
    real(8), intent(in) :: vy(1:ijdim,1:kdim)
    real(8), intent(in) :: vz(1:ijdim,1:kdim)
    real(8), intent(in) :: w(1:ijdim,1:kdim)
    real(8), intent(in) :: pre(1:ijdim,1:kdim)
    real(8), intent(in) :: tem(1:ijdim,1:kdim)
    !
    real(8), intent(out) :: fvx(1:ijdim,1:kdim)
    real(8), intent(out) :: fvy(1:ijdim,1:kdim)
    real(8), intent(out) :: fvz(1:ijdim,1:kdim)
    real(8), intent(out) :: fw(1:ijdim,1:kdim)
    real(8), intent(out) :: fe(1:ijdim,1:kdim)

    real(8) :: T_eq(1:ijdim,1:kdim)

    real(8),parameter :: sigma_b = 0.7D0
    real(8),parameter :: k_f = 1.0D0/(1.0D0*3600.0D0*24.0D0)
    real(8),parameter :: k_a = 1.0D0/(40.0D0*3600.0D0*24.0D0)
    real(8),parameter :: k_s = 1.0D0/(4.0D0*3600.0D0*24.0D0)
    real(8),parameter :: DT_y = 60.0D0   ! K
    real(8),parameter :: Dth_z = 10.0D0  ! K


    real(8) :: sigma
    real(8) :: k_t,k_v

    integer :: l,k,n

    real(8),parameter :: T_eq2=200.0d0
    real(8) :: T_eq1
    real(8) :: acl,asl,ap0
    !
    do k=1,kdim
       do n=1,ijdim
          asl=abs(sin(lat(n)))
          acl=abs(cos(lat(n)))
          ap0=abs(pre(n,k)/CNST_PRE00)
          T_eq1=(315.0d0-DT_y*asl*asl-Dth_z*log(ap0)*acl*acl)*ap0**CNST_KAPPA
          T_eq(n,k)=max(T_eq1,T_eq2)
       enddo
    enddo
    !
    do k = kmin, kmax
       do n=1, ijdim
          sigma = pre(n,k)/( 0.5D0*( pre(n,kmin) + pre(n,kmin-1) ) )
          k_v = k_f * max(0.0D0,(sigma-sigma_b)/(1.0D0-sigma_b))
          fvx(n,k) = - ( k_v * ( vx(n,k)  ) )
          fvy(n,k) = - ( k_v * ( vy(n,k)  ) )
          fvz(n,k) = - ( k_v * ( vz(n,k)  ) )
          k_t = K_a + ( k_s - k_a ) * max(0.0D0,(sigma-sigma_b)/(1.0D0-sigma_b))&
               * dabs(cos(lat(n)))**4.0D0
          fe(n,k) = - ( k_t * CNST_CP * ( tem(n,k) - T_eq(n,k) ) )
       end do
    end do
    fw(:,:) = 0.0D0
    !
  end subroutine af_heldsuarez

  !-----------------------------------------------------------------------------
end module mod_af_heldsuarez
!-------------------------------------------------------------------------------
