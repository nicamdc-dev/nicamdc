!-------------------------------------------------------------------------------
!> Module thermodyanics
!!
!! @par Description
!!          This module calculates thermodyanical variables
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_thrmdyn
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: THRMDYN_qd
  public :: THRMDYN_cv
  public :: THRMDYN_cp
  public :: THRMDYN_rho
  public :: THRMDYN_pre
  public :: THRMDYN_ein
  public :: THRMDYN_tem
  public :: THRMDYN_th
  public :: THRMDYN_eth
  public :: THRMDYN_ent
  public :: THRMDYN_rhoein
  public :: THRMDYN_tempre

  interface THRMDYN_qd
     module procedure THRMDYN_qd_ijk_SP
     module procedure THRMDYN_qd_ijk_DP
     module procedure THRMDYN_qd_ijkl
  end interface THRMDYN_qd

  interface THRMDYN_cv
     module procedure THRMDYN_cv_ijk_SP
     module procedure THRMDYN_cv_ijk_DP
  end interface THRMDYN_cv

  interface THRMDYN_cp
     module procedure THRMDYN_cp_ijk_SP
     module procedure THRMDYN_cp_ijk_DP
  end interface THRMDYN_cp

  interface THRMDYN_rho
     module procedure THRMDYN_rho_ijk_SP
     module procedure THRMDYN_rho_ijk_DP
  end interface THRMDYN_rho

  interface THRMDYN_pre
     module procedure THRMDYN_pre_ijk
  end interface THRMDYN_pre

  interface THRMDYN_ein
     module procedure THRMDYN_ein_ijk_SP
     module procedure THRMDYN_ein_ijk_DP
  end interface THRMDYN_ein

  interface THRMDYN_tem
     module procedure THRMDYN_tem_ijk
  end interface THRMDYN_tem

  interface THRMDYN_th
     module procedure THRMDYN_th_ijk_SP
     module procedure THRMDYN_th_ijk_DP
     module procedure THRMDYN_th_ijkl
  end interface THRMDYN_th

  interface THRMDYN_eth
     module procedure THRMDYN_eth_ijk
     module procedure THRMDYN_eth_ijkl
  end interface THRMDYN_eth

  interface THRMDYN_ent
     module procedure THRMDYN_ent_ijk_SP
     module procedure THRMDYN_ent_ijk_DP
  end interface THRMDYN_ent

  interface THRMDYN_rhoein
     module procedure THRMDYN_rhoein_ijk_SP
     module procedure THRMDYN_rhoein_ijk_DP
     module procedure THRMDYN_rhoein_ijkl_SP
     module procedure THRMDYN_rhoein_ijkl_DP
  end interface THRMDYN_rhoein

  interface THRMDYN_tempre
     module procedure THRMDYN_tempre_ijk_SP
     module procedure THRMDYN_tempre_ijk_DP
     module procedure THRMDYN_tempre_ijkl
  end interface THRMDYN_tempre

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijk_SP( &
       ijdim, &
       kdim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,qd,q)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_SP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(qd) pcopyin(q)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       qd(ij,k) = 1.0_SP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_qd_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijk_DP( &
       ijdim, &
       kdim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,qd,q)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_DP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(qd) pcopyin(q)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       qd(ij,k) = 1.0_DP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_qd_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: q (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: qd(ijdim,kdim,ldim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,qd,q)

!OCL XFILL
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k,l) = 1.0_RP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(qd) pcopyin(q)
!ACC!    do l  = 1, ldim
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       qd(ij,k,l) = 1.0_RP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR,NQW_END
!ACC!          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_qd_ijkl

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cv_ijk_SP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cv     )
    use mod_const, only: &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: cv(ijdim,kdim)       ! specific heat [J/kg/K]

    real(SP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(cv) pcopyin(qd,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = qd(ij,k) * CVdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_cv_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cv_ijk_DP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cv     )
    use mod_const, only: &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: cv(ijdim,kdim)       ! specific heat [J/kg/K]

    real(DP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(cv) pcopyin(qd,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = qd(ij,k) * CVdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_cv_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cp_ijk_SP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cp     )
    use mod_const, only: &
       CONST_CPdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CPW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: cp(ijdim,kdim)       ! specific heat [J/kg/K]

    real(SP) :: CPdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry = CONST_CPdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cp,qd,q,CPW,CPdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(cp) pcopyin(qd,q,CPW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cp(ij,k) = qd(ij,k) * CPdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_cp_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cp_ijk_DP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cp     )
    use mod_const, only: &
       CONST_CPdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CPW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: cp(ijdim,kdim)       ! specific heat [J/kg/K]

    real(DP) :: CPdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry = CONST_CPdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cp,qd,q,CPW,CPdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(cp) pcopyin(qd,q,CPW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cp(ij,k) = qd(ij,k) * CPdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_cp_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate density
  subroutine THRMDYN_rho_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       rho    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       I_QV
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(SP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: rho(ijdim,kdim)       ! density [kg/m3]

    real(SP) :: Rdry, Rvap

    integer  :: ij, k
    !---------------------------------------------------------------------------

    Rdry = CONST_Rdry
    Rvap = CONST_Rvap

!OCL XFILL
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,pre,tem,qd,q,Rdry,Rvap,I_QV)
    !$acc kernels pcopy(rho) pcopyin(pre,tem,qd,q)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
    enddo
    enddo
    !$acc end kernels
    !$omp end parallel do

    return
  end subroutine THRMDYN_rho_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate density
  subroutine THRMDYN_rho_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       rho    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       I_QV
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(DP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: rho(ijdim,kdim)       ! density [kg/m3]

    real(DP) :: Rdry, Rvap

    integer  :: ij, k
    !---------------------------------------------------------------------------

    Rdry = CONST_Rdry
    Rvap = CONST_Rvap

!OCL XFILL
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,pre,tem,qd,q,Rdry,Rvap,I_QV)
    !$acc kernels pcopy(rho) pcopyin(pre,tem,qd,q)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
    enddo
    enddo
    !$acc end kernels
    !$omp end parallel do

    return
  end subroutine THRMDYN_rho_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate pressure
  subroutine THRMDYN_pre_ijk( &
       ijdim, &
       kdim,  &
       rho,   &
       tem,   &
       qd,    &
       q,     &
       pre    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       I_QV
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(RP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(RP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: pre(ijdim,kdim)       ! pressure    [Pa]

    real(RP) :: Rdry, Rvap

    integer  :: ij, k
    !---------------------------------------------------------------------------

    Rdry = CONST_Rdry
    Rvap = CONST_Rvap

!OCL XFILL
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,pre,rho,tem,qd,q,Rdry,Rvap,I_QV)
    !$acc kernels pcopy(pre) pcopyin(rho,tem,qd,q)
    do k  = 1, kdim
    do ij = 1, ijdim
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$acc end kernels
    !$omp end parallel do

    return
  end subroutine THRMDYN_pre_ijk

  !-----------------------------------------------------------------------------
  !> calculate internal energy
  subroutine THRMDYN_ein_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       qd,    &
       q,     &
       ein    )
    use mod_const, only: &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(SP) :: cv(ijdim,kdim)
    real(SP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,ein,tem,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(ein,cv) pcopyin(tem,qd,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = qd(ij,k) * CVdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       ein(ij,k) = tem(ij,k) * cv(ij,k)
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_ein_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate internal energy
  subroutine THRMDYN_ein_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       qd,    &
       q,     &
       ein    )
    use mod_const, only: &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(DP) :: cv(ijdim,kdim)
    real(DP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,ein,tem,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(ein,cv) pcopyin(tem,qd,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = qd(ij,k) * CVdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       ein(ij,k) = tem(ij,k) * cv(ij,k)
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_ein_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate temperature
  subroutine THRMDYN_tem_ijk( &
       ijdim, &
       kdim,  &
       ein,   &
       qd,    &
       q,     &
       tem    )
    use mod_const, only: &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: ein(ijdim,kdim)       ! internal energy [J]
    real(RP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: tem(ijdim,kdim)       ! temperature [K]

    real(RP) :: cv(ijdim,kdim)
    real(RP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,tem,ein,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       tem(ij,k) = ein(ij,k) / cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(tem,cv) pcopyin(ein,qd,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = qd(ij,k) * CVdry
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       tem(ij,k) = ein(ij,k) / cv(ij,k)
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_tem_ijk

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       th     )
    use mod_const, only: &
       CONST_Rdry,  &
       CONST_CPdry, &
       CONST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim) ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim) ! pressure    [Pa]
    real(SP), intent(out) :: th (ijdim,kdim) ! potential temperature [K]

    real(SP) :: RovCP, PRE00

    integer  :: ij, k
    !---------------------------------------------------------------------------

    RovCP = CONST_Rdry / CONST_CPdry
    PRE00 = CONST_PRE00

    !$acc kernels pcopy(th) pcopyin(tem,pre)
!OCL XFILL
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,th,tem,pre,RovCP,PRE00)
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k) = tem(ij,k) * ( PRE00 / pre(ij,k) )**RovCP
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    return
  end subroutine THRMDYN_th_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       th     )
    use mod_const, only: &
       CONST_Rdry,  &
       CONST_CPdry, &
       CONST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim) ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim) ! pressure    [Pa]
    real(DP), intent(out) :: th (ijdim,kdim) ! potential temperature [K]

    real(DP) :: RovCP, PRE00

    integer  :: ij, k
    !---------------------------------------------------------------------------

    RovCP = CONST_Rdry / CONST_CPdry
    PRE00 = CONST_PRE00

    !$acc kernels pcopy(th) pcopyin(tem,pre)
!OCL XFILL
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,th,tem,pre,RovCP,PRE00)
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k) = tem(ij,k) * ( PRE00 / pre(ij,k) )**RovCP
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    return
  end subroutine THRMDYN_th_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       th     )
    use mod_const, only: &
       CONST_Rdry,  &
       CONST_CPdry, &
       CONST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: tem(ijdim,kdim,ldim) ! temperature [K]
    real(RP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure    [Pa]
    real(RP), intent(out) :: th (ijdim,kdim,ldim) ! potential temperature [K]

    real(RP) :: RovCP, PRE00

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    RovCP = CONST_Rdry / CONST_CPdry
    PRE00 = CONST_PRE00

!OCL XFILL
    !$acc kernels pcopy(th) pcopyin(tem,pre)
    !$omp parallel do default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,th,tem,pre,RovCP,PRE00), &
    !$omp collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k,l) = tem(ij,k,l) * ( PRE00 / pre(ij,k,l) )**RovCP
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    return
  end subroutine THRMDYN_th_ijkl

  !-----------------------------------------------------------------------------
  !> calculate enthalpy
  subroutine THRMDYN_eth_ijk( &
       ijdim, &
       kdim,  &
       ein,   &
       pre,   &
       rho,   &
       eth    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: ein(ijdim,kdim) ! internal energy [J]
    real(RP), intent(in)  :: pre(ijdim,kdim) ! pressure        [Pa]
    real(RP), intent(in)  :: rho(ijdim,kdim) ! density         [kg/m3]
    real(RP), intent(out) :: eth(ijdim,kdim) ! enthalpy

    integer  :: ij, k
    !---------------------------------------------------------------------------

!OCL XFILL
    !$acc kernels pcopy(eth) pcopyin(ein,pre,rho)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,eth,ein,pre,rho)
    do k  = 1, kdim
    do ij = 1, ijdim
       eth(ij,k) = ein(ij,k) + pre(ij,k) / rho(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    return
  end subroutine THRMDYN_eth_ijk

  !-----------------------------------------------------------------------------
  !> calculate enthalpy
  subroutine THRMDYN_eth_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       pre,   &
       rho,   &
       eth    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: ein(ijdim,kdim,ldim) ! internal energy [J]
    real(RP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure        [Pa]
    real(RP), intent(in)  :: rho(ijdim,kdim,ldim) ! density         [kg/m3]
    real(RP), intent(out) :: eth(ijdim,kdim,ldim) ! enthalpy

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

!OCL XFILL
    !$acc kernels pcopy(eth) pcopyin(ein,pre,rho)
    !$omp parallel do default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,eth,ein,pre,rho), &
    !$omp collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       eth(ij,k,l) = ein(ij,k,l) + pre(ij,k,l) / rho(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    return
  end subroutine THRMDYN_eth_ijkl

  !-----------------------------------------------------------------------------
  !> calculate entropy
  subroutine THRMDYN_ent_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       ent    )
    use mod_const, only: &
       CONST_Rdry,  &
       CONST_CPdry, &
       CONST_Rvap,  &
       CONST_LHV,   &
       CONST_LHF,   &
       CONST_PRE00, &
       CONST_TEM00, &
       CONST_PSAT0, &
       CONST_EPSvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)
    real(SP), intent(in)  :: pre(ijdim,kdim)
    real(SP), intent(in)  :: qd (ijdim,kdim)
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax)
    real(SP), intent(out) :: ent(ijdim,kdim)

    real(SP) :: Pdry
    real(SP) :: Pvap
    real(SP) :: LH(nqmax)
    real(SP) :: CPdry, Rdry, Rvap, TEM00, PRE00, PSAT0, EPSvap

    real(SP), parameter :: EPS = 1.E-5_SP

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry  = CONST_CPdry
    Rdry   = CONST_Rdry
    Rvap   = CONST_Rvap
    TEM00  = CONST_TEM00
    PRE00  = CONST_PRE00
    PSAT0  = CONST_PSAT0
    EPSvap = CONST_EPSvap

    do nq = NQW_STR, NQW_END
       if ( nq == I_QV ) then
          LH(nq) =  CONST_LHV / TEM00
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          LH(nq) = -CONST_LHF / TEM00
       else
          LH(nq) = 0.0_SP
       endif
    enddo

    !$omp parallel default(none),private(ij,k,Pdry,Pvap), &
    !$omp shared(ijdim,kdim,nq,NQW_STR,NQW_END,ent,tem,pre,qd,q,     &
    !$omp        LH,CVW,CPdry,Rdry,Rvap,TEM00,PRE00,PSAT0,EPSvap,I_QV)

    !$acc kernels pcopy(ent) pcopyin(tem,pre,qd,q)
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSvap*qd(ij,k) / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)    / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    do nq = NQW_STR, NQW_END
       !$acc kernels pcopy(ent) pcopyin(tem,q,CVW,LH)
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
       enddo
       !$omp end do
       !$acc end kernels
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(ent) pcopyin(tem,q,CVW,LH) async(0)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
!ACC!                                + q(ij,k,nq) * LH (nq) / TEM00
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_ent_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate entropy
  subroutine THRMDYN_ent_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       ent    )
    use mod_const, only: &
       CONST_Rdry,  &
       CONST_CPdry, &
       CONST_Rvap,  &
       CONST_LHV,   &
       CONST_LHF,   &
       CONST_PRE00, &
       CONST_TEM00, &
       CONST_PSAT0, &
       CONST_EPSvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)
    real(DP), intent(in)  :: pre(ijdim,kdim)
    real(DP), intent(in)  :: qd (ijdim,kdim)
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax)
    real(DP), intent(out) :: ent(ijdim,kdim)

    real(DP) :: Pdry
    real(DP) :: Pvap
    real(DP) :: LH(nqmax)
    real(DP) :: CPdry, Rdry, Rvap, TEM00, PRE00, PSAT0, EPSvap

    real(DP), parameter :: EPS = 1.E-10_DP

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry  = CONST_CPdry
    Rdry   = CONST_Rdry
    Rvap   = CONST_Rvap
    TEM00  = CONST_TEM00
    PRE00  = CONST_PRE00
    PSAT0  = CONST_PSAT0
    EPSvap = CONST_EPSvap

    do nq = NQW_STR, NQW_END
       if ( nq == I_QV ) then
          LH(nq) =  CONST_LHV / TEM00
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          LH(nq) = -CONST_LHF / TEM00
       else
          LH(nq) = 0.0_DP
       endif
    enddo

    !$omp parallel default(none),private(ij,k,Pdry,Pvap), &
    !$omp shared(ijdim,kdim,nq,NQW_STR,NQW_END,ent,tem,pre,qd,q,     &
    !$omp        LH,CVW,CPdry,Rdry,Rvap,TEM00,PRE00,PSAT0,EPSvap,I_QV)

    !$acc kernels pcopy(ent) pcopyin(tem,pre,qd,q)
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSvap*qd(ij,k) / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)    / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    do nq = NQW_STR, NQW_END
       !$acc kernels pcopy(ent) pcopyin(tem,q,CVW,LH)
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
       enddo
       !$omp end do
       !$acc end kernels
    enddo

    !$omp end parallel

!ACC!    !$acc kernels pcopy(ent) pcopyin(tem,q,CVW,LH) async(0)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
!ACC!                                + q(ij,k,nq) * LH (nq) / TEM00
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_ent_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(SP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(SP) :: cv(ijdim,kdim)
    real(SP) :: qd(ijdim,kdim)
    real(SP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_SP
       qd(ij,k) = 1.0_SP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry

       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = 0.0_SP
!ACC!       qd(ij,k) = 1.0_SP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry
!ACC!
!ACC!       rho(ij,k) = pre(ij,k) / tem(ij,k) / ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )
!ACC!       ein(ij,k) = tem(ij,k) * cv(ij,k)
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_rhoein_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(DP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(DP) :: cv(ijdim,kdim)
    real(DP) :: qd(ijdim,kdim)
    real(DP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_DP
       qd(ij,k) = 1.0_DP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry

       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = 0.0_DP
!ACC!       qd(ij,k) = 1.0_DP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry
!ACC!
!ACC!       rho(ij,k) = pre(ij,k) / tem(ij,k) / ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )
!ACC!       ein(ij,k) = tem(ij,k) * cv(ij,k)
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_rhoein_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijkl_SP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(SP), intent(in)  :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]
    real(SP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(SP), intent(out) :: ein(ijdim,kdim,ldim)       ! internal energy [J]

    real(SP) :: cv(ijdim,kdim,ldim)
    real(SP) :: qd(ijdim,kdim,ldim)
    real(SP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_SP
       qd(ij,k,l) = 1.0_SP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / ( ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l) )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW)
!ACC!    do l  = 1, ldim
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k,l) = 0.0_SP
!ACC!       qd(ij,k,l) = 1.0_SP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
!ACC!          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry
!ACC!
!ACC!       rho(ij,k,l) = pre(ij,k,l) / tem(ij,k,l) / ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap )
!ACC!       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
!ACC!    enddo
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_rhoein_ijkl_SP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijkl_DP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(DP), intent(in)  :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]
    real(DP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(DP), intent(out) :: ein(ijdim,kdim,ldim)       ! internal energy [J]

    real(DP) :: cv(ijdim,kdim,ldim)
    real(DP) :: qd(ijdim,kdim,ldim)
    real(DP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_DP
       qd(ij,k,l) = 1.0_DP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / ( ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l) )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW)
!ACC!    do l  = 1, ldim
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k,l) = 0.0_DP
!ACC!       qd(ij,k,l) = 1.0_DP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
!ACC!          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry
!ACC!
!ACC!       rho(ij,k,l) = pre(ij,k,l) / tem(ij,k,l) / ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap )
!ACC!       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
!ACC!    enddo
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_rhoein_ijkl_DP

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijk_SP( &
       ijdim, &
       kdim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: ein(ijdim,kdim)       ! internal energy [J]
    real(SP), intent(in)  :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(out) :: pre(ijdim,kdim)       ! pressure    [Pa]

    real(SP) :: cv(ijdim,kdim)
    real(SP) :: qd(ijdim,kdim)
    real(SP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,tem,pre,ein,rho,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_SP
       qd(ij,k) = 1.0_SP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv (ij,k) = cv (ij,k) + qd (ij,k) * CVdry
       tem(ij,k) = ein(ij,k) / cv (ij,k)
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = 0.0_SP
!ACC!       qd(ij,k) = 1.0_SP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry
!ACC!
!ACC!       tem(ij,k) = ein(ij,k) / cv(ij,k)
!ACC!       pre(ij,k) = rho(ij,k) * tem(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_tempre_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijk_DP( &
       ijdim, &
       kdim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: ein(ijdim,kdim)       ! internal energy [J]
    real(DP), intent(in)  :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(out) :: pre(ijdim,kdim)       ! pressure    [Pa]

    real(DP) :: cv(ijdim,kdim)
    real(DP) :: qd(ijdim,kdim)
    real(DP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,tem,pre,ein,rho,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_DP
       qd(ij,k) = 1.0_DP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv (ij,k) = cv (ij,k) + qd (ij,k) * CVdry
       tem(ij,k) = ein(ij,k) / cv (ij,k)
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW)
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k) = 0.0_DP
!ACC!       qd(ij,k) = 1.0_DP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
!ACC!          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry
!ACC!
!ACC!       tem(ij,k) = ein(ij,k) / cv(ij,k)
!ACC!       pre(ij,k) = rho(ij,k) * tem(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_tempre_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    use mod_const, only: &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: ein(ijdim,kdim,ldim)       ! internal energy [J]
    real(RP), intent(in)  :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(RP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(RP), intent(out) :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]

    real(RP) :: cv(ijdim,kdim,ldim)
    real(RP) :: qd(ijdim,kdim,ldim)
    real(RP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    CVdry = CONST_CVdry
    Rdry  = CONST_Rdry
    Rvap  = CONST_Rvap

    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,tem,pre,ein,rho,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_RP
       qd(ij,k,l) = 1.0_RP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv (ij,k,l) = cv (ij,k,l) + qd (ij,k,l) * CVdry
       tem(ij,k,l) = ein(ij,k,l) / cv (ij,k,l)
       pre(ij,k,l) = rho(ij,k,l) * ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

!ACC!    !$acc kernels pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW)
!ACC!    do l  = 1, ldim
!ACC!    do k  = 1, kdim
!ACC!    do ij = 1, ijdim
!ACC!       cv(ij,k,l) = 0.0_RP
!ACC!       qd(ij,k,l) = 1.0_RP
!ACC!
!ACC!       !$acc loop seq
!ACC!       do nq = NQW_STR, NQW_END
!ACC!          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
!ACC!          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
!ACC!       enddo
!ACC!       !$acc end loop
!ACC!
!ACC!       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry
!ACC!
!ACC!       tem(ij,k,l) = ein(ij,k,l) / cv(ij,k,l)
!ACC!       pre(ij,k,l) = rho(ij,k,l) * tem(ij,k,l) * ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap )
!ACC!    enddo
!ACC!    enddo
!ACC!    enddo
!ACC!    !$acc end kernels

    return
  end subroutine THRMDYN_tempre_ijkl

end module mod_thrmdyn
