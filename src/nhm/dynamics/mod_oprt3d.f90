!-------------------------------------------------------------------------------
!> Module 3D Operator
!!
!! @par Description
!!          This module contains the subroutines for differential operators using vertical metrics.
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_oprt3d
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
  !++ Public procedures
  !
  public :: OPRT3D_divdamp

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
  !> 3D divergence damping operator
  subroutine OPRT3D_divdamp( &
       ddivdx,    ddivdx_pl,    &
       ddivdy,    ddivdy_pl,    &
       ddivdz,    ddivdz_pl,    &
       rhogvx,    rhogvx_pl,    &
       rhogvy,    rhogvy_pl,    &
       rhogvz,    rhogvz_pl,    &
       rhogw,     rhogw_pl,     &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       TI  => ADM_TI,  &
       TJ  => ADM_TJ,  &
       ADM_nxyz,       &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_vlink,      &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_1d,    &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,    &
       ADM_kmin,       &
       ADM_kmax
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR, &
       GRD_rdgz
    use mod_vmtr, only: &
       VMTR_RGSQRTH,     &
       VMTR_RGSQRTH_pl,  &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_C2WfactGz,   &
       VMTR_C2WfactGz_pl
    implicit none

    real(RP), intent(out) :: ddivdx      (ADM_gall   ,ADM_kall,ADM_lall   ) ! tendency
    real(RP), intent(out) :: ddivdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vx { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vy { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vz { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w  { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_gall   ,1:3,ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_gall_pl,1:3,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_gall,1:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(         1:ADM_vlink,ADM_nxyz,ADM_lall_pl)

#ifdef _OPENACC
    real(RP) :: sclt   (ADM_gall,ADM_kall,TI:TJ) ! scalar on the hexagon vertex
#else
    real(RP) :: sclt   (ADM_gall,TI:TJ)          ! scalar on the hexagon vertex
#endif
    real(RP) :: sclt_pl(ADM_gall_pl)
    real(RP) :: sclt_rhogw
    real(RP) :: sclt_rhogw_pl

#ifdef _OPENACC
    real(RP) :: rhogvx_vm   (ADM_gall,ADM_kall)                ! rho*vx / vertical metrics
    real(RP) :: rhogvy_vm   (ADM_gall,ADM_kall)                ! rho*vy / vertical metrics
    real(RP) :: rhogvz_vm   (ADM_gall,ADM_kall)                ! rho*vz / vertical metrics
#else
    real(RP) :: rhogvx_vm   (ADM_gall)                         ! rho*vx / vertical metrics
    real(RP) :: rhogvy_vm   (ADM_gall)                         ! rho*vy / vertical metrics
    real(RP) :: rhogvz_vm   (ADM_gall)                         ! rho*vz / vertical metrics
#endif
    real(RP) :: rhogvx_vm_pl(ADM_gall_pl)
    real(RP) :: rhogvy_vm_pl(ADM_gall_pl)
    real(RP) :: rhogvz_vm_pl(ADM_gall_pl)
    real(RP) :: rhogw_vm    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  / vertical metrics
    real(RP) :: rhogw_vm_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer  :: gmin, gmax, iall, gall, kall, kmin, kmax, lall, gminm1

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: g, k, l, n, v
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(sclt,rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm) &
    !$acc pcopy(ddivdx,ddivdy,ddivdz) &
    !$acc pcopyin(rhogvx,rhogvy,rhogvz,rhogw,coef_intp,coef_diff) &
    !$acc pcopyin(ADM_have_sgp,GRD_rdgz,VMTR_RGSQRTH,VMTR_RGAM,VMTR_RGAMH,VMTR_C2WfactGz)

    call PROF_rapstart('OPRT3D_divdamp',2)

    gmin = (ADM_gmin-1)*ADM_gall_1d + ADM_gmin
    gmax = (ADM_gmax-1)*ADM_gall_1d + ADM_gmax
    iall = ADM_gall_1d
    gall = ADM_gall
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    gminm1 = (ADM_gmin-1-1)*ADM_gall_1d + ADM_gmin-1

    !$acc kernels present(rhogw_vm) pcopyin(rhogvx,rhogvy,rhogvz,rhogw,VMTR_RGSQRTH,VMTR_RGAMH,VMTR_C2WfactGz)
    !$omp parallel default(none),private(g,k,l), &
    !$omp shared(gall,kmin,kmax,lall,rhogw_vm,rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2WfactGz,VMTR_RGSQRTH,VMTR_RGAMH)
    do l = 1, lall
       !$omp do
       do k = kmin+1, kmax
       do g = 1, gall
          rhogw_vm(g,k,l) = ( VMTR_C2WfactGz(g,k,1,l) * rhogvx(g,k  ,l) &
                            + VMTR_C2WfactGz(g,k,2,l) * rhogvx(g,k-1,l) &
                            + VMTR_C2WfactGz(g,k,3,l) * rhogvy(g,k  ,l) &
                            + VMTR_C2WfactGz(g,k,4,l) * rhogvy(g,k-1,l) &
                            + VMTR_C2WfactGz(g,k,5,l) * rhogvz(g,k  ,l) &
                            + VMTR_C2WfactGz(g,k,6,l) * rhogvz(g,k-1,l) &
                            ) * VMTR_RGAMH(g,k,l)                       & ! horizontal contribution
                          + rhogw(g,k,l) * VMTR_RGSQRTH(g,k,l)            ! vertical   contribution
       enddo
       enddo
       !$omp end do nowait

!OCL XFILL
       !$omp do
       do g = 1, gall
          rhogw_vm(g,kmin  ,l) = 0.0_RP
          rhogw_vm(g,kmax+1,l) = 0.0_RP
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !$omp parallel default(none),private(g,k,l,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1,sclt_rhogw), &
    !$omp shared(ADM_have_sgp,gminm1,gmin,gmax,gall,kmin,kmax,lall,iall,ddivdx,ddivdy,ddivdz,rhogvx,rhogvy,rhogvz, &
    !$omp        rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm,sclt,coef_intp,coef_diff,GRD_rdgz,VMTR_RGAM)
    do l = 1, lall
       !$acc kernels pcopy(ddivdx,ddivdy,ddivdz) present(sclt,rhogvx_vm,rhogvy_vm,rhogvz_vm) &
       !$acc pcopyin(rhogvx,rhogvy,rhogvz,rhogw_vm,coef_intp,coef_diff,ADM_have_sgp,GRD_rdgz,VMTR_RGAM)
       do k = kmin, kmax
!OCL XFILL
          !$omp do
          do g = 1, gall
#ifdef _OPENACC
             rhogvx_vm(g,k) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
             rhogvy_vm(g,k) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
             rhogvz_vm(g,k) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
#else
             rhogvx_vm(g) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
             rhogvy_vm(g) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
             rhogvz_vm(g) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
#endif
          enddo
          !$omp end do

          !$omp do
          do g = gminm1, gmax
             ij     = g
             ip1j   = g + 1
             ip1jp1 = g + iall + 1
             ijp1   = g + iall

             sclt_rhogw = ( ( rhogw_vm(ij,k+1,l) + rhogw_vm(ip1j,k+1,l) + rhogw_vm(ip1jp1,k+1,l) ) &
                          - ( rhogw_vm(ij,k  ,l) + rhogw_vm(ip1j,k  ,l) + rhogw_vm(ip1jp1,k  ,l) ) &
                          ) / 3.0_RP * GRD_rdgz(k)

#ifdef _OPENACC
             sclt(g,k,TI) = coef_intp(g,1,XDIR,TI,l) * rhogvx_vm(ij    ,k) &
                          + coef_intp(g,2,XDIR,TI,l) * rhogvx_vm(ip1j  ,k) &
                          + coef_intp(g,3,XDIR,TI,l) * rhogvx_vm(ip1jp1,k) &
                          + coef_intp(g,1,YDIR,TI,l) * rhogvy_vm(ij    ,k) &
                          + coef_intp(g,2,YDIR,TI,l) * rhogvy_vm(ip1j  ,k) &
                          + coef_intp(g,3,YDIR,TI,l) * rhogvy_vm(ip1jp1,k) &
                          + coef_intp(g,1,ZDIR,TI,l) * rhogvz_vm(ij    ,k) &
                          + coef_intp(g,2,ZDIR,TI,l) * rhogvz_vm(ip1j  ,k) &
                          + coef_intp(g,3,ZDIR,TI,l) * rhogvz_vm(ip1jp1,k) &
                          + sclt_rhogw
#else
             sclt(g,TI)   = coef_intp(g,1,XDIR,TI,l) * rhogvx_vm(ij    ) &
                          + coef_intp(g,2,XDIR,TI,l) * rhogvx_vm(ip1j  ) &
                          + coef_intp(g,3,XDIR,TI,l) * rhogvx_vm(ip1jp1) &
                          + coef_intp(g,1,YDIR,TI,l) * rhogvy_vm(ij    ) &
                          + coef_intp(g,2,YDIR,TI,l) * rhogvy_vm(ip1j  ) &
                          + coef_intp(g,3,YDIR,TI,l) * rhogvy_vm(ip1jp1) &
                          + coef_intp(g,1,ZDIR,TI,l) * rhogvz_vm(ij    ) &
                          + coef_intp(g,2,ZDIR,TI,l) * rhogvz_vm(ip1j  ) &
                          + coef_intp(g,3,ZDIR,TI,l) * rhogvz_vm(ip1jp1) &
                          + sclt_rhogw
#endif
          enddo
          !$omp end do nowait

          !$omp do
          do g = gminm1, gmax
             ij     = g
             ip1j   = g + 1
             ip1jp1 = g + iall + 1
             ijp1   = g + iall

             sclt_rhogw = ( ( rhogw_vm(ij,k+1,l) + rhogw_vm(ip1jp1,k+1,l) + rhogw_vm(ijp1,k+1,l) ) &
                          - ( rhogw_vm(ij,k  ,l) + rhogw_vm(ip1jp1,k  ,l) + rhogw_vm(ijp1,k  ,l) ) &
                          ) / 3.0_RP * GRD_rdgz(k)

#ifdef _OPENACC
             sclt(g,k,TJ) = coef_intp(g,1,XDIR,TJ,l) * rhogvx_vm(ij    ,k) &
                          + coef_intp(g,2,XDIR,TJ,l) * rhogvx_vm(ip1jp1,k) &
                          + coef_intp(g,3,XDIR,TJ,l) * rhogvx_vm(ijp1  ,k) &
                          + coef_intp(g,1,YDIR,TJ,l) * rhogvy_vm(ij    ,k) &
                          + coef_intp(g,2,YDIR,TJ,l) * rhogvy_vm(ip1jp1,k) &
                          + coef_intp(g,3,YDIR,TJ,l) * rhogvy_vm(ijp1  ,k) &
                          + coef_intp(g,1,ZDIR,TJ,l) * rhogvz_vm(ij    ,k) &
                          + coef_intp(g,2,ZDIR,TJ,l) * rhogvz_vm(ip1jp1,k) &
                          + coef_intp(g,3,ZDIR,TJ,l) * rhogvz_vm(ijp1  ,k) &
                          + sclt_rhogw
#else
             sclt(g,TJ)   = coef_intp(g,1,XDIR,TJ,l) * rhogvx_vm(ij    ) &
                          + coef_intp(g,2,XDIR,TJ,l) * rhogvx_vm(ip1jp1) &
                          + coef_intp(g,3,XDIR,TJ,l) * rhogvx_vm(ijp1  ) &
                          + coef_intp(g,1,YDIR,TJ,l) * rhogvy_vm(ij    ) &
                          + coef_intp(g,2,YDIR,TJ,l) * rhogvy_vm(ip1jp1) &
                          + coef_intp(g,3,YDIR,TJ,l) * rhogvy_vm(ijp1  ) &
                          + coef_intp(g,1,ZDIR,TJ,l) * rhogvz_vm(ij    ) &
                          + coef_intp(g,2,ZDIR,TJ,l) * rhogvz_vm(ip1jp1) &
                          + coef_intp(g,3,ZDIR,TJ,l) * rhogvz_vm(ijp1  ) &
                          + sclt_rhogw
#endif
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
#ifdef _OPENACC
             sclt(gminm1,k,TI) = sclt(gminm1+1,k,TJ)
#else
             sclt(gminm1,TI) = sclt(gminm1+1,TJ)
#endif
             !$omp end master
          endif
!OCL XFILL
          !$omp do
          do g = 1, gmin-1
             ddivdx(g,k,l) = 0.0_RP
             ddivdy(g,k,l) = 0.0_RP
             ddivdz(g,k,l) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do g = gmin, gmax
             ij     = g
             im1j   = g - 1
             im1jm1 = g - iall - 1
             ijm1   = g - iall

#ifdef _OPENACC
             ddivdx(g,k,l) = ( coef_diff(g,1,XDIR,l) * ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) &
                             + coef_diff(g,2,XDIR,l) * ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) &
                             + coef_diff(g,3,XDIR,l) * ( sclt(im1j,  k,TI) + sclt(im1jm1,k,TJ) ) &
                             + coef_diff(g,4,XDIR,l) * ( sclt(im1jm1,k,TJ) + sclt(im1jm1,k,TI) ) &
                             + coef_diff(g,5,XDIR,l) * ( sclt(im1jm1,k,TI) + sclt(ijm1  ,k,TJ) ) &
                             + coef_diff(g,6,XDIR,l) * ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) )
#else
             ddivdx(g,k,l) = ( coef_diff(g,1,XDIR,l) * ( sclt(ij,    TI) + sclt(ij,    TJ) ) &
                             + coef_diff(g,2,XDIR,l) * ( sclt(ij,    TJ) + sclt(im1j,  TI) ) &
                             + coef_diff(g,3,XDIR,l) * ( sclt(im1j,  TI) + sclt(im1jm1,TJ) ) &
                             + coef_diff(g,4,XDIR,l) * ( sclt(im1jm1,TJ) + sclt(im1jm1,TI) ) &
                             + coef_diff(g,5,XDIR,l) * ( sclt(im1jm1,TI) + sclt(ijm1  ,TJ) ) &
                             + coef_diff(g,6,XDIR,l) * ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) )
#endif
          enddo
          !$omp end do nowait

          !$omp do
          do g = gmin, gmax
             ij     = g
             im1j   = g - 1
             im1jm1 = g - iall - 1
             ijm1   = g - iall

#ifdef _OPENACC
             ddivdy(g,k,l) = ( coef_diff(g,1,YDIR,l) * ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) &
                             + coef_diff(g,2,YDIR,l) * ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) &
                             + coef_diff(g,3,YDIR,l) * ( sclt(im1j,  k,TI) + sclt(im1jm1,k,TJ) ) &
                             + coef_diff(g,4,YDIR,l) * ( sclt(im1jm1,k,TJ) + sclt(im1jm1,k,TI) ) &
                             + coef_diff(g,5,YDIR,l) * ( sclt(im1jm1,k,TI) + sclt(ijm1  ,k,TJ) ) &
                             + coef_diff(g,6,YDIR,l) * ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) )
#else
             ddivdy(g,k,l) = ( coef_diff(g,1,YDIR,l) * ( sclt(ij,    TI) + sclt(ij,    TJ) ) &
                             + coef_diff(g,2,YDIR,l) * ( sclt(ij,    TJ) + sclt(im1j,  TI) ) &
                             + coef_diff(g,3,YDIR,l) * ( sclt(im1j,  TI) + sclt(im1jm1,TJ) ) &
                             + coef_diff(g,4,YDIR,l) * ( sclt(im1jm1,TJ) + sclt(im1jm1,TI) ) &
                             + coef_diff(g,5,YDIR,l) * ( sclt(im1jm1,TI) + sclt(ijm1  ,TJ) ) &
                             + coef_diff(g,6,YDIR,l) * ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) )
#endif
          enddo
          !$omp end do nowait

          !$omp do
          do g = gmin, gmax
             ij     = g
             im1j   = g - 1
             im1jm1 = g - iall - 1
             ijm1   = g - iall

#ifdef _OPENACC
             ddivdz(g,k,l) = ( coef_diff(g,1,ZDIR,l) * ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) &
                             + coef_diff(g,2,ZDIR,l) * ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) &
                             + coef_diff(g,3,ZDIR,l) * ( sclt(im1j,  k,TI) + sclt(im1jm1,k,TJ) ) &
                             + coef_diff(g,4,ZDIR,l) * ( sclt(im1jm1,k,TJ) + sclt(im1jm1,k,TI) ) &
                             + coef_diff(g,5,ZDIR,l) * ( sclt(im1jm1,k,TI) + sclt(ijm1  ,k,TJ) ) &
                             + coef_diff(g,6,ZDIR,l) * ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) )
#else
             ddivdz(g,k,l) = ( coef_diff(g,1,ZDIR,l) * ( sclt(ij,    TI) + sclt(ij,    TJ) ) &
                             + coef_diff(g,2,ZDIR,l) * ( sclt(ij,    TJ) + sclt(im1j,  TI) ) &
                             + coef_diff(g,3,ZDIR,l) * ( sclt(im1j,  TI) + sclt(im1jm1,TJ) ) &
                             + coef_diff(g,4,ZDIR,l) * ( sclt(im1jm1,TJ) + sclt(im1jm1,TI) ) &
                             + coef_diff(g,5,ZDIR,l) * ( sclt(im1jm1,TI) + sclt(ijm1  ,TJ) ) &
                             + coef_diff(g,6,ZDIR,l) * ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) )
#endif
          enddo
          !$omp end do nowait
!OCL XFILL
          !$omp do
          do g = gmax+1, gall
             ddivdx(g,k,l) = 0.0_RP
             ddivdy(g,k,l) = 0.0_RP
             ddivdz(g,k,l) = 0.0_RP
          enddo
          !$omp end do
       enddo ! loop k
       !$acc end kernels
!OCL XFILL
       !$omp do
       do g = 1, gall
          ddivdx(g,kmin-1,l) = 0.0_RP
          ddivdy(g,kmin-1,l) = 0.0_RP
          ddivdz(g,kmin-1,l) = 0.0_RP
          ddivdx(g,kmax+1,l) = 0.0_RP
          ddivdy(g,kmax+1,l) = 0.0_RP
          ddivdz(g,kmax+1,l) = 0.0_RP
       enddo
       !$omp end do
    enddo ! loop l
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,k,l) = ( VMTR_C2WfactGz_pl(g,k,1,l) * rhogvx_pl(g,k  ,l) &
                                  + VMTR_C2WfactGz_pl(g,k,2,l) * rhogvx_pl(g,k-1,l) &
                                  + VMTR_C2WfactGz_pl(g,k,3,l) * rhogvy_pl(g,k  ,l) &
                                  + VMTR_C2WfactGz_pl(g,k,4,l) * rhogvy_pl(g,k-1,l) &
                                  + VMTR_C2WfactGz_pl(g,k,5,l) * rhogvz_pl(g,k  ,l) &
                                  + VMTR_C2WfactGz_pl(g,k,6,l) * rhogvz_pl(g,k-1,l) &
                                  ) * VMTR_RGAMH_pl(g,k,l)                          & ! horizontal contribution
                                + rhogw_pl(g,k,l) * VMTR_RGSQRTH_pl(g,k,l)            ! vertical   contribution
          enddo
          enddo

          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,ADM_kmin  ,l) = 0.0_RP
             rhogw_vm_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo

       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do v = 1, ADM_gall_pl
                rhogvx_vm_pl(v) = rhogvx_pl(v,k,l) * VMTR_RGAM_pl(v,k,l)
                rhogvy_vm_pl(v) = rhogvy_pl(v,k,l) * VMTR_RGAM_pl(v,k,l)
                rhogvz_vm_pl(v) = rhogvz_pl(v,k,l) * VMTR_RGAM_pl(v,k,l)
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

                sclt_rhogw_pl = ( ( rhogw_vm_pl(n,k+1,l) + rhogw_vm_pl(ij,k+1,l) + rhogw_vm_pl(ijp1,k+1,l) ) &
                                - ( rhogw_vm_pl(n,k  ,l) + rhogw_vm_pl(ij,k  ,l) + rhogw_vm_pl(ijp1,k  ,l) ) &
                                ) / 3.0_RP * GRD_rdgz(k)

                sclt_pl(ij) = coef_intp_pl(v,1,XDIR,l) * rhogvx_vm_pl(n   ) &
                            + coef_intp_pl(v,2,XDIR,l) * rhogvx_vm_pl(ij  ) &
                            + coef_intp_pl(v,3,XDIR,l) * rhogvx_vm_pl(ijp1) &
                            + coef_intp_pl(v,1,YDIR,l) * rhogvy_vm_pl(n   ) &
                            + coef_intp_pl(v,2,YDIR,l) * rhogvy_vm_pl(ij  ) &
                            + coef_intp_pl(v,3,YDIR,l) * rhogvy_vm_pl(ijp1) &
                            + coef_intp_pl(v,1,ZDIR,l) * rhogvz_vm_pl(n   ) &
                            + coef_intp_pl(v,2,ZDIR,l) * rhogvz_vm_pl(ij  ) &
                            + coef_intp_pl(v,3,ZDIR,l) * rhogvz_vm_pl(ijp1) &
                            + sclt_rhogw_pl
             enddo

             ddivdx_pl(:,k,l) = 0.0_RP
             ddivdy_pl(:,k,l) = 0.0_RP
             ddivdz_pl(:,k,l) = 0.0_RP

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijm1 = v - 1
                if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

                ddivdx_pl(n,k,l) = ddivdx_pl(n,k,l) + coef_diff_pl(v-1,XDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
                ddivdy_pl(n,k,l) = ddivdy_pl(n,k,l) + coef_diff_pl(v-1,YDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
                ddivdz_pl(n,k,l) = ddivdz_pl(n,k,l) + coef_diff_pl(v-1,ZDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
             enddo
          enddo

          ddivdx_pl(:,ADM_kmin-1,l) = 0.0_RP
          ddivdx_pl(:,ADM_kmax+1,l) = 0.0_RP
          ddivdy_pl(:,ADM_kmin-1,l) = 0.0_RP
          ddivdy_pl(:,ADM_kmax+1,l) = 0.0_RP
          ddivdz_pl(:,ADM_kmin-1,l) = 0.0_RP
          ddivdz_pl(:,ADM_kmax+1,l) = 0.0_RP
       enddo
    else
       ddivdx_pl(:,:,:) = 0.0_RP
       ddivdy_pl(:,:,:) = 0.0_RP
       ddivdz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT3D_divdamp',2)

    !$acc end data

    return
  end subroutine OPRT3D_divdamp

end module mod_oprt3d
