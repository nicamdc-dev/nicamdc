!-------------------------------------------------------------------------------
!>
!! 3D Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators using vertical metrics.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)    Imported from igdc-4.33
!! @li      2011-09-27 (T.Seiki)     merge optimization by RIST and M.Terai
!!
!<
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
    real(RP), intent(in)  :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_nxyz,ADM_gall,1:6        ,ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(ADM_nxyz,         1:ADM_vlink,ADM_lall_pl)

    real(RP) :: sclt   (ADM_gall   ,TI:TJ) ! scalar on the hexagon vertex
    real(RP) :: sclt_pl(ADM_gall_pl)
    real(RP) :: sclt_rhogw
    real(RP) :: sclt_rhogw_pl
    integer  :: gmin, gmax, gall, kmin, kmax, iall

    real(RP) :: rhogvx_vm   (ADM_gall   )          ! rho*vx / vertical metrics
    real(RP) :: rhogvx_vm_pl(ADM_gall_pl)
    real(RP) :: rhogvy_vm   (ADM_gall   )          ! rho*vy / vertical metrics
    real(RP) :: rhogvy_vm_pl(ADM_gall_pl)
    real(RP) :: rhogvz_vm   (ADM_gall   )          ! rho*vz / vertical metrics
    real(RP) :: rhogvz_vm_pl(ADM_gall_pl)
    real(RP) :: rhogw_vm    (ADM_gall,   ADM_kall) ! rho*w  / vertical metrics
    real(RP) :: rhogw_vm_pl (ADM_gall_pl,ADM_kall)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT3D_divdamp',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    gall = ADM_gall
    kmin = ADM_kmin
    kmax = ADM_kmax
    iall = ADM_gall_1d

    !$omp parallel workshare
    ddivdx(:,:,:) = 0.0_RP
    ddivdy(:,:,:) = 0.0_RP
    ddivdz(:,:,:) = 0.0_RP
    !$omp end parallel workshare

    do l = 1, ADM_lall

       !$omp parallel do default(none), private(ij,k), &
       !$omp shared(l,gall,kmin,kmax,rhogvx,rhogvy,rhogvz,rhogw,rhogw_vm,VMTR_C2WfactGz,VMTR_RGSQRTH,VMTR_RGAMH)
       do k = kmin+1, kmax
          do ij = 1, gall
             rhogw_vm(ij,k) = ( VMTR_C2WfactGz(ij,k,1,l) * rhogvx(ij,k  ,l) &
                              + VMTR_C2WfactGz(ij,k,2,l) * rhogvx(ij,k-1,l) &
                              + VMTR_C2WfactGz(ij,k,3,l) * rhogvy(ij,k  ,l) &
                              + VMTR_C2WfactGz(ij,k,4,l) * rhogvy(ij,k-1,l) &
                              + VMTR_C2WfactGz(ij,k,5,l) * rhogvz(ij,k  ,l) &
                              + VMTR_C2WfactGz(ij,k,6,l) * rhogvz(ij,k-1,l) &
                              ) * VMTR_RGAMH(ij,k,l)                        & ! horizontal contribution
                            + rhogw(ij,k,l) * VMTR_RGSQRTH(ij,k,l)            ! vertical   contribution
          enddo
       enddo

       !$omp parallel workshare
       rhogw_vm(:,kmin  ) = 0.0_RP
       rhogw_vm(:,kmax+1) = 0.0_RP
       !$omp end parallel workshare

       !$omp parallel default(none), private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1,sclt_rhogw), &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,gall,kmin,kmax,iall,ddivdx,ddivdy,ddivdz,rhogvx,rhogvy,rhogvz, &
       !$omp        rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm,sclt,coef_intp,coef_diff,GRD_rdgz,VMTR_RGAM)
       do k = kmin, kmax

          !$omp do
          do ij = 1, gall
             rhogvx_vm(ij) = rhogvx(ij,k,l) * VMTR_RGAM(ij,k,l)
             rhogvy_vm(ij) = rhogvy(ij,k,l) * VMTR_RGAM(ij,k,l)
             rhogvz_vm(ij) = rhogvz(ij,k,l) * VMTR_RGAM(ij,k,l)
          enddo
          !$omp end do

          !$omp do
          do j = gmin-1, gmax
             do i = gmin-1, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1

                sclt_rhogw = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ip1j,k+1) + rhogw_vm(ip1jp1,k+1) ) &
                             - ( rhogw_vm(ij,k  ) + rhogw_vm(ip1j,k  ) + rhogw_vm(ip1jp1,k  ) ) &
                             ) / 3.0_RP * GRD_rdgz(k)

                sclt(ij,TI) = coef_intp(XDIR,ij,1,TI,l) * rhogvx_vm(ij    ) &
                            + coef_intp(XDIR,ij,2,TI,l) * rhogvx_vm(ip1j  ) &
                            + coef_intp(XDIR,ij,3,TI,l) * rhogvx_vm(ip1jp1) &
                            + coef_intp(YDIR,ij,1,TI,l) * rhogvy_vm(ij    ) &
                            + coef_intp(YDIR,ij,2,TI,l) * rhogvy_vm(ip1j  ) &
                            + coef_intp(YDIR,ij,3,TI,l) * rhogvy_vm(ip1jp1) &
                            + coef_intp(ZDIR,ij,1,TI,l) * rhogvz_vm(ij    ) &
                            + coef_intp(ZDIR,ij,2,TI,l) * rhogvz_vm(ip1j  ) &
                            + coef_intp(ZDIR,ij,3,TI,l) * rhogvz_vm(ip1jp1) &
                            + sclt_rhogw
             enddo

             do i = gmin-1, gmax
                ij     = (j-1)*iall + i
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall

                sclt_rhogw = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ip1jp1,k+1) + rhogw_vm(ijp1,k+1) ) &
                             - ( rhogw_vm(ij,k  ) + rhogw_vm(ip1jp1,k  ) + rhogw_vm(ijp1,k  ) ) &
                             ) / 3.0_RP * GRD_rdgz(k)

                sclt(ij,TJ) = coef_intp(XDIR,ij,1,TJ,l) * rhogvx_vm(ij    ) &
                            + coef_intp(XDIR,ij,2,TJ,l) * rhogvx_vm(ip1jp1) &
                            + coef_intp(XDIR,ij,3,TJ,l) * rhogvx_vm(ijp1  ) &
                            + coef_intp(YDIR,ij,1,TJ,l) * rhogvy_vm(ij    ) &
                            + coef_intp(YDIR,ij,2,TJ,l) * rhogvy_vm(ip1jp1) &
                            + coef_intp(YDIR,ij,3,TJ,l) * rhogvy_vm(ijp1  ) &
                            + coef_intp(ZDIR,ij,1,TJ,l) * rhogvz_vm(ij    ) &
                            + coef_intp(ZDIR,ij,2,TJ,l) * rhogvz_vm(ip1jp1) &
                            + coef_intp(ZDIR,ij,3,TJ,l) * rhogvz_vm(ijp1  ) &
                            + sclt_rhogw
             enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
             j = gmin-1
             i = gmin-1

             ij     = (j-1)*iall + i
             ip1j   = ij + 1

             sclt(ij,TI) = sclt(ip1j,TJ)
             !$omp end master
          endif

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                ddivdx(ij,k,l) = ( coef_diff(XDIR,ij,1,l) * ( sclt(ij,    TI) + sclt(ij,    TJ) ) &
                                 + coef_diff(XDIR,ij,2,l) * ( sclt(ij,    TJ) + sclt(im1j,  TI) ) &
                                 + coef_diff(XDIR,ij,3,l) * ( sclt(im1j,  TI) + sclt(im1jm1,TJ) ) &
                                 + coef_diff(XDIR,ij,4,l) * ( sclt(im1jm1,TJ) + sclt(im1jm1,TI) ) &
                                 + coef_diff(XDIR,ij,5,l) * ( sclt(im1jm1,TI) + sclt(ijm1  ,TJ) ) &
                                 + coef_diff(XDIR,ij,6,l) * ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) )
             enddo

             do i = gmin, gmax
                ij     = (j-1)*iall + i
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                ddivdy(ij,k,l) = ( coef_diff(YDIR,ij,1,l) * ( sclt(ij,    TI) + sclt(ij,    TJ) ) &
                                 + coef_diff(YDIR,ij,2,l) * ( sclt(ij,    TJ) + sclt(im1j,  TI) ) &
                                 + coef_diff(YDIR,ij,3,l) * ( sclt(im1j,  TI) + sclt(im1jm1,TJ) ) &
                                 + coef_diff(YDIR,ij,4,l) * ( sclt(im1jm1,TJ) + sclt(im1jm1,TI) ) &
                                 + coef_diff(YDIR,ij,5,l) * ( sclt(im1jm1,TI) + sclt(ijm1  ,TJ) ) &
                                 + coef_diff(YDIR,ij,6,l) * ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) )
             enddo

             do i = gmin, gmax
                ij     = (j-1)*iall + i
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                ddivdz(ij,k,l) = ( coef_diff(ZDIR,ij,1,l) * ( sclt(ij,    TI) + sclt(ij,    TJ) ) &
                                 + coef_diff(ZDIR,ij,2,l) * ( sclt(ij,    TJ) + sclt(im1j,  TI) ) &
                                 + coef_diff(ZDIR,ij,3,l) * ( sclt(im1j,  TI) + sclt(im1jm1,TJ) ) &
                                 + coef_diff(ZDIR,ij,4,l) * ( sclt(im1jm1,TJ) + sclt(im1jm1,TI) ) &
                                 + coef_diff(ZDIR,ij,5,l) * ( sclt(im1jm1,TI) + sclt(ijm1  ,TJ) ) &
                                 + coef_diff(ZDIR,ij,6,l) * ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) )
             enddo
          enddo
          !$omp end do

       enddo
       !$omp end parallel

       !$omp parallel workshare
       ddivdx(:,ADM_kmin-1,l) = 0.0_RP
       ddivdx(:,ADM_kmax+1,l) = 0.0_RP
       ddivdy(:,ADM_kmin-1,l) = 0.0_RP
       ddivdy(:,ADM_kmax+1,l) = 0.0_RP
       ddivdz(:,ADM_kmin-1,l) = 0.0_RP
       ddivdz(:,ADM_kmax+1,l) = 0.0_RP
       !$omp end parallel workshare
    enddo ! l loop

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do ij = 1, ADM_gall_pl
             rhogw_vm_pl(ij,k) = ( VMTR_C2WfactGz_pl(ij,k,1,l) * rhogvx_pl(ij,k  ,l) &
                                 + VMTR_C2WfactGz_pl(ij,k,2,l) * rhogvx_pl(ij,k-1,l) &
                                 + VMTR_C2WfactGz_pl(ij,k,3,l) * rhogvy_pl(ij,k  ,l) &
                                 + VMTR_C2WfactGz_pl(ij,k,4,l) * rhogvy_pl(ij,k-1,l) &
                                 + VMTR_C2WfactGz_pl(ij,k,5,l) * rhogvz_pl(ij,k  ,l) &
                                 + VMTR_C2WfactGz_pl(ij,k,6,l) * rhogvz_pl(ij,k-1,l) &
                                 ) * VMTR_RGAMH_pl(ij,k,l)                           & ! horizontal contribution
                               + rhogw_pl(ij,k,l) * VMTR_RGSQRTH_pl(ij,k,l)            ! vertical   contribution
          enddo
          enddo
          do ij = 1, ADM_gall_pl
             rhogw_vm_pl(ij,ADM_kmin  ) = 0.0_RP
             rhogw_vm_pl(ij,ADM_kmax+1) = 0.0_RP
          enddo

          n = ADM_gslf_pl

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

                sclt_rhogw_pl = ( ( rhogw_vm_pl(n,k+1) + rhogw_vm_pl(ij,k+1) + rhogw_vm_pl(ijp1,k+1) ) &
                                - ( rhogw_vm_pl(n,k  ) + rhogw_vm_pl(ij,k  ) + rhogw_vm_pl(ijp1,k  ) ) &
                                ) / 3.0_RP * GRD_rdgz(k)

                sclt_pl(ij) = coef_intp_pl(XDIR,v,1,l) * rhogvx_vm_pl(n   ) &
                            + coef_intp_pl(XDIR,v,2,l) * rhogvx_vm_pl(ij  ) &
                            + coef_intp_pl(XDIR,v,3,l) * rhogvx_vm_pl(ijp1) &
                            + coef_intp_pl(YDIR,v,1,l) * rhogvy_vm_pl(n   ) &
                            + coef_intp_pl(YDIR,v,2,l) * rhogvy_vm_pl(ij  ) &
                            + coef_intp_pl(YDIR,v,3,l) * rhogvy_vm_pl(ijp1) &
                            + coef_intp_pl(ZDIR,v,1,l) * rhogvz_vm_pl(n   ) &
                            + coef_intp_pl(ZDIR,v,2,l) * rhogvz_vm_pl(ij  ) &
                            + coef_intp_pl(ZDIR,v,3,l) * rhogvz_vm_pl(ijp1) &
                            + sclt_rhogw_pl
             enddo

             ddivdx_pl(:,k,l) = 0.0_RP
             ddivdy_pl(:,k,l) = 0.0_RP
             ddivdz_pl(:,k,l) = 0.0_RP

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijm1 = v - 1
                if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

                ddivdx_pl(n,k,l) = ddivdx_pl(n,k,l) + coef_diff_pl(XDIR,v-1,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
                ddivdy_pl(n,k,l) = ddivdy_pl(n,k,l) + coef_diff_pl(YDIR,v-1,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
                ddivdz_pl(n,k,l) = ddivdz_pl(n,k,l) + coef_diff_pl(ZDIR,v-1,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
             enddo
          enddo

          do ij = 1, ADM_gall_pl
             ddivdx_pl(ij,ADM_kmin-1,l) = 0.0_RP
             ddivdx_pl(ij,ADM_kmax+1,l) = 0.0_RP
             ddivdy_pl(ij,ADM_kmin-1,l) = 0.0_RP
             ddivdy_pl(ij,ADM_kmax+1,l) = 0.0_RP
             ddivdz_pl(ij,ADM_kmin-1,l) = 0.0_RP
             ddivdz_pl(ij,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    else
       ddivdx_pl(:,:,:) = 0.0_RP
       ddivdy_pl(:,:,:) = 0.0_RP
       ddivdz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT3D_divdamp',2)

    return
  end subroutine OPRT3D_divdamp

end module mod_oprt3d
