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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OPRT3D_divdamp
  public :: OPRT3D_gradient_intpl

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
       grdx,     grdx_pl,     &
       grdy,     grdy_pl,     &
       grdz,     grdz_pl,     &
       rhovx_in, rhovx_in_pl, &
       rhovy_in, rhovy_in_pl, &
       rhovz_in, rhovz_in_pl, &
       rhow,     rhow_pl      )
    use mod_adm, only: &
       ADM_prc_me,   &
       ADM_PRC_PL,   &
       ADM_rgn_vnum, &
       ADM_prc_tab,  &
       ADM_W,        &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_AI,       &
       ADM_AIJ,      &
       ADM_AJ,       &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_GSLF_PL,  &
       ADM_GMIN_PL,  &
       ADM_GMAX_PL,  &
       ADM_kmin,     &
       ADM_kmax,     &
       ADM_KNONE
    use mod_grd, only: &
       GRD_dgz,  &
       GRD_afac, &
       GRD_bfac
    use mod_gmtr, only: &
       GMTR_T_RAREA,  &
       GMTR_A_HNX,    &
       GMTR_A_HNY,    &
       GMTR_A_HNZ,    &
       GMTR_A_TNX,    &
       GMTR_A_TNY,    &
       GMTR_A_TNZ,    &
       GMTR_A_TN2X,   &
       GMTR_A_TN2Y,   &
       GMTR_A_TN2Z,   &
       GMTR_P_RAREA,  &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    use mod_vmtr, only: &
       VMTR_RGAM,       &
       VMTR_RGAM_pl,    &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_GSGAMH,     &
       VMTR_GSGAMH_pl,  &
       VMTR_GZXH,       &
       VMTR_GZXH_pl,    &
       VMTR_GZYH,       &
       VMTR_GZYH_pl,    &
       VMTR_GZZH,       &
       VMTR_GZZH_pl,    &
       VMTR_RGSH,       &
       VMTR_RGSH_pl
    implicit none

    real(8), intent(out) :: grdx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grdy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grdz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhovx_in   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhovx_in_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhovy_in   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhovy_in_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhovz_in   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhovz_in_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhow       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhow_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: flx_vz    (ADM_gall,   ADM_kall)
    real(8) :: flx_vz_pl (ADM_gall_pl,ADM_kall)
    real(8) :: sclt      (ADM_gall,   ADM_kall,ADM_TI:ADM_TJ)
    real(8) :: sclt_pl   (ADM_gall_pl,ADM_kall)

    real(8) :: rhovx_k   (ADM_gall   )
    real(8) :: rhovx_k_pl(ADM_gall_pl)
    real(8) :: rhovy_k   (ADM_gall   )
    real(8) :: rhovy_k_pl(ADM_gall_pl)
    real(8) :: rhovz_k   (ADM_gall   )
    real(8) :: rhovz_k_pl(ADM_gall_pl)

    real(8) :: ux1, ux2, ux3
    real(8) :: uy1, uy2, uy3
    real(8) :: uz1, uz2, uz3
    real(8) :: flx_vzt_top, flx_vzt_bot
    real(8) :: dp1, dp2, dp3, dp4, dp5, dp6
    real(8) :: dp_pl

    integer :: nstart, nend
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1
    integer :: ij_pl, ijp1_pl, ijm1_pl

    integer :: k0
    integer :: rgnid
    integer :: ij, k, l, n

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    grdx   (:,1:ADM_kmin-1,:) = 0.D0
    grdx_pl(:,1:ADM_kmin-1,:) = 0.D0
    grdy   (:,1:ADM_kmin-1,:) = 0.D0
    grdy_pl(:,1:ADM_kmin-1,:) = 0.D0
    grdz   (:,1:ADM_kmin-1,:) = 0.D0
    grdz_pl(:,1:ADM_kmin-1,:) = 0.D0

    grdx   (:,ADM_kmax+1:ADM_kall,:) = 0.D0
    grdx_pl(:,ADM_kmax+1:ADM_kall,:) = 0.D0
    grdy   (:,ADM_kmax+1:ADM_kall,:) = 0.D0
    grdy_pl(:,ADM_kmax+1:ADM_kall,:) = 0.D0
    grdz   (:,ADM_kmax+1:ADM_kall,:) = 0.D0
    grdz_pl(:,ADM_kmax+1:ADM_kall,:) = 0.D0

    k0 = ADM_KNONE

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do k = ADM_kmin, ADM_kmax+1
          do ij = 1, ADM_gall
             flx_vz(ij,k) = ( ( GRD_afac(k) * VMTR_RGSGAM2(ij,k,  l) * rhovx_in(ij,k,  l) &
                              + GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l) * rhovx_in(ij,k-1,l) &
                              ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZXH(ij,k,l)         &
                            + ( GRD_afac(k) * VMTR_RGSGAM2(ij,k,  l) * rhovy_in(ij,k,  l) &
                              + GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l) * rhovy_in(ij,k-1,l) &
                              ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZYH(ij,k,l)         &
                            + ( GRD_afac(k) * VMTR_RGSGAM2(ij,k,  l) * rhovz_in(ij,k,  l) &
                              + GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l) * rhovz_in(ij,k-1,l) &
                              ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZZH(ij,k,l)         &
                            ) + rhow(ij,k,l) * VMTR_RGSH(ij,k,l)
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax
          do ij = 1, ADM_gall
             rhovx_k(ij) = rhovx_in(ij,k,l) * VMTR_RGAM(ij,k,l)
             rhovy_k(ij) = rhovy_in(ij,k,l) * VMTR_RGAM(ij,k,l)
             rhovz_k(ij) = rhovz_in(ij,k,l) * VMTR_RGAM(ij,k,l)
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax,  ADM_gmax  )

          do n = nstart, nend
             ij     = n
             ip1j   = ij+1
             ip1jp1 = ij+1 + ADM_gall_1d

             ux1 = -( rhovx_k(ij)     + rhovx_k(ip1j)   )
             uy1 = -( rhovy_k(ij)     + rhovy_k(ip1j)   )
             uz1 = -( rhovz_k(ij)     + rhovz_k(ip1j)   )

             ux2 = -( rhovx_k(ip1j)   + rhovx_k(ip1jp1) )
             uy2 = -( rhovy_k(ip1j)   + rhovy_k(ip1jp1) )
             uz2 = -( rhovz_k(ip1j)   + rhovz_k(ip1jp1) )

             ux3 = +( rhovx_k(ip1jp1) + rhovx_k(ij)     )
             uy3 = +( rhovy_k(ip1jp1) + rhovy_k(ij)     )
             uz3 = +( rhovz_k(ip1jp1) + rhovz_k(ij)     )

             flx_vzt_top = ( flx_vz(ij,    k+1) &
                           + flx_vz(ip1j,  k+1) &
                           + flx_vz(ip1jp1,k+1) ) / 3.D0
             flx_vzt_bot = ( flx_vz(ij,    k  ) &
                           + flx_vz(ip1j,  k  ) &
                           + flx_vz(ip1jp1,k  ) ) / 3.D0

             sclt(ij,k,ADM_TI) = ( ux1 * GMTR_A_var(ij,  k0,l,ADM_AI, GMTR_A_TNX)    &
                                 + uy1 * GMTR_A_var(ij,  k0,l,ADM_AI, GMTR_A_TNY)    &
                                 + uz1 * GMTR_A_var(ij,  k0,l,ADM_AI, GMTR_A_TNZ)    &
                                 + ux2 * GMTR_A_var(ip1j,k0,l,ADM_AJ, GMTR_A_TNX)    &
                                 + uy2 * GMTR_A_var(ip1j,k0,l,ADM_AJ, GMTR_A_TNY)    &
                                 + uz2 * GMTR_A_var(ip1j,k0,l,ADM_AJ, GMTR_A_TNZ)    &
                                 + ux3 * GMTR_A_var(ij,  k0,l,ADM_AIJ,GMTR_A_TNX)    &
                                 + uy3 * GMTR_A_var(ij,  k0,l,ADM_AIJ,GMTR_A_TNY)    &
                                 + uz3 * GMTR_A_var(ij,  k0,l,ADM_AIJ,GMTR_A_TNZ)    &
                                 ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TI,GMTR_T_RAREA) &
                               + ( flx_vzt_top - flx_vzt_bot ) / GRD_dgz(k)

          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax,  ADM_gmax  )

          do n = nstart, nend
             ij     = n
             ijp1   = ij   + ADM_gall_1d
             ip1jp1 = ij+1 + ADM_gall_1d

             ux1 = -( rhovx_k(ij)     + rhovx_k(ip1jp1) )
             uy1 = -( rhovy_k(ij)     + rhovy_k(ip1jp1) )
             uz1 = -( rhovz_k(ij)     + rhovz_k(ip1jp1) )

             ux2 = +( rhovx_k(ip1jp1) + rhovx_k(ijp1)   )
             uy2 = +( rhovy_k(ip1jp1) + rhovy_k(ijp1)   )
             uz2 = +( rhovz_k(ip1jp1) + rhovz_k(ijp1)   )

             ux3 = +( rhovx_k(ijp1)   + rhovx_k(ij)     )
             uy3 = +( rhovy_k(ijp1)   + rhovy_k(ij)     )
             uz3 = +( rhovz_k(ijp1)   + rhovz_k(ij)     )

             flx_vzt_top = ( flx_vz(ij,    k+1) &
                           + flx_vz(ijp1,  k+1) &
                           + flx_vz(ip1jp1,k+1) ) / 3.D0
             flx_vzt_bot = ( flx_vz(ij,    k  ) &
                           + flx_vz(ijp1,  k  ) &
                           + flx_vz(ip1jp1,k  ) ) / 3.D0

             sclt(ij,k,ADM_TJ) = ( ux1 * GMTR_A_var(ij,  k0,l,ADM_AIJ,GMTR_A_TNX)    &
                                 + uy1 * GMTR_A_var(ij,  k0,l,ADM_AIJ,GMTR_A_TNY)    &
                                 + uz1 * GMTR_A_var(ij,  k0,l,ADM_AIJ,GMTR_A_TNZ)    &
                                 + ux2 * GMTR_A_var(ijp1,k0,l,ADM_AI, GMTR_A_TNX)    &
                                 + uy2 * GMTR_A_var(ijp1,k0,l,ADM_AI, GMTR_A_TNY)    &
                                 + uz2 * GMTR_A_var(ijp1,k0,l,ADM_AI, GMTR_A_TNZ)    &
                                 + ux3 * GMTR_A_var(ij,  k0,l,ADM_AJ, GMTR_A_TNX)    &
                                 + uy3 * GMTR_A_var(ij,  k0,l,ADM_AJ, GMTR_A_TNY)    &
                                 + uz3 * GMTR_A_var(ij,  k0,l,ADM_AJ, GMTR_A_TNZ)    &
                                 ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TJ,GMTR_T_RAREA) &
                               + ( flx_vzt_top - flx_vzt_bot ) / GRD_dgz(k)
          enddo

          nstart = suf(ADM_gmin,ADM_gmin)
          nend   = suf(ADM_gmax,ADM_gmax)

          do n = nstart, nend
             ij     = n
             im1j   = ij-1
             im1jm1 = ij-1 - ADM_gall_1d
             ijm1   = ij   - ADM_gall_1d

             dp1 = 0.5D0 * ( sclt(ijm1,  k,ADM_TJ) + sclt(ij,    k,ADM_TI) )
             dp2 = 0.5D0 * ( sclt(ij,    k,ADM_TI) + sclt(ij,    k,ADM_TJ) )
             dp3 = 0.5D0 * ( sclt(ij,    k,ADM_TJ) + sclt(im1j,  k,ADM_TI) )
             dp4 = 0.5D0 * ( sclt(im1jm1,k,ADM_TJ) + sclt(im1j,  k,ADM_TI) )
             dp5 = 0.5D0 * ( sclt(im1jm1,k,ADM_TI) + sclt(im1jm1,k,ADM_TJ) )
             dp6 = 0.5D0 * ( sclt(ijm1  ,k,ADM_TJ) + sclt(im1jm1,k,ADM_TI) )

             grdx(ij,k,l) = ( dp1 * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNX) &
                            + dp2 * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNX) &
                            + dp3 * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNX) &
                            - dp4 * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNX) &
                            - dp5 * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNX) &
                            - dp6 * GMTR_A_var(ijm1,  k0,l,ADM_AJ, GMTR_A_HNX) &
                            ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)

             grdy(ij,k,l) = ( dp1 * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNY) &
                            + dp2 * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNY) &
                            + dp3 * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNY) &
                            - dp4 * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNY) &
                            - dp5 * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNY) &
                            - dp6 * GMTR_A_var(ijm1,  k0,l,ADM_AJ, GMTR_A_HNY) &
                            ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)     

             grdz(ij,k,l) = ( dp1 * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNZ) &
                            + dp2 * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNZ) &
                            + dp3 * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNZ) &
                            - dp4 * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNZ) &
                            - dp5 * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNZ) &
                            - dp6 * GMTR_A_var(ijm1,  k0,l,ADM_AJ, GMTR_A_HNZ) &
                            ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)   
          enddo

          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
             n = suf(ADM_gmin,ADM_gmin)

             ij     = n
             im1j   = ij-1
             im1jm1 = ij-1 - ADM_gall_1d
             ijm1   = ij   - ADM_gall_1d

             sclt(im1jm1,k,ADM_TI) = sclt(ijm1,k,ADM_TJ) ! copy

             dp1 = 0.5D0 * ( sclt(ijm1,  k,ADM_TJ) + sclt(ij,    k,ADM_TI) )
             dp2 = 0.5D0 * ( sclt(ij,    k,ADM_TI) + sclt(ij,    k,ADM_TJ) )
             dp3 = 0.5D0 * ( sclt(ij,    k,ADM_TJ) + sclt(im1j,  k,ADM_TI) )
             dp4 = 0.5D0 * ( sclt(im1jm1,k,ADM_TJ) + sclt(im1j,  k,ADM_TI) )
             dp5 = 0.5D0 * ( sclt(im1jm1,k,ADM_TI) + sclt(im1jm1,k,ADM_TJ) )

             grdx(ij,k,l) = ( dp1 * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNX) &
                            + dp2 * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNX) &
                            + dp3 * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNX) &
                            - dp4 * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNX) &
                            - dp5 * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNX) &
                            ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)

             grdy(ij,k,l) = ( dp1 * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNY) &
                            + dp2 * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNY) &
                            + dp3 * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNY) &
                            - dp4 * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNY) &
                            - dp5 * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNY) &
                            ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)

             grdz(ij,k,l) = ( dp1 * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNZ) &
                            + dp2 * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNZ) &
                            + dp3 * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNZ) &
                            - dp4 * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNZ) &
                            - dp5 * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNZ) &
                            ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)
          endif

       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,k) = ( ( GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k,  l) * rhovx_in_pl(ij,k,  l) &
                                    + GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l) * rhovx_in_pl(ij,k-1,l) &
                                    ) *  0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZXH_pl(ij,k,l)        &
                                  + ( GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k,  l) * rhovy_in_pl(ij,k,  l) &
                                    + GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l) * rhovy_in_pl(ij,k-1,l) &
                                    ) *  0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZYH_pl(ij,k,l)        & 
                                  + ( GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k,  l) * rhovz_in_pl(ij,k,  l) &
                                    + GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l) * rhovz_in_pl(ij,k-1,l) &
                                    ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZZH_pl(ij,k,l)         &
                                  ) + rhow_pl(ij,k,l) * VMTR_RGSH_pl(ij,k,l)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax

             do ij = 1, ADM_gall_pl
                rhovx_k_pl(ij) = rhovx_in_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                rhovy_k_pl(ij) = rhovy_in_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                rhovz_k_pl(ij) = rhovz_in_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
             enddo

             ij = ADM_GSLF_PL
             do n = ADM_GMIN_PL, ADM_GMAX_PL
                ij_pl   = n
                ijp1_pl = n+1

                if( ij_pl == ADM_GMAX_PL ) ijp1_pl = ADM_GMIN_PL ! cyclic condition

                ux1 = +( rhovx_k_pl(ij     ) + rhovx_k_pl(ij_pl  ) )
                uy1 = +( rhovy_k_pl(ij     ) + rhovy_k_pl(ij_pl  ) )
                uz1 = +( rhovz_k_pl(ij     ) + rhovz_k_pl(ij_pl  ) )

                ux2 = +( rhovx_k_pl(ij_pl  ) + rhovx_k_pl(ijp1_pl) )
                uy2 = +( rhovy_k_pl(ij_pl  ) + rhovy_k_pl(ijp1_pl) )
                uz2 = +( rhovz_k_pl(ij_pl  ) + rhovz_k_pl(ijp1_pl) )

                ux3 = -( rhovx_k_pl(ijp1_pl) + rhovx_k_pl(ij     ) )
                uy3 = -( rhovy_k_pl(ijp1_pl) + rhovy_k_pl(ij     ) )
                uz3 = -( rhovz_k_pl(ijp1_pl) + rhovz_k_pl(ij     ) )

                flx_vzt_top = ( flx_vz_pl(ij,     k+1) &
                              + flx_vz_pl(ij_pl,  k+1) &
                              + flx_vz_pl(ijp1_pl,k+1) ) / 3.D0
                flx_vzt_bot = ( flx_vz_pl(ij,     k  ) &
                              + flx_vz_pl(ij_pl,  k  ) &
                              + flx_vz_pl(ijp1_pl,k  ) ) / 3.D0

                sclt_pl(ij_pl,k) = ( ux1 * GMTR_A_var_pl(ij_pl,  k0,l,GMTR_A_TNX )    &
                                   + uy1 * GMTR_A_var_pl(ij_pl,  k0,l,GMTR_A_TNY )    &
                                   + uz1 * GMTR_A_var_pl(ij_pl,  k0,l,GMTR_A_TNZ )    &
                                   + ux2 * GMTR_A_var_pl(ij_pl,  k0,l,GMTR_A_TN2X)    &
                                   + uy2 * GMTR_A_var_pl(ij_pl,  k0,l,GMTR_A_TN2Y)    &
                                   + uz2 * GMTR_A_var_pl(ij_pl,  k0,l,GMTR_A_TN2Z)    &
                                   + ux3 * GMTR_A_var_pl(ijp1_pl,k0,l,GMTR_A_TNX )    &
                                   + uy3 * GMTR_A_var_pl(ijp1_pl,k0,l,GMTR_A_TNY )    &
                                   + uz3 * GMTR_A_var_pl(ijp1_pl,k0,l,GMTR_A_TNZ )    &
                                   ) * 0.5D0 * GMTR_T_var_pl(ij_pl,k0,l,GMTR_T_RAREA) &
                                 + ( flx_vzt_top - flx_vzt_bot ) / GRD_dgz(k)
             enddo

             grdx_pl(ij,k,l) = 0.D0
             grdy_pl(ij,k,l) = 0.D0
             grdz_pl(ij,k,l) = 0.D0
             do n = ADM_GMIN_PL, ADM_GMAX_PL
                ij_pl   = n
                ijm1_pl = n-1

                if( ij_pl == ADM_GMIN_PL ) ijm1_pl = ADM_GMAX_PL ! cyclic condition

                dp_pl = 0.5D0 * ( sclt_pl(ijm1_pl,k) + sclt_pl(ij_pl,k) )

                grdx_pl(ij,k,l) = grdx_pl(ij,k,l) &
                                + dp_pl * GMTR_A_var_pl(ij_pl,k0,l,GMTR_A_HNX)
                grdy_pl(ij,k,l) = grdy_pl(ij,k,l) &
                                + dp_pl * GMTR_A_var_pl(ij_pl,k0,l,GMTR_A_HNY)
                grdz_pl(ij,k,l) = grdz_pl(ij,k,l) &
                                + dp_pl * GMTR_A_var_pl(ij_pl,k0,l,GMTR_A_HNZ)
             enddo

             grdx_pl(ij,k,l) = grdx_pl(ij,k,l) * GMTR_P_var_pl(ij,k0,l,GMTR_P_RAREA)
             grdy_pl(ij,k,l) = grdy_pl(ij,k,l) * GMTR_P_var_pl(ij,k0,l,GMTR_P_RAREA)
             grdz_pl(ij,k,l) = grdz_pl(ij,k,l) * GMTR_P_var_pl(ij,k0,l,GMTR_P_RAREA)

          enddo

       enddo
    endif

    return
  end subroutine OPRT3D_divdamp

  !-----------------------------------------------------------------------------
  subroutine OPRT3D_gradient_intpl( &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl,  &
       scl, scl_pl, &
       mfact        )
    use mod_adm, only: &
       ADM_prc_me,   &
       ADM_PRC_PL,   &
       ADM_rgn_vnum, &
       ADM_prc_tab,  &
       ADM_W,        &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_AI,       &
       ADM_AIJ,      &
       ADM_AJ,       &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_GSLF_PL,  &
       ADM_GMIN_PL,  &
       ADM_GMAX_PL,  &
       ADM_kmin,     &
       ADM_kmax,     &
       ADM_KNONE
    use mod_grd, only: &
       GRD_vz,     &
       GRD_vz_pl, &
       GRD_Z
    use mod_gmtr, only: &
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_A_HNX,      &
         GMTR_A_HNY,      &
         GMTR_A_HNZ,      &
         GMTR_P_RAREA,    &
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    implicit none

    real(8), intent(in)    :: scl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vx    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vy    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vz    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in), optional :: mfact

    real(8) :: flx6    (6)
    real(8) :: flx6_pl
    real(8) :: scl6    (6)
    real(8) :: scl6_pl (ADM_gall_pl)
    real(8) :: sclt6   (ADM_gall,6)
    real(8) :: sclt6_pl(ADM_gall_pl)

    real(8) :: fact

    integer :: nstart, nend
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1
    integer :: ij_pl, ijp1_pl, ijm1_pl

    integer :: k0
    integer :: rgnid
    integer :: ij, k, l, n

    real(8) :: p,z,z1,z2,z3,p1,p2,p3 ! interpolation
    p(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3) ) / ( (z1-z2)*(z1-z3) ) * p1 &
                           + ( (z-z1)*(z-z3) ) / ( (z2-z1)*(z2-z3) ) * p2 &
                           + ( (z-z1)*(z-z2) ) / ( (z3-z1)*(z3-z2) ) * p3
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do k = ADM_kmin, ADM_kmax
          nstart = suf(ADM_gmin,ADM_gmin)
          nend   = suf(ADM_gmax,ADM_gmax)

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ijp1   = n   + ADM_gall_1d
             ip1jp1 = n+1 + ADM_gall_1d
             im1j   = n-1
             ijm1   = n   - ADM_gall_1d
             im1jm1 = n-1 - ADM_gall_1d

             scl6(1) = p( GRD_vz(ij,    k,  l,GRD_Z),                    &
                          GRD_vz(ip1j,  k-1,l,GRD_Z), scl(ip1j,  k-1,l), &
                          GRD_vz(ip1j,  k,  l,GRD_Z), scl(ip1j,  k,  l), &
                          GRD_vz(ip1j,  k+1,l,GRD_Z), scl(ip1j,  k+1,l)  )

             scl6(2) = p( GRD_vz(ij,    k,  l,GRD_Z),                    &
                          GRD_vz(ip1jp1,k-1,l,GRD_Z), scl(ip1jp1,k-1,l), &
                          GRD_vz(ip1jp1,k,  l,GRD_Z), scl(ip1jp1,k,  l), &
                          GRD_vz(ip1jp1,k+1,l,GRD_Z), scl(ip1jp1,k+1,l)  )

             scl6(3) = p( GRD_vz(ij,    k,  l,GRD_Z),                    &
                          GRD_vz(ijp1,  k-1,l,GRD_Z), scl(ijp1,  k-1,l), &
                          GRD_vz(ijp1,  k,  l,GRD_Z), scl(ijp1,  k,  l), &
                          GRD_vz(ijp1,  k+1,l,GRD_Z), scl(ijp1,  k+1,l)  )

             scl6(4) = p( GRD_vz(ij,    k,  l,GRD_Z),                    &
                          GRD_vz(im1j,  k-1,l,GRD_Z), scl(im1j,  k-1,l), &
                          GRD_vz(im1j,  k,  l,GRD_Z), scl(im1j,  k,  l), &
                          GRD_vz(im1j,  k+1,l,GRD_Z), scl(im1j,  k+1,l)  )

             scl6(5) = p( GRD_vz(ij,    k,  l,GRD_Z),                    &
                          GRD_vz(im1jm1,k-1,l,GRD_Z), scl(im1jm1,k-1,l), &
                          GRD_vz(im1jm1,k,  l,GRD_Z), scl(im1jm1,k,  l), &
                          GRD_vz(im1jm1,k+1,l,GRD_Z), scl(im1jm1,k+1,l)  )

             scl6(6) = p( GRD_vz(ij,    k,  l,GRD_Z),                    &
                          GRD_vz(ijm1,  k-1,l,GRD_Z), scl(ijm1,  k-1,l), &
                          GRD_vz(ijm1,  k,  l,GRD_Z), scl(ijm1,  k,  l), &
                          GRD_vz(ijm1,  k+1,l,GRD_Z), scl(ijm1,  k+1,l)  )

             sclt6(ij,1) = GMTR_T_var(ij,    k0,l,ADM_TI,GMTR_T_W1) * scl(ij,k,l) &
                         + GMTR_T_var(ij,    k0,l,ADM_TI,GMTR_T_W2) * scl6(1)     &
                         + GMTR_T_var(ij,    k0,l,ADM_TI,GMTR_T_W3) * scl6(2)

             sclt6(ij,2) = GMTR_T_var(ij,    k0,l,ADM_TJ,GMTR_T_W1) * scl(ij,k,l) &
                         + GMTR_T_var(ij,    k0,l,ADM_TJ,GMTR_T_W2) * scl6(2)     &
                         + GMTR_T_var(ij,    k0,l,ADM_TJ,GMTR_T_W3) * scl6(3)

             sclt6(ij,3) = GMTR_T_var(im1j,  k0,l,ADM_TI,GMTR_T_W1) * scl6(4)     &
                         + GMTR_T_var(im1j,  k0,l,ADM_TI,GMTR_T_W2) * scl(ij,k,l) &
                         + GMTR_T_var(im1j,  k0,l,ADM_TI,GMTR_T_W3) * scl6(3)

             sclt6(ij,4) = GMTR_T_var(im1jm1,k0,l,ADM_TJ,GMTR_T_W1) * scl6(5)     &
                         + GMTR_T_var(im1jm1,k0,l,ADM_TJ,GMTR_T_W2) * scl(ij,k,l) &
                         + GMTR_T_var(im1jm1,k0,l,ADM_TJ,GMTR_T_W3) * scl6(4)

             sclt6(ij,5) = GMTR_T_var(im1jm1,k0,l,ADM_TI,GMTR_T_W1) * scl6(5)     &
                         + GMTR_T_var(im1jm1,k0,l,ADM_TI,GMTR_T_W2) * scl6(6)     &
                         + GMTR_T_var(im1jm1,k0,l,ADM_TI,GMTR_T_W3) * scl(ij,k,l)

             sclt6(ij,6) = GMTR_T_var(ijm1,  k0,l,ADM_TJ,GMTR_T_W1) * scl6(6)     &
                         + GMTR_T_var(ijm1,  k0,l,ADM_TJ,GMTR_T_W2) * scl6(1)     &
                         + GMTR_T_var(ijm1,  k0,l,ADM_TJ,GMTR_T_W3) * scl(ij,k,l)
          enddo

          do n = nstart, nend
             ij     = n
             im1j   = n-1
             ijm1   = n   - ADM_gall_1d
             im1jm1 = n-1 - ADM_gall_1d

             flx6(1) = 0.5D0 * ( sclt6(ij,6) + sclt6(ij,1) ) - scl(ij,k,l)
             flx6(2) = 0.5D0 * ( sclt6(ij,1) + sclt6(ij,2) ) - scl(ij,k,l)
             flx6(3) = 0.5D0 * ( sclt6(ij,2) + sclt6(ij,3) ) - scl(ij,k,l)
             flx6(4) = 0.5D0 * ( sclt6(ij,3) + sclt6(ij,4) ) - scl(ij,k,l)
             flx6(5) = 0.5D0 * ( sclt6(ij,4) + sclt6(ij,5) ) - scl(ij,k,l)
             flx6(6) = 0.5D0 * ( sclt6(ij,5) + sclt6(ij,6) ) - scl(ij,k,l)

             vx(ij,k,l) = ( flx6(1) * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNX) &
                          + flx6(2) * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNX) &
                          + flx6(3) * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNX) &
                          - flx6(4) * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNX) &
                          - flx6(5) * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNX) &
                          - flx6(6) * GMTR_A_var(ijm1,  k0,l,ADM_AJ, GMTR_A_HNX) &
                          ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             vy(ij,k,l) = ( flx6(1) * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNY) &
                          + flx6(2) * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNY) &
                          + flx6(3) * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNY) &
                          - flx6(4) * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNY) &
                          - flx6(5) * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNY) &
                          - flx6(6) * GMTR_A_var(ijm1,  k0,l,ADM_AJ, GMTR_A_HNY) &
                          ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             vz(ij,k,l) = ( flx6(1) * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNZ) &
                          + flx6(2) * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNZ) &
                          + flx6(3) * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNZ) &
                          - flx6(4) * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNZ) &
                          - flx6(5) * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNZ) &
                          - flx6(6) * GMTR_A_var(ijm1,  k0,l,ADM_AJ, GMTR_A_HNZ) &
                          ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact
          enddo

          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
             n = suf(ADM_gmin,ADM_gmin)

             ij     = n
             im1j   = n-1
             ijm1   = n   - ADM_gall_1d
             im1jm1 = n-1 - ADM_gall_1d

             flx6(1) = 0.5D0 * ( sclt6(ij,5) + sclt6(ij,1) ) - scl(ij,k,l)
             flx6(2) = 0.5D0 * ( sclt6(ij,1) + sclt6(ij,2) ) - scl(ij,k,l)
             flx6(3) = 0.5D0 * ( sclt6(ij,2) + sclt6(ij,3) ) - scl(ij,k,l)
             flx6(4) = 0.5D0 * ( sclt6(ij,3) + sclt6(ij,4) ) - scl(ij,k,l)
             flx6(5) = 0.5D0 * ( sclt6(ij,4) + sclt6(ij,5) ) - scl(ij,k,l)

             vx(ij,k,l) = ( flx6(1) * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNX) &
                          + flx6(2) * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNX) &
                          + flx6(3) * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNX) &
                          - flx6(4) * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNX) &
                          - flx6(5) * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNX) &
                          ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             vy(ij,k,l) = ( flx6(1) * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNY) &
                          + flx6(2) * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNY) &
                          + flx6(3) * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNY) &
                          - flx6(4) * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNY) &
                          - flx6(5) * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNY) &
                          ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             vz(ij,k,l) = ( flx6(1) * GMTR_A_var(ij,    k0,l,ADM_AI, GMTR_A_HNZ) &
                          + flx6(2) * GMTR_A_var(ij,    k0,l,ADM_AIJ,GMTR_A_HNZ) &
                          + flx6(3) * GMTR_A_var(ij,    k0,l,ADM_AJ, GMTR_A_HNZ) &
                          - flx6(4) * GMTR_A_var(im1j,  k0,l,ADM_AI, GMTR_A_HNZ) &
                          - flx6(5) * GMTR_A_var(im1jm1,k0,l,ADM_AIJ,GMTR_A_HNZ) &
                          ) * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact
          endif
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax

          ij = ADM_GSLF_PL

          vx_pl(ij,k,l) = 0.D0
          vy_pl(ij,k,l) = 0.D0
          vz_pl(ij,k,l) = 0.D0

          do n = ADM_GSLF_PL, ADM_GMAX_PL
             ij_pl = n

             scl6_pl(ij_pl) = p( GRD_vz_pl(ij,   k,  l,GRD_Z),                      &
                                 GRD_vz_pl(ij_pl,k-1,l,GRD_Z), scl_pl(ij_pl,k-1,l), & 
                                 GRD_vz_pl(ij_pl,k,  l,GRD_Z), scl_pl(ij_pl,k,  l), & 
                                 GRD_vz_pl(ij_pl,k+1,l,GRD_Z), scl_pl(ij_pl,k+1,l)  ) 
          enddo

          do n = ADM_GMIN_PL, ADM_GMAX_PL
             ij_pl   = n
             ijp1_pl = n+1

             if( ij_pl == ADM_GMAX_PL ) ijp1_pl = ADM_GMIN_PL ! cyclic condition

             sclt6_pl(ij_pl)  = GMTR_T_var_pl(ij_pl,k0,l,GMTR_T_W1) * scl6_pl(ij     ) &
                              + GMTR_T_var_pl(ij_pl,k0,l,GMTR_T_W2) * scl6_pl(ij_pl  ) &
                              + GMTR_T_var_pl(ij_pl,k0,l,GMTR_T_W3) * scl6_pl(ijp1_pl)
          enddo

          do n = ADM_GMIN_PL, ADM_GMAX_PL
             ij_pl   = n
             ijm1_pl = n-1

             if( ij_pl == ADM_GMIN_PL ) ijm1_pl = ADM_GMAX_PL ! cyclic condition

             flx6_pl = 0.5D0 * ( sclt6_pl(ijm1_pl) + sclt6_pl(ij_pl) ) - scl6_pl(ij)

             vx_pl(ij,k,l) = vx_pl(ij,k,l) + flx6_pl * GMTR_A_var_pl(n,k0,l,GMTR_A_HNX)
             vy_pl(ij,k,l) = vy_pl(ij,k,l) + flx6_pl * GMTR_A_var_pl(n,k0,l,GMTR_A_HNY)
             vz_pl(ij,k,l) = vz_pl(ij,k,l) + flx6_pl * GMTR_A_var_pl(n,k0,l,GMTR_A_HNZ)
          enddo

          vx_pl(ij,k,l) = vx_pl(ij,k,l) * GMTR_P_var_pl(ij,k0,l,GMTR_P_RAREA) * fact
          vy_pl(ij,k,l) = vy_pl(ij,k,l) * GMTR_P_var_pl(ij,k0,l,GMTR_P_RAREA) * fact
          vz_pl(ij,k,l) = vz_pl(ij,k,l) * GMTR_P_var_pl(ij,k0,l,GMTR_P_RAREA) * fact

       enddo
       enddo
    endif

    return
  end subroutine OPRT3D_gradient_intpl

end module mod_oprt3d
!-------------------------------------------------------------------------------
