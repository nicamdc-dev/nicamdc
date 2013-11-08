!-------------------------------------------------------------------------------
!>
!! variable conversion module
!!
!! @par Description
!!         Conversion tools for prognostic variables
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)   Imported from igdc-4.34
!! @li      2009-07-10 (H.Tomita)   Change the cnvvar_rhokin, cnvvar_kin for the energy conservation.
!! @li      2011-07-22 (T.Ohno)     add subroutines for plane hgrid systems
!!
!<
module mod_cnvvar
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
  public :: cnvvar_rhokin_ijkl

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
  subroutine cnvvar_rhokin_ijkl( &
       rhog,    rhog_pl,   &
       rhogvx,  rhogvx_pl, &
       rhogvy,  rhogvy_pl, &
       rhogvz,  rhogvz_pl, &
       rhogw,   rhogw_pl,  &
       rhogkin, rhogkin_pl )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_cfac, &
       GRD_dfac
    use mod_vmtr, only: &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl,  &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl
    implicit none

    real(8), intent(in)  :: rhog      (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 )
    real(8), intent(in)  :: rhog_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X vx
    real(8), intent(in)  :: rhogvx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X vy
    real(8), intent(in)  :: rhogvy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X vz
    real(8), intent(in)  :: rhogvz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X w
    real(8), intent(in)  :: rhogw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: rhogkin   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho X ( G^{1/2} X gamma2 ) X kin
    real(8), intent(out) :: rhogkin_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhogkin_h   (ADM_gall,   ADM_kall) ! rho X ( G^{1/2} X gamma2 ) X kin (horizontal)
    real(8) :: rhogkin_h_pl(ADM_gall_pl,ADM_kall)
    real(8) :: rhogkin_v   (ADM_gall,   ADM_kall) ! rho X ( G^{1/2} X gamma2 ) X kin (vertical)
    real(8) :: rhogkin_v_pl(ADM_gall_pl,ADM_kall)

    real(8) :: rhog_h   (ADM_gall,   ADM_kall)
    real(8) :: rhog_h_pl(ADM_gall_pl,ADM_kall)

    integer :: ij, k, l
    !---------------------------------------------------------------------------

!del
    do l = 1, ADM_lall
!del
       do k = ADM_kmin+1, ADM_kmax
          !--- rhog at the half level
!del
          do ij = 1, ADM_gall
             rhog_h(ij,k) = 0.5D0 * ( GRD_afac(k) * VMTR_RGSGAM2(ij,k  ,l) * rhog(ij,k  ,l) &
                                    + GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l) * rhog(ij,k-1,l) &
                                    ) * VMTR_GSGAM2H(ij,k,l)
          enddo
       enddo
!del
       do k = ADM_kmin+1, ADM_kmax
          !--- vertical kinetic energy
          do ij = 1, ADM_gall
             rhogkin_v(ij,k) = 0.5D0 * ( rhogw(ij,k,l) * rhogw(ij,k,l) ) / rhog_h(ij,k)
          enddo
       enddo
!del
       do ij = 1, ADM_gall
          rhogkin_v(ij,ADM_kmin  ) = 0.D0
          rhogkin_v(ij,ADM_kmax+1) = 0.D0
       enddo

!del
       do k = ADM_kmin, ADM_kmax
          !--- horizontal kinetic energy
          do ij = 1, ADM_gall
             rhogkin_h(ij,k) = 0.5D0 * ( rhogvx(ij,k,l) * rhogvx(ij,k,l) &
                                       + rhogvy(ij,k,l) * rhogvy(ij,k,l) &
                                       + rhogvz(ij,k,l) * rhogvz(ij,k,l) ) / rhog(ij,k,l)
          enddo
       enddo

!del
       do k = ADM_kmin, ADM_kmax
          !--- total kinetic energy
!del
          do ij = 1, ADM_gall
             rhogkin(ij,k,l) = rhogkin_h(ij,k)                           &
                             + 0.5D0 * ( GRD_dfac(k) * rhogkin_v(ij,k+1) &
                                       + GRD_cfac(k) * rhogkin_v(ij,k  ) )
          enddo
       enddo
!del
       do ij = 1, ADM_gall
          rhogkin(ij,ADM_kmin-1,l) = 0.D0
          rhogkin(ij,ADM_kmax+1,l) = 0.D0
       enddo

    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!del
       do l = 1, ADM_lall_pl
!del
          do k = ADM_kmin+1, ADM_kmax
             !--- rhog at the half level
!del
             do ij = 1, ADM_gall_pl
                rhog_h_pl(ij,k) = 0.5D0 * ( GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k  ,l) * rhog_pl(ij,k  ,l) &
                                          + GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l) * rhog_pl(ij,k-1,l) &
                                          ) * VMTR_GSGAM2H_pl(ij,k,l)
             enddo
          enddo
!del
          do k = ADM_kmin+1, ADM_kmax
             !--- vertical kinetic energy
             do ij = 1, ADM_gall_pl
                rhogkin_v_pl(ij,k) = 0.5D0 * ( rhogw_pl(ij,k,l) * rhogw_pl(ij,k,l) ) / rhog_h_pl(ij,k)
             enddo
          enddo
!del
          do ij = 1, ADM_gall_pl
             rhogkin_v_pl(ij,ADM_kmin  ) = 0.D0
             rhogkin_v_pl(ij,ADM_kmax+1) = 0.D0
          enddo
!del
          do k = ADM_kmin, ADM_kmax
             !--- horizontal kinetic energy
             do ij = 1, ADM_gall_pl
                rhogkin_h_pl(ij,k) = 0.5D0 * ( rhogvx_pl(ij,k,l) * rhogvx_pl(ij,k,l) &
                                             + rhogvy_pl(ij,k,l) * rhogvy_pl(ij,k,l) &
                                             + rhogvz_pl(ij,k,l) * rhogvz_pl(ij,k,l) ) / rhog_pl(ij,k,l)
             enddo
          enddo

!del
          do k = ADM_kmin, ADM_kmax
             !--- total kinetic energy
             do ij = 1, ADM_gall_pl
                rhogkin_pl(ij,k,l) = rhogkin_h_pl(ij,k)                           &
                                   + 0.5D0 * ( GRD_dfac(k) * rhogkin_v_pl(ij,k+1) &
                                             + GRD_cfac(k) * rhogkin_v_pl(ij,k  ) )
             enddo
          enddo
          do ij = 1, ADM_gall_pl
             rhogkin_pl(ij,ADM_kmin-1,l) = 0.D0
             rhogkin_pl(ij,ADM_kmax+1,l) = 0.D0
          enddo

       enddo
    endif

    return
  end subroutine cnvvar_rhokin_ijkl

end module mod_cnvvar
!-------------------------------------------------------------------------------
