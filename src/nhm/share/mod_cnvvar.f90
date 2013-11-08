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
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl, &
       VMTR_W2Cfact,    &
       VMTR_W2Cfact_pl
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

    integer :: n, k, l
    !---------------------------------------------------------------------------
print *, 'OK0'
    do l = 1, ADM_lall
       !--- horizontal kinetic energy
       do k = ADM_kmin, ADM_kmax
       do n = 1, ADM_gall
          rhogkin_h(n,k) = 0.5D0 * ( rhogvx(n,k,l) * rhogvx(n,k,l) &
                                   + rhogvy(n,k,l) * rhogvy(n,k,l) &
                                   + rhogvz(n,k,l) * rhogvz(n,k,l) ) / rhog(n,k,l)
       enddo
       enddo

       !--- vertical kinetic energy
       do k = ADM_kmin+1, ADM_kmax
       do n = 1, ADM_gall
          rhogkin_v(n,k) = 0.5D0 * ( rhogw(n,k,l) * rhogw(n,k,l) ) &
                         / ( VMTR_C2Wfact(1,n,k,l) * rhog(n,k  ,l) &
                           + VMTR_C2Wfact(2,n,k,l) * rhog(n,k-1,l) )
       enddo
       enddo
       rhogkin_v(:,ADM_kmin  ) = 0.D0
       rhogkin_v(:,ADM_kmax+1) = 0.D0

       !--- total kinetic energy
       do k = ADM_kmin, ADM_kmax
       do n = 1, ADM_gall
          rhogkin(n,k,l) = rhogkin_h(n,k)                             & ! horizontal
                         + ( VMTR_W2Cfact(1,n,k,l) * rhogkin_v(n,k+1) & ! vertical
                           + VMTR_W2Cfact(2,n,k,l) * rhogkin_v(n,k  ) )
       enddo
       enddo
       rhogkin(:,ADM_kmin-1,l) = 0.D0
       rhogkin(:,ADM_kmax+1,l) = 0.D0
    enddo
print *, 'OK1'

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          !--- horizontal kinetic energy
          do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall_pl
             rhogkin_h_pl(n,k) = 0.5D0 * ( rhogvx_pl(n,k,l) * rhogvx_pl(n,k,l) &
                                         + rhogvy_pl(n,k,l) * rhogvy_pl(n,k,l) &
                                         + rhogvz_pl(n,k,l) * rhogvz_pl(n,k,l) ) / rhog_pl(n,k,l)
          enddo
          enddo

          !--- vertical kinetic energy
          do k = ADM_kmin+1, ADM_kmax
          do n = 1, ADM_gall_pl
             rhogkin_v_pl(n,k) = 0.5D0 * ( rhogw_pl(n,k,l) * rhogw_pl(n,k,l) ) &
                               / ( VMTR_C2Wfact_pl(1,n,k,l) * rhog_pl(n,k  ,l) &
                                 + VMTR_C2Wfact_pl(2,n,k,l) * rhog_pl(n,k-1,l) )
          enddo
          enddo
          rhogkin_v_pl(:,ADM_kmin  ) = 0.D0
          rhogkin_v_pl(:,ADM_kmax+1) = 0.D0

          !--- total kinetic energy
          do k = ADM_kmin, ADM_kmax
          do n = 1, ADM_gall_pl
             rhogkin_pl(n,k,l) = rhogkin_h_pl(n,k)                                & ! horizontal
                               + ( VMTR_W2Cfact_pl(1,n,k,l) * rhogkin_v_pl(n,k+1) & ! vertical
                                 + VMTR_W2Cfact_pl(2,n,k,l) * rhogkin_v_pl(n,k  ) )
          enddo
          enddo
          rhogkin_pl(:,ADM_kmin-1,l) = 0.D0
          rhogkin_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif
print *, 'OK2'

    return
  end subroutine cnvvar_rhokin_ijkl

end module mod_cnvvar
!-------------------------------------------------------------------------------
