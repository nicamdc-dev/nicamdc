!-------------------------------------------------------------------------------
!
!+  Source calculation module
!
!-------------------------------------------------------------------------------
module mod_src
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the caluculation of source terms
  !       in the nhm-model.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                06-08-11   Add sub[src_update_tracer] for tracer advection.
  !                07-04-25   K.Suzuki: branch for save_memory in src_update_tracer
  !                07-11-13   H.Tomita: Change the vertical advection limiter.
  !                                     apply the limiter from rho q to q.
  !                07-11-14   Y.Niwa: bug fix /rho_h => /rho_h_pl
  !                08-01-24   Y.Niwa: add src_update_tracer
  !                                   save_memory => TRC_SAVE_MEMORY
  !                08-04-28   Y.Niwa: add initialization and avoid zero-dividing. 
  !                09-05-26   Y.Yamada: add directive, derected by T.Asano
  !                11-09-27   T.Seiki: merge optimized routines for K by RIST and M.Terai
  !                12-03-29   T.Yamaura: optimized for K
  !                12-05-30   T.Yashiro: Change arguments from character to index/switch
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
  integer, public, parameter :: I_SRC_horizontal = 1 ! [add] H.Yashiro 20120530
  integer, public, parameter :: I_SRC_vertical   = 2 ! [add] H.Yashiro 20120530
  integer, public, parameter :: I_SRC_default    = 3 ! [add] H.Yashiro 20120530

  real(8), private, allocatable, save :: AFACovGSGAM2   (:,:,:)
  real(8), private, allocatable, save :: AFACovGSGAM2_pl(:,:,:)
  real(8), private, allocatable, save :: BFACovGSGAM2   (:,:,:)
  real(8), private, allocatable, save :: BFACovGSGAM2_pl(:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: src_flux_convergence
  public :: src_advection_convergence
  public :: src_buoyancy
  public :: src_pres_work_uv
  public :: src_advection_convergence_v
  public :: src_gradient
  public :: src_update_tracer

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  logical, private, parameter :: first_layer_remedy = .true.

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  public :: src_flux_convergence_h

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !------ 
  !------ Flux convergence calculation
  !------     1. Horizontal flux convergence is calculated by using
  !------        rhovx, rhovy, and rhovz which are defined at cell
  !------        center (verticaly) and A-grid (horizontally).
  !------     2. Vertical flux convergence is calculated by using
  !------        rhovx, rhovy, rhovz, and rhow.
  !------     3. rhovx, rhovy, and rhovz can be replaced by
  !------        rhovx*h, rhovy*h, and rhovz*h, respectively.
  !------     4. f_type can be set as below.
  !------     5. Calculation region of grho, grho_pl 
  !------             : (:, ADM_kmin:ADM_kmax,:)
  !
  subroutine src_flux_convergence( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       grhog,  grhog_pl,  &
!       f_type ) [mod] H.Yashiro 20120530
       fluxtype           )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_GSGAMH,     &
       VMTR_GSGAMH_pl,  &
       VMTR_RGSH,       &
       VMTR_RGSH_pl,    &
       VMTR_RGAM,       &
       VMTR_RGAM_pl,    &
       VMTR_GZXH,       &
       VMTR_GZXH_pl,    &
       VMTR_GZYH,       &
       VMTR_GZYH_pl,    &
       VMTR_GZZH,       &
       VMTR_GZZH_pl
    use mod_grd, only: &
       GRD_gzh,  &
       GRD_gz,   &
       GRD_afac, &
       GRD_bfac
    use mod_oprt, only: &
       OPRT_divergence
    implicit none

    real(8), intent(in) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz ( gam2 X G^{1/2} )
    real(8), intent(in) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out) :: grhog   (ADM_gall,   ADM_kall,ADM_lall   ) ! source
    real(8), intent(out) :: grhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

!    character(LEN=*),intent(in) :: f_type
!    !<-----                        = 'HORIZONTAL' : horizontal convergence
!    !<-----                        = 'VERTICAL'   : vertical convergence
!    !<-----                        = 'DEFAULT'    : both of them
    integer, intent(in) :: fluxtype ! [mod] H.Yashiro 20120530
    !<-----                = I_SRC_horizontal : horizontal convergence
    !<-----                = I_SRC_vertical   : vertical convergence
    !<-----                = I_SRC_default    : both of them

    real(8) :: flx_vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: hdiv       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: hdiv_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvx_f   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_f_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvy_f   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_f_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvz_f   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_f_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    logical :: iflag = .true.

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    if ( iflag ) then
       iflag = .false.

       allocate( AFACovGSGAM2   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( AFACovGSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( BFACovGSGAM2   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( BFACovGSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )

       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall
                   AFACovGSGAM2(ij,k,l) = GRD_afac(k) * VMTR_RGSGAM2(ij,k  ,l)
                   BFACovGSGAM2(ij,k,l) = GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do ij = 1, ADM_gall_pl
                      AFACovGSGAM2_pl(ij,k,l) = GRD_afac(k) * VMTR_RGSGAM2_pl(ij,k  ,l)
                      BFACovGSGAM2_pl(ij,k,l) = GRD_bfac(k) * VMTR_RGSGAM2_pl(ij,k-1,l)
                enddo
             enddo
          enddo
       endif

    endif ! only once

    !---------------------------------------------------------------------------
    ! Horizontal
    !---------------------------------------------------------------------------
!    if(f_type=='HORIZONTAL') then
    if ( fluxtype == I_SRC_horizontal ) then ! [mod] H.Yashiro 20120530

       !--- horizontal contribution to flx_vz
       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall
                flx_vz(ij,k,l) = ( ( AFACovGSGAM2(ij,k,l) * rhogvx(ij,k  ,l)           &
                                   + BFACovGSGAM2(ij,k,l) * rhogvx(ij,k-1,l)           &
                                   ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZXH(ij,k,l) &
                                 + ( AFACovGSGAM2(ij,k,l) * rhogvy(ij,k  ,l)           &
                                   + BFACovGSGAM2(ij,k,l) * rhogvy(ij,k-1,l)           &
                                   ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZYH(ij,k,l) & 
                                 + ( AFACovGSGAM2(ij,k,l) * rhogvz(ij,k  ,l)           &
                                   + BFACovGSGAM2(ij,k,l) * rhogvz(ij,k-1,l)           &
                                   ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZZH(ij,k,l) &
                                 )
             enddo
          enddo

          !--- vertical flux is zero at the top and bottom boundary.
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmin,l)   = 0.D0
          enddo
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do ij = 1, ADM_gall_pl
                   flx_vz_pl(ij,k,l) = ( ( AFACovGSGAM2_pl(ij,k,l) * rhogvx_pl(ij,k  ,l)           &
                                         + BFACovGSGAM2_pl(ij,k,l) * rhogvx_pl(ij,k-1,l)           &
                                         ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZXH_pl(ij,k,l) &
                                       + ( AFACovGSGAM2_pl(ij,k,l) * rhogvy_pl(ij,k  ,l)           &
                                         + BFACovGSGAM2_pl(ij,k,l) * rhogvy_pl(ij,k-1,l)           &
                                         ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZYH_pl(ij,k,l) &
                                       + ( AFACovGSGAM2_pl(ij,k,l) * rhogvz_pl(ij,k  ,l)           &
                                         + BFACovGSGAM2_pl(ij,k,l) * rhogvz_pl(ij,k-1,l)           &
                                         ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZZH_pl(ij,k,l) &
                                       )
                enddo
             enddo
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmin,l)   = 0.D0
             enddo
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmax+1,l) = 0.D0
             enddo
          enddo
       endif

       !--- Horizontal flux convergence
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do ij = 1, ADM_gall
                rhogvx_f(ij,k,l) = rhogvx(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
             do ij = 1, ADM_gall
                rhogvy_f(ij,k,l) = rhogvy(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
             do ij = 1, ADM_gall
                rhogvz_f(ij,k,l) = rhogvz(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do ij = 1, ADM_gall_pl
                   rhogvx_f_pl(ij,k,l) = rhogvx_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
                do ij = 1, ADM_gall_pl
                   rhogvy_f_pl(ij,k,l) = rhogvy_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
                do ij = 1, ADM_gall_pl
                   rhogvz_f_pl(ij,k,l) = rhogvz_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
             enddo
          enddo
       endif

       call OPRT_divergence( hdiv,     hdiv_pl,     & !--- [OUT]
                             rhogvx_f, rhogvx_f_pl, & !--- [IN]
                             rhogvy_f, rhogvy_f_pl, & !--- [IN]
                             rhogvz_f, rhogvz_f_pl  ) !--- [IN]

    !---------------------------------------------------------------------------
    ! Vertical
    !---------------------------------------------------------------------------
!    else if(f_type=='VERTICAL') then
    elseif( fluxtype == I_SRC_vertical ) then ! [mod] H.Yashiro 20120530

       !--- vertical contribution to flx_vz
       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall
                flx_vz(ij,k,l) = rhogw(ij,k,l) * VMTR_RGSH(ij,k,l)
             enddo
          enddo
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmin  ,l) = 0.D0
          enddo
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do ij = 1, ADM_gall_pl
                   flx_vz_pl(ij,k,l) = rhogw_pl(ij,k,l) * VMTR_RGSH_pl(ij,k,l)
                enddo
             enddo
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmin  ,l) = 0.D0
             enddo
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmax+1,l) = 0.D0
             enddo
          enddo
       endif

       !--- Horizontal flux convergence
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do ij = 1, ADM_gall
                hdiv(ij,k,l) = 0.D0
             enddo
          enddo
       enddo
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                hdiv_pl(ij,k,l) = 0.D0
             enddo
          enddo
       enddo

    !---------------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------------
    else ! f_type=='DEFAULT'

       !--- horizontal + vertical contribution to flx_vz
       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             do ij = 1, ADM_gall
                flx_vz(ij,k,l) = ( ( AFACovGSGAM2(ij,k,l) * rhogvx(ij,k  ,l)           &
                                   + BFACovGSGAM2(ij,k,l) * rhogvx(ij,k-1,l)           &
                                   ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZXH(ij,k,l) &
                                 + ( AFACovGSGAM2(ij,k,l) * rhogvy(ij,k  ,l)           &
                                   + BFACovGSGAM2(ij,k,l) * rhogvy(ij,k-1,l)           &
                                   ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZYH(ij,k,l) & 
                                 + ( AFACovGSGAM2(ij,k,l) * rhogvz(ij,k  ,l)           &
                                   + BFACovGSGAM2(ij,k,l) * rhogvz(ij,k-1,l)           &
                                   ) * 0.5D0 * VMTR_GSGAMH(ij,k,l) * VMTR_GZZH(ij,k,l) &
                                 ) + rhogw(ij,k,l) * VMTR_RGSH(ij,k,l) ! vertical contribution
             enddo
          enddo
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmin  ,l) = 0.D0
          enddo
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do ij = 1, ADM_gall_pl
                   flx_vz_pl(ij,k,l) = ( ( AFACovGSGAM2_pl(ij,k,l) * rhogvx_pl(ij,k  ,l)           &
                                         + BFACovGSGAM2_pl(ij,k,l) * rhogvx_pl(ij,k-1,l)           &
                                         ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZXH_pl(ij,k,l) &
                                       + ( AFACovGSGAM2_pl(ij,k,l) * rhogvy_pl(ij,k  ,l)           &
                                         + BFACovGSGAM2_pl(ij,k,l) * rhogvy_pl(ij,k-1,l)           &
                                         ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZYH_pl(ij,k,l) &
                                       + ( AFACovGSGAM2_pl(ij,k,l) * rhogvz_pl(ij,k  ,l)           &
                                         + BFACovGSGAM2_pl(ij,k,l) * rhogvz_pl(ij,k-1,l)           &
                                         ) * 0.5D0 * VMTR_GSGAMH_pl(ij,k,l) * VMTR_GZZH_pl(ij,k,l) &
                                       ) + rhogw_pl(ij,k,l) * VMTR_RGSH_pl(ij,k,l) ! vertical contribution
                enddo
             enddo
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmin  ,l) = 0.D0
             enddo
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmax+1,l) = 0.D0
             enddo
          enddo
       endif

       !--- Horizontal flux convergence
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do ij = 1, ADM_gall
                rhogvx_f(ij,k,l) = rhogvx(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
             do ij = 1, ADM_gall
                rhogvy_f(ij,k,l) = rhogvy(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
             do ij = 1, ADM_gall
                rhogvz_f(ij,k,l) = rhogvz(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do ij = 1, ADM_gall_pl
                   rhogvx_f_pl(ij,k,l) = rhogvx_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
                do ij = 1, ADM_gall_pl
                   rhogvy_f_pl(ij,k,l) = rhogvy_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
                do ij = 1, ADM_gall_pl
                   rhogvz_f_pl(ij,k,l) = rhogvz_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
             enddo
          enddo
       endif

       call OPRT_divergence( hdiv,     hdiv_pl,     & !--- [OUT]
                             rhogvx_f, rhogvx_f_pl, & !--- [IN]
                             rhogvy_f, rhogvy_f_pl, & !--- [IN]
                             rhogvz_f, rhogvz_f_pl  ) !--- [IN]

    endif

    !--- Total flux convergence
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          do ij = 1, ADM_gall
             grhog(ij,k,l) = - ( flx_vz(ij,k+1,l) - flx_vz(ij,k,l) ) &
                             / ( GRD_gzh(k+1)     - GRD_gzh(k)     ) - hdiv(ij,k,l)
          enddo
       enddo
       do ij = 1, ADM_gall
          grhog(ij,ADM_kmin-1,l) = 0.D0
       enddo
       do ij = 1, ADM_gall
          grhog(ij,ADM_kmax+1,l) = 0.D0
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do ij = 1, ADM_gall_pl
                grhog_pl(ij,k,l) = - ( flx_vz_pl(ij,k+1,l) - flx_vz_pl(ij,k,l) ) &
                                   / ( GRD_gzh(k+1)        - GRD_gzh(k)        ) - hdiv_pl(ij,k,l)
             enddo
          enddo
          do ij = 1, ADM_gall_pl
             grhog_pl(ij,ADM_kmin-1,l) = 0.D0
          enddo
          do ij = 1, ADM_gall_pl
             grhog_pl(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo
    endif

    return
  end subroutine src_flux_convergence

  !-----------------------------------------------------------------------------
  subroutine src_flux_convergence_h(&
         rhogvx, rhogvx_pl,        & !--- IN : rho*Vx  ( gam2 X G^{1/2} )
         rhogvy, rhogvy_pl,        & !--- IN : rho*Vy  ( gam2 X G^{1/2} )
         rhogvz, rhogvz_pl,        & !--- IN : rho*Vz  ( gam2 X G^{1/2} )
         rhogw,  rhogw_pl,         & !--- IN : rho*w   ( gam2 X G^{1/2} )
         grhog,  grhog_pl,         & !--- OUT: source
!         f_type )                    !--- IN : flux type
         fluxtype                  ) !--- IN : flux type [mod] H.Yashiro 20120530
    !------ 
    !------ Flux convergence calculation
    !------     1. Horizontal flux convergence is calculated by using
    !------        rhovx, rhovy, and rhovz which are defined at cell
    !------        center (verticaly) and A-grid (horizontally).
    !------     2. Vertical flux convergence is calculated by using
    !------        rhovx, rhovy, rhovz, and rhow.
    !------     3. rhovx, rhovy, and rhovz can be replaced by
    !------        rhovx*h, rhovy*h, and rhovz*h, respectively.
    !------     4. f_type can be set as below.
    !------     5. Calculation region of grho, grho_pl 
    !------             : (:, ADM_kmin:ADM_kmax,:)
    !
    use mod_adm, only :  &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_GALL_PL,    &
         ADM_LALL_PL,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_vmtr, only : &
         VMTR_RGSGAM2H,   &
         VMTR_RGSGAM2H_pl,&
         VMTR_GZX,       &
         VMTR_GZX_pl,    &
         VMTR_GZY,       &
         VMTR_GZY_pl,    &
         VMTR_GZZ,       &
         VMTR_GZZ_pl,    &
         VMTR_GSGAMH,    &
         VMTR_GSGAMH_pl, &
         VMTR_RGSH,      &
         VMTR_RGSH_pl,   &
         VMTR_RGAM,      &
         VMTR_RGAM_pl
    use mod_grd, only :  &
         GRD_gz,         &
         GRD_cfac,       &
         GRD_dfac
    use mod_oprt, only :        &
         OPRT_divergence
    !
    implicit none
    !
    real(8), intent(in) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(out) :: grhog(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
!    character(LEN=*),intent(in) :: f_type
    !<-----                    = 'HORIZONTAL' : horizontal convergence
    !<-----                    = 'VERTICAL'   : vertical convergence
    !<-----                    = 'DEFAULT'    : both of them
    integer, intent(in) :: fluxtype ! [mod] H.Yashiro 20120530
    !<-----                    = I_SRC_horizontal : horizontal convergence
    !<-----                    = I_SRC_vertical   : vertical convergence
    !<-----                    = I_SRC_default    : both of them

    real(8) :: flx_vz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: flx_vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8) :: hdiv(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: hdiv_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhogvx_f(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvx_f_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhogvy_f(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvy_f_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhogvz_f(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvz_f_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    integer :: l,k
    !
    !
    !--- horizontal contribution to flx_vz
    flx_vz = 0.0D0
    flx_vz_pl = 0.0D0
!    if(f_type=='VERTICAL') then
    if ( fluxtype == I_SRC_vertical ) then ! [mod] H.Yashiro 20120530
       flx_vz(:,:,:) = 0.0D0
       flx_vz_pl(:,:,:) = 0.0D0
    else !--- default or 'HORIZONTAL'
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax
             flx_vz(:,k,l) = (                                            &
                  ( GRD_cfac(k)*VMTR_RGSGAM2H(:,k+1,l)*rhogvx(:,k+1,  l)  &
                   +GRD_dfac(k)*VMTR_RGSGAM2H(:,k  ,l)*rhogvx(:,k  ,  l) )&
                  * 0.5D0*VMTR_GSGAMH(:,k,l)*VMTR_GZX(:,k,l)               &
                 +( GRD_cfac(k)*VMTR_RGSGAM2H(:,k+1,l)*rhogvy(:,k+1,  l)  &
                   +GRD_dfac(k)*VMTR_RGSGAM2H(:,k  ,l)*rhogvy(:,k  ,  l) )&
                  * 0.5D0*VMTR_GSGAMH(:,k,l)*VMTR_GZY(:,k,l)               &
                 +( GRD_cfac(k)*VMTR_RGSGAM2H(:,k+1,l)*rhogvz(:,k+1,  l)  &
                   +GRD_dfac(k)*VMTR_RGSGAM2H(:,k  ,l)*rhogvz(:,k  ,  l) )&
                  * 0.5D0*VMTR_GSGAMH(:,k,l)*VMTR_GZZ(:,k,l)               &
                  )
          end do
       end do
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_LALL_PL
             do k = ADM_kmin, ADM_kmax
                flx_vz_pl(:,k,l) = (                                            &
                     ( GRD_cfac(k)*VMTR_RGSGAM2H_pl(:,k+1,l)*rhogvx_pl(:,k+1,  l)  &
                     +GRD_dfac(k)*VMTR_RGSGAM2H_pl(:,k  ,l)*rhogvx_pl(:,k  ,  l) )&
                     * 0.5D0*VMTR_GSGAMH_pl(:,k,l)*VMTR_GZX_pl(:,k,l)               &
                    +( GRD_cfac(k)*VMTR_RGSGAM2H_pl(:,k+1,l)*rhogvy_pl(:,k+1,  l)  &
                     +GRD_dfac(k)*VMTR_RGSGAM2H_pl(:,k  ,l)*rhogvy_pl(:,k  ,  l) )&
                     * 0.5D0*VMTR_GSGAMH_pl(:,k,l)*VMTR_GZY_pl(:,k,l)               &
                    +( GRD_cfac(k)*VMTR_RGSGAM2H_pl(:,k+1,l)*rhogvz_pl(:,k+1,  l)  &
                     +GRD_dfac(k)*VMTR_RGSGAM2H_pl(:,k  ,l)*rhogvz_pl(:,k  ,  l) )&
                     * 0.5D0*VMTR_GSGAMH_pl(:,k,l)*VMTR_GZZ_pl(:,k,l)               &
                     )
             end do
          end do
       end if
    end if
    !
    !--- vertical contribution to flx_vz
!    if(f_type/='HORIZONTAL') then
    if ( fluxtype /= I_SRC_horizontal ) then ! [mod] H.Yashiro 20120530
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax
             flx_vz(:,k,l)                                    &
                  = flx_vz(:,k,l)                             &
                  + rhogw(:,k,l) * VMTR_RGSH(:,k,l)
          end do
       end do
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_LALL_PL
             do k = ADM_kmin, ADM_kmax
                flx_vz_pl(:,k,l)                             &
                  = flx_vz_pl(:,k,l)                      &
                  + rhogw_pl(:,k,l) * VMTR_RGSH_pl(:,k,l)
             end do
          end do
       end if
    end if
    !
    !--- vertical flux is zero at the top and bottom boundary.
    flx_vz(:,ADM_kmax+1,:) = 0.0D0
    if(ADM_prc_me==ADM_prc_pl) then
       flx_vz_pl(:,ADM_kmax+1,:) = 0.0D0
    end if
    !
    !--- Horizontal flux convergence
!    if(f_type/='VERTICAL') then
    if ( fluxtype /= I_SRC_vertical ) then ! [mod] H.Yashiro 20120530
       hdiv(:,:,:) = 0.0D0
       hdiv_pl(:,:,:) = 0.0D0
       !
       rhogvx_f(:,:,:) = rhogvx(:,:,:)*VMTR_RGAM(:,:,:)
       rhogvy_f(:,:,:) = rhogvy(:,:,:)*VMTR_RGAM(:,:,:)
       rhogvz_f(:,:,:) = rhogvz(:,:,:)*VMTR_RGAM(:,:,:)
       !
       if(ADM_prc_me==ADM_prc_pl) then
          rhogvx_f_pl(:,:,:) = rhogvx_pl(:,:,:)*VMTR_RGAM_pl(:,:,:)
          rhogvy_f_pl(:,:,:) = rhogvy_pl(:,:,:)*VMTR_RGAM_pl(:,:,:)
          rhogvz_f_pl(:,:,:) = rhogvz_pl(:,:,:)*VMTR_RGAM_pl(:,:,:)
       end if
       
       call OPRT_divergence(      &
            hdiv,hdiv_pl,         &  !--- out
            rhogvx_f, rhogvx_f_pl,&  !--- in
            rhogvy_f, rhogvy_f_pl,&  !--- in
            rhogvz_f, rhogvz_f_pl &  !--- in
            )
    else
       hdiv(:,:,:) = 0.0D0
       hdiv_pl(:,:,:) = 0.0D0
    end if
    !
    !--- Total flux convergence
    do l=1,ADM_lall
       do k = ADM_kmin+1, ADM_kmax
          grhog(:,k,l)                                  &
               = - ( flx_vz(:,k,l) - flx_vz(:,k-1,l) ) &
                 / ( GRD_gz(k) - GRD_gz(k-1) )       &
                 - hdiv(:,k,l)
       end do
    end do
    grhog(:,ADM_kmin,:) = 0.0D0
    grhog(:,ADM_kmin-1,:) = 0.0D0
    grhog(:,ADM_kmax+1,:) = 0.0D0
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
             grhog_pl(:,k,l)                                  &
                  = - ( flx_vz_pl(:,k,l) - flx_vz_pl(:,k-1,l) ) &
                  / ( GRD_gz(k) - GRD_gz(k-1) )       &
                  - hdiv_pl(:,k,l)
          end do
       end do
       grhog_pl(:,ADM_kmin,:) = 0.0D0
       grhog_pl(:,ADM_kmin-1,:) = 0.0D0
       grhog_pl(:,ADM_kmax+1,:) = 0.0D0
    end if
    !
    !
  end subroutine src_flux_convergence_h

  !-----------------------------------------------------------------------------
  subroutine src_buoyancy(&
       rhog, rhog_pl,     & !--- IN : perturb density ( gam2 X G^{1/2} )
       gbz, gbz_pl )        !--- OUT : buoyancy force  at half level
    !------ 
    !------ Calculation of buoyacy force
    !------     1. Calculation region of gbz, gbz_pl
    !------                    : (:, ADM_kmin:ADM_kmax+1,:)
    !------ 
    !
    use mod_adm, only :  &
         ADM_VMISS,      &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_GALL_PL,    &
         ADM_LALL_PL,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_vmtr, only : &
         VMTR_RGSGAM2,   &
         VMTR_RGSGAM2_pl,&
         VMTR_GSGAM2H,   &
         VMTR_GSGAM2H_pl
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac
    use mod_cnst, only : &
         CNST_EGRAV
    !
    implicit none
    !
    real(8), intent(in) :: rhog(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gbz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gbz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    integer :: g,k,l





    !---- NOTICE! : Upward direction is positive for gbz.
    do l = 1, ADM_lall
       do g = 1, ADM_gall
          gbz(g,ADM_kmin-1,l) = ADM_VMISS
       end do
    end do
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall
             gbz(g,k,l) =                               &
                  - ( GRD_afac(k)*rhog(g,k,l)*VMTR_RGSGAM2(g,k,l)     &
                  + GRD_bfac(k)*rhog(g,k-1,l)*VMTR_RGSGAM2(g,k-1,l)  ) &
                  * VMTR_GSGAM2H(g,k,l) &
                  * CNST_EGRAV * 0.5D0
          end do
       end do
    end do
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do g = 1, ADM_gall_pl
             gbz_pl(g,ADM_kmin-1,l) = ADM_VMISS
          end do
       end do
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall_pl
                gbz_pl(g,k,l) =                               &
                     - ( GRD_afac(k)*rhog_pl(g,k,l)*VMTR_RGSGAM2_pl(g,k,l)     &
                     + GRD_bfac(k)*rhog_pl(g,k-1,l)*VMTR_RGSGAM2_pl(g,k-1,l)  ) &
                     * VMTR_GSGAM2H_pl(g,k,l)&
                     * CNST_EGRAV * 0.5D0
             end do
          end do
       end do
    end if





  end subroutine src_buoyancy

  !-----------------------------------------------------------------------------
  !------ Advection convergence
  !------     1. f_type is for sub[flux_convergence]
  !------     2. vint_type is only for vertical direction.
  !------        (current version)
  !------     3. If e_w is present, it is used in the vertical flux
  !------        calculation.
  !------     4. Calculation region of grhoe, grhoe_pl
  !------                    : (:,ADM_kmin:ADM_kmax,:)
  !------ 
  subroutine src_advection_convergence( &  
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       e,      e_pl,      &
       grhoge, grhoge_pl, &
!       f_type  [mod] H.Yashiro 20120530
       fluxtype           )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    implicit none

    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: e        (ADM_gall,   ADM_kall,ADM_lall   ) ! arbitrary scalar
    real(8), intent(in)  :: e_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhoge   (ADM_gall,   ADM_kall,ADM_lall   ) ! tendency
    real(8), intent(out) :: grhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

!    character(LEN=*),intent(in) :: f_type
!    !<-----                        = 'HORIZONTAL' : horizontal convergence
!    !<-----                        = 'VERTICAL'   : vertical convergence
!    !<-----                        = 'DEFAULT'    : both of them
    integer, intent(in) :: fluxtype ! [mod] H.Yashiro 20120530
    !<-----                = I_SRC_horizontal : horizontal convergence
    !<-----                = I_SRC_vertical   : vertical convergence
    !<-----                = I_SRC_default    : both of them

    real(8) :: fez   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: fez_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8) :: rhogvxe   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvxe_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvye   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvye_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvze   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvze_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: ij, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do ij = 1, ADM_gall
             rhogvxe(ij,k,l) = rhogvx(ij,k,l) * e(ij,k,l)
          enddo
          do ij = 1, ADM_gall
             rhogvye(ij,k,l) = rhogvy(ij,k,l) * e(ij,k,l)
          enddo
          do ij = 1, ADM_gall
             rhogvze(ij,k,l) = rhogvz(ij,k,l) * e(ij,k,l)
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do ij = 1, ADM_gall_pl
                rhogvxe_pl(ij,k,l) = rhogvx_pl(ij,k,l) * e_pl(ij,k,l)
             enddo
             do ij = 1, ADM_gall_pl
                rhogvye_pl(ij,k,l) = rhogvy_pl(ij,k,l) * e_pl(ij,k,l)
             enddo
             do ij = 1, ADM_gall_pl
                rhogvze_pl(ij,k,l) = rhogvz_pl(ij,k,l) * e_pl(ij,k,l)
             enddo
          enddo
       enddo
    endif

    !--- fez = rhow * e_w ( at half level ).
    !------ if 'HORIZONTAL', skip.
!    if ( f_type /= 'HORIZONTAL' ) then
    if ( fluxtype /= I_SRC_horizontal ) then ! [mod] H.Yashiro 20120530
       do l = 1, ADM_lall
          do k =ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                fez(g,k,l) = 0.5D0 * ( GRD_afac(k) * e(g,k,  l) &
                                     + GRD_bfac(k) * e(g,k-1,l) ) * rhogw(g,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   fez_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * e_pl(g,k,  l) &
                                           + GRD_bfac(k) * e_pl(g,k-1,l) ) * rhogw_pl(g,k,l)
                enddo
             enddo
          enddo
       endif

    else !--- = 'HORIZONTAL'

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                fez(g,k,l) = 0.D0
             enddo
          enddo
       enddo
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                fez_pl(g,k,l)=0.D0
             enddo
          enddo
       enddo

    endif

    !--- flux convergence step
    call src_flux_convergence( rhogvxe, rhogvxe_pl, & !--- [IN]
                               rhogvye, rhogvye_pl, & !--- [IN]
                               rhogvze, rhogvze_pl, & !--- [IN]
                               fez,     fez_pl,     & !--- [IN]
                               grhoge,  grhoge_pl,  & !--- [OUT]
!                               f_type=f_type )        ! [mod] H.Yashiro 20120530
                               fluxtype             ) !--- [IN] 

    return
  end subroutine src_advection_convergence


  subroutine src_advection_convergence_h(&  
       rhogvx, rhogvx_pl,              & !--- IN  : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy, rhogvy_pl,              & !--- IN  : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz, rhogvz_pl,              & !--- IN  : rho*Vz  ( gam2 X G^{1/2} )
       rhogw,  rhogw_pl,               & !--- IN  : rho*w   ( gam2 X G^{1/2} )
       e, e_pl,                        & !--- IN  : arbitrary scalar
       grhoge, grhoge_pl,              & !--- OUT : tendency
!       f_type                          & !--- IN  : type of flux conv.
       fluxtype                        ) !--- IN : flux type [mod] H.Yashiro 20120530
    !------ 
    !------ Advection convergence
    !------     1. f_type is for sub[flux_convergence]
    !------     2. vint_type is only for vertical direction.
    !------        (current version)
    !------     3. If e_w is present, it is used in the vertical flux
    !------        calculation.
    !------     4. Calculation region of grhoe, grhoe_pl
    !------                    : (:,ADM_kmin:ADM_kmax,:)
    !------ 
    !
    use mod_adm, only :  &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_GALL_PL,    &
         ADM_LALL_PL,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_grd, only :  &
         GRD_cfac,       &
         GRD_dfac
    !
    implicit none
    !
    real(8), intent(in) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: e(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: e_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(out) :: grhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhoge_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

!    character(LEN=*),intent(in) :: f_type
    !<-----                    = 'HORIZONTAL' : horizontal convergence
    !<-----                    = 'VERTICAL'   : vertical convergence
    !<-----                    = 'DEFAULT'    : both of them
    integer, intent(in) :: fluxtype ! [mod] H.Yashiro 20120530
    !<-----                    = I_SRC_horizontal : horizontal convergence
    !<-----                    = I_SRC_vertical   : vertical convergence
    !<-----                    = I_SRC_default    : both of them

    real(8) :: fez(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: fez_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8) :: rhogvxe(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvxe_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhogvye(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvye_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rhogvze(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhogvze_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    integer :: k
    !
    !
    rhogvxe(:,:,:)=rhogvx(:,:,:)*e(:,:,:)
    rhogvye(:,:,:)=rhogvy(:,:,:)*e(:,:,:)
    rhogvze(:,:,:)=rhogvz(:,:,:)*e(:,:,:)
    !
    if(ADM_prc_me==ADM_prc_pl) then
       rhogvxe_pl(:,:,:)=rhogvx_pl(:,:,:)*e_pl(:,:,:)
       rhogvye_pl(:,:,:)=rhogvy_pl(:,:,:)*e_pl(:,:,:)
       rhogvze_pl(:,:,:)=rhogvz_pl(:,:,:)*e_pl(:,:,:)
    end if
    !
    !--- fez = rhow * e_w ( at half level ).
    !------ if 'HORIZONTAL', skip.
    fez = 0.0D0
    fez_pl = 0.0D0
!    if ( f_type /= 'HORIZONTAL' ) then
    if ( fluxtype /= I_SRC_horizontal ) then ! [mod] H.Yashiro 20120530
       do k =ADM_kmin,ADM_kmax
          fez(:,k,:) = 0.5D0                  &
               * ( GRD_cfac(k)  * e(:,k+1,:)  &
                 + GRD_dfac(k)  * e(:,k  ,:) )&
               * rhogw(:,k,:)
       end do
       if(ADM_prc_me==ADM_prc_pl) then
          do k =ADM_kmin,ADM_kmax
             fez_pl(:,k,:) = 0.5D0                  &
                  * ( GRD_cfac(k)  * e_pl(:,k+1,:)    &
                    + GRD_cfac(k)  * e_pl(:,k  ,:) )&
                  * rhogw_pl(:,k,:)
          end do
       end if
    else !--- = 'HORIZONTAL'
       fez(:,:,:)=0.0D0
       fez_pl(:,:,:)=0.0D0
    end if
    !
    !--- flux convergence step
    call src_flux_convergence_h(&
         rhogvxe, rhogvxe_pl,  & !--- in
         rhogvye, rhogvye_pl,  & !--- in
         rhogvze, rhogvze_pl,  & !--- in
         fez,    fez_pl,       & !--- in
         grhoge,  grhoge_pl,   & !--- out
!         f_type=f_type )         !--- in
         fluxtype             ) !--- in [mod] H.Yashiro 20120530

  end subroutine src_advection_convergence_h
  !-----------------------------------------------------------------------------
  subroutine src_pres_work_uv(&
       vx, vx_pl,             &  !--- IN : Vx
       vy, vy_pl,             &  !--- IN : Vy
       vz, vz_pl,             &  !--- IN : Vz
       pgx, pgx_pl,           &  !--- IN : pres. grad. force ( x comp. )
       pgy, pgy_pl,           &  !--- IN : pres. grad. force ( y comp. )
       pgz, pgz_pl,           &  !--- IN : pres. grad. force ( z comp. )
       grhoge, grhoge_pl)        !--- OUT : pressure work
    !------ 
    !------ Pressure work term (horizontal)
    !------     1. vx*pgx + vy*pgy + vz*pgz
    !------     2.    pgx = dp/dx + d (G3X p )/d(gz)
    !------           pgy = dp/dy + d (G3Y p )/d(gz)
    !------           pgz = dp/dz + d (G3Z p )/d(gz)
    !------     3. Calculation region of grhoge, grhoge_pl
    !------                    : (:,:,:)
    !------ 
    !
    use mod_adm, only :  &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_GALL_PL,    &
         ADM_LALL_PL,    &
         ADM_prc_me,     &
         ADM_prc_pl
    !
    implicit none
    !
    real(8), intent(in) :: vx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: pgx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: pgx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: pgy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: pgy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: pgz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: pgz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(out) :: grhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhoge_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    integer :: g,k,l





    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do g = 1, ADM_gall
             grhoge(g,k,l)                &
                  = vx(g,k,l) * pgx(g,k,l)&
                  + vy(g,k,l) * pgy(g,k,l)&
                  + vz(g,k,l) * pgz(g,k,l)
          end do
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                grhoge_pl(g,k,l)                   &
                     = vx_pl(g,k,l) * pgx_pl(g,k,l)&
                     + vy_pl(g,k,l) * pgy_pl(g,k,l)&
                     + vz_pl(g,k,l) * pgz_pl(g,k,l)
             end do
          end do
       end do
    end if





  end subroutine src_pres_work_uv
  !-----------------------------------------------------------------------------
  subroutine src_advection_convergence_v( &  
       vx, vx_pl,                         & !--- IN  : Vx
       vy, vy_pl,                         & !--- IN  : Vy
       vz, vz_pl,                         & !--- IN  : Vz
       w,  w_pl,                          & !--- IN  : w 
       rhog, rhog_pl,                     & !--- IN  : rho( X G^{1/2} )
       rhogvx, rhogvx_pl,                 & !--- IN  : rho*Vx ( X G^{1/2} )
       rhogvy, rhogvy_pl,                 & !--- IN  : rho*Vy ( X G^{1/2} )
       rhogvz, rhogvz_pl,                 & !--- IN  : rho*Vz ( X G^{1/2} )
       rhogw,  rhogw_pl,                  & !--- IN  : rho*w  ( X G^{1/2} )
       grhogvx, grhogvx_pl,               & !--- OUT : tendency  of rhovx
       grhogvy, grhogvy_pl,               & !--- OUT : tendency  of rhovy
       grhogvz, grhogvz_pl,               & !--- OUT : tendency  of rhovz
       grhogw,  grhogw_pl                 & !--- OUT : tendency  of rhovw
       )
    !------ 
    !------ Advection convergence for velocity
    !------ 
    !
    use mod_adm, only :  &
         ADM_KNONE,      &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_GALL_PL,    &
         ADM_LALL_PL,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_grd, only :  &
         GRD_x,          &
         GRD_x_pl,       &
         GRD_XDIR,       &
         GRD_YDIR,       &
         GRD_ZDIR,       &
         GRD_rscale,     &
         GRD_afac,       &
         GRD_bfac,       &
         GRD_cfac,       &
         GRD_dfac
    use mod_cnst, only : &
         CNST_EOHM
    !
    use mod_runconf, only : &
         NON_HYDRO_ALPHA
    !
    implicit none
    !
    real(8), intent(in) :: vx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: w(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: w_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(in) :: rhog(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: rhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(out) :: grhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: grhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: grhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: grhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: grhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8) :: vvx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: vvy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: vvz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: dvvx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: dvvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: dvvy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: dvvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: dvvz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: dvvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8) :: grhogwc(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: grhogwc_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8) :: f_zcom(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: f_zcom_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    integer :: g,k,l









    do l = 1, ADM_lall
       do k=ADM_kmin,ADM_kmax
          do g = 1, ADM_gall
             vvx(g,k,l) = vx(g,k,l) &
                  + ( GRD_cfac(k) * w(g,k+1,l) + GRD_dfac(k) * w(g,k,l) )  * 0.5D0&
                  * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
          end do
          do g = 1, ADM_gall
             vvy(g,k,l) = vy(g,k,l) &
                  + ( GRD_cfac(k) * w(g,k+1,l) + GRD_dfac(k) * w(g,k,l) )  * 0.5D0&
                  * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
          end do
          do g = 1, ADM_gall
             vvz(g,k,l) = vz(g,k,l) &
                  + ( GRD_cfac(k) * w(g,k+1,l) + GRD_dfac(k) * w(g,k,l) )  * 0.5D0&
                  * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
          end do
       end do
    end do
    do l = 1, ADM_lall
       do g = 1, ADM_gall
          vvx(g,ADM_kmin-1,l) = 0.0D0
       end do
       do g = 1, ADM_gall
          vvx(g,ADM_kmax+1,l) = 0.0D0
       end do
       do g = 1, ADM_gall
          vvy(g,ADM_kmin-1,l) = 0.0D0
       end do
       do g = 1, ADM_gall
          vvy(g,ADM_kmax+1,l) = 0.0D0
       end do
       do g = 1, ADM_gall
          vvz(g,ADM_kmin-1,l) = 0.0D0
       end do
       do g = 1, ADM_gall
          vvz(g,ADM_kmax+1,l) = 0.0D0
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do k=ADM_kmin,ADM_kmax
             do g = 1, ADM_gall_pl
                vvx_pl(g,k,l) = vx_pl(g,k,l) &
                     + ( GRD_cfac(k) * w_pl(g,k+1,l) + GRD_dfac(k) * w_pl(g,k,l) )  * 0.5D0&
                     * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
             end do
             do g = 1, ADM_gall_pl
                vvy_pl(g,k,l) = vy_pl(g,k,l) &
                     + ( GRD_cfac(k) * w_pl(g,k+1,l) + GRD_dfac(k) * w_pl(g,k,l) )  * 0.5D0&
                     * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
             end do
             do g = 1, ADM_gall_pl
                vvz_pl(g,k,l) = vz_pl(g,k,l) &
                     + ( GRD_cfac(k) * w_pl(g,k+1,l) + GRD_dfac(k) * w_pl(g,k,l) )  * 0.5D0&
                     * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
             end do
          end do
       end do
       do l = 1, ADM_lall_pl
          do g = 1, ADM_gall_pl
             vvx_pl(g,ADM_kmin-1,l) = 0.0D0
          end do
          do g = 1, ADM_gall_pl
             vvx_pl(g,ADM_kmax+1,l) = 0.0D0
          end do
          do g = 1, ADM_gall_pl
             vvy_pl(g,ADM_kmin-1,l) = 0.0D0
          end do
          do g = 1, ADM_gall_pl
             vvy_pl(g,ADM_kmax+1,l) = 0.0D0
          end do
          do g = 1, ADM_gall_pl
             vvz_pl(g,ADM_kmin-1,l) = 0.0D0
          end do
          do g = 1, ADM_gall_pl
             vvz_pl(g,ADM_kmax+1,l) = 0.0D0
          end do
       end do
    end if
    !
    !
    !
    call src_advection_convergence(&  
         rhogvx, rhogvx_pl,              & !--- in
         rhogvy, rhogvy_pl,              & !--- in
         rhogvz, rhogvz_pl,              & !--- in
         rhogw,  rhogw_pl,               & !--- in
         vvx, vvx_pl,                  & !--- in
         dvvx, dvvx_pl,                & !--- out
!         f_type = 'DEFAULT'            ) !--- in
         I_SRC_default                 ) !--- in [mod] H.Yashiro 20120530
    call src_advection_convergence(&  
         rhogvx, rhogvx_pl,              & !--- in
         rhogvy, rhogvy_pl,              & !--- in
         rhogvz, rhogvz_pl,              & !--- in
         rhogw,  rhogw_pl,               & !--- in
         vvy, vvy_pl,                  & !--- in
         dvvy, dvvy_pl,                & !--- out
!         f_type = 'DEFAULT'            ) !--- in
         I_SRC_default                 ) !--- in [mod] H.Yashiro 20120530
    call src_advection_convergence(&  
         rhogvx, rhogvx_pl,              & !--- in
         rhogvy, rhogvy_pl,              & !--- in
         rhogvz, rhogvz_pl,              & !--- in
         rhogw,  rhogw_pl,               & !--- in
         vvz, vvz_pl,                  & !--- in
         dvvz, dvvz_pl,                & !--- out
!         f_type = 'DEFAULT'            ) !--- in
         I_SRC_default                 ) !--- in [mod] H.Yashiro 20120530

!    Omega(GRD_XDIR) = 0.0D0
!    Omega(GRD_YDIR) = 0.0D0
!    Omega(GRD_ZDIR) = CNST_EOHM

    do l = 1, ADM_lall
       do k=ADM_kmin,ADM_kmax
          do g = 1, ADM_gall
             dvvx(g,k,l) = dvvx(g,k,l) &
                  - 2.0D0*rhog(g,k,l) * ( -CNST_EOHM * vvy(g,k,l) ) !( Omega(GRD_YDIR)*vvz(g,k,l) - Omega(GRD_ZDIR)*vvy(g,k,l) )
             dvvy(g,k,l) = dvvy(g,k,l) &
                  - 2.0D0*rhog(g,k,l) * (  CNST_EOHM * vvx(g,k,l) ) !( Omega(GRD_ZDIR)*vvx(g,k,l) - Omega(GRD_XDIR)*vvz(g,k,l) )
             !dvvz(g,k,l) = dvvz(g,k,l) &
             !     - 2.0D0*rhog(g,k,l) * ( Omega(GRD_XDIR)*vvy(g,k,l) - Omega(GRD_YDIR)*vvx(g,k,l) )
             f_zcom(g,k,l) = dvvx(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale  &
                           + dvvy(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale  &
                           + dvvz(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
             grhogvx(g,k,l) = dvvx(g,k,l) &
                  - f_zcom(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
             grhogvy(g,k,l) = dvvy(g,k,l) &
                  - f_zcom(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
             grhogvz(g,k,l) = dvvz(g,k,l) &
                  - f_zcom(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
             grhogwc(g,k,l) = f_zcom(g,k,l) * dble(NON_HYDRO_ALPHA)
          end do
       end do
    end do
    !
    do l = 1, ADM_lall
       do g = 1, ADM_gall
          grhogvx(g,ADM_kmin-1,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogvy(g,ADM_kmin-1,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogvz(g,ADM_kmin-1,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogvx(g,ADM_kmax+1,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogvy(g,ADM_kmax+1,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogvz(g,ADM_kmax+1,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogw(g,ADM_kmin-1,l) = 0.0D0
       end do
       do g = 1, ADM_gall
          grhogw(g,ADM_kmin  ,l)=0.0D0
       end do
       do g = 1, ADM_gall
          grhogw(g,ADM_kmax+1,l)=0.0D0
       end do
    end do
    !
    do l = 1, ADM_lall
       do k=ADM_kmin+1,ADM_kmax
          do g = 1, ADM_gall
             grhogw(g,k,l) &
                  = ( GRD_afac(k) * grhogwc(g,k,l)  &
                    + GRD_bfac(k) * grhogwc(g,k-1,l) ) * 0.5D0
          end do
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do k=ADM_kmin,ADM_kmax
             do g = 1, ADM_gall_pl
                dvvx_pl(g,k,l) = dvvx_pl(g,k,l) &
                     - 2.0D0*rhog_pl(g,k,l) * ( -CNST_EOHM * vvy_pl(g,k,l) ) !( Omega(GRD_YDIR)*vvz_pl(g,k,l) - Omega(GRD_ZDIR)*vvy_pl(g,k,l) )
                dvvy_pl(g,k,l) = dvvy_pl(g,k,l) &
                     - 2.0D0*rhog_pl(g,k,l) * (  CNST_EOHM * vvx_pl(g,k,l) ) !( Omega(GRD_ZDIR)*vvx_pl(g,k,l) - Omega(GRD_XDIR)*vvz_pl(g,k,l) )
                !dvvz_pl(g,k,l) = dvvz_pl(g,k,l) &
                !     - 2.0D0*rhog_pl(g,k,l) * ( Omega(GRD_XDIR)*vvy_pl(g,k,l) - Omega(GRD_YDIR)*vvx_pl(g,k,l) )
                f_zcom_pl(g,k,l) = dvvx_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale  &
                                 + dvvy_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale  &
                                 + dvvz_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
                grhogvx_pl(g,k,l) = dvvx_pl(g,k,l) &
                     - f_zcom_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
                grhogvy_pl(g,k,l) = dvvy_pl(g,k,l) &
                     - f_zcom_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
                grhogvz_pl(g,k,l) = dvvz_pl(g,k,l) &
                     - f_zcom_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
                grhogwc_pl(g,k,l) = f_zcom_pl(g,k,l) * dble(NON_HYDRO_ALPHA)
             end do
          end do
       end do
       do l = 1, ADM_lall_pl
          do g = 1, ADM_gall_pl
             grhogvx_pl(g,ADM_kmin-1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogvy_pl(g,ADM_kmin-1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogvz_pl(g,ADM_kmin-1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogvx_pl(g,ADM_kmax+1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogvy_pl(g,ADM_kmax+1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogvz_pl(g,ADM_kmax+1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogw_pl(g,ADM_kmin-1,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogw_pl(g,ADM_kmin  ,l)=0.0D0
          end do
          do g = 1, ADM_gall_pl
             grhogw_pl(g,ADM_kmax+1,l)=0.0D0
          end do
       end do
       do l = 1, ADM_lall_pl
          do k=ADM_kmin+1,ADM_kmax
             do g = 1, ADM_gall_pl
                grhogw_pl(g,k,l) &
                      = ( GRD_afac(k) * grhogwc_pl(g,k,l)  &
                        + GRD_bfac(k) * grhogwc_pl(g,k-1,l) ) * 0.5D0
             end do
          end do
       end do
    end if





  end subroutine src_advection_convergence_v
  !-----------------------------------------------------------------------------  
  subroutine src_gradient( &
       e, e_pl,            &  !--- IN  : phi * G^{1/2} * gamma^2
       gex, gex_pl,        &  !--- OUT : horizontal grad. ( x-comp. )
       gey, gey_pl,        &  !--- OUT : horizontal grad. ( y-comp. )
       gez, gez_pl,        &  !--- OUT : horizontal grad. ( z-comp. )
       gevz, gevz_pl,      &  !--- OUT : vertical grad.   ( half point )
!       f_type )               !--- IN  : calculation type
       grad_type            ) !--- IN : calculation type [mod] H.Yashiro 20120530
    !------ 
    !------ Gradient operator for 3D
    !------     1. Calculation region of gex, gey, gez
    !------                    : (:, ADM_kmin:ADM_kmax,:)
    !------     2. Calculation region of gevz
    !------                    : (:, ADM_kmin:ADM_kmax+1,:)
    !------     3. Return values are
    !------ 
    !
    use mod_adm, only :  &
         ADM_VMISS,      &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_GALL_PL,    &
         ADM_LALL_PL,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_vmtr, only : &
         VMTR_RGSGAM2,   &
         VMTR_RGSGAM2_pl,&
         VMTR_GZXH,      &
         VMTR_GZXH_pl,   &
         VMTR_GZYH,      &
         VMTR_GZYH_pl,   &
         VMTR_GZZH,      &
         VMTR_GZZH_pl,   &
         VMTR_GSGAMH,    &
         VMTR_GSGAMH_pl, &
         VMTR_RGAM,      &
         VMTR_RGAM_pl,   &
         VMTR_GAM2H,     &
         VMTR_GAM2H_pl
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac,       &
         GRD_rdgz,       &
         GRD_rdgzh
    use mod_oprt, only :        &
         OPRT_horizontalize_vec,&
         OPRT_gradient
    implicit none

    real(8), intent(in) :: e(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: e_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8), intent(out) :: gex(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gex_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gey(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gey_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gez(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gez_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gevz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gevz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8) :: fex_h(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: fex_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: fey_h(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: fey_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: fez_h(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: fez_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
!    character(LEN=*),intent(in) :: f_type
    !<-----            HORIZONTAL          : horizontal gradient
    !<-----            VERTICAL            : vertical gradient ( do not used! )
    !<-----            DEFAULT             : both of them
    integer, intent(in) :: grad_type ! [mod] H.Yashiro 20120530
    !<-----                    = I_SRC_horizontal : horizontal gradient
    !<-----                    = I_SRC_vertical   : vertical gradient ( no use )
    !<-----                    = I_SRC_default    : both of them

    !
    real(8) :: gee(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: gee_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    integer :: k,n,l






    !------ horizontal gradient without mountain
    do l=1,ADM_lall
       do k=1,ADM_kall
          do n=1,ADM_gall
             gex(n,k,l) = ADM_VMISS
          end do
          do n=1,ADM_gall
             gey(n,k,l) = ADM_VMISS
          end do
          do n=1,ADM_gall
             gez(n,k,l) = ADM_VMISS
          end do
       end do
    end do

    do l=1,ADM_lall_pl
       do k=1,ADM_kall
          do n=1,ADM_gall_pl
             gex_pl(n,k,l) = ADM_VMISS
          end do
          do n=1,ADM_gall_pl
             gey_pl(n,k,l) = ADM_VMISS
          end do
          do n=1,ADM_gall_pl
             gez_pl(n,k,l) = ADM_VMISS
          end do
       end do
    end do
    !
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=1,ADM_gall
          gee(n,k,l) = e(n,k,l)*VMTR_RGAM(n,k,l)
        end do
      end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
      do l=1,ADM_lall_pl
        do k=1,ADM_kall
          do n=1,ADM_gall_pl
            gee_pl(n,k,l) = e_pl(n,k,l) *VMTR_RGAM_pl(n,k,l) 
          end do
        end do
      end do
    end if






    call OPRT_gradient(&
         gex, gex_pl,  &  !--- out
         gey, gey_pl,  &  !--- out
         gez, gez_pl,  &  !--- out
         gee, gee_pl )    !--- in
    




    !--- horizontal gradient
    do l = 1, ADM_lall





       do k = ADM_kmin, ADM_kmax+1
          do n=1, ADM_gall
             fex_h(n,k,l)                                                  &
                  = VMTR_GZXH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          end do
          do n=1, ADM_gall
             fey_h(n,k,l)                                                  &
                  = VMTR_GZYH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          end do
          do n=1, ADM_gall
             fez_h(n,k,l)                                                  &
                  = VMTR_GZZH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          end do
       end do
       !




       !
       do k = ADM_kmin, ADM_kmax
          do n=1, ADM_gall
             gex(n,k,l) = gex(n,k,l) &
                  + (  fex_h(n,k+1,l) - fex_h(n,k,l) ) * GRD_rdgz(k)
             gey(n,k,l) = gey(n,k,l) &
                  + (  fey_h(n,k+1,l) - fey_h(n,k,l) ) * GRD_rdgz(k)
             gez(n,k,l) = gez(n,k,l) &
                  + (  fez_h(n,k+1,l) - fez_h(n,k,l) ) * GRD_rdgz(k)
          end do
       end do
       !




       !
       !--- At the lowest layer, do not use the extrapolation value!
       if(first_layer_remedy) then
         do n=1,ADM_gall
           gex(n,ADM_kmin,l) = gex(n,ADM_kmin+1,l)
           gey(n,ADM_kmin,l) = gey(n,ADM_kmin+1,l)
           gez(n,ADM_kmin,l) = gez(n,ADM_kmin+1,l)
         end do
       end if





    end do
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl





          do k = ADM_kmin, ADM_kmax+1
             do n=1, ADM_gall_pl
                fex_h_pl(n,k,l) &
                     = VMTR_GZXH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)
             end do
             do n=1, ADM_gall_pl
                fey_h_pl(n,k,l)                                                     &
                     = VMTR_GZYH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)
             end do
             do n=1, ADM_gall_pl
                fez_h_pl(n,k,l)                                                     &
                     = VMTR_GZZH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)

             end do
          end do






          do k = ADM_kmin, ADM_kmax
             do n=1, ADM_gall_pl
                gex_pl(n,k,l) = gex_pl(n,k,l) &
                     + (  fex_h_pl(n,k+1,l) - fex_h_pl(n,k,l) ) * GRD_rdgz(k)
                gey_pl(n,k,l) = gey_pl(n,k,l) &
                     + (  fey_h_pl(n,k+1,l) - fey_h_pl(n,k,l) ) * GRD_rdgz(k)
                gez_pl(n,k,l) = gez_pl(n,k,l) &
                     + (  fez_h_pl(n,k+1,l) - fez_h_pl(n,k,l) ) * GRD_rdgz(k)
             end do
          end do
          !




          !
          !--- At the lowest layer, do not use the extrapolation value!
          if(first_layer_remedy) then
            do n=1,ADM_gall_pl
              gex_pl(n,ADM_kmin,l) = gex_pl(n,ADM_kmin+1,l)
              gey_pl(n,ADM_kmin,l) = gey_pl(n,ADM_kmin+1,l)
              gez_pl(n,ADM_kmin,l) = gez_pl(n,ADM_kmin+1,l)
            end do
          end if





       end do
    end if






    !--- horizontalize
    call OPRT_horizontalize_vec( &
          gex, gex_pl,           & !--- inout
          gey, gey_pl,           & !--- inout
          gez, gez_pl )            !--- inout





    !--- vertical gradient ( note : half points )
!    if ( f_type == 'HORIZONTAL' ) then
    if ( grad_type == I_SRC_horizontal ) then ! [mod] H.Yashiro 20120530





       do l=1,ADM_lall
         do k=1,ADM_kall
           do n=1,ADM_gall
             gevz(n,k,l) = 0.0D0
           end do
         end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
         do l=1,ADM_lall_pl
           do k=1,ADM_kall
             do n=1,ADM_gall_pl 
               gevz_pl(n,k,l) = 0.0D0
             end do
           end do
         end do
       end if





    else





       do l=1,ADM_lall
         do n=1, ADM_gall
           gevz(n,ADM_kmin-1,l) = ADM_VMISS
         end do
       end do
       do l=1,ADM_lall
         do k = ADM_kmin, ADM_kmax+1
           do n=1,ADM_gall
             gevz(n,k,l)                                &
               = ( e(n,k,l) * VMTR_RGSGAM2(n,k,l)       &
                -  e(n,k-1,l) * VMTR_RGSGAM2(n,k-1,l) ) &
                 * GRD_rdgzh(k) * VMTR_GAM2H(n,k,l)
           end do
         end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
            do n=1,ADM_gall_pl
              gevz_pl(n,ADM_kmin-1,l) = ADM_VMISS
            end do
            do k = ADM_kmin, ADM_kmax+1
              do n=1,ADM_gall_pl
                gevz_pl(n,k,l) &
                  = ( e_pl(n,k,l) * VMTR_RGSGAM2_pl(n,k,l)       &
                    - e_pl(n,k-1,l) * VMTR_RGSGAM2_pl(n,k-1,l) ) &
                    * GRD_rdgzh(k) * VMTR_GAM2H_pl(n,k,l)
              end do
            end do
          end do
       end if





    end if





  end subroutine src_gradient

  !----------------------------------------------------------------------------------
  ! Y.Niwa add 080124
  ! This routine is revised version of src_update_tracer
  subroutine src_update_tracer( &
       nqmax,                       & !--- IN    : number of tracers
       rhogq,       rhogq_pl,       & !--- INOUT : rhogq   ( gam2 X G^{1/2} )
       rhog_in,     rhog_in_pl,     & !--- IN    : rho(old)( gam2 X G^{1/2} )
       rhog_mean,   rhog_mean_pl,   & !--- IN    : rho     ( gam2 X G^{1/2} )
       rhogvx_mean, rhogvx_mean_pl, & !--- IN    : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy_mean, rhogvy_mean_pl, & !--- IN    : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz_mean, rhogvz_mean_pl, & !--- IN    : rho*Vz  ( gam2 X G^{1/2} )
       rhogw_mean,  rhogw_mean_pl,  & !--- IN    : rho*w   ( gam2 X G^{1/2} )
       frhog,       frhog_pl,       & !--- IN    : hyperviscosity tendency for rhog
       dt,                          & !--- IN    : delta t
       thubern_lim                  ) !--- IN    : switch of thubern limiter
    use mod_adm, only :  &
         ADM_gall,       &
         ADM_gmin,       &
         ADM_gmax,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl,     &
         ADM_AI,ADM_AJ,  &
         ADM_GSLF_PL,    &
         ADM_GMIN_PL,    &
         ADM_GMAX_PL,    &
         ADM_gall_1d
    use mod_vmtr, only : &
         VMTR_RGSGAM2,   &
         VMTR_RGSGAM2_pl,&
         VMTR_GZXH,      &
         VMTR_GZXH_pl,   &
         VMTR_GZYH,      &
         VMTR_GZYH_pl,   &
         VMTR_GZZH,      &
         VMTR_GZZH_pl,   &
         VMTR_GSGAMH,    &
         VMTR_GSGAMH_pl, &
         VMTR_RGSH,      &
         VMTR_RGSH_pl,   &
         VMTR_RGAM,      &
         VMTR_RGAM_pl
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_XDIR, &
       GRD_ZDIR, &
       GRD_rdgz
    use mod_oprt, only: &
       oprt_divergence2_prep, &
       oprt_divergence2
    use mod_advlim_thuburn, only: &
       advlim_thuburn_v
    implicit none

    integer, intent(in)    :: nqmax
    real(8), intent(inout) :: rhogq         (ADM_gall,   ADM_kall,ADM_lall,   nqmax)
    real(8), intent(inout) :: rhogq_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,nqmax)
    real(8), intent(in)    :: rhog_in       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhog_in_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhog_mean     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhog_mean_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogvx_mean   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogvx_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogvy_mean   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogvy_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogvz_mean   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogvz_mean_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: rhogw_mean    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rhogw_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: frhog         (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: frhog_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)    :: dt

    logical, intent(in)    :: thubern_lim  ![add] 20130613 R.Yoshida

    real(8) :: rhog    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rrhog   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rrhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: q       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: q_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: q_h     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: q_h_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: d       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: d_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: flx_v   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: ck      (ADM_gall,   ADM_kall,ADM_lall,   2)
    real(8) :: ck_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(8) :: flx_h   (6,ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_h_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: c       (6,ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: c_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: vx_r    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vx_r_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vy_r    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vy_r_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vz_r    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vz_r_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: hdiv    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: hdiv_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: cp      (ADM_gall,   ADM_kall,ADM_lall,   ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8) :: cp_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR)

    real(8), parameter :: b1 = 0.D0
    real(8), parameter :: b2 = 1.D0
    real(8), parameter :: b3 = 1.D0 - (b1+b2)

    integer :: nstart, nend
    integer :: g, k, l, n, nq

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 1st
    !---------------------------------------------------------------------------
    do l = 1, ADM_lall

       do k = 1, ADM_kall
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog_in(g,k,l)
             d(g,k,l) = b1 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall
             flx_v(g,k,l) = ( ( ( GRD_afac(k) * VMTR_RGSGAM2(g,k  ,l) * rhogvx_mean(g,k  ,l) &
                                + GRD_bfac(k) * VMTR_RGSGAM2(g,k-1,l) * rhogvx_mean(g,k-1,l) &
                                ) * 0.5D0 * VMTR_GSGAMH(g,k,l) * VMTR_GZXH(g,k,l)            &
                              + ( GRD_afac(k) * VMTR_RGSGAM2(g,k  ,l) * rhogvy_mean(g,k  ,l) &
                                + GRD_bfac(k) * VMTR_RGSGAM2(g,k-1,l) * rhogvy_mean(g,k-1,l) &
                                ) * 0.5D0 * VMTR_GSGAMH(g,k,l) * VMTR_GZYH(g,k,l)            & 
                              + ( GRD_afac(k) * VMTR_RGSGAM2(g,k  ,l) * rhogvz_mean(g,k  ,l) &
                                + GRD_bfac(k) * VMTR_RGSGAM2(g,k-1,l) * rhogvz_mean(g,k-1,l) &
                                ) * 0.5D0 * VMTR_GSGAMH(g,k,l) * VMTR_GZZH(g,k,l)            &
                              ) + rhogw_mean(g,k,l) * VMTR_RGSH(g,k,l) &
                            ) * (0.5D0*dt)
          enddo
       enddo

       do g = 1, ADM_gall
          flx_v(g,ADM_kmin,  l) = 0.D0
          flx_v(g,ADM_kmax+1,l) = 0.D0
       enddo

       !---- Courant numbers at cell boundary
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             ck(g,k,l,1) = -flx_v(g,k  ,l) * rrhog(g,k,l) * GRD_rdgz(k)
             ck(g,k,l,2) =  flx_v(g,k+1,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
       enddo

       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,1) = 0.D0 
          ck(g,ADM_kmin-1,l,2) = 0.D0
          ck(g,ADM_kmax+1,l,1) = 0.D0
          ck(g,ADM_kmax+1,l,2) = 0.D0
       enddo

    enddo ! l LOOP

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_in_pl(g,k,l)
                d_pl(g,k,l) = b1 * frhog_pl(g,k,l) * dt * rrhog_pl(g,k,l)
             enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
             do g = 1, ADM_gall_pl
                flx_v_pl(g,k,l) = ( ( ( GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * rhogvx_mean_pl(g,k  ,l) &
                                      + GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * rhogvx_mean_pl(g,k-1,l) &
                                      ) * 0.5D0 * VMTR_GSGAMH_pl(g,k,l) * VMTR_GZXH_pl(g,k,l)            &
                                    + ( GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * rhogvy_mean_pl(g,k  ,l) &
                                      + GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * rhogvy_mean_pl(g,k-1,l) &
                                      ) * 0.5D0 * VMTR_GSGAMH_pl(g,k,l) * VMTR_GZYH_pl(g,k,l)            & 
                                    + ( GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * rhogvz_mean_pl(g,k  ,l) &
                                      + GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * rhogvz_mean_pl(g,k-1,l) &
                                      ) * 0.5D0 * VMTR_GSGAMH_pl(g,k,l) * VMTR_GZZH_pl(g,k,l)            &
                                    ) + rhogw_mean_pl(g,k,l) * VMTR_RGSH_pl(g,k,l)                       &
                                  ) * (0.5d0*dt)
             enddo
          enddo

          do g = 1, ADM_gall_pl
             flx_v_pl(g,ADM_kmin,  l) = 0.D0
             flx_v_pl(g,ADM_kmax+1,l) = 0.D0
          enddo

          !---- Courant numbers at cell boundary
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,1) = -flx_v_pl(g,k  ,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
                ck_pl(g,k,l,2) =  flx_v_pl(g,k+1,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
          enddo

          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,1) = 0.D0
             ck_pl(g,ADM_kmin-1,l,2) = 0.D0
             ck_pl(g,ADM_kmax+1,l,1) = 0.D0
             ck_pl(g,ADM_kmax+1,l,2) = 0.D0
          enddo

       enddo
    endif

    !--- basic scheme ( 2nd-order centered difference )
    do nq = 1, nqmax

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                q_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * q(g,k,  l) &
                                     + GRD_bfac(k) * q(g,k-1,l) )
             enddo
          enddo

          do g = 1, ADM_gall
             q_h(g,ADM_kmin-1,l) = 0.D0
          enddo
       enddo


       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo

             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   q_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(g,k,  l) &
                                           + GRD_bfac(k) * q_pl(g,k-1,l) )
                enddo
             enddo

             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin-1,l) = 0.D0
             enddo
          enddo
       endif

       ! [mod] 20130613 R.Yoshida
       if (thubern_lim) call advlim_thuburn_v( q_h, q_h_pl, & !--- [INOUT]
                                             q,   q_pl,   & !--- [IN]
                                             ck,  ck_pl,  & !--- [IN]
                                             d,   d_pl    ) !--- [IN]

       !--- update rhogq
       do l = 1, ADM_lall
          do g = 1, ADM_gall
             q_h(g,ADM_kmin  ,l) = 0.D0
             q_h(g,ADM_kmax+1,l) = 0.D0
          enddo

          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) &
                                - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                  - flx_v(g,k,  l) * q_h(g,k,  l) ) * GRD_rdgz(k)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin  ,l) = 0.D0
                q_h_pl(g,ADM_kmax+1,l) = 0.D0
             enddo

             do k = ADM_kmin, ADM_kmax
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) &
                                      - ( flx_v_pl(g,k+1,l) * q_h_pl(g,k+1,l) &
                                        - flx_v_pl(g,k,  l) * q_h_pl(g,k,  l) ) * GRD_rdgz(k)
                enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             rhog(g,k,l) = rhog_in(g,k,l)                                  &
                         - ( flx_v(g,k+1,l) - flx_v(g,k,l) ) * GRD_rdgz(k) &
                         + b1 * frhog(g,k,l) * dt
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo

       do g = 1, ADM_gall
          rhog(g,ADM_kmin-1,l) = rhog_in(g,ADM_kmin,l)
          rhog(g,ADM_kmax+1,l) = rhog_in(g,ADM_kmax,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                rhog_pl(g,k,l) = rhog_in_pl(g,k,l)                                     &
                               - ( flx_v_pl(g,k+1,l) - flx_v_pl(g,k,l) ) * GRD_rdgz(k) &
                               + b1 * frhog_pl(g,k,l) * dt
             enddo
          enddo

          do g = 1, ADM_gall_pl
             rhog_pl(g,ADM_kmin-1,l) = rhog_in_pl(g,ADM_kmin,l)
             rhog_pl(g,ADM_kmax+1,l) = rhog_in_pl(g,ADM_kmax,l)
          enddo
       enddo
    endif


    !---------------------------------------------------------------------------
    ! Horizontal advection by MIURA 2004 scheme
    !---------------------------------------------------------------------------
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do g = 1, ADM_gall
             vx_r(g,k,l) = rhogvx_mean(g,k,l) * VMTR_RGAM(g,k,l)
             vy_r(g,k,l) = rhogvy_mean(g,k,l) * VMTR_RGAM(g,k,l) 
             vz_r(g,k,l) = rhogvz_mean(g,k,l) * VMTR_RGAM(g,k,l) 
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                vx_r_pl(g,k,l) = rhogvx_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)  
                vy_r_pl(g,k,l) = rhogvy_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)   
                vz_r_pl(g,k,l) = rhogvz_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)   
             enddo
          enddo
       enddo
    endif

    call oprt_divergence2_prep( flx_h,     flx_h_pl,     & !--- [OUT]
                                cp,        cp_pl,        & !--- [OUT]
                                vx_r,      vx_r_pl,      & !--- [IN]
                                vy_r,      vy_r_pl,      & !--- [IN]
                                vz_r,      vz_r_pl,      & !--- [IN]
                                rhog_mean, rhog_mean_pl, & !--- [IN]
                                dt                       ) !--- [IN]

    do l = 1, ADM_lall
       do k = 1, ADM_kall

          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog(g,k,l)
             d(g,k,l) = b2 * frhog(g,k,l) * dt * rrhog(g,k,l)
             !--- Courant number
             c(1,g,k,l) = flx_h(1,g,k,l) * rrhog(g,k,l)
             c(2,g,k,l) = flx_h(2,g,k,l) * rrhog(g,k,l)
             c(3,g,k,l) = flx_h(3,g,k,l) * rrhog(g,k,l)
             c(4,g,k,l) = flx_h(4,g,k,l) * rrhog(g,k,l)
             c(5,g,k,l) = flx_h(5,g,k,l) * rrhog(g,k,l)
             c(6,g,k,l) = flx_h(6,g,k,l) * rrhog(g,k,l)
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall

             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_pl(g,k,l)
             enddo

             do n = ADM_GMIN_PL, ADM_GMAX_PL
                c_pl(n,k,l) = flx_h_pl(n,k,l) * rrhog_pl(ADM_GSLF_PL,k,l)
             enddo

             d_pl(ADM_GSLF_PL,k,l) = b2 * frhog_pl(ADM_GSLF_PL,k,l) * dt * rrhog_pl(ADM_GSLF_PL,k,l)
          enddo
       enddo
    endif

    do nq = 1, nqmax

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo
          enddo
       endif


       call oprt_divergence2( hdiv,  hdiv_pl,  & !--- [OUT]
                              q,     q_pl,     & !--- [IN]
                              flx_h, flx_h_pl, & !--- [IN]
                              c,     c_pl,     & !--- [IN]
                              cp,    cp_pl,    & !--- [IN]
                              d,     d_pl,     & !--- [IN]
                              dt               ) !--- [IN]

       !--- update rhogq
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) - hdiv(g,k,l) * dt
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) - hdiv_pl(g,k,l) * dt
                enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    !--- update rhog
    do l = 1, ADM_lall
       nstart = suf(ADM_gmin,ADM_gmin)
       nend   = suf(ADM_gmax,ADM_gmax)

       do k = 1, ADM_kall
          do n = nstart, nend
             rhog(n,k,l)= rhog(n,k,l) - ( flx_h(1,n,k,l) &
                                        + flx_h(2,n,k,l) &
                                        + flx_h(3,n,k,l) &
                                        + flx_h(4,n,k,l) &
                                        + flx_h(5,n,k,l) &
                                        + flx_h(6,n,k,l) &
                                        ) + b2 * frhog(n,k,l) * dt
          enddo
       enddo

       do k = 1, ADM_kall
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_GSLF_PL

       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             rhog_pl(n,k,l)= rhog_pl(n,k,l) - ( flx_h_pl(ADM_GMIN_PL  ,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+1,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+2,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+3,k,l) &
                                              + flx_h_pl(ADM_GMIN_PL+4,k,l) &
                                              ) + b2 * frhog_pl(n,k,l) * dt
          enddo
       enddo
    endif

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 2nd
    !---------------------------------------------------------------------------
    do l = 1, ADM_lall

       do k = 1, ADM_kall
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog(g,k,l)
          enddo
       enddo

       !---- Courant numbers at cell boundary
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             ck(g,k,l,1) = -flx_v(g,k  ,l) * rrhog(g,k,l) * GRD_rdgz(k)
             ck(g,k,l,2) =  flx_v(g,k+1,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
       enddo

       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,1) = 0.D0 
          ck(g,ADM_kmin-1,l,2) = 0.D0
          ck(g,ADM_kmax+1,l,1) = 0.D0
          ck(g,ADM_kmax+1,l,2) = 0.D0
       enddo

       do k = 1, ADM_kall
          do g = 1, ADM_gall
             d(g,k,l) = b3 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo
       enddo

    enddo ! l LOOP

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_pl(g,k,l)
             enddo
          enddo

          !---- Courant numbers at cell boundary
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,1) = -flx_v_pl(g,k  ,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
                ck_pl(g,k,l,2) =  flx_v_pl(g,k+1,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
          enddo

          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,1) = 0.D0
             ck_pl(g,ADM_kmin-1,l,2) = 0.D0
             ck_pl(g,ADM_kmax+1,l,1) = 0.D0
             ck_pl(g,ADM_kmax+1,l,2) = 0.D0
          enddo

          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                d_pl(g,k,l) = b3 * frhog_pl(g,k,l) * dt * rrhog_pl(g,k,l)
             enddo
          enddo

       enddo
    endif

    !--- basic scheme ( 2nd-order centered difference )
    do nq = 1, nqmax

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                q_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * q(g,k,  l) &
                                     + GRD_bfac(k) * q(g,k-1,l) )
             enddo
          enddo

          do g = 1, ADM_gall
             q_h(g,ADM_kmin-1,l) = 0.D0
          enddo
       enddo


       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo

             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   q_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(g,k,  l) &
                                           + GRD_bfac(k) * q_pl(g,k-1,l) )
                enddo
             enddo

             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin-1,l) = 0.D0
             enddo
          enddo
       endif

       ! [mod] 20130613 R.Yoshida
       if (thubern_lim) call advlim_thuburn_v( q_h, q_h_pl, & !--- [INOUT]
                                             q,   q_pl,   & !--- [IN]
                                             ck,  ck_pl,  & !--- [IN]
                                             d,   d_pl    ) !--- [IN]

       !--- update rhogq
       do l = 1, ADM_lall
          do g = 1, ADM_gall
             q_h(g,ADM_kmin  ,l) = 0.D0
             q_h(g,ADM_kmax+1,l) = 0.D0
          enddo

          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) &
                                - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                  - flx_v(g,k,  l) * q_h(g,k,  l) ) * GRD_rdgz(k)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin  ,l) = 0.D0
                q_h_pl(g,ADM_kmax+1,l) = 0.D0
             enddo

             do k = ADM_kmin, ADM_kmax
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) &
                                      - ( flx_v_pl(g,k+1,l) * q_h_pl(g,k+1,l) &
                                        - flx_v_pl(g,k,  l) * q_h_pl(g,k,  l) ) * GRD_rdgz(k)
                enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

    return
  end subroutine src_update_tracer

end module mod_src
