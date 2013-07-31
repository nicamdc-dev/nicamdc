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
  !++ Public procedure
  !
  public :: src_advection_convergence_momentum
  public :: src_advection_convergence
  public :: src_flux_convergence

  public :: src_buoyancy
  public :: src_pres_work_uv
  public :: src_gradient

  public :: src_update_tracer

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: I_SRC_horizontal = 1 ! [add] H.Yashiro 20120530
  integer, public, parameter :: I_SRC_vertical   = 2 ! [add] H.Yashiro 20120530
  integer, public, parameter :: I_SRC_default    = 3 ! [add] H.Yashiro 20120530

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, parameter :: first_layer_remedy = .true.

  real(8), private, allocatable, save :: AFACovGSGAM2   (:,:,:)
  real(8), private, allocatable, save :: AFACovGSGAM2_pl(:,:,:)
  real(8), private, allocatable, save :: BFACovGSGAM2   (:,:,:)
  real(8), private, allocatable, save :: BFACovGSGAM2_pl(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Advection convergence for momentum
  subroutine src_advection_convergence_momentum( &  
       vx,      vx_pl,      &
       vy,      vy_pl,      &
       vz,      vz_pl,      &
       w,       w_pl,       &
       rhog,    rhog_pl,    &
       rhogvx,  rhogvx_pl,  &
       rhogvy,  rhogvy_pl,  &
       rhogvz,  rhogvz_pl,  &
       rhogw,   rhogw_pl,   &
       grhogvx, grhogvx_pl, &
       grhogvy, grhogvy_pl, &
       grhogvz, grhogvz_pl, &
       grhogw,  grhogw_pl   )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_EOHM
    use mod_grd, only: &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR,   &
       GRD_rscale, &
       GRD_e,      &
       GRD_e_pl,   &
       GRD_x,      &
       GRD_x_pl,   &
       GRD_afac,   &
       GRD_bfac,   &
       GRD_cfac,   &
       GRD_dfac
    use mod_runconf, only : &
       NON_HYDRO_ALPHA
    implicit none

    real(8), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: w    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: w_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out) :: grhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: vvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: vvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: vvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: vvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
 
    real(8) :: dvvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dvvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dvvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dvvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dvvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dvvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: grhogwc   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhogwc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: prd, wc

    integer :: g, k, l
    !---------------------------------------------------------------------------

    !---< merge horizontal velocity & vertical velocity >

    do l = 1, ADM_lall
       do k = ADM_kmin,ADM_kmax
       do g = 1, ADM_gall
          wc = 0.5D0 * ( GRD_cfac(k) * w(g,k+1,l) &
                       + GRD_dfac(k) * w(g,k  ,l) )

          vvx(g,k,l) = vx(g,k,l) + wc * GRD_e(g,l,GRD_XDIR)
          vvy(g,k,l) = vy(g,k,l) + wc * GRD_e(g,l,GRD_YDIR)
          vvz(g,k,l) = vz(g,k,l) + wc * GRD_e(g,l,GRD_ZDIR)
       enddo
       enddo

       vvx(:,ADM_kmin-1,l) = 0.D0
       vvx(:,ADM_kmax+1,l) = 0.D0
       vvy(:,ADM_kmin-1,l) = 0.D0
       vvy(:,ADM_kmax+1,l) = 0.D0
       vvz(:,ADM_kmin-1,l) = 0.D0
       vvz(:,ADM_kmax+1,l) = 0.D0
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             wc = 0.5D0 * ( GRD_cfac(k) * w_pl(g,k+1,l) &
                          + GRD_dfac(k) * w_pl(g,k  ,l) )

             vvx_pl(g,k,l) = vx_pl(g,k,l) + wc * GRD_e_pl(g,l,GRD_XDIR)
             vvy_pl(g,k,l) = vy_pl(g,k,l) + wc * GRD_e_pl(g,l,GRD_YDIR)
             vvz_pl(g,k,l) = vz_pl(g,k,l) + wc * GRD_e_pl(g,l,GRD_ZDIR)
          enddo
          enddo

          vvx_pl(:,ADM_kmin-1,l) = 0.D0
          vvx_pl(:,ADM_kmax+1,l) = 0.D0
          vvy_pl(:,ADM_kmin-1,l) = 0.D0
          vvy_pl(:,ADM_kmax+1,l) = 0.D0
          vvz_pl(:,ADM_kmin-1,l) = 0.D0
          vvz_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    !---< advection term for momentum >

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvx,    vvx_pl,    & ! [IN]
                                    dvvx,   dvvx_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvy,    vvy_pl,    & ! [IN]
                                    dvvy,   dvvy_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvz,    vvz_pl,    & ! [IN]
                                    dvvz,   dvvz_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]


    do l = 1, ADM_lall
       !---< coriolis force >
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          dvvx(g,k,l) = dvvx(g,k,l) - 2.D0 * rhog(g,k,l) * ( -CNST_EOHM * vvy(g,k,l) )
          dvvy(g,k,l) = dvvy(g,k,l) - 2.D0 * rhog(g,k,l) * (  CNST_EOHM * vvx(g,k,l) )
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          prd = dvvx(g,k,l) * GRD_e(g,l,GRD_XDIR) &
              + dvvy(g,k,l) * GRD_e(g,l,GRD_YDIR) &
              + dvvz(g,k,l) * GRD_e(g,l,GRD_ZDIR)

          grhogvx(g,k,l) = dvvx(g,k,l) - prd * GRD_e(g,l,GRD_XDIR)
          grhogvy(g,k,l) = dvvy(g,k,l) - prd * GRD_e(g,l,GRD_YDIR)
          grhogvz(g,k,l) = dvvz(g,k,l) - prd * GRD_e(g,l,GRD_ZDIR)

          grhogwc(g,k,l) = prd * real(NON_HYDRO_ALPHA,kind=8)
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          grhogw(g,k,l) = 0.5D0 * ( GRD_afac(k) * grhogwc(g,k  ,l) &
                                  + GRD_bfac(k) * grhogwc(g,k-1,l) )
       enddo
       enddo

       grhogvx(:,ADM_kmin-1,l) = 0.D0
       grhogvx(:,ADM_kmax+1,l) = 0.D0
       grhogvy(:,ADM_kmin-1,l) = 0.D0
       grhogvy(:,ADM_kmax+1,l) = 0.D0
       grhogvz(:,ADM_kmin-1,l) = 0.D0
       grhogvz(:,ADM_kmax+1,l) = 0.D0
       grhogw (:,ADM_kmin-1,l) = 0.D0
       grhogw (:,ADM_kmin  ,l) = 0.D0
       grhogw (:,ADM_kmax+1,l) = 0.D0
    enddo
 
    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             dvvx_pl(g,k,l) = dvvx_pl(g,k,l) - 2.D0 * rhog_pl(g,k,l) * ( -CNST_EOHM * vvy_pl(g,k,l) )
             dvvy_pl(g,k,l) = dvvy_pl(g,k,l) - 2.D0 * rhog_pl(g,k,l) * (  CNST_EOHM * vvx_pl(g,k,l) )
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             prd = dvvx_pl(g,k,l) * GRD_e_pl(g,l,GRD_XDIR)  &
                 + dvvy_pl(g,k,l) * GRD_e_pl(g,l,GRD_YDIR)  &
                 + dvvz_pl(g,k,l) * GRD_e_pl(g,l,GRD_ZDIR)

             grhogvx_pl(g,k,l) = dvvx_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_XDIR)
             grhogvy_pl(g,k,l) = dvvy_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_YDIR)
             grhogvz_pl(g,k,l) = dvvz_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_ZDIR)

             grhogwc_pl(g,k,l) = prd * real(NON_HYDRO_ALPHA,kind=8)
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             grhogw_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * grhogwc_pl(g,k  ,l) &
                                        + GRD_bfac(k) * grhogwc_pl(g,k-1,l) )
          enddo
          enddo

          grhogvx_pl(:,ADM_kmin-1,l) = 0.D0
          grhogvy_pl(:,ADM_kmin-1,l) = 0.D0
          grhogvz_pl(:,ADM_kmin-1,l) = 0.D0
          grhogvx_pl(:,ADM_kmax+1,l) = 0.D0
          grhogvy_pl(:,ADM_kmax+1,l) = 0.D0
          grhogvz_pl(:,ADM_kmax+1,l) = 0.D0
          grhogw_pl (:,ADM_kmin-1,l) = 0.D0
          grhogw_pl (:,ADM_kmin  ,l) = 0.D0
          grhogw_pl (:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    return
  end subroutine src_advection_convergence_momentum

  !-----------------------------------------------------------------------------
  !> Advection convergence
  subroutine src_advection_convergence( &
       rhogvx,   rhogvx_pl,   &
       rhogvy,   rhogvy_pl,   &
       rhogvz,   rhogvz_pl,   &
       rhogw,    rhogw_pl,    &
       scl,      scl_pl,      &
       grhogscl, grhogscl_pl, &
       fluxtype               )
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

    real(8), intent(in)  :: rhogvx     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw      (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: scl        (ADM_gall,   ADM_kall,ADM_lall   ) ! scalar
    real(8), intent(in)  :: scl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogscl   (ADM_gall,   ADM_kall,ADM_lall   ) ! scalar tendency
    real(8), intent(out) :: grhogscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in)  :: fluxtype                                    ! scheme type
                                                                        ! I_SRC_horizontal : horizontal convergence
                                                                        ! I_SRC_vertical   : vertical convergence
                                                                        ! I_SRC_default    : both of them

    real(8) :: rhogvxscl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvxscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvyscl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvyscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvzscl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvzscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogwscl    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogwscl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: k, l
    !---------------------------------------------------------------------------

    rhogvxscl(:,:,:) = rhogvx(:,:,:) * scl(:,:,:)
    rhogvyscl(:,:,:) = rhogvy(:,:,:) * scl(:,:,:)
    rhogvzscl(:,:,:) = rhogvz(:,:,:) * scl(:,:,:)

    if ( ADM_prc_me == ADM_prc_pl ) then
       rhogvxscl_pl(:,:,:) = rhogvx_pl(:,:,:) * scl_pl(:,:,:)
       rhogvyscl_pl(:,:,:) = rhogvy_pl(:,:,:) * scl_pl(:,:,:)
       rhogvzscl_pl(:,:,:) = rhogvz_pl(:,:,:) * scl_pl(:,:,:)
    endif

    ! rhogwscl = rhow * e_w ( at half level ).
    if ( fluxtype == I_SRC_horizontal ) then

       rhogwscl(:,:,:) = 0.D0
       if ( ADM_prc_me == ADM_prc_pl ) then
          rhogwscl_pl(:,:,:) = 0.D0
       endif

    else

       do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
          rhogwscl(:,k,l) = rhogw(:,k,l) * 0.5D0 * ( GRD_afac(k) * scl(:,k,  l) &
                                                   + GRD_bfac(k) * scl(:,k-1,l) )
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             rhogwscl_pl(:,k,l) = rhogw_pl(:,k,l) * 0.5D0 * ( GRD_afac(k) * scl_pl(:,k  ,l) &
                                                            + GRD_bfac(k) * scl_pl(:,k-1,l) )
          enddo
          enddo
       endif

    endif

    !--- flux convergence step
    call src_flux_convergence( rhogvxscl, rhogvxscl_pl, & !--- [IN]
                               rhogvyscl, rhogvyscl_pl, & !--- [IN]
                               rhogvzscl, rhogvzscl_pl, & !--- [IN]
                               rhogwscl,  rhogwscl_pl,  & !--- [IN]
                               grhogscl,  grhogscl_pl,  & !--- [OUT]
                               fluxtype                 ) !--- [IN]

    return
  end subroutine src_advection_convergence

  !-----------------------------------------------------------------------------
  !> Flux convergence calculation
  !! 1. Horizontal flux convergence is calculated by using rhovx, rhovy, and 
  !!    rhovz which are defined at cell center (vertical) and A-grid (horizontal).
  !! 2. Vertical flux convergence is calculated by using rhovx, rhovy, rhovz, and rhow.
  !! 3. rhovx, rhovy, and rhovz can be replaced by rhovx*h, rhovy*h, and rhovz*h, respectively.
  !! 4. fluxtype can be set as below.
  subroutine src_flux_convergence( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       grhog,  grhog_pl,  &
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
       GRD_rdgz, &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_RGAM,       &
       VMTR_RGAM_pl,    &
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSH,       &
       VMTR_RGSH_pl,    &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_oprt, only: &
       OPRT_divergence
    implicit none

    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhog    (ADM_gall,   ADM_kall,ADM_lall   ) ! source
    real(8), intent(out) :: grhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in)  :: fluxtype                                    ! scheme type
                                                                        ! I_SRC_horizontal : horizontal convergence
                                                                        ! I_SRC_vertical   : vertical convergence
                                                                        ! I_SRC_default    : both of them

    real(8) :: div_flx_rhogvh   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: div_flx_rhogvh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: flx_rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: flx_rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: flx_rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: flx_rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: flx_rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    ! boundary condition
    do l = 1, ADM_lall
       flx_rhogw(:,ADM_kmin  ,l) = 0.D0
       flx_rhogw(:,ADM_kmax+1,l) = 0.D0
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          flx_rhogw_pl(:,ADM_kmin  ,l) = 0.D0
          flx_rhogw_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    if ( fluxtype == I_SRC_horizontal ) then ! Horizontal

       do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          flx_rhogw(g,k,l) = ( VMTR_C2Wfact(1,g,k,l) * rhogvx(g,k  ,l) &
                             + VMTR_C2Wfact(2,g,k,l) * rhogvx(g,k-1,l) &
                             + VMTR_C2Wfact(3,g,k,l) * rhogvy(g,k  ,l) &
                             + VMTR_C2Wfact(4,g,k,l) * rhogvy(g,k-1,l) &
                             + VMTR_C2Wfact(5,g,k,l) * rhogvz(g,k  ,l) &
                             + VMTR_C2Wfact(6,g,k,l) * rhogvz(g,k-1,l) &
                             ) * VMTR_RGAMH(g,k,l) ! horizontal contribution
       enddo
       enddo
       enddo

       flx_rhogvx(:,:,:) = rhogvx(:,:,:) * VMTR_RGAM(:,:,:)
       flx_rhogvy(:,:,:) = rhogvy(:,:,:) * VMTR_RGAM(:,:,:)
       flx_rhogvz(:,:,:) = rhogvz(:,:,:) * VMTR_RGAM(:,:,:)

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             flx_rhogw_pl(g,k,l) = ( VMTR_C2Wfact_pl(1,g,k,l) * rhogvx_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(2,g,k,l) * rhogvx_pl(g,k-1,l) &
                                   + VMTR_C2Wfact_pl(3,g,k,l) * rhogvy_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(4,g,k,l) * rhogvy_pl(g,k-1,l) &
                                   + VMTR_C2Wfact_pl(5,g,k,l) * rhogvz_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(6,g,k,l) * rhogvz_pl(g,k-1,l) &
                                   ) / VMTR_RGAMH_pl(g,k,l) ! horizontal contribution
          enddo
          enddo
          enddo

          flx_rhogvx_pl(:,:,:) = rhogvx_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          flx_rhogvy_pl(:,:,:) = rhogvy_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          flx_rhogvz_pl(:,:,:) = rhogvz_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
       endif

       !--- Horizontal flux convergence
       call OPRT_divergence( div_flx_rhogvh, div_flx_rhogvh_pl, & !--- [OUT]
                             flx_rhogvx,     flx_rhogvx_pl,     & !--- [IN]
                             flx_rhogvy,     flx_rhogvy_pl,     & !--- [IN]
                             flx_rhogvz,     flx_rhogvz_pl      ) !--- [IN]

    elseif( fluxtype == I_SRC_vertical ) then ! Vertical

       do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          flx_rhogw(g,k,l) = rhogw(g,k,l) * VMTR_RGSH(g,k,l) ! vertical contribution
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             flx_rhogw_pl(g,k,l) = rhogw_pl(g,k,l) * VMTR_RGSH_pl(g,k,l) ! vertical contribution
          enddo
          enddo
          enddo
       endif

       !--- Horizontal flux convergence
       div_flx_rhogvh   (:,:,:) = 0.D0
       div_flx_rhogvh_pl(:,:,:) = 0.D0

    else ! Default

       !--- horizontal + vertical contribution to flx_rhogw
       do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          flx_rhogw(g,k,l) = ( VMTR_C2Wfact(1,g,k,l) * rhogvx(g,k  ,l) &
                             + VMTR_C2Wfact(2,g,k,l) * rhogvx(g,k-1,l) &
                             + VMTR_C2Wfact(3,g,k,l) * rhogvy(g,k  ,l) &
                             + VMTR_C2Wfact(4,g,k,l) * rhogvy(g,k-1,l) &
                             + VMTR_C2Wfact(5,g,k,l) * rhogvz(g,k  ,l) &
                             + VMTR_C2Wfact(6,g,k,l) * rhogvz(g,k-1,l) &
                             ) * VMTR_RGAMH(g,k,l)                     & ! horizontal contribution
                           + rhogw(g,k,l) * VMTR_RGSH(g,k,l)             ! vertical   contribution
       enddo
       enddo
       enddo

       flx_rhogvx(:,:,:) = rhogvx(:,:,:) * VMTR_RGAM(:,:,:)
       flx_rhogvy(:,:,:) = rhogvy(:,:,:) * VMTR_RGAM(:,:,:)
       flx_rhogvz(:,:,:) = rhogvz(:,:,:) * VMTR_RGAM(:,:,:)

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             flx_rhogw_pl(g,k,l) = ( VMTR_C2Wfact_pl(1,g,k,l) * rhogvx_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(2,g,k,l) * rhogvx_pl(g,k-1,l) &
                                   + VMTR_C2Wfact_pl(3,g,k,l) * rhogvy_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(4,g,k,l) * rhogvy_pl(g,k-1,l) &
                                   + VMTR_C2Wfact_pl(5,g,k,l) * rhogvz_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(6,g,k,l) * rhogvz_pl(g,k-1,l) &
                                   ) / VMTR_RGAMH_pl(g,k,l)                        & ! horizontal contribution
                                 + rhogw_pl(g,k,l) * VMTR_RGSH_pl(g,k,l)             ! vertical   contribution
          enddo
          enddo
          enddo

          flx_rhogvx_pl(:,:,:) = rhogvx_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          flx_rhogvy_pl(:,:,:) = rhogvy_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          flx_rhogvz_pl(:,:,:) = rhogvz_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
       endif

       !--- Horizontal flux convergence
       call OPRT_divergence( div_flx_rhogvh, div_flx_rhogvh_pl, & !--- [OUT]
                             flx_rhogvx,     flx_rhogvx_pl,     & !--- [IN]
                             flx_rhogvy,     flx_rhogvy_pl,     & !--- [IN]
                             flx_rhogvz,     flx_rhogvz_pl      ) !--- [IN]

    endif

    !--- Total flux convergence
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          grhog(:,k,l) = - ( div_flx_rhogvh(:,k,l)                                 &
                           + ( flx_rhogw(:,k+1,l)-flx_rhogw(:,k,l) ) * GRD_rdgz(k) )
       enddo

       grhog(:,ADM_kmin-1,l) = 0.D0
       grhog(:,ADM_kmax+1,l) = 0.D0
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             grhog_pl(:,k,l) = - ( div_flx_rhogvh_pl(:,k,l)                                    &
                                 + ( flx_rhogw_pl(:,k+1,l)-flx_rhogw_pl(:,k,l) ) * GRD_rdgz(k) )
          enddo

          grhog_pl(:,ADM_kmin-1,l) = 0.D0
          grhog_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    return
  end subroutine src_flux_convergence

  !-----------------------------------------------------------------------------
  subroutine src_buoyancy(&
       rhog, rhog_pl,     & ! [IN] : perturb density ( gam2 X G^{1/2} )
       gbz, gbz_pl )        ! [OUT] : buoyancy force  at half level
    !------ 
    !------ Calculation of buoyacy force
    !------     1. Calculation region of gbz, gbz_pl
    !------                    : (:, ADM_kmin:ADM_kmax+1,:)
    !------ 
    !
    use mod_adm, only:  &
         ADM_VMISS,      &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_vmtr, only: &
         VMTR_RGSGAM2,   &
         VMTR_RGSGAM2_pl,&
         VMTR_GSGAM2H,   &
         VMTR_GSGAM2H_pl
    use mod_grd, only:  &
         GRD_afac,       &
         GRD_bfac
    use mod_cnst, only: &
         CNST_EGRAV
    implicit none

    real(8), intent(in)  :: rhog   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gbz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: gbz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    integer :: g,k,l

    !---- NOTICE! : Upward direction is positive for gbz.
    do l = 1, ADM_lall
       do g = 1, ADM_gall
          gbz(g,ADM_kmin-1,l) = ADM_VMISS
       enddo
    enddo
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall
             gbz(g,k,l) =                               &
                  - ( GRD_afac(k)*rhog(g,k,l)*VMTR_RGSGAM2(g,k,l)     &
                  + GRD_bfac(k)*rhog(g,k-1,l)*VMTR_RGSGAM2(g,k-1,l)  ) &
                  * VMTR_GSGAM2H(g,k,l) &
                  * CNST_EGRAV * 0.5D0
          enddo
       enddo
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do g = 1, ADM_gall_pl
             gbz_pl(g,ADM_kmin-1,l) = ADM_VMISS
          enddo
       enddo
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall_pl
                gbz_pl(g,k,l) =                               &
                     - ( GRD_afac(k)*rhog_pl(g,k,l)*VMTR_RGSGAM2_pl(g,k,l)     &
                     + GRD_bfac(k)*rhog_pl(g,k-1,l)*VMTR_RGSGAM2_pl(g,k-1,l)  ) &
                     * VMTR_GSGAM2H_pl(g,k,l)&
                     * CNST_EGRAV * 0.5D0
             enddo
          enddo
       enddo
    end if

  end subroutine src_buoyancy

  !-----------------------------------------------------------------------------
  subroutine src_pres_work_uv(&
       vx, vx_pl,             &  ! [IN] : Vx
       vy, vy_pl,             &  ! [IN] : Vy
       vz, vz_pl,             &  ! [IN] : Vz
       pgx, pgx_pl,           &  ! [IN] : pres. grad. force ( x comp. )
       pgy, pgy_pl,           &  ! [IN] : pres. grad. force ( y comp. )
       pgz, pgz_pl,           &  ! [IN] : pres. grad. force ( z comp. )
       grhoge, grhoge_pl)        ! [OUT] : pressure work
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
    use mod_adm, only:  &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_prc_me,     &
         ADM_prc_pl
    !
    implicit none

    real(8), intent(in) :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: pgx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: pgx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: pgy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: pgy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: pgz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: pgz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(out) :: grhoge   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhoge_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g,k,l





    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do g = 1, ADM_gall
             grhoge(g,k,l)                &
                  = vx(g,k,l) * pgx(g,k,l)&
                  + vy(g,k,l) * pgy(g,k,l)&
                  + vz(g,k,l) * pgz(g,k,l)
          enddo
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                grhoge_pl(g,k,l)                   &
                     = vx_pl(g,k,l) * pgx_pl(g,k,l)&
                     + vy_pl(g,k,l) * pgy_pl(g,k,l)&
                     + vz_pl(g,k,l) * pgz_pl(g,k,l)
             enddo
          enddo
       enddo
    endif

  end subroutine src_pres_work_uv

  !-----------------------------------------------------------------------------  
  subroutine src_gradient( &
       e, e_pl,            &  ! [IN]  : phi * G^{1/2} * gamma^2
       gex, gex_pl,        &  ! [OUT] : horizontal grad. ( x-comp. )
       gey, gey_pl,        &  ! [OUT] : horizontal grad. ( y-comp. )
       gez, gez_pl,        &  ! [OUT] : horizontal grad. ( z-comp. )
       gevz, gevz_pl,      &  ! [OUT] : vertical grad.   ( half point )
!       f_type )               ! [IN]  : calculation type
       grad_type            ) ! [IN] : calculation type [mod] H.Yashiro 20120530
    !------ 
    !------ Gradient operator for 3D
    !------     1. Calculation region of gex, gey, gez
    !------                    : (:, ADM_kmin:ADM_kmax,:)
    !------     2. Calculation region of gevz
    !------                    : (:, ADM_kmin:ADM_kmax+1,:)
    !------     3. Return values are
    !------ 
    !
    use mod_adm, only:  &
         ADM_VMISS,      &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_vmtr, only: &
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
    use mod_grd, only:  &
         GRD_afac,       &
         GRD_bfac,       &
         GRD_rdgz,       &
         GRD_rdgzh
    use mod_oprt, only:        &
         OPRT_horizontalize_vec,&
         OPRT_gradient
    implicit none

    real(8), intent(in) :: e   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: e_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out) :: gex   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: gex_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gey   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: gey_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gez   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: gez_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gevz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: gevz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: fex_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: fex_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: fey_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: fey_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: fez_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: fez_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
!    character(LEN=*),intent(in) :: f_type
    !<-----            HORIZONTAL          : horizontal gradient
    !<-----            VERTICAL            : vertical gradient ( do not used! )
    !<-----            DEFAULT             : both of them
    integer, intent(in) :: grad_type ! [mod] H.Yashiro 20120530
    !<-----                    = I_SRC_horizontal : horizontal gradient
    !<-----                    = I_SRC_vertical   : vertical gradient ( no use )
    !<-----                    = I_SRC_default    : both of them

    real(8) :: gee   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: gee_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: k,n,l






    !------ horizontal gradient without mountain
    do l=1,ADM_lall
       do k=1,ADM_kall
          do n=1,ADM_gall
             gex(n,k,l) = ADM_VMISS
          enddo
          do n=1,ADM_gall
             gey(n,k,l) = ADM_VMISS
          enddo
          do n=1,ADM_gall
             gez(n,k,l) = ADM_VMISS
          enddo
       enddo
    enddo

    do l=1,ADM_lall_pl
       do k=1,ADM_kall
          do n=1,ADM_gall_pl
             gex_pl(n,k,l) = ADM_VMISS
          enddo
          do n=1,ADM_gall_pl
             gey_pl(n,k,l) = ADM_VMISS
          enddo
          do n=1,ADM_gall_pl
             gez_pl(n,k,l) = ADM_VMISS
          enddo
       enddo
    enddo
    !
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=1,ADM_gall
          gee(n,k,l) = e(n,k,l)*VMTR_RGAM(n,k,l)
        enddo
      enddo
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
      do l=1,ADM_lall_pl
        do k=1,ADM_kall
          do n=1,ADM_gall_pl
            gee_pl(n,k,l) = e_pl(n,k,l) *VMTR_RGAM_pl(n,k,l) 
          enddo
        enddo
      enddo
    endif






    call OPRT_gradient(&
         gex, gex_pl,  &  ! [OUT]
         gey, gey_pl,  &  ! [OUT]
         gez, gez_pl,  &  ! [OUT]
         gee, gee_pl )    ! [IN]
    




    !--- horizontal gradient
    do l = 1, ADM_lall





       do k = ADM_kmin, ADM_kmax+1
          do n=1, ADM_gall
             fex_h(n,k,l)                                                  &
                  = VMTR_GZXH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          enddo
          do n=1, ADM_gall
             fey_h(n,k,l)                                                  &
                  = VMTR_GZYH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          enddo
          do n=1, ADM_gall
             fez_h(n,k,l)                                                  &
                  = VMTR_GZZH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          enddo
       enddo
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
          enddo
       enddo
       !




       !
       !--- At the lowest layer, do not use the extrapolation value!
       if(first_layer_remedy) then
         do n=1,ADM_gall
           gex(n,ADM_kmin,l) = gex(n,ADM_kmin+1,l)
           gey(n,ADM_kmin,l) = gey(n,ADM_kmin+1,l)
           gez(n,ADM_kmin,l) = gez(n,ADM_kmin+1,l)
         enddo
       endif





    enddo
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl





          do k = ADM_kmin, ADM_kmax+1
             do n=1, ADM_gall_pl
                fex_h_pl(n,k,l) &
                     = VMTR_GZXH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)
             enddo
             do n=1, ADM_gall_pl
                fey_h_pl(n,k,l)                                                     &
                     = VMTR_GZYH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)
             enddo
             do n=1, ADM_gall_pl
                fez_h_pl(n,k,l)                                                     &
                     = VMTR_GZZH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)

             enddo
          enddo






          do k = ADM_kmin, ADM_kmax
             do n=1, ADM_gall_pl
                gex_pl(n,k,l) = gex_pl(n,k,l) &
                     + (  fex_h_pl(n,k+1,l) - fex_h_pl(n,k,l) ) * GRD_rdgz(k)
                gey_pl(n,k,l) = gey_pl(n,k,l) &
                     + (  fey_h_pl(n,k+1,l) - fey_h_pl(n,k,l) ) * GRD_rdgz(k)
                gez_pl(n,k,l) = gez_pl(n,k,l) &
                     + (  fez_h_pl(n,k+1,l) - fez_h_pl(n,k,l) ) * GRD_rdgz(k)
             enddo
          enddo

          !--- At the lowest layer, do not use the extrapolation value!
          if(first_layer_remedy) then
            do n=1,ADM_gall_pl
              gex_pl(n,ADM_kmin,l) = gex_pl(n,ADM_kmin+1,l)
              gey_pl(n,ADM_kmin,l) = gey_pl(n,ADM_kmin+1,l)
              gez_pl(n,ADM_kmin,l) = gez_pl(n,ADM_kmin+1,l)
            enddo
          endif

       enddo
    endif

    !--- horizontalize
    call OPRT_horizontalize_vec( &
          gex, gex_pl,           & ! [IN]out
          gey, gey_pl,           & ! [IN]out
          gez, gez_pl )            ! [IN]out





    !--- vertical gradient ( note : half points )
!    if ( f_type == 'HORIZONTAL' ) then
    if ( grad_type == I_SRC_horizontal ) then ! [mod] H.Yashiro 20120530





       do l=1,ADM_lall
         do k=1,ADM_kall
           do n=1,ADM_gall
             gevz(n,k,l) = 0.D0
           enddo
         enddo
       enddo

       if(ADM_prc_me==ADM_prc_pl) then
         do l=1,ADM_lall_pl
           do k=1,ADM_kall
             do n=1,ADM_gall_pl 
               gevz_pl(n,k,l) = 0.D0
             enddo
           enddo
         enddo
       endif

    else

       do l=1,ADM_lall
         do n=1, ADM_gall
           gevz(n,ADM_kmin-1,l) = ADM_VMISS
         enddo
       enddo
       do l=1,ADM_lall
         do k = ADM_kmin, ADM_kmax+1
           do n=1,ADM_gall
             gevz(n,k,l)                                &
               = ( e(n,k,l) * VMTR_RGSGAM2(n,k,l)       &
                -  e(n,k-1,l) * VMTR_RGSGAM2(n,k-1,l) ) &
                 * GRD_rdgzh(k) * VMTR_GAM2H(n,k,l)
           enddo
         enddo
       enddo

       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
            do n=1,ADM_gall_pl
              gevz_pl(n,ADM_kmin-1,l) = ADM_VMISS
            enddo
            do k = ADM_kmin, ADM_kmax+1
              do n=1,ADM_gall_pl
                gevz_pl(n,k,l) &
                  = ( e_pl(n,k,l) * VMTR_RGSGAM2_pl(n,k,l)       &
                    - e_pl(n,k-1,l) * VMTR_RGSGAM2_pl(n,k-1,l) ) &
                    * GRD_rdgzh(k) * VMTR_GAM2H_pl(n,k,l)
              enddo
            enddo
          enddo
       endif
    endif

    return
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

    integer :: suf, i, j
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
