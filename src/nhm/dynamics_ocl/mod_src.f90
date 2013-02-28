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
  !                08-01-24   Y.Niwa: add src_update_tracer_rev
  !                                   save_memory => TRC_SAVE_MEMORY
  !                08-04-28   Y.Niwa: add initialization and avoid zero-dividing. 
  !                09-05-26   Y.Yamada: add directive, derected by T.Asano
  !                11-09-27   T.Seiki: merge optimized routines for K by RIST and M.Terai
  !                                    (see directive !OCL)
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
  public :: src_update_tracer_rev ! Y.Niwa add 080124

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

#ifdef _FJTIMER_
call timer_sta(2620)
#endif

    if ( iflag ) then
       iflag = .false.

       allocate( AFACovGSGAM2   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( AFACovGSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
       allocate( BFACovGSGAM2   (ADM_gall,   ADM_kall,ADM_lall   ) )
       allocate( BFACovGSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                   AFACovGSGAM2(ij,k,l) = GRD_afac(k) * VMTR_RGSGAM2(ij,k  ,l)
                   BFACovGSGAM2(ij,k,l) = GRD_bfac(k) * VMTR_RGSGAM2(ij,k-1,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
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
!OCL SERIAL
       do l = 1, ADM_lall
!OCL SERIAL
          do k = ADM_kmin, ADM_kmax+1
!OCL PARALLEL
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
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmin,l)   = 0.D0
          enddo
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL SERIAL
             do k = ADM_kmin, ADM_kmax+1
!OCL PARALLEL
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
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmin,l)   = 0.D0
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmax+1,l) = 0.D0
             enddo
          enddo
       endif

       !--- Horizontal flux convergence
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogvx_f(ij,k,l) = rhogvx(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogvy_f(ij,k,l) = rhogvy(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogvz_f(ij,k,l) = rhogvz(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rhogvx_f_pl(ij,k,l) = rhogvx_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rhogvy_f_pl(ij,k,l) = rhogvy_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
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
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                flx_vz(ij,k,l) = rhogw(ij,k,l) * VMTR_RGSH(ij,k,l)
             enddo
          enddo
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmin  ,l) = 0.D0
          enddo
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   flx_vz_pl(ij,k,l) = rhogw_pl(ij,k,l) * VMTR_RGSH_pl(ij,k,l)
                enddo
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmin  ,l) = 0.D0
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmax+1,l) = 0.D0
             enddo
          enddo
       endif

       !--- Horizontal flux convergence
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
             do ij = 1, ADM_gall
                hdiv(ij,k,l) = 0.D0
             enddo
          enddo
       enddo
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
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
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
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
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmin  ,l) = 0.D0
          enddo
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             flx_vz(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
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
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmin  ,l) = 0.D0
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                flx_vz_pl(ij,ADM_kmax+1,l) = 0.D0
             enddo
          enddo
       endif

       !--- Horizontal flux convergence
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogvx_f(ij,k,l) = rhogvx(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogvy_f(ij,k,l) = rhogvy(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do ij = 1, ADM_gall
                rhogvz_f(ij,k,l) = rhogvz(ij,k,l) * VMTR_RGAM(ij,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rhogvx_f_pl(ij,k,l) = rhogvx_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
                do ij = 1, ADM_gall_pl
                   rhogvy_f_pl(ij,k,l) = rhogvy_pl(ij,k,l) * VMTR_RGAM_pl(ij,k,l)
                enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
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
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax
!OCL XFILL, PREFETCH_STRONG
          do ij = 1, ADM_gall
             grhog(ij,k,l) = - ( flx_vz(ij,k+1,l) - flx_vz(ij,k,l) ) &
                             / ( GRD_gzh(k+1)     - GRD_gzh(k)     ) - hdiv(ij,k,l)
          enddo
       enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do ij = 1, ADM_gall
          grhog(ij,ADM_kmin-1,l) = 0.D0
       enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do ij = 1, ADM_gall
          grhog(ij,ADM_kmax+1,l) = 0.D0
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
             do ij = 1, ADM_gall_pl
                grhog_pl(ij,k,l) = - ( flx_vz_pl(ij,k+1,l) - flx_vz_pl(ij,k,l) ) &
                                   / ( GRD_gzh(k+1)        - GRD_gzh(k)        ) - hdiv_pl(ij,k,l)
             enddo
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             grhog_pl(ij,ADM_kmin-1,l) = 0.D0
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
          do ij = 1, ADM_gall_pl
             grhog_pl(ij,ADM_kmax+1,l) = 0.D0
          enddo
       enddo
    endif

#ifdef _FJTIMER_
call timer_end(2620)
#endif

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
         VMTR_RGSGAM2,   &
         VMTR_RGSGAM2_pl,&
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
         GRD_gzh,        &
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

    !--- horizontal contribution to flx_vz
    flx_vz = 0.0D0
    flx_vz_pl = 0.0D0
!    if(f_type=='VERTICAL') then
    if ( fluxtype == I_SRC_vertical ) then ! [mod] H.Yashiro 20120530
       flx_vz(:,:,:) = 0.0D0
       flx_vz_pl(:,:,:) = 0.0D0
    else !--- default or 'HORIZONTAL'
!OCL SERIAL
       do l=1,ADM_lall
!OCL PARALLEL
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
!OCL SERIAL
          do l=1,ADM_LALL_PL
!OCL PARALLEL
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
!OCL SERIAL
       do l=1,ADM_lall
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
             flx_vz(:,k,l)                                    &
                  = flx_vz(:,k,l)                             &
                  + rhogw(:,k,l) * VMTR_RGSH(:,k,l)
          end do
       end do
       if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
          do l=1,ADM_LALL_PL
!OCL PARALLEL
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
!OCL SERIAL
    do l=1,ADM_lall
!OCL PARALLEL
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
!OCL SERIAL
       do l=1,ADM_lall_pl
!OCL PARALLEL
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

#ifdef _FJTIMER_
call timer_sta(2640)
#endif

    !---- NOTICE! : Upward direction is positive for gbz.
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG
       do g = 1, ADM_gall
          gbz(g,ADM_kmin-1,l) = ADM_VMISS
       end do
    end do
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG
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
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
          do g = 1, ADM_gall_pl
             gbz_pl(g,ADM_kmin-1,l) = ADM_VMISS
          end do
       end do
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
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

#ifdef _FJTIMER_
call timer_end(2640)
#endif

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

#ifdef _FJTIMER_
call timer_sta(2650)
#endif

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogvxe(ij,k,l) = rhogvx(ij,k,l) * e(ij,k,l)
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogvye(ij,k,l) = rhogvy(ij,k,l) * e(ij,k,l)
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do ij = 1, ADM_gall
             rhogvze(ij,k,l) = rhogvz(ij,k,l) * e(ij,k,l)
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                rhogvxe_pl(ij,k,l) = rhogvx_pl(ij,k,l) * e_pl(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do ij = 1, ADM_gall_pl
                rhogvye_pl(ij,k,l) = rhogvy_pl(ij,k,l) * e_pl(ij,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
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
!OCL SERIAL
       do l = 1, ADM_lall
!OCL SERIAL
          do k =ADM_kmin, ADM_kmax+1
!OCL PARALLEL
             do g = 1, ADM_gall
                fez(g,k,l) = 0.5D0 * ( GRD_afac(k) * e(g,k,  l) &
                                     + GRD_bfac(k) * e(g,k-1,l) ) * rhogw(g,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL SERIAL
             do k = ADM_kmin, ADM_kmax+1
!OCL PARALLEL
                do g = 1, ADM_gall_pl
                   fez_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * e_pl(g,k,  l) &
                                           + GRD_bfac(k) * e_pl(g,k-1,l) ) * rhogw_pl(g,k,l)
                enddo
             enddo
          enddo
       endif

    else !--- = 'HORIZONTAL'

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
             do g = 1, ADM_gall
                fez(g,k,l) = 0.D0
             enddo
          enddo
       enddo
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
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


#ifdef _FJTIMER_
call timer_end(2650)
#endif

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

    return
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

#ifdef _FJTIMER_
call timer_sta(2670)
#endif

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG
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
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, UNROLL('full')
             do g = 1, ADM_gall_pl
                grhoge_pl(g,k,l)                   &
                     = vx_pl(g,k,l) * pgx_pl(g,k,l)&
                     + vy_pl(g,k,l) * pgy_pl(g,k,l)&
                     + vz_pl(g,k,l) * pgz_pl(g,k,l)
             end do
          end do
       end do
    end if

#ifdef _FJTIMER_
call timer_end(2670)
#endif

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
    real(8) :: Omega(GRD_XDIR:GRD_ZDIR)
    !
    real(8) :: alpha = 0.0D0

    integer :: g,k,l

#ifdef _FJTIMER_
call timer_sta(2680)
#endif

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k=ADM_kmin,ADM_kmax
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             vvx(g,k,l) = vx(g,k,l) &
                  + ( GRD_cfac(k) * w(g,k+1,l) + GRD_dfac(k) * w(g,k,l) )  * 0.5D0&
                  * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
          end do
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             vvy(g,k,l) = vy(g,k,l) &
                  + ( GRD_cfac(k) * w(g,k+1,l) + GRD_dfac(k) * w(g,k,l) )  * 0.5D0&
                  * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
          end do
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             vvz(g,k,l) = vz(g,k,l) &
                  + ( GRD_cfac(k) * w(g,k+1,l) + GRD_dfac(k) * w(g,k,l) )  * 0.5D0&
                  * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
          end do
       end do
    end do
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          vvx(g,ADM_kmin-1,l) = 0.0D0
       end do
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          vvx(g,ADM_kmax+1,l) = 0.0D0
       end do
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          vvy(g,ADM_kmin-1,l) = 0.0D0
       end do
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          vvy(g,ADM_kmax+1,l) = 0.0D0
       end do
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          vvz(g,ADM_kmin-1,l) = 0.0D0
       end do
!OCL PARALLEL, XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          vvz(g,ADM_kmax+1,l) = 0.0D0
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k=ADM_kmin,ADM_kmax
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do g = 1, ADM_gall_pl
                vvx_pl(g,k,l) = vx_pl(g,k,l) &
                     + ( GRD_cfac(k) * w_pl(g,k+1,l) + GRD_dfac(k) * w_pl(g,k,l) )  * 0.5D0&
                     * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
             end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
             do g = 1, ADM_gall_pl
                vvy_pl(g,k,l) = vy_pl(g,k,l) &
                     + ( GRD_cfac(k) * w_pl(g,k+1,l) + GRD_dfac(k) * w_pl(g,k,l) )  * 0.5D0&
                     * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
             end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
             do g = 1, ADM_gall_pl
                vvz_pl(g,k,l) = vz_pl(g,k,l) &
                     + ( GRD_cfac(k) * w_pl(g,k+1,l) + GRD_dfac(k) * w_pl(g,k,l) )  * 0.5D0&
                     * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
             end do
          end do
       end do
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             vvx_pl(g,ADM_kmin-1,l) = 0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             vvx_pl(g,ADM_kmax+1,l) = 0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             vvy_pl(g,ADM_kmin-1,l) = 0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             vvy_pl(g,ADM_kmax+1,l) = 0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             vvz_pl(g,ADM_kmin-1,l) = 0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
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

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
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
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogvx(g,ADM_kmin-1,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogvy(g,ADM_kmin-1,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogvz(g,ADM_kmin-1,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogvx(g,ADM_kmax+1,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogvy(g,ADM_kmax+1,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogvz(g,ADM_kmax+1,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogw(g,ADM_kmin-1,l) = 0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogw(g,ADM_kmin  ,l)=0.0D0
       end do
!OCL PARALLEL,XFILL, PREFETCH_STRONG,LOOP_NOFUSION
       do g = 1, ADM_gall
          grhogw(g,ADM_kmax+1,l)=0.0D0
       end do
    end do
    !
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k=ADM_kmin+1,ADM_kmax
!OCL XFILL, PREFETCH_STRONG
          do g = 1, ADM_gall
             grhogw(g,k,l) &
                  = ( GRD_afac(k) * grhogwc(g,k,l)  &
                    + GRD_bfac(k) * grhogwc(g,k-1,l) ) * 0.5D0
          end do
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
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
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogvx_pl(g,ADM_kmin-1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogvy_pl(g,ADM_kmin-1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogvz_pl(g,ADM_kmin-1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogvx_pl(g,ADM_kmax+1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogvy_pl(g,ADM_kmax+1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogvz_pl(g,ADM_kmax+1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogw_pl(g,ADM_kmin-1,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogw_pl(g,ADM_kmin  ,l)=0.0D0
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do g = 1, ADM_gall_pl
             grhogw_pl(g,ADM_kmax+1,l)=0.0D0
          end do
       end do
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k=ADM_kmin+1,ADM_kmax
!OCL XFILL, PREFETCH_STRONG,UNROLL('full')
             do g = 1, ADM_gall_pl
                grhogw_pl(g,k,l) &
                      = ( GRD_afac(k) * grhogwc_pl(g,k,l)  &
                        + GRD_bfac(k) * grhogwc_pl(g,k-1,l) ) * 0.5D0
             end do
          end do
       end do
    end if

#ifdef _FJTIMER_
call timer_end(2680)
#endif

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
         VMTR_GZZ,       &
         VMTR_GZZ_pl,    &
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
    !
    !
    implicit none
    !
    real(8), intent(in) :: e(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: e_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8), intent(out) :: gex(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gex_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gey(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gey_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gez(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gez_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: gevz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: gevz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
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
!    real(8) :: ge(ADM_gall,ADM_kall,ADM_lall)
!    real(8) :: ge_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: gee(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: gee_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    integer :: k,n,l

#ifdef _FJTIMER_
call timer_sta(2690)
#endif

!    ge(:,:,:)    = e(:,:,:)
!    if(ADM_prc_me==ADM_prc_pl) then
!       ge_pl(:,:,:) = e_pl(:,:,:)
!    end if
    !
    !------ horizontal gradient without mountain
!OCL SERIAL
    do l=1,ADM_lall
!OCL PARALLEL
       do k=1,ADM_kall
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION
          do n=1,ADM_gall
             gex(n,k,l) = ADM_VMISS
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION
          do n=1,ADM_gall
             gey(n,k,l) = ADM_VMISS
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION
          do n=1,ADM_gall
             gez(n,k,l) = ADM_VMISS
          end do
       end do
    end do

!OCL SERIAL
    do l=1,ADM_lall_pl
!OCL PARALLEL
       do k=1,ADM_kall
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do n=1,ADM_gall_pl
             gex_pl(n,k,l) = ADM_VMISS
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do n=1,ADM_gall_pl
             gey_pl(n,k,l) = ADM_VMISS
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
          do n=1,ADM_gall_pl
             gez_pl(n,k,l) = ADM_VMISS
          end do
       end do
    end do
    !
!OCL SERIAL
    do l=1,ADM_lall
!OCL PARALLEL
      do k=1,ADM_kall
!OCL XFILL, PREFETCH_STRONG
        do n=1,ADM_gall
          gee(n,k,l) = e(n,k,l)*VMTR_RGAM(n,k,l)
        end do
      end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
      do l=1,ADM_lall_pl
!OCL PARALLEL
        do k=1,ADM_kall
!OCL XFILL, PREFETCH_STRONG,UNROLL('full')
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
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION
          do n=1, ADM_gall
             fex_h(n,k,l)                                                  &
                  = VMTR_GZXH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION
          do n=1, ADM_gall
             fey_h(n,k,l)                                                  &
                  = VMTR_GZYH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION
          do n=1, ADM_gall
             fez_h(n,k,l)                                                  &
                  = VMTR_GZZH(n,k,l)                                       &
                  * ( GRD_afac(k)*e(n,k,l)*VMTR_RGSGAM2(n,k,l)             &
                     +GRD_bfac(k)*e(n,k-1,l)*VMTR_RGSGAM2(n,k-1,l) )*0.5D0 &
                  * VMTR_GSGAMH(n,k,l)
          end do
       end do

!OCL PARALLEL
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

       !--- At the lowest layer, do not use the extrapolation value!
       if(first_layer_remedy) then
!OCL PARALLEL
         do n=1,ADM_gall
           gex(n,ADM_kmin,l) = gex(n,ADM_kmin+1,l)
           gey(n,ADM_kmin,l) = gey(n,ADM_kmin+1,l)
           gez(n,ADM_kmin,l) = gez(n,ADM_kmin+1,l)
         end do
       end if
!       else
!          k = ADM_kmin
!          do n=1, ADM_gall
!             fex_h(n,k,l)                                     &
!                  = VMTR_GZX(n,k,l) * gee(n,k,l)
!             fey_h(n,k,l)                                     &
!                  = VMTR_GZY(n,k,l) * gee(n,k,l)
!             fez_h(n,k,l)                                     &
!                  = VMTR_GZZ(n,k,l) * gee(n,k,l)
!          end do

!          do n=1, ADM_gall
!             gex(n,k,l) &
!                  = ( 0.5D0*(GRD_cfac(k)*gex(n,k+1,l)+GRD_dfac(k)*gex(n,k,l))&
!                  + gex(n,k,l) ) * 0.5D0
!             gey(n,k,l) &
!                  = ( 0.5D0*(GRD_cfac(k)*gey(n,k+1,l)+GRD_dfac(k)*gey(n,k,l))&
!                  + gey(n,k,l) ) * 0.5D0
!             gez(n,k,l) &
!                  = ( 0.5D0*(GRD_cfac(k)*gez(n,k+1,l)+GRD_dfac(k)*gez(n,k,l))&
!                  + gez(n,k,l) ) * 0.5D0
!          end do
!          !
!          do n=1, ADM_gall
!             gex(n,k,l)                                             &
!                  = gex(n,k,l)                                      &
!                  + (  fex_h(n,k+1,l) - fex_h(n,k,l) ) * GRD_rdgz(k) * 2.0D0
!             gey(n,k,l)                                             &
!                  = gey(n,k,l)                                      &
!                  + (  fey_h(n,k+1,l) - fey_h(n,k,l) ) * GRD_rdgz(k) * 2.0D0
!             gez(n,k,l)                                             &
!                  = gez(n,k,l)                                      &
!                  + (  fez_h(n,k+1,l) - fez_h(n,k,l) ) * GRD_rdgz(k) * 2.0D0
!          end do
!       end if
    end do

    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
             do n=1, ADM_gall_pl
                fex_h_pl(n,k,l) &
                     = VMTR_GZXH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)
             end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
             do n=1, ADM_gall_pl
                fey_h_pl(n,k,l)                                                     &
                     = VMTR_GZYH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)
             end do
!OCL XFILL, PREFETCH_STRONG,LOOP_NOFUSION,UNROLL('full')
             do n=1, ADM_gall_pl
                fez_h_pl(n,k,l)                                                     &
                     = VMTR_GZZH_pl(n,k,l)                                          &
                     * ( GRD_afac(k)*e_pl(n,k,l)*VMTR_RGSGAM2_pl(n,k,l)             &
                        +GRD_bfac(k)*e_pl(n,k-1,l)*VMTR_RGSGAM2_pl(n,k-1,l) )*0.5D0 &
                     * VMTR_GSGAMH_pl(n,k,l)

             end do
          end do

!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
!OCL UNROLL('full')
             do n=1, ADM_gall_pl
                gex_pl(n,k,l) = gex_pl(n,k,l) &
                     + (  fex_h_pl(n,k+1,l) - fex_h_pl(n,k,l) ) * GRD_rdgz(k)
                gey_pl(n,k,l) = gey_pl(n,k,l) &
                     + (  fey_h_pl(n,k+1,l) - fey_h_pl(n,k,l) ) * GRD_rdgz(k)
                gez_pl(n,k,l) = gez_pl(n,k,l) &
                     + (  fez_h_pl(n,k+1,l) - fez_h_pl(n,k,l) ) * GRD_rdgz(k)
             end do
          end do

          !--- At the lowest layer, do not use the extrapolation value!
          if(first_layer_remedy) then
!OCL SERIAL,UNROLL('full')
            do n=1,ADM_gall_pl
              gex_pl(n,ADM_kmin,l) = gex_pl(n,ADM_kmin+1,l)
              gey_pl(n,ADM_kmin,l) = gey_pl(n,ADM_kmin+1,l)
              gez_pl(n,ADM_kmin,l) = gez_pl(n,ADM_kmin+1,l)
            end do
          end if
!          else
!             k = ADM_kmin
!             do n=1, ADM_gall_pl
!                fex_h_pl(n,k,l)                                  &
!                     = VMTR_GZX_pl(n,k,l) * gee_pl(n,k,l)
!                fey_h_pl(n,k,l)                                  &
!                     = VMTR_GZY_pl(n,k,l) * gee_pl(n,k,l)
!                fez_h_pl(n,k,l)                                  &
!                     = VMTR_GZZ_pl(n,k,l) * gee_pl(n,k,l)
!             end do
!             do n=1, ADM_gall_pl
!                gex_pl(n,k,l) &
!                     = ( 0.5D0*(GRD_cfac(k)*gex_pl(n,k+1,l)+GRD_dfac(k)*gex_pl(n,k,l))&
!                     + gex_pl(n,k,l) ) * 0.5D0
!                gey_pl(n,k,l) &
!                     = ( 0.5D0*(GRD_cfac(k)*gey_pl(n,k+1,l)+GRD_dfac(k)*gey_pl(n,k,l))&
!                     + gey_pl(n,k,l) ) * 0.5D0
!                gez_pl(n,k,l) &
!                     = ( 0.5D0*(GRD_cfac(k)*gez_pl(n,k+1,l)+GRD_dfac(k)*gez_pl(n,k,l))&
!                     + gez_pl(n,k,l) ) * 0.5D0
!             end do
             
!             do n=1, ADM_gall_pl
!                gex_pl(n,k,l)                                             &
!                     = gex_pl(n,k,l)                                      &
!                     + (  fex_h_pl(n,k+1,l) - fex_h_pl(n,k,l) ) * GRD_rdgz(k) * 2.0D0
!                gey_pl(n,k,l)                                             &
!                     = gey_pl(n,k,l)                                      &
!                     + (  fey_h_pl(n,k+1,l) - fey_h_pl(n,k,l) ) * GRD_rdgz(k) * 2.0D0
!                gez_pl(n,k,l)                                             &
!                     = gez_pl(n,k,l)                                      &
!                     + (  fez_h_pl(n,k+1,l) - fez_h_pl(n,k,l) ) * GRD_rdgz(k) * 2.0D0
!             end do
!          end if
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
!OCL SERIAL
       do l=1,ADM_lall
!OCL PARALLEL
         do k=1,ADM_kall
!OCL XFILL, PREFETCH_STRONG
           do n=1,ADM_gall
             gevz(n,k,l) = 0.0D0
           end do
         end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
         do l=1,ADM_lall_pl
!OCL PARALLEL
           do k=1,ADM_kall
!OCL XFILL, PREFETCH_STRONG,UNROLL('full')
             do n=1,ADM_gall_pl 
               gevz_pl(n,k,l) = 0.0D0
             end do
           end do
         end do
       end if
    else
!OCL SERIAL
       do l=1,ADM_lall
!OCL PARALLEL,XFILL, PREFETCH_STRONG
         do n=1, ADM_gall
           gevz(n,ADM_kmin-1,l) = ADM_VMISS
         end do
       end do
!OCL SERIAL
       do l=1,ADM_lall
!OCL PARALLEL
         do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG
           do n=1,ADM_gall
             gevz(n,k,l)                                &
               = ( e(n,k,l) * VMTR_RGSGAM2(n,k,l)       &
                -  e(n,k-1,l) * VMTR_RGSGAM2(n,k-1,l) ) &
                 * GRD_rdgzh(k) * VMTR_GAM2H(n,k,l)
           end do
         end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
          do l=1,ADM_lall_pl
!OCL XFILL, PREFETCH_STRONG,UNROLL('full')
            do n=1,ADM_gall_pl
              gevz_pl(n,ADM_kmin-1,l) = ADM_VMISS
            end do
!OCL PARALLEL
            do k = ADM_kmin, ADM_kmax+1
!OCL XFILL, PREFETCH_STRONG,UNROLL('full')
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

#ifdef _FJTIMER_
call timer_end(2690)
#endif

  end subroutine src_gradient
  !-----------------------------------------------------------------------------  
  subroutine src_update_tracer(        &
       rhogq,  rhogq_pl,               & !--- INOUT rhoq    ( gam2 X G^{1/2} )
       rhog,   rhog_pl,                & !--- IN  : rho     ( gam2 X G^{1/2} )
       rhogvx, rhogvx_pl,              & !--- IN  : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy, rhogvy_pl,              & !--- IN  : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz, rhogvz_pl,              & !--- IN  : rho*Vz  ( gam2 X G^{1/2} )
       rhogw,  rhogw_pl,               & !--- IN  : rho*w   ( gam2 X G^{1/2} )
       dt                              & !--- IN  : delta t
       )
    !------ 
    !------ Update of tracer by advection
    !------   07/11/13 H.Tomita: Change the vertical limiter.
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
         ADM_prc_pl,     &
         ADM_AI,ADM_AJ,ADM_AIJ,&
         ADM_prc_run_master, &
         ADM_LOG_FID,&
         ADM_KNONE
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
         VMTR_RGAM_pl,   &
         VMTR_GSGAM2H,   &
         VMTR_GSGAM2H_pl
         
    use mod_grd, only :  &
         GRD_gzh,        &
         GRD_gz,         &
         GRD_afac,       &
         GRD_bfac,       &
         GRD_XDIR,       &
         GRD_YDIR,       &
         GRD_ZDIR,       &
         GRD_dgz
    use mod_oprt, only :       &
         oprt_divergence2_prep,&
         oprt_divergence2_all, &
         oprt_divergence2      ! K.Suzuki [add] 07/04/25
    use mod_gtl, only : &
         GTL_max_k,GTL_min_k,&
         GTL_global_sum_srf
    use mod_cnst, only : &
         CNST_MAX_REAL,  &
         CNST_EPS_ZERO
    use mod_runconf, only : &
         nqmax =>TRC_VMAX,  &
         I_QV, I_QC,        &
         TRC_SAVE_MEMORY
    !
    implicit none
    !
    real(8), intent(inout) :: rhogq(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8), intent(inout) :: rhogq_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)
    real(8) :: q(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: q_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)
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
    real(8), intent(in) :: dt
    !
    real(8) :: flx_vz(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: flx_vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)
    !
    real(8) :: vx(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: vy(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: vz(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8) :: vx_r(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vx_r_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: vy_r(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vy_r_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: vz_r(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vz_r_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    !
    real(8) :: c(ADM_gall,ADM_kall,ADM_lall,6)
    real(8) :: c_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: cp(ADM_gall,ADM_kall,ADM_lall,ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8) :: cp_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,GRD_XDIR:GRD_ZDIR)

    real(8) :: hdiv(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: hdiv_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)
    real(8) :: wh(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: wh_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)

    real(8) :: rhogqv_max(ADM_kall)
    real(8) :: rhogqv_min(ADM_kall)
    real(8) :: rhogqc_max(ADM_kall)
    real(8) :: rhogqc_min(ADM_kall)

    real(8) :: ck(ADM_gall,ADM_kall,ADM_lall,2)
    real(8) :: ck_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,2)

    real(8) :: q_in_min(ADM_gall,ADM_kall,ADM_lall,2,nqmax)
    real(8) :: q_in_min_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,2,nqmax)
    real(8) :: q_in_max(ADM_gall,ADM_kall,ADM_lall,2,nqmax)
    real(8) :: q_in_max_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,2,nqmax)

    real(8) ::  q_m1_k_min(ADM_gall)
    real(8) ::  q_m1_k_max(ADM_gall)
    real(8) ::  q_m1_k_min_pl(ADM_gall_pl)
    real(8) ::  q_m1_k_max_pl(ADM_gall_pl)

    real(8) ::  c_in_sum(ADM_gall)
    real(8) ::  c_out_sum(ADM_gall)
    real(8) ::  c_in_sum_pl(ADM_gall_pl)
    real(8) ::  c_out_sum_pl(ADM_gall_pl)
    
    real(8) ::  c_qin_sum_max(ADM_gall)
    real(8) ::  c_qin_sum_min(ADM_gall)
    real(8) ::  c_qin_sum_max_pl(ADM_gall_pl)
    real(8) ::  c_qin_sum_min_pl(ADM_gall_pl)

    real(8) :: q_out_k_min(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: q_out_k_min_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)
    real(8) :: q_out_k_max(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: q_out_k_max_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)

    real(8) :: q_h(ADM_gall,ADM_kall,ADM_lall,nqmax)
    real(8) :: q_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nqmax)

    real(8) :: rho_h(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rho_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    
    real(8) :: SMALL_ZERO = 0.0D0

!del Y.Niwa 080124    logical :: save_memory = .true.  ! K.Suzuki [add] 07/04/25
    !
    integer :: k,l,n,nq
    !
    !
    !-->
!    do k=ADM_kmin,ADM_kmax
!       rhogqv_max(k) = GTL_max_k(rhogq(:,:,:,I_QV),rhogq_pl(:,:,:,I_QV),k)
!       rhogqv_min(k) = GTL_min_k(rhogq(:,:,:,I_QV),rhogq_pl(:,:,:,I_QV),k)
!       rhogqc_max(k) = GTL_max_k(rhogq(:,:,:,I_QC),rhogq_pl(:,:,:,I_QC),k)
!       rhogqc_min(k) = GTL_min_k(rhogq(:,:,:,I_QC),rhogq_pl(:,:,:,I_QC),k)
!    end do
!    if(ADM_prc_me==ADM_prc_run_master) then
!       write(ADM_LOG_FID,*) '------------- start of src_update_tracer'
!       do k = ADM_kmin, ADM_kmax
!          write(ADM_LOG_FID,'(I4,4E14.6)') k, rhogqv_max(k), rhogqv_min(k), rhogqc_max(k), rhogqc_min(k)
!       end do
!    end if
    !<---
    !
    !
    !--- horizontal velocity
    vx(:,:,:)=rhogvx(:,:,:)/rhog(:,:,:)
    vy(:,:,:)=rhogvy(:,:,:)/rhog(:,:,:)
    vz(:,:,:)=rhogvz(:,:,:)/rhog(:,:,:)
    !
    if(ADM_prc_me==ADM_prc_pl) then
       vx_pl(:,:,:)=rhogvx_pl(:,:,:)/rhog_pl(:,:,:)
       vy_pl(:,:,:)=rhogvy_pl(:,:,:)/rhog_pl(:,:,:)
       vz_pl(:,:,:)=rhogvz_pl(:,:,:)/rhog_pl(:,:,:)
    end if
    !
    !--- rho at half level
!OCL SERIAL
    do l=1, ADM_lall 
!OCL PARALLEL
       do k =ADM_kmin,ADM_kmax+1
          rho_h(:,k,l) = 0.5D0 *         &
               ( GRD_afac(k) * rhog(:,k  ,l)* VMTR_RGSGAM2(:,k  ,l)&
               + GRD_bfac(k) * rhog(:,k-1,l)* VMTR_RGSGAM2(:,k-1,l)&
               ) 
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1, ADM_lall_pl
!OCL PARALLEL
          do k =ADM_kmin,ADM_kmax+1
             rho_h_pl(:,k,l) = 0.5D0 *         &
                  ( GRD_afac(k) * rhog_pl(:,k  ,l)* VMTR_RGSGAM2_pl(:,k  ,l)&
                  + GRD_bfac(k) * rhog_pl(:,k-1,l)* VMTR_RGSGAM2_pl(:,k-1,l)&
                  )
          end do
       end do
    end if
    !
    !
!OCL SERIAL
    do l=1,ADM_lall
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax+1
          wh(:,k,l) =( (                                                 &
               ( GRD_afac(k)*VMTR_RGSGAM2(:,k  ,l)*rhogvx(:,k  ,  l)    &
                +GRD_bfac(k)*VMTR_RGSGAM2(:,k-1,l)*rhogvx(:,k-1,  l) )  &
                * 0.5D0*VMTR_GSGAMH(:,k,l)*VMTR_GZXH(:,k,l)             &
                +( GRD_afac(k)*VMTR_RGSGAM2(:,k  ,l)*rhogvy(:,k  ,  l)  &
                +GRD_bfac(k)*VMTR_RGSGAM2(:,k-1,l)*rhogvy(:,k-1,  l) )  &
                * 0.5D0*VMTR_GSGAMH(:,k,l)*VMTR_GZYH(:,k,l)             & 
                +( GRD_afac(k)*VMTR_RGSGAM2(:,k  ,l)*rhogvz(:,k  ,  l)  &
                +GRD_bfac(k)*VMTR_RGSGAM2(:,k-1,l)*rhogvz(:,k-1,  l) )  &
                * 0.5D0*VMTR_GSGAMH(:,k,l)*VMTR_GZZH(:,k,l)             &
                )&
                + rhogw(:,k,l) * VMTR_RGSH(:,k,l) )/rho_h(:,k,l)

       end do
       wh(:,ADM_kmin,l) = 0.0D0
       wh(:,ADM_kmax+1,l) = 0.0D0
    end do
    !
    !
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1,ADM_LALL_PL
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
             wh_pl(:,k,l) = ((                                               &
                  ( GRD_afac(k)*VMTR_RGSGAM2_pl(:,k  ,l)*rhogvx_pl(:,k  ,  l)   &
                  +GRD_bfac(k)*VMTR_RGSGAM2_pl(:,k-1,l)*rhogvx_pl(:,k-1,  l) )  &
                  * 0.5D0*VMTR_GSGAMH_pl(:,k,l)*VMTR_GZXH_pl(:,k,l)             &
                  + ( GRD_afac(k)*VMTR_RGSGAM2_pl(:,k  ,l)*rhogvy_pl(:,k  ,  l) &
                  +GRD_bfac(k)*VMTR_RGSGAM2_pl(:,k-1,l)*rhogvy_pl(:,k-1,  l) )  &
                  * 0.5D0*VMTR_GSGAMH_pl(:,k,l)*VMTR_GZYH_pl(:,k,l)             & 
                  +( GRD_afac(k)*VMTR_RGSGAM2_pl(:,k  ,l)*rhogvz_pl(:,k  ,  l)  &
                  +GRD_bfac(k)*VMTR_RGSGAM2_pl(:,k-1,l)*rhogvz_pl(:,k-1,  l) )  &
                  * 0.5D0*VMTR_GSGAMH_pl(:,k,l)*VMTR_GZZH_pl(:,k,l)             &
                  )&
                  + rhogw_pl(:,k,l) * VMTR_RGSH_pl(:,k,l))/rho_h_pl(:,k,l)
          end do
          wh_pl(:,ADM_kmin,l) = 0.0D0
          wh_pl(:,ADM_kmax+1,l) = 0.0D0
       end do
    end if
    !
    !
    !---- Courant number
!OCl SERIAL
    do l=1,ADM_lall
       ck(:,ADM_kmin-1,l,1) = 0.0D0
       ck(:,ADM_kmin-1,l,2) = 0.0D0
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax
          ck(:,k,l,1) =-wh(:,k,l)  &
               /GRD_dgz(k)*(dt*0.5D0) * rho_h(:,k,l)/(rhog(:,k,l)* VMTR_RGSGAM2(:,k,l))
          ck(:,k,l,2) = wh(:,k+1,l)  &
               /GRD_dgz(k)*(dt*0.5D0) * rho_h(:,k+1,l)/(rhog(:,k,l)* VMTR_RGSGAM2(:,k,l))
       end do
       ck(:,ADM_kmax+1,l,1) = 0.0D0
       ck(:,ADM_kmax+1,l,2) = 0.0D0
    end do
    !
    !
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1,ADM_lall_pl
          ck_pl(:,ADM_kmin-1,l,1) = 0.0D0
          ck_pl(:,ADM_kmin-1,l,2) = 0.0D0
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
             ck_pl(:,k,l,1) =-wh_pl(:,k,l)   &
                  /GRD_dgz(k)*(0.5D0*dt) * rho_h_pl(:,k,l)/(rhog_pl(:,k,l)* VMTR_RGSGAM2_pl(:,k,l))
             ck_pl(:,k,l,2) = wh_pl(:,k+1,l) &
                  /GRD_dgz(k)*(0.5D0*dt) * rho_h_pl(:,k+1,l)/(rhog_pl(:,k,l)* VMTR_RGSGAM2_pl(:,k,l))
          end do
          ck_pl(:,ADM_kmax+1,l,1) = 0.0D0
          ck_pl(:,ADM_kmax+1,l,2) = 0.0D0
       end do
    end if
    !
    do nq = 1, nqmax
       q(:,:,:,nq) = rhogq(:,:,:,nq)/rhog(:,:,:)
       if(ADM_prc_me==ADM_prc_pl) then
          q_pl(:,:,:,nq) = rhogq_pl(:,:,:,nq)/rhog_pl(:,:,:)
       end if
    end do
    !
    !--- inflow & outflow limiter
    do nq = 1, nqmax
       q_in_min(:,:,:,:,nq) = CNST_MAX_REAL
       q_in_max(:,:,:,:,nq) =-CNST_MAX_REAL

       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax+1
!OCL SIMD
             do n=1,ADM_gall
                if(ck(n,k,l,1)<=0.0D0) then
                   q_in_min(n,k,l,1,nq) =min(q(n,k,l,nq),q(n,k-1,l,nq))
                   q_in_max(n,k,l,1,nq) =max(q(n,k,l,nq),q(n,k-1,l,nq))
                else
                   q_in_min(n,k-1,l,2,nq) =min(q(n,k,l,nq),q(n,k-1,l,nq))
                   q_in_max(n,k-1,l,2,nq) =max(q(n,k,l,nq),q(n,k-1,l,nq))
                end if
             end do
             !
          end do
       end do
       !
       q_in_min_pl(:,:,:,:,nq) = CNST_MAX_REAL
       q_in_max_pl(:,:,:,:,nq) =-CNST_MAX_REAL
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
!OCL SIMD
                do n = 1, ADM_gall_pl
                   if(ck_pl(n,k,l,1)<=0.0D0) then
                      q_in_min_pl(n,k,l,1,nq) =min(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                      q_in_max_pl(n,k,l,1,nq) =max(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                   else
                      q_in_min_pl(n,k-1,l,2,nq) =min(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                      q_in_max_pl(n,k-1,l,2,nq) =max(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                   end if
                end do
                !
             end do
          end do
       end if
       !
       do l=1,ADM_lall
          do k = 1,ADM_kall
!OCL SIMD
             do n=1,ADM_gall
                !
                q_m1_k_min(n) = min(q_in_min(n,k,l,1,nq),q_in_min(n,k,l,2,nq))
                if(q_m1_k_min(n) == CNST_MAX_REAL) q_m1_k_min(n) = q(n,k,l,nq)
                q_m1_k_min(n) = max(SMALL_ZERO,q_m1_k_min(n))
                q_m1_k_max(n) = max(q_in_max(n,k,l,1,nq),q_in_max(n,k,l,2,nq))
                if(q_m1_k_max(n) ==-CNST_MAX_REAL) q_m1_k_max(n) = q(n,k,l,nq)
                !
                c_in_sum(n) &
                     = (0.5D0-sign(0.5D0,ck(n,k,l,1)))*ck(n,k,l,1)&
                     + (0.5D0-sign(0.5D0,ck(n,k,l,2)))*ck(n,k,l,2)
                c_out_sum(n) &
                     = (0.5D0+sign(0.5D0,ck(n,k,l,1)))*ck(n,k,l,1)&
                     + (0.5D0+sign(0.5D0,ck(n,k,l,2)))*ck(n,k,l,2)
                c_qin_sum_max(n) &
                     = (0.5D0-sign(0.5D0,ck(n,k,l,1)))*(ck(n,k,l,1)*q_in_max(n,k,l,1,nq))&
                     + (0.5D0-sign(0.5D0,ck(n,k,l,2)))*(ck(n,k,l,2)*q_in_max(n,k,l,2,nq))
                c_qin_sum_min(n) &
                     = (0.5D0-sign(0.5D0,ck(n,k,l,1)))*(ck(n,k,l,1)*q_in_min(n,k,l,1,nq))&
                     + (0.5D0-sign(0.5D0,ck(n,k,l,2)))*(ck(n,k,l,2)*q_in_min(n,k,l,2,nq))
                if(abs(c_out_sum(n))<CNST_EPS_ZERO) then
                   q_out_k_min(n,k,l,nq) = q(n,k,l,nq)
                   q_out_k_max(n,k,l,nq) = q(n,k,l,nq)
                else
                   q_out_k_min(n,k,l,nq) = (q(n,k,l,nq)-c_qin_sum_max(n)&
                        -q_m1_k_max(n)*(1.0D0-c_in_sum(n)-c_out_sum(n)))&
                        /c_out_sum(n)
                   q_out_k_max(n,k,l,nq) = (q(n,k,l,nq)-c_qin_sum_min(n)&
                        -q_m1_k_min(n)*(1.0D0-c_in_sum(n)-c_out_sum(n)))&
                        /c_out_sum(n)
                end if
             end do
          end do
       end do
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
             do k = 1,ADM_kall
!OCL SIMD
                do n=1,ADM_gall_pl
                   !
                   q_m1_k_min_pl(n) = min(q_in_min_pl(n,k,l,1,nq),q_in_min_pl(n,k,l,2,nq))
                   if(q_m1_k_min_pl(n) == CNST_MAX_REAL) q_m1_k_min_pl(n) = q_pl(n,k,l,nq)
                   q_m1_k_min_pl(n) = max(SMALL_ZERO,q_m1_k_min_pl(n))
                   q_m1_k_max_pl(n) = max(q_in_max_pl(n,k,l,1,nq),q_in_max_pl(n,k,l,2,nq))
                   if(q_m1_k_max_pl(n) ==-CNST_MAX_REAL) q_m1_k_max_pl(n) = q_pl(n,k,l,nq)
                   !
                   c_in_sum_pl(n) &
                        = (0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))*ck_pl(n,k,l,1)&
                        + (0.5D0-sign(0.5D0,ck_pl(n,k,l,2)))*ck_pl(n,k,l,2)
                   c_out_sum_pl(n) &
                        = (0.5D0+sign(0.5D0,ck_pl(n,k,l,1)))*ck_pl(n,k,l,1)&
                        + (0.5D0+sign(0.5D0,ck_pl(n,k,l,2)))*ck_pl(n,k,l,2)
                   c_qin_sum_max_pl(n) &
                        = (0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))*(ck_pl(n,k,l,1)*q_in_max_pl(n,k,l,1,nq))&
                        + (0.5D0-sign(0.5D0,ck_pl(n,k,l,2)))*(ck_pl(n,k,l,2)*q_in_max_pl(n,k,l,2,nq))
                   c_qin_sum_min_pl(n) &
                        = (0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))*(ck_pl(n,k,l,1)*q_in_min_pl(n,k,l,1,nq))&
                        + (0.5D0-sign(0.5D0,ck_pl(n,k,l,2)))*(ck_pl(n,k,l,2)*q_in_min_pl(n,k,l,2,nq))
                   if(abs(c_out_sum_pl(n))<CNST_EPS_ZERO) then
                      q_out_k_min_pl(n,k,l,nq) = q_pl(n,k,l,nq)
                      q_out_k_max_pl(n,k,l,nq) = q_pl(n,k,l,nq)
                   else
                      q_out_k_min_pl(n,k,l,nq) = (q_pl(n,k,l,nq)-c_qin_sum_max_pl(n)&
                           -q_m1_k_max_pl(n)*(1.0D0-c_in_sum_pl(n)-c_out_sum_pl(n)))&
                           /c_out_sum_pl(n)
                      q_out_k_max_pl(n,k,l,nq) = (q_pl(n,k,l,nq)-c_qin_sum_min_pl(n)&
                           -q_m1_k_min_pl(n)*(1.0D0-c_in_sum_pl(n)-c_out_sum_pl(n)))&
                           /c_out_sum_pl(n)
                   end if
                end do
             end do
          end do
       end if
    end do
    !
    !--- basic scheme
    do nq = 1, nqmax
       do l=1, ADM_lall 
          do k =ADM_kmin,ADM_kmax+1
             q_h(:,k,l,nq) = 0.5D0 *           &
                  ( GRD_afac(k) * q(:,k  ,l,nq)&
                  + GRD_bfac(k) * q(:,k-1,l,nq)&
                  )
          end do
       end do
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1, ADM_lall_pl
             do k =ADM_kmin,ADM_kmax+1
                q_h_pl(:,k,l,nq) = 0.5D0 *         &
                     ( GRD_afac(k) * q_pl(:,k  ,l,nq)&
                     + GRD_bfac(k) * q_pl(:,k-1,l,nq)&
                     )
             end do
          end do
       end if
    end do
    !
    !--- apply limiter
    !
    do nq = 1, nqmax
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             do n=1,ADM_gall 
                q_h(n,k,l,nq) &
                     =(0.5D0-sign(0.5D0,ck(n,k,l,1)))&
                     *min(max(q_h(n,k,l,nq),q_in_min(n,k,l,1,nq)),q_in_max(n,k,l,1,nq))&
                     +(0.5D0+sign(0.5D0,ck(n,k,l,1)))&
                     *min(max(q_h(n,k,l,nq),q_in_min(n,k-1,l,2,nq)),q_in_max(n,k-1,l,2,nq))
             end do
             do n=1,ADM_gall 
                q_h(n,k,l,nq) &
                     =(0.5D0-sign(0.5D0,ck(n,k,l,1)))&
                     *max(min(q_h(n,k,l,nq),q_out_k_max(n,k-1,l,nq)),q_out_k_min(n,k-1,l,nq))&
                     +(0.5D0+sign(0.5D0,ck(n,k,l,1)))&
                     *max(min(q_h(n,k,l,nq),q_out_k_max(n,k,l,nq)),q_out_k_min(n,k,l,nq))
             end do
          end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do n=1,ADM_gall_pl
                   q_h_pl(n,k,l,nq) &
                        =(0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))&
                        *min(max(q_h_pl(n,k,l,nq),q_in_min_pl(n,k,l,1,nq)),q_in_max_pl(n,k,l,1,nq))&
                        +(0.5D0+sign(0.5D0,ck_pl(n,k,l,1)))&
                        *min(max(q_h_pl(n,k,l,nq),q_in_min_pl(n,k-1,l,2,nq)),q_in_max_pl(n,k-1,l,2,nq))
                end do
                do n=1,ADM_gall_pl
                   q_h_pl(n,k,l,nq) &
                        =(0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))&
                        *max(min(q_h_pl(n,k,l,nq),q_out_k_max_pl(n,k-1,l,nq)),q_out_k_min_pl(n,k-1,l,nq))&
                        +(0.5D0+sign(0.5D0,ck_pl(n,k,l,1)))&
                        *max(min(q_h_pl(n,k,l,nq),q_out_k_max_pl(n,k,l,nq)),q_out_k_min_pl(n,k,l,nq))
                end do
             end do
          end do
       end if
    end do
    !
    !--- determine the flux and update
    do nq = 1, nqmax
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             flx_vz(:,k,l,nq) = rho_h(:,k,l)*wh(:,k,l) * q_h(:,k,l,nq)
          end do
          flx_vz(:,ADM_kmin,l,nq)   = 0.0D0
          flx_vz(:,ADM_kmax+1,l,nq) = 0.0D0
          do k = ADM_kmin, ADM_kmax
             rhogq(:,k,l,nq) = rhogq(:,k,l,nq) &
                  - ( flx_vz(:,k+1,l,nq) - flx_vz(:,k,l,nq) ) &
                  / ( GRD_gzh(k+1) - GRD_gzh(k) )*(dt*0.5D0)
          end do
       end do
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_LALL_PL
             do k = ADM_kmin, ADM_kmax+1
                flx_vz_pl(:,k,l,nq) = rho_h_pl(:,k,l)*wh_pl(:,k,l) * q_h_pl(:,k,l,nq)
             end do
             flx_vz_pl(:,ADM_kmin,l,nq)   = 0.0D0
             flx_vz_pl(:,ADM_kmax+1,l,nq) = 0.0D0
             do k = ADM_kmin, ADM_kmax
                rhogq_pl(:,k,l,nq) = rhogq_pl(:,k,l,nq) &
                     - ( flx_vz_pl(:,k+1,l,nq) - flx_vz_pl(:,k,l,nq) ) &
                     / ( GRD_gzh(k+1) - GRD_gzh(k) )*(dt*0.5D0)
             end do
          end do
       end if
    end do
    !
    !
    !--- Horizontal flux convergence
    hdiv(:,:,:,:) = 0.0D0
    hdiv_pl(:,:,:,:) = 0.0D0
    vx_r(:,:,:)=vx(:,:,:)*VMTR_RGAM(:,:,:)
    vy_r(:,:,:)=vy(:,:,:)*VMTR_RGAM(:,:,:)
    vz_r(:,:,:)=vz(:,:,:)*VMTR_RGAM(:,:,:)
    !
    if(ADM_prc_me==ADM_prc_pl) then
       vx_r_pl(:,:,:)=vx_pl(:,:,:)*VMTR_RGAM_pl(:,:,:)
       vy_r_pl(:,:,:)=vy_pl(:,:,:)*VMTR_RGAM_pl(:,:,:)
       vz_r_pl(:,:,:)=vz_pl(:,:,:)*VMTR_RGAM_pl(:,:,:)
    end if
    !
    call OPRT_divergence2_prep(&
         c,  c_pl,             &
         cp, cp_pl,            &
         vx_r, vx_r_pl,        &
         vy_r, vy_r_pl,        &
         vz_r, vz_r_pl,        &
         dt              &
         )
    !
!del Y.Niwa 080124    if ( save_memory ) then  ! K.Suzuki 07/04/25
    if (TRC_SAVE_MEMORY) then    ! Y.Niwa 080124
      do nq = 1, nqmax
        call OPRT_divergence2( &
             hdiv(:,:,:,nq),hdiv_pl(:,:,:,nq),     &
             rhogq(:,:,:,nq),rhogq_pl(:,:,:,nq),   &
             c, c_pl,                              &
             cp, cp_pl,                            &
             dt  )
      end do
    else
      call OPRT_divergence2_all( &
           hdiv,hdiv_pl,     &
           rhogq,rhogq_pl,   &
           c,  c_pl,         &
           cp, cp_pl,        &
           dt, nqmax )
    end if
    !
    !--- update rhogq by horizotal flux convergence
    do nq = 1, nqmax
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax
             rhogq(:,k,l,nq) = rhogq(:,k,l,nq) - hdiv(:,k,l,nq)*dt
          end do
       end do
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_LALL_PL
             do k = ADM_kmin, ADM_kmax
                rhogq_pl(:,k,l,nq) = rhogq_pl(:,k,l,nq) - hdiv_pl(:,k,l,nq)*dt
             end do
          end do
       end if
    end do
    !
    !--- vertical transport (fractional step) : 2nd 
    !
    do nq = 1, nqmax
       q(:,:,:,nq) = rhogq(:,:,:,nq)/rhog(:,:,:)
       if(ADM_prc_me==ADM_prc_pl) then
          q_pl(:,:,:,nq) = rhogq_pl(:,:,:,nq)/rhog_pl(:,:,:)
       end if
    end do
    !
    !--- inflow & outflow limiter
    do nq = 1, nqmax
       q_in_min(:,:,:,:,nq) = CNST_MAX_REAL
       q_in_max(:,:,:,:,nq) =-CNST_MAX_REAL
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax+1
!OCL SIMD
             do n=1,ADM_gall
                if(ck(n,k,l,1)<=0.0D0) then
                   q_in_min(n,k,l,1,nq) =min(q(n,k,l,nq),q(n,k-1,l,nq))
                   q_in_max(n,k,l,1,nq) =max(q(n,k,l,nq),q(n,k-1,l,nq))
                else
                   q_in_min(n,k-1,l,2,nq) =min(q(n,k,l,nq),q(n,k-1,l,nq))
                   q_in_max(n,k-1,l,2,nq) =max(q(n,k,l,nq),q(n,k-1,l,nq))
                end if
             end do
             !
          end do
       end do
       q_in_min_pl(:,:,:,:,nq) = CNST_MAX_REAL
       q_in_max_pl(:,:,:,:,nq) =-CNST_MAX_REAL
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
!OCL SIMD
                do n = 1, ADM_gall_pl
                   if(ck_pl(n,k,l,1)<=0.0D0) then
                      q_in_min_pl(n,k,l,1,nq) =min(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                      q_in_max_pl(n,k,l,1,nq) =max(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                   else
                      q_in_min_pl(n,k-1,l,2,nq) =min(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                      q_in_max_pl(n,k-1,l,2,nq) =max(q_pl(n,k,l,nq),q_pl(n,k-1,l,nq))
                   end if
                end do
                !
             end do
          end do
       end if

       do l=1,ADM_lall
          do k = 1,ADM_kall
!OCL SIMD
             do n=1,ADM_gall
                !
                q_m1_k_min(n) = min(q_in_min(n,k,l,1,nq),q_in_min(n,k,l,2,nq))
                if(q_m1_k_min(n) == CNST_MAX_REAL) q_m1_k_min(n) = q(n,k,l,nq)
                q_m1_k_min(n) = max(SMALL_ZERO,q_m1_k_min(n))
                q_m1_k_max(n) = max(q_in_max(n,k,l,1,nq),q_in_max(n,k,l,2,nq))
                if(q_m1_k_max(n) ==-CNST_MAX_REAL) q_m1_k_max(n) = q(n,k,l,nq)
                !
                c_in_sum(n) &
                     = (0.5D0-sign(0.5D0,ck(n,k,l,1)))*ck(n,k,l,1)&
                     + (0.5D0-sign(0.5D0,ck(n,k,l,2)))*ck(n,k,l,2)
                c_out_sum(n) &
                     = (0.5D0+sign(0.5D0,ck(n,k,l,1)))*ck(n,k,l,1)&
                     + (0.5D0+sign(0.5D0,ck(n,k,l,2)))*ck(n,k,l,2)
                c_qin_sum_max(n) &
                     = (0.5D0-sign(0.5D0,ck(n,k,l,1)))*(ck(n,k,l,1)*q_in_max(n,k,l,1,nq))&
                     + (0.5D0-sign(0.5D0,ck(n,k,l,2)))*(ck(n,k,l,2)*q_in_max(n,k,l,2,nq))
                c_qin_sum_min(n) &
                     = (0.5D0-sign(0.5D0,ck(n,k,l,1)))*(ck(n,k,l,1)*q_in_min(n,k,l,1,nq))&
                     + (0.5D0-sign(0.5D0,ck(n,k,l,2)))*(ck(n,k,l,2)*q_in_min(n,k,l,2,nq))
                if(abs(c_out_sum(n))<CNST_EPS_ZERO) then
                   q_out_k_min(n,k,l,nq) = q(n,k,l,nq)
                   q_out_k_max(n,k,l,nq) = q(n,k,l,nq)
                else
                   q_out_k_min(n,k,l,nq) = (q(n,k,l,nq)-c_qin_sum_max(n)&
                        -q_m1_k_max(n)*(1.0D0-c_in_sum(n)-c_out_sum(n)))&
                        /c_out_sum(n)
                   q_out_k_max(n,k,l,nq) = (q(n,k,l,nq)-c_qin_sum_min(n)&
                        -q_m1_k_min(n)*(1.0D0-c_in_sum(n)-c_out_sum(n)))&
                        /c_out_sum(n)
                end if
             end do
          end do
       end do
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
             do k = 1,ADM_kall
!OCL SIMD
                do n=1,ADM_gall_pl
                   !
                   q_m1_k_min_pl(n) = min(q_in_min_pl(n,k,l,1,nq),q_in_min_pl(n,k,l,2,nq))
                   if(q_m1_k_min_pl(n) == CNST_MAX_REAL) q_m1_k_min_pl(n) = q_pl(n,k,l,nq)
                   q_m1_k_min_pl(n) = max(SMALL_ZERO,q_m1_k_min_pl(n))
                   q_m1_k_max_pl(n) = max(q_in_max_pl(n,k,l,1,nq),q_in_max_pl(n,k,l,2,nq))
                   if(q_m1_k_max_pl(n) ==-CNST_MAX_REAL) q_m1_k_max_pl(n) = q_pl(n,k,l,nq)
                   !
                   c_in_sum_pl(n) &
                        = (0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))*ck_pl(n,k,l,1)&
                        + (0.5D0-sign(0.5D0,ck_pl(n,k,l,2)))*ck_pl(n,k,l,2)
                   c_out_sum_pl(n) &
                        = (0.5D0+sign(0.5D0,ck_pl(n,k,l,1)))*ck_pl(n,k,l,1)&
                        + (0.5D0+sign(0.5D0,ck_pl(n,k,l,2)))*ck_pl(n,k,l,2)
                   c_qin_sum_max_pl(n) &
                        = (0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))*(ck_pl(n,k,l,1)*q_in_max_pl(n,k,l,1,nq))&
                        + (0.5D0-sign(0.5D0,ck_pl(n,k,l,2)))*(ck_pl(n,k,l,2)*q_in_max_pl(n,k,l,2,nq))
                   c_qin_sum_min_pl(n) &
                        = (0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))*(ck_pl(n,k,l,1)*q_in_min_pl(n,k,l,1,nq))&
                        + (0.5D0-sign(0.5D0,ck_pl(n,k,l,2)))*(ck_pl(n,k,l,2)*q_in_min_pl(n,k,l,2,nq))
                   if(abs(c_out_sum_pl(n))<CNST_EPS_ZERO) then
                      q_out_k_min_pl(n,k,l,nq) = q_pl(n,k,l,nq)
                      q_out_k_max_pl(n,k,l,nq) = q_pl(n,k,l,nq)
                   else
                      q_out_k_min_pl(n,k,l,nq) = (q_pl(n,k,l,nq)-c_qin_sum_max_pl(n)&
                           -q_m1_k_max_pl(n)*(1.0D0-c_in_sum_pl(n)-c_out_sum_pl(n)))&
                           /c_out_sum_pl(n)
                      q_out_k_max_pl(n,k,l,nq) = (q_pl(n,k,l,nq)-c_qin_sum_min_pl(n)&
                           -q_m1_k_min_pl(n)*(1.0D0-c_in_sum_pl(n)-c_out_sum_pl(n)))&
                           /c_out_sum_pl(n)
                   end if
                end do
             end do
          end do
       end if
    end do
    !
    !--- basic scheme
    do nq = 1, nqmax
       do l=1, ADM_lall 
          do k =ADM_kmin,ADM_kmax+1
             q_h(:,k,l,nq) = 0.5D0 *           &
                  ( GRD_afac(k) * q(:,k  ,l,nq)&
                  + GRD_bfac(k) * q(:,k-1,l,nq)&
                  )
          end do
       end do
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1, ADM_lall_pl
             do k =ADM_kmin,ADM_kmax+1
                q_h_pl(:,k,l,nq) = 0.5D0 *         &
                     ( GRD_afac(k) * q_pl(:,k  ,l,nq)&
                     + GRD_bfac(k) * q_pl(:,k-1,l,nq)&
                     )
             end do
          end do
       end if
    end do
    !
    !--- apply limiter
    !
    !
    do nq = 1, nqmax
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             do n=1,ADM_gall 
                q_h(n,k,l,nq) &
                     =(0.5D0-sign(0.5D0,ck(n,k,l,1)))&
                     *min(max(q_h(n,k,l,nq),q_in_min(n,k,l,1,nq)),q_in_max(n,k,l,1,nq))&
                     +(0.5D0+sign(0.5D0,ck(n,k,l,1)))&
                     *min(max(q_h(n,k,l,nq),q_in_min(n,k-1,l,2,nq)),q_in_max(n,k-1,l,2,nq))
             end do
             do n=1,ADM_gall 
                q_h(n,k,l,nq) &
                     =(0.5D0-sign(0.5D0,ck(n,k,l,1)))&
                     *max(min(q_h(n,k,l,nq),q_out_k_max(n,k-1,l,nq)),q_out_k_min(n,k-1,l,nq))&
                     +(0.5D0+sign(0.5D0,ck(n,k,l,1)))&
                     *max(min(q_h(n,k,l,nq),q_out_k_max(n,k,l,nq)),q_out_k_min(n,k,l,nq))
             end do
          end do
       end do

       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
                do n=1,ADM_gall_pl
                   q_h_pl(n,k,l,nq) &
                        =(0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))&
                        *min(max(q_h_pl(n,k,l,nq),q_in_min_pl(n,k,l,1,nq)),q_in_max_pl(n,k,l,1,nq))&
                        +(0.5D0+sign(0.5D0,ck_pl(n,k,l,1)))&
                        *min(max(q_h_pl(n,k,l,nq),q_in_min_pl(n,k-1,l,2,nq)),q_in_max_pl(n,k-1,l,2,nq))
                end do
                do n=1,ADM_gall_pl
                   q_h_pl(n,k,l,nq) &
                        =(0.5D0-sign(0.5D0,ck_pl(n,k,l,1)))&
                        *max(min(q_h_pl(n,k,l,nq),q_out_k_max_pl(n,k-1,l,nq)),q_out_k_min_pl(n,k-1,l,nq))&
                        +(0.5D0+sign(0.5D0,ck_pl(n,k,l,1)))&
                        *max(min(q_h_pl(n,k,l,nq),q_out_k_max_pl(n,k,l,nq)),q_out_k_min_pl(n,k,l,nq))
                end do
             end do
          end do
       end if
    end do
    !
    !
    !--- determine the flux and update
    do nq = 1, nqmax
       do l=1,ADM_lall
          do k = ADM_kmin, ADM_kmax+1
             flx_vz(:,k,l,nq) = rho_h(:,k,l)*wh(:,k,l) * q_h(:,k,l,nq)
          end do
          flx_vz(:,ADM_kmin,l,nq)   = 0.0D0
          flx_vz(:,ADM_kmax+1,l,nq) = 0.0D0
          do k = ADM_kmin, ADM_kmax
             rhogq(:,k,l,nq) = rhogq(:,k,l,nq) &
                  - ( flx_vz(:,k+1,l,nq) - flx_vz(:,k,l,nq) ) &
                  / ( GRD_gzh(k+1) - GRD_gzh(k) )*(dt*0.5D0)
          end do
       end do
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1,ADM_LALL_PL
             do k = ADM_kmin, ADM_kmax+1
                flx_vz_pl(:,k,l,nq) = rho_h_pl(:,k,l)*wh_pl(:,k,l) * q_h_pl(:,k,l,nq)
             end do
             flx_vz_pl(:,ADM_kmin,l,nq)   = 0.0D0
             flx_vz_pl(:,ADM_kmax+1,l,nq) = 0.0D0
             do k = ADM_kmin, ADM_kmax
                rhogq_pl(:,k,l,nq) = rhogq_pl(:,k,l,nq) &
                     - ( flx_vz_pl(:,k+1,l,nq) - flx_vz_pl(:,k,l,nq) ) &
                     / ( GRD_gzh(k+1) - GRD_gzh(k) )*(dt*0.5D0)
             end do
          end do
       end if
    end do
    !
    !
   !-->
!    do k=ADM_kmin,ADM_kmax
!       rhogqv_max(k) = GTL_max_k(rhogq(:,:,:,I_QV),rhogq_pl(:,:,:,I_QV),k)
!       rhogqv_min(k) = GTL_min_k(rhogq(:,:,:,I_QV),rhogq_pl(:,:,:,I_QV),k)
!       rhogqc_max(k) = GTL_max_k(rhogq(:,:,:,I_QC),rhogq_pl(:,:,:,I_QC),k)
!       rhogqc_min(k) = GTL_min_k(rhogq(:,:,:,I_QC),rhogq_pl(:,:,:,I_QC),k)
!    end do
!    if(ADM_prc_me==ADM_prc_run_master) then
!       write(ADM_LOG_FID,*) '------------- end of src_update_tracer'
!       do k = ADM_kmin, ADM_kmax
!          write(ADM_LOG_FID,'(I4,4E14.6)') k, rhogqv_max(k), rhogqv_min(k), rhogqc_max(k), rhogqc_min(k)
!       end do
!    end if
!    !<---


  end subroutine src_update_tracer
  !----------------------------------------------------------------------------------
  ! Y.Niwa add 080124
  ! This routine is revised version of src_update_tracer
  subroutine src_update_tracer_rev( &
       nqmax,                       & !--- IN    : number of tracers
       rhogq,       rhogq_pl,       & !--- INOUT : rhogq   ( gam2 X G^{1/2} )
       rhog_in,     rhog_in_pl,     & !--- IN    : rho(old)( gam2 X G^{1/2} )
       rhog_mean,   rhog_mean_pl,   & !--- IN    : rho     ( gam2 X G^{1/2} )
       rhogvx_mean, rhogvx_mean_pl, & !--- IN    : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy_mean, rhogvy_mean_pl, & !--- IN    : rho*Vy  ( gam2 X G^{1/2} )
       rhogvz_mean, rhogvz_mean_pl, & !--- IN    : rho*Vz  ( gam2 X G^{1/2} )
       rhogw_mean,  rhogw_mean_pl,  & !--- IN    : rho*w   ( gam2 X G^{1/2} )
       frhog,       frhog_pl,       & !--- IN    : hyperviscosity tendency for rhog
       dt                           ) !--- IN    : delta t
    use mod_adm, only :  &
         ADM_LOG_FID,    &
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
         ADM_AI,ADM_AJ,ADM_AIJ,&
         ADM_KNONE,      &
         ADM_GSLF_PL,    &
         ADM_GMIN_PL,    &
         ADM_GMAX_PL,    &
         ADM_gall_1d,    &
         ADM_prc_tab
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
         VMTR_RGAM_pl,   &
         VMTR_GSGAM2H,   &
         VMTR_GSGAM2H_pl
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_rdgz
    use mod_oprt, only: &
       oprt_divergence2_prep_rev, &
       oprt_divergence2_rev
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

#ifdef _FJTIMER_
call timer_sta(2600)
call timer_sta(2601)
#endif

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 1st
    !---------------------------------------------------------------------------
!OCL SERIAL
    do l = 1, ADM_lall

!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog_in(g,k,l)
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             d(g,k,l) = b1 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo
       enddo

!OCL SERIAL
       do k = ADM_kmin+1, ADM_kmax
!OCL PARALLEL
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

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          flx_v(g,ADM_kmin,  l) = 0.D0
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          flx_v(g,ADM_kmax+1,l) = 0.D0
       enddo

       !---- Courant numbers at cell boundary
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             ck(g,k,l,1) = -flx_v(g,k  ,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             ck(g,k,l,2) =  flx_v(g,k+1,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
       enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,1) = 0.D0 
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,2) = 0.D0
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmax+1,l,1) = 0.D0
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmax+1,l,2) = 0.D0
       enddo

    enddo ! l LOOP

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl

!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_in_pl(g,k,l)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                d_pl(g,k,l) = b1 * frhog_pl(g,k,l) * dt * rrhog_pl(g,k,l)
             enddo
          enddo

!OCL SERIAL
          do k = ADM_kmin+1, ADM_kmax
!OCL PARALLEL
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

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             flx_v_pl(g,ADM_kmin,  l) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             flx_v_pl(g,ADM_kmax+1,l) = 0.D0
          enddo

          !---- Courant numbers at cell boundary
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,1) = -flx_v_pl(g,k  ,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,2) =  flx_v_pl(g,k+1,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
          enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,1) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,2) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmax+1,l,1) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmax+1,l,2) = 0.D0
          enddo

       enddo
    endif

#ifdef _FJTIMER_
call timer_end(2601)
call timer_sta(2602)
#endif

    !--- basic scheme ( 2nd-order centered difference )
!OCL SERIAL
    do nq = 1, nqmax

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo

!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                q_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * q(g,k,  l) &
                                     + GRD_bfac(k) * q(g,k-1,l) )
             enddo
          enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             q_h(g,ADM_kmin-1,l) = 0.D0
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo

!OCL PARALLEL
             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   q_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(g,k,  l) &
                                           + GRD_bfac(k) * q_pl(g,k-1,l) )
                enddo
             enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin-1,l) = 0.D0
             enddo
          enddo
       endif

#ifdef _FJTIMER_
call timer_sta(2609)
#endif
       call advlim_thuburn_v( q_h, q_h_pl, & !--- [INOUT]
                              q,   q_pl,   & !--- [IN]
                              ck,  ck_pl,  & !--- [IN]
                              d,   d_pl    ) !--- [IN]
#ifdef _FJTIMER_
call timer_end(2609)
#endif

       !--- update rhogq
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             q_h(g,ADM_kmin  ,l) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             q_h(g,ADM_kmax+1,l) = 0.D0
          enddo

!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) &
                                - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                  - flx_v(g,k,  l) * q_h(g,k,  l) ) * GRD_rdgz(k)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin  ,l) = 0.D0
             enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmax+1,l) = 0.D0
             enddo

!OCL PARALLEL
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

#ifdef _FJTIMER_
call timer_end(2602)
call timer_sta(2603)
#endif

    !--- update rhog
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             rhog(g,k,l) = rhog_in(g,k,l)                                  &
                         - ( flx_v(g,k+1,l) - flx_v(g,k,l) ) * GRD_rdgz(k) &
                         + b1 * frhog(g,k,l) * dt
          enddo
       enddo

!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          rhog(g,ADM_kmin-1,l) = rhog_in(g,ADM_kmin,l)
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          rhog(g,ADM_kmax+1,l) = rhog_in(g,ADM_kmax,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall_pl
                rhog_pl(g,k,l) = rhog_in_pl(g,k,l)                                     &
                               - ( flx_v_pl(g,k+1,l) - flx_v_pl(g,k,l) ) * GRD_rdgz(k) &
                               + b1 * frhog_pl(g,k,l) * dt
             enddo
          enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             rhog_pl(g,ADM_kmin-1,l) = rhog_in_pl(g,ADM_kmin,l)
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             rhog_pl(g,ADM_kmax+1,l) = rhog_in_pl(g,ADM_kmax,l)
          enddo
       enddo
    endif

#ifdef _FJTIMER_
call timer_end(2603)
call timer_sta(2604)
#endif

    !---------------------------------------------------------------------------
    ! Horizontal advection by MIURA 2004 scheme
    !---------------------------------------------------------------------------
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             vx_r(g,k,l) = rhogvx_mean(g,k,l) * VMTR_RGAM(g,k,l)
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             vy_r(g,k,l) = rhogvy_mean(g,k,l) * VMTR_RGAM(g,k,l) 
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             vz_r(g,k,l) = rhogvz_mean(g,k,l) * VMTR_RGAM(g,k,l) 
          enddo
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                vx_r_pl(g,k,l) = rhogvx_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)  
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                vy_r_pl(g,k,l) = rhogvy_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)   
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                vz_r_pl(g,k,l) = rhogvz_mean_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)   
             enddo
          enddo
       enddo
    endif

#ifdef _FJTIMER_
call timer_sta(2610)
#endif
    call OPRT_divergence2_prep_rev( flx_h,     flx_h_pl,     & !--- [OUT]
                                    cp,        cp_pl,        & !--- [OUT]
                                    vx_r,      vx_r_pl,      & !--- [IN]
                                    vy_r,      vy_r_pl,      & !--- [IN]
                                    vz_r,      vz_r_pl,      & !--- [IN]
                                    rhog_mean, rhog_mean_pl, & !--- [IN]
                                    dt                       ) !--- [IN]
#ifdef _FJTIMER_
call timer_end(2610)
#endif

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall

!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog(g,k,l)
          enddo

          !--- Courant number
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             c(1,g,k,l) = flx_h(1,g,k,l) * rrhog(g,k,l)
             c(2,g,k,l) = flx_h(2,g,k,l) * rrhog(g,k,l)
             c(3,g,k,l) = flx_h(3,g,k,l) * rrhog(g,k,l)
             c(4,g,k,l) = flx_h(4,g,k,l) * rrhog(g,k,l)
             c(5,g,k,l) = flx_h(5,g,k,l) * rrhog(g,k,l)
             c(6,g,k,l) = flx_h(6,g,k,l) * rrhog(g,k,l)
          enddo

!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             d(g,k,l) = b2 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo

       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall

!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_pl(g,k,l)
             enddo

!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do n = ADM_GMIN_PL, ADM_GMAX_PL
                c_pl(n,k,l) = flx_h_pl(n,k,l) * rrhog_pl(ADM_GSLF_PL,k,l)
             enddo

             d_pl(ADM_GSLF_PL,k,l) = b2 * frhog_pl(ADM_GSLF_PL,k,l) * dt * rrhog_pl(ADM_GSLF_PL,k,l)
          enddo
       enddo
    endif

#ifdef _FJTIMER_
call timer_end(2604)
call timer_sta(2605)
#endif

!OCL SERIAL
    do nq = 1, nqmax

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo
          enddo
       endif

#ifdef _FJTIMER_
call timer_sta(2611)
#endif
       call OPRT_divergence2_rev( hdiv,  hdiv_pl,  & !--- [OUT]
                                  q,     q_pl,     & !--- [IN]
                                  flx_h, flx_h_pl, & !--- [IN]
                                  c,     c_pl,     & !--- [IN]
                                  cp,    cp_pl,    & !--- [IN]
                                  d,     d_pl,     & !--- [IN]
                                  dt               ) !--- [IN]
#ifdef _FJTIMER_
call timer_end(2611)
#endif

       !--- update rhogq
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) - hdiv(g,k,l) * dt
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   rhogq_pl(g,k,l,nq) = rhogq_pl(g,k,l,nq) - hdiv_pl(g,k,l) * dt
                enddo
             enddo
          enddo
       endif

    enddo ! tracer q LOOP

#ifdef _FJTIMER_
call timer_end(2605)
call timer_sta(2606)
#endif

    !--- update rhog
!OCL SERIAL
    do l = 1, ADM_lall
       nstart = suf(ADM_gmin,ADM_gmin)
       nend   = suf(ADM_gmax,ADM_gmax)

!OCL PARALLEL
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

!OCL PARALLEL
       do k = 1, ADM_kall
          rhog(suf(ADM_gall_1d,1),k,l) = rhog(suf(ADM_gmax+1,ADM_gmin),k,l)
          rhog(suf(1,ADM_gall_1d),k,l) = rhog(suf(ADM_gmin,ADM_gmax+1),k,l)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_GSLF_PL

!OCL SERIAL
       do l = 1, ADM_lall_pl
!OCL PARALLEL
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

#ifdef _FJTIMER_
call timer_end(2606)
call timer_sta(2607)
#endif

    !---------------------------------------------------------------------------
    ! Vertical Advection (fractioanl step) : 2nd
    !---------------------------------------------------------------------------
!OCL SERIAL
    do l = 1, ADM_lall

!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             rrhog(g,k,l) = 1.D0 / rhog(g,k,l)
          enddo
       enddo

       !---- Courant numbers at cell boundary
!OCL PARALLEL
       do k = ADM_kmin, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             ck(g,k,l,1) = -flx_v(g,k  ,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             ck(g,k,l,2) =  flx_v(g,k+1,l) * rrhog(g,k,l) * GRD_rdgz(k)
          enddo
       enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,1) = 0.D0 
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmin-1,l,2) = 0.D0
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmax+1,l,1) = 0.D0
       enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
       do g = 1, ADM_gall
          ck(g,ADM_kmax+1,l,2) = 0.D0
       enddo

!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             d(g,k,l) = b3 * frhog(g,k,l) * dt * rrhog(g,k,l)
          enddo
       enddo

    enddo ! l LOOP

    if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
       do l = 1, ADM_lall_pl

!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                rrhog_pl(g,k,l) = 1.D0 / rhog_pl(g,k,l)
             enddo
          enddo

          !---- Courant numbers at cell boundary
!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,1) = -flx_v_pl(g,k  ,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                ck_pl(g,k,l,2) =  flx_v_pl(g,k+1,l) * rrhog_pl(g,k,l) * GRD_rdgz(k)
             enddo
          enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,1) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmin-1,l,2) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmax+1,l,1) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall_pl
             ck_pl(g,ADM_kmax+1,l,2) = 0.D0
          enddo

!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                d_pl(g,k,l) = b3 * frhog_pl(g,k,l) * dt * rrhog_pl(g,k,l)
             enddo
          enddo

       enddo
    endif

#ifdef _FJTIMER_
call timer_end(2607)
call timer_sta(2608)
#endif

    !--- basic scheme ( 2nd-order centered difference )
!OCL SERIAL
    do nq = 1, nqmax

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall
                q(g,k,l) = rhogq(g,k,l,nq) * rrhog(g,k,l)
             enddo
          enddo

!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall
                q_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * q(g,k,  l) &
                                     + GRD_bfac(k) * q(g,k-1,l) )
             enddo
          enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             q_h(g,ADM_kmin-1,l) = 0.D0
          enddo
       enddo


       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
             do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
                do g = 1, ADM_gall_pl
                   q_pl(g,k,l) = rhogq_pl(g,k,l,nq) * rrhog_pl(g,k,l)
                enddo
             enddo

!OCL PARALLEL
             do k = ADM_kmin, ADM_kmax+1
                do g = 1, ADM_gall_pl
                   q_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * q_pl(g,k,  l) &
                                           + GRD_bfac(k) * q_pl(g,k-1,l) )
                enddo
             enddo

!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin-1,l) = 0.D0
             enddo
          enddo
       endif

#ifdef _FJTIMER_
call timer_sta(2609)
#endif
       call advlim_thuburn_v( q_h, q_h_pl, & !--- [INOUT]
                              q,   q_pl,   & !--- [IN]
                              ck,  ck_pl,  & !--- [IN]
                              d,   d_pl    ) !--- [IN]
#ifdef _FJTIMER_
call timer_end(2609)
#endif

       !--- update rhogq
!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             q_h(g,ADM_kmin  ,l) = 0.D0
          enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do g = 1, ADM_gall
             q_h(g,ADM_kmax+1,l) = 0.D0
          enddo

!OCL PARALLEL
          do k = ADM_kmin, ADM_kmax
             do g = 1, ADM_gall
                rhogq(g,k,l,nq) = rhogq(g,k,l,nq) &
                                - ( flx_v(g,k+1,l) * q_h(g,k+1,l) &
                                  - flx_v(g,k,  l) * q_h(g,k,  l) ) * GRD_rdgz(k)
             enddo
          enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
!OCL SERIAL
          do l = 1, ADM_lall_pl
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmin  ,l) = 0.D0
             enddo
!OCL PARALLEL
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
             do g = 1, ADM_gall_pl
                q_h_pl(g,ADM_kmax+1,l) = 0.D0
             enddo

!OCL PARALLEL
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

#ifdef _FJTIMER_
call timer_end(2608)
call timer_end(2600)
#endif

    return
  end subroutine src_update_tracer_rev

end module mod_src
