!-------------------------------------------------------------------------------
!>
!! Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)    Imported from igdc-4.33
!! @li      2006-04-17 (H.Tomita)    Add sub[OPRT_divergence]
!! @li      2006-08-11 (H.Tomita)    Implementatio of miura scheme(2004) with Thurbuen(1996)'s limiter.
!!                                   Add sub[OPRT_divergence2]
!!                                       sub[OPRT_divergence2_prep]
!!                                       sub[OPRT_divergence2_all].
!! @li      2006-08-22 (Y.Niwa)      divide the rows for the calc. of clap and clap_pl due to the rule of XLF.
!! @li      2007-01-26 (H.Tomita)    Optimization of sub[oprt_diffusion].
!! @li      2007-11-28 (T.Mitsui)    bugfix in oprt_divergence2, _all
!! @li      2008-01-24 (Y.Niwa)      add OPRT_divergence2{_prep,,_all}_rev
!! @li      2008-04-28 (T.Mitsui)    bug fix in OPRT_divergence2{_all}_rev
!! @li      2009-09-04 (H.Taniguchi) bug fix in OPRT_divergence2{_all}_rev Zero clear of wrk[_pl] is needed.
!! @li      2010-06-08 (S.Iga)       new grid is implemented (see, string XTMS)
!! @li      2011-09-27 (T.Seiki)     merge optimization by RIST and M.Terai
!! @li      2012-01-17 (M.Terai)     update optimization(case6) in div2rev
!! @li      2012-05-01 (T.Yamaura)   bug fix in div2rev
!! @li      2012-06-28 (M.Terai)     Removed wrapper subroutine to invoked the operators directly macro var.
!!
!<
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,    &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
     K0  => ADM_KNONE
  use mod_gmtr, only: &
     P_RAREA => GMTR_P_RAREA, &
     T_RAREA => GMTR_T_RAREA, &
     W1      => GMTR_T_W1,    &
     W2      => GMTR_T_W2,    &
     W3      => GMTR_T_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
     HTX     => GMTR_A_HTX,   &
     HTY     => GMTR_A_HTY,   &
     HTZ     => GMTR_A_HTZ,   &
     TNX     => GMTR_A_TNX,   &
     TNY     => GMTR_A_TNY,   &
     TNZ     => GMTR_A_TNZ,   &
     TN2X    => GMTR_A_TN2X,  &
     TN2Y    => GMTR_A_TN2Y,  &
     TN2Z    => GMTR_A_TN2Z
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OPRT_setup
  public :: OPRT_divergence
  public :: OPRT_gradient
  public :: OPRT_laplacian
  public :: OPRT_diffusion
  public :: OPRT_horizontalize_vec
  public :: OPRT_vorticity
  public :: OPRT_divdamp

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, save :: OPRT_nstart
  integer, public, save :: OPRT_nend

  ! < for divergence operator >
  real(8), public, allocatable, save :: cdiv   (:,:,:,:)
  real(8), public, allocatable, save :: cdiv_pl(:,:,:,:)

  ! < for gradient operator >
  real(8), public, allocatable, save :: cgrad   (:,:,:,:)
  real(8), public, allocatable, save :: cgrad_pl(:,:,:,:)

  ! < for laplacian operator >
  real(8), public, allocatable, save :: clap   (:,:,:)
  real(8), public, allocatable, save :: clap_pl(:,:,:)

  ! < for diffusion operator >
  real(8), public, allocatable, save :: cinterp_TN (:,:,:,:)
  real(8), public, allocatable, save :: cinterp_HN (:,:,:,:)
  real(8), public, allocatable, save :: cinterp_TRA(:,:,:)
  real(8), public, allocatable, save :: cinterp_PRA(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OPRT_setup
    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_vlink_nmax, &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl
    use mod_gmtr, only: &
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    implicit none

    integer :: n0,n1,n2,n3,n4

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, l, m, md, v

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[oprt]/Category[common share]'

    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)

    !---< setup coefficient of divergence operator >
    write(ADM_LOG_FID,*) '*** setup coefficient of divergence operator'

    allocate( cdiv   (0:6,             ADM_gall   ,ADM_lall   ,3) )
    allocate( cdiv_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl,3) )

    do l = 1, ADM_lall

       do m = 1, 3
          md = m + HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cdiv(0,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1j
             cdiv(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1jp1
             cdiv(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijp1
             cdiv(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0*GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1j
             cdiv(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1jm1
             cdiv(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijm1
             cdiv(6,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + HNX - 1

             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cdiv(0,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
             ! ip1j
             cdiv(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
             ! ip1jp1
             cdiv(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
             ! ijp1
             cdiv(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
             ! im1j
             cdiv(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
             ! im1jm1
             cdiv(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
             ! ijm1
             cdiv(6,n,l,m) = ( + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,P_RAREA)
          enddo
       endif

    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + HNX - 1

             cdiv_pl(0,n,l,m) = 0.D0
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cdiv_pl(0,n,l,m) = cdiv_pl(0,n,l,m) + ( GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                                      + GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ijp1,k0,l,md) )
             enddo
             cdiv_pl(0,n,l,m) = cdiv_pl(0,n,l,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cdiv_pl(v-1,n,l,m) = ( + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                       + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                     ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif



    !---< setup coefficient of gradient operator >

    write(ADM_LOG_FID,*) '*** setup coefficient of gradient operator'

    allocate( cgrad   (0:6,             ADM_gall   ,ADM_lall   ,3) )
    allocate( cgrad_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl,3) )

    do l = 1, ADM_lall

       do m = 1, 3
          md = m + HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cgrad(0,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AI ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AIJ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AJ ,md)                          &
                                + 2.D0 * GMTR_A_var(im1j  ,k0,l,AI ,md)                          &
                                + 2.D0 * GMTR_A_var(im1jm1,k0,l,AIJ,md)                          &
                                + 2.D0 * GMTR_A_var(ijm1  ,k0,l,AJ ,md)                          &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1j
             cgrad(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1jp1
             cgrad(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijp1
             cgrad(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1j
             cgrad(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1jm1
             cgrad(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijm1
             cgrad(6,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + HNX - 1

             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cgrad(0,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AI ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AIJ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AJ ,md)                          &
                                + 2.D0 * GMTR_A_var(im1j  ,k0,l,AI ,md)                          &
                                + 2.D0 * GMTR_A_var(im1jm1,k0,l,AIJ,md)                          &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1j
             cgrad(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1jp1
             cgrad(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijp1
             cgrad(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1j
             cgrad(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1jm1
             cgrad(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijm1
             cgrad(6,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,P_RAREA)
          enddo
       endif
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + HNX - 1

             cgrad_pl(0,n,l,m) = 0.D0
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cgrad_pl(0,n,l,m) = cgrad_pl(0,n,l,m) &
                                  + 2.D0 * ( GMTR_T_var_pl(ij,k0,l,W1) - 1.D0 ) * GMTR_A_var_pl(ijp1,k0,l,md)
             enddo
             cgrad_pl(0,n,l,m) = cdiv_pl(0,n,l,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cgrad_pl(v-1,n,l,m) = ( + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                        + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                        + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                        + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                      ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif

    ! ---- setup coefficient of laplacian operator

    write(ADM_LOG_FID,*) '*** setup coefficient of laplacian operator'

    allocate( clap   (0:6,             ADM_gall   ,ADM_lall   ) )
    allocate( clap_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl) )

    do l = 1, ADM_lall

       do n = OPRT_nstart, OPRT_nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          ! ij
          clap(0,ij,l) = ( &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          clap(0,ij,l) = clap(0,ij,l) + ( &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ip1j
          clap(1,ij,l) = ( &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ip1jp1
          clap(2,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ijp1
          clap(3,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! im1j
          clap(4,ij,l) = ( &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! im1jm1
          clap(5,ij,l) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ijm1
          clap(6,ij,l) = ( &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +2*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +2*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +2*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          ! 0: ij
          clap(0,ij,l) = ( &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0  ! Y.Niwa add 06/08/22

          clap(0,ij,l) = clap(0,ij,l) + ( &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ip1j
          clap(1,ij,l) = ( &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ip1jp1
          clap(2,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ijp1
          clap(3,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! im1j
          clap(4,ij,l) = ( &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! im1jm1
          clap(5,ij,l) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

          ! ijm1
          clap(6,ij,l) = ( &
      -1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
      +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
      -2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
      -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
      +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
      -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
      +1*GMTR_A_var(ijm1,k0,l,AJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
      -2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
      -2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
         )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0d0

      endif
    enddo

    if ( ADM_have_pl ) then

      n =ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4

      do l = 1,ADM_lall_pl
          clap_pl(0,n,l)=( &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0   ! Y.Niwa add 060822

          clap_pl(0,n,l)= clap_pl(0,n,l) + ( &                        ! Y.Niwa add 060822
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0
          !
          ! n0
          clap_pl(1,n,l)=( &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0
          !
          ! n1
          clap_pl(2,n,l)=( &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0
          !
          ! n2
          clap_pl(3,n,l)=( &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0
          !
          ! n3
          clap_pl(4,n,l)=( &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0
          !
          ! n4
          clap_pl(5,n,l)=( &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0d0
      enddo
    endif




    ! ---- setup coefficient of diffusion operator

    write(ADM_LOG_FID,*) '*** setup coefficient of diffusion operator'

    allocate( cinterp_TN(AI:AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cinterp_HN(AI:AJ,1:3,ADM_gall,ADM_lall) )

    allocate( cinterp_TRA(TI:TJ,ADM_gall,ADM_lall) )
    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )

    do m = AI, AJ
    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cinterp_TN(m,1,n,l) = GMTR_A_var(n,k0,l,m,TNX)
       cinterp_TN(m,2,n,l) = GMTR_A_var(n,k0,l,m,TNY)
       cinterp_TN(m,3,n,l) = GMTR_A_var(n,k0,l,m,TNZ)

       cinterp_HN(m,1,n,l) = GMTR_A_var(n,k0,l,m,HNX)
       cinterp_HN(m,2,n,l) = GMTR_A_var(n,k0,l,m,HNY)
       cinterp_HN(m,3,n,l) = GMTR_A_var(n,k0,l,m,HNZ)
    enddo
    enddo
    enddo

    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cinterp_TRA(TI,n,l) = GMTR_T_var(n,k0,l,TI,T_RAREA)
       cinterp_TRA(TJ,n,l) = GMTR_T_var(n,k0,l,TJ,T_RAREA)

       cinterp_PRA(n,l) = GMTR_P_var(n,k0,l,P_RAREA)
    enddo
    enddo

    return
  end subroutine OPRT_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence( &
       scl, scl_pl, &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl,  &
       mfact        )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(8), intent(out) :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: mfact

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_divergence')

    !$acc data present_or_copyin(cdiv,vx,vy,vz) copyout(scl)
    !$acc kernels
    !$acc loop seq
    do l = 1, ADM_lall
    !$acc loop gang vector(32)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d
       im1j   = n - 1
       ijm1   = n     - ADM_gall_1d
       im1jm1 = n - 1 - ADM_gall_1d

       scl(n,k,l) = ( cdiv(0,n,l,1) * vx(ij    ,k,l) &
                    + cdiv(1,n,l,1) * vx(ip1j  ,k,l) &
                    + cdiv(2,n,l,1) * vx(ip1jp1,k,l) &
                    + cdiv(3,n,l,1) * vx(ijp1  ,k,l) &
                    + cdiv(4,n,l,1) * vx(im1j  ,k,l) &
                    + cdiv(5,n,l,1) * vx(im1jm1,k,l) &
                    + cdiv(6,n,l,1) * vx(ijm1  ,k,l) &
                    + cdiv(0,n,l,2) * vy(ij    ,k,l) &
                    + cdiv(1,n,l,2) * vy(ip1j  ,k,l) &
                    + cdiv(2,n,l,2) * vy(ip1jp1,k,l) &
                    + cdiv(3,n,l,2) * vy(ijp1  ,k,l) &
                    + cdiv(4,n,l,2) * vy(im1j  ,k,l) &
                    + cdiv(5,n,l,2) * vy(im1jm1,k,l) &
                    + cdiv(6,n,l,2) * vy(ijm1  ,k,l) &
                    + cdiv(0,n,l,3) * vz(ij    ,k,l) &
                    + cdiv(1,n,l,3) * vz(ip1j  ,k,l) &
                    + cdiv(2,n,l,3) * vz(ip1jp1,k,l) &
                    + cdiv(3,n,l,3) * vz(ijp1  ,k,l) &
                    + cdiv(4,n,l,3) * vz(im1j  ,k,l) &
                    + cdiv(5,n,l,3) * vz(im1jm1,k,l) &
                    + cdiv(6,n,l,3) * vz(ijm1  ,k,l) ) * mfact
    enddo
    enddo
    enddo
    !$acc end kernels
    !$acc end data

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          scl_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( cdiv_pl(v-1,n,l,1) * vx_pl(v,k,l) &
                                             + cdiv_pl(v-1,n,l,2) * vy_pl(v,k,l) &
                                             + cdiv_pl(v-1,n,l,3) * vz_pl(v,k,l) )
          enddo

          scl_pl(n,k,l) = scl_pl(n,k,l) * mfact
       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_divergence')

    return
  end subroutine OPRT_divergence

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient( &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl,  &
       scl, scl_pl, &
       mfact        )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(8), intent(in)  :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: mfact

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_gradient')

    !!$acc data present_or_copyin(cgrad,scl) copyout(vx,vy,vz)
    !!$acc kernels
    !!$acc loop seq
    do l = 1, ADM_lall
    !!$acc loop gang vector(32)
    do k = ADM_kmin, ADM_kmax
    !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d
       im1j   = n - 1
       ijm1   = n     - ADM_gall_1d
       im1jm1 = n - 1 - ADM_gall_1d

       vx(n,k,l) = ( cgrad(0,ij,l,1) * scl(ij    ,k,l) &
                   + cgrad(1,ij,l,1) * scl(ip1j  ,k,l) &
                   + cgrad(2,ij,l,1) * scl(ip1jp1,k,l) &
                   + cgrad(3,ij,l,1) * scl(ijp1  ,k,l) &
                   + cgrad(4,ij,l,1) * scl(im1j  ,k,l) &
                   + cgrad(5,ij,l,1) * scl(im1jm1,k,l) &
                   + cgrad(6,ij,l,1) * scl(ijm1  ,k,l) ) * mfact

       vy(n,k,l) = ( cgrad(0,ij,l,2) * scl(ij    ,k,l) &
                   + cgrad(1,ij,l,2) * scl(ip1j  ,k,l) &
                   + cgrad(2,ij,l,2) * scl(ip1jp1,k,l) &
                   + cgrad(3,ij,l,2) * scl(ijp1  ,k,l) &
                   + cgrad(4,ij,l,2) * scl(im1j  ,k,l) &
                   + cgrad(5,ij,l,2) * scl(im1jm1,k,l) &
                   + cgrad(6,ij,l,2) * scl(ijm1  ,k,l) ) * mfact

       vz(n,k,l) = ( cgrad(0,ij,l,3) * scl(ij    ,k,l) &
                   + cgrad(1,ij,l,3) * scl(ip1j  ,k,l) &
                   + cgrad(2,ij,l,3) * scl(ip1jp1,k,l) &
                   + cgrad(3,ij,l,3) * scl(ijp1  ,k,l) &
                   + cgrad(4,ij,l,3) * scl(im1j  ,k,l) &
                   + cgrad(5,ij,l,3) * scl(im1jm1,k,l) &
                   + cgrad(6,ij,l,3) * scl(ijm1  ,k,l) ) * mfact
    enddo
    enddo
    enddo
    !!$acc end kernels
    !!$acc end data

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          vx_pl(n,k,l) = 0.D0
          vy_pl(n,k,l) = 0.D0
          vz_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             vx_pl(n,k,l) = vx_pl(n,k,l) + cgrad_pl(v-1,n,l,1) * scl_pl(v,k,l)
             vy_pl(n,k,l) = vy_pl(n,k,l) + cgrad_pl(v-1,n,l,2) * scl_pl(v,k,l)
             vz_pl(n,k,l) = vz_pl(n,k,l) + cgrad_pl(v-1,n,l,3) * scl_pl(v,k,l)
          enddo

          vx_pl(n,k,l) = vx_pl(n,k,l) * mfact
          vy_pl(n,k,l) = vy_pl(n,k,l) * mfact
          vz_pl(n,k,l) = vz_pl(n,k,l) * mfact
       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_gradient')

    return
  end subroutine OPRT_gradient

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       mfact          )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(8), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: mfact

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_laplacian')

    !!$acc data present_or_copyin(clap,scl) copyout(dscl)
    !!$acc kernels
    !!$acc loop seq
    do l = 1, ADM_lall
    !!$acc loop gang vector(32)
    do k = ADM_kmin, ADM_kmax
    !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d
       im1j   = n - 1
       ijm1   = n     - ADM_gall_1d
       im1jm1 = n - 1 - ADM_gall_1d

       dscl(n,k,l) = ( clap(0,ij,l) * scl(ij    ,k,l) &
                     + clap(1,ij,l) * scl(ip1j  ,k,l) &
                     + clap(2,ij,l) * scl(ip1jp1,k,l) &
                     + clap(3,ij,l) * scl(ijp1  ,k,l) &
                     + clap(4,ij,l) * scl(im1j  ,k,l) &
                     + clap(5,ij,l) * scl(im1jm1,k,l) &
                     + clap(6,ij,l) * scl(ijm1  ,k,l) ) * mfact
    enddo
    enddo
    enddo
    !!$acc end kernels
    !!$acc end data

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          dscl_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + clap_pl(v-1,n,l) * scl_pl(v,k,l)
          enddo

          dscl_pl(n,k,l) = dscl_pl(n,k,l) * mfact
       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_laplacian')

    return
  end subroutine OPRT_laplacian

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       kh,   kh_pl,   &
       mfact          )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    implicit none

    real(8), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: mfact

    real(8)  :: vxt    (ADM_gall,TI:TJ)
    real(8)  :: vyt    (ADM_gall,TI:TJ)
    real(8)  :: vzt    (ADM_gall,TI:TJ)
    real(8)  :: flux   (ADM_gall,AI:AJ)
    real(8)  :: vxt_pl (ADM_gall_pl)
    real(8)  :: vyt_pl (ADM_gall_pl)
    real(8)  :: vzt_pl (ADM_gall_pl)
    real(8)  :: flux_pl(ADM_gall_pl)

    real(8) :: u1, u2, u3, smean

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_diffusion')

    !!$acc data present_or_copyin(cinterp_TN,cinterp_TRA,cinterp_HN,cinterp_PRA,scl,kh) copyout(dscl) &
    !!$acc&     create(vxt,vyt,vzt,flux)
    !!$acc kernels
    !!$acc loop seq
    do l = 1, ADM_lall
    !!$acc loop gang vector(32)
    do k = 1, ADM_kall

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          smean = ( scl(ij,k,l) + scl(ip1j,k,l) + scl(ip1jp1,k,l) ) / 3.D0

          u1 = 0.5D0 * (scl(ij    ,k,l)+scl(ip1j  ,k,l)) - smean
          u2 = 0.5D0 * (scl(ip1j  ,k,l)+scl(ip1jp1,k,l)) - smean
          u3 = 0.5D0 * (scl(ip1jp1,k,l)+scl(ij    ,k,l)) - smean

          vxt(n,TI) = ( - u1 * cinterp_TN(AI ,1,ij  ,l) &
                        - u2 * cinterp_TN(AJ ,1,ip1j,l) &
                        + u3 * cinterp_TN(AIJ,1,ij  ,l) ) * cinterp_TRA(TI,ij,l)
          vyt(n,TI) = ( - u1 * cinterp_TN(AI ,2,ij  ,l) &
                        - u2 * cinterp_TN(AJ ,2,ip1j,l) &
                        + u3 * cinterp_TN(AIJ,2,ij  ,l) ) * cinterp_TRA(TI,ij,l)
          vzt(n,TI) = ( - u1 * cinterp_TN(AI ,3,ij  ,l) &
                        - u2 * cinterp_TN(AJ ,3,ip1j,l) &
                        + u3 * cinterp_TN(AIJ,3,ij  ,l) ) * cinterp_TRA(TI,ij,l)
       enddo

       !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          smean = ( scl(ij,k,l) + scl(ip1jp1,k,l) + scl(ijp1,k,l) ) / 3.D0

          u1 = 0.5D0 * (scl(ij    ,k,l)+scl(ip1jp1,k,l)) - smean
          u2 = 0.5D0 * (scl(ip1jp1,k,l)+scl(ijp1  ,k,l)) - smean
          u3 = 0.5D0 * (scl(ijp1  ,k,l)+scl(ij    ,k,l)) - smean

          vxt(n,TJ) = ( - u1 * cinterp_TN(AIJ,1,ij  ,l) &
                        + u2 * cinterp_TN(AI ,1,ijp1,l) &
                        + u3 * cinterp_TN(AJ ,1,ij  ,l) ) * cinterp_TRA(TJ,ij,l)
          vyt(n,TJ) = ( - u1 * cinterp_TN(AIJ,2,ij  ,l) &
                        + u2 * cinterp_TN(AI ,2,ijp1,l) &
                        + u3 * cinterp_TN(AJ ,2,ij  ,l) ) * cinterp_TRA(TJ,ij,l)
          vzt(n,TJ) = ( - u1 * cinterp_TN(AIJ,3,ij  ,l) &
                        + u2 * cinterp_TN(AI ,3,ijp1,l) &
                        + u3 * cinterp_TN(AJ ,3,ij  ,l) ) * cinterp_TRA(TJ,ij,l)
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          vxt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vxt(suf(ADM_gmin,ADM_gmin-1),TJ)
          vyt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vyt(suf(ADM_gmin,ADM_gmin-1),TJ)
          vzt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vzt(suf(ADM_gmin,ADM_gmin-1),TJ)
       endif

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          flux(n,AI ) = 0.25D0 * ( (vxt(ijm1,TJ)+vxt(ij  ,TI)) * cinterp_HN(AI ,1,ij,l) &
                                 + (vyt(ijm1,TJ)+vyt(ij  ,TI)) * cinterp_HN(AI ,2,ij,l) &
                                 + (vzt(ijm1,TJ)+vzt(ij  ,TI)) * cinterp_HN(AI ,3,ij,l) &
                                 ) * (kh(ij,k,l)+kh(ip1j  ,k,l))
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
       do n = nstart, nend
          ij     = n
          ip1jp1 = n + 1 + ADM_gall_1d

          flux(n,AIJ) = 0.25D0 * ( (vxt(ij  ,TI)+vxt(ij  ,TJ)) * cinterp_HN(AIJ,1,ij,l) &
                                 + (vyt(ij  ,TI)+vyt(ij  ,TJ)) * cinterp_HN(AIJ,2,ij,l) &
                                 + (vzt(ij  ,TI)+vzt(ij  ,TJ)) * cinterp_HN(AIJ,3,ij,l) &
                                 ) * (kh(ij,k,l)+kh(ip1jp1,k,l))
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1

          flux(n,AJ ) = 0.25D0 * ( (vxt(ij  ,TJ)+vxt(im1j,TI)) * cinterp_HN(AJ ,1,ij,l) &
                                 + (vyt(ij  ,TJ)+vyt(im1j,TI)) * cinterp_HN(AJ ,2,ij,l) &
                                 + (vzt(ij  ,TJ)+vzt(im1j,TI)) * cinterp_HN(AJ ,3,ij,l) &
                                 ) * (kh(ij,k,l)+kh(ijp1  ,k,l))
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          flux(suf(ADM_gmin,ADM_gmin-1),AJ) = 0.D0
       endif

       !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
       do n = OPRT_nstart, OPRT_nend
          ij     = n
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          dscl(n,k,l) = ( flux(ij,AI ) - flux(im1j  ,AI ) &
                        + flux(ij,AIJ) - flux(im1jm1,AIJ) &
                        + flux(ij,AJ ) - flux(ijm1  ,AJ ) ) * cinterp_PRA(ij,l) * mfact
       enddo

    enddo
    enddo
    !!$acc end kernels
    !!$acc end data

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             smean = ( scl_pl(n,k,l) + scl_pl(ij,k,l) + scl_pl(ijp1,k,l) ) / 3.D0

             u1 = 0.5D0 * (scl_pl(n   ,k,l)+scl_pl(ij  ,k,l)) - smean
             u2 = 0.5D0 * (scl_pl(ij  ,k,l)+scl_pl(ijp1,k,l)) - smean
             u3 = 0.5D0 * (scl_pl(ijp1,k,l)+scl_pl(n   ,k,l)) - smean

             vxt_pl(v) = ( + u1 * GMTR_A_var_pl(ij  ,k0,l,TNX ) &
                           + u2 * GMTR_A_var_pl(ij  ,k0,l,TN2X) &
                           - u3 * GMTR_A_var_pl(ijp1,k0,l,TNX ) ) * GMTR_T_var_pl(ij,k0,l,T_RAREA)
             vyt_pl(v) = ( + u1 * GMTR_A_var_pl(ij  ,k0,l,TNY ) &
                           + u2 * GMTR_A_var_pl(ij  ,k0,l,TN2Y) &
                           - u3 * GMTR_A_var_pl(ijp1,k0,l,TNY ) ) * GMTR_T_var_pl(ij,k0,l,T_RAREA)
             vzt_pl(v) = ( + u1 * GMTR_A_var_pl(ij  ,k0,l,TNZ ) &
                           + u2 * GMTR_A_var_pl(ij  ,k0,l,TN2Z) &
                           - u3 * GMTR_A_var_pl(ijp1,k0,l,TNZ ) ) * GMTR_T_var_pl(ij,k0,l,T_RAREA)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             flux_pl(v) = 0.25D0 * ( (vxt_pl(ij)+vxt_pl(ijm1)) * GMTR_A_var_pl(ij,k0,l,HNX) &
                                   + (vyt_pl(ij)+vyt_pl(ijm1)) * GMTR_A_var_pl(ij,k0,l,HNY) &
                                   + (vzt_pl(ij)+vzt_pl(ijm1)) * GMTR_A_var_pl(ij,k0,l,HNZ) &
                                   ) * (kh_pl(n,k,l)+kh_pl(ij,k,l))
          enddo

          dscl_pl(n,k,l) = 0.D0
          do v = ADM_gmin_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + flux_pl(v)
          enddo
          dscl_pl(n,k,l) = dscl_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_RAREA) * mfact

       enddo
       enddo

    endif

    call DEBUG_rapend('++++OPRT_diffusion')

    return
  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------
  subroutine OPRT_horizontalize_vec( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_e,    &
       GRD_e_pl
    implicit none

    real(8), intent(inout) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: prd
    integer :: n, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_horizontalize_vec')

    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax
    do n = 1, ADM_gall
       prd = ( vx(n,k,l) * GRD_e(n,l,GRD_XDIR) &
             + vy(n,k,l) * GRD_e(n,l,GRD_YDIR) &
             + vz(n,k,l) * GRD_e(n,l,GRD_ZDIR) )

       vx(n,k,l) = vx(n,k,l) - prd * GRD_e(n,l,GRD_XDIR)
       vy(n,k,l) = vy(n,k,l) - prd * GRD_e(n,l,GRD_YDIR)
       vz(n,k,l) = vz(n,k,l) - prd * GRD_e(n,l,GRD_ZDIR)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
       do n = 1, ADM_gall_pl
          prd = ( vx_pl(n,k,l) * GRD_e_pl(n,l,GRD_XDIR) &
                + vy_pl(n,k,l) * GRD_e_pl(n,l,GRD_YDIR) &
                + vz_pl(n,k,l) * GRD_e_pl(n,l,GRD_ZDIR) )

          vx_pl(n,k,l) = vx_pl(n,k,l) - prd * GRD_e_pl(n,l,GRD_XDIR)
          vy_pl(n,k,l) = vy_pl(n,k,l) - prd * GRD_e_pl(n,l,GRD_YDIR)
          vz_pl(n,k,l) = vz_pl(n,k,l) - prd * GRD_e_pl(n,l,GRD_ZDIR)
       enddo
       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_horizontalize_vec')

    return
  end subroutine OPRT_horizontalize_vec

  !-----------------------------------------------------------------------------
  subroutine OPRT_vorticity( &
       scl, scl_pl, &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl,  &
       mfact        )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    implicit none

    real(8), intent(out) :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: mfact

    real(8)  :: vxt    (ADM_gall   ,TI:TJ)
    real(8)  :: vxt_pl (ADM_gall_pl)
    real(8)  :: vyt    (ADM_gall   ,TI:TJ)
    real(8)  :: vyt_pl (ADM_gall_pl)
    real(8)  :: vzt    (ADM_gall   ,TI:TJ)
    real(8)  :: vzt_pl (ADM_gall_pl)
    real(8)  :: flux   (ADM_gall   ,AI:AJ)
    real(8)  :: flux_pl(ADM_gall_pl)

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_vorticity')

    do l = 1, ADM_lall
    do k = 1, ADM_kall

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          vxt(n,TI) = vx(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                    + vx(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                    + vx(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
          vyt(n,TI) = vy(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                    + vy(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                    + vy(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
          vzt(n,TI) = vz(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                    + vz(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                    + vz(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
       enddo

       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          vxt(n,TJ) = vx(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                    + vx(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                    + vx(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
          vyt(n,TJ) = vy(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                    + vy(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                    + vy(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
          vzt(n,TJ) = vz(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                    + vz(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                    + vz(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
       enddo

       if ( ADM_have_sgp(l) ) then
          vxt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vxt(suf(ADM_gmin,ADM_gmin-1),TJ)
          vyt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vyt(suf(ADM_gmin,ADM_gmin-1),TJ)
          vzt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vzt(suf(ADM_gmin,ADM_gmin-1),TJ)
       endif

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d

          flux(n,AI) = 0.5D0 * ( (vxt(ijm1,TJ)+vxt(ij,TI)) * cinterp_HN(AI ,1,ij,l) &
                               + (vyt(ijm1,TJ)+vyt(ij,TI)) * cinterp_HN(AI ,2,ij,l) &
                               + (vzt(ijm1,TJ)+vzt(ij,TI)) * cinterp_HN(AI ,3,ij,l) )
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          flux(n,AIJ) = 0.5D0 * ( (vxt(ij,TI)+vxt(ij,TJ)) * cinterp_HN(AIJ,1,ij,l) &
                                + (vyt(ij,TI)+vyt(ij,TJ)) * cinterp_HN(AIJ,2,ij,l) &
                                + (vzt(ij,TI)+vzt(ij,TJ)) * cinterp_HN(AIJ,3,ij,l) )
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart,nend
          ij     = n
          im1j   = n - 1

          flux(n,AJ) = 0.5D0 * ( ( vxt(ij,TJ)+vxt(im1j,TI)) * cinterp_HN(AJ ,1,ij,l) &
                               + ( vyt(ij,TJ)+vyt(im1j,TI)) * cinterp_HN(AJ ,2,ij,l) &
                               + ( vzt(ij,TJ)+vzt(im1j,TI)) * cinterp_HN(AJ ,3,ij,l) )
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          flux(suf(ADM_gmin,ADM_gmin-1),AJ) = 0.D0
       endif

       do n = OPRT_nstart, OPRT_nend
          ij     = n
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          scl(n,k,l) = - ( flux(ij,AI ) - flux(im1j  ,AI ) &
                         + flux(ij,AIJ) - flux(im1jm1,AIJ) &
                         + flux(ij,AJ ) - flux(ijm1  ,AJ ) ) * cinterp_PRA(ij,l) * mfact
       enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             vxt_pl(v) = vx_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                       + vx_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                       + vx_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)

             vyt_pl(v) = vy_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                       + vy_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                       + vy_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)

             vzt_pl(v) = vz_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                       + vz_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                       + vz_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             flux_pl(v) = 0.5D0 * ( (vxt_pl(ijm1)+vxt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HTX) &
                                  + (vyt_pl(ijm1)+vyt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HTY) &
                                  + (vzt_pl(ijm1)+vzt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HTZ) )
          enddo

          scl_pl(n,k,l) = 0.D0
          do v = ADM_gmin_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + flux_pl(v)
          enddo
          scl_pl(n,k,l) = scl_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_RAREA) * mfact

       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_vorticity')

    return
  end subroutine OPRT_vorticity

  !-----------------------------------------------------------------------------
  subroutine OPRT_divdamp( &
       grdx, grdx_pl, &
       grdy, grdy_pl, &
       grdz, grdz_pl, &
       vx,   vx_pl,   &
       vy,   vy_pl,   &
       vz,   vz_pl,   &
       mfact          )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl,  &
       ADM_kmin,     &
       ADM_kmax
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    implicit none

    real(8), intent(out) :: grdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vx     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: mfact

    real(8) :: sclt   (ADM_gall   ,TI:TJ)
    real(8) :: sclt_pl(ADM_gall_pl,ADM_kall)

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_divdamp')

    !!$acc data present_or_copyin(cinterp_TN,cinterp_TRA,cinterp_HN,cinterp_PRA,vx,vy,vz) &
    !!$acc&     copyout(grdx,grdy,grdz) create(sclt)
    !!$acc kernels
    !!$acc loop seq
    do l = 1, ADM_lall
       !!$acc loop gang vector(32)
       do k = ADM_kmin, ADM_kmax

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )

          !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
          do n = nstart, nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d

             sclt(n,TI) = ( - ( vx(ij    ,k,l) + vx(ip1j  ,k,l) ) * cinterp_TN(AI ,1,ij  ,l) &
                            - ( vx(ip1j  ,k,l) + vx(ip1jp1,k,l) ) * cinterp_TN(AJ ,1,ip1j,l) &
                            + ( vx(ip1jp1,k,l) + vx(ij    ,k,l) ) * cinterp_TN(AIJ,1,ij  ,l) &
                            - ( vy(ij    ,k,l) + vy(ip1j  ,k,l) ) * cinterp_TN(AI ,2,ij  ,l) &
                            - ( vy(ip1j  ,k,l) + vy(ip1jp1,k,l) ) * cinterp_TN(AJ ,2,ip1j,l) &
                            + ( vy(ip1jp1,k,l) + vy(ij    ,k,l) ) * cinterp_TN(AIJ,2,ij  ,l) &
                            - ( vz(ij    ,k,l) + vz(ip1j  ,k,l) ) * cinterp_TN(AI ,3,ij  ,l) &
                            - ( vz(ip1j  ,k,l) + vz(ip1jp1,k,l) ) * cinterp_TN(AJ ,3,ip1j,l) &
                            + ( vz(ip1jp1,k,l) + vz(ij    ,k,l) ) * cinterp_TN(AIJ,3,ij  ,l) &
                          ) * 0.5D0 * cinterp_TRA(TI,ij,l)

             sclt(n,TJ) = ( - ( vx(ij    ,k,l) + vx(ip1jp1,k,l) ) * cinterp_TN(AIJ,1,ij  ,l) &
                            + ( vx(ip1jp1,k,l) + vx(ijp1  ,k,l) ) * cinterp_TN(AI ,1,ijp1,l) &
                            + ( vx(ijp1  ,k,l) + vx(ij    ,k,l) ) * cinterp_TN(AJ ,1,ij  ,l) &
                            - ( vy(ij    ,k,l) + vy(ip1jp1,k,l) ) * cinterp_TN(AIJ,2,ij  ,l) &
                            + ( vy(ip1jp1,k,l) + vy(ijp1  ,k,l) ) * cinterp_TN(AI ,2,ijp1,l) &
                            + ( vy(ijp1  ,k,l) + vy(ij    ,k,l) ) * cinterp_TN(AJ ,2,ij  ,l) &
                            - ( vz(ij    ,k,l) + vz(ip1jp1,k,l) ) * cinterp_TN(AIJ,3,ij  ,l) &
                            + ( vz(ip1jp1,k,l) + vz(ijp1  ,k,l) ) * cinterp_TN(AI ,3,ijp1,l) &
                            + ( vz(ijp1  ,k,l) + vz(ij    ,k,l) ) * cinterp_TN(AJ ,3,ij  ,l) &
                          ) * 0.5D0 * cinterp_TRA(TJ,ij,l)
          enddo

          !!$acc loop gang vector(32) private(ij,ip1j,ijp1,ip1jp1,im1j,ijm1,im1jm1)
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             grdx(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,1,ij,    l) &
                             + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,1,ij,    l) &
                             + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,1,ij,    l) &
                             - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,1,im1j,  l) &
                             - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,1,im1jm1,l) &
                             - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(AJ ,1,ijm1,  l) &
                           ) * 0.5D0 * cinterp_PRA(ij,l) * mfact

             grdy(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,2,ij,    l) &
                             + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,2,ij,    l) &
                             + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,2,ij,    l) &
                             - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,2,im1j,  l) &
                             - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,2,im1jm1,l) &
                             - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(AJ ,2,ijm1,  l) &
                           ) * 0.5D0 * cinterp_PRA(ij,l) * mfact

             grdz(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,3,ij,    l) &
                             + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,3,ij,    l) &
                             + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,3,ij,    l) &
                             - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,3,im1j,  l) &
                             - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,3,im1jm1,l) &
                             - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(AJ ,3,ijm1,  l) &
                           ) * 0.5D0 * cinterp_PRA(ij,l) * mfact
          enddo

          if ( ADM_have_sgp(l) ) then
             n = suf(ADM_gmin,ADM_gmin)

             ij     = n
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             sclt(im1jm1,TI) = sclt(ijm1,TJ) ! copy

             grdx(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,1,ij,    l) &
                             + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,1,ij,    l) &
                             + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,1,ij,    l) &
                             - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,1,im1j,  l) &
                             - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,1,im1jm1,l) &
                           ) * 0.5D0 * cinterp_PRA(ij,l) * mfact

             grdy(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,2,ij,    l) &
                             + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,2,ij,    l) &
                             + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,2,ij,    l) &
                             - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,2,im1j,  l) &
                             - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,2,im1jm1,l) &
                           ) * 0.5D0 * cinterp_PRA(ij,l) * mfact

             grdz(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,3,ij,    l) &
                             + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,3,ij,    l) &
                             + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,3,ij,    l) &
                             - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,3,im1j,  l) &
                             - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,3,im1jm1,l) &
                           ) * 0.5D0 * cinterp_PRA(ij,l) * mfact
          endif
       enddo

       grdx   (:,ADM_kmin-1,l) = 0.D0
       grdx   (:,ADM_kmax+1,l) = 0.D0
       grdy   (:,ADM_kmin-1,l) = 0.D0
       grdy   (:,ADM_kmax+1,l) = 0.D0
       grdz   (:,ADM_kmin-1,l) = 0.D0
       grdz   (:,ADM_kmax+1,l) = 0.D0
    enddo
    !!$acc end kernels
    !!$acc end data

    if ( ADM_have_pl ) then
       n = ADM_GSLF_PL

       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             sclt_pl(v,k) = ( + ( vx_pl(n   ,k,l) + vx_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNX ) &
                              + ( vy_pl(n   ,k,l) + vy_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNY ) &
                              + ( vz_pl(n   ,k,l) + vz_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNZ ) &
                              + ( vx_pl(ij  ,k,l) + vx_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2X) &
                              + ( vy_pl(ij  ,k,l) + vy_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2Y) &
                              + ( vz_pl(ij  ,k,l) + vz_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2Z) &
                              - ( vx_pl(ijp1,k,l) + vx_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNX ) &
                              - ( vy_pl(ijp1,k,l) + vy_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNY ) &
                              - ( vz_pl(ijp1,k,l) + vz_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNZ ) &
                            ) * 0.5D0 * GMTR_T_var_pl(ij,k0,l,T_RAREA)
          enddo

          grdx_pl(n,k,l) = 0.D0
          grdy_pl(n,k,l) = 0.D0
          grdz_pl(n,k,l) = 0.D0

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 < ADM_gmin_pl ) ijm1 = ADM_gmax_pl ! cyclic condition

             grdx_pl(n,k,l) = grdx_pl(n,k,l) + ( sclt_pl(ijm1,k) + sclt_pl(ij,k) ) * GMTR_A_var_pl(ij,k0,l,HNX)
             grdy_pl(n,k,l) = grdy_pl(n,k,l) + ( sclt_pl(ijm1,k) + sclt_pl(ij,k) ) * GMTR_A_var_pl(ij,k0,l,HNY)
             grdz_pl(n,k,l) = grdz_pl(n,k,l) + ( sclt_pl(ijm1,k) + sclt_pl(ij,k) ) * GMTR_A_var_pl(ij,k0,l,HNZ)
          enddo

          grdx_pl(n,k,l) = grdx_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA) * mfact
          grdy_pl(n,k,l) = grdy_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA) * mfact
          grdz_pl(n,k,l) = grdz_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,P_RAREA) * mfact

       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_divdamp')

    return
  end subroutine OPRT_divdamp

end module mod_oprt
!-------------------------------------------------------------------------------
