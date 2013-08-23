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
     ADM_LOG_FID
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
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, save :: OPRT_nstart
  integer, private, save :: OPRT_nend

  ! < for divergence operator >
  real(8), private, allocatable, save :: cdiv   (:,:,:,:)
  real(8), private, allocatable, save :: cdiv_pl(:,:,:,:)

  ! < for gradient operator >
  real(8), private, allocatable, save :: cgrad   (:,:,:,:)
  real(8), private, allocatable, save :: cgrad_pl(:,:,:,:)

  ! < for laplacian operator >
  real(8), private, allocatable, save :: clap   (:,:,:)
  real(8), private, allocatable, save :: clap_pl(:,:,:)

  ! < for diffusion operator >
  real(8), private, allocatable, save :: cmdif_P (:,:)     !(n,l) <- GMTR_P_VAR(n,1,l,P_RAREA)
  real(8), private, allocatable, save :: cmdif_T (:,:,:)   !(TI:TJ,n,l) <- GMTR_T_VAR(n,1,l,TI:TJ,T_RAREA)
  real(8), private, allocatable, save :: cmdif_AH(:,:,:,:) !(AI:AJ,1:3,n,l) <-GMTR_A_VAR(n,1,l,TI:TJ,HN[XYZ])
  real(8), private, allocatable, save :: cmdif_AT(:,:,:,:) !(AI:AJ,1:3,n,l) <-GMTR_A_VAR(n,1,l,TI:TJ,TN[XYZ])

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OPRT_setup
    use mod_adm, only: &
       ADM_W,          &
       ADM_TI,         &
       ADM_TJ,         &
       ADM_AI,         &
       ADM_AIJ,        &
       ADM_AJ,         &
       ADM_prc_tab,    &
       ADM_prc_me,     &
       ADM_prc_pl,     &
       ADM_rgn_vnum,   &
       ADM_vlink_nmax, &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,    &
       ADM_KNONE,      &
       ADM_kmin,       &
       ADM_kmax
    use mod_gmtr, only: &
       GMTR_P_rarea,  &
       GMTR_T_W1,     &
       GMTR_T_W2,     &
       GMTR_T_W3,     &
       GMTR_T_rarea,  &
       GMTR_A_hnx,    &
       GMTR_A_hny,    &
       GMTR_A_hnz,    &
       GMTR_A_tnx,    &
       GMTR_A_tny,    &
       GMTR_A_tnz,    &
       GMTR_A_tn2x,   &
       GMTR_A_tn2y,   &
       GMTR_A_tn2z,   &
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    implicit none

    integer :: n0,n1,n2,n3,n4

    integer :: ij
    integer :: im1j,ijm1,im1jm1
    integer :: ip1j,ijp1,ip1jp1

    integer :: rgnid
    integer :: n, l, m, md, v

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: TI,TJ,AI,AIJ,AJ,W1,W2,W3
    integer :: k0,a0
    integer :: tx1,ty1,tz1
    integer :: tx2,ty2,tz2
    integer :: hx1,hy1,hz1
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[oprt]/Category[common share]'

    k0  = ADM_KNONE

    TI  = ADM_TI
    TJ  = ADM_TJ
    AI  = ADM_AI
    AIJ = ADM_AIJ
    AJ  = ADM_AJ
    W1  = GMTR_T_W1
    W2  = GMTR_T_W2
    W3  = GMTR_T_W3

    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)

    !---< setup coefficient of divergence operator >
    write(ADM_LOG_FID,*) '*** setup coefficient of divergence operator'

    allocate( cdiv   (0:6,             ADM_gall   ,ADM_lall   ,3) )
    allocate( cdiv_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl,3) )

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do m = 1, 3
          md = m + GMTR_A_HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

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
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cdiv(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cdiv(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cdiv(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0*GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cdiv(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cdiv(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cdiv(6,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
          enddo
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

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
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ip1j
             cdiv(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cdiv(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ijp1
             cdiv(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! im1j
             cdiv(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! im1jm1
             cdiv(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ijm1
             cdiv(6,n,l,m) = ( + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
          enddo
       endif

    enddo ! loop l

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             cdiv_pl(0,n,l,m) = 0.D0
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cdiv_pl(0,n,l,m) = cdiv_pl(0,n,l,m) + ( GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                                      + GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ijp1,k0,l,md) )
             enddo
             cdiv_pl(0,n,l,m) = cdiv_pl(0,n,l,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)

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
                                     ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif



    !---< setup coefficient of gradient operator >

    write(ADM_LOG_FID,*) '*** setup coefficient of gradient operator'

    allocate( cgrad   (0:6,             ADM_gall   ,ADM_lall   ,3) )
    allocate( cgrad_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl,3) )

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do m = 1, 3
          md = m + GMTR_A_HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

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
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cgrad(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cgrad(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cgrad(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cgrad(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cgrad(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cgrad(6,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
          enddo
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

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
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cgrad(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cgrad(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cgrad(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cgrad(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cgrad(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cgrad(6,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
          enddo
       endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             cgrad_pl(0,n,l,m) = 0.D0
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cgrad_pl(0,n,l,m) = cgrad_pl(0,n,l,m) &
                                  + 2.D0 * ( GMTR_T_var_pl(ij,k0,l,W1) - 1.D0 ) * GMTR_A_var_pl(ijp1,k0,l,md)
             enddo
             cgrad_pl(0,n,l,m) = cdiv_pl(0,n,l,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)

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
                                      ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif

    ! ---- setup coefficient of laplacian operator

    write(ADM_LOG_FID,*) '*** setup coefficient of laplacian operator'

    allocate( clap   (0:6,             ADM_gall   ,ADM_lall   ) )
    allocate( clap_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl) )

    a0  = GMTR_T_rarea
    tx1 = GMTR_A_tnx
    ty1 = GMTR_A_tny
    tz1 = GMTR_A_tnz
    hx1 = GMTR_A_hnx
    hy1 = GMTR_A_hny
    hz1 = GMTR_A_hnz

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = OPRT_nstart, OPRT_nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          ! ij
          clap(0,ij,l) = ( &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          clap(0,ij,l) = clap(0,ij,l) + ( &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1j
          clap(1,ij,l) = ( &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1jp1
          clap(2,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijp1
          clap(3,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &     
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1j
          clap(4,ij,l) = ( &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1jm1
          clap(5,ij,l) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &  
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijm1
          clap(6,ij,l) = ( &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +2*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +2*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +2*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          ! 0: ij
          clap(0,ij,l) = ( &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0  ! Y.Niwa add 06/08/22

          clap(0,ij,l) = clap(0,ij,l) + ( &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1j
          clap(1,ij,l) = ( &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1jp1
          clap(2,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijp1
          clap(3,ij,l) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &     
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1j
          clap(4,ij,l) = ( &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1jm1
          clap(5,ij,l) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijm1
          clap(6,ij,l) = ( &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &  
      +1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
      -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
      -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
      -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
         )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

      endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

      n =ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4
      k0=ADM_KNONE
      a0=GMTR_T_rarea
      tx1=gMtr_a_tnx
      tx2=GMTR_A_tn2x
      ty1=GMTR_A_tny
      ty2=GMTR_A_tn2y
      tz1=GMTR_A_tnz
      tz2=GMTR_A_tn2z
      hx1=GMTR_A_hnx
      hy1=GMTR_A_hny
      hz1=GMTR_A_hnz

      do l = 1,ADM_lall_pl
          clap_pl(0,n,l)=( &  
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0   ! Y.Niwa add 060822

          clap_pl(0,n,l)= clap_pl(0,n,l) + ( &                        ! Y.Niwa add 060822 
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n0
          clap_pl(1,n,l)=( &  
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n1
          clap_pl(2,n,l)=( &  
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n2
          clap_pl(3,n,l)=( &  
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n3
          clap_pl(4,n,l)=( &  
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n4
          clap_pl(5,n,l)=( &  
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
      enddo
    endif



    allocate( cmdif_P (                  ADM_gall,ADM_lall) )
    allocate( cmdif_T (ADM_TI:ADM_TJ,    ADM_gall,ADM_lall) )
    allocate( cmdif_AT(ADM_AI:ADM_AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cmdif_AH(ADM_AI:ADM_AJ,1:3,ADM_gall,ADM_lall) )

    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cmdif_P(n,l) = GMTR_P_var(n,k0,l,GMTR_P_RAREA)
    enddo
    enddo

    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cmdif_T(TI,n,l) = GMTR_T_var(n,k0,l,TI,GMTR_T_RAREA)
       cmdif_T(TJ,n,l) = GMTR_T_var(n,k0,l,TJ,GMTR_T_RAREA)
    enddo
    enddo

    do m=ADM_AI,ADM_AJ
       do l=1,ADM_lall
          do n=1,ADM_gall
             cmdif_AT(m, 1 , n, l) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_TNX)
             cmdif_AT(m, 2 , n, l) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_TNY)
             cmdif_AT(m, 3 , n, l) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_TNZ)

             cmdif_AH(m, 1 , n, l) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_HNX)
             cmdif_AH(m, 2 , n, l) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_HNY)
             cmdif_AH(m, 3 , n, l) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_HNZ)
          enddo
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
       ADM_prc_me,  &
       ADM_prc_pl,  &
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

    real(8), intent(inout) :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: fact

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_divergence')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ip1jp1 = n + 1 + ADM_gall_1d
       ijp1   = n     + ADM_gall_1d
       im1j   = n - 1
       im1jm1 = n - 1 - ADM_gall_1d
       ijm1   = n     - ADM_gall_1d

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
                    + cdiv(6,n,l,3) * vz(ijm1  ,k,l) ) * fact
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          scl_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( cdiv_pl(v-1,n,l,1) * vx_pl(v,k,l) &
                                             + cdiv_pl(v-1,n,l,2) * vy_pl(v,k,l) &
                                             + cdiv_pl(v-1,n,l,3) * vz_pl(v,k,l) )
          enddo

          scl_pl(n,k,l) = scl_pl(n,k,l) * fact
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
       ADM_prc_me,  &
       ADM_prc_pl,  &
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

    real(8), intent(in)    :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: fact

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_gradient')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ip1jp1 = n + 1 + ADM_gall_1d
       ijp1   = n     + ADM_gall_1d
       im1j   = n - 1
       im1jm1 = n - 1 - ADM_gall_1d
       ijm1   = n     - ADM_gall_1d

       vx(n,k,l) = ( cgrad(0,ij,l,1) * scl(ij    ,k,l) &
                   + cgrad(1,ij,l,1) * scl(ip1j  ,k,l) &
                   + cgrad(2,ij,l,1) * scl(ip1jp1,k,l) &
                   + cgrad(3,ij,l,1) * scl(ijp1  ,k,l) &
                   + cgrad(4,ij,l,1) * scl(im1j  ,k,l) &
                   + cgrad(5,ij,l,1) * scl(im1jm1,k,l) &
                   + cgrad(6,ij,l,1) * scl(ijm1  ,k,l) ) * fact

       vy(n,k,l) = ( cgrad(0,ij,l,2) * scl(ij    ,k,l) &
                   + cgrad(1,ij,l,2) * scl(ip1j  ,k,l) &
                   + cgrad(2,ij,l,2) * scl(ip1jp1,k,l) &
                   + cgrad(3,ij,l,2) * scl(ijp1  ,k,l) &
                   + cgrad(4,ij,l,2) * scl(im1j  ,k,l) &
                   + cgrad(5,ij,l,2) * scl(im1jm1,k,l) &
                   + cgrad(6,ij,l,2) * scl(ijm1  ,k,l) ) * fact

       vz(n,k,l) = ( cgrad(0,ij,l,3) * scl(ij    ,k,l) &
                   + cgrad(1,ij,l,3) * scl(ip1j  ,k,l) &
                   + cgrad(2,ij,l,3) * scl(ip1jp1,k,l) &
                   + cgrad(3,ij,l,3) * scl(ijp1  ,k,l) &
                   + cgrad(4,ij,l,3) * scl(im1j  ,k,l) &
                   + cgrad(5,ij,l,3) * scl(im1jm1,k,l) &
                   + cgrad(6,ij,l,3) * scl(ijm1  ,k,l) ) * fact
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          vx_pl(n,k,l) = 0.D0
          vy_pl(n,k,l) = 0.D0
          vz_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             vx_pl(n,k,l) = vx_pl(n,k,l) + cgrad_pl(v-1,n,l,1) * scl_pl(v,k,l)
             vy_pl(n,k,l) = vy_pl(n,k,l) + cgrad_pl(v-1,n,l,2) * scl_pl(v,k,l)
             vz_pl(n,k,l) = vz_pl(n,k,l) + cgrad_pl(v-1,n,l,3) * scl_pl(v,k,l)
          enddo

          vx_pl(n,k,l) = vx_pl(n,k,l) * fact
          vy_pl(n,k,l) = vy_pl(n,k,l) * fact
          vz_pl(n,k,l) = vz_pl(n,k,l) * fact
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
       ADM_prc_me,  &
       ADM_prc_pl,  &
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

    real(8), intent(inout) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: fact

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_laplacian')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall
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
                     + clap(6,ij,l) * scl(ijm1  ,k,l) ) * fact
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          dscl_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + clap_pl(v-1,n,l) * scl_pl(v,k,l)
          enddo

          dscl_pl(n,k,l) = dscl_pl(n,k,l) * fact
       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_laplacian')

    return
  end subroutine OPRT_laplacian

  !-----------------------------------------------------------------------------  
  subroutine OPRT_diffusion( &
       dscl, dscl_pl,        &
       scl,scl_pl,           &
       kh, kh_pl,            &
       mfact )
    use mod_adm, only :   &
         !--- public parameters
         ADM_prc_pl,      &
         ADM_KNONE,       &
         ADM_TI,          &
         ADM_TJ,          &
         adm_ai,          &
         adm_aij,         &
         adm_aj,          &
         adm_w,           &
         ADM_lall_pl,     &
         ADM_GMAX_PL,     &
         ADM_gmin_pl,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         ADM_rgn_vnum,    &
         !--- public variables
         ADM_prc_me,      &
         ADM_prc_tab,     &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall
    use mod_gmtr, only : &
         !--- public parameters
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_T_rarea,    &
         GMTR_A_hnx,      &
         GMTR_A_hny,      &
         GMTR_A_hnz,      &
         GMTR_A_htx,      &
         GMTR_A_hty,      &
         GMTR_A_htz,      &
         GMTR_A_tnx,      &
         GMTR_A_tny,      &
         GMTR_A_tnz,      &
         GMTR_A_tn2x,     &
         GMTR_A_tn2y,     &
         GMTR_A_tn2z,     &
         GMTR_P_rarea,    &
         !--- public variables
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    implicit none

    real(8), intent(inout) :: dscl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: scl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: kh(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: kh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in), optional :: mfact
    !
    real(8)  :: vxt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vxt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: vyt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vyt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: vzt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vzt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: flux(        &
         ADM_gall,           &
         ADM_AI:ADM_AJ)
    real(8)  :: flux_pl(     &
         ADM_gall_pl)
    !
    integer :: l,n,k
    integer :: rgnid
    real(8) :: fact
    real(8) :: u1,u2,u3
    real(8) :: smean
    !
    integer :: nstart, nend
    integer :: nstart2, nstart3 !2011/02/25 for K by RIST 
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_diffusion')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact=1.0D0
    end if
    !
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k= 1, ADM_kall
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart, nend
             !
             smean = (scl(n,k,l)+scl(n+1,k,l)+scl(n+1+ADM_gall_1d,k,l))/3.0D0
             u1=-cmdif_T(ADM_TI, n, l) &
                    *(0.5D0*(scl(n,k,l)+scl(n+1,k,l))-smean)
             u2=-cmdif_T(ADM_TI, n, l) &
                    *(0.5D0*(scl(n+1,k,l)+scl(n+1+ADM_gall_1d,k,l))-smean)
             u3=+cmdif_T(ADM_TI, n, l) &
                    *(0.5D0*(scl(n+1+ADM_gall_1d,k,l)+scl(n,k,l))-smean)
             vxt( n, k, l, ADM_TI ) &
                  =u1*cmdif_AT(ADM_AI , 1 , n  , l)&
                  +u2*cmdif_AT(ADM_AJ , 1 , n+1, l)&
                  +u3*cmdif_AT(ADM_AIJ, 1 , n  , l)
             vyt( n, k, l, ADM_TI ) &
                  =u1*cmdif_AT(ADM_AI , 2 , n  , l)&
                  +u2*cmdif_AT(ADM_AJ , 2 , n+1, l)&
                  +u3*cmdif_AT(ADM_AIJ, 2 , n  , l)
             vzt( n, k, l, ADM_TI ) &
                  =u1*cmdif_AT(ADM_AI , 3 , n  , l)&
                  +u2*cmdif_AT(ADM_AJ , 3 , n+1, l)&
                  +u3*cmdif_AT(ADM_AIJ, 3 , n  , l)

             smean = (scl(n,k,l)+scl(n+1+ADM_gall_1d,k,l)+scl(n+ADM_gall_1d,k,l))/3.0D0
             u1=-cmdif_T(ADM_TJ, n, l) &
                    *(0.5D0*(scl(n              ,k,l)+scl(n+1+ADM_gall_1d,k,l))-smean)
             u2=+cmdif_T(ADM_TJ, n, l) &
                    *(0.5D0*(scl(n+1+ADM_gall_1d,k,l)+scl(n  +ADM_gall_1d,k,l))-smean)
             u3=+cmdif_T(ADM_TJ, n, l) &
                    *(0.5D0*(scl(n  +ADM_gall_1d,k,l)+scl(n              ,k,l))-smean)
             vxt( n, k, l, ADM_TJ ) &
                  =u1*cmdif_AT(ADM_AIJ, 1 , n            , l)&
                  +u2*cmdif_AT(ADM_AI , 1 , n+ADM_gall_1d, l)&
                  +u3*cmdif_AT(ADM_AJ , 1 , n            , l)
             vyt( n, k, l, ADM_TJ ) &
                  =u1*cmdif_AT(ADM_AIJ, 2 , n            , l)&
                  +u2*cmdif_AT(ADM_AI , 2 , n+ADM_gall_1d, l)&
                  +u3*cmdif_AT(ADM_AJ , 2 , n            , l)
             vzt( n, k, l, ADM_TJ ) &
                  =u1*cmdif_AT(ADM_AIJ, 3 , n            , l)&
                  +u2*cmdif_AT(ADM_AI , 3 , n+ADM_gall_1d, l)&
                  +u3*cmdif_AT(ADM_AJ , 3 , n            , l)
             !
          end do
          !
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             vxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
          !
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do k= 1, ADM_kall
          do l=1,ADM_lall_pl
             do n=ADM_gmin_pl,ADM_GMAX_PL-1
                smean = (scl_pl(ADM_gslf_pl,k,l)+scl_pl(n,k,l)+scl_pl(n+1,k,l))/3.0D0
                u1=+GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)&
                     *(0.5D0*(scl_pl(ADM_gslf_pl,k,l)+scl_pl(n,k,l))-smean)
                u2=+GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)&
                     *(0.5D0*(scl_pl(n,k,l)+scl_pl(n+1,k,l))-smean)
                u3=-GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)&
                     *(0.5D0*(scl_pl(n+1,k,l)+scl_pl(ADM_gslf_pl,k,l))-smean)
                vxt_pl(n,k,l)&
                     =u1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNX)&
                     +u2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2X)&
                     +u3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNX)
                vyt_pl(n,k,l)&
                     =u1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNY)&
                     +u2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2Y)&
                     +u3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNY)
                vzt_pl(n,k,l)&
                     =u1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNZ)&
                     +u2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2Z)&
                     +u3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNZ)
             enddo
             smean = (scl_pl(ADM_gslf_pl,k,l)+scl_pl(ADM_GMAX_PL,k,l)+scl_pl(ADM_gmin_pl,k,l))/3.0D0
             u1=+GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)&
                  *(0.5D0*(scl_pl(ADM_gslf_pl,k,l)+scl_pl(ADM_GMAX_PL,k,l))-smean)
             u2=+GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)&
                  *(0.5D0*(scl_pl(ADM_GMAX_PL,k,l)+scl_pl(ADM_gmin_pl,k,l))-smean)
             u3=-GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)&
                  *(0.5D0*(scl_pl(ADM_gmin_pl,k,l)+scl_pl(ADM_gslf_pl,k,l))-smean)
             vxt_pl(ADM_GMAX_PL,k,l)&
                  =u1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNX)&
                  +u2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2X)&
                  +u3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNX)
             vyt_pl(ADM_GMAX_PL,k,l)&
                  =u1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNY)&
                  +u2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2Y)&
                  +u3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNY)
             vzt_pl(ADM_GMAX_PL,k,l)&
                  =u1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNZ)&
                  +u2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2Z)&
                  +u3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNZ)
          enddo
       enddo
    end if
    !
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       !Do k=ADM_kmin,ADM_kmax
       do k= 1, ADM_kall
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux( n, ADM_AI )  =(&
                  &(vxt( n-ADM_gall_1d, k, l, ADM_TJ ) +vxt( n, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AI, 1 , n, l)&
                  +(vyt( n-ADM_gall_1d, k, l, ADM_TJ ) +vyt( n, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AI, 2 , n, l)&
                  +(vzt( n-ADM_gall_1d, k, l, ADM_TJ ) +vzt( n, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AI, 3 , n, l)&
                  &)*0.5D0
             flux( n, ADM_AI )  = flux( n, ADM_AI )  * (kh(n,k,l)+kh(n+1,k,l))*0.5D0
             flux( n, ADM_AIJ )  =(&
                  &(vxt( n, k, l, ADM_TI ) +vxt( n, k, l, ADM_TJ ) ) &
                  *cmdif_AH(ADM_AIJ, 1 , n, l)&
                  +(vyt( n, k, l, ADM_TI ) +vyt( n, k, l, ADM_TJ ) ) &
                  *cmdif_AH(ADM_AIJ, 2 , n, l)&
                  +(vzt( n, k, l, ADM_TI ) +vzt( n, k, l, ADM_TJ ) ) &
                  *cmdif_AH(ADM_AIJ, 3 , n, l)&
                  &)*0.5D0
             flux( n, ADM_AIJ )  = flux( n, ADM_AIJ )  * (kh(n,k,l)+kh(n+ADM_gall_1d+1,k,l))*0.5D0
             flux( n, ADM_AJ )  =(&
                  &(vxt( n, k, l, ADM_TJ ) +vxt( n-1, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AJ, 1 , n, l)&
                  +(vyt( n, k, l, ADM_TJ ) +vyt( n-1, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AJ, 2 , n, l)&
                  +(vzt( n, k, l, ADM_TJ ) +vzt( n-1, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AJ, 3 , n, l)&
                  &)*0.5D0
             flux( n, ADM_AJ )  = flux( n, ADM_AJ )  * (kh(n,k,l)+kh(n+ADM_gall_1d,k,l))*0.5D0
          end do
          nstart2 = suf(ADM_gmin-1,ADM_gmin-1)
          do n = nstart2, nstart-1
             flux( n, ADM_AIJ )  =(&
                  &(vxt( n, k, l, ADM_TI ) +vxt( n, k, l, ADM_TJ ) ) &
                  *cmdif_AH(ADM_AIJ, 1 , n, l)&
                  +(vyt( n, k, l, ADM_TI ) +vyt( n, k, l, ADM_TJ ) ) &
                  *cmdif_AH(ADM_AIJ, 2 , n, l)&
                  +(vzt( n, k, l, ADM_TI ) +vzt( n, k, l, ADM_TJ ) ) &
                  *cmdif_AH(ADM_AIJ, 3 , n, l)&
                  &)*0.5D0
             flux( n, ADM_AIJ )  = flux( n, ADM_AIJ )  * (kh(n,k,l)+kh(n+ADM_gall_1d+1,k,l))*0.5D0
          end do
          nstart3 = suf(ADM_gmin  ,ADM_gmin-1)
          do n = nstart3, nstart-1
             flux( n, ADM_AJ )  =(&
                  &(vxt( n, k, l, ADM_TJ ) +vxt( n-1, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AJ, 1 , n, l)&
                  +(vyt( n, k, l, ADM_TJ ) +vyt( n-1, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AJ, 2 , n, l)&
                  +(vzt( n, k, l, ADM_TJ ) +vzt( n-1, k, l, ADM_TI ) ) &
                  *cmdif_AH(ADM_AJ, 3 , n, l)&
                  &)*0.5D0
             flux( n, ADM_AJ )  = flux( n, ADM_AJ )  * (kh(n,k,l)+kh(n+ADM_gall_1d,k,l))*0.5D0
          end do
          !
          nstart = suf(ADM_gmin  ,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             dscl(n,k,l)=(&
                  +flux( n              , ADM_AI  ) &
                  +flux( n              , ADM_AIJ ) &
                  +flux( n              , ADM_AJ  ) &
                  -flux( n-1            , ADM_AI  ) &
                  -flux( n-1-ADM_gall_1d, ADM_AIJ ) &
                  -flux( n  -ADM_gall_1d, ADM_AJ ) &
                  ) * cmdif_P( n, l)&
                  * fact
          end do
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             dscl(suf(ADM_gmin,ADM_gmin),k,l)=(&
                  +flux(suf(ADM_gmin,ADM_gmin),ADM_AI)&
                  +flux(suf(ADM_gmin,ADM_gmin),ADM_AIJ)&
                  +flux(suf(ADM_gmin,ADM_gmin),ADM_AJ)&
                  -flux(suf(ADM_gmin-1,ADM_gmin),ADM_AI)&
                  -flux(suf(ADM_gmin-1,ADM_gmin-1),ADM_AIJ)&
                  ) * GMTR_P_var(suf(ADM_gmin,ADM_gmin),ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          end if
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do k= 1, ADM_kall
          do l=1,ADM_lall_pl
             flux_pl(ADM_gmin_pl)&
                  =((vxt_pl(ADM_GMAX_PL,k,l)+vxt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HNX)&
                  +(vyt_pl(ADM_GMAX_PL,k,l)+vyt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HNY)&
                  +(vzt_pl(ADM_GMAX_PL,k,l)+vzt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HNZ) )*0.5D0
             flux_pl(ADM_gmin_pl) = flux_pl(ADM_gmin_pl)&
                  * ( kh_pl(ADM_gslf_pl,k,l)+kh_pl(ADM_gmin_pl,k,l) )*0.5D0
             do n=ADM_gmin_pl+1,ADM_GMAX_PL
                flux_pl(n)&
                     =((vxt_pl(n-1,k,l)+vxt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNX)&
                     +(vyt_pl(n-1,k,l)+vyt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNY)&
                     +(vzt_pl(n-1,k,l)+vzt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNZ) )*0.5D0
                flux_pl(n) = flux_pl(n)&
                     * ( kh_pl(ADM_gslf_pl,k,l)+kh_pl(n,k,l) )*0.5D0
             enddo
             !
             dscl_pl(ADM_gslf_pl,k,l)=(&
                  +flux_pl(ADM_gmin_pl  )&
                  +flux_pl(ADM_gmin_pl+1)&
                  +flux_pl(ADM_gmin_pl+2)&
                  +flux_pl(ADM_gmin_pl+3)&
                  +flux_pl(ADM_gmin_pl+4)&
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          enddo
       enddo
    end if



    call DEBUG_rapend('++++OPRT_diffusion')


  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------  
  subroutine OPRT_horizontalize_vec( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
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
    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_horizontalize_vec')

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       prd = ( vx(g,k,l) * GRD_e(g,l,GRD_XDIR) &
             + vy(g,k,l) * GRD_e(g,l,GRD_YDIR) &
             + vz(g,k,l) * GRD_e(g,l,GRD_ZDIR) )

       vx(g,k,l) = vx(g,k,l) - prd * GRD_e(g,l,GRD_XDIR)
       vy(g,k,l) = vy(g,k,l) - prd * GRD_e(g,l,GRD_YDIR)
       vz(g,k,l) = vz(g,k,l) - prd * GRD_e(g,l,GRD_ZDIR)
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          prd = ( vx_pl(g,k,l) * GRD_e_pl(g,l,GRD_XDIR) &
                + vy_pl(g,k,l) * GRD_e_pl(g,l,GRD_YDIR) &
                + vz_pl(g,k,l) * GRD_e_pl(g,l,GRD_ZDIR) )

          vx_pl(g,k,l) = vx_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_XDIR)
          vy_pl(g,k,l) = vy_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_YDIR)
          vz_pl(g,k,l) = vz_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_ZDIR)
       enddo
       enddo
       enddo
    endif

    call DEBUG_rapend('++++OPRT_horizontalize_vec')

    return
  end subroutine OPRT_horizontalize_vec

  !-----------------------------------------------------------------------------
  subroutine OPRT_vorticity(       &
       scl, scl_pl,                &
       vx, vx_pl,                  &
       vy, vy_pl,                  &
       vz, vz_pl,                  &
       mfact )
    !
    use mod_adm, only :   &
         !--- public parameters
         ADM_W,           &
         ADM_prc_pl,      &
         ADM_TI,          &
         ADM_TJ,          &
         ADM_AI,          &
         ADM_AIJ,         &
         ADM_AJ,          &
         ADM_KNONE,       &
         ADM_lall_pl,     &
         ADM_gmin_pl,     &
         ADM_GMAX_PL,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         !--- public variables
         ADM_prc_me,      &
         ADM_prc_tab,     &
         ADM_rgn_vnum,    &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall
    use mod_gmtr, only :  &
         !--- public parameters
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_A_htx,      &
         GMTR_A_hty,      &
         GMTR_A_htz,      &
         GMTR_P_rarea,    &
         !--- public variables
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    !
    implicit none
    !
    real(8), intent(inout) :: scl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in)  :: vx(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in), optional :: mfact
    !
    real(8)  :: vxt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vxt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: vyt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vyt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: vzt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vzt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: flux(        &
         ADM_gall,           &
         ADM_AI:ADM_AJ)
    real(8)  :: flux_pl(     &
         ADM_gall_pl)
    !
    integer :: l,n,k
    integer :: rgnid
    real(8) :: fact
    !
    integer :: nstart,nend
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    call DEBUG_rapstart('++++OPRT_vorticity')

    if(present(mfact)) then
       fact=mfact
    else
       fact=1.0D0
    end if

    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k=ADM_kmin,ADM_kmax
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend = suf(ADM_gmax,ADM_gmax)
          do n = nstart,nend
             vxt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *vx(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *vx(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *vx(n+1+ADM_gall_1d,k,l)
             vxt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *vx(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *vx(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *vx(n+ADM_gall_1d,k,l)
             vyt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *vy(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *vy(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *vy(n+1+ADM_gall_1d,k,l)
             vyt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *vy(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *vy(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *vy(n+ADM_gall_1d,k,l)
             vzt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *vz(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *vz(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *vz(n+1+ADM_gall_1d,k,l)
             vzt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *vz(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *vz(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *vz(n+ADM_gall_1d,k,l)
          enddo
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             vxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do k=ADM_kmin,ADM_kmax
          do l=1,ADM_lall_pl
             do n=ADM_gmin_pl,ADM_GMAX_PL-1
                vxt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *vx_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *vx_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *vx_pl(n+1,k,l)
                vyt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *vy_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *vy_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *vy_pl(n+1,k,l)
                vzt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *vz_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *vz_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *vz_pl(n+1,k,l)
             enddo
             vxt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *vx_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *vx_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *vx_pl(ADM_gmin_pl,k,l)
             vyt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *vy_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *vy_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *vy_pl(ADM_gmin_pl,k,l)
             vzt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *vz_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *vz_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *vz_pl(ADM_gmin_pl,k,l)
          enddo
       enddo
    end if
    !
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k=ADM_kmin,ADM_kmax
          !
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,ADM_AI)&
                  =((vxt(n-ADM_gall_1d,k,l,ADM_TJ)+vxt(n,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HTX)&
                  +(vyt(n-ADM_gall_1d,k,l,ADM_TJ)+vyt(n,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HTY)&
                  +(vzt(n-ADM_gall_1d,k,l,ADM_TJ)+vzt(n,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HTZ))*0.5D0
          enddo
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,ADM_AIJ)&
                  =((vxt(n,k,l,ADM_TI)+vxt(n,k,l,ADM_TJ))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HTX)&
                  +(vyt(n,k,l,ADM_TI)+vyt(n,k,l,ADM_TJ))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HTY)&
                  +(vzt(n,k,l,ADM_TI)+vzt(n,k,l,ADM_TJ))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HTZ))*0.5D0
          enddo
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,ADM_AJ)&
                  =((vxt(n,k,l,ADM_TJ)+vxt(n-1,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HTX)&
                  +(vyt(n,k,l,ADM_TJ)+vyt(n-1,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HTY)&
                  +(vzt(n,k,l,ADM_TJ)+vzt(n-1,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HTZ))*0.5D0
          enddo
          !
          nstart = suf(ADM_gmin  ,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             scl(n,k,l)=-(&
                  +flux(n,ADM_AI)&
                  +flux(n,ADM_AIJ)&
                  +flux(n,ADM_AJ)&
                  -flux(n-1,ADM_AI)&
                  -flux(n-1-ADM_gall_1d,ADM_AIJ)&
                  -flux(n-ADM_gall_1d,ADM_AJ)&
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          enddo
          !
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             scl(suf(ADM_gmin,ADM_gmin),k,l)=-(&
                  +flux(suf(ADM_gmin,ADM_gmin),ADM_AI)&
                  +flux(suf(ADM_gmin,ADM_gmin),ADM_AIJ)&
                  +flux(suf(ADM_gmin,ADM_gmin),ADM_AJ)&
                  -flux(suf(ADM_gmin-1,ADM_gmin),ADM_AI)&
                  -flux(suf(ADM_gmin-1,ADM_gmin-1),ADM_AIJ)&
                  ) * GMTR_P_var(suf(ADM_gmin,ADM_gmin),ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          end if
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          do l=1,ADM_lall_pl
             !
             flux_pl(ADM_gmin_pl)&
                  =((vxt_pl(ADM_GMAX_PL,k,l)+vxt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HTX)&
                  +(vyt_pl(ADM_GMAX_PL,k,l)+vyt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HTY)&
                  +(vzt_pl(ADM_GMAX_PL,k,l)+vzt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HTZ) )*0.5D0
             do n=ADM_gmin_pl+1,ADM_GMAX_PL
                flux_pl(n)&
                     =((vxt_pl(n-1,k,l)+vxt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HTX)&
                     +(vyt_pl(n-1,k,l)+vyt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HTY)&
                     +(vzt_pl(n-1,k,l)+vzt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HTZ) )*0.5D0
             enddo
             !
             scl_pl(ADM_gslf_pl,k,l)=-(&
                  +flux_pl(ADM_gmin_pl  )&
                  +flux_pl(ADM_gmin_pl+1)&
                  +flux_pl(ADM_gmin_pl+2)&
                  +flux_pl(ADM_gmin_pl+3)&
                  +flux_pl(ADM_gmin_pl+4)&
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          enddo
       enddo
    end if

    call DEBUG_rapend('++++OPRT_vorticity')

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
       ADM_W,        &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_AI,       &
       ADM_AIJ,      &
       ADM_AJ,       &
       ADM_prc_tab,  &
       ADM_prc_me,   &
       ADM_prc_pl,   &
       ADM_rgn_vnum, &
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
       ADM_KNONE,    &
       ADM_kmin,     &
       ADM_kmax
    use mod_gmtr, only: &
       GMTR_P_RAREA,  &
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
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    implicit none

    real(8), intent(inout) :: grdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: grdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: grdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: sclt   (ADM_gall,   ADM_kall,ADM_TI:ADM_TJ)
    real(8) :: sclt_pl(ADM_gall_pl,ADM_kall)
    real(8) :: fact

    integer :: rgnid
    integer :: nstart, nend

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: k, l, n, v, k0

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: TI,TJ,AI,AIJ,AJ,TNX,TNY,TNZ,TN2X,TN2Y,TN2Z,HNX,HNY,HNZ
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_divdamp')

    k0   = ADM_KNONE
    TI   = ADM_TI
    TJ   = ADM_TJ
    AI   = ADM_AI
    AIJ  = ADM_AIJ
    AJ   = ADM_AJ
    TNX  = GMTR_A_TNX
    TNY  = GMTR_A_TNY
    TNZ  = GMTR_A_TNZ
    HNX  = GMTR_A_HNX
    HNY  = GMTR_A_HNY
    HNZ  = GMTR_A_HNZ
    TN2X = GMTR_A_TN2X
    TN2Y = GMTR_A_TN2Y
    TN2Z = GMTR_A_TN2Z

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax,  ADM_gmax  )

       do k = ADM_kmin, ADM_kmax
          do n = nstart, nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d

             sclt(n,k,ADM_TI) = ( - ( vx(ij    ,k,l) + vx(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNX) &
                                  - ( vy(ij    ,k,l) + vy(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNY) &
                                  - ( vz(ij    ,k,l) + vz(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNZ) &
                                  - ( vx(ip1j  ,k,l) + vx(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNX) &
                                  - ( vy(ip1j  ,k,l) + vy(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNY) &
                                  - ( vz(ip1j  ,k,l) + vz(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNZ) &
                                  + ( vx(ip1jp1,k,l) + vx(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNX) &
                                  + ( vy(ip1jp1,k,l) + vy(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNY) &
                                  + ( vz(ip1jp1,k,l) + vz(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNZ) &
                                ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TI,GMTR_T_RAREA)

             sclt(n,k,ADM_TJ) = ( - ( vx(ij    ,k,l) + vx(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNX) &
                                  - ( vy(ij    ,k,l) + vy(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNY) &
                                  - ( vz(ij    ,k,l) + vz(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNZ) &
                                  + ( vx(ip1jp1,k,l) + vx(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNX) &
                                  + ( vy(ip1jp1,k,l) + vy(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNY) &
                                  + ( vz(ip1jp1,k,l) + vz(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNZ) &
                                  + ( vx(ijp1  ,k,l) + vx(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNX) &
                                  + ( vy(ijp1  ,k,l) + vy(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNY) &
                                  + ( vz(ijp1  ,k,l) + vz(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNZ) &
                                ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TJ,GMTR_T_RAREA)
          enddo
       enddo

       nstart = suf(ADM_gmin,ADM_gmin)

       do k = ADM_kmin, ADM_kmax
          do n = nstart, nend
             ij     = n
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

             grdx(n,k,l) = ( + ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) * GMTR_A_var(ij,    k0,l,AI, HNX) &
                             + ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) * GMTR_A_var(ij,    k0,l,AIJ,HNX) &
                             + ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(ij,    k0,l,AJ, HNX) &
                             - ( sclt(im1jm1,k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(im1j,  k0,l,AI, HNX) &
                             - ( sclt(im1jm1,k,TI) + sclt(im1jm1,k,TJ) ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
                             - ( sclt(ijm1  ,k,TJ) + sclt(im1jm1,k,TI) ) * GMTR_A_var(ijm1,  k0,l,AJ, HNX) &
                           ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             grdy(n,k,l) = ( + ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) * GMTR_A_var(ij,    k0,l,AI, HNY) &
                             + ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) * GMTR_A_var(ij,    k0,l,AIJ,HNY) &
                             + ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(ij,    k0,l,AJ, HNY) &
                             - ( sclt(im1jm1,k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(im1j,  k0,l,AI, HNY) &
                             - ( sclt(im1jm1,k,TI) + sclt(im1jm1,k,TJ) ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
                             - ( sclt(ijm1  ,k,TJ) + sclt(im1jm1,k,TI) ) * GMTR_A_var(ijm1,  k0,l,AJ, HNY) &
                           ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             grdz(n,k,l) = ( + ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) * GMTR_A_var(ij,    k0,l,AI, HNZ) &
                             + ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) * GMTR_A_var(ij,    k0,l,AIJ,HNZ) &
                             + ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(ij,    k0,l,AJ, HNZ) &
                             - ( sclt(im1jm1,k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(im1j,  k0,l,AI, HNZ) &
                             - ( sclt(im1jm1,k,TI) + sclt(im1jm1,k,TJ) ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
                             - ( sclt(ijm1  ,k,TJ) + sclt(im1jm1,k,TI) ) * GMTR_A_var(ijm1,  k0,l,AJ, HNZ) &
                           ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact
          enddo
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
          n = suf(ADM_gmin,ADM_gmin)

          ij     = n
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          do k = ADM_kmin, ADM_kmax
             sclt(im1jm1,k,ADM_TI) = sclt(ijm1,k,ADM_TJ) ! copy

             grdx(n,k,l) = ( + ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) * GMTR_A_var(ij,    k0,l,AI, HNX) &
                             + ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) * GMTR_A_var(ij,    k0,l,AIJ,HNX) &
                             + ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(ij,    k0,l,AJ, HNX) &
                             - ( sclt(im1jm1,k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(im1j,  k0,l,AI, HNX) &
                             - ( sclt(im1jm1,k,TI) + sclt(im1jm1,k,TJ) ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
                           ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             grdy(n,k,l) = ( + ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) * GMTR_A_var(ij,    k0,l,AI, HNY) &
                             + ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) * GMTR_A_var(ij,    k0,l,AIJ,HNY) &
                             + ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(ij,    k0,l,AJ, HNY) &
                             - ( sclt(im1jm1,k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(im1j,  k0,l,AI, HNY) &
                             - ( sclt(im1jm1,k,TI) + sclt(im1jm1,k,TJ) ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
                           ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

             grdz(n,k,l) = ( + ( sclt(ijm1,  k,TJ) + sclt(ij,    k,TI) ) * GMTR_A_var(ij,    k0,l,AI, HNZ) &
                             + ( sclt(ij,    k,TI) + sclt(ij,    k,TJ) ) * GMTR_A_var(ij,    k0,l,AIJ,HNZ) &
                             + ( sclt(ij,    k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(ij,    k0,l,AJ, HNZ) &
                             - ( sclt(im1jm1,k,TJ) + sclt(im1j,  k,TI) ) * GMTR_A_var(im1j,  k0,l,AI, HNZ) &
                             - ( sclt(im1jm1,k,TI) + sclt(im1jm1,k,TJ) ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
                           ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact
          enddo
       endif

       grdx   (:,ADM_kmin-1,l) = 0.D0
       grdx   (:,ADM_kmax+1,l) = 0.D0
       grdy   (:,ADM_kmin-1,l) = 0.D0
       grdy   (:,ADM_kmax+1,l) = 0.D0
       grdz   (:,ADM_kmin-1,l) = 0.D0
       grdz   (:,ADM_kmax+1,l) = 0.D0
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
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
                               ) * 0.5D0 * GMTR_T_var_pl(ij,k0,l,GMTR_T_RAREA)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax
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

             grdx_pl(n,k,l) = grdx_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA) * fact
             grdy_pl(n,k,l) = grdy_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA) * fact
             grdz_pl(n,k,l) = grdz_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA) * fact
          enddo

       enddo
    endif

    call DEBUG_rapend('++++OPRT_divdamp')

    return
  end subroutine OPRT_divdamp

end module mod_oprt
!-------------------------------------------------------------------------------
