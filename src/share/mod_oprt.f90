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
  public :: OPRT_divergence2_prep
  public :: OPRT_divergence2

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
  integer,allocatable,private,save::n021(:)

  ! < for divergence operator >
  real(8), private, allocatable, save :: cdiv   (:,:,:,:)
  real(8), private, allocatable, save :: cdiv_pl(:,:,:,:)

  ! < for gradient operator >
  real(8), private, allocatable, save :: cgrad(:,:,:,:)!(0:6,gall,lall,1:3)
  real(8), private, allocatable, save :: cgrad0_pl(:,:,:)
  real(8), private, allocatable, save :: cgrad1_pl(:,:,:)
  real(8), private, allocatable, save :: cgrad2_pl(:,:,:)
  real(8), private, allocatable, save :: cgrad3_pl(:,:,:)
  real(8), private, allocatable, save :: cgrad4_pl(:,:,:)
  real(8), private, allocatable, save :: cgrad5_pl(:,:,:)

  ! < for laplacian operator >
  real(8), private, allocatable, save :: clap0(:,:)
  real(8), private, allocatable, save :: clap1(:,:)
  real(8), private, allocatable, save :: clap2(:,:)
  real(8), private, allocatable, save :: clap3(:,:)
  real(8), private, allocatable, save :: clap4(:,:)
  real(8), private, allocatable, save :: clap5(:,:)
  real(8), private, allocatable, save :: clap6(:,:)

  real(8), private, allocatable, save :: clap0_pl(:,:)
  real(8), private, allocatable, save :: clap1_pl(:,:)
  real(8), private, allocatable, save :: clap2_pl(:,:)
  real(8), private, allocatable, save :: clap3_pl(:,:)
  real(8), private, allocatable, save :: clap4_pl(:,:)
  real(8), private, allocatable, save :: clap5_pl(:,:)

  ! < for diffusion operator >
  real(8), public,  allocatable, save :: cmdif_T(:,:,:)! (TI:TJ,n,l) <- GMTR_T_VAR(n,1,l,TI:TJ,T_RAREA)
  real(8), allocatable, public,  save  :: cmdif_AH(:,:,:,:)!(AI:AJ,1:3,n,l) <-GMTR_A_VAR(n,1,l,TI:TJ,HN[XYZ])
  real(8), allocatable, public,  save  :: cmdif_AT(:,:,:,:)!(AI:AJ,1:3,n,l) <-GMTR_A_VAR(n,1,l,TI:TJ,TN[XYZ])
  real(8), allocatable, public,  save  :: cmdif_P(:,:) !(n,l) <- GMTR_P_VAR(n,1,l,P_RAREA)

  !
  ! <  For XTMS (when it is used, above clap0-5_pl,cgrad0-5_pl,cdiv0-5_pl are 
  !            not needed)
  real(8), private, allocatable, save :: cdivN_pl(:,:,:,:) !S.iga100608
  real(8), private, allocatable, save :: cgradN_pl(:,:,:,:)!S.iga100608
  real(8), private, allocatable, save :: clapN_pl(:,:,:)   !S.iga100608
  real(8),private::tmp1 !S.iga100608
  integer,private::itmp !S.iga100608

  ! < for OPRT_divergence2_prep >  ! Y.Niwa add 080130
  real(8), private, allocatable, save :: local_t_var(:,:,:,:,:)
  real(8), private, allocatable, save :: local_t_var_pl(:,:,:,:)
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
    integer :: k0,a0
    integer :: tx1,ty1,tz1
    integer :: tx2,ty2,tz2
    integer :: hx1,hy1,hz1

    integer :: ij
    integer :: im1j,ijm1,im1jm1
    integer :: ip1j,ijp1,ip1jp1

    integer :: rgnid
    integer :: nstart, nend
    integer :: n, k, l, m, md, v

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: TI,TJ,AI,AIJ,AJ,W1,W2,W3
    !---------------------------------------------------------------------------

    TI  = ADM_TI
    TJ  = ADM_TJ
    AI  = ADM_AI
    AIJ = ADM_AIJ
    AJ  = ADM_AJ
    W1  = GMTR_T_W1
    W2  = GMTR_T_W2
    W3  = GMTR_T_W3

    k0 = ADM_KNONE

    allocate( n021((ADM_gall_1d-2)**2) )

    do n = 1, (ADM_gall_1d-2)**2
      n021(n) = ADM_gall_1d * ( (n-1)/(ADM_gall_1d-2)+1) + mod(n-1,ADM_gall_1d-2)+2
    enddo

    !---< setup coefficient of divergence operator >
    allocate( cdiv   (0:6,             ADM_gall   ,ADM_lall   ,3) )
    allocate( cdiv_pl(0:ADM_vlink_nmax,ADM_gall_pl,ADM_lall_pl,3) )

    nstart = suf(ADM_gmin,ADM_gmin)
    nend   = suf(ADM_gmax,ADM_gmax)

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do m = 1, 3
          md = m + GMTR_A_HNX - 1
          do n = nstart, nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

             ! ij
             cdiv(0,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,GMTR_T_W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cdiv(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cdiv(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cdiv(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0*GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cdiv(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cdiv(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,GMTR_T_W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cdiv(6,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,GMTR_T_W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
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
             cdiv(0,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ip1j
             cdiv(1,n,l,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cdiv(2,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ijp1
             cdiv(3,n,l,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! im1j
             cdiv(4,n,l,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,GMTR_T_W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! im1jm1
             cdiv(5,n,l,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ijm1
             cdiv(6,n,l,m) = ( + GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,GMTR_T_W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
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

                cdiv_pl(0,n,l,m) = cdiv_pl(0,n,l,m) + ( GMTR_T_var_pl(ij,k0,l,GMTR_T_W1) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                                      + GMTR_T_var_pl(ij,k0,l,GMTR_T_W1) * GMTR_A_var_pl(ijp1,k0,l,md) )
             enddo
             cdiv_pl(0,n,l,m) = cdiv_pl(0,n,l,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cdiv_pl(v-1,n,l,m) = ( + GMTR_T_var_pl(ijm1,k0,l,GMTR_T_W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                       + GMTR_T_var_pl(ijm1,k0,l,GMTR_T_W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,GMTR_T_W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,GMTR_T_W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                     ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif



    ! ---- setup coefficient of gradient operator
    allocate(cgrad(0:6,ADM_gall,ADM_lall,1:3))
    !
    allocate(cgrad0_pl(ADM_gall_pl,ADM_lall_pl,1:3))
    allocate(cgrad1_pl(ADM_gall_pl,ADM_lall_pl,1:3))
    allocate(cgrad2_pl(ADM_gall_pl,ADM_lall_pl,1:3))
    allocate(cgrad3_pl(ADM_gall_pl,ADM_lall_pl,1:3))
    allocate(cgrad4_pl(ADM_gall_pl,ADM_lall_pl,1:3))
    allocate(cgrad5_pl(ADM_gall_pl,ADM_lall_pl,1:3))
    !
    nstart=suf(ADM_gmin,ADM_gmin)
    nend  =suf(ADM_gmax,ADM_gmax)
    do l=1,ADM_lall
      do m=1,3
        md=m+GMTR_A_hnx-1
        do n = nstart,nend
          ij=n
          ip1j=ij+1
          ip1jp1=ij+1+ADM_gall_1d
          ijp1=ij+ADM_gall_1d
          im1j=ij-1
          im1jm1=ij-1-ADM_gall_1d
          ijm1=ij-ADM_gall_1d
          ! 
          ! ij
          cgrad(0,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W1) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W1) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W2) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W2) &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W2) &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TI,GMTR_T_W3) &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W2) &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3) &
                        *GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TI,GMTR_T_W3) &
                        *GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      -2*GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      -2*GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      -2*GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      +2*GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      +2*GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      +2*GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ip1j
          cgrad(1,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W2) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W2) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2) &
                        *GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ip1jp1
          cgrad(2,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W3) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W3) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ijp1
          cgrad(3,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W3) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W3) &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! im1j
          cgrad(4,n,l,m)=(GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W1) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W1) &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W3) &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W3) &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! im1jm1
          cgrad(5,n,l,m)=(-GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W1) &
                         *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                         -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TI,GMTR_T_W1) &
                         *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                         -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TI,GMTR_T_W1) &
                         *GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                         -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W1) &
                         *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ijm1
          cgrad(6,n,l,m)=(GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1) &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TI,GMTR_T_W2) &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1) &
                        *GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TI,GMTR_T_W2) &
                        *GMTR_A_var(ijm1  ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0    
        enddo !loop n
      enddo !loop m
    enddo !loop l
    !
    do l=1,ADM_lall
      rgnid=ADM_prc_tab(l,ADM_prc_me)
      if (ADM_rgn_vnum(ADM_W,rgnid)==3) then
        n=suf(ADM_gmin,ADM_gmin)
        ij=n
        ip1j=ij+1
        ip1jp1=ij+1+ADM_gall_1d
        ijp1=ij+ADM_gall_1d
        im1j=ij-1
        im1jm1=ij-1-ADM_gall_1d
        ijm1=ij-ADM_gall_1d
        do m=1,3
          md=m+GMTR_A_hnx-1
          ! ij
          cgrad(0,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W1)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W1)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W2)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W2)                                   &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)                                   &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)                                   &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)                                   &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      -2*GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      -2*GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      -2*GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      +2*GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                      +2*GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ip1j
          cgrad(1,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W2)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W2)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)                                   &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ip1jp1
          cgrad(2,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W3)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TI,GMTR_T_W3)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ijp1
          cgrad(3,n,l,m)=(GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(ij    ,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        +GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W3)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W3)                                   &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! im1j
          cgrad(4,n,l,m)=(GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W1)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AJ ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1j  ,ADM_KNONE,l,ADM_TI,GMTR_T_W1)                                   &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)                                   &
                        *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)                                   &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! im1jm1
          cgrad(5,n,l,m)=(-GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)                                   &
                         *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                         -GMTR_T_var(im1jm1,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)                                   &
                         *GMTR_A_var(im1j  ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
          ! ijm1
          cgrad(6,n,l,m)=(GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)                                   &
                        *GMTR_A_var(ij    ,ADM_KNONE,l,ADM_AI ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea) &
                        -GMTR_T_var(ijm1  ,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)                                   &
                        *GMTR_A_var(im1jm1,ADM_KNONE,l,ADM_AIJ,md)*GMTR_P_var(ij,ADM_KNONE,l,GMTR_P_rarea))*0.5d0
        enddo
      endif
    enddo
    !
    if (ADM_prc_me==ADM_prc_pl) then
      n=ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4
      do l=1,ADM_lall_pl
        do m=1,3
          md=m+GMTR_A_hnx-1
          !n
          cgrad0_pl(n,l,m)=(                                                   &
                            +2*(GMTR_T_var_pl(n4,ADM_KNONE,l,GMTR_T_W1)-1.0d0) &
                            *GMTR_A_var_pl(n0,ADM_KNONE,l,md)                  &
                            +2*(GMTR_T_var_pl(n0,ADM_KNONE,l,GMTR_T_W1)-1.0d0) &
                            *GMTR_A_var_pl(n1,ADM_KNONE,l,md)                  &
                            +2*(GMTR_T_var_pl(n1,ADM_KNONE,l,GMTR_T_W1)-1.0d0) &
                            *GMTR_A_var_pl(n2,ADM_KNONE,l,md)                  &
                            +2*(GMTR_T_var_pl(n2,ADM_KNONE,l,GMTR_T_W1)-1.0d0) &
                            *GMTR_A_var_pl(n3,ADM_KNONE,l,md)                  &
                            +2*(GMTR_T_var_pl(n3,ADM_KNONE,l,GMTR_T_W1)-1.0d0) &
                            *GMTR_A_var_pl(n4,ADM_KNONE,l,md)                  &
                           )*GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_rarea)*0.5d0                                           
          !n0
          cgrad1_pl(n,l,m)=(                                         &
                            +GMTR_T_var_pl(n4,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n4,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n4,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n0,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n0,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n0,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n0,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n1,ADM_KNONE,l,md)        &
                           )*GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_rarea)*0.5d0                                 
          !n1
          cgrad2_pl(n,l,m)=(                                         &
                            +GMTR_T_var_pl(n0,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n0,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n0,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n1,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n1,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n1,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n1,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n2,ADM_KNONE,l,md)        &
                           )*GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_rarea)*0.5d0                                 
          !n2
          cgrad3_pl(n,l,m)=(                                         &
                            +GMTR_T_var_pl(n1,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n1,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n1,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n2,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n2,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n2,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n2,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n3,ADM_KNONE,l,md)        &
                           )*GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_rarea)*0.5d0                                 
          !n3
          cgrad4_pl(n,l,m)=(                                         &
                            +GMTR_T_var_pl(n2,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n2,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n2,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n3,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n3,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n3,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n3,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n4,ADM_KNONE,l,md)        &
                           )*GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_rarea)*0.5d0                                 
          !n4
          cgrad5_pl(n,l,m)=(                                         &
                            +GMTR_T_var_pl(n3,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n3,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n3,ADM_KNONE,l,GMTR_T_W3) &
                            *GMTR_A_var_pl(n4,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n4,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n4,ADM_KNONE,l,md)        &
                            +GMTR_T_var_pl(n4,ADM_KNONE,l,GMTR_T_W2) &
                            *GMTR_A_var_pl(n0,ADM_KNONE,l,md)        &
                           )*GMTR_P_var_pl(n,ADM_KNONE,l,GMTR_P_rarea)*0.5d0                                 
        enddo
      enddo
    endif
    !
    ! ---- setup coefficient of laplacian operator
    allocate(clap0(ADM_gall,ADM_lall))
    allocate(clap1(ADM_gall,ADM_lall))
    allocate(clap2(ADM_gall,ADM_lall))
    allocate(clap3(ADM_gall,ADM_lall))
    allocate(clap4(ADM_gall,ADM_lall))
    allocate(clap5(ADM_gall,ADM_lall))
    allocate(clap6(ADM_gall,ADM_lall))
    !
    k0=ADM_KNONE
    a0=GMTR_T_rarea
    tx1=GMTR_A_tnx
    ty1=GMTR_A_tny
    tz1=GMTR_A_tnz
    hx1=GMTR_A_hnx
    hy1=GMTR_A_hny
    hz1=GMTR_A_hnz
    nstart=suf(ADM_gmin,ADM_gmin)
    nend  =suf(ADM_gmax,ADM_gmax)
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=nstart,nend
          ij    =n
          ip1j  =ij+1
          ip1jp1=ij+1+ADM_gall_1d
          ijp1  =ij+ADM_gall_1d
          im1j  =ij-1
          im1jm1=ij-1-ADM_gall_1d
          ijm1  =ij-ADM_gall_1d
          !
     ! 0: ij
     clap0(ij,l)=( &
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
     clap0(ij,l)=clap0(ij,l)+( &
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
     ! 1: ip1j
     clap1(ij,l)=( &
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
     ! 2: ip1jp1
     clap2(ij,l)=( &
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
     ! 3: ijp1
     clap3(ij,l)=( &
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
     ! 4: im1j
     clap4(ij,l)=( &
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
     ! 5: im1jm1
     clap5(ij,l)=( &
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
     ! 6: ijm1
     clap6(ij,l)=( &
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
      enddo
    enddo
    !
    do l=1,ADM_lall
      rgnid=ADM_prc_tab(l,ADM_prc_me)
      if (ADM_rgn_vnum(ADM_W,rgnid)==3) then
        n=suf(ADM_gmin,ADM_gmin)
        ij    =n
        ip1j  =ij+1
        ip1jp1=ij+1+ADM_gall_1d
        ijp1  =ij+ADM_gall_1d
        im1j  =ij-1
        im1jm1=ij-1-ADM_gall_1d
        ijm1  =ij-ADM_gall_1d
        !
        do k= 1,ADM_kall
     ! 0: ij
     clap0(ij,l)=( &
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
     clap0(ij,l)=clap0(ij,l)+( &                 ! Y.Niwa add 06/08/22 
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
     ! 1: ip1j
     clap1(ij,l)=( &
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
     ! 2: ip1jp1
     clap2(ij,l)=( &
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
     ! 3: ijp1
     clap3(ij,l)=( &
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
     ! 4: im1j
     clap4(ij,l)=( &
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
     ! 5: im1jm1
     clap5(ij,l)=( &
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
     ! 6: ijm1
     clap6(ij,l)=( &
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
        !
        enddo
      endif
    enddo
    !
    allocate(clap0_pl(ADM_gall_pl,ADM_lall_pl))
    allocate(clap1_pl(ADM_gall_pl,ADM_lall_pl))
    allocate(clap2_pl(ADM_gall_pl,ADM_lall_pl))
    allocate(clap3_pl(ADM_gall_pl,ADM_lall_pl))
    allocate(clap4_pl(ADM_gall_pl,ADM_lall_pl))
    allocate(clap5_pl(ADM_gall_pl,ADM_lall_pl))
    !
    if (ADM_prc_me==ADM_prc_pl) then
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
      do l=1,ADM_lall_pl
          !
          ! n
          clap0_pl(n,l)=( &  
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
          clap0_pl(n,l)= clap0_pl(n,l) + ( &                        ! Y.Niwa add 060822 
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
          clap1_pl(n,l)=( &  
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
          clap2_pl(n,l)=( &  
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
          clap3_pl(n,l)=( &  
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
          clap4_pl(n,l)=( &  
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
          clap5_pl(n,l)=( &  
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
    allocate(cmdif_T(ADM_TI:ADM_TJ,ADM_gall,ADM_lall))
    allocate(cmdif_AT(ADM_AI:ADM_AJ,1:3,ADM_gall,ADM_lall))
    allocate(cmdif_AH(ADM_AI:ADM_AJ,1:3,ADM_gall,ADM_lall))
    allocate(cmdif_P(ADM_gall,ADM_lall))

    do m=ADM_TI, ADM_TJ
       do l=1,ADM_lall
          do n=1,ADM_gall
             cmdif_T(m, n, l) = GMTR_T_var(n, ADM_KNONE, l, m, GMTR_T_RAREA)
          enddo
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

    do l=1,ADM_lall
       do n=1,ADM_gall
          cmdif_P( n, l ) = GMTR_P_var( n,  ADM_KNONE,  l,    GMTR_P_RAREA  )
       enddo
    enddo
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
       ADM_gmin,    &
       ADM_gmax,    &
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

    integer :: nstart, nend
    integer :: n, k, l, v

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    nstart = suf(ADM_gmin,ADM_gmin)
    nend   = suf(ADM_gmax,ADM_gmax)

    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax
    do n = nstart, nend
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
          scl_pl(n,k,l) = 0

          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( cdiv_pl(v-1,n,l,1) * vx_pl(v,k,l) &
                                             + cdiv_pl(v-1,n,l,2) * vy_pl(v,k,l) &
                                             + cdiv_pl(v-1,n,l,3) * vz_pl(v,k,l) )
          enddo

          scl_pl(n,k,l) = scl_pl(n,k,l) * fact
       enddo
       enddo
    endif

    return
  end subroutine OPRT_divergence

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient(        &
       vx, vx_pl,                  &
       vy, vy_pl,                  &
       vz, vz_pl,                  &
       scl, scl_pl,                &
       mfact )
    use mod_adm, only :   &
         !--- public parameters
         ADM_prc_pl,      &
         ADM_lall_pl,     &
         ADM_gmin_pl,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         !--- public variables
         ADM_prc_me,      &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall
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

    integer :: im1j,ijm1,im1jm1
    integer :: ip1j,ijp1,ip1jp1
    integer :: ij
    integer :: n0,n1,n2,n3,n4
    !
    integer :: l,n,k
    integer :: nstart,nend
    integer :: rgnid
    real(8)::fact
    !
    integer :: suf,i,j
    suf(i,j)=ADM_gall_1d*((j)-1)+(i)

    if(present(mfact)) then
       fact=mfact
    else
       fact=1.0d0
    end if
    !
    nstart=suf(ADM_gmin,ADM_gmin)
    nend  =suf(ADM_gmax,ADM_gmax)
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=nstart,nend
          ij=n
          ip1j=ij+1
          ip1jp1=ij+1+ADM_gall_1d
          ijp1=ij+ADM_gall_1d
          im1j=ij-1
          im1jm1=ij-1-ADM_gall_1d
          ijm1=ij-ADM_gall_1d
          vx(n,k,l)=(                                &
                     +cgrad(0,ij,l,1)*scl(ij    ,k,l) &
                     +cgrad(1,ij,l,1)*scl(ip1j  ,k,l) &
                     +cgrad(2,ij,l,1)*scl(ip1jp1,k,l) &
                     +cgrad(3,ij,l,1)*scl(ijp1  ,k,l) &
                     +cgrad(4,ij,l,1)*scl(im1j  ,k,l) &
                     +cgrad(5,ij,l,1)*scl(im1jm1,k,l) &
                     +cgrad(6,ij,l,1)*scl(ijm1  ,k,l) &
                    )*fact
          vy(n,k,l)=(                                &
                     +cgrad(0,ij,l,2)*scl(ij    ,k,l) &
                     +cgrad(1,ij,l,2)*scl(ip1j  ,k,l) &
                     +cgrad(2,ij,l,2)*scl(ip1jp1,k,l) &
                     +cgrad(3,ij,l,2)*scl(ijp1  ,k,l) &
                     +cgrad(4,ij,l,2)*scl(im1j  ,k,l) &
                     +cgrad(5,ij,l,2)*scl(im1jm1,k,l) &
                     +cgrad(6,ij,l,2)*scl(ijm1  ,k,l) &
                    )*fact
          vz(n,k,l)=(                                &
                     +cgrad(0,ij,l,3)*scl(ij    ,k,l) &
                     +cgrad(1,ij,l,3)*scl(ip1j  ,k,l) &
                     +cgrad(2,ij,l,3)*scl(ip1jp1,k,l) &
                     +cgrad(3,ij,l,3)*scl(ijp1  ,k,l) &
                     +cgrad(4,ij,l,3)*scl(im1j  ,k,l) &
                     +cgrad(5,ij,l,3)*scl(im1jm1,k,l) &
                     +cgrad(6,ij,l,3)*scl(ijm1  ,k,l) &
                    )*fact
        enddo
      enddo
    enddo
    !
    if (ADM_prc_me==ADM_prc_pl) then
      n=ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4
      do l=1,ADM_lall_pl
        do k=1,ADM_kall
          vx_pl(n,k,l)=(                                &
                       +cgrad0_pl(n,l,1)*scl_pl(n ,k,l) & 
                       +cgrad1_pl(n,l,1)*scl_pl(n0,k,l) & 
                       +cgrad2_pl(n,l,1)*scl_pl(n1,k,l) & 
                       +cgrad3_pl(n,l,1)*scl_pl(n2,k,l) & 
                       +cgrad4_pl(n,l,1)*scl_pl(n3,k,l) & 
                       +cgrad5_pl(n,l,1)*scl_pl(n4,k,l) & 
                       )*fact
          vy_pl(n,k,l)=(                                &
                       +cgrad0_pl(n,l,2)*scl_pl(n ,k,l) & 
                       +cgrad1_pl(n,l,2)*scl_pl(n0,k,l) & 
                       +cgrad2_pl(n,l,2)*scl_pl(n1,k,l) & 
                       +cgrad3_pl(n,l,2)*scl_pl(n2,k,l) & 
                       +cgrad4_pl(n,l,2)*scl_pl(n3,k,l) & 
                       +cgrad5_pl(n,l,2)*scl_pl(n4,k,l) & 
                       )*fact
          vz_pl(n,k,l)=(                                &
                       +cgrad0_pl(n,l,3)*scl_pl(n ,k,l) & 
                       +cgrad1_pl(n,l,3)*scl_pl(n0,k,l) & 
                       +cgrad2_pl(n,l,3)*scl_pl(n1,k,l) & 
                       +cgrad3_pl(n,l,3)*scl_pl(n2,k,l) & 
                       +cgrad4_pl(n,l,3)*scl_pl(n3,k,l) & 
                       +cgrad5_pl(n,l,3)*scl_pl(n4,k,l) & 
                       )*fact
        enddo
      enddo
    endif

  end subroutine OPRT_gradient
  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian(dscl,dscl_pl,scl, scl_pl,mfact)
    use mod_adm, only :   &
         !--- public parameters
         ADM_prc_pl,      &
         ADM_TI,          &
         ADM_lall_pl,     &
         ADM_gmin_pl,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         !--- public variables
         ADM_prc_me,      &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall
    !
    implicit none
    !
    real(8),intent(inout)::dscl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8),intent(inout)::dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8),intent(in)   :: scl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8),intent(in)   :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8),intent(in),optional::mfact
    !
    integer :: im1j,ijm1,im1jm1
    integer :: ip1j,ijp1,ip1jp1
    integer :: ij
    integer :: n0,n1,n2,n3,n4
    integer :: nstart,nend
    integer :: l,n,k
    !
    real(8)::fact
    !
    integer ::suf,i,j
    suf(i,j)=ADM_gall_1d*((j)-1)+(i)

    if(present(mfact)) then
       fact=mfact
    else
       fact=1.0d0
    end if
    !
    nstart=suf(ADM_gmin,ADM_gmin)
    nend  =suf(ADM_gmax,ADM_gmax)
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=nstart,nend
          !
          ij    =n
          ip1j  =ij+1
          ip1jp1=ij+1+ADM_gall_1d
          ijp1  =ij+ADM_gall_1d
          im1j  =ij-1
          im1jm1=ij-1-ADM_gall_1d
          ijm1  =ij-ADM_gall_1d
          !
          dscl(n,k,l)=( &
                       +clap0(ij,l)*scl(ij    ,k,l) &
                       +clap1(ij,l)*scl(ip1j  ,k,l) &
                       +clap2(ij,l)*scl(ip1jp1,k,l) &
                       +clap3(ij,l)*scl(ijp1  ,k,l) &
                       +clap4(ij,l)*scl(im1j  ,k,l) &
                       +clap5(ij,l)*scl(im1jm1,k,l) &
                       +clap6(ij,l)*scl(ijm1  ,k,l) &
                      )*fact
        enddo
      enddo
    enddo
    !
    if (ADM_prc_me==ADM_prc_pl) then
      n =ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4
      do l=1,ADM_lall_pl
        do k=1,ADM_kall
          dscl_pl(n,k,l)=( &
                          +clap0_pl(n,l)*scl_pl(n ,k,l) &  
                          +clap1_pl(n,l)*scl_pl(n0,k,l) &  
                          +clap2_pl(n,l)*scl_pl(n1,k,l) &  
                          +clap3_pl(n,l)*scl_pl(n2,k,l) &  
                          +clap4_pl(n,l)*scl_pl(n3,k,l) &  
                          +clap5_pl(n,l)*scl_pl(n4,k,l) &  
                         )*fact
        enddo
      enddo
    endif

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
 
    !

    implicit none
    !
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

    if(present(mfact)) then
       fact=mfact
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
             end do
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
          end do
       end do
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
             end do
             !
             dscl_pl(ADM_gslf_pl,k,l)=(&
                  +flux_pl(ADM_gmin_pl  )&
                  +flux_pl(ADM_gmin_pl+1)&
                  +flux_pl(ADM_gmin_pl+2)&
                  +flux_pl(ADM_gmin_pl+3)&
                  +flux_pl(ADM_gmin_pl+4)&
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          end do
       end do
    end if





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

    real(8), intent(inout) :: vx   (ADM_gall,   ADM_kall,   ADM_lall)
    real(8), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: prd
    integer :: g, k, l
    !---------------------------------------------------------------------------

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
    !
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
          end do
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             vxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
       end do
    end do
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
             end do
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
          end do
       end do
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
          end do
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
          end do
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
          end do
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
          end do
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
       end do
    end do
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
             end do
             !
             scl_pl(ADM_gslf_pl,k,l)=-(&
                  +flux_pl(ADM_gmin_pl  )&
                  +flux_pl(ADM_gmin_pl+1)&
                  +flux_pl(ADM_gmin_pl+2)&
                  +flux_pl(ADM_gmin_pl+3)&
                  +flux_pl(ADM_gmin_pl+4)&
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          end do
       end do
    end if
  end subroutine OPRT_vorticity
  !-----------------------------------------------------------------------------
  subroutine OPRT_divdamp(       &
       grdx, grdx_pl,            &
       grdy, grdy_pl,            &
       grdz, grdz_pl,            &
       vx, vx_pl,                &
       vy, vy_pl,                &
       vz, vz_pl,                &
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
         GMTR_T_rarea,    &
         GMTR_A_hnx,      &
         GMTR_A_hny,      &
         GMTR_A_hnz,      &
         GMTR_A_tnx,      &
         GMTR_A_tny,      &
         GMTR_A_tnz,      &
         GMTR_A_tn2x,     &
         GMTR_A_tn2y,     &
         GMTR_A_tn2z,     &
         GMTR_P_rarea,    &
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    implicit none

    real(8), intent(inout) :: grdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: grdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: grdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
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
    real(8)  :: sclt(             &
         ADM_gall,                &
         ADM_kall,                &
         ADM_lall,                &
         ADM_TI:ADM_TJ )
    real(8)  :: sclt_pl(          &
         ADM_gall_pl,             &
         ADM_kall,                &
         ADM_lall_pl )
    real(8)  :: flux_pl(          &
         ADM_gall_pl)
    !
    integer :: l,n,k
    integer :: nstart, nend
    integer :: rgnid
    real(8) :: fact
    real(8) :: dp1,dp2,dp3,dp4,dp5,dp6
    !
    real(8) :: ux1,ux2,ux3
    real(8) :: uy1,uy2,uy3
    real(8) :: uz1,uz2,uz3
    real(8) :: f
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    if(present(mfact)) then
       fact=mfact
    else
       fact=1.0D0
    end if

    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k=ADM_kmin,ADM_kmax
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax,  ADM_gmax  )
          do n = nstart,nend
             !
             f=GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_RAREA)*0.5D0
             !
             ux1=-(vx(n,k,l)+vx(n+1,k,l))
             uy1=-(vy(n,k,l)+vy(n+1,k,l))
             uz1=-(vz(n,k,l)+vz(n+1,k,l))
             !
             ux2=-(vx(n+1,k,l)+vx(n+1+ADM_gall_1d,k,l))
             uy2=-(vy(n+1,k,l)+vy(n+1+ADM_gall_1d,k,l))
             uz2=-(vz(n+1,k,l)+vz(n+1+ADM_gall_1d,k,l))
             !
             ux3=+(vx(n+1+ADM_gall_1d,k,l)+vx(n,k,l))
             uy3=+(vy(n+1+ADM_gall_1d,k,l)+vy(n,k,l))
             uz3=+(vz(n+1+ADM_gall_1d,k,l)+vz(n,k,l))
             !
             sclt(n,k,l,ADM_TI)=f*(&
                  +ux1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_TNX)  &
                  +uy1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_TNY)  &
                  +uz1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_TNZ)  &
                  +ux2*GMTR_A_var(n+1,ADM_KNONE,l,ADM_AJ,GMTR_A_TNX)&
                  +uy2*GMTR_A_var(n+1,ADM_KNONE,l,ADM_AJ,GMTR_A_TNY)&
                  +uz2*GMTR_A_var(n+1,ADM_KNONE,l,ADM_AJ,GMTR_A_TNZ)&
                  +ux3*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_TNX) &
                  +uy3*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_TNY) &
                  +uz3*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_TNZ))
             !
             f=GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_RAREA)*0.5D0
             !
             ux1=-(vx(n,k,l)+vx(n+1+ADM_gall_1d,k,l))
             uy1=-(vy(n,k,l)+vy(n+1+ADM_gall_1d,k,l))
             uz1=-(vz(n,k,l)+vz(n+1+ADM_gall_1d,k,l))
             !
             ux2=+(vx(n+1+ADM_gall_1d,k,l)+vx(n+ADM_gall_1d,k,l))
             uy2=+(vy(n+1+ADM_gall_1d,k,l)+vy(n+ADM_gall_1d,k,l))
             uz2=+(vz(n+1+ADM_gall_1d,k,l)+vz(n+ADM_gall_1d,k,l))
             !
             ux3=+(vx(n+ADM_gall_1d,k,l)+vx(n,k,l))
             uy3=+(vy(n+ADM_gall_1d,k,l)+vy(n,k,l))
             uz3=+(vz(n+ADM_gall_1d,k,l)+vz(n,k,l))
             !
             sclt(n,k,l,ADM_TJ)=f*(&
                  +ux1*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_TNX)&
                  +uy1*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_TNY)&
                  +uz1*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_TNZ)&
                  +ux2*GMTR_A_var(n+ADM_gall_1d,ADM_KNONE,l,ADM_AI,GMTR_A_TNX)&
                  +uy2*GMTR_A_var(n+ADM_gall_1d,ADM_KNONE,l,ADM_AI,GMTR_A_TNY)&
                  +uz2*GMTR_A_var(n+ADM_gall_1d,ADM_KNONE,l,ADM_AI,GMTR_A_TNZ)&
                  +ux3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_TNX)&
                  +uy3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_TNY)&
                  +uz3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_TNZ))
          end do
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             sclt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =sclt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          do l=1,ADM_lall_pl
             do n=ADM_gmin_pl,ADM_GMAX_PL-1
                !
                f=GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)*0.5D0
                !
                ux1=+(vx_pl(ADM_gslf_pl,k,l)+vx_pl(n,k,l))
                uy1=+(vy_pl(ADM_gslf_pl,k,l)+vy_pl(n,k,l))
                uz1=+(vz_pl(ADM_gslf_pl,k,l)+vz_pl(n,k,l))
                !
                ux2=+(vx_pl(n,k,l)+vx_pl(n+1,k,l))
                uy2=+(vy_pl(n,k,l)+vy_pl(n+1,k,l))
                uz2=+(vz_pl(n,k,l)+vz_pl(n+1,k,l))
                !
                ux3=-(vx_pl(n+1,k,l)+vx_pl(ADM_gslf_pl,k,l))
                uy3=-(vy_pl(n+1,k,l)+vy_pl(ADM_gslf_pl,k,l))
                uz3=-(vz_pl(n+1,k,l)+vz_pl(ADM_gslf_pl,k,l))
                sclt_pl(n,k,l)=f*(&
                     +ux1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNX)&
                     +uy1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNY)&
                     +uz1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNZ)&
                     +ux2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2X)&
                     +uy2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2Y)&
                     +uz2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2Z)&
                     +ux3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNX)&
                     +uy3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNY)&
                     +uz3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNZ))
             end do
             !
             f=GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)*0.5D0
             !
             ux1=+(vx_pl(ADM_gslf_pl,k,l)+vx_pl(ADM_GMAX_PL,k,l))
             uy1=+(vy_pl(ADM_gslf_pl,k,l)+vy_pl(ADM_GMAX_PL,k,l))
             uz1=+(vz_pl(ADM_gslf_pl,k,l)+vz_pl(ADM_GMAX_PL,k,l))
             !
             ux2=+(vx_pl(ADM_GMAX_PL,k,l)+vx_pl(ADM_gmin_pl,k,l))
             uy2=+(vy_pl(ADM_GMAX_PL,k,l)+vy_pl(ADM_gmin_pl,k,l))
             uz2=+(vz_pl(ADM_GMAX_PL,k,l)+vz_pl(ADM_gmin_pl,k,l))
             !
             ux3=-(vx_pl(ADM_gmin_pl,k,l)+vx_pl(ADM_gslf_pl,k,l))
             uy3=-(vy_pl(ADM_gmin_pl,k,l)+vy_pl(ADM_gslf_pl,k,l))
             uz3=-(vz_pl(ADM_gmin_pl,k,l)+vz_pl(ADM_gslf_pl,k,l))
             !   
             sclt_pl(ADM_GMAX_PL,k,l)=f*(&
                  +ux1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNX)&
                  +uy1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNY)&
                  +uz1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNZ)&
                  +ux2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2X)&
                  +uy2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2Y)&
                  +uz2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2Z)&
                  +ux3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNX)&
                  +uy3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNY)&
                  +uz3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNZ))
             !
          end do
       end do
    end if

    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k=ADM_kmin,ADM_kmax
          nstart = suf(ADM_gmin  ,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             dp1=(sclt(n-ADM_gall_1d,   k,l,ADM_TJ)+sclt(n,               k,l,ADM_TI))*0.5D0
             dp2=(sclt(n,               k,l,ADM_TI)+sclt(n,               k,l,ADM_TJ))*0.5D0
             dp3=(sclt(n,               k,l,ADM_TJ)+sclt(n-1,             k,l,ADM_TI))*0.5D0
             dp4=(sclt(n-1-ADM_gall_1d, k,l,ADM_TJ)+sclt(n-1,             k,l,ADM_TI))*0.5D0
             dp5=(sclt(n-1-ADM_gall_1d, k,l,ADM_TI)+sclt(n-1-ADM_gall_1d, k,l,ADM_TJ))*0.5D0
             dp6=(sclt(n-ADM_gall_1d,   k,l,ADM_TJ)+sclt(n-1-ADM_gall_1d, k,l,ADM_TI))*0.5D0
             !
             grdx(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNX) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  -dp6*GMTR_A_var(n-ADM_gall_1d,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNX) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             grdy(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNY) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  -dp6*GMTR_A_var(n-ADM_gall_1d,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNY) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             grdz(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNZ) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  -dp6*GMTR_A_var(n-ADM_gall_1d,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNZ) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
          end do
          !
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             i=ADM_gmin
             j=ADM_gmin
             dp1=(sclt(suf(i,j)-ADM_gall_1d,     k,l,ADM_TJ)+sclt(suf(i,j),        k,l,ADM_TI))*0.5D0
             dp2=(sclt(suf(i,j),                 k,l,ADM_TI)+sclt(suf(i,j),        k,l,ADM_TJ))*0.5D0
             dp3=(sclt(suf(i,j),                 k,l,ADM_TJ)+sclt(suf(i,j)-1,      k,l,ADM_TI))*0.5D0
             dp4=(sclt(suf(i-1,j)-ADM_gall_1d,   k,l,ADM_TJ)+sclt(suf(i-1,j),      k,l,ADM_TI))*0.5D0
             dp5=(sclt(suf(i-1,j-1),             k,l,ADM_TI)+sclt(suf(i-1,j-1),    k,l,ADM_TJ))*0.5D0
             !
             grdx(suf(i,j),k,l)=(                                     &
                  +dp1*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  +dp2*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  +dp3*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AJ ,GMTR_A_HNX) &
                  -dp4*GMTR_A_var(suf(i-1,j),ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  -dp5*GMTR_A_var(suf(i-1,j-1),ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  ) * GMTR_P_var(suf(i,j),ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             grdy(suf(i,j),k,l)=(                                     &
                  +dp1*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  +dp2*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  +dp3*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AJ ,GMTR_A_HNY) &
                  -dp4*GMTR_A_var(suf(i-1,j),ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  -dp5*GMTR_A_var(suf(i-1,j-1),ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  ) * GMTR_P_var(suf(i,j),ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             grdz(suf(i,j),k,l)=(                                     &
                  +dp1*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  +dp2*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  +dp3*GMTR_A_var(suf(i,j),ADM_KNONE,l,ADM_AJ ,GMTR_A_HNZ) &
                  -dp4*GMTR_A_var(suf(i-1,j),ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  -dp5*GMTR_A_var(suf(i-1,j-1),ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  ) * GMTR_P_var(suf(i,j),ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
          end if
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          do l=1,ADM_lall_pl
             !
             flux_pl(ADM_gmin_pl)&
                  =(sclt_pl(ADM_GMAX_PL,k,l)+sclt_pl(ADM_gmin_pl,k,l))*0.5D0
             do n=ADM_gmin_pl+1,ADM_GMAX_PL
                flux_pl(n)&
                     =(sclt_pl(n-1,k,l)+sclt_pl(n,k,l))*0.5D0
             end do
             !
             dp1=flux_pl(ADM_gmin_pl  )
             dp2=flux_pl(ADM_gmin_pl+1)
             dp3=flux_pl(ADM_gmin_pl+2)
             dp4=flux_pl(ADM_gmin_pl+3)
             dp5=flux_pl(ADM_gmin_pl+4)
             !
             grdx_pl(ADM_gslf_pl,k,l)=(                                 &
                  +dp1*GMTR_A_var_pl(ADM_gmin_pl  ,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp2*GMTR_A_var_pl(ADM_gmin_pl+1,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp3*GMTR_A_var_pl(ADM_gmin_pl+2,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp4*GMTR_A_var_pl(ADM_gmin_pl+3,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp5*GMTR_A_var_pl(ADM_gmin_pl+4,ADM_KNONE,l,GMTR_A_HNX) &
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)  &
                  * fact
             !
             grdy_pl(ADM_gslf_pl,k,l)=(                                      &
                  +dp1*GMTR_A_var_pl(ADM_gmin_pl  ,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp2*GMTR_A_var_pl(ADM_gmin_pl+1,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp3*GMTR_A_var_pl(ADM_gmin_pl+2,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp4*GMTR_A_var_pl(ADM_gmin_pl+3,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp5*GMTR_A_var_pl(ADM_gmin_pl+4,ADM_KNONE,l,GMTR_A_HNY) &
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)  &
                  * fact
             !
             grdz_pl(ADM_gslf_pl,k,l)=(                                      &
                  +dp1*GMTR_A_var_pl(ADM_gmin_pl  ,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp2*GMTR_A_var_pl(ADM_gmin_pl+1,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp3*GMTR_A_var_pl(ADM_gmin_pl+2,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp4*GMTR_A_var_pl(ADM_gmin_pl+3,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp5*GMTR_A_var_pl(ADM_gmin_pl+4,ADM_KNONE,l,GMTR_A_HNZ) &
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)  &
                  * fact
             !
          end do
       end do
    end if

  end subroutine OPRT_divdamp

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient_nobase(        &
       vx, vx_pl,                  &
       vy, vy_pl,                  &
       vz, vz_pl,                  &
       scl, scl_pl,                &
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
         GMTR_A_hnx,      &
         GMTR_A_hny,      &
         GMTR_A_hnz,      &
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
    real(8), intent(in) :: scl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(inout)  :: vx(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout)  :: vy(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout)  :: vz(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in), optional :: mfact
    !
    real(8)  :: sclt(             &
         ADM_gall,                &
         ADM_kall,                &
         ADM_lall,                &
         ADM_TI:ADM_TJ )
    real(8)  :: sclt_pl(          &
         ADM_gall_pl,             &
         ADM_kall,                &
         ADM_lall_pl )
    real(8)  :: flux(             &
         ADM_gall,                &
         ADM_AI:ADM_AJ)
    real(8)  :: flux_pl(          &
         ADM_gall_pl)
    !
    integer :: l,n,k
    integer :: nstart,nend
    integer :: rgnid
    real(8) :: fact
    real(8) :: dp1,dp2,dp3,dp4,dp5,dp6
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    if(present(mfact)) then
       fact=mfact
    else
       fact=1.0D0
    end if

    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k=1, ADM_kall
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend = suf(ADM_gmax,ADM_gmax)
          do n = nstart,nend
             sclt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *scl(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *scl(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *scl(n+1+ADM_gall_1d,k,l)
             sclt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *scl(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *scl(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *scl(n+ADM_gall_1d,k,l)
          end do
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             sclt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =sclt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do k=1, ADM_kall
          do l=1,ADM_lall_pl
             do n=ADM_gmin_pl,ADM_GMAX_PL-1
                sclt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *scl_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *scl_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *scl_pl(n+1,k,l)
             end do
             sclt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *scl_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *scl_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *scl_pl(ADM_gmin_pl,k,l)
          end do
       end do
    end if
    !
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k=1, ADM_kall
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax,  ADM_gmax  )
          do n = nstart,nend
             flux(n,ADM_AI)&
                  =(sclt(n-ADM_gall_1d,k,l,ADM_TJ)+sclt(n,k,l,ADM_TI))*0.5D0
          end do
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,ADM_AIJ)&
                  =(sclt(n,k,l,ADM_TI)+sclt(n,k,l,ADM_TJ))*0.5D0
          end do
          !
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,ADM_AJ)&
                  =(sclt(n,k,l,ADM_TJ)+sclt(n-1,k,l,ADM_TI))*0.5D0
          end do
          !
          nstart = suf(ADM_gmin  ,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             dp1=flux(n,ADM_AI )
             dp2=flux(n,ADM_AIJ)
             dp3=flux(n,ADM_AJ )
             dp4=flux(n-1,ADM_AI )
             dp5=flux(n-1-ADM_gall_1d,ADM_AIJ)
             dp6=flux(n-ADM_gall_1d,ADM_AJ )
             !
             vx(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNX) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  -dp6*GMTR_A_var(n-ADM_gall_1d,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNX) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             vy(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNY) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  -dp6*GMTR_A_var(n-ADM_gall_1d,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNY) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             vz(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNZ) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  -dp6*GMTR_A_var(n-ADM_gall_1d,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNZ) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
          end do
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             n=suf(ADM_gmin,ADM_gmin)
             dp1=flux(n,ADM_AI )
             dp2=flux(n,ADM_AIJ)
             dp3=flux(n,ADM_AJ )
             dp4=flux(n-1,ADM_AI )
             dp5=flux(n-1-ADM_gall_1d,ADM_AIJ)
             !
             vx(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNX) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNX) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             vy(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNY) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNY) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
             !
             vz(n,k,l)=(                                     &
                  +dp1*GMTR_A_var(n,ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  +dp2*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  +dp3*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ ,GMTR_A_HNZ) &
                  -dp4*GMTR_A_var(n-1,ADM_KNONE,l,ADM_AI ,GMTR_A_HNZ) &
                  -dp5*GMTR_A_var(n-1-ADM_gall_1d,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ) &
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)      &
                  * fact
          end if
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do k=1, ADM_kall
          do l=1,ADM_lall_pl
             !
             flux_pl(ADM_gmin_pl)&
                  =(sclt_pl(ADM_GMAX_PL,k,l)+sclt_pl(ADM_gmin_pl,k,l))*0.5D0
             do n=ADM_gmin_pl+1,ADM_GMAX_PL
                flux_pl(n)&
                     =(sclt_pl(n-1,k,l)+sclt_pl(n,k,l))*0.5D0
             end do
             !
             dp1=flux_pl(ADM_gmin_pl  )
             dp2=flux_pl(ADM_gmin_pl+1)
             dp3=flux_pl(ADM_gmin_pl+2)
             dp4=flux_pl(ADM_gmin_pl+3)
             dp5=flux_pl(ADM_gmin_pl+4)
             !
             vx_pl(ADM_gslf_pl,k,l)=(                                 &
                  +dp1*GMTR_A_var_pl(ADM_gmin_pl  ,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp2*GMTR_A_var_pl(ADM_gmin_pl+1,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp3*GMTR_A_var_pl(ADM_gmin_pl+2,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp4*GMTR_A_var_pl(ADM_gmin_pl+3,ADM_KNONE,l,GMTR_A_HNX) &
                  +dp5*GMTR_A_var_pl(ADM_gmin_pl+4,ADM_KNONE,l,GMTR_A_HNX) &
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)  &
                  * fact
             !
             vy_pl(ADM_gslf_pl,k,l)=(                                      &
                  +dp1*GMTR_A_var_pl(ADM_gmin_pl  ,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp2*GMTR_A_var_pl(ADM_gmin_pl+1,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp3*GMTR_A_var_pl(ADM_gmin_pl+2,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp4*GMTR_A_var_pl(ADM_gmin_pl+3,ADM_KNONE,l,GMTR_A_HNY) &
                  +dp5*GMTR_A_var_pl(ADM_gmin_pl+4,ADM_KNONE,l,GMTR_A_HNY) &
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)  &
                  * fact
             !
             vz_pl(ADM_gslf_pl,k,l)=(                                      &
                  +dp1*GMTR_A_var_pl(ADM_gmin_pl  ,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp2*GMTR_A_var_pl(ADM_gmin_pl+1,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp3*GMTR_A_var_pl(ADM_gmin_pl+2,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp4*GMTR_A_var_pl(ADM_gmin_pl+3,ADM_KNONE,l,GMTR_A_HNZ) &
                  +dp5*GMTR_A_var_pl(ADM_gmin_pl+4,ADM_KNONE,l,GMTR_A_HNZ) &
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)  &
                  * fact
             !
          end do
       end do
    end if
  end subroutine OPRT_gradient_nobase

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence2_prep( &
       flx_h,   flx_h_pl,  &
       cp,      cp_pl,     &
       rhovx,  rhovx_pl, &
       rhovy,  rhovy_pl, &
       rhovz,  rhovz_pl, &
       rho,    rho_pl,   &
       dt                  &
       )
    !
    !--- Miura(2004)'s scheme : path1
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
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         ADM_GMAX_PL,     &
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
         ADM_gall,        &
         ADM_ImoJmo_nmax, &  ! Y.Niwa add 080130
         ADM_ImoJmo,      &  ! Y.Niwa add 080130
         ADM_GIoJo,       &  ! Y.Niwa add 080130
         ADM_VMISS           ! Y.Niwa add 080130
    use mod_gmtr, only :  &
         !--- public parameters
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_T_rarea,    &
         GMTR_A_hnx,      &
         GMTR_A_hny,      &
         GMTR_A_hnz,      &
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
    use mod_grd, only :    &
         GRD_x, GRD_x_pl,  &
         GRD_xt, GRD_xt_pl,&
         GRD_XDIR,         &  
         GRD_YDIR,         &
         GRD_ZDIR
    use mod_comm, only : &
         COMM_data_transfer
    use mod_cnst, only : &
         CNST_MAX_REAL,  &
         CNST_EPS_ZERO
    implicit none

    real(8), intent(out) :: flx_h(6,ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(out) :: flx_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(out) :: cp(ADM_gall,ADM_kall,ADM_lall,ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8), intent(out) :: cp_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,GRD_XDIR:GRD_ZDIR)
    real(8), intent(in)  :: rhovx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rhovx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: rhovy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rhovy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: rhovz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rhovz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: rho(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: rho_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in)  :: dt

    real(8)  :: rhovxt(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8)  :: rhovxt_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8)  :: rhovyt(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8)  :: rhovyt_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8)  :: rhovzt(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8)  :: rhovzt_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8)  :: rhot(ADM_gall,ADM_kall,ADM_lall,ADM_TI:ADM_TJ)
    real(8) ::  rhot_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    real(8)  :: ccc
    real(8)  :: ccc_pl
    !
    real(8)  :: rpx(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rpx_pl(ADM_gall_pl)
    real(8)  :: rpy(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rpy_pl(ADM_gall_pl)
    real(8)  :: rpz(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rpz_pl(ADM_gall_pl)
    !
    real(8)  :: rhovxa(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhovxa_pl(ADM_GALL_PL)
    real(8)  :: rhovya(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhovya_pl(ADM_GALL_PL)
    real(8)  :: rhovza(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhovza_pl(ADM_GALL_PL)
    real(8)  :: rhoa(ADM_gall,ADM_AI:ADM_AJ)
    real(8)  :: rhoa_pl(ADM_GALL_PL)
    !
    integer :: l,n,k
    integer :: rgnid
    !
    integer :: np1(ADM_gall_pl)
    integer :: nm1(ADM_gall_pl)
    !
    integer :: nstart,nend
    logical, save :: first = .true.  ! Y.Niwa add 080130
    integer :: t  ! Y.Niwa add 080130
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)





    ! Y.Niwa add 080130 =>
    if(first) then
      allocate(local_t_var(ADM_gall,ADM_KNONE,ADM_lall,ADM_TI:ADM_TJ,GMTR_T_W1:GMTR_T_W3))
      allocate(local_t_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,GMTR_T_W1:GMTR_T_W3))

      do l=1, ADM_lall
         do t=ADM_TI, ADM_TJ
            do n=1, ADM_ImoJmo_nmax
               local_t_var(ADM_ImoJmo(n,ADM_GIoJo),ADM_KNONE,l,t,GMTR_T_W1) &
                  =  GMTR_t_var(ADM_ImoJmo(n,ADM_GIoJo),ADM_KNONE,l,t,GMTR_T_W1)
               local_t_var(ADM_ImoJmo(n,ADM_GIoJo),ADM_KNONE,l,t,GMTR_T_W2) &
                  =  GMTR_t_var(ADM_ImoJmo(n,ADM_GIoJo),ADM_KNONE,l,t,GMTR_T_W2)
               local_t_var(ADM_ImoJmo(n,ADM_GIoJo),ADM_KNONE,l,t,GMTR_T_W3) &
                  =  GMTR_t_var(ADM_ImoJmo(n,ADM_GIoJo),ADM_KNONE,l,t,GMTR_T_W3)
            end do
         end do
      end do
      !
      if(ADM_prc_me==ADM_prc_pl) then
         do l=1, ADM_lall_pl
            do n=ADM_gmin_pl, ADM_GMAX_PL
               local_t_var_pl(n,ADM_KNONE,l,GMTR_T_W1) &
                  = GMTR_t_var_pl(n,ADM_KNONE,l,GMTR_T_W1)
               local_t_var_pl(n,ADM_KNONE,l,GMTR_T_W2) &
                  = GMTR_t_var_pl(n,ADM_KNONE,l,GMTR_T_W2)
               local_t_var_pl(n,ADM_KNONE,l,GMTR_T_W3) &
                  = GMTR_t_var_pl(n,ADM_KNONE,l,GMTR_T_W3)
            end do
         end do
      end if
      first = .false.
    end if
    !<= Y.Niwa add 080130
    ! 
    flx_h(:,:,:,:) = 0.0d0
    flx_h_pl(:,:,:) = 0.0d0
    cp(:,:,:,:,:) = 0.0d0
    cp_pl(:,:,:,:) = 0.0d0
    !
    nm1(ADM_gmin_pl) = ADM_GMAX_PL
    do n=ADM_gmin_pl+1,ADM_GMAX_PL
       nm1(n) = n-1
    end do
    !
    do n=ADM_gmin_pl,ADM_GMAX_PL-1
       np1(n) = n+1
    end do
    np1(ADM_GMAX_PL) = ADM_gmin_pl
    !
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend = suf(ADM_gmax,ADM_gmax)
          do n = nstart,nend
             rhovxt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *rhovx(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *rhovx(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *rhovx(n+1+ADM_gall_1d,k,l)
             rhovxt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *rhovx(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *rhovx(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *rhovx(n+ADM_gall_1d,k,l)
             rhovyt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *rhovy(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *rhovy(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *rhovy(n+1+ADM_gall_1d,k,l)
             rhovyt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *rhovy(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *rhovy(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *rhovy(n+ADM_gall_1d,k,l)
             rhovzt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *rhovz(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *rhovz(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *rhovz(n+1+ADM_gall_1d,k,l)
             rhovzt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *rhovz(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *rhovz(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *rhovz(n+ADM_gall_1d,k,l)
             rhot(n,k,l,ADM_TI)                            &
                  =local_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *rho(n,k,l)                         &
                  +local_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *rho(n+1,k,l)                       &
                  +local_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *rho(n+1+ADM_gall_1d,k,l)
!                 =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&  ! A.Noda 110728
!                 *rho(n,k,l)                         &
!                 +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
!                 *rho(n+1,k,l)                       &
!                 +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
!                 *rho(n+1+ADM_gall_1d,k,l)
             rhot(n,k,l,ADM_TJ)                            &
                 =local_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *rho(n,k,l)                         &
                  +local_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *rho(n+1+ADM_gall_1d,k,l)                     &
                  +local_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *rho(n+ADM_gall_1d,k,l)
!!$                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&   ! Y.Niwa 080130
!!$                  *rho(n,k,l)                         &
!!$                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
!!$                  *rho(n+1+ADM_gall_1d,k,l)                     &
!!$                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
!!$                  *rho(n+ADM_gall_1d,k,l)
          end do
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             rhovxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhovxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             rhovyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhovyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             rhovzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhovzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             rhot(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =rhot(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
          !
       end do
       !
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_GMAX_PL
                rhovxt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *rhovx_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *rhovx_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *rhovx_pl(np1(n),k,l)
                rhovyt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *rhovy_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *rhovy_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *rhovy_pl(np1(n),k,l)
                rhovzt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *rhovz_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *rhovz_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *rhovz_pl(np1(n),k,l)
                rhot_pl(n,k,l)                            &
                     =local_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *rho_pl(ADM_gslf_pl,k,l)              &
                     +local_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *rho_pl(n,k,l)                        &
                     +local_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *rho_pl(np1(n),k,l)
!!$                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&  ! Y.Niwa 080130
!!$                     *rho_pl(ADM_gslf_pl,k,l)              &
!!$                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
!!$                     *rho_pl(n,k,l)                        &
!!$                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
!!$                     *rho_pl(np1(n),k,l)
             end do
          end do
       end do
    end if
    !
    !--- calculation of courant number
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k = 1, ADM_kall
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax,  ADM_gmax  )
          do n = nstart,nend
             rpx(n,ADM_AI) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TI,GRD_XDIR)+GRD_xt(n-ADM_gall_1d,ADM_KNONE,l,ADM_TJ,GRD_XDIR))
             rpy(n,ADM_AI) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TI,GRD_YDIR)+GRD_xt(n-ADM_gall_1d,ADM_KNONE,l,ADM_TJ,GRD_YDIR))
             rpz(n,ADM_AI) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TI,GRD_ZDIR)+GRD_xt(n-ADM_gall_1d,ADM_KNONE,l,ADM_TJ,GRD_ZDIR))
             rhovxa(n,ADM_AI) = (rhovxt(n-ADM_gall_1d,k,l,ADM_TJ)+rhovxt(n,k,l,ADM_TI))*0.5D0
             rhovya(n,ADM_AI) = (rhovyt(n-ADM_gall_1d,k,l,ADM_TJ)+rhovyt(n,k,l,ADM_TI))*0.5D0
             rhovza(n,ADM_AI) = (rhovzt(n-ADM_gall_1d,k,l,ADM_TJ)+rhovzt(n,k,l,ADM_TI))*0.5D0
             rhoa(n,ADM_AI) = (rhot(n-ADM_gall_1d,k,l,ADM_TJ)+rhot(n,k,l,ADM_TI))*0.5D0
             !
             cp(n,k,l,ADM_AI,GRD_XDIR) = rpx(n,ADM_AI) - rhovxa(n,ADM_AI)/rhoa(n,ADM_AI)*dt*0.5D0
             cp(n,k,l,ADM_AI,GRD_YDIR) = rpy(n,ADM_AI) - rhovya(n,ADM_AI)/rhoa(n,ADM_AI)*dt*0.5D0
             cp(n,k,l,ADM_AI,GRD_ZDIR) = rpz(n,ADM_AI) - rhovza(n,ADM_AI)/rhoa(n,ADM_AI)*dt*0.5D0
             !
             ccc = &
                  (rhovxa(n,ADM_AI)*GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HNX)&
                  +rhovya(n,ADM_AI)*GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HNY)&
                  +rhovza(n,ADM_AI)*GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HNZ))
             flx_h(1,n  ,k,l) = ccc*dt*GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)
             flx_h(4,n+1,k,l) =-ccc*dt*GMTR_P_var(n+1,ADM_KNONE,l,GMTR_P_RAREA)
             !
          end do
          !
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax,  ADM_gmax  )
          do n = nstart,nend
             rpx(n,ADM_AIJ) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TI,GRD_XDIR)+GRD_xt(n,ADM_KNONE,l,ADM_TJ,GRD_XDIR))
             rpy(n,ADM_AIJ) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TI,GRD_YDIR)+GRD_xt(n,ADM_KNONE,l,ADM_TJ,GRD_YDIR))
             rpz(n,ADM_AIJ) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TI,GRD_ZDIR)+GRD_xt(n,ADM_KNONE,l,ADM_TJ,GRD_ZDIR))
             rhovxa(n,ADM_AIJ) = (rhovxt(n,k,l,ADM_TI)+rhovxt(n,k,l,ADM_TJ))*0.5D0
             rhovya(n,ADM_AIJ) = (rhovyt(n,k,l,ADM_TI)+rhovyt(n,k,l,ADM_TJ))*0.5D0
             rhovza(n,ADM_AIJ) = (rhovzt(n,k,l,ADM_TI)+rhovzt(n,k,l,ADM_TJ))*0.5D0
             rhoa(n,ADM_AIJ) = (rhot(n,k,l,ADM_TI)+rhot(n,k,l,ADM_TJ))*0.5D0
             !
             cp(n,k,l,ADM_AIJ,GRD_XDIR) = rpx(n,ADM_AIJ) - rhovxa(n,ADM_AIJ)/rhoa(n,ADM_AIJ)*dt*0.5D0
             cp(n,k,l,ADM_AIJ,GRD_YDIR) = rpy(n,ADM_AIJ) - rhovya(n,ADM_AIJ)/rhoa(n,ADM_AIJ)*dt*0.5D0
             cp(n,k,l,ADM_AIJ,GRD_ZDIR) = rpz(n,ADM_AIJ) - rhovza(n,ADM_AIJ)/rhoa(n,ADM_AIJ)*dt*0.5D0
             !
             ccc = &
                  (rhovxa(n,ADM_AIJ)*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNX)&
                  +rhovya(n,ADM_AIJ)*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNY)&
                  +rhovza(n,ADM_AIJ)*GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HNZ))
             flx_h(2,n              ,k,l) = ccc*dt*GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)
             flx_h(5,n+1+ADM_gall_1d,k,l) =-ccc*dt*GMTR_P_var(n+1+ADM_gall_1d,ADM_KNONE,l,GMTR_P_RAREA)
             !
          end do
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             rpx(n,ADM_AJ) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TJ,GRD_XDIR)+GRD_xt(n-1,ADM_KNONE,l,ADM_TI,GRD_XDIR))
             rpy(n,ADM_AJ) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TJ,GRD_YDIR)+GRD_xt(n-1,ADM_KNONE,l,ADM_TI,GRD_YDIR))
             rpz(n,ADM_AJ) &
                  = 0.5D0*(GRD_xt(n,ADM_KNONE,l,ADM_TJ,GRD_ZDIR)+GRD_xt(n-1,ADM_KNONE,l,ADM_TI,GRD_ZDIR))
             rhovxa(n,ADM_AJ)=  (rhovxt(n,k,l,ADM_TJ)+rhovxt(n-1,k,l,ADM_TI))*0.5D0
             rhovya(n,ADM_AJ)=  (rhovyt(n,k,l,ADM_TJ)+rhovyt(n-1,k,l,ADM_TI))*0.5D0
             rhovza(n,ADM_AJ)=  (rhovzt(n,k,l,ADM_TJ)+rhovzt(n-1,k,l,ADM_TI))*0.5D0
             rhoa(n,ADM_AJ)=  (rhot(n,k,l,ADM_TJ)+rhot(n-1,k,l,ADM_TI))*0.5D0
             !
             cp(n,k,l,ADM_AJ,GRD_XDIR) = rpx(n,ADM_AJ) - rhovxa(n,ADM_AJ)/rhoa(n,ADM_AJ)*dt*0.5D0
             cp(n,k,l,ADM_AJ,GRD_YDIR) = rpy(n,ADM_AJ) - rhovya(n,ADM_AJ)/rhoa(n,ADM_AJ)*dt*0.5D0
             cp(n,k,l,ADM_AJ,GRD_ZDIR) = rpz(n,ADM_AJ) - rhovza(n,ADM_AJ)/rhoa(n,ADM_AJ)*dt*0.5D0
             !
             ccc = &
                  (rhovxa(n,ADM_AJ)*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HNX)&
                  +rhovya(n,ADM_AJ)*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HNY)&
                  +rhovza(n,ADM_AJ)*GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HNZ))
             flx_h(3,n            ,k,l) = ccc*dt*GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)
             flx_h(6,n+ADM_gall_1d,k,l) =-ccc*dt*GMTR_P_var(n+ADM_gall_1d,ADM_KNONE,l,GMTR_P_RAREA)
          end do
          !
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             flx_h(6,suf(ADM_gmin,ADM_gmin),k,l) = 0.0D0
          end if
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_GMAX_PL
                rpx_pl(n) &
                     = 0.5D0*(GRD_xt_pl(nm1(n),ADM_KNONE,l,GRD_XDIR)+GRD_xt_pl(n,ADM_KNONE,l,GRD_XDIR))
                rpy_pl(n) &
                     = 0.5D0*(GRD_xt_pl(nm1(n),ADM_KNONE,l,GRD_YDIR)+GRD_xt_pl(n,ADM_KNONE,l,GRD_YDIR))
                rpz_pl(n) &
                     = 0.5D0*(GRD_xt_pl(nm1(n),ADM_KNONE,l,GRD_ZDIR)+GRD_xt_pl(n,ADM_KNONE,l,GRD_ZDIR))
                !
                rhovxa_pl(n)=  (rhovxt_pl(nm1(n),k,l)+rhovxt_pl(n,k,l))*0.5D0
                rhovya_pl(n)=  (rhovyt_pl(nm1(n),k,l)+rhovyt_pl(n,k,l))*0.5D0
                rhovza_pl(n)=  (rhovzt_pl(nm1(n),k,l)+rhovzt_pl(n,k,l))*0.5D0
                rhoa_pl(n)=  (rhot_pl(nm1(n),k,l)+rhot_pl(n,k,l))*0.5D0
                !
                cp_pl(n,k,l,GRD_XDIR) = rpx_pl(n) - rhovxa_pl(n)/rhoa_pl(n)*dt*0.5D0
                cp_pl(n,k,l,GRD_YDIR) = rpy_pl(n) - rhovya_pl(n)/rhoa_pl(n)*dt*0.5D0
                cp_pl(n,k,l,GRD_ZDIR) = rpz_pl(n) - rhovza_pl(n)/rhoa_pl(n)*dt*0.5D0
                !
                ccc_pl = &
                     (rhovxa_pl(n)*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNX)&
                     +rhovya_pl(n)*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNY)&
                     +rhovza_pl(n)*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNZ))
                flx_h_pl(n,k,l) = ccc_pl*dt*GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)
             end do
          end do
       end do
    end if

    return
  end subroutine OPRT_divergence2_prep

  !-----------------------------------------------------------------------------
  ! Miura(2004)'s scheme with Thuburn(1996) limiter : path2 
  subroutine OPRT_divergence2( &
       scl,   scl_pl,   &
       s,     s_pl,     &
       flx_h, flx_h_pl, &
       c,     c_pl,     &
       cp,    cp_pl,    &
       d,     d_pl,     &
       dt,              &
       limiter,         &
       mfact            )
    use mod_adm, only: &
       ADM_prc_pl,      &
       ADM_AI,          &
       ADM_AIJ,         &
       ADM_AJ,          &
       ADM_KNONE,       &
       ADM_lall_pl,     &
       ADM_gmin_pl,     &
       ADM_gslf_pl,     &
       ADM_gall_pl,     &
       ADM_GMAX_PL,     &
       ADM_prc_me,      &
       ADM_prc_tab,     &
       ADM_gall_1d,     &
       ADM_kall,        &
       ADM_lall,        &
       ADM_kmin,        &
       ADM_kmax,        &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_gall,        &
       ADM_W,           &
       ADM_rgn_vnum
    use mod_grd, only: &
       GRD_XDIR, &  
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_x,    &
       GRD_x_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_cnst, only: &
       CNST_EPS_ZERO, & 
       CNST_MAX_REAL
    implicit none

    real(8), intent(out) :: scl     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: scl_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: s       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: s_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: flx_h   (6,ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: flx_h_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: c       (6,ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in)  :: c_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: cp      (ADM_gall,ADM_kall,ADM_lall,ADM_AI:ADM_AJ,GRD_XDIR:GRD_ZDIR)
    real(8), intent(in)  :: cp_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR)
    real(8), intent(in)  :: d       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: d_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: dt

    character(len=*), intent(in), optional :: limiter
    real(8),          intent(in), optional :: mfact

    real(8)  :: sa(ADM_AI:ADM_AJ,ADM_gall,ADM_kall,ADM_lall)
    real(8)  :: sa_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8)  :: sa_p,sa_m

    real(8)  :: sw(ADM_AI:ADM_AJ,ADM_gall)

    real(8)  :: s_in_min(6,ADM_gall,ADM_kall,ADM_lall)
    real(8)  :: s_in_min_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(8)  :: s_in_max(6,ADM_gall,ADM_kall,ADM_lall)
    real(8)  :: s_in_max_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    !
    real(8)  :: s_m1_k_min(ADM_gall)
    real(8)  :: s_m1_k_min_pl
    real(8)  :: s_m1_k_max(ADM_gall)
    real(8)  :: s_m1_k_max_pl
    real(8)  :: c_in_sum(ADM_gall)
    real(8)  :: c_in_sum_pl
    real(8)  :: c_out_sum(ADM_gall)
    real(8)  :: c_out_sum_pl
    real(8)  :: c_qin_sum_max(ADM_gall)
    real(8)  :: c_qin_sum_max_pl
    real(8)  :: c_qin_sum_min(ADM_gall)
    real(8)  :: c_qin_sum_min_pl
    !
    real(8)  :: s_m1_k_min_n
    real(8)  :: s_m1_k_max_n
    real(8)  :: c_in_sum_n
    real(8)  :: c_out_sum_n
    real(8)  :: c_qin_sum_max_n
    real(8)  :: c_qin_sum_min_n
    !
    integer :: np1(ADM_gall_pl)
    integer :: nm1(ADM_gall_pl)
    !
    real(8) :: SMALL_ZERO = 0.0d0
    !
    integer, parameter :: dsx=1
    integer, parameter :: dsy=2
    integer, parameter :: dsz=3
    integer, parameter :: s_out_k_min=4
    integer, parameter :: s_out_k_max=5
    real(8)  :: wrk(ADM_gall,ADM_kall,ADM_lall,dsx:s_out_k_max)
    real(8)  :: wrk_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,dsx:s_out_k_max)
    !
!    logical :: NON_NEG=.false.
    !
    integer :: l,n,k,m
    integer :: rgnid
    real(8) :: fact
    !
    integer :: nstart,nend
    integer :: nstart2,nstart3

    real(8) :: AI_min, AIJ_min, AJ_min
    real(8) :: AI_max, AIJ_max, AJ_max

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1, ip2jp1
    integer :: im1j, ijm1

    integer :: ierr

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( present(mfact) ) then
       fact = mfact / dt
    else
       fact = 1.D0 / dt
    endif

    nm1(ADM_gmin_pl) = ADM_GMAX_PL
    do n = ADM_gmin_pl+1, ADM_GMAX_PL
       nm1(n) = n-1
    enddo

    do n = ADM_gmin_pl, ADM_GMAX_PL-1
       np1(n) = n+1
    enddo
    np1(ADM_GMAX_PL) = ADM_gmin_pl

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       scl(n,k,l) = 0.D0
    enddo
    enddo
    enddo

    do l = 1, ADM_lall_pl
    do k = 1, ADM_kall
    do n = 1, ADM_gall_pl
       scl_pl(n,k,l) = 0.D0
    enddo
    enddo
    enddo

    ! [Fix] 08/04/28 T.Mitsui, to avoid undefined reference in outflow limiter
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do n = 1, ADM_gall
             wrk(n,k,l,s_out_k_min) = s(n,k,l)
          enddo
          do n = 1, ADM_gall
             wrk(n,k,l,s_out_k_max) = s(n,k,l)
          enddo
       enddo
    enddo

    do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          do n = 1,ADM_gall_pl
             wrk_pl(n,k,l,s_out_k_min) = s_pl(n,k,l)
          enddo
          do n = 1,ADM_gall_pl
             wrk_pl(n,k,l,s_out_k_max) = s_pl(n,k,l)
          enddo
       enddo
    enddo

    if(present(limiter)) then
       if(trim(limiter)=='NON_LIM') go to 1000
    endif

    !--- calculation of inflow limiter
!    do l = 1, ADM_lall
!    do k = 1, ADM_kall
!    do n = 1, ADM_gall
!    do m = 1, 6
!       s_in_min(m,n,k,l) = CNST_MAX_REAL
!    enddo
!    enddo
!    enddo
!    enddo
!
!    do l = 1, ADM_lall
!    do k = 1, ADM_kall
!    do n = 1, ADM_gall
!    do m = 1, 6
!       s_in_max(m,n,k,l) = -CNST_MAX_REAL
!    enddo
!    enddo
!    enddo
!    enddo

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do k = 1, ADM_kall

          do n = 1, ADM_gall
             do m = ADM_AI, ADM_AJ
                sw(m,n) = 0.D0
             enddo
          enddo

          do n = 1, ADM_gall
             if ( c(1,n,k,l) <= 0.D0 ) sw(ADM_AI, n) = 1.D0
             if ( c(2,n,k,l) <= 0.D0 ) sw(ADM_AIJ,n) = 1.D0
             if ( c(3,n,k,l) <= 0.D0 ) sw(ADM_AJ, n) = 1.D0
          enddo

          nstart = suf(ADM_gmin,ADM_gmin)
          nend   = suf(ADM_gmax,ADM_gmax)

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ijp1   = n  +ADM_gall_1d
             ip1jp1 = n+1+ADM_gall_1d
             im1j   = n-1
             ijm1   = n  -ADM_gall_1d

             AI_min  = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
             AI_max  = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
             AIJ_min = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
             AIJ_max = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
             AJ_min  = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
             AJ_max  = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )

!             if ( c(1,n,k,l) <= 0.D0 ) then
!                s_in_min(1,ij,  k,l) = AI_min * 1.D0 + 0.D0
!                s_in_max(1,ij,  k,l) = AI_max * 1.D0 + 0.D0
!                s_in_min(4,ip1j,k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(4,ip1j,k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!             else
!                s_in_min(1,ij,  k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(1,ij,  k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_min(4,ip1j,k,l) = AI_min * 1.D0 + 0.D0
!                s_in_max(4,ip1j,k,l) = AI_max * 1.D0 + 0.D0
!             endif
!             if ( c(2,n,k,l) <= 0.D0 ) then
!                s_in_min(2,ij,    k,l) = AIJ_min * 1.D0 + 0.D0
!                s_in_max(2,ij,    k,l) = AIJ_max * 1.D0 + 0.D0
!                s_in_min(5,ip1jp1,k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(5,ip1jp1,k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!             else
!                s_in_min(2,ij,    k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(2,ij,    k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_min(5,ip1jp1,k,l) = AIJ_min * 1.D0 + 0.D0
!                s_in_max(5,ip1jp1,k,l) = AIJ_max * 1.D0 + 0.D0
!             endif
!             if ( c(3,n,k,l) <= 0.D0 ) then
!                s_in_min(3,ij,  k,l) = AJ_min * 1.D0 + 0.D0
!                s_in_max(3,ij,  k,l) = AJ_max * 1.D0 + 0.D0
!                s_in_min(6,ijp1,k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(6,ijp1,k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!             else
!                s_in_min(3,ij,  k,l) = CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_max(3,ij,  k,l) = -CNST_MAX_REAL * 1.D0 + 0.D0
!                s_in_min(6,ijp1,k,l) = AJ_min * 1.D0 + 0.D0
!                s_in_max(6,ijp1,k,l) = AJ_max * 1.D0 + 0.D0
!             endif

             s_in_min(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_min        &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * CNST_MAX_REAL
             s_in_min(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_min       &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * CNST_MAX_REAL
             s_in_min(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_min        &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * CNST_MAX_REAL
             s_in_min(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_min
             s_in_min(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_min
             s_in_min(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_min

             s_in_max(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_max         &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * (-CNST_MAX_REAL)
             s_in_max(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_max        &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL)
             s_in_max(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_max         &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * (-CNST_MAX_REAL)
             s_in_max(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_max
             s_in_max(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_max
             s_in_max(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_max
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmin-1,ADM_gmin  )

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ip1jp1 = n+1+ADM_gall_1d
             ijm1   = n  -ADM_gall_1d

             AI_min  = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
             AI_max  = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )

             s_in_min(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_min        &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * CNST_MAX_REAL
             s_in_min(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_min

             s_in_max(1,ij,    k,l) = (        sw(ADM_AI, n) ) * AI_max         &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * (-CNST_MAX_REAL)
             s_in_max(4,ip1j,  k,l) = (        sw(ADM_AI, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AI, n) ) * AI_max

!             if ( c(1,n,k,l) <= 0.D0 ) then
!                s_in_min(1,ij,  k,l) = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!                s_in_max(1,ij,  k,l) = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!             else
!                s_in_min(4,ip1j,k,l) = min( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!                s_in_max(4,ip1j,k,l) = max( s(ij,k,l), s(ip1j,k,l), s(ip1jp1,k,l), s(ijm1,k,l) )
!             endif
          enddo

          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmin-1,ADM_gmin  )

          do n = nstart, nend
             ij     = n
             ip1j   = n+1
             ip1jp1 = n+1+ADM_gall_1d
             ijp1   = n  +ADM_gall_1d

             AIJ_min = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
             AIJ_max = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )

             s_in_min(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_min       &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * CNST_MAX_REAL
             s_in_min(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_min
             s_in_max(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_max        &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL)
             s_in_max(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_max

!             if ( c(2,n,k,l) <= 0.D0 ) then
!                s_in_min(2,ij,    k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!                s_in_max(2,ij,    k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!             else
!                s_in_min(5,ip1jp1,k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!                s_in_max(5,ip1jp1,k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip1j,k,l), s(ijp1,k,l) )
!             endif
          enddo

          nstart = suf(ADM_gmin,  ADM_gmin-1)
          nend   = suf(ADM_gmin-1,ADM_gmin  )

          do n = nstart, nend
             ij     = n
             ip1jp1 = n+1+ADM_gall_1d
             ijp1   = n  +ADM_gall_1d
             im1j   = n-1

             AJ_min  = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
             AJ_max  = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )

             s_in_min(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_min        &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * CNST_MAX_REAL
             s_in_min(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_min
             s_in_max(3,ij,    k,l) = (        sw(ADM_AJ, n) ) * AJ_max         &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * (-CNST_MAX_REAL)
             s_in_max(6,ijp1,  k,l) = (        sw(ADM_AJ, n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AJ, n) ) * AJ_max

!             if ( c(3,n,k,l) <= 0.D0 ) then
!                s_in_min(3,ij,  k,l) = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!                s_in_max(3,ij,  k,l) = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!             else
!                s_in_min(6,ijp1,k,l) = min( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!                s_in_max(6,ijp1,k,l) = max( s(ij,k,l), s(ijp1,k,l), s(ip1jp1,k,l), s(im1j,k,l) )
!             endif
          enddo

          if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
             n = suf(ADM_gmin-1,ADM_gmin-1)
             ij     = n
             ip1jp1 = n+1+ADM_gall_1d
             ip2jp1 = n+2+ADM_gall_1d
             ijp1   = n  +ADM_gall_1d

             AIJ_min = min( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
             AIJ_max = max( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )

             s_in_min(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_min       &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * CNST_MAX_REAL
             s_in_min(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * CNST_MAX_REAL &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_min
             s_in_max(2,ij,    k,l) = (        sw(ADM_AIJ,n) ) * AIJ_max        &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL)
             s_in_max(5,ip1jp1,k,l) = (        sw(ADM_AIJ,n) ) * (-CNST_MAX_REAL) &
                                    + ( 1.D0 - sw(ADM_AIJ,n) ) * AIJ_max

!             if ( c(2,n,k,l) <= 0.D0 ) then
!                s_in_min(2,ij,    k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!                s_in_max(2,ij,    k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!             else
!                s_in_min(5,ip1jp1,k,l) = min( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!                s_in_max(5,ip1jp1,k,l) = max( s(ij,k,l), s(ip1jp1,k,l), s(ip2jp1,k,l), s(ijp1,k,l) )
!             endif
          endif

       enddo
    enddo

    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             s_in_min_pl(:,k,l,:) = CNST_MAX_REAL
             s_in_max_pl(:,k,l,:) =-CNST_MAX_REAL
             do n=ADM_gmin_pl,ADM_GMAX_PL
                if(c_pl(n,k,l)<=0.0D0) then
                   s_in_min_pl(n,k,l,1) = min(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                   s_in_max_pl(n,k,l,1) = max(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                else
                   s_in_min_pl(n,k,l,2) = min(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                   s_in_max_pl(n,k,l,2) = max(s_pl(ADM_gslf_pl,k,l),s_pl(n,k,l),s_pl(nm1(n),k,l),s_pl(np1(n),k,l))
                end if
             end do

          end do
          !
       end do
       !
    end if
    !
    !--- calcluation outflow limiter
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
       nstart = suf(ADM_gmin, ADM_gmin )
       nend  = suf(ADM_gmax ,ADM_gmax )
       do n = nstart,nend
          s_m1_k_min_n = min(s_in_min(1,n,k,l),s_in_min(2,n,k,l),s_in_min(3,n,k,l),&
                             s_in_min(4,n,k,l),s_in_min(5,n,k,l),s_in_min(6,n,k,l))
          if(s_m1_k_min_n==CNST_MAX_REAL) s_m1_k_min_n = s(n,k,l)
          s_m1_k_max_n = max(s_in_max(1,n,k,l),s_in_max(2,n,k,l),s_in_max(3,n,k,l),&
                             s_in_max(4,n,k,l),s_in_max(5,n,k,l),s_in_max(6,n,k,l))
          if(s_m1_k_max_n==-CNST_MAX_REAL) s_m1_k_max_n = s(n,k,l)

          c_in_sum_n &
            = (0.5D0-sign(0.5D0,c(1,n,k,l)))*c(1,n,k,l)&
            + (0.5D0-sign(0.5D0,c(2,n,k,l)))*c(2,n,k,l)&
            + (0.5D0-sign(0.5D0,c(3,n,k,l)))*c(3,n,k,l)&
            + (0.5D0-sign(0.5D0,c(4,n,k,l)))*c(4,n,k,l)&
            + (0.5D0-sign(0.5D0,c(5,n,k,l)))*c(5,n,k,l)&
            + (0.5D0-sign(0.5D0,c(6,n,k,l)))*c(6,n,k,l)

          c_out_sum_n &
            = (0.5D0+sign(0.5D0,c(1,n,k,l)))*c(1,n,k,l)&
            + (0.5D0+sign(0.5D0,c(2,n,k,l)))*c(2,n,k,l)&
            + (0.5D0+sign(0.5D0,c(3,n,k,l)))*c(3,n,k,l)&
            + (0.5D0+sign(0.5D0,c(4,n,k,l)))*c(4,n,k,l)&
            + (0.5D0+sign(0.5D0,c(5,n,k,l)))*c(5,n,k,l)&
            + (0.5D0+sign(0.5D0,c(6,n,k,l)))*c(6,n,k,l)
          
          c_qin_sum_max_n &
            = (0.5D0-sign(0.5D0,c(1,n,k,l)))*(c(1,n,k,l)*s_in_max(1,n,k,l))&
            + (0.5D0-sign(0.5D0,c(2,n,k,l)))*(c(2,n,k,l)*s_in_max(2,n,k,l))&
            + (0.5D0-sign(0.5D0,c(3,n,k,l)))*(c(3,n,k,l)*s_in_max(3,n,k,l))&
            + (0.5D0-sign(0.5D0,c(4,n,k,l)))*(c(4,n,k,l)*s_in_max(4,n,k,l))&
            + (0.5D0-sign(0.5D0,c(5,n,k,l)))*(c(5,n,k,l)*s_in_max(5,n,k,l))&
            + (0.5D0-sign(0.5D0,c(6,n,k,l)))*(c(6,n,k,l)*s_in_max(6,n,k,l))

          c_qin_sum_min_n &
            = (0.5D0-sign(0.5D0,c(1,n,k,l)))*(c(1,n,k,l)*s_in_min(1,n,k,l))&
            + (0.5D0-sign(0.5D0,c(2,n,k,l)))*(c(2,n,k,l)*s_in_min(2,n,k,l))&
            + (0.5D0-sign(0.5D0,c(3,n,k,l)))*(c(3,n,k,l)*s_in_min(3,n,k,l))&
            + (0.5D0-sign(0.5D0,c(4,n,k,l)))*(c(4,n,k,l)*s_in_min(4,n,k,l))&
            + (0.5D0-sign(0.5D0,c(5,n,k,l)))*(c(5,n,k,l)*s_in_min(5,n,k,l))&
            + (0.5D0-sign(0.5D0,c(6,n,k,l)))*(c(6,n,k,l)*s_in_min(6,n,k,l))

          if(abs(c_out_sum_n)<CNST_EPS_ZERO) then
            wrk(n,k,l,s_out_k_min) = s(n,k,l)
            wrk(n,k,l,s_out_k_max) = s(n,k,l)
          else
            wrk(n,k,l,s_out_k_min) = ( &
              s(n,k,l)-c_qin_sum_max_n&
              -s_m1_k_max_n*(1.0D0-c_in_sum_n-c_out_sum_n+d(n,k,l)) &
              )/c_out_sum_n
            wrk(n,k,l,s_out_k_max) = ( &
              s(n,k,l)-c_qin_sum_min_n&
              -s_m1_k_min_n*(1.0D0-c_in_sum_n-c_out_sum_n+d(n,k,l)) &
              )/c_out_sum_n
          end if
       end do !N
    end do !K
    end do !L
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             s_m1_k_min_pl &
                  = min(s_in_min_pl(ADM_gmin_pl,k,l,1),s_in_min_pl(ADM_gmin_pl+1,k,l,1),&
                  s_in_min_pl(ADM_gmin_pl+2,k,l,1),s_in_min_pl(ADM_gmin_pl+3,k,l,1),s_in_min_pl(ADM_gmin_pl+4,k,l,1))
             if(s_m1_k_min_pl== CNST_MAX_REAL) s_m1_k_min_pl=s_pl(ADM_gslf_pl,k,l)
             !          s_m1_k_min_pl = max(SMALL_ZERO,s_m1_k_min_pl)
             s_m1_k_max_pl &
                  = max(s_in_min_pl(ADM_gmin_pl,k,l,1),s_in_min_pl(ADM_gmin_pl+1,k,l,1),&
                  s_in_min_pl(ADM_gmin_pl+2,k,l,1),s_in_min_pl(ADM_gmin_pl+3,k,l,1),s_in_min_pl(ADM_gmin_pl+4,k,l,1))
             if(s_m1_k_max_pl==-CNST_MAX_REAL) s_m1_k_max_pl=s_pl(ADM_gslf_pl,k,l)
             !
             c_in_sum_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*c_pl(ADM_gmin_pl  ,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*c_pl(ADM_gmin_pl+1,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*c_pl(ADM_gmin_pl+2,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*c_pl(ADM_gmin_pl+3,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*c_pl(ADM_gmin_pl+4,k,l)
             c_out_sum_pl &
                  = (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*c_pl(ADM_gmin_pl  ,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*c_pl(ADM_gmin_pl+1,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*c_pl(ADM_gmin_pl+2,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*c_pl(ADM_gmin_pl+3,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*c_pl(ADM_gmin_pl+4,k,l)
             c_qin_sum_max_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*(c_pl(ADM_gmin_pl  ,k,l)*s_in_max_pl(ADM_gmin_pl  ,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*(c_pl(ADM_gmin_pl+1,k,l)*s_in_max_pl(ADM_gmin_pl+1,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*(c_pl(ADM_gmin_pl+2,k,l)*s_in_max_pl(ADM_gmin_pl+2,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*(c_pl(ADM_gmin_pl+3,k,l)*s_in_max_pl(ADM_gmin_pl+3,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*(c_pl(ADM_gmin_pl+4,k,l)*s_in_max_pl(ADM_gmin_pl+4,k,l,1))
             c_qin_sum_min_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl  ,k,l)))*(c_pl(ADM_gmin_pl  ,k,l)*s_in_min_pl(ADM_gmin_pl  ,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+1,k,l)))*(c_pl(ADM_gmin_pl+1,k,l)*s_in_min_pl(ADM_gmin_pl+1,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+2,k,l)))*(c_pl(ADM_gmin_pl+2,k,l)*s_in_min_pl(ADM_gmin_pl+2,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+3,k,l)))*(c_pl(ADM_gmin_pl+3,k,l)*s_in_min_pl(ADM_gmin_pl+3,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_gmin_pl+4,k,l)))*(c_pl(ADM_gmin_pl+4,k,l)*s_in_min_pl(ADM_gmin_pl+4,k,l,1))
             !
             if(abs(c_out_sum_pl)<CNST_EPS_ZERO) then
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_min) = s_pl(ADM_gslf_pl,k,l)
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_max) = s_pl(ADM_gslf_pl,k,l)
             else
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_min) = ( &
                     s_pl(ADM_gslf_pl,k,l)-c_qin_sum_max_pl&
                     -s_m1_k_max_pl*(1.0D0-c_in_sum_pl-c_out_sum_pl+d_pl(ADM_gslf_pl,k,l)) &
                     )/c_out_sum_pl
                wrk_pl(ADM_gslf_pl,k,l,s_out_k_max) = ( &
                     s_pl(ADM_gslf_pl,k,l)-c_qin_sum_min_pl&
                     -s_m1_k_min_pl*(1.0D0-c_in_sum_pl-c_out_sum_pl+d_pl(ADM_gslf_pl,k,l)) &
                     )/c_out_sum_pl
             end if
          end do
       end do
    endif
    !
    !
1000 continue
    !
    !--- H.Tomita 090414 
    do l=1,ADM_lall
      do k=1,ADM_kall
        do n=1,ADM_gall
          wrk(n,k,l,dsx:dsz)=0.0D0
        end do
      end do
    end do

    do l=1,ADM_lall_pl
      do k=1,ADM_kall
        do n=1,ADM_gall_pl
          wrk_pl(n,k,l,dsx:dsz)=0.0D0
        end do
      end do
    end do

    call OPRT_gradient(                    &
         wrk(:,:,:,dsx), wrk_pl(:,:,:,dsx),&
         wrk(:,:,:,dsy), wrk_pl(:,:,:,dsy),&
         wrk(:,:,:,dsz), wrk_pl(:,:,:,dsz),&
         s, s_pl )

    call COMM_data_transfer(wrk,wrk_pl)

    !--- basic scheme
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k = 1, ADM_kall
        nstart = suf(ADM_gmin-1,ADM_gmin-1 )
        nend  = suf(ADM_gmax, ADM_gmax )
        do n = nstart, nend
          sa_p = s(n,k,l)  &
            +wrk(n,k,l,dsx)*(cp(n,k,l,ADM_AI,GRD_XDIR)-GRD_x(n,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n,k,l,dsy)*(cp(n,k,l,ADM_AI,GRD_YDIR)-GRD_x(n,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n,k,l,dsz)*(cp(n,k,l,ADM_AI,GRD_ZDIR)-GRD_x(n,ADM_KNONE,l,GRD_ZDIR))
          sa_m = s(n+1,k,l) &
            +wrk(n+1,k,l,dsx)*(cp(n,k,l,ADM_AI,GRD_XDIR)-GRD_x(n+1,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n+1,k,l,dsy)*(cp(n,k,l,ADM_AI,GRD_YDIR)-GRD_x(n+1,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n+1,k,l,dsz)*(cp(n,k,l,ADM_AI,GRD_ZDIR)-GRD_x(n+1,ADM_KNONE,l,GRD_ZDIR))
          sa(ADM_AI,n,k,l) &
            =(0.5D0+sign(0.5D0,c(1,n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c(1,n,k,l)))*sa_m

          sa_p = s(n,k,l) &
            +wrk(n,k,l,dsx)*(cp(n,k,l,ADM_AIJ,GRD_XDIR)-GRD_x(n,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n,k,l,dsy)*(cp(n,k,l,ADM_AIJ,GRD_YDIR)-GRD_x(n,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n,k,l,dsz)*(cp(n,k,l,ADM_AIJ,GRD_ZDIR)-GRD_x(n,ADM_KNONE,l,GRD_ZDIR))
          sa_m = s(n+1+ADM_gall_1d,k,l) &
            +wrk(n+1+ADM_gall_1d,k,l,dsx)*(cp(n,k,l,ADM_AIJ,GRD_XDIR)- &
            GRD_x(n+1+ADM_gall_1d,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n+1+ADM_gall_1d,k,l,dsy)*(cp(n,k,l,ADM_AIJ,GRD_YDIR)- &
            GRD_x(n+1+ADM_gall_1d,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n+1+ADM_gall_1d,k,l,dsz)*(cp(n,k,l,ADM_AIJ,GRD_ZDIR)- &
            GRD_x(n+1+ADM_gall_1d,ADM_KNONE,l,GRD_ZDIR))
          sa(ADM_AIJ,n,k,l) &
            =(0.5D0+sign(0.5D0,c(2,n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c(2,n,k,l)))*sa_m

          sa_p = s(n,k,l) &
            +wrk(n,k,l,dsx)*(cp(n,k,l,ADM_AJ,GRD_XDIR)-GRD_x(n,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n,k,l,dsy)*(cp(n,k,l,ADM_AJ,GRD_YDIR)-GRD_x(n,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n,k,l,dsz)*(cp(n,k,l,ADM_AJ,GRD_ZDIR)-GRD_x(n,ADM_KNONE,l,GRD_ZDIR))
          sa_m = s(n+ADM_gall_1d,k,l) &
            +wrk(n+ADM_gall_1d,k,l,dsx)*(cp(n,k,l,ADM_AJ,GRD_XDIR)- &
            GRD_x(n+ADM_gall_1d,ADM_KNONE,l,GRD_XDIR)) &
            +wrk(n+ADM_gall_1d,k,l,dsy)*(cp(n,k,l,ADM_AJ,GRD_YDIR)- GRD_x(n+ADM_gall_1d,ADM_KNONE,l,GRD_YDIR)) &
            +wrk(n+ADM_gall_1d,k,l,dsz)*(cp(n,k,l,ADM_AJ,GRD_ZDIR)- GRD_x(n+ADM_gall_1d,ADM_KNONE,l,GRD_ZDIR))
          sa(ADM_AJ,n,k,l) &
            =(0.5D0+sign(0.5D0,c(3,n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c(3,n,k,l)))*sa_m
        end do !N
      end do !K
    end do !L
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_GMAX_PL
                sa_p = s_pl(ADM_gslf_pl,k,l) &
                     +wrk_pl(ADM_gslf_pl,k,l,dsx)*(cp_pl(n,k,l,GRD_XDIR)-GRD_x_pl(ADM_gslf_pl,ADM_KNONE,l,GRD_XDIR)) &
                     +wrk_pl(ADM_gslf_pl,k,l,dsy)*(cp_pl(n,k,l,GRD_YDIR)-GRD_x_pl(ADM_gslf_pl,ADM_KNONE,l,GRD_YDIR)) &
                     +wrk_pl(ADM_gslf_pl,k,l,dsz)*(cp_pl(n,k,l,GRD_ZDIR)-GRD_x_pl(ADM_gslf_pl,ADM_KNONE,l,GRD_ZDIR))
                sa_m = s_pl(n,k,l) &
                     +wrk_pl(n,k,l,dsx)*(cp_pl(n,k,l,GRD_XDIR)-GRD_x_pl(n,ADM_KNONE,l,GRD_XDIR)) &
                     +wrk_pl(n,k,l,dsy)*(cp_pl(n,k,l,GRD_YDIR)-GRD_x_pl(n,ADM_KNONE,l,GRD_YDIR)) &
                     +wrk_pl(n,k,l,dsz)*(cp_pl(n,k,l,GRD_ZDIR)-GRD_x_pl(n,ADM_KNONE,l,GRD_ZDIR))
                sa_pl(n,k,l) &
                     =(0.5D0+sign(0.5D0,c_pl(n,k,l)))*sa_p+(0.5D0-sign(0.5D0,c_pl(n,k,l)))*sa_m
             end do
          end do
       end do
       !
    end if
    !
    if(present(limiter)) then
       if(trim(limiter)=='NON_LIM') goto 2000
    end if
    !
    !---- apply inflow limiter
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
        nstart = suf(ADM_gmin-1,ADM_gmin-1 )
        nend = suf(ADM_gmax ,ADM_gmax )
        do n=nstart,nend
          sa(ADM_AI,n,k,l) &
            =(0.5D0-sign(0.5D0,c(1,n,k,l)))&
            *min(max(sa(ADM_AI,n,k,l),s_in_min(1,n,k,l)),s_in_max(1,n,k,l))&
            +(0.5D0+sign(0.5D0,c(1,n,k,l)))&
            *min(max(sa(ADM_AI,n,k,l),s_in_min(4,n+1,k,l)),s_in_max(4,n+1,k,l))

          sa(ADM_AIJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(2,n,k,l)))&
            *min(max(sa(ADM_AIJ,n,k,l),s_in_min(2,n,k,l)),s_in_max(2,n,k,l))&
            +(0.5D0+sign(0.5D0,c(2,n,k,l)))&
            *min(max(sa(ADM_AIJ,n,k,l),s_in_min(5,n+1+ADM_gall_1d,k,l)),s_in_max(5,n+1+ADM_gall_1d,k,l))

          sa(ADM_AJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(3,n,k,l)))&
            *min(max(sa(ADM_AJ,n,k,l),s_in_min(3,n,k,l)),s_in_max(3,n,k,l))&
            +(0.5D0+sign(0.5D0,c(3,n,k,l)))&
            *min(max(sa(ADM_AJ,n,k,l),s_in_min(6,n+ADM_gall_1d,k,l)),s_in_max(6,n+ADM_gall_1d,k,l))
        end do
      end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n=ADM_gmin_pl,ADM_GMAX_PL
                sa_pl(n,k,l) &
                  =(0.5D0-sign(0.5D0,c_pl(n,k,l)))&
                  *min(max(sa_pl(n,k,l),s_in_min_pl(n,k,l,1)),s_in_max_pl(n,k,l,1))&
                  +(0.5D0+sign(0.5D0,c_pl(n,k,l)))&
                  *min(max(sa_pl(n,k,l),s_in_min_pl(n,k,l,2)),s_in_max_pl(n,k,l,2))
             end do
          end do
       end do
    end if
    !
    !---- apply outflow limitter
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k = 1, ADM_kall
        nstart = suf(ADM_gmin-1,ADM_gmin-1 )
        nend = suf(ADM_gmax ,ADM_gmax )
        do n = nstart,nend
          sa(ADM_AI,n,k,l) &
            =(0.5D0-sign(0.5D0,c(1,n,k,l)))&
            *max(min(sa(ADM_AI,n,k,l),wrk(n+1,k,l,s_out_k_max)),wrk(n+1,k,l,s_out_k_min))&
            +(0.5D0+sign(0.5D0,c(1,n,k,l)))&
            *max(min(sa(ADM_AI,n,k,l),wrk(n,k,l,s_out_k_max)),wrk(n,k,l,s_out_k_min))

          sa(ADM_AIJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(2,n,k,l)))&
            *max(min(sa(ADM_AIJ,n,k,l),wrk(n+1+ADM_gall_1d,k,l,s_out_k_max)),wrk(n+1+ADM_gall_1d,k,l,s_out_k_min))&
            +(0.5D0+sign(0.5D0,c(2,n,k,l)))&
            *max(min(sa(ADM_AIJ,n,k,l),wrk(n,k,l,s_out_k_max)),wrk(n,k,l,s_out_k_min))

          sa(ADM_AJ,n,k,l) &
            =(0.5D0-sign(0.5D0,c(3,n,k,l)))&
            *max(min(sa(ADM_AJ,n,k,l),wrk(n+ADM_gall_1d,k,l,s_out_k_max)),wrk(n+ADM_gall_1d,k,l,s_out_k_min))&
            +(0.5D0+sign(0.5D0,c(3,n,k,l)))&
            *max(min(sa(ADM_AJ,n,k,l),wrk(n,k,l,s_out_k_max)),wrk(n,k,l,s_out_k_min))
        end do
       end do
    end do
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             do n = ADM_gmin_pl,ADM_GMAX_PL
                sa_pl(n,k,l)&
                     =(0.5D0-sign(0.5D0,c_pl(n,k,l)))&
                     *max(min(sa_pl(n,k,l),wrk_pl(n,k,l,s_out_k_max)),wrk_pl(n,k,l,s_out_k_min))&
                     +(0.5D0+sign(0.5D0,c_pl(n,k,l)))&
                     *max(min(sa_pl(n,k,l),wrk_pl(ADM_gslf_pl,k,l,s_out_k_max)),wrk_pl(ADM_gslf_pl,k,l,s_out_k_min))
             end do
          end do
       end do
    end if
    !
    !
2000 continue
    !
    !--- update
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       nstart = suf(ADM_gmin  ,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )
       do k = 1, ADM_kall
          do n = nstart,nend
             scl(n,k,l) = &
                  ( flx_h(1,n,k,l)*sa(ADM_AI,n,k,l)   &
                  + flx_h(2,n,k,l)*sa(ADM_AIJ,n,k,l)  &
                  + flx_h(3,n,k,l)*sa(ADM_AJ,n,k,l)   &
                  + flx_h(4,n,k,l)*sa(ADM_AI,n-1,k,l) &
                  + flx_h(5,n,k,l)*sa(ADM_AIJ,n-1-ADM_gall_1d,k,l) &
                  + flx_h(6,n,k,l)*sa(ADM_AJ,n-ADM_gall_1d,k,l)    &
                  ) * fact
          end do
       end do
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do l=1,ADM_lall_pl
          do k = 1, ADM_kall
             scl_pl(ADM_gslf_pl,k,l)=  &
                  ( flx_h_pl(ADM_gmin_pl  ,k,l)*sa_pl(ADM_gmin_pl  ,k,l) &
                  + flx_h_pl(ADM_gmin_pl+1,k,l)*sa_pl(ADM_gmin_pl+1,k,l) &
                  + flx_h_pl(ADM_gmin_pl+2,k,l)*sa_pl(ADM_gmin_pl+2,k,l) &
                  + flx_h_pl(ADM_gmin_pl+3,k,l)*sa_pl(ADM_gmin_pl+3,k,l) &
                  + flx_h_pl(ADM_gmin_pl+4,k,l)*sa_pl(ADM_gmin_pl+4,k,l) &
                  ) * fact
          end do
       end do
    end if

    return
  end subroutine OPRT_divergence2

end module mod_oprt
!-------------------------------------------------------------------------------
