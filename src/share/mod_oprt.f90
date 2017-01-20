!-------------------------------------------------------------------------------
!> Module operator
!!
!! @par Description
!!          This module contains the subroutines for differential operators
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof

  use mod_adm, only: &
     ADM_nxyz,           &
     TI    => ADM_TI,    &
     TJ    => ADM_TJ,    &
     AI    => ADM_AI,    &
     AIJ   => ADM_AIJ,   &
     AJ    => ADM_AJ,    &
     K0    => ADM_KNONE, &
     vlink => ADM_vlink, &
     ADM_lall,           &
     ADM_lall_pl,        &
     ADM_gall,           &
     ADM_gall_pl,        &
     ADM_kall
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: OPRT_setup
  public :: OPRT_divergence
  public :: OPRT_gradient
  public :: OPRT_laplacian
  public :: OPRT_diffusion
  public :: OPRT_horizontalize_vec
  public :: OPRT_rotation
  public :: OPRT_divdamp

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
#ifdef _FIXEDINDEX_
  real(RP), public              :: OPRT_coef_div    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
  real(RP), public              :: OPRT_coef_div_pl (ADM_nxyz,         0:vlink,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_rot    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
  real(RP), public              :: OPRT_coef_rot_pl (ADM_nxyz,         0:vlink,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
  real(RP), public              :: OPRT_coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_lap    (         ADM_gall,0:6    ,ADM_lall   )
  real(RP), public              :: OPRT_coef_lap_pl (                  0:vlink,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
  real(RP), public              :: OPRT_coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
  real(RP), public              :: OPRT_coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   )
  real(RP), public              :: OPRT_coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl)
#else
  real(RP), public, allocatable :: OPRT_coef_div    (:,:,:,:)   ! coefficient for divergence operator
  real(RP), public, allocatable :: OPRT_coef_div_pl (:,:,:)
  real(RP), public, allocatable :: OPRT_coef_rot    (:,:,:,:)   ! coefficient for rotation operator
  real(RP), public, allocatable :: OPRT_coef_rot_pl (:,:,:)
  real(RP), public, allocatable :: OPRT_coef_grad   (:,:,:,:)   ! coefficient for gradient operator
  real(RP), public, allocatable :: OPRT_coef_grad_pl(:,:,:)
  real(RP), public, allocatable :: OPRT_coef_lap    (:,:,:)     ! coefficient for laplacian operator
  real(RP), public, allocatable :: OPRT_coef_lap_pl (:,:)
  real(RP), public, allocatable :: OPRT_coef_intp   (:,:,:,:,:) ! coefficient for interpolation operator
  real(RP), public, allocatable :: OPRT_coef_intp_pl(:,:,:,:)
  real(RP), public, allocatable :: OPRT_coef_diff   (:,:,:,:)
  real(RP), public, allocatable :: OPRT_coef_diff_pl(:,:,:)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: OPRT_fname   = ''
  character(len=H_SHORT), private :: OPRT_io_mode = 'ADVANCED'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine OPRT_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_gmtr, only: &
       GMTR_p,    &
       GMTR_p_pl, &
       GMTR_t,    &
       GMTR_t_pl, &
       GMTR_a,    &
       GMTR_a_pl
    implicit none

    namelist / OPRTPARAM / &
       OPRT_io_mode, &
       OPRT_fname

    integer :: gall

    integer :: ierr
    integer :: ij, l, d, v
    !---------------------------------------------------------------------------

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[oprt]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=OPRTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** OPRTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist OPRTPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist OPRTPARAM. STOP.'
       call PRC_MPIstop
    endif
    write(IO_FID_LOG,nml=OPRTPARAM)

#ifndef _FIXEDINDEX_
    allocate( OPRT_coef_div    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   ) )
    allocate( OPRT_coef_div_pl (ADM_nxyz,         0:vlink,ADM_lall_pl) )

    allocate( OPRT_coef_rot    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   ) )
    allocate( OPRT_coef_rot_pl (ADM_nxyz,         0:vlink,ADM_lall_pl) )

    allocate( OPRT_coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   ) )
    allocate( OPRT_coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl) )

    allocate( OPRT_coef_lap    (         ADM_gall,0:6    ,ADM_lall   ) )
    allocate( OPRT_coef_lap_pl (                  0:vlink,ADM_lall_pl) )

    allocate( OPRT_coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   ) )
    allocate( OPRT_coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl) )

    allocate( OPRT_coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   ) )
    allocate( OPRT_coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl) )
#endif

    gall = ADM_gall

    do l = 1, ADM_lall
       !$omp parallel do default(none),private(ij,d,v), &
       !$omp shared(l,gall,OPRT_coef_div)
       do ij = 1, gall
       do v  = 0, 6
       do d  = 1, ADM_nxyz
          OPRT_coef_div(d,ij,v,l) = 0.0_RP
       enddo
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,d,v), &
       !$omp shared(l,gall,OPRT_coef_rot)
       do ij = 1, gall
       do v  = 0, 6
       do d  = 1, ADM_nxyz
          OPRT_coef_rot(d,ij,v,l) = 0.0_RP
       enddo
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,d,v), &
       !$omp shared(l,gall,OPRT_coef_grad)
       do ij = 1, gall
       do v  = 0, 6
       do d  = 1, ADM_nxyz
          OPRT_coef_grad(d,ij,v,l) = 0.0_RP
       enddo
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,v), &
       !$omp shared(l,gall,OPRT_coef_lap)
       do ij = 1, gall
       do v = 0, 6
          OPRT_coef_lap(ij,v,l) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel default(none),private(ij,d), &
       !$omp shared(l,gall,OPRT_coef_intp,OPRT_coef_diff)

       !$omp do
       do ij = 1, gall
          do d = 1, ADM_nxyz
             OPRT_coef_intp(d,ij,1,TI,l) = 0.0_RP
             OPRT_coef_intp(d,ij,2,TI,l) = 0.0_RP
             OPRT_coef_intp(d,ij,3,TI,l) = 0.0_RP
          enddo

          do d = 1, ADM_nxyz
             OPRT_coef_intp(d,ij,1,TJ,l) = 0.0_RP
             OPRT_coef_intp(d,ij,2,TJ,l) = 0.0_RP
             OPRT_coef_intp(d,ij,3,TJ,l) = 0.0_RP
          enddo
       enddo
       !$omp end do

       !$omp do
       do ij = 1, gall
       do d  = 1, ADM_nxyz
          OPRT_coef_diff(d,ij,1,l) = 0.0_RP
          OPRT_coef_diff(d,ij,2,l) = 0.0_RP
          OPRT_coef_diff(d,ij,3,l) = 0.0_RP
          OPRT_coef_diff(d,ij,4,l) = 0.0_RP
          OPRT_coef_diff(d,ij,5,l) = 0.0_RP
          OPRT_coef_diff(d,ij,6,l) = 0.0_RP
       enddo
       enddo
       !$omp end do

       !$omp end parallel
    enddo

    call OPRT_divergence_setup( GMTR_p        (:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                GMTR_t        (:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                GMTR_a        (:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_div (:,:,:,:),   OPRT_coef_div_pl (:,:,:)    ) ! [OUT]

    call OPRT_rotation_setup  ( GMTR_p        (:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                GMTR_t        (:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                GMTR_a        (:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_rot (:,:,:,:),   OPRT_coef_rot_pl (:,:,:)    ) ! [OUT]

    call OPRT_gradient_setup  ( GMTR_p        (:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                GMTR_t        (:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                GMTR_a        (:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_grad(:,:,:,:),   OPRT_coef_grad_pl(:,:,:)    ) ! [OUT]

    call OPRT_laplacian_setup ( GMTR_p        (:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                GMTR_t        (:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                GMTR_a        (:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_lap (:,:,:),     OPRT_coef_lap_pl (:,:)      ) ! [OUT]

    call OPRT_diffusion_setup ( GMTR_p        (:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                GMTR_t        (:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                GMTR_a        (:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_intp(:,:,:,:,:), OPRT_coef_intp_pl(:,:,:,:), & ! [OUT]
                                OPRT_coef_diff(:,:,:,:),   OPRT_coef_diff_pl(:,:,:)    ) ! [OUT]

    if ( OPRT_fname /= '' ) then
       call OPRT_output_coef( OPRT_fname )
    endif

    return
  end subroutine OPRT_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_div, coef_div_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       W1      => GMTR_t_W1,    &
       W2      => GMTR_t_W2,    &
       W3      => GMTR_t_W3,    &
       HNX     => GMTR_a_HNX,   &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_div   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(out) :: coef_div_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, hn
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** setup coefficient of divergence operator'

    coef_div   (:,:,:,:) = 0.0_RP
    coef_div_pl(:,  :,:) = 0.0_RP

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       hn = d + HNX - 1

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          ! ij
          coef_div(d,ij,0,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_div(d,ij,1,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_div(d,ij,2,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_div(d,ij,3,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                               ) * 0.5_RP*GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_div(d,ij,4,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_div(d,ij,5,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_div(d,ij,6,l) = ( - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_gmin
          i = ADM_gmin

          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          ! ij
          coef_div(d,ij,0,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_div(d,ij,1,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_div(d,ij,2,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_div(d,ij,3,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_div(d,ij,4,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_div(d,ij,5,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_div(d,ij,6,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       endif

    enddo ! loop d
    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,hn) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,hn) )
          enddo
          coef_div_pl(d,0,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             coef_div_pl(d,v-1,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_divergence_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_rotation_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_rot, coef_rot_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       W1      => GMTR_t_W1,    &
       W2      => GMTR_t_W2,    &
       W3      => GMTR_t_W3,    &
       HTX     => GMTR_a_HTX,   &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_rot   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(out) :: coef_rot_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, ht
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** setup coefficient of rotation operator'

    coef_rot   (:,:,:,:) = 0.0_RP
    coef_rot_pl(:,  :,:) = 0.0_RP

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       ht = d + HTX - 1

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          ! ij
          coef_rot(d,ij,0,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_rot(d,ij,1,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_rot(d,ij,2,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_rot(d,ij,3,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                               ) * 0.5_RP*GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_rot(d,ij,4,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_rot(d,ij,5,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q5 * b5
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_rot(d,ij,6,l) = ( - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q5 * b4
                                 - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q5 * b5
                                 - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ijm1  ,k0,l,AJ ,ht) & ! Q6 * b5
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_gmin
          i = ADM_gmin

          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          ! ij
          coef_rot(d,ij,0,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                                 - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_rot(d,ij,1,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_rot(d,ij,2,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q1 * b6
                                 + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q1 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_rot(d,ij,3,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,ht) & ! Q2 * b1
                                 + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q2 * b2
                                 + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_rot(d,ij,4,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,ht) & ! Q3 * b2
                                 - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q3 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_rot(d,ij,5,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,ht) & ! Q4 * b3
                                 - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q4 * b4
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_rot(d,ij,6,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,ht) & ! Q6 * b4
                                 + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,ht) & ! Q6 * b6
                               ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       endif

    enddo ! loop d
    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          ht = d + HTX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,ht) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,ht) )
          enddo
          coef_rot_pl(d,0,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             coef_rot_pl(d,v-1,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,ht) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,ht) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_rotation_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient_setup( &
       GMTR_p,    GMTR_p_pl,   &
       GMTR_t,    GMTR_t_pl,   &
       GMTR_a,    GMTR_a_pl,   &
       coef_grad, coef_grad_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       W1      => GMTR_t_W1,    &
       W2      => GMTR_t_W2,    &
       W3      => GMTR_t_W3,    &
       HNX     => GMTR_a_HNX,   &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(out) :: coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, hn
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** setup coefficient of gradient operator'

    coef_grad   (:,:,:,:) = 0.0_RP
    coef_grad_pl(:,  :,:) = 0.0_RP

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       hn = d + HNX - 1

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          ! ij
          coef_grad(d,ij,0,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,hn)                    & ! P0 * b1
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,hn)                    & ! P0 * b2
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,hn)                    & ! P0 * b3
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,hn)                    & ! P0 * b4
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,hn)                    & ! P0 * b5
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,hn)                    & ! P0 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_grad(d,ij,1,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_grad(d,ij,2,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_grad(d,ij,3,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_grad(d,ij,4,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_grad(d,ij,5,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_grad(d,ij,6,l) = ( - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(im1jm1,k0,l,TI,W2) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ijm1  ,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_gmin
          i = ADM_gmin

          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          ! ij
          coef_grad(d,ij,0,l) = ( + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(ijm1  ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,hn)                    & ! P0 * b1
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,hn)                    & ! P0 * b2
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,hn)                    & ! P0 * b3
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,hn)                    & ! P0 * b4
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,hn)                    & ! P0 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1j
          coef_grad(d,ij,1,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(ijm1  ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ip1jp1
          coef_grad(d,ij,2,l) = ( + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(ij    ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W2) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijp1
          coef_grad(d,ij,3,l) = ( + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(ij    ,k0,l,TJ,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1j
          coef_grad(d,ij,4,l) = ( + GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(ij    ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(im1j  ,k0,l,TI,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W3) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! im1jm1
          coef_grad(d,ij,5,l) = ( - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(im1jm1,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          ! ijm1
          coef_grad(d,ij,6,l) = ( - GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(im1jm1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(ijm1  ,k0,l,TJ,W1) * GMTR_a(ij    ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
       endif
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             coef = coef + 2.0_RP * ( GMTR_t_pl(ij,k0,l,W1) - 1.0_RP ) * GMTR_a_pl(ijp1,k0,l,hn)
          enddo
          coef_grad_pl(d,0,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             coef_grad_pl(d,v-1,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                     ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_gradient_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_lap, coef_lap_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       T_RAREA => GMTR_t_RAREA, &
       HNX     => GMTR_a_HNX,   &
       TNX     => GMTR_a_TNX,   &
       TN2X    => GMTR_a_TN2X,  &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_lap   (ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(out) :: coef_lap_pl(         0:vlink,ADM_lall_pl)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer  :: i, j, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** setup coefficient of laplacian operator'

    coef_lap   (:,:,:) = 0.0_RP
    coef_lap_pl(  :,:) = 0.0_RP

    do l = 1, ADM_lall

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ! ij
             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1jm1,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(ij,5,l) = coef_lap(ij,5,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             coef_lap(ij,5,l) = coef_lap(ij,5,l) &
                              + GMTR_t(im1jm1,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) )

             ! ijm1
             coef_lap(ij,6,l) = coef_lap(ij,6,l) &
                              + GMTR_t(im1jm1,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) )

             coef_lap(ij,6,l) = coef_lap(ij,6,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ijm1  ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )
          enddo
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_gmin
          i = ADM_gmin

          ij     = (j-1)*ADM_gall_1d + i
          ip1j   = ij + 1
          ip1jp1 = ij + ADM_gall_1d + 1
          ijp1   = ij + ADM_gall_1d
          im1j   = ij - 1
          im1jm1 = ij - ADM_gall_1d - 1
          ijm1   = ij - ADM_gall_1d

          coef_lap(ij,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ! ij
             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             coef_lap(ij,0,l) = coef_lap(ij,0,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )

             coef_lap(ij,1,l) = coef_lap(ij,1,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ip1j  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) )

             coef_lap(ij,2,l) = coef_lap(ij,2,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(ij    ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(ijp1  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) )

             coef_lap(ij,3,l) = coef_lap(ij,3,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1j  ,k0,l,TI,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AJ ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) )

             coef_lap(ij,4,l) = coef_lap(ij,4,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 2.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(ij,5,l) = coef_lap(ij,5,l) &
                              + GMTR_t(im1jm1,k0,l,TJ,T_RAREA) &
                              * ( - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1j  ,k0,l,AI ,hn) &
                                  - 1.0_RP * GMTR_a(im1jm1,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(im1j  ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 1.0_RP * GMTR_a(im1jm1,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) )

             ! ijm1
             coef_lap(ij,6,l) = coef_lap(ij,6,l) &
                              + GMTR_t(ijm1  ,k0,l,TJ,T_RAREA) &
                              * ( + 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  + 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(im1jm1,k0,l,AIJ,hn) &
                                  - 1.0_RP * GMTR_a(ijm1  ,k0,l,AIJ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  - 2.0_RP * GMTR_a(ij    ,k0,l,AI ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) &
                                  + 1.0_RP * GMTR_a(ijm1  ,k0,l,AJ ,tn) * GMTR_a(ij    ,k0,l,AI ,hn) )
          enddo
       endif

       coef_lap(:,0,l) = coef_lap(:,0,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,1,l) = coef_lap(:,1,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,2,l) = coef_lap(:,2,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,3,l) = coef_lap(:,3,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,4,l) = coef_lap(:,4,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,5,l) = coef_lap(:,5,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,6,l) = coef_lap(:,6,l) * GMTR_p(:,k0,l,P_RAREA) / 12.0_RP
    enddo

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1,ADM_lall_pl

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                   * ( - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) )

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                   * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) )
             enddo
          enddo ! d loop

          do v = ADM_gslf_pl, ADM_gmax_pl
             coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) * GMTR_p_pl(n,k0,l,P_RAREA) / 12.0_RP
          enddo

       enddo ! l loop
    endif

    return
  end subroutine OPRT_laplacian_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion_setup( &
       GMTR_p,    GMTR_p_pl,    &
       GMTR_t,    GMTR_t_pl,    &
       GMTR_a,    GMTR_a_pl,    &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       T_RAREA => GMTR_t_RAREA, &
       HNX     => GMTR_a_HNX,   &
       TNX     => GMTR_a_TNX,   &
       TN2X    => GMTR_a_TN2X,  &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
    real(RP), intent(out) :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
    real(RP), intent(out) :: coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   )
    real(RP), intent(out) :: coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl)

    integer :: ij
    integer :: ip1j, ijp1
    integer :: im1j, ijm1, im1jm1

    integer :: i, j, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** setup coefficient of diffusion operator'

    coef_intp   (:,:,:,:,:) = 0.0_RP
    coef_intp_pl(:,:,:,  :) = 0.0_RP
    coef_diff   (:,:,:,:)   = 0.0_RP
    coef_diff_pl(:,  :,:)   = 0.0_RP

    do l = 1, ADM_lall
       do d = 1, ADM_nxyz
          hn = d + HNX - 1
          tn = d + TNX - 1

          do j = ADM_gmin-1, ADM_gmax
          do i = ADM_gmin-1, ADM_gmax
             ij     = (j-1)*ADM_gall_1d + i
             ip1j   = ij + 1
             ijp1   = ij + ADM_gall_1d

             coef_intp(d,ij,1,TI,l) = + GMTR_a(ij  ,k0,l,AIJ,tn) - GMTR_a(ij  ,k0,l,AI ,tn)
             coef_intp(d,ij,2,TI,l) = - GMTR_a(ij  ,k0,l,AI ,tn) - GMTR_a(ip1j,k0,l,AJ ,tn)
             coef_intp(d,ij,3,TI,l) = - GMTR_a(ip1j,k0,l,AJ ,tn) + GMTR_a(ij  ,k0,l,AIJ,tn)

             coef_intp(d,ij,1,TJ,l) = + GMTR_a(ij  ,k0,l,AJ ,tn) - GMTR_a(ij  ,k0,l,AIJ,tn)
             coef_intp(d,ij,2,TJ,l) = - GMTR_a(ij  ,k0,l,AIJ,tn) + GMTR_a(ijp1,k0,l,AI ,tn)
             coef_intp(d,ij,3,TJ,l) = + GMTR_a(ijp1,k0,l,AI ,tn) + GMTR_a(ij  ,k0,l,AJ ,tn)

             coef_intp(d,ij,:,TI,l) = coef_intp(d,ij,:,TI,l) * 0.5_RP * GMTR_t(ij,k0,l,TI,T_RAREA)
             coef_intp(d,ij,:,TJ,l) = coef_intp(d,ij,:,TJ,l) * 0.5_RP * GMTR_t(ij,k0,l,TJ,T_RAREA)
          enddo
          enddo

          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij     = (j-1)*ADM_gall_1d + i
             im1j   = ij - 1
             im1jm1 = ij - ADM_gall_1d - 1
             ijm1   = ij - ADM_gall_1d

             coef_diff(d,ij,1,l) = + GMTR_a(ij    ,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(d,ij,2,l) = + GMTR_a(ij    ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(d,ij,3,l) = - GMTR_a(im1j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(d,ij,4,l) = - GMTR_a(im1jm1,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(d,ij,5,l) = - GMTR_a(ijm1  ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(d,ij,6,l) = + GMTR_a(ij    ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
          enddo
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_gmin
          i = ADM_gmin

          ij = (j-1)*ADM_gall_1d + i

          coef_diff(:,ij,5,l) = 0.0_RP
       endif

    enddo ! l loop

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1, ADM_lall_pl

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

                coef_intp_pl(d,v,1,l) = - GMTR_a_pl(ijp1,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn )
                coef_intp_pl(d,v,2,l) = + GMTR_a_pl(ij  ,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn2)
                coef_intp_pl(d,v,3,l) = + GMTR_a_pl(ij  ,k0,l,tn2) - GMTR_a_pl(ijp1,k0,l,tn )

                coef_intp_pl(d,v,:,l) = coef_intp_pl(d,v,:,l) * 0.5_RP * GMTR_t_pl(v,k0,l,T_RAREA)

                coef_diff_pl(d,v-1,l) = GMTR_a_pl(v,k0,l,hn) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
             enddo
          enddo

       enddo ! l loop
    endif

    return
  end subroutine OPRT_diffusion_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence( &
       scl,      scl_pl,     &
       vx,       vx_pl,      &
       vy,       vy_pl,      &
       vz,       vz_pl,      &
       coef_div, coef_div_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: scl        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: scl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_div   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_div_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)

    integer  :: gmin, gmax, kall, iall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_divergence',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    !$omp parallel workshare
    scl(:,:,:) = 0.0_RP
    !$omp end parallel workshare

    do l = 1, ADM_lall
       !$omp parallel do default(none),private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(l,gmin,gmax,kall,iall,scl,vx,vy,vz,coef_div)
       do k = 1, kall
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                scl(ij,k,l) = scl(ij,k,l) + ( coef_div(XDIR,ij,0,l) * vx(ij    ,k,l) &
                                            + coef_div(XDIR,ij,1,l) * vx(ip1j  ,k,l) &
                                            + coef_div(XDIR,ij,2,l) * vx(ip1jp1,k,l) &
                                            + coef_div(XDIR,ij,3,l) * vx(ijp1  ,k,l) &
                                            + coef_div(XDIR,ij,4,l) * vx(im1j  ,k,l) &
                                            + coef_div(XDIR,ij,5,l) * vx(im1jm1,k,l) &
                                            + coef_div(XDIR,ij,6,l) * vx(ijm1  ,k,l) )
             enddo

             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                scl(ij,k,l) = scl(ij,k,l) + ( coef_div(YDIR,ij,0,l) * vy(ij    ,k,l) &
                                            + coef_div(YDIR,ij,1,l) * vy(ip1j  ,k,l) &
                                            + coef_div(YDIR,ij,2,l) * vy(ip1jp1,k,l) &
                                            + coef_div(YDIR,ij,3,l) * vy(ijp1  ,k,l) &
                                            + coef_div(YDIR,ij,4,l) * vy(im1j  ,k,l) &
                                            + coef_div(YDIR,ij,5,l) * vy(im1jm1,k,l) &
                                            + coef_div(YDIR,ij,6,l) * vy(ijm1  ,k,l) )
             enddo

             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                scl(ij,k,l) = scl(ij,k,l) + ( coef_div(ZDIR,ij,0,l) * vz(ij    ,k,l) &
                                            + coef_div(ZDIR,ij,1,l) * vz(ip1j  ,k,l) &
                                            + coef_div(ZDIR,ij,2,l) * vz(ip1jp1,k,l) &
                                            + coef_div(ZDIR,ij,3,l) * vz(ijp1  ,k,l) &
                                            + coef_div(ZDIR,ij,4,l) * vz(im1j  ,k,l) &
                                            + coef_div(ZDIR,ij,5,l) * vz(im1jm1,k,l) &
                                            + coef_div(ZDIR,ij,6,l) * vz(ijm1  ,k,l) )
             enddo
          enddo
       enddo
       !$omp end parallel do
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          scl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( coef_div_pl(XDIR,v-1,l) * vx_pl(v,k,l) &
                                             + coef_div_pl(YDIR,v-1,l) * vy_pl(v,k,l) &
                                             + coef_div_pl(ZDIR,v-1,l) * vz_pl(v,k,l) )
          enddo
       enddo
       enddo
    else
       scl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_divergence',2)

    return
  end subroutine OPRT_divergence

  !-----------------------------------------------------------------------------
  subroutine OPRT_rotation( &
       scl,      scl_pl,     &
       vx,       vx_pl,      &
       vy,       vy_pl,      &
       vz,       vz_pl,      &
       coef_rot, coef_rot_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: scl        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: scl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_rot   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_rot_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)

    integer  :: gmin, gmax, kall, iall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_rotation',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    do l = 1, ADM_lall
       !$omp parallel default(none),private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(l,gmin,gmax,kall,iall,scl,vx,vy,vz,coef_rot)
       do k = 1, kall

          !$omp workshare
          scl(:,k,l) = 0.0_RP
          !$omp end workshare

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,0,l) * vx(ij    ,k,l) &
                                            + coef_rot(YDIR,ij,0,l) * vy(ij    ,k,l) &
                                            + coef_rot(ZDIR,ij,0,l) * vz(ij    ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,1,l) * vx(ip1j  ,k,l) &
                                            + coef_rot(YDIR,ij,1,l) * vy(ip1j  ,k,l) &
                                            + coef_rot(ZDIR,ij,1,l) * vz(ip1j  ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1jp1 = ij + iall + 1

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,2,l) * vx(ip1jp1,k,l) &
                                            + coef_rot(YDIR,ij,2,l) * vy(ip1jp1,k,l) &
                                            + coef_rot(ZDIR,ij,2,l) * vz(ip1jp1,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ijp1   = ij + iall

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,3,l) * vx(ijp1  ,k,l) &
                                            + coef_rot(YDIR,ij,3,l) * vy(ijp1  ,k,l) &
                                            + coef_rot(ZDIR,ij,3,l) * vz(ijp1  ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                im1j   = ij - 1

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,4,l) * vx(im1j  ,k,l) &
                                            + coef_rot(YDIR,ij,4,l) * vy(im1j  ,k,l) &
                                            + coef_rot(ZDIR,ij,4,l) * vz(im1j  ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                im1jm1 = ij - iall -1

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,5,l) * vx(im1jm1,k,l) &
                                            + coef_rot(YDIR,ij,5,l) * vy(im1jm1,k,l) &
                                            + coef_rot(ZDIR,ij,5,l) * vz(im1jm1,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ijm1   = ij - iall

                scl(ij,k,l) = scl(ij,k,l) + ( coef_rot(XDIR,ij,6,l) * vx(ijm1  ,k,l) &
                                            + coef_rot(YDIR,ij,6,l) * vy(ijm1  ,k,l) &
                                            + coef_rot(ZDIR,ij,6,l) * vz(ijm1  ,k,l) )
             enddo
          enddo
          !$omp end do

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          scl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( coef_rot_pl(XDIR,v-1,l) * vx_pl(v,k,l) &
                                             + coef_rot_pl(YDIR,v-1,l) * vy_pl(v,k,l) &
                                             + coef_rot_pl(ZDIR,v-1,l) * vz_pl(v,k,l) )
          enddo
       enddo
       enddo
    else
       scl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_rotation',2)

    return
  end subroutine OPRT_rotation

  !-----------------------------------------------------------------------------
  !> horizontal gradient operator
  subroutine OPRT_gradient( &
       grad,      grad_pl,     &
       scl,       scl_pl,      &
       coef_grad, coef_grad_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: grad        (ADM_gall   ,ADM_kall,ADM_lall   ,ADM_nxyz)
    real(RP), intent(out) :: grad_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,ADM_nxyz)
    real(RP), intent(in)  :: scl         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)

    integer  :: gmin, gmax, kall, iall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_gradient',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    !$omp parallel workshare
    grad(:,:,:,:) = 0.0_RP
    !$omp end parallel workshare

    do l = 1, ADM_lall
       !$omp parallel do default(none),private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(l,gmin,gmax,kall,iall,grad,scl,coef_grad)
       do k = 1, kall
          do j = gmin, gmax
             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                grad(ij,k,l,XDIR) = ( coef_grad(XDIR,ij,0,l) * scl(ij    ,k,l) &
                                    + coef_grad(XDIR,ij,1,l) * scl(ip1j  ,k,l) &
                                    + coef_grad(XDIR,ij,2,l) * scl(ip1jp1,k,l) &
                                    + coef_grad(XDIR,ij,3,l) * scl(ijp1  ,k,l) &
                                    + coef_grad(XDIR,ij,4,l) * scl(im1j  ,k,l) &
                                    + coef_grad(XDIR,ij,5,l) * scl(im1jm1,k,l) &
                                    + coef_grad(XDIR,ij,6,l) * scl(ijm1  ,k,l) )
             enddo

             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                grad(ij,k,l,YDIR) = ( coef_grad(YDIR,ij,0,l) * scl(ij    ,k,l) &
                                    + coef_grad(YDIR,ij,1,l) * scl(ip1j  ,k,l) &
                                    + coef_grad(YDIR,ij,2,l) * scl(ip1jp1,k,l) &
                                    + coef_grad(YDIR,ij,3,l) * scl(ijp1  ,k,l) &
                                    + coef_grad(YDIR,ij,4,l) * scl(im1j  ,k,l) &
                                    + coef_grad(YDIR,ij,5,l) * scl(im1jm1,k,l) &
                                    + coef_grad(YDIR,ij,6,l) * scl(ijm1  ,k,l) )
             enddo

             do i = gmin, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall
                im1j   = ij - 1
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                grad(ij,k,l,ZDIR) = ( coef_grad(ZDIR,ij,0,l) * scl(ij    ,k,l) &
                                    + coef_grad(ZDIR,ij,1,l) * scl(ip1j  ,k,l) &
                                    + coef_grad(ZDIR,ij,2,l) * scl(ip1jp1,k,l) &
                                    + coef_grad(ZDIR,ij,3,l) * scl(ijp1  ,k,l) &
                                    + coef_grad(ZDIR,ij,4,l) * scl(im1j  ,k,l) &
                                    + coef_grad(ZDIR,ij,5,l) * scl(im1jm1,k,l) &
                                    + coef_grad(ZDIR,ij,6,l) * scl(ijm1  ,k,l) )
             enddo
          enddo
       enddo
       !$omp end parallel do
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          grad_pl(:,k,l,XDIR) = 0.0_RP
          grad_pl(:,k,l,YDIR) = 0.0_RP
          grad_pl(:,k,l,ZDIR) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             grad_pl(n,k,l,XDIR) = grad_pl(n,k,l,XDIR) + coef_grad_pl(XDIR,v-1,l) * scl_pl(v,k,l)
             grad_pl(n,k,l,YDIR) = grad_pl(n,k,l,YDIR) + coef_grad_pl(YDIR,v-1,l) * scl_pl(v,k,l)
             grad_pl(n,k,l,ZDIR) = grad_pl(n,k,l,ZDIR) + coef_grad_pl(ZDIR,v-1,l) * scl_pl(v,k,l)
          enddo
       enddo
       enddo
    else
       grad_pl(:,:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_gradient',2)

    return
  end subroutine OPRT_gradient

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian( &
       dscl,     dscl_pl,    &
       scl,      scl_pl,     &
       coef_lap, coef_lap_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall_1d, &
       ADM_gmin,    &
       ADM_gmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    implicit none

    real(RP), intent(out) :: dscl       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_lap   (ADM_gall   ,0:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_lap_pl(            0:vlink,ADM_lall_pl)

    integer  :: gmin, gmax, kall, iall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_laplacian',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    !$omp parallel workshare
    dscl(:,:,:) = 0.0_RP
    !$omp end parallel workshare

    do l = 1, ADM_lall
       !$omp parallel do default(none),private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(l,gmin,gmax,kall,iall,dscl,scl,coef_lap)
       do k = 1, kall
          do j = gmin, gmax
          do i = gmin, gmax
             ij     = (j-1)*iall + i
             ip1j   = ij + 1
             ip1jp1 = ij + iall + 1
             ijp1   = ij + iall
             im1j   = ij - 1
             im1jm1 = ij - iall -1
             ijm1   = ij - iall

             dscl(ij,k,l) = ( coef_lap(ij,0,l) * scl(ij    ,k,l) &
                            + coef_lap(ij,1,l) * scl(ip1j  ,k,l) &
                            + coef_lap(ij,2,l) * scl(ip1jp1,k,l) &
                            + coef_lap(ij,3,l) * scl(ijp1  ,k,l) &
                            + coef_lap(ij,4,l) * scl(im1j  ,k,l) &
                            + coef_lap(ij,5,l) * scl(im1jm1,k,l) &
                            + coef_lap(ij,6,l) * scl(ijm1  ,k,l) )
          enddo
          enddo
       enddo
       !$omp end parallel do
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          dscl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + coef_lap_pl(v-1,l) * scl_pl(v,k,l)
          enddo
       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_laplacian',2)

    return
  end subroutine OPRT_laplacian

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion( &
       dscl,      dscl_pl,      &
       scl,       scl_pl,       &
       kh,        kh_pl,        &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: dscl        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl)

    real(RP) :: vt   (ADM_nxyz,ADM_gall   ,TI:TJ)
    real(RP) :: vt_pl(ADM_nxyz,ADM_gall_pl)

    integer  :: gmin, gmax, kall, iall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: i, j, k, l, d, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_diffusion',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    !$omp parallel workshare
    dscl(:,:,:) = 0.0_RP
    !$omp end parallel workshare

    do l = 1, ADM_lall
       !$omp parallel default(none),private(i,j,k,d,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,kall,iall,dscl,scl,kh,vt,coef_intp,coef_diff)
       do k = 1, kall

          !$omp do
          do j = gmin-1, gmax
             do i = gmin-1, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1

                do d = 1, 3
                   vt(d,ij,TI) = ( ( + 2.0_RP * coef_intp(d,ij,1,TI,l) &
                                     - 1.0_RP * coef_intp(d,ij,2,TI,l) &
                                     - 1.0_RP * coef_intp(d,ij,3,TI,l) ) * scl(ij    ,k,l) &
                                 + ( - 1.0_RP * coef_intp(d,ij,1,TI,l) &
                                     + 2.0_RP * coef_intp(d,ij,2,TI,l) &
                                     - 1.0_RP * coef_intp(d,ij,3,TI,l) ) * scl(ip1j  ,k,l) &
                                 + ( - 1.0_RP * coef_intp(d,ij,1,TI,l) &
                                     - 1.0_RP * coef_intp(d,ij,2,TI,l) &
                                     + 2.0_RP * coef_intp(d,ij,3,TI,l) ) * scl(ip1jp1,k,l) &
                                 ) / 3.0_RP
                enddo
             enddo

             do i = gmin-1, gmax
                ij     = (j-1)*iall + i
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall

                do d = 1, 3
                   vt(d,ij,TJ) = ( ( + 2.0_RP * coef_intp(d,ij,1,TJ,l) &
                                     - 1.0_RP * coef_intp(d,ij,2,TJ,l) &
                                     - 1.0_RP * coef_intp(d,ij,3,TJ,l) ) * scl(ij    ,k,l) &
                                 + ( - 1.0_RP * coef_intp(d,ij,1,TJ,l) &
                                     + 2.0_RP * coef_intp(d,ij,2,TJ,l) &
                                     - 1.0_RP * coef_intp(d,ij,3,TJ,l) ) * scl(ip1jp1,k,l) &
                                 + ( - 1.0_RP * coef_intp(d,ij,1,TJ,l) &
                                     - 1.0_RP * coef_intp(d,ij,2,TJ,l) &
                                     + 2.0_RP * coef_intp(d,ij,3,TJ,l) ) * scl(ijp1  ,k,l) &
                                 ) / 3.0_RP
                enddo
             enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
             j = gmin-1
             i = gmin-1

             ij   = (j-1)*iall + i
             ip1j = ij + 1

             vt(:,ij,TI) = vt(:,ip1j,TJ)
             !$omp end master
          endif

          !$omp do
          do j = gmin  , gmax
             do i = gmin  , gmax
                ij     = (j-1)*iall + i
                ip1jp1 = ij + iall + 1

                dscl(ij,k,l) = dscl(ij,k,l) &
                             + ( coef_diff(XDIR,ij,1,l) * ( vt(XDIR,ij    ,TI) + vt(XDIR,ij    ,TJ) ) &
                               + coef_diff(YDIR,ij,1,l) * ( vt(YDIR,ij    ,TI) + vt(YDIR,ij    ,TJ) ) &
                               + coef_diff(ZDIR,ij,1,l) * ( vt(ZDIR,ij    ,TI) + vt(ZDIR,ij    ,TJ) ) &
                               ) * 0.5_RP * ( kh(ij    ,k,l) + kh(ip1jp1,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
             do i = gmin  , gmax
                ij     = (j-1)*iall + i
                ijp1   = ij + iall
                im1j   = ij - 1

                dscl(ij,k,l) = dscl(ij,k,l) &
                             + ( coef_diff(XDIR,ij,2,l) * ( vt(XDIR,ij    ,TJ) + vt(XDIR,im1j  ,TI) ) &
                               + coef_diff(YDIR,ij,2,l) * ( vt(YDIR,ij    ,TJ) + vt(YDIR,im1j  ,TI) ) &
                               + coef_diff(ZDIR,ij,2,l) * ( vt(ZDIR,ij    ,TJ) + vt(ZDIR,im1j  ,TI) ) &
                               ) * 0.5_RP * ( kh(ij    ,k,l) + kh(ijp1  ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
             do i = gmin  , gmax
                ij     = (j-1)*iall + i
                im1j   = ij - 1
                im1jm1 = ij - iall -1

                dscl(ij,k,l) = dscl(ij,k,l) &
                             + ( coef_diff(XDIR,ij,3,l) * ( vt(XDIR,im1j  ,TI) + vt(XDIR,im1jm1,TJ) ) &
                               + coef_diff(YDIR,ij,3,l) * ( vt(YDIR,im1j  ,TI) + vt(YDIR,im1jm1,TJ) ) &
                               + coef_diff(ZDIR,ij,3,l) * ( vt(ZDIR,im1j  ,TI) + vt(ZDIR,im1jm1,TJ) ) &
                               ) * 0.5_RP * ( kh(im1j  ,k,l) + kh(ij    ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
             do i = gmin  , gmax
                ij     = (j-1)*iall + i
                im1jm1 = ij - iall -1

                dscl(ij,k,l) = dscl(ij,k,l) &
                             + ( coef_diff(XDIR,ij,4,l) * ( vt(XDIR,im1jm1,TJ) + vt(XDIR,im1jm1,TI) ) &
                               + coef_diff(YDIR,ij,4,l) * ( vt(YDIR,im1jm1,TJ) + vt(YDIR,im1jm1,TI) ) &
                               + coef_diff(ZDIR,ij,4,l) * ( vt(ZDIR,im1jm1,TJ) + vt(ZDIR,im1jm1,TI) ) &
                               ) * 0.5_RP * ( kh(im1jm1,k,l) + kh(ij    ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
             do i = gmin  , gmax
                ij     = (j-1)*iall + i
                im1jm1 = ij - iall -1
                ijm1   = ij - iall

                dscl(ij,k,l) = dscl(ij,k,l) &
                             + ( coef_diff(XDIR,ij,5,l) * ( vt(XDIR,im1jm1,TI) + vt(XDIR,ijm1  ,TJ) ) &
                               + coef_diff(YDIR,ij,5,l) * ( vt(YDIR,im1jm1,TI) + vt(YDIR,ijm1  ,TJ) ) &
                               + coef_diff(ZDIR,ij,5,l) * ( vt(ZDIR,im1jm1,TI) + vt(ZDIR,ijm1  ,TJ) ) &
                               ) * 0.5_RP * ( kh(ijm1  ,k,l) + kh(ij    ,k,l) )
             enddo
          enddo
          !$omp end do

          !$omp do
          do j = gmin  , gmax
             do i = gmin  , gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ijm1   = ij - iall

                dscl(ij,k,l) = dscl(ij,k,l)  &
                             + ( coef_diff(XDIR,ij,6,l) * ( vt(XDIR,ijm1  ,TJ) + vt(XDIR,ij    ,TI) ) &
                               + coef_diff(YDIR,ij,6,l) * ( vt(YDIR,ijm1  ,TJ) + vt(YDIR,ij    ,TI) ) &
                               + coef_diff(ZDIR,ij,6,l) * ( vt(ZDIR,ijm1  ,TJ) + vt(ZDIR,ij    ,TI) ) &
                               ) * 0.5_RP * ( kh(ij    ,k,l) + kh(ip1j  ,k,l) )
             enddo
          enddo
          !$omp end do

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                vt_pl(d,ij) = ( ( + 2.0_RP * coef_intp_pl(d,v,1,l) &
                                  - 1.0_RP * coef_intp_pl(d,v,2,l) &
                                  - 1.0_RP * coef_intp_pl(d,v,3,l) ) * scl_pl(n   ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(d,v,1,l) &
                                  + 2.0_RP * coef_intp_pl(d,v,2,l) &
                                  - 1.0_RP * coef_intp_pl(d,v,3,l) ) * scl_pl(ij  ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(d,v,1,l) &
                                  - 1.0_RP * coef_intp_pl(d,v,2,l) &
                                  + 2.0_RP * coef_intp_pl(d,v,3,l) ) * scl_pl(ijp1,k,l) &
                              ) / 3.0_RP
             enddo
          enddo

          dscl_pl(:,k,l) = 0.0_RP

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

             dscl_pl(n,k,l) = dscl_pl(n,k,l) &
                            + ( coef_diff_pl(XDIR,v-1,l) * ( vt_pl(XDIR,ijm1) + vt_pl(XDIR,ij) ) &
                              + coef_diff_pl(YDIR,v-1,l) * ( vt_pl(YDIR,ijm1) + vt_pl(YDIR,ij) ) &
                              + coef_diff_pl(ZDIR,v-1,l) * ( vt_pl(ZDIR,ijm1) + vt_pl(ZDIR,ij) ) &
                              ) * 0.5_RP * ( kh_pl(n,k,l) + kh_pl(ij,k,l) )
          enddo

       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_diffusion',2)

    return
  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------
  subroutine OPRT_divdamp( &
       ddivdx,    ddivdx_pl,    &
       ddivdy,    ddivdy_pl,    &
       ddivdz,    ddivdz_pl,    &
       vx,        vx_pl,        &
       vy,        vy_pl,        &
       vz,        vz_pl,        &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: ddivdx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl)

    real(RP) :: sclt   (ADM_gall   ,TI:TJ)
    real(RP) :: sclt_pl(ADM_gall_pl)

    integer  :: gmin, gmax, kall, iall

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_divdamp',2)

    gmin = ADM_gmin
    gmax = ADM_gmax
    kall = ADM_kall
    iall = ADM_gall_1d

    !$omp parallel workshare
    ddivdx(:,:,:) = 0.0_RP
    ddivdy(:,:,:) = 0.0_RP
    ddivdz(:,:,:) = 0.0_RP
    !$omp end parallel workshare

    do l = 1, ADM_lall
       !$omp parallel default(none),private(i,j,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(l,ADM_have_sgp,gmin,gmax,kall,iall,ddivdx,ddivdy,ddivdz,vx,vy,vz,sclt,coef_intp,coef_diff)
       do k = 1, kall

          !$omp do
          do j = gmin-1, gmax
             do i = gmin-1, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall

                sclt(ij,TI) = coef_intp(XDIR,ij,1,TI,l) * vx(ij    ,k,l) &
                            + coef_intp(XDIR,ij,2,TI,l) * vx(ip1j  ,k,l) &
                            + coef_intp(XDIR,ij,3,TI,l) * vx(ip1jp1,k,l) &
                            + coef_intp(YDIR,ij,1,TI,l) * vy(ij    ,k,l) &
                            + coef_intp(YDIR,ij,2,TI,l) * vy(ip1j  ,k,l) &
                            + coef_intp(YDIR,ij,3,TI,l) * vy(ip1jp1,k,l) &
                            + coef_intp(ZDIR,ij,1,TI,l) * vz(ij    ,k,l) &
                            + coef_intp(ZDIR,ij,2,TI,l) * vz(ip1j  ,k,l) &
                            + coef_intp(ZDIR,ij,3,TI,l) * vz(ip1jp1,k,l)
             enddo

             do i = gmin-1, gmax
                ij     = (j-1)*iall + i
                ip1j   = ij + 1
                ip1jp1 = ij + iall + 1
                ijp1   = ij + iall

                sclt(ij,TJ) = coef_intp(XDIR,ij,1,TJ,l) * vx(ij    ,k,l) &
                            + coef_intp(XDIR,ij,2,TJ,l) * vx(ip1jp1,k,l) &
                            + coef_intp(XDIR,ij,3,TJ,l) * vx(ijp1  ,k,l) &
                            + coef_intp(YDIR,ij,1,TJ,l) * vy(ij    ,k,l) &
                            + coef_intp(YDIR,ij,2,TJ,l) * vy(ip1jp1,k,l) &
                            + coef_intp(YDIR,ij,3,TJ,l) * vy(ijp1  ,k,l) &
                            + coef_intp(ZDIR,ij,1,TJ,l) * vz(ij    ,k,l) &
                            + coef_intp(ZDIR,ij,2,TJ,l) * vz(ip1jp1,k,l) &
                            + coef_intp(ZDIR,ij,3,TJ,l) * vz(ijp1  ,k,l)
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
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

                sclt_pl(ij) = coef_intp_pl(XDIR,v,1,l) * vx_pl(n   ,k,l) &
                            + coef_intp_pl(XDIR,v,2,l) * vx_pl(ij  ,k,l) &
                            + coef_intp_pl(XDIR,v,3,l) * vx_pl(ijp1,k,l) &
                            + coef_intp_pl(YDIR,v,1,l) * vy_pl(n   ,k,l) &
                            + coef_intp_pl(YDIR,v,2,l) * vy_pl(ij  ,k,l) &
                            + coef_intp_pl(YDIR,v,3,l) * vy_pl(ijp1,k,l) &
                            + coef_intp_pl(ZDIR,v,1,l) * vz_pl(n   ,k,l) &
                            + coef_intp_pl(ZDIR,v,2,l) * vz_pl(ij  ,k,l) &
                            + coef_intp_pl(ZDIR,v,3,l) * vz_pl(ijp1,k,l)
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
       enddo
    else
       ddivdx_pl(:,:,:) = 0.0_RP
       ddivdy_pl(:,:,:) = 0.0_RP
       ddivdz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_divdamp',2)

    return
  end subroutine OPRT_divdamp

  !-----------------------------------------------------------------------------
  subroutine OPRT_horizontalize_vec( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_have_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR, &
       GRD_rscale,       &
       GRD_x,            &
       GRD_x_pl
    implicit none

    real(RP), intent(inout) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer  :: gall, kall

    real(RP) :: rscale
    real(RP) :: prd
    integer  :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_horizontalize_vec',2)

    rscale = GRD_rscale
    kall   = ADM_kall
    gall   = ADM_gall

    do l = 1, ADM_lall
       !$omp parallel do default(none), private(g,k,prd), &
       !$omp shared(l,gall,kall,vx,vy,vz,GRD_x,rscale)
       do k = 1, kall
       do g = 1, gall
          prd = vx(g,k,l) * GRD_x(g,1,l,XDIR) / rscale &
              + vy(g,k,l) * GRD_x(g,1,l,YDIR) / rscale &
              + vz(g,k,l) * GRD_x(g,1,l,ZDIR) / rscale

          vx(g,k,l) = vx(g,k,l) - prd * GRD_x(g,1,l,XDIR) / rscale
          vy(g,k,l) = vy(g,k,l) - prd * GRD_x(g,1,l,YDIR) / rscale
          vz(g,k,l) = vz(g,k,l) - prd * GRD_x(g,1,l,ZDIR) / rscale
       enddo
       enddo
       !$omp end parallel do
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             prd = vx_pl(g,k,l) * GRD_x_pl(g,1,l,XDIR) / rscale &
                 + vy_pl(g,k,l) * GRD_x_pl(g,1,l,YDIR) / rscale &
                 + vz_pl(g,k,l) * GRD_x_pl(g,1,l,ZDIR) / rscale

             vx_pl(g,k,l) = vx_pl(g,k,l) - prd * GRD_x_pl(g,1,l,XDIR) / rscale
             vy_pl(g,k,l) = vy_pl(g,k,l) - prd * GRD_x_pl(g,1,l,YDIR) / rscale
             vz_pl(g,k,l) = vz_pl(g,k,l) - prd * GRD_x_pl(g,1,l,ZDIR) / rscale
          enddo
          enddo
       enddo
    else
       vx_pl(:,:,:) = 0.0_RP
       vy_pl(:,:,:) = 0.0_RP
       vz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_horizontalize_vec',2)

    return
  end subroutine OPRT_horizontalize_vec

  !-----------------------------------------------------------------------------
  subroutine OPRT_output_coef( &
       basename )
    use mod_process, only: &
       PRC_MPIstop
    use mod_io_param, only: &
       IO_REAL8, &
       IO_REAL4
    use mod_adm, only: &
       ADM_have_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_fio, only: &
       FIO_output
    use mod_hio, only: &
       HIO_output
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_MID) :: desc = 'Coefficients info'

    real(RP) :: tmp   (ADM_gall   ,106,ADM_lall   ,1)
    real(RP) :: tmp_pl(ADM_gall_pl,106,ADM_lall_pl,1)

    integer :: dtype
    integer :: g, l
    !---------------------------------------------------------------------------

    if    ( RP == SP ) then
       dtype = IO_REAL4
    elseif( RP == DP ) then
       dtype = IO_REAL8
    endif

    do l = 1, ADM_lall
       do g = 1, ADM_gall
          tmp(g,  1,l,1) = abs(OPRT_coef_div (1,g,0,l))
          tmp(g,  2,l,1) = abs(OPRT_coef_div (2,g,0,l))
          tmp(g,  3,l,1) = abs(OPRT_coef_div (3,g,0,l))
          tmp(g,  4,l,1) = abs(OPRT_coef_div (1,g,1,l))
          tmp(g,  5,l,1) = abs(OPRT_coef_div (2,g,1,l))
          tmp(g,  6,l,1) = abs(OPRT_coef_div (3,g,1,l))
          tmp(g,  7,l,1) = abs(OPRT_coef_div (1,g,2,l))
          tmp(g,  8,l,1) = abs(OPRT_coef_div (2,g,2,l))
          tmp(g,  9,l,1) = abs(OPRT_coef_div (3,g,2,l))
          tmp(g, 10,l,1) = abs(OPRT_coef_div (1,g,3,l))
          tmp(g, 11,l,1) = abs(OPRT_coef_div (2,g,3,l))
          tmp(g, 12,l,1) = abs(OPRT_coef_div (3,g,3,l))
          tmp(g, 13,l,1) = abs(OPRT_coef_div (1,g,4,l))
          tmp(g, 14,l,1) = abs(OPRT_coef_div (2,g,4,l))
          tmp(g, 15,l,1) = abs(OPRT_coef_div (3,g,4,l))
          tmp(g, 16,l,1) = abs(OPRT_coef_div (1,g,5,l))
          tmp(g, 17,l,1) = abs(OPRT_coef_div (2,g,5,l))
          tmp(g, 18,l,1) = abs(OPRT_coef_div (3,g,5,l))
          tmp(g, 19,l,1) = abs(OPRT_coef_div (1,g,6,l))
          tmp(g, 20,l,1) = abs(OPRT_coef_div (2,g,6,l))
          tmp(g, 21,l,1) = abs(OPRT_coef_div (3,g,6,l))
          tmp(g, 22,l,1) = abs(OPRT_coef_rot (1,g,0,l))
          tmp(g, 23,l,1) = abs(OPRT_coef_rot (2,g,0,l))
          tmp(g, 24,l,1) = abs(OPRT_coef_rot (3,g,0,l))
          tmp(g, 25,l,1) = abs(OPRT_coef_rot (1,g,1,l))
          tmp(g, 26,l,1) = abs(OPRT_coef_rot (2,g,1,l))
          tmp(g, 27,l,1) = abs(OPRT_coef_rot (3,g,1,l))
          tmp(g, 28,l,1) = abs(OPRT_coef_rot (1,g,2,l))
          tmp(g, 29,l,1) = abs(OPRT_coef_rot (2,g,2,l))
          tmp(g, 30,l,1) = abs(OPRT_coef_rot (3,g,2,l))
          tmp(g, 31,l,1) = abs(OPRT_coef_rot (1,g,3,l))
          tmp(g, 32,l,1) = abs(OPRT_coef_rot (2,g,3,l))
          tmp(g, 33,l,1) = abs(OPRT_coef_rot (3,g,3,l))
          tmp(g, 34,l,1) = abs(OPRT_coef_rot (1,g,4,l))
          tmp(g, 35,l,1) = abs(OPRT_coef_rot (2,g,4,l))
          tmp(g, 36,l,1) = abs(OPRT_coef_rot (3,g,4,l))
          tmp(g, 37,l,1) = abs(OPRT_coef_rot (1,g,5,l))
          tmp(g, 38,l,1) = abs(OPRT_coef_rot (2,g,5,l))
          tmp(g, 39,l,1) = abs(OPRT_coef_rot (3,g,5,l))
          tmp(g, 40,l,1) = abs(OPRT_coef_rot (1,g,6,l))
          tmp(g, 41,l,1) = abs(OPRT_coef_rot (2,g,6,l))
          tmp(g, 42,l,1) = abs(OPRT_coef_rot (3,g,6,l))
          tmp(g, 43,l,1) = abs(OPRT_coef_grad(1,g,0,l))
          tmp(g, 44,l,1) = abs(OPRT_coef_grad(2,g,0,l))
          tmp(g, 45,l,1) = abs(OPRT_coef_grad(3,g,0,l))
          tmp(g, 46,l,1) = abs(OPRT_coef_grad(1,g,1,l))
          tmp(g, 47,l,1) = abs(OPRT_coef_grad(2,g,1,l))
          tmp(g, 48,l,1) = abs(OPRT_coef_grad(3,g,1,l))
          tmp(g, 49,l,1) = abs(OPRT_coef_grad(1,g,2,l))
          tmp(g, 50,l,1) = abs(OPRT_coef_grad(2,g,2,l))
          tmp(g, 51,l,1) = abs(OPRT_coef_grad(3,g,2,l))
          tmp(g, 52,l,1) = abs(OPRT_coef_grad(1,g,3,l))
          tmp(g, 53,l,1) = abs(OPRT_coef_grad(2,g,3,l))
          tmp(g, 54,l,1) = abs(OPRT_coef_grad(3,g,3,l))
          tmp(g, 55,l,1) = abs(OPRT_coef_grad(1,g,4,l))
          tmp(g, 56,l,1) = abs(OPRT_coef_grad(2,g,4,l))
          tmp(g, 57,l,1) = abs(OPRT_coef_grad(3,g,4,l))
          tmp(g, 58,l,1) = abs(OPRT_coef_grad(1,g,5,l))
          tmp(g, 59,l,1) = abs(OPRT_coef_grad(2,g,5,l))
          tmp(g, 60,l,1) = abs(OPRT_coef_grad(3,g,5,l))
          tmp(g, 61,l,1) = abs(OPRT_coef_grad(1,g,6,l))
          tmp(g, 62,l,1) = abs(OPRT_coef_grad(2,g,6,l))
          tmp(g, 63,l,1) = abs(OPRT_coef_grad(3,g,6,l))
          tmp(g, 64,l,1) = abs(OPRT_coef_lap (  g,0,l))
          tmp(g, 65,l,1) = abs(OPRT_coef_lap (  g,1,l))
          tmp(g, 66,l,1) = abs(OPRT_coef_lap (  g,2,l))
          tmp(g, 67,l,1) = abs(OPRT_coef_lap (  g,3,l))
          tmp(g, 68,l,1) = abs(OPRT_coef_lap (  g,4,l))
          tmp(g, 69,l,1) = abs(OPRT_coef_lap (  g,5,l))
          tmp(g, 70,l,1) = abs(OPRT_coef_lap (  g,6,l))
          tmp(g, 71,l,1) = abs(OPRT_coef_intp(1,g,1,1,l))
          tmp(g, 72,l,1) = abs(OPRT_coef_intp(2,g,1,1,l))
          tmp(g, 73,l,1) = abs(OPRT_coef_intp(3,g,1,1,l))
          tmp(g, 74,l,1) = abs(OPRT_coef_intp(1,g,2,1,l))
          tmp(g, 75,l,1) = abs(OPRT_coef_intp(2,g,2,1,l))
          tmp(g, 76,l,1) = abs(OPRT_coef_intp(3,g,2,1,l))
          tmp(g, 77,l,1) = abs(OPRT_coef_intp(1,g,3,1,l))
          tmp(g, 78,l,1) = abs(OPRT_coef_intp(2,g,3,1,l))
          tmp(g, 79,l,1) = abs(OPRT_coef_intp(3,g,3,1,l))
          tmp(g, 80,l,1) = abs(OPRT_coef_intp(1,g,1,2,l))
          tmp(g, 81,l,1) = abs(OPRT_coef_intp(2,g,1,2,l))
          tmp(g, 82,l,1) = abs(OPRT_coef_intp(3,g,1,2,l))
          tmp(g, 83,l,1) = abs(OPRT_coef_intp(1,g,2,2,l))
          tmp(g, 84,l,1) = abs(OPRT_coef_intp(2,g,2,2,l))
          tmp(g, 85,l,1) = abs(OPRT_coef_intp(3,g,2,2,l))
          tmp(g, 86,l,1) = abs(OPRT_coef_intp(1,g,3,2,l))
          tmp(g, 87,l,1) = abs(OPRT_coef_intp(2,g,3,2,l))
          tmp(g, 88,l,1) = abs(OPRT_coef_intp(3,g,3,2,l))
          tmp(g, 89,l,1) = abs(OPRT_coef_diff(1,g,1,l))
          tmp(g, 90,l,1) = abs(OPRT_coef_diff(2,g,1,l))
          tmp(g, 91,l,1) = abs(OPRT_coef_diff(3,g,1,l))
          tmp(g, 92,l,1) = abs(OPRT_coef_diff(1,g,2,l))
          tmp(g, 93,l,1) = abs(OPRT_coef_diff(2,g,2,l))
          tmp(g, 94,l,1) = abs(OPRT_coef_diff(3,g,2,l))
          tmp(g, 95,l,1) = abs(OPRT_coef_diff(1,g,3,l))
          tmp(g, 96,l,1) = abs(OPRT_coef_diff(2,g,3,l))
          tmp(g, 97,l,1) = abs(OPRT_coef_diff(3,g,3,l))
          tmp(g, 98,l,1) = abs(OPRT_coef_diff(1,g,4,l))
          tmp(g, 99,l,1) = abs(OPRT_coef_diff(2,g,4,l))
          tmp(g,100,l,1) = abs(OPRT_coef_diff(3,g,4,l))
          tmp(g,101,l,1) = abs(OPRT_coef_diff(1,g,5,l))
          tmp(g,102,l,1) = abs(OPRT_coef_diff(2,g,5,l))
          tmp(g,103,l,1) = abs(OPRT_coef_diff(3,g,5,l))
          tmp(g,104,l,1) = abs(OPRT_coef_diff(1,g,6,l))
          tmp(g,105,l,1) = abs(OPRT_coef_diff(2,g,6,l))
          tmp(g,106,l,1) = abs(OPRT_coef_diff(3,g,6,l))
       enddo
    enddo

    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          tmp_pl(g,  1,l,1) = abs(OPRT_coef_div_pl (1,0,l))
          tmp_pl(g,  2,l,1) = abs(OPRT_coef_div_pl (2,0,l))
          tmp_pl(g,  3,l,1) = abs(OPRT_coef_div_pl (3,0,l))
          tmp_pl(g,  4,l,1) = abs(OPRT_coef_div_pl (1,1,l))
          tmp_pl(g,  5,l,1) = abs(OPRT_coef_div_pl (2,1,l))
          tmp_pl(g,  6,l,1) = abs(OPRT_coef_div_pl (3,1,l))
          tmp_pl(g,  7,l,1) = abs(OPRT_coef_div_pl (1,2,l))
          tmp_pl(g,  8,l,1) = abs(OPRT_coef_div_pl (2,2,l))
          tmp_pl(g,  9,l,1) = abs(OPRT_coef_div_pl (3,2,l))
          tmp_pl(g, 10,l,1) = abs(OPRT_coef_div_pl (1,3,l))
          tmp_pl(g, 11,l,1) = abs(OPRT_coef_div_pl (2,3,l))
          tmp_pl(g, 12,l,1) = abs(OPRT_coef_div_pl (3,3,l))
          tmp_pl(g, 13,l,1) = abs(OPRT_coef_div_pl (1,4,l))
          tmp_pl(g, 14,l,1) = abs(OPRT_coef_div_pl (2,4,l))
          tmp_pl(g, 15,l,1) = abs(OPRT_coef_div_pl (3,4,l))
          tmp_pl(g, 16,l,1) = abs(OPRT_coef_div_pl (1,5,l))
          tmp_pl(g, 17,l,1) = abs(OPRT_coef_div_pl (2,5,l))
          tmp_pl(g, 18,l,1) = abs(OPRT_coef_div_pl (3,5,l))
          tmp_pl(g, 19,l,1) = abs(OPRT_coef_div_pl (1,1,l))
          tmp_pl(g, 20,l,1) = abs(OPRT_coef_div_pl (2,1,l))
          tmp_pl(g, 21,l,1) = abs(OPRT_coef_div_pl (3,1,l))
          tmp_pl(g, 22,l,1) = abs(OPRT_coef_rot_pl (1,0,l))
          tmp_pl(g, 23,l,1) = abs(OPRT_coef_rot_pl (2,0,l))
          tmp_pl(g, 24,l,1) = abs(OPRT_coef_rot_pl (3,0,l))
          tmp_pl(g, 25,l,1) = abs(OPRT_coef_rot_pl (1,1,l))
          tmp_pl(g, 26,l,1) = abs(OPRT_coef_rot_pl (2,1,l))
          tmp_pl(g, 27,l,1) = abs(OPRT_coef_rot_pl (3,1,l))
          tmp_pl(g, 28,l,1) = abs(OPRT_coef_rot_pl (1,2,l))
          tmp_pl(g, 29,l,1) = abs(OPRT_coef_rot_pl (2,2,l))
          tmp_pl(g, 30,l,1) = abs(OPRT_coef_rot_pl (3,2,l))
          tmp_pl(g, 31,l,1) = abs(OPRT_coef_rot_pl (1,3,l))
          tmp_pl(g, 32,l,1) = abs(OPRT_coef_rot_pl (2,3,l))
          tmp_pl(g, 33,l,1) = abs(OPRT_coef_rot_pl (3,3,l))
          tmp_pl(g, 34,l,1) = abs(OPRT_coef_rot_pl (1,4,l))
          tmp_pl(g, 35,l,1) = abs(OPRT_coef_rot_pl (2,4,l))
          tmp_pl(g, 36,l,1) = abs(OPRT_coef_rot_pl (3,4,l))
          tmp_pl(g, 37,l,1) = abs(OPRT_coef_rot_pl (1,5,l))
          tmp_pl(g, 38,l,1) = abs(OPRT_coef_rot_pl (2,5,l))
          tmp_pl(g, 39,l,1) = abs(OPRT_coef_rot_pl (3,5,l))
          tmp_pl(g, 40,l,1) = abs(OPRT_coef_rot_pl (1,1,l))
          tmp_pl(g, 41,l,1) = abs(OPRT_coef_rot_pl (2,1,l))
          tmp_pl(g, 42,l,1) = abs(OPRT_coef_rot_pl (3,1,l))
          tmp_pl(g, 43,l,1) = abs(OPRT_coef_grad_pl(1,0,l))
          tmp_pl(g, 44,l,1) = abs(OPRT_coef_grad_pl(2,0,l))
          tmp_pl(g, 45,l,1) = abs(OPRT_coef_grad_pl(3,0,l))
          tmp_pl(g, 46,l,1) = abs(OPRT_coef_grad_pl(1,1,l))
          tmp_pl(g, 47,l,1) = abs(OPRT_coef_grad_pl(2,1,l))
          tmp_pl(g, 48,l,1) = abs(OPRT_coef_grad_pl(3,1,l))
          tmp_pl(g, 49,l,1) = abs(OPRT_coef_grad_pl(1,2,l))
          tmp_pl(g, 50,l,1) = abs(OPRT_coef_grad_pl(2,2,l))
          tmp_pl(g, 51,l,1) = abs(OPRT_coef_grad_pl(3,2,l))
          tmp_pl(g, 52,l,1) = abs(OPRT_coef_grad_pl(1,3,l))
          tmp_pl(g, 53,l,1) = abs(OPRT_coef_grad_pl(2,3,l))
          tmp_pl(g, 54,l,1) = abs(OPRT_coef_grad_pl(3,3,l))
          tmp_pl(g, 55,l,1) = abs(OPRT_coef_grad_pl(1,4,l))
          tmp_pl(g, 56,l,1) = abs(OPRT_coef_grad_pl(2,4,l))
          tmp_pl(g, 57,l,1) = abs(OPRT_coef_grad_pl(3,4,l))
          tmp_pl(g, 58,l,1) = abs(OPRT_coef_grad_pl(1,5,l))
          tmp_pl(g, 59,l,1) = abs(OPRT_coef_grad_pl(2,5,l))
          tmp_pl(g, 60,l,1) = abs(OPRT_coef_grad_pl(3,5,l))
          tmp_pl(g, 61,l,1) = abs(OPRT_coef_grad_pl(1,1,l))
          tmp_pl(g, 62,l,1) = abs(OPRT_coef_grad_pl(2,1,l))
          tmp_pl(g, 63,l,1) = abs(OPRT_coef_grad_pl(3,1,l))
          tmp_pl(g, 64,l,1) = abs(OPRT_coef_lap_pl (  0,l))
          tmp_pl(g, 65,l,1) = abs(OPRT_coef_lap_pl (  1,l))
          tmp_pl(g, 66,l,1) = abs(OPRT_coef_lap_pl (  2,l))
          tmp_pl(g, 67,l,1) = abs(OPRT_coef_lap_pl (  3,l))
          tmp_pl(g, 68,l,1) = abs(OPRT_coef_lap_pl (  4,l))
          tmp_pl(g, 69,l,1) = abs(OPRT_coef_lap_pl (  5,l))
          tmp_pl(g, 70,l,1) = abs(OPRT_coef_lap_pl (  1,l))
          tmp_pl(g, 71,l,1) = abs(OPRT_coef_intp_pl(1,g,1,l))
          tmp_pl(g, 72,l,1) = abs(OPRT_coef_intp_pl(2,g,1,l))
          tmp_pl(g, 73,l,1) = abs(OPRT_coef_intp_pl(3,g,1,l))
          tmp_pl(g, 74,l,1) = abs(OPRT_coef_intp_pl(1,g,2,l))
          tmp_pl(g, 75,l,1) = abs(OPRT_coef_intp_pl(2,g,2,l))
          tmp_pl(g, 76,l,1) = abs(OPRT_coef_intp_pl(3,g,2,l))
          tmp_pl(g, 77,l,1) = abs(OPRT_coef_intp_pl(1,g,3,l))
          tmp_pl(g, 78,l,1) = abs(OPRT_coef_intp_pl(2,g,3,l))
          tmp_pl(g, 79,l,1) = abs(OPRT_coef_intp_pl(3,g,3,l))
          tmp_pl(g, 80,l,1) = abs(OPRT_coef_intp_pl(1,g,1,l))
          tmp_pl(g, 81,l,1) = abs(OPRT_coef_intp_pl(2,g,1,l))
          tmp_pl(g, 82,l,1) = abs(OPRT_coef_intp_pl(3,g,1,l))
          tmp_pl(g, 83,l,1) = abs(OPRT_coef_intp_pl(1,g,2,l))
          tmp_pl(g, 84,l,1) = abs(OPRT_coef_intp_pl(2,g,2,l))
          tmp_pl(g, 85,l,1) = abs(OPRT_coef_intp_pl(3,g,2,l))
          tmp_pl(g, 86,l,1) = abs(OPRT_coef_intp_pl(1,g,3,l))
          tmp_pl(g, 87,l,1) = abs(OPRT_coef_intp_pl(2,g,3,l))
          tmp_pl(g, 88,l,1) = abs(OPRT_coef_intp_pl(3,g,3,l))
          tmp_pl(g, 89,l,1) = abs(OPRT_coef_diff_pl(1,1,l))
          tmp_pl(g, 90,l,1) = abs(OPRT_coef_diff_pl(2,1,l))
          tmp_pl(g, 91,l,1) = abs(OPRT_coef_diff_pl(3,1,l))
          tmp_pl(g, 92,l,1) = abs(OPRT_coef_diff_pl(1,2,l))
          tmp_pl(g, 93,l,1) = abs(OPRT_coef_diff_pl(2,2,l))
          tmp_pl(g, 94,l,1) = abs(OPRT_coef_diff_pl(3,2,l))
          tmp_pl(g, 95,l,1) = abs(OPRT_coef_diff_pl(1,3,l))
          tmp_pl(g, 96,l,1) = abs(OPRT_coef_diff_pl(2,3,l))
          tmp_pl(g, 97,l,1) = abs(OPRT_coef_diff_pl(3,3,l))
          tmp_pl(g, 98,l,1) = abs(OPRT_coef_diff_pl(1,4,l))
          tmp_pl(g, 99,l,1) = abs(OPRT_coef_diff_pl(2,4,l))
          tmp_pl(g,100,l,1) = abs(OPRT_coef_diff_pl(3,4,l))
          tmp_pl(g,101,l,1) = abs(OPRT_coef_diff_pl(1,5,l))
          tmp_pl(g,102,l,1) = abs(OPRT_coef_diff_pl(2,5,l))
          tmp_pl(g,103,l,1) = abs(OPRT_coef_diff_pl(3,5,l))
          tmp_pl(g,104,l,1) = abs(OPRT_coef_diff_pl(1,1,l))
          tmp_pl(g,105,l,1) = abs(OPRT_coef_diff_pl(2,1,l))
          tmp_pl(g,106,l,1) = abs(OPRT_coef_diff_pl(3,1,l))
       enddo
       enddo
    endif

    call COMM_data_transfer( tmp, tmp_pl )

    if ( OPRT_io_mode == 'ADVANCED' ) then

       call FIO_output( tmp(:,:,:,1), basename, desc, '',               & ! [IN]
                        'oprtcoef', 'oprt coef', '',                    & ! [IN]
                        '', dtype, 'LAYERNM', 1, 106, 1, 0.0_DP, 0.0_DP ) ! [IN]

    elseif( OPRT_io_mode == 'POH5' ) then

       call HIO_output( tmp(:,:,:,1), basename, desc, '',               & ! [IN]
                        'oprtcoef', 'oprt coef', '',                    & ! [IN]
                        '', dtype, 'LAYERNM', 1, 106, 1, 0.0_DP, 0.0_DP ) ! [IN]

    else
       write(IO_FID_LOG,*) 'Invalid io_mode!'
       call PRC_MPIstop
    endif

    return
  end subroutine OPRT_output_coef

end module mod_oprt
