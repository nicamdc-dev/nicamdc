!-------------------------------------------------------------------------------
!
!+  Numerical smoothing module
!
!-------------------------------------------------------------------------------
module mod_numfilter
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for numerical smoothings or filters.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                05-11-02   add 2 dimensional damping (N. Hirota)
  !                05-11-10   add z dependent horizontal diffusion (N. Hirota)
  !                05-12-09   M.Satoh bug fix
  !                05-12-10   1/gamma -> tau when type==E_FOLD_TIME (N. Hirota)
  !                05-12-17   M.Satoh: read namelist twice for compatibility
  !                06-01-10   S.Iga: add 'Kh_coef_lap1' 
  !                06-04-17   H.Tomita : Add nonlinear diffusion
  !                06-10-20   K.Suzuki : not calling tracer diffusion in using
  !                                      Miura(2004) advection scheme
  !                07-01-26   H.Tomita : Add an option [rayleigh_damp_only_w].
  !                                      Some change in 
  !                                      sub[numfilter_rayleigh_damping].
  !                07-08-07   T.Mitsui : trivial fix
  !                08-01-24   Y.Niwa : add MIURA2004OLD in numerical_hdiff
  !                                    hdiff_fact_q = 0.D0
  !                08-01-30   Y.Niwa : bug fix
  !                08-04-12   T.Mitsui: omit needless calculation in hdiff
  !                09-04-14   H.Tomita: Add the initilization ( zero clear )
  !                                     for Kh_coef_lap1.
  !                11-11-28   Y.Yamada: Merge Terai-san code
  !                                        into the original code.
  !                12-02-10   T.Yamaura: Optimized numfilter_divdamp and
  !                                        numfilter_numerical_hdiff for K-Computer (comitted by Iga 12-03-09)
  !                12-03-28   T.Seiki : fix undefined reference,
  !                                     change the method to calculate horiz_dx2
  !                                     collective communication => analytic diagnosis
  !                12-05-30   T.Yashiro: Change arguments from character to index/switch
  !                12-07-21   S.Iga: add switch 'logical:smooth_1var' which control 
  !                                 'call smooth_1var'
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  logical, public, save :: NUMFILTER_DOrayleigh     = .false.
  logical, public, save :: NUMFILTER_DOverticaldiff = .false.
  logical, public, save :: NUMFILTER_DOdivdamp      = .false.
  logical, public, save :: NUMFILTER_DOdivdamp_2d   = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: numfilter_setup
  private :: numfilter_rayleigh_damping_setup
! private :: numfilter_horizontal_diffusion_setup
  private :: numfilter_vertical_diffusion_setup
  private :: numfilter_divdamp_setup
  private :: numfilter_divdamp_2d_setup

  public :: numfilter_rayleigh_damping
!  public :: numfilter_horizontal_diffusion
  public :: numfilter_numerical_hdiff
  public :: numfilter_vertical_diffusion
  public :: numfilter_divdamp
  public :: numfilter_divdamp_2d

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  real(8), private, allocatable, save :: rayleigh_coef  (:) ! Rayleigh damping coefficient at cell center
  real(8), private, allocatable, save :: rayleigh_coef_h(:) ! Rayleigh damping coefficient at cell wall
  logical, private, save :: rayleigh_damp_only_w = .false. 

  real(8), private, allocatable, save :: Kv_coef  (:)       ! vertical diffusion coefficient at cell center
  real(8), private, allocatable, save :: Kv_coef_h(:)       ! vertical diffusion coefficient at cell wall

  real(8), private, allocatable, save :: divdamp_coef   (:,:,:) ! divergence damping coefficient at cell center
  real(8), private, allocatable, save :: divdamp_coef_pl(:,:,:)
  integer, private, save :: lap_order_divdamp = 1               ! laplacian order
  real(8), private, save :: divdamp_coef_v

  real(8), private, allocatable, save :: divdamp_2d_coef   (:,:,:) ! divergence damping coefficient at cell center
  real(8), private, allocatable, save :: divdamp_2d_coef_pl(:,:,:)
  integer, private, save :: lap_order_divdamp_2d = 1               ! laplacian order

  !------ Dependency of horizontal grid size 
  logical, private, save :: dep_hgrid = .false.

  !------ Horizontal diffusion type 
  character(ADM_NSYS), private, save :: hdiff_type = 'NONDIM_COEF'
  !                                                  'DIRECT'
  !                                                  'E_FOLD_TIME'
  !
  !------ coefficient for diffusion (horizontal)
  integer, private, save :: lap_order_hdiff = 2
  real(8), private, save :: gamma_h = 1.0D0/16.D0/10.D0
  !<-----        Durran and Klemp(1983) 0.015 horizontal
  real(8), private, save :: gamma_h_zdep = 0.D0
  real(8), private, save :: tau_h = 160000.d0   ! N. Hirota 051210
  real(8), private, save :: tau_h_top = 40000.d0   ! N. Hirota 051210
  real(8), private, save :: ZD_hdiff     = 25000.D0 !(add 'save'  S.Iga060110)
  real(8), private, save :: hdiff_fact_rho = 1.0D0
  real(8), private, save :: hdiff_fact_q = 0.D0  ! Y.Niwa 080124 1.0D0 => 0.D0
  !------ coefficient for diffusion (horizontal) with laplacian order 1 
  !  (independent from gamma_h.  e-fold time is not available)  S.Iga060110
  real(8), private, save :: gamma_h_lap1 = 0.D0
  real(8), private, save :: ZD_hdiff_lap1     = 25000.D0
  !
  !------ hight for decay of nonlinear diffusion
  real(8), private, save :: ZD_hdiff_nl = 20000.D0

  !<-----        Durran and Klemp(1983) 0.001 vertical
  !
  !------ max coefficient for diffusion (horizontal)
  real(8), private, save :: gamma_hmax = 1.0D0/16.D0/10.D0
  !<-----        Durran and Klemp(1983) 1.d0/16.d0
  !
  !------ horizontal diffusion coefficient
  real(8), private, allocatable, save :: Kh_coef(:,:,:)   !<--- at cell center
  real(8), private, allocatable, save :: Kh_coef_h(:,:,:) !<--- at cell wall
  real(8), private, allocatable, save :: Kh_coef_pl(:,:,:)   !<--- at cell center
  real(8), private, allocatable, save :: Kh_coef_h_pl(:,:,:) !<--- at cell wall
  real(8), private, allocatable, save :: Kh_coef_zdep(:,:,:)   !<--- at cell center
  real(8), private, allocatable, save :: Kh_coef_zdep_pl(:,:,:)   !<--- at cell center
  !------ horizontal diffusion coefficient for laplacian-order 1 
  !------ (additional to Kh_coef) (S.Iga 0601010)
  real(8), private, allocatable, save :: Kh_coef_lap1(:,:,:)   !<--- at cell center
  real(8), private, allocatable, save :: Kh_coef_lap1_h(:,:,:) !<--- at cell wall
  real(8), private, allocatable, save :: Kh_coef_lap1_pl(:,:,:)   !<--- at cell center
  real(8), private, allocatable, save :: Kh_coef_lap1_h_pl(:,:,:) !<--- at cell wall
  !
  !

  real(8), private, save :: Kh_coef_minlim = 3.0D+9
  real(8), private, save :: Kh_coef_maxlim = 1.0D+99
  real(8), private, save :: cfact = 2.0D0

  logical, private, save :: DEEP_EFFECT = .true.

  real(8), private, save :: horiz_dx2

  logical, private, save :: debug_output = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine numfilter_setup
    use mod_adm, only :  &
         ADM_CTL_FID,    &
         ADM_LOG_FID,    &
         ADM_proc_stop, &
         ADM_GLEVEL,     & ! [Add] 12/03/28 T.Seiki
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_prc_me,     &
         ADM_prc_pl,     &
         ADM_kmax,       &
         ADM_KNONE
    use mod_cnst, only : &
         CNST_PI,        &
         CNST_RAIR,      &
         CNST_GAMMA,     &
         CNST_ERADIUS
    use mod_grd, only :  &
         GRD_gz,         &
         GRD_dgz,        &
         GRD_gzh,        &
         GRD_dgzh
    use mod_gmtr, only: &
         GMTR_area,    &
         GMTR_area_pl
    use mod_time, only: &
         TIME_DTS, &
         TIME_DTL
    use mod_gtl, only: &
         GTL_min_k
    implicit none

    ! rayleigh damping
    real(8) :: alpha_r = 1.D0 / ( 3600D0*3.D0 ) / 10.D0 ! coefficient for rayleigh damping
    real(8) :: ZD      = 5000.D0                        ! lower limit of rayleigh damping [m]

    ! vertical diffusion
    real(8) :: gamma_v = 0.001D0 ! coefficient of vertical diffusion

    ! 3D divergence damping
    character(len=ADM_NSYS) :: divdamp_type      = 'NONDIM_COEF' ! damping type
    real(8)                 :: alpha_d           = 1.D0          ! coefficient    for divergence damping
    real(8)                 :: tau_d             = 132800.D0     ! e-folding time for divergence damping
    real(8)                 :: alpha_dv          = 0.D0          ! vertical coefficient

    ! 2D divergence damping
    character(len=ADM_NSYS) :: divdamp_2d_type      = 'NONDIM_COEF' ! damping type
    real(8)                 :: alpha_d_2d           = 0.D0          ! coefficient    for divergence damping
    real(8)                 :: tau_d_2d             = 1328000.D0    ! e-folding time for divergence damping
    real(8)                 :: ZD_d_2d              = 25000.D0      ! lower limit of divergence damping [m]

    logical,save :: smooth_1var = .true. !(S.Iga 20120721.  It sould be false for stretched grid)
                                  
    namelist / NUMFILTERPARAM / & ! add tau (N. Hirota 051210)
         debug_output,         &
         alpha_r,              &
         ZD,                   &
         rayleigh_damp_only_w, & !--- rayleigh damping for only w?
         gamma_v,              &
         divdamp_type,         &
         lap_order_divdamp,    &
         alpha_d,              &
         tau_d,                &
         alpha_dv,             &
         divdamp_2d_type,      &
         lap_order_divdamp_2d, &
         alpha_d_2d,           &
         tau_d_2d,             &
         ZD_d_2d,              &
         dep_hgrid,            & !--- depnedency of horzontal grid size
         hdiff_type,           & !--- horizontal diffusion type
         lap_order_hdiff,      & !--- laplacian order of horz. diff.
         hdiff_fact_rho,       & !--- factor of hdiff for density
         hdiff_fact_q,         & !--- factor of hdiff for tracer
         gamma_h,              & !--- coefficient of horizontal diffusion.
         gamma_h_zdep,         & !--- coefficient of z dep h diffusion.
         tau_h,                & !--- E_FOLD_TIME of h diff.
         tau_h_top,            & !--- E_FOLD_TIME of h diff at top.
         ZD_hdiff,             & !--- at z>ZD, gamma_h_zdep (tau_h_top) is effective.
         gamma_h_lap1,         & !--- coefficient of 1-order laplacian type !  h-diffusion.  additional to gamma_h
         ZD_hdiff_lap1,        & !--- at z>ZD, gamma_h_zdep (tau_h_top) is effective.            
         Kh_coef_minlim,       & !--- Kh_coef (minimum limit)
         Kh_coef_maxlim,       & !--- Kh_coef (maximum limit)
         cfact,                & !--- coefficient for nonlinear diffusion
         ZD_hdiff_nl,          &
         DEEP_EFFECT,          &
         smooth_1var             !--- smooth Kh and divdamp coef. (default true. In fact, I recommend false)

    real(8) :: global_area, global_grid

    integer :: ierr
    integer :: k, l, n
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[numfilter]/Category[nhm dynamics]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=NUMFILTERPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** NUMFILTERPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist NUMFILTERPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist NUMFILTERPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,NUMFILTERPARAM)

    if ( hdiff_type == 'E_FOLD_TIME' ) then
       gamma_h      = -999.D0
       gamma_h_zdep = -999.D0
    endif

    if ( gamma_h > 0.D0 ) then
       tau_h     = gamma_h
       tau_h_top = gamma_h
    endif

    if ( gamma_h_zdep > 0.D0 ) tau_h_top = gamma_h_zdep

    if (      divdamp_type /= "DIRECT" &
         .OR. hdiff_type   /= "DIRECT" ) then
       global_area = 4.D0 * CNST_PI * CNST_ERADIUS * CNST_ERADIUS
       global_grid = 10.D0 * 4.D0**ADM_GLEVEL
       horiz_dx2   = global_area / global_grid
    endif



    call numfilter_rayleigh_damping_setup( alpha_r, & ! [IN]
                                           ZD       ) ! [IN]

    call numfilter_vertical_diffusion_setup( gamma_v ) ! [IN]

    call numfilter_divdamp_setup( divdamp_type,      & ! [IN]
                                  dep_hgrid,         & ! [IN]
                                  lap_order_divdamp, & ! [IN]
                                  alpha_d,           & ! [IN]
                                  tau_d,             & ! [IN]
                                  alpha_dv           ) ! [IN]

    call numfilter_divdamp_2d_setup( divdamp_2d_type,      & ! [IN]
                                     dep_hgrid,            & ! [IN]
                                     lap_order_divdamp_2d, & ! [IN]
                                     alpha_d_2d,           & ! [IN]
                                     tau_d_2d,             & ! [IN]
                                     ZD_d_2d               ) ! [IN]

!    call numfilter_divdamp_2d_setup

    !--- coefficient for numerical diffusion

    allocate( Kh_coef        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Kh_coef_h      (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_h_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Kh_coef_zdep   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_zdep_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( Kh_coef_lap1     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_lap1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Kh_coef_lap1_h   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_lap1_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    Kh_coef_lap1      = 0.D0
    Kh_coef_lap1_pl   = 0.D0
    Kh_coef_lap1_h    = 0.D0
    Kh_coef_lap1_h_pl = 0.D0

    do l = 1, ADM_lall
       do k = 1, ADM_kall
          ! if hdiff_type = 'DIRECT'
          !   k=gamma
          !   k_zdep=gamma_zdep
          ! if hdiff_type = 'nondim_coef'
          !   k=gamma_h
          !     *(GMTR_area(:,l)**lap_order_hdiff)/TIME_DTL)
          !   k_zdep=gamma_h_zdep
          !     *(GMTR_area(:,l)**lap_order_hdiff)/TIME_DTL)
          ! if hdiff_type = 'E_FOLD_TIME'
          !   k=1/tau_h
          !   k_zdep=1/tau_h_top-1/tau_h
          !   
          !   k=k                                          for z<zd
          !   k=k+k_zdep*[1-cos(pi(z-zd)/(ztop-zd))]/2}    for z>=zd
          !         
          if (trim(hdiff_type)=="DIRECT") then
             !--- in this case, gamma_h is an absolute value.
             do n = 1, ADM_gall
                Kh_coef(n,k,l)   = gamma_h
                Kh_coef_h(n,k,l) = gamma_h
                Kh_coef_zdep(n,k,l)   = gamma_h_zdep
             enddo
          else if (trim(hdiff_type)=="NONDIM_COEF") then
             !--- in this case, gamma_h is a non-dimensional number.
             if(dep_hgrid) then
                do n = 1, ADM_gall
                   Kh_coef(n,k,l) = gamma_h/TIME_DTL&
                        *(GMTR_area(n,l)**lap_order_hdiff)
                   Kh_coef_h(n,k,l) = Kh_coef(n,k,l)
                   Kh_coef_zdep(n,k,l) = gamma_h_zdep/TIME_DTL&
                        *(GMTR_area(n,l)**lap_order_hdiff)
                enddo
             else
                do n = 1, ADM_gall
                   Kh_coef(n,k,l) = gamma_h/TIME_DTL&
                        *(horiz_dx2**lap_order_hdiff)
                   Kh_coef_h(n,k,l) = Kh_coef(n,k,l)
                   Kh_coef_zdep(n,k,l) = gamma_h_zdep/TIME_DTL&
                        *(horiz_dx2**lap_order_hdiff)
                enddo
             end if
          else if (trim(hdiff_type)=="E_FOLD_TIME") then
             !--- in this case, tau_h is e-folding time for 2DX waves.
             if(dep_hgrid) then
                do n = 1, ADM_gall
                   Kh_coef(n,k,l) &
                        = ((sqrt(GMTR_area(n,l))/CNST_PI)&
                        **(2*lap_order_hdiff))&
                        *(1.d0/tau_h)   ! 051210 N. Hirota
                   Kh_coef_h(n,k,l) = Kh_coef(n,k,l)
                   Kh_coef_zdep(n,k,l) &
                        = ((sqrt(GMTR_area(n,l))/CNST_PI)&
                        **(2*lap_order_hdiff))&
                        *(1.d0/tau_h_top-1.d0/tau_h)   ! 051210 N. Hirota
!                   Kh_coef_zdep(n,k,l) = 0.d0 ! 05/12/09 M.Satoh
                enddo
             else
                do n = 1, ADM_gall
                   Kh_coef(n,k,l)  = ((sqrt(horiz_dx2)/CNST_PI)&
                        **(2*lap_order_hdiff))&
                        *(1.d0/tau_h)   ! 051210 N. Hirota
                   Kh_coef_h(n,k,l) = Kh_coef(n,k,l)
                   Kh_coef_zdep(n,k,l)  = ((sqrt(horiz_dx2)/CNST_PI)&
                        **(2*lap_order_hdiff))&
                        *(1.d0/tau_h_top-1.d0/tau_h)   ! 051210 N. Hirota
!                   Kh_coef_zdep(n,k,l) = 0.d0 ! 05/12/09 M.Satoh
                enddo
             end if
          else !!! 07.08.07 [fix] T.Mitsui for udefined value
             Kh_coef(:,k,l) = 0.d0
             Kh_coef_h(:,k,l) = 0.d0
             Kh_coef_zdep(:,k,l) = 0.d0
          end if
          if ( GRD_gz(k) < ZD_hdiff ) then
             Kh_coef_zdep(:,k,l) = 0.D0
          else
             Kh_coef_zdep(:,k,l) = 0.5D0 * Kh_coef_zdep(:,k,l)  &
                  * ( 1.0D0 - cos ( CNST_PI                           &
                  * ( GRD_gz(k) - ZD_hdiff )                          &
                  / ( GRD_gzh(ADM_kmax+1) - ZD_hdiff ) ) )
          end if
          Kh_coef(:,k,l)=Kh_coef(:,k,l)+Kh_coef_zdep(:,k,l)
          Kh_coef_h(:,k,l) = Kh_coef(:,k,l)
          !---- Kh-coef for order 1 laplacian (S.Iga060110) =>
                  !---  gamma_h_lap1 is an absolute value.
          do n = 1, ADM_gall
             Kh_coef_lap1(n,k,l)   = gamma_h_lap1
!             Kh_coef_lap1_h(n,k,l) = gamma_h_lap1
          enddo
          if ( GRD_gz(k) < ZD_hdiff_lap1 ) then
             Kh_coef_lap1(:,k,l) = 0.D0
          else
             Kh_coef_lap1(:,k,l) = 0.5D0 * Kh_coef_lap1(:,k,l)  &
                  * ( 1.0D0 - cos ( CNST_PI                           &
                  * ( GRD_gz(k) - ZD_hdiff_lap1 )                          &
                  / ( GRD_gzh(ADM_kmax+1) - ZD_hdiff_lap1 ) ) )
          end if
          Kh_coef_lap1_h(:,k,l) = Kh_coef_lap1(:,k,l)
          !----                               (S.Iga060110) >=
          !
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             !
             !--- p-th order horizontal numerical diffusion
             if (trim(hdiff_type)=="DIRECT") then
                !--- in this case, gamma_h is an absolute value.
                do n = 1, ADM_gall_pl
                   Kh_coef_pl(n,k,l)   = gamma_h
                   Kh_coef_h_pl(n,k,l) = gamma_h
                   Kh_coef_zdep_pl(n,k,l) = gamma_h_zdep
                enddo
             else if (trim(hdiff_type)=="NONDIM_COEF") then
                !--- in this case, gamma_h is a non-dimensional number.
                if(dep_hgrid) then
                   do n = 1, ADM_gall_pl
                      Kh_coef_pl(n,k,l) = gamma_h/TIME_DTL&
                           *(GMTR_area_pl(n,l)**lap_order_hdiff)
                      Kh_coef_h_pl(n,k,l) = Kh_coef_pl(n,k,l)
                      Kh_coef_zdep_pl(n,k,l) = gamma_h_zdep/TIME_DTL&
                           *(GMTR_area_pl(n,l)**lap_order_hdiff)
                   enddo
                else
                   do n = 1, ADM_gall_pl
                      Kh_coef_pl(n,k,l) = gamma_h/TIME_DTL&
                           *(horiz_dx2**lap_order_hdiff)
                      Kh_coef_h_pl(n,k,l) = Kh_coef_pl(n,k,l)
                      Kh_coef_zdep_pl(n,k,l) = gamma_h_zdep/TIME_DTL&
                           *(horiz_dx2**lap_order_hdiff)
                   enddo
                end if
             else if (trim(hdiff_type)=="E_FOLD_TIME") then
                !--- in this case, gamma_h is e-folding time for 2DX waves.
                if(dep_hgrid) then
                   do n = 1, ADM_gall_pl
                      Kh_coef_pl(n,k,l) &
                           = ((sqrt(GMTR_area_pl(n,l))/CNST_PI)&
                           **(2*lap_order_hdiff))&
                           *(1.d0/tau_h)   ! 051210 N. Hirota
                      Kh_coef_h_pl(n,k,l) = Kh_coef_pl(n,k,l)
                      Kh_coef_zdep_pl(n,k,l) &
                           = ((sqrt(GMTR_area_pl(n,l))/CNST_PI)&
                           **(2*lap_order_hdiff))&
                           *(1.d0/tau_h_top-1.d0/tau_h)   ! 051210 N. Hirota
!                      Kh_coef_zdep_pl(n,k,l) = 0.d0 ! 05/12/09 M.Satoh
                   enddo
                else
                   do n = 1, ADM_gall_pl
                      Kh_coef_pl(n,k,l)  = ((sqrt(horiz_dx2)/CNST_PI)&
                           **(2*lap_order_hdiff))&
                           *(1.d0/tau_h)   ! 051210 N. Hirota
                      Kh_coef_h_pl(n,k,l) = Kh_coef_pl(n,k,l)
                      Kh_coef_zdep_pl(n,k,l)  = ((sqrt(horiz_dx2)/CNST_PI)&
                           **(2*lap_order_hdiff))&
                           *(1.d0/tau_h_top-1.d0/tau_h)   ! 051210 N. Hirota
!                      Kh_coef_zdep_pl(n,k,l) = 0.d0 ! 05/12/09 M.Satoh
                   enddo
                end if
             else !!! 07.08.07 [fix] T.Mitsui for udefined value
                Kh_coef_pl(:,k,l) = 0.d0
                Kh_coef_h_pl(:,k,l) = 0.d0
                Kh_coef_zdep_pl(:,k,l) = 0.d0
             end if
             if ( GRD_gz(k) < ZD_hdiff ) then
                Kh_coef_zdep_pl(:,k,l) = 0.D0
             else
                Kh_coef_zdep_pl(:,k,l) = 0.5D0 * Kh_coef_zdep_pl(:,k,l)  &
                     * ( 1.0D0 - cos ( CNST_PI                           &
                     * ( GRD_gz(k) - ZD_hdiff )                          &
                     / ( GRD_gzh(ADM_kmax+1) - ZD_hdiff ) ) )
             end if
             Kh_coef_pl(:,k,l)=Kh_coef_pl(:,k,l)+Kh_coef_zdep_pl(:,k,l)
             Kh_coef_h_pl(:,k,l) = Kh_coef_pl(:,k,l)

             !---- Kh-coef for order 1 laplacian (S.Iga060110) =>
                !---  gamma_h_lap1 is an absolute value.
                !--- in this case, gamma_h is an absolute value.
             do n = 1, ADM_gall_pl
                Kh_coef_lap1_pl(n,k,l)   = gamma_h_lap1
!                Kh_coef_lap1_h_pl(n,k,l) = gamma_h_lap1
             enddo
             if ( GRD_gz(k) < ZD_hdiff_lap1 ) then
                Kh_coef_lap1_pl(:,k,l) = 0.D0
             else
                Kh_coef_lap1_pl(:,k,l) = 0.5D0 * Kh_coef_lap1_pl(:,k,l)  &
                     * ( 1.0D0 - cos ( CNST_PI                           &
                     * ( GRD_gz(k) - ZD_hdiff_lap1 )                          &
                     / ( GRD_gzh(ADM_kmax+1) - ZD_hdiff_lap1 ) ) )
             end if
             Kh_coef_lap1_h_pl(:,k,l) = Kh_coef_lap1_pl(:,k,l)
             !----                               (S.Iga060110) >=

          enddo
       enddo
    end if


    if ( hdiff_type /= "DIRECT" ) then
       if (smooth_1var) then ! iga 20120721 (add if)
          call numfilter_smooth_1var(               Kh_coef, Kh_coef_pl                  )
          call numfilter_smooth_1var(               Kh_coef_h, Kh_coef_h_pl                 )
       endif
       Kh_coef  (:,:,:) = max( Kh_coef  (:,:,:), Kh_coef_minlim )
       Kh_coef_h(:,:,:) = max( Kh_coef_h(:,:,:), Kh_coef_minlim )
    endif

    if (trim(divdamp_type)/="DIRECT") then
       if (smooth_1var) then ! iga 20120721 (add if)
          call numfilter_smooth_1var(&
               divdamp_coef, divdamp_coef_pl   &  !--- [INOUT]
               )
       endif ! iga 20120721 (add if)
       divdamp_coef(:,:,:) = max(divdamp_coef(:,:,:),Kh_coef_minlim)
    end if

    if(DEEP_EFFECT) then
       do k=1, ADM_kall
          !
          Kh_coef(:,k,:) = Kh_coef(:,k,:)* ((GRD_gz(k)+CNST_ERADIUS)/CNST_ERADIUS)**(2*lap_order_hdiff)
          Kh_coef_pl(:,k,:) = Kh_coef_pl(:,k,:)* ((GRD_gz(k)+CNST_ERADIUS)/CNST_ERADIUS)**(2*lap_order_hdiff)
          !
          Kh_coef_h(:,k,:) = Kh_coef_h(:,k,:)*((GRD_gzh(k)+CNST_ERADIUS)/CNST_ERADIUS)**(2*lap_order_hdiff)
          Kh_coef_h_pl(:,k,:) = Kh_coef_h_pl(:,k,:)*((GRD_gzh(k)+CNST_ERADIUS)/CNST_ERADIUS)**(2*lap_order_hdiff)
          !
          divdamp_coef(:,k,:) = divdamp_coef(:,k,:)*((GRD_gz(k)+CNST_ERADIUS)/CNST_ERADIUS)**(2*lap_order_divdamp)
          divdamp_coef_pl(:,k,:) = divdamp_coef_pl(:,k,:)*((GRD_gz(k)+CNST_ERADIUS)/CNST_ERADIUS)**(2*lap_order_divdamp)

          !S.Iga060110=>
          Kh_coef_lap1(:,k,:) = Kh_coef_lap1(:,k,:)* ((GRD_gz(k)+CNST_ERADIUS)/CNST_ERADIUS)**2
          Kh_coef_lap1_pl(:,k,:) = Kh_coef_lap1_pl(:,k,:)* ((GRD_gz(k)+CNST_ERADIUS)/CNST_ERADIUS)**2
          !
          Kh_coef_lap1_h(:,k,:) = Kh_coef_lap1_h(:,k,:)*((GRD_gzh(k)+CNST_ERADIUS)/CNST_ERADIUS)**2
          Kh_coef_lap1_h_pl(:,k,:) = Kh_coef_lap1_h_pl(:,k,:)*((GRD_gzh(k)+CNST_ERADIUS)/CNST_ERADIUS)**2
          !S.Iga060110>=
       enddo
    end if

    if ( debug_output ) then
       call output_info
    endif

    return
  end subroutine numfilter_setup

  !-----------------------------------------------------------------------------
  subroutine output_info
    use mod_adm, only :  &
         ADM_LOG_FID,    &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_prc_me,     &
         ADM_prc_pl,     &
         ADM_kmin,       &
         ADM_kmax
    use mod_grd, only :  &
         GRD_gz,         &
         GRD_dgz,        &
         GRD_gzh,        &
         GRD_dgzh
    use mod_cnst, only : &
         CNST_PI
    use mod_gtl, only :  &
         GTL_max_k,      &
         GTL_min_k
    use mod_gmtr, only : &
         GMTR_area,      &
         GMTR_area_pl
    implicit none

    integer :: k,l,n
    real(8) :: eft(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: eft_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: eft_h(ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: eft_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: eft_max,eft_min
    real(8) :: coef_max,coef_min
    real(8) :: eft_h_max,eft_h_min
    real(8) :: coef_h_max,coef_h_min
    !
    !
    !--- Horizontal numerical diffusion output
    eft_h = 0.D0
    eft_h_pl = 0.D0
    eft=0.D0
    eft_pl=0.D0
    do l = 1, ADM_lall
       do k = 1, ADM_kall
          do n=1,ADM_gall
             eft(n,k,l) = ((sqrt(GMTR_area(n,l))/CNST_PI)&
                  **(2*lap_order_hdiff))&
                  /(Kh_coef(n,k,l)+1.0D-99)
             eft_h(n,k,l) = ((sqrt(GMTR_area(n,l))/CNST_PI)&
                  **(2*lap_order_hdiff))&
                  /(Kh_coef_h(n,k,l)+1.0D-99)
          enddo
       enddo
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             do n=1,ADM_gall_pl
                eft_pl(n,k,l) = ((sqrt(GMTR_area_pl(n,l))/CNST_PI)&
                     **(2*lap_order_hdiff))&
                     /(Kh_coef_pl(n,k,l)+1.0D-99)
                eft_h_pl(n,k,l) = ((sqrt(GMTR_area_pl(n,l))/CNST_PI)&
                     **(2*lap_order_hdiff))&
                     /(Kh_coef_h_pl(n,k,l)+1.0D-99)
             enddo
          enddo
       enddo
    end if
    !
    write(ADM_LOG_FID,*) '   --- Horizontal numerical diffusion ---'    
    write(ADM_LOG_FID,*) '      z      max coef      min coef  max eft(2DX)  min eft(2DX)'
    do k = ADM_kmin, ADM_kmax+1
       eft_max = GTL_max_k(eft,eft_pl,k)
       eft_min = GTL_min_k(eft,eft_pl,k)
       eft_h_max = GTL_max_k(eft,eft_pl,k)
       eft_h_min = GTL_min_k(eft,eft_pl,k)
       coef_max = GTL_max_k(Kh_coef,Kh_coef_pl,k)
       coef_min = GTL_min_k(Kh_coef,Kh_coef_pl,k)
       coef_h_max = GTL_max_k(Kh_coef_h,Kh_coef_h_pl,k)
       coef_h_min = GTL_min_k(Kh_coef_h,Kh_coef_h_pl,k)
       write(ADM_LOG_FID,'(F8.2,4E14.6)') GRD_gzh(k),  &
            coef_h_min, coef_h_max,                    &
            eft_h_max,eft_h_min
       if(k/=ADM_kmax+1) then
          write(ADM_LOG_FID,'(F8.2,4E14.6)') GRD_gz(k),  &
               coef_min, coef_max,                       &
               eft_max,eft_min
       end if
    enddo

    return
  end subroutine output_info

  !-----------------------------------------------------------------------------
  !> setup coefficient for rayleigh damping
  subroutine numfilter_rayleigh_damping_setup( &
       alpha, &
       zlimit )
    use mod_adm, only: &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       PI => CNST_PI
    use mod_grd, only: &
       GRD_gz,  &
       GRD_gzh
    implicit none

    real(8), intent(in) :: alpha  ! coefficient for rayleigh damping
    real(8), intent(in) :: zlimit ! lower limit of rayleigh damping [m]

    real(8) :: Htop

    integer :: k
    !---------------------------------------------------------------------------

    if ( alpha == 0.D0 ) NUMFILTER_DOrayleigh = .false.

    allocate( rayleigh_coef  (ADM_kall) )
    allocate( rayleigh_coef_h(ADM_kall) )
    rayleigh_coef  (:) = 0.D0
    rayleigh_coef_h(:) = 0.D0

    Htop = GRD_gzh(ADM_kmax+1)

    do k = 1, ADM_kall
       if ( GRD_gz(k) >= zlimit ) then
          rayleigh_coef(k) = 0.5D0 * alpha * ( 1.D0 - cos( PI * ( GRD_gz(k)-zlimit ) / ( Htop-zlimit ) ) )
       endif
    enddo

    do k = 1, ADM_kall
       if ( GRD_gzh(k) >= zlimit ) then
          rayleigh_coef_h(k) = 0.5D0 * alpha * ( 1.D0 - cos( PI * ( GRD_gzh(k)-zlimit ) / ( Htop-zlimit ) ) )
       endif
    enddo

    if ( debug_output ) then
       write(ADM_LOG_FID,*) '   --- Rayleigh damping ---'
       write(ADM_LOG_FID,*) '       z      ray.coef   e-time(2DX)'
       do k = ADM_kmin, ADM_kmax
          write(ADM_LOG_FID,'(F8.2,3E14.6)') GRD_gzh(k), rayleigh_coef_h(k), 1.D0 / ( rayleigh_coef_h(k) + 1.D-99 )
          write(ADM_LOG_FID,'(F8.2,3E14.6)') GRD_gz (k), rayleigh_coef  (k), 1.D0 / ( rayleigh_coef  (k) + 1.D-99 )
       enddo
       k = ADM_kmax + 1
       write(ADM_LOG_FID,'(F8.2,3E14.6)') GRD_gzh(k), rayleigh_coef_h(k), 1.D0 / ( rayleigh_coef_h(k) + 1.D-99 )
    endif

    return
  end subroutine numfilter_rayleigh_damping_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for vertical numerical diffusion
  subroutine numfilter_vertical_diffusion_setup( &
       gamma )
    use mod_adm, only: &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_PI
    use mod_grd, only: &
       GRD_gz,   &
       GRD_gzh,  &
       GRD_dgz,  &
       GRD_dgzh
    use mod_time, only: &
       TIME_DTL
    implicit none

    real(8), intent(in) :: gamma ! coefficient for vertical diffusion

    integer :: k
    !---------------------------------------------------------------------------

    if ( gamma == 0.D0 ) NUMFILTER_DOverticaldiff = .false.

    allocate( Kv_coef  (ADM_kall) )
    allocate( Kv_coef_h(ADM_kall) )
    Kv_coef  (k) = 0.D0
    Kv_coef_h(k) = 0.D0

    ! 6th order vertical numerical diffusion
    do k = 1, ADM_kall
       Kv_coef  (k) = gamma * GRD_dgz (k)**6 / TIME_DTL
       Kv_coef_h(k) = gamma * GRD_dgzh(k)**6 / TIME_DTL
    enddo

    if ( debug_output ) then
       write(ADM_LOG_FID,*) '   --- Vertical numerical diffusion ---'
       write(ADM_LOG_FID,*) '       z     vdiff.coef   e-time(2DX)'
       do k = ADM_kmin, ADM_kmax
          write(ADM_LOG_FID,'(F8.2,3E14.6)') GRD_gzh(k), Kv_coef_h(k), (GRD_dgzh(k)/CNST_PI)**6 / ( Kv_coef_h(k) + 1.D-99 )
          write(ADM_LOG_FID,'(F8.2,3E14.6)') GRD_gz (k), Kv_coef  (k), (GRD_dgz (k)/CNST_PI)**6 / ( Kv_coef  (k) + 1.D-99 )
       enddo
       k = ADM_kmax + 1
       write(ADM_LOG_FID,'(F8.2,3E14.6)') GRD_gzh(k), Kv_coef_h(k), (GRD_dgzh(k)/CNST_PI)**6 / ( Kv_coef_h(k) + 1.D-99 )
    endif

    return
  end subroutine numfilter_vertical_diffusion_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for 3D divergence damping
  subroutine numfilter_divdamp_setup( &
       divdamp_type, &
       dep_hgrid,    &
       lap_order,    &
       alpha,        &
       tau,          &
       alpha_v       )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       PI    => CNST_PI,    &
       RAIR  => CNST_RAIR,  &
       GAMMA => CNST_GAMMA
    use mod_grd, only: &
       GRD_gz
    use mod_gmtr, only: &
       GMTR_area,    &
       GMTR_area_pl
    use mod_time, only: &
       TIME_DTS
    use mod_gtl, only: &
       GTL_max_k, &
       GTL_min_k
    implicit none

    character(len=*), intent(in) :: divdamp_type ! type of divergence damping
    logical,          intent(in) :: dep_hgrid    ! depend on each horizontal grid?
    integer,          intent(in) :: lap_order    ! laplacian order
    real(8),          intent(in) :: alpha        ! coefficient    for divergence damping
    real(8),          intent(in) :: tau          ! e-folding time for divergence damping
    real(8),          intent(in) :: alpha_v      ! coefficient    for divergence damping

    real(8) :: e_fold_time   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: e_fold_time_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: coef
    real(8) :: coef_max, coef_min
    real(8) :: eft_max,  eft_min

    integer :: k, l
    !---------------------------------------------------------------------------

    allocate( divdamp_coef   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( divdamp_coef_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    divdamp_coef    = 0.D0
    divdamp_coef_pl = 0.D0

    if ( divdamp_type == "DIRECT") then
       if ( alpha == 0.D0 ) NUMFILTER_DOdivdamp = .false.

       ! alpha_d is an absolute value.
       coef = alpha

       divdamp_coef   (:,:,:) = coef
       divdamp_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "NONDIM_COEF" ) then
       if ( alpha == 0.D0 ) NUMFILTER_DOdivdamp = .false.

       ! alpha_d is a non-dimensional number.
       ! alpha_d * (c_s)^p * dt^{2p-1}
       coef = alpha * ( GAMMA * RAIR * 273.D0 )**lap_order * TIME_DTS**(2*lap_order-1)

       divdamp_coef   (:,:,:) = coef
       divdamp_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "E_FOLD_TIME" ) then
       if ( tau == 0.D0 ) then
          NUMFILTER_DOdivdamp = .false.
       endif

       ! tau_d is e-folding time for 2*dx.
       if ( dep_hgrid ) then

          do l = 1, ADM_lall
          do k = 1, ADM_kall
             divdamp_coef(:,k,l) = ( sqrt(GMTR_area(:,l)) / PI )**(2*lap_order) / ( tau + 1.D-99 )
          enddo
          enddo
          if ( ADM_prc_me == ADM_prc_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                divdamp_coef_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l)) / PI )**(2*lap_order) / ( tau + 1.D-99 )
             enddo
             enddo
          endif

       else

          coef = ( sqrt(horiz_dx2) / PI )**(2*lap_order) / ( tau + 1.D-99 )

          divdamp_coef   (:,:,:) = coef
          divdamp_coef_pl(:,:,:) = coef

       endif
    endif

    if ( debug_output ) then
       e_fold_time   (:,:,:) = 0.D0
       e_fold_time_pl(:,:,:) = 0.D0

       do l = 1, ADM_lall
       do k = 1, ADM_kall
          e_fold_time(:,k,l) = ( sqrt(GMTR_area(:,l)) / PI )**(2*lap_order_divdamp) &
                             / ( divdamp_coef(:,k,l) + 1.D-99 )
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             e_fold_time_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l)) / PI )**(2*lap_order_divdamp) &
                                   / ( divdamp_coef_pl(:,k,l) + 1.D-99 )
          enddo
          enddo
       endif

       write(ADM_LOG_FID,*) '   --- 3D divergence damping ---'    
       write(ADM_LOG_FID,*) '      z      max coef      min coef  max eft(2DX)  min eft(2DX)'
       do k = ADM_kmin, ADM_kmax
          eft_max  = GTL_max_k( e_fold_time,  e_fold_time_pl,  k )
          eft_min  = GTL_min_k( e_fold_time,  e_fold_time_pl,  k )
          coef_max = GTL_max_k( divdamp_coef, divdamp_coef_pl, k )
          coef_min = GTL_min_k( divdamp_coef, divdamp_coef_pl, k )
          write(ADM_LOG_FID,'(F8.2,4E14.6)') GRD_gz(k), coef_min, coef_max, eft_max, eft_min
       enddo
    endif

    divdamp_coef_v = -alpha_v * GAMMA * RAIR * 273.D0 * TIME_DTS

    return
  end subroutine numfilter_divdamp_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for vertical numerical diffusion
  subroutine numfilter_divdamp_2d_setup( &
       divdamp_type, &
       dep_hgrid,    &
       lap_order,    &
       alpha,        &
       tau,          &
       zlimit        )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       PI    => CNST_PI,    &
       RAIR  => CNST_RAIR,  &
       GAMMA => CNST_GAMMA
    use mod_grd, only: &
       GRD_gz,  &
       GRD_gzh
    use mod_gmtr, only: &
       GMTR_area,    &
       GMTR_area_pl
    use mod_time, only: &
       TIME_DTS
    use mod_gtl, only: &
       GTL_max_k, &
       GTL_min_k
    implicit none

    character(len=*), intent(in) :: divdamp_type ! type of divergence damping
    logical,          intent(in) :: dep_hgrid    ! depend on each horizontal grid?
    integer,          intent(in) :: lap_order    ! laplacian order
    real(8),          intent(in) :: alpha        ! coefficient    for divergence damping
    real(8),          intent(in) :: tau          ! e-folding time for divergence damping
    real(8),          intent(in) :: zlimit       ! lower limit of divergence damping [m]

    real(8) :: Htop

    real(8) :: e_fold_time   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: e_fold_time_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: coef
    real(8) :: coef_max, coef_min
    real(8) :: eft_max,  eft_min

    integer :: k, l
    !---------------------------------------------------------------------------

    allocate( divdamp_2d_coef   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( divdamp_2d_coef_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    divdamp_2d_coef    = 0.D0
    divdamp_2d_coef_pl = 0.D0

    if ( divdamp_type == "DIRECT") then
       if ( alpha == 0.D0 ) NUMFILTER_DOdivdamp_2d = .false.

       ! alpha is the absolute value.
       coef = alpha

       divdamp_2d_coef   (:,:,:) = coef
       divdamp_2d_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "NONDIM_COEF" ) then
       if ( alpha == 0.D0 ) NUMFILTER_DOdivdamp_2d = .false.

       ! alpha is the non-dimensional number.
       ! alpha * (c_s)^p * dt^{2p-1}
       coef = alpha * ( GAMMA * RAIR * 273.D0 )**lap_order * TIME_DTS**(2*lap_order-1)

       divdamp_2d_coef   (:,:,:) = coef
       divdamp_2d_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "E_FOLD_TIME" ) then
       if ( tau == 0.D0 ) then
          NUMFILTER_DOdivdamp_2d = .false.
       endif

       ! tau is e-folding time for 2*dx.
       if ( dep_hgrid ) then

          do l = 1, ADM_lall
          do k = 1, ADM_kall
             divdamp_2d_coef(:,k,l) = ( sqrt(GMTR_area(:,l)) / PI )**(2*lap_order) / ( tau + 1.D-99 )
          enddo
          enddo
          if ( ADM_prc_me == ADM_prc_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                divdamp_2d_coef_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l)) / PI )**(2*lap_order) / ( tau + 1.D-99 )
             enddo
             enddo
          endif

       else

          coef = ( sqrt(horiz_dx2) / PI )**(2*lap_order) / ( tau + 1.D-99 )

          divdamp_2d_coef   (:,:,:) = coef
          divdamp_2d_coef_pl(:,:,:) = coef

       endif
    endif

    Htop = GRD_gzh(ADM_kmax+1)

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       if ( GRD_gz(k) < zlimit ) then
          divdamp_2d_coef(:,k,l) = 0.D0
       else
          divdamp_2d_coef(:,k,l) = 0.5D0 * divdamp_2d_coef(:,k,l) &
                                 * ( 1.D0 - cos( PI * ( GRD_gz(k)-zlimit ) / ( Htop-zlimit ) ) )
       endif
    enddo
    enddo

   if ( ADM_prc_me == ADM_prc_pl ) then
      do l = 1, ADM_lall_pl
      do k = 1, ADM_kall
         if ( GRD_gz(k) < zlimit ) then
            divdamp_2d_coef_pl(:,k,l) = 0.D0
         else
            divdamp_2d_coef_pl(:,k,l) = 0.5D0 * divdamp_2d_coef_pl(:,k,l) &
                                      * ( 1.D0 - cos( PI * ( GRD_gz(k)-zlimit ) / ( Htop-zlimit ) ) )
         endif
      enddo
      enddo
   endif

    if ( debug_output ) then

       e_fold_time   (:,:,:) = 0.D0
       e_fold_time_pl(:,:,:) = 0.D0

       do l = 1, ADM_lall
       do k = 1, ADM_kall
          e_fold_time(:,k,l) = ( sqrt(GMTR_area(:,l)) / PI )**(2*lap_order_divdamp) &
                             / ( divdamp_2d_coef(:,k,l) + 1.D-99 )
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             e_fold_time_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l)) / PI )**(2*lap_order_divdamp) &
                                   / ( divdamp_2d_coef_pl(:,k,l) + 1.D-99 )
          enddo
          enddo
       endif

       write(ADM_LOG_FID,*) '   --- 2D divergence damping ---'    
       write(ADM_LOG_FID,*) '      z      max coef      min coef  max eft(2DX)  min eft(2DX)'
       do k = ADM_kmin, ADM_kmax
          eft_max  = GTL_max_k( e_fold_time,  e_fold_time_pl,  k )
          eft_min  = GTL_min_k( e_fold_time,  e_fold_time_pl,  k )
          coef_max = GTL_max_k( divdamp_coef, divdamp_coef_pl, k )
          coef_min = GTL_min_k( divdamp_coef, divdamp_coef_pl, k )
          write(ADM_LOG_FID,'(F8.2,4E14.6)') GRD_gz(k), coef_min, coef_max, eft_max, eft_min
       enddo

    endif

    return
  end subroutine numfilter_divdamp_2d_setup

  !-----------------------------------------------------------------------------
  !> Rayleigh damping
  subroutine numfilter_rayleigh_damping( &
       rho,     rho_pl,     &
       vx,      vx_pl,      &
       vy,      vy_pl,      &
       vz,      vz_pl,      &
       w,       w_pl,       &
       frhogvx, frhogvx_pl, &
       frhogvy, frhogvy_pl, &
       frhogvz, frhogvz_pl, &
       frhogw,  frhogw_pl   )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall_pl, &
       ADM_lall,    &
       ADM_gall_pl, &
       ADM_gall,    &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_afac,  &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_GSGAM2,     &
       VMTR_GSGAM2_pl,  &
       VMTR_GSGAM2H,    &
       VMTR_GSGAM2H_pl
    implicit none

    real(8), intent(in)    :: rho       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rho_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: w         (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: w_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: coef

    integer :: g, k, l
    !---------------------------------------------------------------------------

    if( .NOT. NUMFILTER_DOrayleigh ) return

    if ( .NOT. rayleigh_damp_only_w ) then
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          coef = rayleigh_coef(k) * rho(g,k,l) * VMTR_GSGAM2(g,k,l)

          frhogvx(g,k,l) = frhogvx(g,k,l) - coef * vx(g,k,l)
          frhogvy(g,k,l) = frhogvy(g,k,l) - coef * vy(g,k,l)
          frhogvz(g,k,l) = frhogvz(g,k,l) - coef * vz(g,k,l)
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             coef = rayleigh_coef(k) * rho_pl(g,k,l) * VMTR_GSGAM2_pl(g,k,l)

             frhogvx_pl(g,k,l) = frhogvx_pl(g,k,l) - coef * vx_pl(g,k,l)
             frhogvy_pl(g,k,l) = frhogvy_pl(g,k,l) - coef * vy_pl(g,k,l)
             frhogvz_pl(g,k,l) = frhogvz_pl(g,k,l) - coef * vz_pl(g,k,l)
          enddo
          enddo
          enddo
       endif
    endif

    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax+1
    do g = 1, ADM_gall
       frhogw(g,k,l) = frhogw(g,k,l) - rayleigh_coef_h(k) * w(g,k,l) * VMTR_GSGAM2H(g,k,l) &
                                                          * 0.5D0 * ( GRD_afac(k) * rho(g,k  ,l) &
                                                                    + GRD_bfac(k) * rho(g,k-1,l) )
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall_pl
          frhogw_pl(g,k,l) = frhogw_pl(g,k,l) - rayleigh_coef_h(k) * w_pl(g,k,l) * VMTR_GSGAM2H_pl(g,k,l) &
                                                                   * 0.5D0 * ( GRD_afac(k) * rho_pl(g,k  ,l) &
                                                                             + GRD_bfac(k) * rho_pl(g,k-1,l) )
       enddo
       enddo
       enddo
    endif

    return
  end subroutine numfilter_rayleigh_damping

  !-----------------------------------------------------------------------------
  subroutine numfilter_vertical_diffusion( &
       rho,       rho_pl,       &
       vx,        vx_pl,        &
       vy,        vy_pl,        &
       vz,        vz_pl,        &
       w,         w_pl,         &
       tem,       tem_pl,       &
       q,         q_pl,         &
       frhog,     frhog_pl,     &
       frhogvx,   frhogvx_pl,   &
       frhogvy,   frhogvy_pl,   &
       frhogvz,   frhogvz_pl,   &
       frhogw,    frhogw_pl,    &
       frhoge,    frhoge_pl,    &
       frhogetot, frhogetot_pl, &
       frhogq,    frhogq_pl     )
    use mod_adm, only: &
       ADM_gall,    &
       ADM_kall,    &
       ADM_lall,    &
       ADM_gall_pl, &
       ADM_lall_pl, &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_CV
    use mod_grd, only: &
       GRD_rdgz,  &
       GRD_rdgzh, &
       GRD_afac,  &
       GRD_bfac
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_GSGAM2,     &
       VMTR_GSGAM2_pl,  &
       VMTR_GSGAM2H,    &
       VMTR_GSGAM2H_pl
    use mod_runconf, only: &
       TRC_VMAX
    implicit none

    real(8), intent(in)    :: rho         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rho_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: w           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: w_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: tem         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: tem_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: q           (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8), intent(in)    :: q_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(8), intent(inout) :: frhog       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhog_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvx     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvy     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvz     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogw      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhoge      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhoge_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogetot   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogq      (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8), intent(inout) :: frhogq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    integer, parameter :: vmax  = 6
    integer, parameter :: I_RHO = 1
    integer, parameter :: I_VX  = 2
    integer, parameter :: I_VY  = 3
    integer, parameter :: I_VZ  = 4
    integer, parameter :: I_W   = 5
    integer, parameter :: I_TEM = 6

    real(8) :: rhog_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: rhog_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: flux    (ADM_gall   ,ADM_kall,ADM_lall   ,vmax+TRC_VMAX)
    real(8) :: flux_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax+TRC_VMAX)
    real(8) :: vtmp0   (ADM_gall   ,ADM_kall,ADM_lall   ,vmax+TRC_VMAX)
    real(8) :: vtmp0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax+TRC_VMAX)
    real(8) :: vtmp2   (ADM_gall   ,ADM_kall,ADM_lall   ,vmax+TRC_VMAX)
    real(8) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax+TRC_VMAX)

    real(8) :: coef

    integer :: k, l, nq, p
    !---------------------------------------------------------------------------

    if( .NOT. NUMFILTER_DOverticaldiff ) return

    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
          rhog_h(:,k,l) = 0.5D0 * ( GRD_afac(k) * VMTR_GSGAM2(:,k,  l) * rho(:,k,  l) &
                                  + GRD_bfac(k) * VMTR_GSGAM2(:,k-1,l) * rho(:,k-1,l) )
       enddo
       rhog_h(:,ADM_kmin-1,l) = rhog_h(:,ADM_kmin,l)
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             rhog_h_pl(:,k,l) = 0.5D0 * ( GRD_afac(k) * VMTR_GSGAM2_pl(:,k,  l) * rho_pl(:,k  ,l) &
                                        + GRD_bfac(k) * VMTR_GSGAM2_pl(:,k-1,l) * rho_pl(:,k-1,l) )
          enddo
          rhog_h_pl(:,ADM_kmin-1,l) = rhog_h_pl(:,ADM_kmin,l)
       enddo
    endif

    vtmp0(:,:,:,I_RHO) = rho(:,:,:)
    vtmp0(:,:,:,I_VX ) = vx (:,:,:)
    vtmp0(:,:,:,I_VY ) = vy (:,:,:)
    vtmp0(:,:,:,I_VZ ) = vz (:,:,:)
    vtmp0(:,:,:,I_TEM) = tem(:,:,:)
    do nq = 1, TRC_VMAX
       vtmp0(:,:,:,vmax+nq) = rho(:,:,:) * q(:,:,:,nq)
    enddo
    vtmp0(:,:,:,I_W) = w(:,:,:)

    !--- bottom boundary
    vtmp0(:,ADM_kmin-1,:,I_RHO) = 3.D0 * vtmp0(:,ADM_kmin  ,:,I_RHO) &
                                - 3.D0 * vtmp0(:,ADM_kmin+1,:,I_RHO) &
                                + 1.D0 * vtmp0(:,ADM_kmin+2,:,I_RHO)
    vtmp0(:,ADM_kmin-1,:,I_VX ) = vtmp0(:,ADM_kmin,:,I_VX)
    vtmp0(:,ADM_kmin-1,:,I_VY ) = vtmp0(:,ADM_kmin,:,I_VY)
    vtmp0(:,ADM_kmin-1,:,I_VZ ) = vtmp0(:,ADM_kmin,:,I_VZ)
    vtmp0(:,ADM_kmin-1,:,I_TEM) = 3.D0 * vtmp0(:,ADM_kmin  ,:,I_TEM) &
                                - 3.D0 * vtmp0(:,ADM_kmin+1,:,I_TEM) &
                                + 1.D0 * vtmp0(:,ADM_kmin+2,:,I_TEM)
    do nq = 1, TRC_VMAX
       vtmp0(:,ADM_kmin-1,:,vmax+nq) = 3.D0 * vtmp0(:,ADM_kmin  ,:,vmax+nq) &
                                     - 3.D0 * vtmp0(:,ADM_kmin+1,:,vmax+nq) &
                                     + 1.D0 * vtmp0(:,ADM_kmin+2,:,vmax+nq)
    enddo
    vtmp0(:,ADM_kmin  ,:,I_W  ) = 0.D0

    !--- top boundary
    vtmp0(:,ADM_kmax+1,:,I_RHO) = 3.D0 * vtmp0(:,ADM_kmax  ,:,I_RHO) &
                                - 3.D0 * vtmp0(:,ADM_kmax-1,:,I_RHO) &
                                + 1.D0 * vtmp0(:,ADM_kmax-2,:,I_RHO)
    vtmp0(:,ADM_kmax+1,:,I_VX ) = vtmp0(:,ADM_kmax,:,I_VX)
    vtmp0(:,ADM_kmax+1,:,I_VY ) = vtmp0(:,ADM_kmax,:,I_VY)
    vtmp0(:,ADM_kmax+1,:,I_VZ ) = vtmp0(:,ADM_kmax,:,I_VZ)
    vtmp0(:,ADM_kmax+1,:,I_TEM) = 3.D0 * vtmp0(:,ADM_kmax  ,:,I_TEM) &
                                - 3.D0 * vtmp0(:,ADM_kmax-1,:,I_TEM) &
                                + 1.D0 * vtmp0(:,ADM_kmax-2,:,I_TEM)
    do nq = 1, TRC_VMAX
       vtmp0(:,ADM_kmax+1,:,nq+vmax) = 3.D0 * vtmp0(:,ADM_kmax  ,:,nq+vmax) &
                                     - 3.D0 * vtmp0(:,ADM_kmax-1,:,nq+vmax) &
                                     + 1.D0 * vtmp0(:,ADM_kmax-2,:,nq+vmax)
    enddo
    vtmp0(:,ADM_kmax+1,:,I_W  ) = 0.D0

    do l = 1, ADM_lall
       do p = 1, 2
          do k = ADM_kmin, ADM_kmax
             vtmp2(:,k,l,I_RHO) = ( ( vtmp0(:,k+1,l,I_RHO) - vtmp0(:,k  ,l,I_RHO) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_RHO) - vtmp0(:,k-1,l,I_RHO) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp2(:,k,l,I_VX ) = ( ( vtmp0(:,k+1,l,I_VX ) - vtmp0(:,k  ,l,I_VX ) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_VX ) - vtmp0(:,k-1,l,I_VX ) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp2(:,k,l,I_VY ) = ( ( vtmp0(:,k+1,l,I_VY ) - vtmp0(:,k  ,l,I_VY ) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_VY ) - vtmp0(:,k-1,l,I_VY ) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp2(:,k,l,I_VZ ) = ( ( vtmp0(:,k+1,l,I_VZ ) - vtmp0(:,k  ,l,I_VZ ) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_VZ ) - vtmp0(:,k-1,l,I_VZ ) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp2(:,k,l,I_TEM) = ( ( vtmp0(:,k+1,l,I_TEM) - vtmp0(:,k  ,l,I_TEM) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_TEM) - vtmp0(:,k-1,l,I_TEM) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             do nq = 1, TRC_VMAX
                vtmp2(:,k,l,nq+vmax) = ( ( vtmp0(:,k+1,l,nq+vmax) - vtmp0(:,k  ,l,nq+vmax) ) * GRD_rdgzh(k+1) &
                                       - ( vtmp0(:,k  ,l,nq+vmax) - vtmp0(:,k-1,l,nq+vmax) ) * GRD_rdgzh(k)   &
                                       )  * GRD_rdgz(k)
             enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
             vtmp2(:,k,l,I_W) = ( ( vtmp0(:,k+1,l,I_W) - vtmp0(:,k  ,l,I_W) ) * GRD_rdgz(k)   &
                                - ( vtmp0(:,k  ,l,I_W) - vtmp0(:,k-1,l,I_W) ) * GRD_rdgz(k-1) &
                                ) * GRD_rdgzh(k)
          enddo

          if ( p == 1 ) then
             !--- bottom boundary
             vtmp2(:,ADM_kmin-1,l,I_RHO) = vtmp2(:,ADM_kmin,l,I_RHO) * 2.D0 - vtmp2(:,ADM_kmin+1,l,I_RHO)
             vtmp2(:,ADM_kmin-1,l,I_VX ) = vtmp2(:,ADM_kmin,l,I_VX )
             vtmp2(:,ADM_kmin-1,l,I_VY ) = vtmp2(:,ADM_kmin,l,I_VY )
             vtmp2(:,ADM_kmin-1,l,I_VZ ) = vtmp2(:,ADM_kmin,l,I_VZ )
             vtmp2(:,ADM_kmin-1,l,I_TEM) = vtmp2(:,ADM_kmin,l,I_TEM) * 2.D0 - vtmp2(:,ADM_kmin+1,l,I_TEM)
             do nq = 1, TRC_VMAX
                vtmp2(:,ADM_kmin-1,l,nq+vmax) = 2.D0 * vtmp2(:,ADM_kmin,l,nq+vmax) - vtmp2(:,ADM_kmin+1,l,nq+vmax)
             enddo
             vtmp2(:,ADM_kmin,l,I_W) = vtmp2(:,ADM_kmin+1,l,I_W)

             !--- top boundary
             vtmp2(:,ADM_kmax+1,l,I_RHO) = vtmp2(:,ADM_kmax,l,I_RHO) * 2.D0 - vtmp2(:,ADM_kmax-1,l,I_RHO)
             vtmp2(:,ADM_kmax+1,l,I_VX ) = vtmp2(:,ADM_kmax,l,I_VX )
             vtmp2(:,ADM_kmax+1,l,I_VY ) = vtmp2(:,ADM_kmax,l,I_VY )
             vtmp2(:,ADM_kmax+1,l,I_VZ ) = vtmp2(:,ADM_kmax,l,I_VZ )
             vtmp2(:,ADM_kmax+1,l,I_TEM) = vtmp2(:,ADM_kmax,l,I_TEM) * 2.D0 - vtmp2(:,ADM_kmax-1,l,I_TEM)
             do nq = 1, TRC_VMAX
                vtmp2(:,ADM_kmax+1,l,nq+vmax) = 2.D0 * vtmp2(:,ADM_kmax,l,nq+vmax) - vtmp2(:,ADM_kmax-1,l,nq+vmax)
             enddo
             vtmp2(:,ADM_kmax+1,l,I_W) = vtmp2(:,ADM_kmax,l,I_W)

             vtmp0(:,:,l,:) = vtmp2(:,:,l,:)
          elseif( p == 2 ) then
             !--- bottom boundary
             vtmp2(:,ADM_kmin-1,l,I_RHO) = vtmp2(:,ADM_kmin,l,I_RHO)
             vtmp2(:,ADM_kmin-1,l,I_VX ) = vtmp2(:,ADM_kmin,l,I_VX )
             vtmp2(:,ADM_kmin-1,l,I_VY ) = vtmp2(:,ADM_kmin,l,I_VY )
             vtmp2(:,ADM_kmin-1,l,I_VZ ) = vtmp2(:,ADM_kmin,l,I_VZ )
             vtmp2(:,ADM_kmin-1,l,I_TEM) = vtmp2(:,ADM_kmin,l,I_TEM)
             do nq = 1, TRC_VMAX
                vtmp2(:,ADM_kmin-1,l,nq+vmax) = vtmp2(:,ADM_kmin,l,nq+vmax)
             enddo
             vtmp2(:,ADM_kmin,l,I_W) = vtmp2(:,ADM_kmin+1,l,I_W)

             !--- top boundary
             vtmp2(:,ADM_kmax+1,l,I_RHO) = vtmp2(:,ADM_kmax,l,I_RHO)
             vtmp2(:,ADM_kmax+1,l,I_VX ) = vtmp2(:,ADM_kmax,l,I_VX )
             vtmp2(:,ADM_kmax+1,l,I_VY ) = vtmp2(:,ADM_kmax,l,I_VY )
             vtmp2(:,ADM_kmax+1,l,I_VZ ) = vtmp2(:,ADM_kmax,l,I_VZ )
             vtmp2(:,ADM_kmax+1,l,I_TEM) = vtmp2(:,ADM_kmax,l,I_TEM)
             do nq = 1, TRC_VMAX
                vtmp2(:,ADM_kmax+1,l,nq+vmax) = vtmp2(:,ADM_kmax,l,nq+vmax)
             enddo
             vtmp2(:,ADM_kmax+1,l,I_W) = vtmp2(:,ADM_kmax,l,I_W)

             vtmp0(:,:,l,:) = vtmp2(:,:,l,:)
          endif
       enddo

       do k = ADM_kmin, ADM_kmax+1
          coef = Kv_coef_h(k) * GRD_rdgzh(k)

          flux(:,k,l,I_RHO) = coef * ( vtmp0(:,k,l,I_RHO)-vtmp0(:,k-1,l,I_RHO) ) * VMTR_GSGAM2H(:,k,l)
          flux(:,k,l,I_VX ) = coef * ( vtmp0(:,k,l,I_VX )-vtmp0(:,k-1,l,I_VX ) ) * rhog_h(:,k,l)
          flux(:,k,l,I_VY ) = coef * ( vtmp0(:,k,l,I_VY )-vtmp0(:,k-1,l,I_VY ) ) * rhog_h(:,k,l)
          flux(:,k,l,I_VZ ) = coef * ( vtmp0(:,k,l,I_VZ )-vtmp0(:,k-1,l,I_VZ ) ) * rhog_h(:,k,l)
          flux(:,k,l,I_TEM) = coef * ( vtmp0(:,k,l,I_TEM)-vtmp0(:,k-1,l,I_TEM) ) * rhog_h(:,k,l) * CNST_CV
          do nq = 1, TRC_VMAX
             flux(:,k,l,nq+vmax) = coef * ( vtmp0(:,k,l,nq+vmax) - vtmp0(:,k-1,l,nq+vmax) )
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax
          coef = Kv_coef(k) * GRD_rdgz(k)

          flux(:,k,l,I_W) = ( vtmp0(:,k+1,l,I_W)-vtmp0(:,k,l,I_W) ) * rho(:,k,l) * VMTR_GSGAM2(:,k,l)
       enddo

       !--- update tendency
       do k = ADM_kmin, ADM_kmax
          frhog    (:,k,l) = frhog    (:,k,l) + ( flux(:,k+1,l,I_RHO) - flux(:,k,l,I_RHO) ) * GRD_rdgz(k)
          frhogvx  (:,k,l) = frhogvx  (:,k,l) + ( flux(:,k+1,l,I_VX ) - flux(:,k,l,I_VX ) ) * GRD_rdgz(k)
          frhogvy  (:,k,l) = frhogvy  (:,k,l) + ( flux(:,k+1,l,I_VY ) - flux(:,k,l,I_VY ) ) * GRD_rdgz(k)
          frhogvz  (:,k,l) = frhogvz  (:,k,l) + ( flux(:,k+1,l,I_VZ ) - flux(:,k,l,I_VZ ) ) * GRD_rdgz(k)
          frhoge   (:,k,l) = frhoge   (:,k,l) + ( flux(:,k+1,l,I_TEM) - flux(:,k,l,I_TEM) ) * GRD_rdgz(k)
          frhogetot(:,k,l) = frhogetot(:,k,l) + ( flux(:,k+1,l,I_TEM) - flux(:,k,l,I_TEM) ) * GRD_rdgz(k)
          do nq = 1, TRC_VMAX
             frhogq(:,k,l,nq) = frhogq(:,k,l,nq) + ( flux(:,k+1,l,nq+vmax) - flux(:,k,l,nq+vmax) ) * GRD_rdgz(k)
          enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
          frhogw(:,k,l) = frhogw(:,k,l) + ( flux(:,k,l,I_W) - flux(:,k-1,l,I_W) ) * GRD_rdgzh(k)
       enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

       vtmp0_pl(:,:,:,I_RHO) = rho_pl(:,:,:)
       vtmp0_pl(:,:,:,I_VX ) = vx_pl (:,:,:)
       vtmp0_pl(:,:,:,I_VY ) = vy_pl (:,:,:)
       vtmp0_pl(:,:,:,I_VZ ) = vz_pl (:,:,:)
       vtmp0_pl(:,:,:,I_TEM) = tem_pl(:,:,:)
       do nq = 1, TRC_VMAX
          vtmp0_pl(:,:,:,vmax+nq) = rho_pl(:,:,:) * q_pl(:,:,:,nq)
       enddo
       vtmp0_pl(:,:,:,I_W) = w_pl(:,:,:)

       !--- bottom boundary
       vtmp0_pl(:,ADM_kmin-1,:,I_RHO) = 3.D0 * vtmp0_pl(:,ADM_kmin  ,:,I_RHO) &
                                      - 3.D0 * vtmp0_pl(:,ADM_kmin+1,:,I_RHO) &
                                      + 1.D0 * vtmp0_pl(:,ADM_kmin+2,:,I_RHO)
       vtmp0_pl(:,ADM_kmin-1,:,I_VX ) = vtmp0_pl(:,ADM_kmin,:,I_VX)
       vtmp0_pl(:,ADM_kmin-1,:,I_VY ) = vtmp0_pl(:,ADM_kmin,:,I_VY)
       vtmp0_pl(:,ADM_kmin-1,:,I_VZ ) = vtmp0_pl(:,ADM_kmin,:,I_VZ)
       vtmp0_pl(:,ADM_kmin-1,:,I_TEM) = 3.D0 * vtmp0_pl(:,ADM_kmin  ,:,I_TEM) &
                                      - 3.D0 * vtmp0_pl(:,ADM_kmin+1,:,I_TEM) &
                                      + 1.D0 * vtmp0_pl(:,ADM_kmin+2,:,I_TEM)
       do nq = 1, TRC_VMAX
          vtmp0_pl(:,ADM_kmin-1,:,vmax+nq) = 3.D0 * vtmp0_pl(:,ADM_kmin  ,:,vmax+nq) &
                                           - 3.D0 * vtmp0_pl(:,ADM_kmin+1,:,vmax+nq) &
                                           + 1.D0 * vtmp0_pl(:,ADM_kmin+2,:,vmax+nq)
       enddo
       vtmp0_pl(:,ADM_kmin,:,I_W) = 0.D0

       !--- top boundary
       vtmp0_pl(:,ADM_kmax+1,:,I_RHO) = 3.D0 * vtmp0_pl(:,ADM_kmax  ,:,I_RHO) &
                                      - 3.D0 * vtmp0_pl(:,ADM_kmax-1,:,I_RHO) &
                                      + 1.D0 * vtmp0_pl(:,ADM_kmax-2,:,I_RHO)
       vtmp0_pl(:,ADM_kmax+1,:,I_VX) = vtmp0_pl(:,ADM_kmax,:,I_VX)
       vtmp0_pl(:,ADM_kmax+1,:,I_VY) = vtmp0_pl(:,ADM_kmax,:,I_VY)
       vtmp0_pl(:,ADM_kmax+1,:,I_VZ) = vtmp0_pl(:,ADM_kmax,:,I_VZ)
       vtmp0_pl(:,ADM_kmax+1,:,I_TEM) = 3.D0 * vtmp0_pl(:,ADM_kmax  ,:,I_TEM) &
                                      - 3.D0 * vtmp0_pl(:,ADM_kmax-1,:,I_TEM) &
                                      + 1.D0 * vtmp0_pl(:,ADM_kmax-2,:,I_TEM)
       do nq = 1, TRC_VMAX
          vtmp0_pl(:,ADM_kmax+1,:,nq+vmax) = 3.D0 * vtmp0_pl(:,ADM_kmax  ,:,nq+vmax) &
                                           - 3.D0 * vtmp0_pl(:,ADM_kmax-1,:,nq+vmax) &
                                           + 1.D0 * vtmp0_pl(:,ADM_kmax-2,:,nq+vmax)
       enddo
       vtmp0_pl(:,ADM_kmax+1,:,I_W) = 0.D0

       do l = 1, ADM_lall
          do p = 1, 2
             do k = ADM_kmin, ADM_kmax
                vtmp2_pl(:,k,l,I_RHO) = ( ( vtmp0_pl(:,k+1,l,I_RHO)-vtmp0_pl(:,k,l,I_RHO) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k,l,I_RHO)-vtmp0_pl(:,k-1,l,I_RHO) ) * GRD_rdgzh(k) &
                                        ) * GRD_rdgz(k)
                vtmp2_pl(:,k,l,I_VX ) = ( ( vtmp0_pl(:,k+1,l,I_VX)-vtmp0_pl(:,k,l,I_VX) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k,l,I_VX)-vtmp0_pl(:,k-1,l,I_VX) ) * GRD_rdgzh(k) &
                                        ) * GRD_rdgz(k)
                vtmp2_pl(:,k,l,I_VY ) = ( ( vtmp0_pl(:,k+1,l,I_VY)-vtmp0_pl(:,k,l,I_VY) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k,l,I_VY)-vtmp0_pl(:,k-1,l,I_VY) ) * GRD_rdgzh(k) &
                                        ) * GRD_rdgz(k)
                vtmp2_pl(:,k,l,I_VZ ) = ( ( vtmp0_pl(:,k+1,l,I_VZ)-vtmp0_pl(:,k,l,I_VZ) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k,l,I_VZ)-vtmp0_pl(:,k-1,l,I_VZ) ) * GRD_rdgzh(k) &
                                        ) * GRD_rdgz(k)
                vtmp2_pl(:,k,l,I_TEM) = ( ( vtmp0_pl(:,k+1,l,I_TEM)-vtmp0_pl(:,k,l,I_TEM) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k,l,I_TEM)-vtmp0_pl(:,k-1,l,I_TEM) ) * GRD_rdgzh(k) &
                                        ) * GRD_rdgz(k)
                do nq = 1, TRC_VMAX
                   vtmp2_pl(:,k,l,nq+vmax) = ( ( vtmp0_pl(:,k+1,l,nq+vmax)-vtmp0_pl(:,k,l,nq+vmax) ) * GRD_rdgzh(k+1) &
                                             - ( vtmp0_pl(:,k,l,nq+vmax)-vtmp0_pl(:,k-1,l,nq+vmax) ) * GRD_rdgzh(k) &
                                             )  * GRD_rdgz(k)
                enddo
             enddo

             do k = ADM_kmin+1, ADM_kmax
                vtmp2_pl(:,k,l,I_W) = ( ( vtmp0_pl(:,k+1,l,I_W)-vtmp0_pl(:,k,l,I_W) ) * GRD_rdgz(k) &
                                      - ( vtmp0_pl(:,k,l,I_W)-vtmp0_pl(:,k-1,l,I_W) ) * GRD_rdgz(k-1) &
                                      ) * GRD_rdgzh(k)
             enddo

             if(p==1) then
                !--- bottom boundary
                vtmp2_pl(:,ADM_kmin-1,l,I_RHO)             &
                     = 2.0D0*vtmp2_pl(:,ADM_kmin,l,I_RHO)  &
                     - 1.D0 * vtmp2_pl(:,ADM_kmin+1,l,I_RHO)
                vtmp2_pl(:,ADM_kmin-1,l,I_VX) = vtmp2_pl(:,ADM_kmin,l,I_VX)
                vtmp2_pl(:,ADM_kmin-1,l,I_VY) = vtmp2_pl(:,ADM_kmin,l,I_VY)
                vtmp2_pl(:,ADM_kmin-1,l,I_VZ) = vtmp2_pl(:,ADM_kmin,l,I_VZ)
                vtmp2_pl(:,ADM_kmin-1,l,I_TEM)             &
                     = 2.0D0*vtmp2_pl(:,ADM_kmin,l,I_TEM)  &
                     - 1.D0 * vtmp2_pl(:,ADM_kmin+1,l,I_TEM)
                vtmp2_pl(:,ADM_kmin,l,I_W)= vtmp2_pl(:,ADM_kmin+1,l,I_W)
                do nq = 1, TRC_VMAX
                   vtmp2_pl(:,ADM_kmin-1,l,nq+vmax)             &
                        = 2.0D0*vtmp2_pl(:,ADM_kmin,l,nq+vmax)  &
                        - 1.D0 * vtmp2_pl(:,ADM_kmin+1,l,nq+vmax)
                enddo
                !--- top boundary
                vtmp2_pl(:,ADM_kmax+1,l,I_RHO)             &
                     = 2.0D0*vtmp2_pl(:,ADM_kmax,l,I_RHO)  &
                     - 1.D0 * vtmp2_pl(:,ADM_kmax-1,l,I_RHO)
                vtmp2_pl(:,ADM_kmax+1,l,I_VX) = vtmp2_pl(:,ADM_kmax,l,I_VX)
                vtmp2_pl(:,ADM_kmax+1,l,I_VY) = vtmp2_pl(:,ADM_kmax,l,I_VY)
                vtmp2_pl(:,ADM_kmax+1,l,I_VZ) = vtmp2_pl(:,ADM_kmax,l,I_VZ)
                vtmp2_pl(:,ADM_kmax+1,l,I_TEM)             &
                     = 2.0D0*vtmp2_pl(:,ADM_kmax,l,I_TEM)  &
                     - 1.D0 * vtmp2_pl(:,ADM_kmax-1,l,I_TEM)
                vtmp2_pl(:,ADM_kmax+1,l,I_W) = vtmp2_pl(:,ADM_kmax,l,I_W)
                do nq = 1, TRC_VMAX
                   vtmp2_pl(:,ADM_kmax+1,l,nq+vmax)             &
                        = 2.0D0*vtmp2_pl(:,ADM_kmax,l,nq+vmax)  &
                        - 1.D0 * vtmp2_pl(:,ADM_kmax-1,l,nq+vmax)
                enddo
                !
                vtmp0_pl(:,:,l,:)=vtmp2_pl(:,:,l,:)
                !
             else if(p==2) then
                !
                vtmp2_pl(:,ADM_kmin-1,l,I_RHO) = vtmp2_pl(:,ADM_kmin,l,I_RHO)
                vtmp2_pl(:,ADM_kmin-1,l,I_VX)= vtmp2_pl(:,ADM_kmin,l,I_VX)
                vtmp2_pl(:,ADM_kmin-1,l,I_VY)= vtmp2_pl(:,ADM_kmin,l,I_VY)
                vtmp2_pl(:,ADM_kmin-1,l,I_VZ)= vtmp2_pl(:,ADM_kmin,l,I_VZ)
                vtmp2_pl(:,ADM_kmin-1,l,I_TEM) = vtmp2_pl(:,ADM_kmin,l,I_TEM)
                vtmp2_pl(:,ADM_kmin,l,I_W)= vtmp2_pl(:,ADM_kmin+1,l,I_W)
                do nq = 1, TRC_VMAX
                   vtmp2_pl(:,ADM_kmin-1,l,nq+vmax) = vtmp2_pl(:,ADM_kmin,l,nq+vmax)
                enddo
                !
                vtmp2_pl(:,ADM_kmax+1,l,I_RHO) = vtmp2_pl(:,ADM_kmax,l,I_RHO)
                vtmp2_pl(:,ADM_kmax+1,l,I_VX) = vtmp2_pl(:,ADM_kmax,l,I_VX)
                vtmp2_pl(:,ADM_kmax+1,l,I_VY) = vtmp2_pl(:,ADM_kmax,l,I_VY)
                vtmp2_pl(:,ADM_kmax+1,l,I_VZ) = vtmp2_pl(:,ADM_kmax,l,I_VZ)
                vtmp2_pl(:,ADM_kmax+1,l,I_TEM) = vtmp2_pl(:,ADM_kmax,l,I_TEM)
                vtmp2_pl(:,ADM_kmax+1,l,I_W) = vtmp2_pl(:,ADM_kmax,l,I_W)
                do nq = 1, TRC_VMAX
                   vtmp2_pl(:,ADM_kmax+1,l,nq+vmax) = vtmp2_pl(:,ADM_kmax,l,nq+vmax)
                enddo
                !
                vtmp0_pl(:,:,l,:)=vtmp2_pl(:,:,l,:)
                !
             end if
          enddo
          !
          do k=ADM_kmin, ADM_kmax+1
             flux_pl(:,k,l,I_RHO) = ( vtmp0_pl(:,k,l,I_RHO) - vtmp0_pl(:,k-1,l,I_RHO))&
                  *GRD_rdgzh(k) * VMTR_GSGAM2H_pl(:,k,l)&
                  * Kv_coef_h(k)
             flux_pl(:,k,l,I_VX) = (vtmp0_pl(:,k,l,I_VX) - vtmp0_pl(:,k-1,l,I_VX)) &
                  *GRD_rdgzh(k)&
                  *rhog_h_pl(:,k,l) * Kv_coef_h(k)
             flux_pl(:,k,l,I_VY) = (vtmp0_pl(:,k,l,I_VY) - vtmp0_pl(:,k-1,l,I_VY)) &
                  *GRD_rdgzh(k)&
                  *rhog_h_pl(:,k,l) * Kv_coef_h(k)
             flux_pl(:,k,l,I_VZ) = (vtmp0_pl(:,k,l,I_VZ) - vtmp0_pl(:,k-1,l,I_VZ)) &
                  *GRD_rdgzh(k)&
                  *rhog_h_pl(:,k,l) * Kv_coef_h(k)
             flux_pl(:,k,l,I_TEM) = (vtmp0_pl(:,k,l,I_TEM) - vtmp0_pl(:,k-1,l,I_TEM))     &
                  *GRD_rdgzh(k)&
                  *rhog_h_pl(:,k,l) * Kv_coef_h(k) * CNST_CV
             do nq = 1, TRC_VMAX
                flux_pl(:,k,l,nq+vmax) = ( vtmp0_pl(:,k,l,nq+vmax) - vtmp0_pl(:,k-1,l,nq+vmax))&
                     *GRD_rdgzh(k)&
                     * Kv_coef_h(k)
             enddo
          enddo
          !
          do k=ADM_kmin, ADM_kmax
             flux_pl(:,k,l,I_W) = (vtmp0_pl(:,k+1,l,I_W) - vtmp0_pl(:,k,l,I_W))    &
                  * GRD_rdgz(k)&
                  *rho_pl(:,k,l) * VMTR_GSGAM2_pl(:,k,l) * Kv_coef(k)
          enddo
          !
          !--- update tendency
          do k =ADM_kmin,ADM_kmax
             frhog_pl(:,k,l) = frhog_pl(:,k,l) &
                  + (flux_pl(:,k+1,l,I_RHO)-flux_pl(:,k,l,I_RHO))&
                  * GRD_rdgz(k)
             frhogvx_pl(:,k,l) = frhogvx_pl(:,k,l) &
                  + (flux_pl(:,k+1,l,I_VX)-flux_pl(:,k,l,I_VX))&
                  * GRD_rdgz(k)
             frhogvy_pl(:,k,l) = frhogvy_pl(:,k,l) &
                  + (flux_pl(:,k+1,l,I_VY)-flux_pl(:,k,l,I_VY))&
                  * GRD_rdgz(k)
             frhogvz_pl(:,k,l) = frhogvz_pl(:,k,l) &
                  + (flux_pl(:,k+1,l,I_VZ)-flux_pl(:,k,l,I_VZ))&
                  * GRD_rdgz(k)
             frhoge_pl(:,k,l) = frhoge_pl(:,k,l) &
                  + (flux_pl(:,k+1,l,I_TEM)-flux_pl(:,k,l,I_TEM))&
                  * GRD_rdgz(k)
             frhogetot_pl(:,k,l) = frhogetot_pl(:,k,l) &
                  + (flux_pl(:,k+1,l,I_TEM)-flux_pl(:,k,l,I_TEM))&
                  * GRD_rdgz(k)
             do nq = 1, TRC_VMAX
                frhogq_pl(:,k,l,nq) = frhogq_pl(:,k,l,nq) &
                     + (flux_pl(:,k+1,l,nq+vmax)-flux_pl(:,k,l,nq+vmax))&
                     * GRD_rdgz(k)
             enddo
          enddo
          do k =ADM_kmin+1,ADM_kmax
             frhogw_pl(:,k,l) = frhogw_pl(:,k,l) &
                  + (flux_pl(:,k,l,I_W)-flux_pl(:,k-1,l,I_W))&
                  * GRD_rdgzh(k)
          enddo
       enddo
    end if

    call OPRT_horizontalize_vec( frhogvx, frhogvx_pl, & !--- [INOUT]
                                 frhogvy, frhogvy_pl, & !--- [INOUT]
                                 frhogvz, frhogvz_pl  ) !--- [INOUT]

    return
  end subroutine numfilter_vertical_diffusion

  !-----------------------------------------------------------------------------
  !------ Numerical diffusion
  !------     1. Calculation region of frhogvx,...., frhoge
  !------                    : (:,:,:)
  subroutine numfilter_numerical_hdiff( &
       rho,  rho_pl,            & !--- [IN]    : density
       vx,   vx_pl,             & !--- [IN]    : Vx
       vy,   vy_pl,             & !--- [IN]    : Vy
       vz,   vz_pl,             & !--- [IN]    : Vz
       w,    w_pl,              & !--- [IN]    : w
       temd, temd_pl,           & !--- [IN]    : temperature
       q,    q_pl,              & !--- [IN]    : q
       frhog,     frhog_pl,     & !--- [INOUT] : tend. of rhog
       frhogvx,   frhogvx_pl,   & !--- [INOUT] : tend. of rhogvx
       frhogvy,   frhogvy_pl,   & !--- [INOUT] : tend. of rhogvy
       frhogvz,   frhogvz_pl,   & !--- [INOUT] : tend. of rhogvz
       frhogw,    frhogw_pl,    & !--- [INOUT] : tend. of rhogw
       frhoge,    frhoge_pl,    & !--- [INOUT] : tend. of rhoge
       frhogetot, frhogetot_pl, & !--- [INOUT] : tend. of rhoge
       frhogq,    frhogq_pl     ) !--- [INOUT] : tend. of rhogq
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_cnst, only: &
       CNST_CV, &
       CNST_PI
    use mod_comm, only: &
       COMM_data_transfer
    use mod_time, only: &
       TIME_DTL
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac, &
       GRD_gz,   &
       GRD_gzh
    use mod_vmtr, only: &
       VMTR_GSGAM2,    &
       VMTR_GSGAM2_pl, &
       VMTR_GSGAM2H,   &
       VMTR_GSGAM2H_pl
    use mod_runconf, only : &
       TRC_VMAX,    &
       TRC_ADV_TYPE
    use mod_oprt, only: &
       OPRT_horizontalize_vec, &
       OPRT_laplacian,         &
       OPRT_diffusion
    use mod_bsstate, only: &
       rho_bs,   &
       rho_bs_pl
    implicit none

    real(8), intent(in)    :: rho    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: rho_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: w      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: w_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: temd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: temd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: q      (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8), intent(in)    :: q_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(8), intent(inout) :: frhog       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhog_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvx     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvy     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogvz     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogw      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhoge      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhoge_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogetot   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: frhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: frhogq      (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8), intent(inout) :: frhogq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(8) :: vtmp        (ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(8) :: vtmp_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(8) :: vtmp2       (ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(8) :: vtmp2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(8) :: vtmp_lap1   (ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(8) :: vtmp_lap1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(8) :: qtmp        (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8) :: qtmp_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(8) :: qtmp2       (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8) :: qtmp2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(8) :: qtmp_lap1   (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(8) :: qtmp_lap1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !=> [add] H.Yashiro 20120530
    real(8) :: wk   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: wk_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !<= [add] H.Yashiro 20120530
    real(8) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhog_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), parameter :: T0 = 300.D0
    real(8) :: kh_max(ADM_kall)
    real(8) :: d2T_dx2, coef

    integer :: g, k, l, nq, p
    !---------------------------------------------------------------------------

    if ( trim(hdiff_type) == "NONLINEAR1" ) then
       do k = 1, ADM_kall

          kh_max(k) = Kh_coef_maxlim

          if ( GRD_gz(k) >= ZD_hdiff_nl ) then
             kh_max(k) = Kh_coef_maxlim &
                       - 0.5D0 * ( Kh_coef_maxlim - Kh_coef_minlim ) &
                       * ( 1.D0 - cos( CNST_PI * (GRD_gz(k)-ZD_hdiff_nl) / (GRD_gzh(ADM_kmax+1)-ZD_hdiff_nl) ) )
          endif
       enddo
    endif

    rhog(:,:,:) = rho(:,:,:) * VMTR_GSGAM2(:,:,:)

    rhog_h(:,ADM_kmin,:) = 0.D0
    do l = 1, ADM_lall
    do k = ADM_kmin+1, ADM_kmax
    do g = 1, ADM_gall
       rhog_h(g,k,l) = 0.5D0 * ( GRD_afac(k) * rho(g,k,  l) &
                               + GRD_bfac(k) * rho(g,k-1,l) ) * VMTR_GSGAM2H(g,k,l)
    enddo
    enddo
    enddo

    rhog_pl(:,:,:) = rho_pl(:,:,:) * VMTR_GSGAM2_pl(:,:,:)

    rhog_h_pl(:,ADM_kmin,:) = 0.D0

    do l = 1, ADM_lall_pl
    do k = ADM_kmin+1, ADM_kmax
    do g = 1, ADM_gall_pl
       rhog_h_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * rho_pl(g,k,  l) &
                                  + GRD_bfac(k) * rho_pl(g,k-1,l) ) * VMTR_GSGAM2H_pl(g,k,l)
    enddo
    enddo
    enddo

    vtmp   (:,:,:,1) = vx     (:,:,:)
    vtmp   (:,:,:,2) = vy     (:,:,:)
    vtmp   (:,:,:,3) = vz     (:,:,:)
    vtmp   (:,:,:,4) = w      (:,:,:)
    vtmp   (:,:,:,5) = temd   (:,:,:)
    vtmp   (:,:,:,6) = rho    (:,:,:) - rho_bs   (:,:,:)

    vtmp_pl(:,:,:,1) = vx_pl  (:,:,:)
    vtmp_pl(:,:,:,2) = vy_pl  (:,:,:)
    vtmp_pl(:,:,:,3) = vz_pl  (:,:,:)
    vtmp_pl(:,:,:,4) = w_pl   (:,:,:)
    vtmp_pl(:,:,:,5) = temd_pl(:,:,:)
    vtmp_pl(:,:,:,6) = rho_pl (:,:,:) - rho_bs_pl(:,:,:)

    ! copy beforehand
    if ( gamma_h_lap1 /= 0.D0 ) then
       vtmp_lap1   (:,:,:,:) = vtmp   (:,:,:,:)
       vtmp_lap1_pl(:,:,:,:) = vtmp_pl(:,:,:,:)
    else
       vtmp_lap1   (:,:,:,:) = 0.D0
       vtmp_lap1_pl(:,:,:,:) = 0.D0
    endif

    ! high order laplacian
    do p = 1, lap_order_hdiff
       ! for momentum
       call OPRT_laplacian( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & !--- [OUT]
                            vtmp (:,:,:,1), vtmp_pl (:,:,:,1)  ) !--- [IN]
 
       call OPRT_laplacian( vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & !--- [OUT]
                            vtmp (:,:,:,2), vtmp_pl (:,:,:,2)  ) !--- [IN]
 
       call OPRT_laplacian( vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & !--- [OUT]
                            vtmp (:,:,:,3), vtmp_pl (:,:,:,3)  ) !--- [IN]
 
       call OPRT_laplacian( vtmp2(:,:,:,4), vtmp2_pl(:,:,:,4), & !--- [OUT]
                            vtmp (:,:,:,4), vtmp_pl (:,:,:,4)  ) !--- [IN]
       ! for scalar
       if ( p == lap_order_hdiff ) then

          if ( trim(hdiff_type) == "NONLINEAR1" ) then
             do l = 1, ADM_lall
                do k = 1, ADM_kall
                do g = 1, ADM_gall
                   d2T_dx2 = abs(vtmp(g,k,l,5)) / T0 * horiz_dx2
                   coef    = cfact * ( horiz_dx2 * horiz_dx2 ) / TIME_DTL * d2T_dx2

                   KH_coef(g,k,l) = max( min( coef, Kh_max(k) ), Kh_coef_minlim )
                enddo
                enddo

                do k = ADM_kmin+1, ADM_kmax
                do g = 1, ADM_gall
                   KH_coef_h(g,k,l) = 0.5D0 * ( KH_coef(g,k,l) + KH_coef(g,k-1,l) )
                enddo
                enddo

                do g = 1, ADM_gall
                   KH_coef_h(g,ADM_kmin,l) = 0.D0
                enddo
             enddo

             do l = 1, ADM_lall_pl
                do k = 1, ADM_kall
                do g = 1, ADM_gall_pl
                   d2T_dx2 = abs(vtmp_pl(g,k,l,5)) / T0 * horiz_dx2
                   coef    = cfact * ( horiz_dx2 * horiz_dx2 ) / TIME_DTL * d2T_dx2

                   KH_coef_pl(g,k,l) = max( min( coef, Kh_max(k) ), Kh_coef_minlim )
                enddo
                enddo

                do k = ADM_kmin+1, ADM_kmax
                do g = 1, ADM_gall_pl
                   KH_coef_h_pl(g,k,l) = 0.5D0 * ( KH_coef_pl(g,k,l) + KH_coef_pl(g,k-1,l) )
                enddo
                enddo

                do g = 1, ADM_gall_pl
                   KH_coef_h_pl(g,ADM_kmin,l) = 0.D0
                enddo
             enddo
          endif ! nonlinear1

          wk   (:,:,:) = rhog   (:,:,:) * CNST_CV * KH_coef   (:,:,:)
          wk_pl(:,:,:) = rhog_pl(:,:,:) * CNST_CV * KH_coef_pl(:,:,:)

          call OPRT_diffusion( vtmp2(:,:,:,5), vtmp2_pl(:,:,:,5), &
                               vtmp (:,:,:,5), vtmp_pl (:,:,:,5), &
                               wk   (:,:,:),   wk_pl   (:,:,:)    )

          wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_rho * KH_coef   (:,:,:)
          wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_rho * KH_coef_pl(:,:,:)

          call OPRT_diffusion( vtmp2(:,:,:,6), vtmp2_pl(:,:,:,6), &
                               vtmp (:,:,:,6), vtmp_pl (:,:,:,6), &
                               wk   (:,:,:),   wk_pl   (:,:,:)    )
       else
          call OPRT_laplacian( vtmp2(:,:,:,5), vtmp2_pl(:,:,:,5), & !--- [OUT]
                               vtmp (:,:,:,5), vtmp_pl (:,:,:,5)  ) !--- [IN]

          call OPRT_laplacian( vtmp2(:,:,:,6), vtmp2_pl(:,:,:,6), & !--- [OUT]
                               vtmp (:,:,:,6), vtmp_pl (:,:,:,6)  ) !--- [IN]
       endif

       vtmp   (:,:,:,:) = -vtmp2   (:,:,:,:)
       vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

       call COMM_data_transfer( vtmp, vtmp_pl )

    enddo ! laplacian order loop

    !--- 1st order laplacian filter
    if ( gamma_h_lap1 /= 0.D0 ) then
       call OPRT_laplacian( vtmp2    (:,:,:,1), vtmp2_pl    (:,:,:,1), &
                            vtmp_lap1(:,:,:,1), vtmp_lap1_pl(:,:,:,1)  )

       call OPRT_laplacian( vtmp2    (:,:,:,2), vtmp2_pl    (:,:,:,2), &
                            vtmp_lap1(:,:,:,2), vtmp_lap1_pl(:,:,:,2)  )

       call OPRT_laplacian( vtmp2    (:,:,:,3), vtmp2_pl    (:,:,:,3), &
                            vtmp_lap1(:,:,:,3), vtmp_lap1_pl(:,:,:,3)  )

       call OPRT_laplacian( vtmp2    (:,:,:,4), vtmp2_pl    (:,:,:,4), &
                            vtmp_lap1(:,:,:,4), vtmp_lap1_pl(:,:,:,4)  )

       wk   (:,:,:) = rhog   (:,:,:) * CNST_CV * KH_coef_lap1   (:,:,:)
       wk_pl(:,:,:) = rhog_pl(:,:,:) * CNST_CV * KH_coef_lap1_pl(:,:,:)

       call OPRT_diffusion( vtmp2    (:,:,:,5), vtmp2_pl    (:,:,:,5), &
                            vtmp_lap1(:,:,:,5), vtmp_lap1_pl(:,:,:,5), &
                            wk       (:,:,:),   wk_pl       (:,:,:)    )

       wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_rho * KH_coef_lap1   (:,:,:)
       wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_rho * KH_coef_lap1_pl(:,:,:)

       call OPRT_diffusion( vtmp2    (:,:,:,6), vtmp2_pl    (:,:,:,6), &
                            vtmp_lap1(:,:,:,6), vtmp_lap1_pl(:,:,:,6), &
                            wk       (:,:,:),   wk_pl       (:,:,:)    )

       vtmp_lap1   (:,:,:,:) = -vtmp2   (:,:,:,:)
       vtmp_lap1_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

       call COMM_data_transfer( vtmp_lap1, vtmp_lap1_pl )
    endif

    !--- Update tendency
    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax
    do g = 1, ADM_gall
       frhogvx  (g,k,l) = frhogvx  (g,k,l) - ( vtmp     (g,k,l,1) * KH_coef     (g,k,l) &
                                             + vtmp_lap1(g,k,l,1) * KH_coef_lap1(g,k,l) &
                                             ) * rhog(g,k,l)
       frhogvy  (g,k,l) = frhogvy  (g,k,l) - ( vtmp     (g,k,l,2) * KH_coef     (g,k,l) &
                                             + vtmp_lap1(g,k,l,2) * KH_coef_lap1(g,k,l) &
                                             ) * rhog(g,k,l)
       frhogvz  (g,k,l) = frhogvz  (g,k,l) - ( vtmp     (g,k,l,3) * KH_coef     (g,k,l) &
                                             + vtmp_lap1(g,k,l,3) * KH_coef_lap1(g,k,l) &
                                             ) * rhog(g,k,l)
       frhogw   (g,k,l) = frhogw   (g,k,l) - ( vtmp     (g,k,l,4) * KH_coef_h     (g,k,l) &
                                             + vtmp_lap1(g,k,l,4) * KH_coef_lap1_h(g,k,l) &
                                             ) * rhog_h(g,k,l)
       frhoge   (g,k,l) = frhoge   (g,k,l) - ( vtmp(g,k,l,5) + vtmp_lap1(g,k,l,5) )
       frhogetot(g,k,l) = frhogetot(g,k,l) - ( vtmp(g,k,l,5) + vtmp_lap1(g,k,l,5) )
       frhog    (g,k,l) = frhog    (g,k,l) - ( vtmp(g,k,l,6) + vtmp_lap1(g,k,l,6) )
    enddo
    enddo
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall_pl
          frhogvx_pl  (g,k,l) = frhogvx_pl  (g,k,l) - ( vtmp_pl     (g,k,l,1) * KH_coef_pl     (g,k,l) &
                                                      + vtmp_lap1_pl(g,k,l,1) * KH_coef_lap1_pl(g,k,l) &
                                                      ) * rhog_pl(g,k,l)
          frhogvy_pl  (g,k,l) = frhogvy_pl  (g,k,l) - ( vtmp_pl     (g,k,l,2) * KH_coef_pl     (g,k,l) &
                                                      + vtmp_lap1_pl(g,k,l,2) * KH_coef_lap1_pl(g,k,l) &
                                                      ) * rhog_pl(g,k,l)
          frhogvz_pl  (g,k,l) = frhogvz_pl  (g,k,l) - ( vtmp_pl     (g,k,l,3) * KH_coef_pl     (g,k,l) &
                                                      + vtmp_lap1_pl(g,k,l,3) * KH_coef_lap1_pl(g,k,l) &
                                                      ) * rhog_pl(g,k,l)
          frhogw_pl   (g,k,l) = frhogw_pl   (g,k,l) - ( vtmp_pl     (g,k,l,4) * KH_coef_h_pl     (g,k,l) &
                                                      + vtmp_lap1_pl(g,k,l,4) * KH_coef_lap1_h_pl(g,k,l) &
                                                      ) * rhog_h_pl(g,k,l)
          frhoge_pl   (g,k,l) = frhoge_pl   (g,k,l) - ( vtmp_pl(g,k,l,5) + vtmp_lap1_pl(g,k,l,5) )
          frhogetot_pl(g,k,l) = frhogetot_pl(g,k,l) - ( vtmp_pl(g,k,l,5) + vtmp_lap1_pl(g,k,l,5) )
          frhog_pl    (g,k,l) = frhog_pl    (g,k,l) - ( vtmp_pl(g,k,l,6) + vtmp_lap1_pl(g,k,l,6) )
       enddo
       enddo
       enddo
    endif

    call OPRT_horizontalize_vec( frhogvx(:,:,:), frhogvx_pl(:,:,:), & !--- [INOUT]
                                 frhogvy(:,:,:), frhogvy_pl(:,:,:), & !--- [INOUT]
                                 frhogvz(:,:,:), frhogvz_pl(:,:,:)  ) !--- [INOUT]

    !---------------------------------------------------------------------------
    ! For tracer
    !---------------------------------------------------------------------------
    ! 08/04/12 [Mod] T.Mitsui, hyper diffusion is needless for tracer if MIURA2004
    !                          because that is upwind-type advection(already diffusive)
    if ( TRC_ADV_TYPE /= 'MIURA2004' ) then 

       qtmp   (:,:,:,:) = q   (:,:,:,:)
       qtmp_pl(:,:,:,:) = q_pl(:,:,:,:)

       ! copy beforehand
       if ( gamma_h_lap1 /= 0.D0 ) then
          qtmp_lap1   (:,:,:,:) = qtmp   (:,:,:,:)
          qtmp_lap1_pl(:,:,:,:) = qtmp_pl(:,:,:,:)
       else
          qtmp_lap1   (:,:,:,:) = 0.D0
          qtmp_lap1_pl(:,:,:,:) = 0.D0
       endif

       ! high order laplacian filter
       do p = 1, lap_order_hdiff

          if ( p == lap_order_hdiff ) then

             wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_q * KH_coef   (:,:,:)
             wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_q * KH_coef_pl(:,:,:)

             do nq = 1, TRC_VMAX
                call OPRT_diffusion( qtmp2(:,:,:,nq), qtmp2_pl(:,:,:,nq), &
                                     qtmp (:,:,:,nq), qtmp_pl (:,:,:,nq), &
                                     wk   (:,:,:),    wk_pl   (:,:,:)     )
             enddo
          else
             do nq = 1, TRC_VMAX
                call OPRT_laplacian( qtmp2(:,:,:,nq), qtmp2_pl(:,:,:,nq), & !--- [OUT]
                                     qtmp (:,:,:,nq), qtmp_pl (:,:,:,nq)  ) !--- [IN
             enddo
          endif

          qtmp   (:,:,:,:) = -qtmp2   (:,:,:,:)
          qtmp_pl(:,:,:,:) = -qtmp2_pl(:,:,:,:)

          call COMM_data_transfer( qtmp, qtmp_pl )

       enddo ! laplacian order loop

       !--- 1st order laplacian filter
       if ( gamma_h_lap1 /= 0.D0 ) then

          wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_q * KH_coef_lap1   (:,:,:)
          wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_q * KH_coef_lap1_pl(:,:,:)

          do nq = 1, TRC_VMAX
             call OPRT_diffusion( qtmp2    (:,:,:,nq), qtmp2_pl    (:,:,:,nq), &
                                  qtmp_lap1(:,:,:,nq), qtmp_lap1_pl(:,:,:,nq), &
                                  wk       (:,:,:),    wk_pl       (:,:,:)     )
          enddo

          qtmp_lap1   (:,:,:,:) = -qtmp2   (:,:,:,:)
          qtmp_lap1_pl(:,:,:,:) = -qtmp2_pl(:,:,:,:)

          call COMM_data_transfer( qtmp_lap1(:,:,:,:), qtmp_lap1_pl(:,:,:,:) )
       endif

       do nq = 1, TRC_VMAX
       do l  = 1, ADM_lall
       do k  = ADM_kmin, ADM_kmax
          frhogq(:,k,l,nq)=frhogq(:,k,l,nq) - ( qtmp(:,k,l,nq) + qtmp_lap1(:,k,l,nq) )
       enddo
       enddo
       enddo

       if ( ADM_prc_me == ADM_prc_pl ) then
          do nq = 1, TRC_VMAX
          do l  = 1, ADM_lall_pl
          do k  = ADM_kmin, ADM_kmax
             frhogq_pl(:,k,l,nq) = frhogq_pl(:,k,l,nq) - ( qtmp_pl(:,k,l,nq) + qtmp_lap1_pl(:,k,l,nq) )
          enddo
          enddo
          enddo
       endif

    endif ! apply filter to tracer?

    return
  end subroutine numfilter_numerical_hdiff


  !-----------------------------------------------------------------------------
  !> 3D divergence damping
  subroutine numfilter_divdamp( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       gdx,    gdx_pl,    &
       gdy,    gdy_pl,    &
       gdz,    gdz_pl,    &
       gdvz,   gdvz_pl    )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_cnst, only: &
       CNST_RAIR, &
       CNST_GAMMA
    use mod_grd, only:  &
       GRD_rdgzh
    use mod_time, only: &
       TIME_DTS
    use mod_comm, only: &
       COMM_data_transfer
    use mod_oprt, only: &
       OPRT_horizontalize_vec, &
       OPRT_divdamp
    use mod_oprt3d, only: &
       OPRT3D_divdamp
    use mod_src, only: &
       src_flux_convergence, &
       I_SRC_default
    implicit none

    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gdx      (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( gam2 X G^{1/2} )
    real(8), intent(out) :: gdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gdy      (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( gam2 X G^{1/2} )
    real(8), intent(out) :: gdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gdz      (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( gam2 X G^{1/2} )
    real(8), intent(out) :: gdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gdvz     (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( gam2 X G^{1/2} )
    real(8), intent(out) :: gdvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: vtmp    (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(8) :: vtmp_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(8) :: vtmp2   (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(8) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    real(8) :: cnv     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: cnv_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l, p
    !---------------------------------------------------------------------------

    if ( .NOT. NUMFILTER_DOdivdamp ) then
       gdx    (:,:,:) = 0.D0
       gdx_pl (:,:,:) = 0.D0
       gdy    (:,:,:) = 0.D0
       gdy_pl (:,:,:) = 0.D0
       gdz    (:,:,:) = 0.D0
       gdz_pl (:,:,:) = 0.D0
       gdvz   (:,:,:) = 0.D0
       gdvz_pl(:,:,:) = 0.D0
       return
    endif

    !--- 3D divergence divdamp
    call OPRT3D_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & !--- [OUT]
                         vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & !--- [OUT]
                         vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & !--- [OUT]
                         rhogvx(:,:,:),  rhogvx_pl(:,:,:),  & !--- [IN]
                         rhogvy(:,:,:),  rhogvy_pl(:,:,:),  & !--- [IN]
                         rhogvz(:,:,:),  rhogvz_pl(:,:,:),  & !--- [IN]
                         rhogw (:,:,:),  rhogw_pl (:,:,:)   ) !--- [IN]

    if ( lap_order_divdamp > 1 ) then
       do p = 1, lap_order_divdamp-1

          call COMM_data_transfer(vtmp2,vtmp2_pl)

          !--- note : sign changes
          vtmp   (:,:,:,:) = -vtmp2(:,:,:,:)
          vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

          !--- 2D dinvergence divdamp
          call OPRT_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & !--- [OUT]
                             vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & !--- [OUT]
                             vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & !--- [OUT]
                             vtmp (:,:,:,1), vtmp_pl (:,:,:,1), & !--- [IN]
                             vtmp (:,:,:,2), vtmp_pl (:,:,:,2), & !--- [IN]
                             vtmp (:,:,:,3), vtmp_pl (:,:,:,3)  ) !--- [IN]

       enddo ! lap_order
    endif

    !--- X coeffcient
    gdx(:,:,:) = divdamp_coef(:,:,:) * vtmp2(:,:,:,1)
    gdy(:,:,:) = divdamp_coef(:,:,:) * vtmp2(:,:,:,2)
    gdz(:,:,:) = divdamp_coef(:,:,:) * vtmp2(:,:,:,3)

    if ( ADM_prc_me == ADM_prc_pl ) then
       gdx_pl(:,:,:) = divdamp_coef_pl(:,:,:) * vtmp2_pl(:,:,:,1)
       gdy_pl(:,:,:) = divdamp_coef_pl(:,:,:) * vtmp2_pl(:,:,:,2)
       gdz_pl(:,:,:) = divdamp_coef_pl(:,:,:) * vtmp2_pl(:,:,:,3)
    endif

    call OPRT_horizontalize_vec( gdx(:,:,:), gdx_pl(:,:,:), & !--- [INOUT]
                                 gdy(:,:,:), gdy_pl(:,:,:), & !--- [INOUT]
                                 gdz(:,:,:), gdz_pl(:,:,:)  ) !--- [INOUT]

    call src_flux_convergence( rhogvx(:,:,:), rhogvx_pl(:,:,:), & !--- [IN]
                               rhogvy(:,:,:), rhogvy_pl(:,:,:), & !--- [IN]
                               rhogvz(:,:,:), rhogvz_pl(:,:,:), & !--- [IN]
                               rhogw (:,:,:), rhogw_pl (:,:,:), & !--- [IN]
                               cnv   (:,:,:), cnv_pl   (:,:,:), & !--- [OUT]
                               I_SRC_default                    ) !--- [IN]

    do l = 1, ADM_lall
    do k = ADM_kmin+1, ADM_kmax
    do g = 1, ADM_gall
       gdvz(g,k,l) = divdamp_coef_v * ( cnv(g,k,l) - cnv(g,k-1,l) ) * GRD_rdgzh(k)
    enddo
    enddo
    enddo
    gdvz(:,ADM_kmin  ,:) = 0.D0
    gdvz(:,ADM_kmax+1,:) = 0.D0

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall_pl
          gdvz_pl(g,k,l) = divdamp_coef_v * ( cnv_pl(g,k,l) - cnv_pl(g,k-1,l) ) * GRD_rdgzh(k)
       enddo
       enddo
       enddo
       gdvz_pl(:,ADM_kmin  ,:) = 0.D0
       gdvz_pl(:,ADM_kmax+1,:) = 0.D0
    endif

    return
  end subroutine numfilter_divdamp

  !-----------------------------------------------------------------------------
  !> 2D dinvergence divdamp
  subroutine numfilter_divdamp_2d( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       gdx,    gdx_pl,    &
       gdy,    gdy_pl,    &
       gdz,    gdz_pl     )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl
    use mod_comm, only: &
       COMM_data_transfer
    use mod_oprt, only: &
       OPRT_horizontalize_vec, &
       OPRT_divdamp
    implicit none

    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out) :: gdx      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: gdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gdy      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: gdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gdz      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: gdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: vtmp    (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(8) :: vtmp_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(8) :: vtmp2   (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(8) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    integer :: p
    !---------------------------------------------------------------------------

    if ( .NOT. NUMFILTER_DOdivdamp_2d ) then
       gdx   (:,:,:) = 0.D0
       gdx_pl(:,:,:) = 0.D0
       gdy   (:,:,:) = 0.D0
       gdy_pl(:,:,:) = 0.D0
       gdz   (:,:,:) = 0.D0
       gdz_pl(:,:,:) = 0.D0
       return
    endif

    !--- 2D dinvergence divdamp
    call OPRT_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & !--- [OUT]
                       vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & !--- [OUT]
                       vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & !--- [OUT]
                       rhogvx(:,:,:),  rhogvx_pl(:,:,:),  & !--- [IN]
                       rhogvy(:,:,:),  rhogvy_pl(:,:,:),  & !--- [IN]
                       rhogvz(:,:,:),  rhogvz_pl(:,:,:)   ) !--- [IN]

    if ( lap_order_divdamp_2d > 1 ) then
       do p = 1, lap_order_divdamp_2d-1

          call COMM_data_transfer(vtmp2,vtmp2_pl)

          !--- note : sign changes
          vtmp   (:,:,:,:) = -vtmp2   (:,:,:,:)
          vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

          !--- 2D dinvergence divdamp
          call OPRT_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & !--- [OUT]
                             vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & !--- [OUT]
                             vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & !--- [OUT]
                             vtmp (:,:,:,1), vtmp_pl (:,:,:,1), & !--- [IN]
                             vtmp (:,:,:,2), vtmp_pl (:,:,:,2), & !--- [IN]
                             vtmp (:,:,:,3), vtmp_pl (:,:,:,3)  ) !--- [IN]

       enddo ! lap_order
    endif

    !--- X coeffcient
    gdx(:,:,:) = divdamp_2d_coef(:,:,:) * vtmp2(:,:,:,1)
    gdy(:,:,:) = divdamp_2d_coef(:,:,:) * vtmp2(:,:,:,2)
    gdz(:,:,:) = divdamp_2d_coef(:,:,:) * vtmp2(:,:,:,3)

    if ( ADM_prc_me == ADM_prc_pl ) then
       gdx_pl(:,:,:) = divdamp_2d_coef_pl(:,:,:) * vtmp2_pl(:,:,:,1)
       gdy_pl(:,:,:) = divdamp_2d_coef_pl(:,:,:) * vtmp2_pl(:,:,:,2)
       gdz_pl(:,:,:) = divdamp_2d_coef_pl(:,:,:) * vtmp2_pl(:,:,:,3)
    endif

    call OPRT_horizontalize_vec( gdx(:,:,:), gdx_pl(:,:,:), & !--- [INOUT]
                                 gdy(:,:,:), gdy_pl(:,:,:), & !--- [INOUT]
                                 gdz(:,:,:), gdz_pl(:,:,:)  ) !--- [INOUT]

    return
  end subroutine numfilter_divdamp_2d


  !------ Numerical diffusion
  !------     1. Calculation region of frhogvx,...., frhoge
  subroutine numfilter_smooth_1var(&
       s, s_pl                     &  !--- [INOUT]
       )
    use mod_adm, only :  &
         ADM_gall,       &
         ADM_kall,       &
         ADM_lall,       &
         ADM_gall_pl,    &
         ADM_lall_pl,    &
         ADM_prc_me,     &
         ADM_prc_pl
    use mod_oprt, only : &
         OPRT_laplacian
    use mod_gmtr, only : &
         GMTR_area,      &
         GMTR_area_pl
    use mod_comm, only :    &
         COMM_data_transfer
    !
    implicit none
    !
    real(8), intent(inout) :: s   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: s_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    
    real(8), parameter :: ggamma_h = 1.0D0/16.0D0/10.D0
    integer, parameter :: pmax = 10*8
    !
    integer :: p, pp
    integer :: k
    !
    real(8) :: vtmp   (ADM_gall   ,ADM_kall,ADM_lall   ,1)
    real(8) :: vtmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,1)
    real(8) :: vtmp2   (ADM_gall   ,ADM_kall,ADM_lall   ,1)
    real(8) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,1)

    do pp = 1,pmax

       vtmp=0.D0
       vtmp_pl=0.D0
       do k = 1, ADM_kall
          vtmp(:,k,:,1) = s(:,k,:)
       enddo
       if(ADM_prc_me==ADM_prc_pl) then
          do k = 1, ADM_kall
             vtmp_pl(:,k,:,1) = s_pl(:,k,:)
          enddo
       end if
       !
       call COMM_data_transfer(vtmp,vtmp_pl)
       !
       vtmp2= 0.D0
       vtmp2_pl= 0.D0
       !
       do p=1,2
          call OPRT_laplacian(                   &
               vtmp2(:,:,:,1),vtmp2_pl(:,:,:,1), &  !--- [OUT]
               vtmp (:,:,:,1),vtmp_pl (:,:,:,1))    !--- [IN]
          !
          vtmp(:,:,:,:)=-vtmp2(:,:,:,:)
          if (ADM_prc_me==ADM_prc_pl) then
             vtmp_pl(:,:,:,:)=-vtmp2_pl(:,:,:,:)
          endif
          !
          call COMM_data_transfer(vtmp,vtmp_pl)
       enddo
       !
       do k=1, ADM_kall
          s(:,k,:)=s(:,k,:) &
               - ggamma_h*(GMTR_area(:,:)**2)*vtmp(:,k,:,1)
       enddo
       !
       if(ADM_prc_me==ADM_prc_pl) then
          do k=1, ADM_kall
             s_pl(:,k,:)=s_pl(:,k,:) &
                  - ggamma_h*(GMTR_area_pl(:,:)**2)*vtmp_pl(:,k,:,1)
          enddo
       end if
    enddo

    do k = 1, ADM_kall
       vtmp(:,k,:,1) = s(:,k,:)
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
       do k = 1, ADM_kall
          vtmp_pl(:,k,:,1) = s_pl(:,k,:)
       enddo
    end if
    !
    call COMM_data_transfer(vtmp,vtmp_pl)

    s(:,:,:) = vtmp(:,:,:,1)
    s_pl(:,:,:) = vtmp_pl(:,:,:,1)
    !
    return
    !
  end subroutine numfilter_smooth_1var

end module mod_numfilter
!-------------------------------------------------------------------------------
