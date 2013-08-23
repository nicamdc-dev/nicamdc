!-------------------------------------------------------------------------------
!
!+  Basic state module
!
!-------------------------------------------------------------------------------
module mod_bsstate
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the set of basic state for non-hydrostatic
  !       model.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                11-08-13   A.Noda : add twp-ice exp
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only :  &
       ADM_MAXFNAME,   &
       ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !
  !--- density
  real(8),allocatable, public, save :: rho_bs(:,:,:)
  real(8),allocatable, public, save :: rho_bs_pl(:,:,:)
  !
  !--- pressure
  real(8),allocatable, public, save :: pre_bs(:,:,:)
  real(8),allocatable, public, save :: pre_bs_pl(:,:,:)
  !
  !--- temperature
  real(8),allocatable, public, save :: tem_bs(:,:,:)
  real(8),allocatable, public, save :: tem_bs_pl(:,:,:)
  !
  !--- pot temperature
  real(8),allocatable, public, save :: th_bs(:,:,:)
  real(8),allocatable, public, save :: th_bs_pl(:,:,:)
  !
  !--- water vap.
  real(8),allocatable, public, save :: qv_bs(:,:,:)
  real(8),allocatable, public, save :: qv_bs_pl(:,:,:)
  !
  !--- geo-potential ( g X z )
  real(8),allocatable, public, save :: phi(:,:,:)
  real(8),allocatable, public, save :: phi_pl(:,:,:)
  !
  !--- Basic state type
  character(ADM_NSYS), public, save :: ref_type = 'TEM'
  !                                  ='TEM': temerature is given.
  !                                  ='TH' : potential temperature is given.
  !                                  ='NOBASE' : no basic state
  !                                  ='INPUT'  : input
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: bsstate_setup
  public :: bsstate_input_ref
  public :: bsstate_output_ref
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !--- reference pressure at the ground
  real(8), private, save :: pre_g = 101325.0D0
  !
  !--- reference temperature at the ground
  real(8), private, save :: tem_g = 300.0D0
  !
  !--- reference pot. temperature at the ground
  real(8), private, save :: th_g = 300.0D0
  !
  !--- reference density at the ground ( calculated by using pre_g & tem_g )
  real(8), private, save :: rho_g
  !
  !--- reference Brunt-Vaisala frequency ( used if ref_type = 'TH'. )
  real(8), private, save :: BV_freq = 0.0D0
  !
  !--- lapse rate ( used if ref_type = 'TEM'. )
  real(8), private, save :: TGAMMA = 0.0D0
  !
  !--- lower boundary of constant (potential) temperature
  real(8), private, save :: ZT = 0.0D0
  !
  !--- geopotential at the ground
  real(8), parameter, public :: PHI0=0.0D0
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !--- reference phi
  real(8), allocatable, private, save :: phi_ref(:)
  !
  !--- reference density
  real(8), allocatable, private, save :: rho_ref(:)
  !
  !--- reference pressure
  real(8), allocatable, private, save :: pre_ref(:)
  !
  !--- reference temperature
  real(8), allocatable, private, save :: tem_ref(:)
  !
  !--- water vapor
  real(8), allocatable, private, save :: qv_ref(:)
  !
  !--- reference potential temperature
  real(8), allocatable, private, save :: th_ref(:)
  !
  character(ADM_MAXFNAME), private, save :: ref_fname = 'ref.dat'
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: set_referencestate
  private :: set_basicstate
  private :: output_info
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine bsstate_setup
    !------
    !------ Setup routine of this module
    !------    1. read the parameters.
    !------    2. allocate the memory for basic state variables.
    !------    3. set the reference state.
    !------    4. set the basic state.
    !------
    !
    use mod_adm, only :  &
         ADM_LOG_FID,    &
         ADM_CTL_FID,    &
         ADM_LALL_PL,    &
         ADM_GALL_PL,    &
         ADM_lall,       &
         ADM_kall,       &
         ADM_gall
    !
    implicit none
    !
    integer :: ierr
    !
    namelist / BSSTATEPARAM / &
         ref_type,            & !--- type of basic state
         ZT,                  & !--- if z>ZT, equi-temperature
         pre_g,               & !--- reference pressure
         tem_g,               & !--- reference temperature
         TGAMMA,              & !--- lapse rate when ref_type='TEM'.
         th_g,                & !--- reference potential temp. for when ref_type='TH'
         BV_freq,             & !--- Vaisala freq when ref_type='TH'
         ref_fname              !--- in/output base name if necessary.
    !
    !--- reading the parameters
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=BSSTATEPARAM,iostat=ierr)
    if(ierr<0) then
       write(ADM_LOG_FID,*) &
            'Msg : Sub[bsstate_setup_setup]/Mod[bsstate]'
       write(ADM_LOG_FID,*) &
            ' *** Not found namelist.'
       write(ADM_LOG_FID,*) &
            ' *** Use default values.'
    else if(ierr>0) then
       write(*,*) &
            'Msg : Sub[bsstate_setup_setup]/Mod[bsstate]'
       write(*,*) &
            ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
    end if
    !
    !--- allocation of reference variables
    allocate(phi_ref(ADM_kall))
    allocate(rho_ref(ADM_kall))
    allocate(pre_ref(ADM_kall))
    allocate(tem_ref(ADM_kall))
    allocate(th_ref(ADM_kall))
    allocate(qv_ref(ADM_kall))
    !
    ! add by kgoto
    ! initialize
    phi_ref=0.0d0
    rho_ref=0.0d0
    pre_ref=0.0d0
    tem_ref=0.0d0
    th_ref=0.0d0
    !--- allocation of the basic variables
    allocate(rho_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(rho_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(pre_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(pre_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(tem_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(tem_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(th_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(th_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(qv_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(qv_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(phi(ADM_gall,ADM_kall,ADM_lall))
    allocate(phi_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    !
    ! add by kgoto
    ! initialize
    rho_bs    = 0.0D0
    rho_bs_pl = 0.0D0
    pre_bs    = 0.0D0
    pre_bs_pl = 0.0D0
    tem_bs    = 0.0D0
    tem_bs_pl = 0.0D0
    th_bs    = 0.0D0
    th_bs_pl = 0.0D0
    qv_bs    = 0.0D0
    qv_bs_pl = 0.0D0
    !
    if(ref_type=='INPUT') then
       !
       !---- input reference state
       call bsstate_input_ref(ref_fname)
       !
       !--- calculation of basic state
       call set_basicstate
    else !--- other type
       !
       !--- calculation of reference state
       call set_referencestate
       !
       !--- calculation of basic state
       call set_basicstate
    end if
    !
    !--- output the information
    call output_info
    !
  end subroutine bsstate_setup
  !-----------------------------------------------------------------------------
  subroutine output_info
    !
    use mod_adm, only :  &
         ADM_LOG_FID,    &
         ADM_kall,       &
         ADM_kmin,       &
         ADM_kmax
    !
    implicit none
    integer :: k
    write(ADM_LOG_FID,*) &
         'Msg : Sub[nhm_bs_output_info]/Mod[basicstate]'
    write(ADM_LOG_FID,*) &
         ' --- Basic state type                : ',&
            trim(ref_type)
    write(ADM_LOG_FID,*) &
         ' --- Reference pressure              : ',&
         pre_g
    write(ADM_LOG_FID,*) &
         ' --- Reference temperature           : ',&
         tem_g
    write(ADM_LOG_FID,*) &
         ' --- Reference density               : ',&
         rho_g
    write(ADM_LOG_FID,*) &
         ' --- Vaisala frequency               : ',&
         BV_freq
    write(ADM_LOG_FID,*) &
         ' --- Lapse rate of temperature       : ',&
         TGAMMA
    write(ADM_LOG_FID,*) &
         ' --- Effective height                : ',&
         ZT
    write(ADM_LOG_FID,*) &
         '-------------------------------------------------------'
    write(ADM_LOG_FID,*) &
         'Level   Density  Pressure     Temp. Pot. Tem.        qv'
    do k=1,ADM_kall
       write(ADM_LOG_FID,'(I4,F12.4,F10.2,F10.2,F10.2,F10.7)') k,&
            rho_ref(k), pre_ref(k), tem_ref(k),th_ref(k), qv_ref(k)
       if(k==ADM_kmin-1) then
          write(ADM_LOG_FID,*) &
               '-------------------------------------------------------'
       end if
       if(k==ADM_kmax) then
          write(ADM_LOG_FID,*) &
               '-------------------------------------------------------'
       end if
    end do
  end subroutine output_info
  !-----------------------------------------------------------------------------
  subroutine set_referencestate
    !------
    !------ Calculate the reference state.
    !------
    !
    use mod_cnst, only : &
         CNST_EGRAV,     &
         CNST_RAIR,      &
         CNST_RVAP,      &
         CNST_PRE00,     &
         CNST_KAPPA
    use mod_adm, only :  &
         ADM_kall,       &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_prc_pl,     &
         ADM_prc_me
    use mod_grd, only :  &
         GRD_gz,         &
         GRD_dgz,        &
         GRD_gzh,        &
         GRD_afac,       &
         GRD_bfac
    !
    implicit none
    integer :: k
    real(8) :: dpre_ref_k, pre_ref_k, tem_ref_k
    real(8) :: pre_s, rho_s, total_mass0, total_mass,mass_diff_ratio
    !
    !
    !--- calculation of reference geopotential
    do k=1,ADM_kall
       phi_ref(k) = CNST_EGRAV * GRD_gz(k)
    end do
    !
    if(ref_type=='NOBASE') then
       !
       phi_ref = 0.0D0
       rho_ref = 0.0D0
       pre_ref = 0.0D0 
       tem_ref = 0.0D0
       th_ref  = 0.0D0
       qv_ref  = 0.0D0
       !
       pre_bs    = 0.0D0
       pre_bs_pl = 0.0D0
       tem_bs    = 0.0D0
       tem_bs_pl = 0.0D0
       th_bs    = 0.0D0
       th_bs_pl = 0.0D0
       qv_bs    = 0.0D0
       qv_bs_pl = 0.0D0
       rho_bs    = 0.0D0
       rho_bs_pl = 0.0D0
       !
       ! 04/12/25 M.Satoh add
       if ( tem_g /= 0.0d0 ) then
          rho_g = pre_g / CNST_RAIR / tem_g
       else
          rho_g = 0.0d0
       end if
       !
    else if ( ref_type == 'TEM' ) then
       qv_ref = 0.0D0
       !
       !---  calculation of reference temperature
       do k=1,ADM_kall
          if(GRD_gz(k) <=ZT ) then
             tem_ref(k) = tem_g - TGAMMA * GRD_gz(k)
          else
             tem_ref(k) = tem_g - TGAMMA * ZT
          end if
       end do
       !
       !--- calculation of density at the surface
       rho_g = pre_g / CNST_RAIR / tem_g
       !
       !--- calculation of reference pressure and density 
       !--- just below the ground level 
       pre_ref(ADM_kmin-1)                                &
            = pre_g                                        &
            + 0.5D0                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_g
       rho_ref(ADM_kmin-1)        &
            = pre_ref(ADM_kmin-1) &
            / CNST_RAIR           &
            / tem_ref(ADM_kmin-1)
       !
       !--- Reference pressure and density at the first level
       pre_ref(ADM_kmin)                                        &
            = pre_ref(ADM_kmin-1)                               &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_g
       rho_ref(ADM_kmin)        &
            = pre_ref(ADM_kmin) &
            / CNST_RAIR         &
            / tem_ref(ADM_kmin)
       !
       !--- Reference pressure and density at k level
       !---    In this caluculation, the hydrostatic balance equation 
       !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer 
       !---    level ( k-1/2 ). RHO is obtained by extrapolation from
       !---    the values at lower level.
       !--- 
       !--- Consistent way(?) in the scheme.
       !--- 
       !--- pre_ref(k) - pre_ref(k-1) 
       !--- = - 0.5D0 * ( GRD_afac(k) * rho_ref(k) 
       !---              +GRD_bfac(k) * rho_ref(k-1)  )
       !---   * ( phi_ref(k) - phi_ref(k-1) )
       !--- 
       !--- rho_ref(k)*CNST_RAIR*tem_ref(k)  - pre_ref(k-1) 
       !--- = - 0.5D0 * ( GRD_afac(k) * rho_ref(k) 
       !---              +GRD_bfac(k) * rho_ref(k-1)  )
       !---   * ( phi_ref(k) - phi_ref(k-1) )
       !
       !--- rho_ref(k)*( CNST_RAIR*tem_ref(k) 
       !---            + 0.5D0 * GRD_afac(k)
       !---            * ( phi_ref(k) - phi_ref(k-1) ) 
       !---            )
       !--- = pre_ref(k-1) 
       !---            - 0.5D0 * GRD_bfac(k) * rho_ref(k-1)
       !---            * ( phi_ref(k) - phi_ref(k-1) ) 
       !--- 
       do k = ADM_kmin+1, ADM_kmax+1
          rho_ref(k) = &
               ( pre_ref(k-1) &
               - 0.5D0 * GRD_bfac(k) * rho_ref(k-1) &
               * ( phi_ref(k) - phi_ref(k-1) )&
               ) / &
               ( CNST_RAIR*tem_ref(k) &
               + 0.5D0 * GRD_afac(k)  &
               * ( phi_ref(k) - phi_ref(k-1) )  &
               )
          pre_ref(k) = rho_ref(k) * CNST_RAIR * tem_ref(k)
       end do
       !
       !--- calculation of reference potential temperature
       do k = 1, ADM_kall
          th_ref(k)                                        &
               = tem_ref(k)                                &
               * ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
       end do
       !
    else if ( ref_type == 'RHO' ) then
       qv_ref = 0.0D0
       !
       !---  calculation of reference density
       total_mass0 = CNST_PRE00/CNST_EGRAV
       pre_s = CNST_PRE00
       do 
          do k = ADM_kmin-1, ADM_kmax+1
             rho_ref(k) = pre_s/CNST_RAIR/tem_g / exp(CNST_EGRAV*GRD_gz(k)/CNST_RAIR/tem_g)
          end do
          total_mass = 0.0D0
          do k = ADM_kmin, ADM_kmax
             total_mass = total_mass + GRD_dgz(k)*rho_ref(k)
          end do
          mass_diff_ratio = total_mass0/total_mass
          if(abs(mass_diff_ratio-1.0D0)<1.0D-8) then
             exit
          else
             pre_s = pre_s * mass_diff_ratio
          end if
       end do
       !
       !--- calculation of density at the surface
       rho_s = pre_s /CNST_RAIR / tem_g
       !
       !--- calculation of reference pressure and density 
       !--- just below the ground level 
       pre_ref(ADM_kmin-1)                                &
            = pre_s                                       &
            + 0.5D0                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_s
       !
       !--- Reference pressure and density at the first level
       pre_ref(ADM_kmin)                                        &
            = pre_ref(ADM_kmin-1)                               &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_s
       !
       !--- Reference pressure and density at k level
       do k = ADM_kmin+1, ADM_kmax+1
          pre_ref(k) = pre_ref(k-1) &
               - 0.5D0 * ( GRD_afac(k) * rho_ref(k)     &
                          +GRD_bfac(k) * rho_ref(k-1)  )&
               * ( phi_ref(k) - phi_ref(k-1) )
       end do
       !
       !--- calculation of reference temperature & pot temperature
       do k = 1, ADM_kall
          tem_ref(k) = pre_ref(k) / rho_ref(k) / CNST_RAIR
          th_ref(k)                                        &
               = tem_ref(k)                                &
               * ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
       end do
       !
    else if ( ref_type == 'TH' ) then
       qv_ref = 0.0D0
       !
       !--- calculation of reference pot. temp. 
       !--- just below the surface level.
       th_ref(ADM_kmin-1)                                 &
            = th_g                                         &
            / exp ( BV_freq**2 / CNST_EGRAV                     &
                   * ( GRD_gzh(ADM_kmin) - GRD_gz(ADM_kmin-1) ) )
       !
       !--- calculation of reference pot. temp. at the first level.
       th_ref(ADM_kmin)                                   &
            = th_g                                         &
            * exp ( BV_freq**2 / CNST_EGRAV                     &
                   * ( GRD_gz(ADM_kmin) - GRD_gzh(ADM_kmin) ) )
       !
       !--- calculation of reference pot. temp. at k level
       do k = ADM_kmin+1, ADM_kmax+1
          if ( GRD_gz(k) <= ZT ) then
             th_ref(k) = th_ref(k-1)                     &
                  * exp ( BV_freq**2 / CNST_EGRAV              &
                         * ( GRD_gz(k) - GRD_gz(k-1) ) )
          else
             th_ref(k) = th_ref(k-1)
          end if
       end do
       !
       !--- calculation of density at the surface
       tem_g = th_g / ( CNST_PRE00 / pre_g ) ** CNST_KAPPA
       rho_g = pre_g / CNST_RAIR / tem_g
       !
       !--- calculation of reference pressure, temperature, density  
       !--- just below the ground level.
       pre_ref(ADM_kmin-1)                                &
            = pre_g                                        &
            + 0.5D0                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_g
       tem_ref(ADM_kmin-1)                                        &
            = th_ref(ADM_kmin-1)                                  &
            / ( CNST_PRE00 / pre_ref(ADM_kmin-1) ) ** CNST_KAPPA
       rho_ref(ADM_kmin-1)        &
            = pre_ref(ADM_kmin-1) &
            / CNST_RAIR           &
            / tem_ref(ADM_kmin-1)
       !
       !--- calculation of reference pressure and density 
       !--- at the first level
       pre_ref(ADM_kmin)                                         &
            = pre_ref(ADM_kmin-1)                                &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_g
       tem_ref(ADM_kmin)                                       &
            = th_ref(ADM_kmin)                                 &
            / ( CNST_PRE00 / pre_ref(ADM_kmin) ) ** CNST_KAPPA
       rho_ref(ADM_kmin)                                       &
            = pre_ref(ADM_kmin) / CNST_RAIR / tem_ref(ADM_kmin)
       !
       !--- Reference pressure and density at k level
       !---    In this caluculation, the hydrostatic balance equation 
       !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer 
       !---    level ( k-1/2 ). RHO is obtained by extrapolation from
       !---    the values at lower level.
       !
       !--- fist guess
       do k = ADM_kmin+1, ADM_kmax+1
          pre_ref(k) = pre_ref(k-1)              &
               - ( phi_ref(k) - phi_ref(k-1) )   &
               * ( rho_ref(k-1)                  &
               + ( rho_ref(k-1) - rho_ref(k-2) ) &
               / ( GRD_gz(k-1)  - GRD_gz(k-2) )   &
               * ( GRD_gz(k)    - GRD_gz(k-1) ) * 0.5D0 )
          tem_ref(k)                                      &
               = th_ref(k)                                &
               / ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
          rho_ref(k)        &
               = pre_ref(k) &
               / CNST_RAIR  &
               / tem_ref(k)
       end do
       !--- hydro static balance adjustment
       do k = ADM_kmin+1, ADM_kmax+1
          do 
             tem_ref(k) = th_ref(k)                                &
                  / ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
             rho_ref(k) = &
                  ( pre_ref(k-1) &
                  - 0.5D0 * GRD_bfac(k) * rho_ref(k-1) &
                  * ( phi_ref(k) - phi_ref(k-1) )&
                  ) / &
                  ( CNST_RAIR*tem_ref(k) &
                  + 0.5D0 * GRD_afac(k)  &
                  * ( phi_ref(k) - phi_ref(k-1) )  &
                  )
             dpre_ref_k = rho_ref(k) * CNST_RAIR * tem_ref(k) - pre_ref(k)
             pre_ref(k) = pre_ref(k)+ dpre_ref_k
             if(abs(dpre_ref_k) < 1.0D-10 ) exit
             if(ADM_prc_me==ADM_prc_pl) then
                write(*,*) k,abs(dpre_ref_k)
             end if
          end do
       end do
    else if ( ref_type == 'TH-SP' ) then
       qv_ref = 0.0D0
       !
       !--- calculation of reference pot. temp. 
       !--- just below the surface level.
       do k=ADM_kmin-1,ADM_kmax+1
          if(GRD_gz(k)<10000.0D0) then
             th_ref(k) = th_g
          else
             th_ref(k) = th_g+10.0D0/1000.0D0*(GRD_gz(k)-10000.0D0)
          end if
       end do
       !
       !--- calculation of density at the surface
       tem_g = th_g / ( CNST_PRE00 / pre_g ) ** CNST_KAPPA
       rho_g = pre_g / CNST_RAIR / tem_g
       !
       !--- calculation of reference pressure, temperature, density  
       !--- just below the ground level.
       pre_ref(ADM_kmin-1)                                &
            = pre_g                                        &
            + 0.5D0                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_g
       tem_ref(ADM_kmin-1)                                        &
            = th_ref(ADM_kmin-1)                                  &
            / ( CNST_PRE00 / pre_ref(ADM_kmin-1) ) ** CNST_KAPPA
       rho_ref(ADM_kmin-1)        &
            = pre_ref(ADM_kmin-1) &
            / CNST_RAIR           &
            / tem_ref(ADM_kmin-1)
       !
       !--- calculation of reference pressure and density 
       !--- at the first level
       pre_ref(ADM_kmin)                                         &
            = pre_ref(ADM_kmin-1)                                &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_g
       tem_ref(ADM_kmin)                                       &
            = th_ref(ADM_kmin)                                 &
            / ( CNST_PRE00 / pre_ref(ADM_kmin) ) ** CNST_KAPPA
       rho_ref(ADM_kmin)                                       &
            = pre_ref(ADM_kmin) / CNST_RAIR / tem_ref(ADM_kmin)
       !
       !--- Reference pressure and density at k level
       !---    In this caluculation, the hydrostatic balance equation 
       !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer 
       !---    level ( k-1/2 ). RHO is obtained by extrapolation from
       !---    the values at lower level.
       !
       !--- fist guess
       do k = ADM_kmin+1, ADM_kmax+1
          pre_ref(k) = pre_ref(k-1)              &
               - ( phi_ref(k) - phi_ref(k-1) )   &
               * ( rho_ref(k-1)                  &
               + ( rho_ref(k-1) - rho_ref(k-2) ) &
               / ( GRD_gz(k-1)  - GRD_gz(k-2) )   &
               * ( GRD_gz(k)    - GRD_gz(k-1) ) * 0.5D0 )
          tem_ref(k)                                      &
               = th_ref(k)                                &
               / ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
          rho_ref(k)        &
               = pre_ref(k) &
               / CNST_RAIR  &
               / tem_ref(k)
       end do
       !--- hydro static balance adjustment
       do k = ADM_kmin+1, ADM_kmax+1
          do 
             tem_ref(k) = th_ref(k)                                &
                  / ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
             rho_ref(k) = &
                  ( pre_ref(k-1) &
                  - 0.5D0 * GRD_bfac(k) * rho_ref(k-1) &
                  * ( phi_ref(k) - phi_ref(k-1) )&
                  ) / &
                  ( CNST_RAIR*tem_ref(k) &
                  + 0.5D0 * GRD_afac(k)  &
                  * ( phi_ref(k) - phi_ref(k-1) )  &
                  )
             dpre_ref_k = rho_ref(k) * CNST_RAIR * tem_ref(k) - pre_ref(k)
             pre_ref(k) = pre_ref(k)+ dpre_ref_k
             if(abs(dpre_ref_k) < 1.0D-10 ) exit
             if(ADM_prc_me==ADM_prc_pl) then
                write(*,*) k,abs(dpre_ref_k)
             end if
          end do
       end do
    else if ( ref_type == 'OOYAMA' ) then
       call ooyama_reference
       do k = 1, ADM_kall
          pre_ref(k) = rho_ref(k) * tem_ref(k) &
               * ( (1.0D0-qv_ref(k))*CNST_RAIR+qv_ref(k)*CNST_RVAP )
          th_ref(k)                                        &
               = tem_ref(k)                                &
               * ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
       end do
    end if

!   if ( ref_type == 'GCSS_CASE1' ) then
    if ( ref_type == 'GCSS_CASE1' .or. ref_type=='TWP-ICE') then ! 11/08/13 A.Noda [add]
       call gcss_reference
       do k = 1, ADM_kall
          tem_ref(k) = th_ref(k)                                &
               / ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
          rho_ref(k) = pre_ref(k)/tem_ref(k) &
               / ( (1.0D0-qv_ref(k))*CNST_RAIR+qv_ref(k)*CNST_RVAP )
       end do
    end if


  end subroutine set_referencestate
  !-----------------------------------------------------------------------------
  subroutine gcss_reference
    !
    use mod_misc, only : &
         MISC_get_available_fid
    use mod_cnst, only : &
         CNST_RVAP
    !
    implicit none
    !
    integer :: fid
    !
!   character(128) :: fname = 'gcss_profile.dat' ! 11/08/13 A.Noda
    !
    fid = MISC_get_available_fid()
!   Open(fid,file=Trim(fname),status='old',form='unformatted') ! 11/08/03 A.Noda
    Open(fid,file=Trim(ref_fname),status='old',form='unformatted')
    read(fid) th_ref(:)
    read(fid) pre_ref(:)
    read(fid) qv_ref(:)
    close(fid)
  end subroutine gcss_reference
  !-----------------------------------------------------------------------------
  subroutine ooyama_reference
    !
    use mod_misc, only : &
         MISC_get_available_fid
    use mod_cnst, only : &
         CNST_RVAP
    use mod_adm, only :  &
         ADM_kall,       &
         ADM_kmin,       &
         ADM_kmax
    use mod_grd, only :  &
         GRD_gz
    use mod_runconf, only : &
         TRC_VMAX
    !
    !
    implicit none
    !
    integer :: fid
    !
    character(128) :: fname = 'ooyama_profile.dat'
    !
    integer :: kmax
    real(8), allocatable :: z_s(:), rho_s(:), tem_s(:), qv_s(:)
    real(8) :: rho_s_intpl(ADM_kall), tem_s_intpl(ADM_kall), qv_s_intpl(ADM_kall)
    !
    real(8) :: lag_intpl
    real(8) :: z,z1,p1,z2,p2,z3,p3
    !
    integer :: l,n,k,kk,kp
    !
    lag_intpl(z,z1,p1,z2,p2,z3,p3)             &
         = ((z-z2)*(z-z3))/((z1-z2)*(z1-z3))*p1&
         + ((z-z1)*(z-z3))/((z2-z1)*(z2-z3))*p2&
         + ((z-z1)*(z-z2))/((z3-z1)*(z3-z2))*p3

    !
    !--- read sounding data ( ooyama(2001) )
    fid = MISC_get_available_fid()
    Open(fid,file=Trim(fname),status='old',form='unformatted')
    !
    read(fid) kmax
    !
    allocate(z_s(kmax))
    allocate(rho_s(kmax))
    allocate(tem_s(kmax))
    allocate(qv_s(kmax))
    !
    read(fid) z_s(:)
    read(fid) rho_s(:)
    read(fid) tem_s(:)
    read(fid) qv_s(:)
    !
    close(fid)
    !
    !--- initialization
    do k = ADM_kmin,ADM_kmax
       do kk = 1, kmax-1
          if ( ( z_s(kk)   <= GRD_gz(k) ) .and. &
               ( z_s(kk+1) >=  GRD_gz(k) ) ) then
             if(kk==1) then
                kp=2
             else
                kp=kk
             end if
             rho_ref(k)                            &
                  = lag_intpl(GRD_gz(k), &
                  z_s(kp+1),rho_s(kp+1),           &
                  z_s(kp  ),rho_s(kp  ),           &
                  z_s(kp-1),rho_s(kp-1))
             tem_ref(k)                            &
                  = lag_intpl(GRD_gz(k), &
                  z_s(kp+1),tem_s(kp+1),           &
                  z_s(kp  ),tem_s(kp  ),           &
                  z_s(kp-1),tem_s(kp-1))
             qv_ref(k)                             &
                  = lag_intpl(GRD_gz(k), &
                  z_s(kp+1),qv_s(kp+1),            &
                  z_s(kp  ),qv_s(kp  ),            &
                  z_s(kp-1),qv_s(kp-1))
             exit
          end if
       end do
    end do
    !
    k = ADM_kmin-1
    kp = 2
    rho_ref(k)                            &
         = lag_intpl(GRD_gz(k), &
         z_s(kp+1),rho_s(kp+1),           &
         z_s(kp  ),rho_s(kp  ),           &
         z_s(kp-1),rho_s(kp-1))
    tem_ref(k)                            &
         = lag_intpl(GRD_gz(k), &
         z_s(kp+1),tem_s(kp+1),           &
         z_s(kp  ),tem_s(kp  ),           &
         z_s(kp-1),tem_s(kp-1))
    qv_ref(k)                             &
         = lag_intpl(GRD_gz(k), &
         z_s(kp+1),qv_s(kp+1),            &
         z_s(kp  ),qv_s(kp  ),            &
         z_s(kp-1),qv_s(kp-1))
    !
    k = ADM_kmax+1
    kp = kmax-1
    rho_ref(k)                            &
         = lag_intpl(GRD_gz(k), &
         z_s(kp+1),rho_s(kp+1),           &
         z_s(kp  ),rho_s(kp  ),           &
         z_s(kp-1),rho_s(kp-1))
    tem_ref(k)                            &
         = lag_intpl(GRD_gz(k), &
         z_s(kp+1),tem_s(kp+1),           &
         z_s(kp  ),tem_s(kp  ),           &
         z_s(kp-1),tem_s(kp-1))
    qv_ref(k)                             &
         = lag_intpl(GRD_gz(k), &
         z_s(kp+1),qv_s(kp+1),            &
         z_s(kp  ),qv_s(kp  ),            &
         z_s(kp-1),qv_s(kp-1))
    !
    deallocate(rho_s)
    deallocate(tem_s)
    deallocate(qv_s)
    !
  end subroutine ooyama_reference
  !-----------------------------------------------------------------------------
  subroutine set_basicstate
    !------
    !------ generation of basic state from reference state
    !------
    !
    use mod_cnst, only : &
         CNST_EGRAV,     &
         CNST_RVAP
    use mod_adm, only :  &
         ADM_LALL_PL,    &
         ADM_GALL_PL,    &
         ADM_lall,       &
         ADM_kall,       &
         ADM_gall,       &
         ADM_prc_pl,     &
         ADM_prc_me
    use mod_grd, only :  &
         GRD_Z,          &
         GRD_vz,         &
         GRD_vz_pl
    use mod_vintrpl, only : &
         VINTRPL_zstar_level
    use mod_bndcnd, only :  &
         bndcnd_thermo
    use mod_thrmdyn, only : &
         thrmdyn_th,        &
         thrmdyn_rho,       &
         thrmdyn_qd
    use mod_runconf, only : &
         TRC_VMAX,          &
         I_QV
    !
    implicit none
    !
    integer :: k,l
    !
    !---  pot. temp. ( dummy )
    real(8) :: qd(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: qd_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: q(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: q_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,TRC_VMAX)
    !
    !--- calculation of geo-potential
    phi(:,:,:)    = CNST_EGRAV * GRD_vz(:,:,:,GRD_Z)
    if(ADM_prc_me==ADM_prc_pl) then
       phi_pl(:,:,:) = CNST_EGRAV * GRD_vz_pl(:,:,:,GRD_Z)
    end if
    !
    if(ref_type=='NOBASE') then
       return
    end if
    !
    do k=1, ADM_kall
       pre_bs(:,k,:) = pre_ref(k)
       pre_bs_pl(:,k,:) = pre_ref(k)
       tem_bs(:,k,:) = tem_ref(k)
       tem_bs_pl(:,k,:) = tem_ref(k)
       qv_bs(:,k,:) = qv_ref(k)
       qv_bs_pl(:,k,:) = qv_ref(k)
    end do
    !
    !-- from z-level to zstar-level
    call VINTRPL_zstar_level( pre_bs, pre_bs_pl, .false. )
    call VINTRPL_zstar_level( tem_bs, tem_bs_pl, .false. )
    call VINTRPL_zstar_level( qv_bs, qv_bs_pl, .false. )
    !
    !--- Setting of mass concentration
    !--- Note :: The basic state is "dry" and TKE=0
    q = 0.0D0
    q_pl = 0.0D0
    q(:,:,:,I_QV) = qv_bs(:,:,:)
    q_pl(:,:,:,I_QV) = qv_bs_pl(:,:,:)
    do l=1, ADM_lall
       call thrmdyn_qd( &
            ADM_gall,   & !--- IN
            qd(:,:,l),  & !--- OUT  : dry mass concentration
            q(:,:,l,:)  & !--- IN  : mass concentration
            )
    end do
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1, ADM_lall_pl
          call thrmdyn_qd(   &
               ADM_gall_pl,  & !--- in
               qd_pl(:,:,l), & !--- OUT  : dry mass concentration
               q_pl(:,:,l,:) & !--- IN  : mass concentration
               )
       end do
    end if
    !
    !--- calculation of density
    do l=1, ADM_lall
       call thrmdyn_rho(   &
            ADM_gall,      & !--- in
            rho_bs(:,:,l), & !--- out
            pre_bs(:,:,l), & !--- in
            tem_bs(:,:,l), & !--- in
            qd(:,:,l),     & !--- in
            q(:,:,l,:)     & !--- in
            )
    end do
    if(ADM_prc_me==ADM_prc_pl) then
       do l=1, ADM_lall_pl
          call thrmdyn_rho(      &
               ADM_gall_pl,      & !--- in
               rho_bs_pl(:,:,l), & !--- out
               pre_bs_pl(:,:,l), & !--- in
               tem_bs_pl(:,:,l), & !--- in
               qd_pl(:,:,l),     & !--- in
               q_pl(:,:,l,:)     & !--- in
               )
       end do
    end if
    !
    !--- set boundary conditions of basic state
    do l=1, ADM_lall
       call bndcnd_thermo( &
            ADM_gall,      & !--- in
            tem_bs(:,:,l), & !--- inout
            rho_bs(:,:,l), & !--- inout
            pre_bs(:,:,l), & !--- inout
            phi(:,:,l)     & !--- in
            )
       call thrmdyn_th(    &
            ADM_gall,      &
            th_bs(:,:,l),  &  !--- OUT : potential temperature
            tem_bs(:,:,l), &  !--- IN  : temperature       
            pre_bs(:,:,l)  &  !--- IN  : pressure
            )
    end do

    if(ADM_prc_me==ADM_prc_pl) then
       do l=1, ADM_lall_pl
          call bndcnd_thermo(    &
               ADM_gall_pl,      & !--- in
               tem_bs_pl(:,:,l), & !--- inout
               rho_bs_pl(:,:,l), & !--- inout
               pre_bs_pl(:,:,l), & !--- inout
               phi_pl(:,:,l)     & !--- in
               )
          call thrmdyn_th(       &
               ADM_gall_pl,      &
               th_bs_pl(:,:,l),  &  !--- OUT : potential temperature
               tem_bs_pl(:,:,l), &  !--- IN  : temperature       
               pre_bs_pl(:,:,l)  &  !--- IN  : pressure
               )
       end do
    end if
    !
    return
    !
  end subroutine set_basicstate
  !-----------------------------------------------------------------------------
  subroutine bsstate_output_ref( basename )
    !
    use mod_misc, only : &
         MISC_get_available_fid
    use mod_cnst, only : &
         CNST_RVAP
    use mod_adm, only :  &
         ADM_prc_me,     &
         ADM_prc_run_master
    !
    implicit none
    Character(*), Intent(in) :: basename
    integer :: fid
    !
    !--- output
    if(ADM_prc_me==ADM_prc_run_master) then
       fid = MISC_get_available_fid()
       Open(fid,file=Trim(basename),form='unformatted')
       Write(fid) pre_ref(:)
       Write(fid) tem_ref(:)
       Write(fid) qv_ref(:)
       Close(fid)
       !
    End if
    !
  end subroutine bsstate_output_ref
  !-----------------------------------------------------------------------------
  subroutine bsstate_input_ref( basename )
    !
    use mod_misc, only : &
         MISC_get_available_fid
    use mod_cnst, only : &
         CNST_EGRAV,     &
         CNST_RAIR,      &
         CNST_RVAP,      &
         CNST_KAPPA,     &
         CNST_PRE00
    use mod_adm, only :  &
         ADM_kall
    use mod_grd, only :  &
         GRD_gz
    !
    implicit none
    Character(*), Intent(in) :: basename
    integer :: fid
    !
    integer :: k
    !
    !--- input
    fid = MISC_get_available_fid()
    Open(fid,file=Trim(basename),status='old',form='unformatted')
    read(fid) pre_ref(:)
    read(fid) tem_ref(:)
    read(fid) qv_ref(:)
    Close(fid)
    !
    !--- additional reference state.
    do k = 1, ADM_kall
       th_ref(k)                                        &
            = tem_ref(k)                                &
            * ( CNST_PRE00 / pre_ref(k) ) ** CNST_KAPPA
       phi_ref(k) = CNST_EGRAV * GRD_gz(k)
       rho_ref(k) = pre_ref(k)/ ( (1.0D0-qv_ref(k))*CNST_RAIR+qv_ref(k)*CNST_RVAP )/ tem_ref(k)
    end do
    !
  end subroutine bsstate_input_ref
  !-----------------------------------------------------------------------------
end module mod_bsstate
!-------------------------------------------------------------------------------
