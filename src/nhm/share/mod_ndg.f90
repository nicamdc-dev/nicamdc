!-----------------------------------------------------------------------------------
!
!+ Nudging module
!
!-----------------------------------------------------------------------------------
module mod_ndg
  !---------------------------------------------------------------------------------
  !
  !++ Current Corresponding Author: Y.Niwa
  !
  !++ History: 
  !      Version   Date       Comment 
  !      ---------------------------------------------------------------------------
  !      1.00      08-09-09   Y.Niwa renew
  !                08-09-10   M.Hara : modified nudging routines to calculate
  !                                    nudging coeffient depending on the distance
  !                                    from a point.
  !                09-09-08   S.Iga  frhog and frhog_pl in ndg are deleted ( suggested by ES staff)
  !                12-02-07   T.Seiki: add option to use yearly data
  !      ---------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------------
  !
  !++ Used modules
  use mod_adm, only :    &
       ADM_MAXFNAME,     &
       ADM_NSYS
  !
  !---------------------------------------------------------------------------------
  implicit none
  private
  !
  !++ Public variable
  !
  integer, parameter, private :: allowable_mts = 60
  integer, parameter, private :: allowable_vars = 10
  !
  !++ Private variables
  integer, private, save :: varmax = 5
  integer, private, save :: nU
  integer, private, save :: nV
  integer, private, save :: nW
  integer, private, save :: nT
  integer, private, save :: nP
  integer, private, save :: nQV
  integer, private, save :: nVX
  integer, private, save :: nVY
  integer, private, save :: nVZ
  !
  !--- namelist variables
  logical, private, save :: periodic = .false.
  character(ADM_NSYS), private, save :: VEL_TYPE = 'UV'   ! or VXVYVZ
  character(ADM_MAXFNAME), private, save :: inbasedir(allowable_mts) = '.'
  character(ADM_MAXFNAME), private, save :: U_basename = 'u'
  character(ADM_MAXFNAME), private, save :: V_basename = 'v'
  character(ADM_MAXFNAME), private, save :: VX_basename = 'vx'
  character(ADM_MAXFNAME), private, save :: VY_basename = 'vy'
  character(ADM_MAXFNAME), private, save :: VZ_basename = 'vz'
  character(ADM_MAXFNAME), private, save :: PVEL_basename = 'pvel'
  character(ADM_MAXFNAME), private, save :: T_basename = 'tem'
  character(ADM_MAXFNAME), private, save :: PRES_basename = 'pres'
  character(ADM_MAXFNAME), private, save :: QV_basename = 'qv'
  logical, private, save :: input_short_data = .false. 
  real(8), private, save :: tau_vel  = 1.0d10
  real(8), private, save :: tau_tem  = 1.0d10
  real(8), private, save :: tau_pres = 1.0d10 
  real(8), private, save :: tau_qv   = 1.0d10 
  integer, private, save :: ndg_kmin0 = -1
  integer, private, save :: ndg_kmax0 = 1000
  integer, private, save :: ndg_kmin1 = -1
  integer, private, save :: ndg_kmax1 = 1000
  integer, private, save :: ndg_months = 1
  integer, private, save :: ndg_start_yr = 1
  integer, private, save :: ndg_start_mt = 1
  integer, private, save :: ndg_start_dy = 1 
  integer, private, save :: ndg_start_hr = 0 
  integer, private, save :: ndg_start_mn = 0 
  integer, private, save :: ndg_start_sc = 0 
  integer, private, save :: ndg_times_dy = 4
  logical, private, save :: tem_ndg = .true.
  logical, private, save :: pres_ndg = .false.
  logical, private, save :: w_ndg = .false.
  logical, private, save :: qv_ndg = .false.
  integer, private, save :: kmax_in
  !
  logical, private, save :: ndg_flag = .false.
  real(8), allocatable, private, save :: scale_factor(:,:)
  real(8), allocatable, private, save :: offset(:,:)
  integer, private, save :: ndg_intv
  integer, private, save :: ndg_start_st
  real(8), allocatable, private, save :: rtau_vel(:,:,:)
  real(8), allocatable, private, save :: rtau_vel_pl(:,:,:)
  real(8), allocatable, private, save :: rtau_tem(:,:,:)
  real(8), allocatable, private, save :: rtau_tem_pl(:,:,:)
  real(8), allocatable, private, save :: rtau_pres(:,:,:) 
  real(8), allocatable, private, save :: rtau_pres_pl(:,:,:) 
  real(8), allocatable, private, save :: rtau_qv(:,:,:) 
  real(8), allocatable, private, save :: rtau_qv_pl(:,:,:) 
  character(ADM_MAXFNAME), private, save :: basenames(allowable_mts,allowable_vars) = ''
  integer, private, save :: rl 
  integer, private, save :: cmonth
  integer, private, save :: cyear
  integer, private, save :: rlinmt
  integer, private, save :: rlinyr ! [Add] 12/02/07 T.Seiki
  integer, private, save :: mtnum 
  logical, private, save :: rd1st = .true.
  logical, private, save :: opt_yearly_data = .false. ! [Add] 12/02/07 T.Seiki
  !
  real(8), allocatable, private, save :: var1(:,:,:,:)
  real(8), allocatable, private, save :: var2(:,:,:,:)
  real(8), allocatable, private, save :: var_bs(:,:,:,:)
  real(8), allocatable, private, save :: var_bs_pl(:,:,:,:)
  !
  ! 2008/09/10 [Add] M.Hara
  logical, private, save :: wt_ndg       = .false.      ! weighted nudging option 
                                                        ! depending on the distance from the pole
  real(8), private, save :: wt_ndg_lon_c = 135.0d0      ! lon. of the pole (-180<=v<=180)
  real(8), private, save :: wt_ndg_lat_c = 35.0d0       ! lat. of the pole (-90<=v<=90)
  real(8), private, save :: wt_ndg_min   = 0.0d0        ! min. coefficient (0<=v<=wt_ngd_max)
  real(8), private, save :: wt_ndg_max   = 1.0d0        ! max. coefficient (wt_ngd_min<=v<=1)
  real(8), private, save :: wt_ndg_halo1 = 0.0d0        ! distance from the pole to the halo1 
                                                        ! in m (0<=v<=wt_ndg_halo2)
  real(8), private, save :: wt_ndg_halo2 = 2.0015778D+7 ! distance from the pole to the halo2 
                                                        ! in m (wt_ngd_halo1<=v<=pi*r_e)
  !
  integer, parameter :: nvel_max=4
  integer, parameter :: nvel_vx=1
  integer, parameter :: nvel_vy=2
  integer, parameter :: nvel_vz=3
  integer, parameter :: nvel_w=4
  real(8), allocatable, private, save :: vel_bs(:,:,:,:)
  real(8), allocatable, private, save :: vel_bs_pl(:,:,:,:)
  !
  integer, private, save :: stepv1
  integer, private, save :: stepv2
  !
  real(8), allocatable, private, save :: vitpl_fact(:,:,:)
  integer, allocatable, private, save :: vitpl_tab(:,:)
  real(8), allocatable, private, save :: vspln_mat(:,:,:)
  real(8), allocatable, private, save :: vspln_h(:,:)
  real(8), allocatable, private, save :: vspln_dh(:,:)
  integer, allocatable, private, save :: vspln_idx(:,:)
  real(8), allocatable, private, save :: vspln_x0(:,:)
  !
  ! 2008/09/10 [Add] M.Hara
  real(8), allocatable, private, save ::wt_ndg_weight(:,:)
  real(8), allocatable, private, save ::wt_ndg_weight_pl(:,:)
  !
  !
  !++ Public procedures
  public :: ndg_setup
  public :: ndg_nudging_uvtp
  public :: ndg_nudging_q
  public :: ndg_update_var
  !
  !++ Private procedures
  private :: read_var
  private :: set_rl_mtnum
  !
  ! 2008/09/10 [Add] M.Hara
  private :: get_wt_ndg_weight
  private :: get_great_circle_dist
  !
contains
  !----------------------------------------------------------------------------------
  subroutine ndg_setup(ctime, dtime)
    !
    use mod_adm, only :      &
         ADM_MAXFNAME,       &
         ADM_LOG_FID,        &
         ADM_CTL_FID,        &
         ADM_prc_me,         &
         ADM_prc_run_master, &
         ADM_prc_tab,        &
         ADM_gall,           &
         ADM_lall,           &
         ADM_GALL_PL,        &
         ADM_LALL_PL,        &
         ADM_kall,           &
         ADM_kmin,           &
         ADM_kmax,           &
         ADM_proc_stop,      &
         ADM_vlayer
    use mod_misc, only :           &
         MISC_get_available_fid,   &
         MISC_make_idstr
    use mod_time, only :  &
         TIME_DTL
    use mod_calendar, only : &
         calendar_yh2ss,     &
         calendar_daymo,     &
         calendar_dayyr       ! [Add] 12/02/07 T.Seiki
    use mod_grd, only : &
         GRD_gz
    use mod_cnst, only : &
         CNST_PI,   &
         CNST_CV,   &
         CNST_RAIR
    use mod_vmtr, only : &
         VMTR_GSGAM2,    &
         VMTR_GSGAM2_pl
    !
    implicit none
    !
    real(8), intent(in) :: ctime
    real(8), intent(in) :: dtime
    !
    namelist /NUDGEPARAM/ &
         periodic,     &
         VEL_TYPE,     &
         inbasedir,    &
         U_basename,   &
         V_basename,   &
         VX_basename,  &
         VY_basename,  &
         VZ_basename,  &
         PVEL_basename,&
         T_basename,   &
         PRES_basename, &
         QV_basename,  &
         input_short_data, &
         tau_vel,      &
         tau_tem,      &
         tau_pres,     &
         tau_qv,       &
         ndg_kmin0,    &
         ndg_kmax0,    &
         ndg_kmin1,    &  ! 1 must be greater than 0
         ndg_kmax1,    &
         ndg_months,   &
         ndg_start_yr, &
         ndg_start_mt, &
         ndg_start_dy, &
         ndg_start_hr, &
         ndg_start_mn, &
         ndg_start_sc, &
         ndg_times_dy, &
         tem_ndg,      &
         pres_ndg,     &
         w_ndg,        &
         qv_ndg,       &
         kmax_in,      &
         ! 2012/02/07 [Add] T.Seiki
         opt_yearly_data, &
         !
         ! 2008/09/10 [Add] M.Hara
         wt_ndg,       &
         wt_ndg_lon_c, &
         wt_ndg_lat_c, &
         wt_ndg_max,   &
         wt_ndg_min,   &
         wt_ndg_halo1, &
         wt_ndg_halo2
    !
    integer :: n, l
    integer :: fid
    integer :: rgnid
    character(ADM_MAXFNAME) :: fname
    integer :: lsteps_dy
    integer :: idate(6) !" yymmddhhmmss
    real(8) :: ss
    real(8) :: ss0
    integer :: mt, dyinmt
    integer :: dyinyr ![add] 12/02/07 T.Seiki
    integer :: ierr
    real(8) :: fact
    integer :: k0, k1, k
    integer :: input_size = 4
    !
    !
    rewind(ADM_CTL_FID)  
    read(ADM_CTL_FID, nml=NUDGEPARAM, iostat=ierr)
    if(ierr < 0) then
       write(ADM_LOG_FID,*) &
            'MSG : Sub[NUDGE_setup]/Mod[NUDGE]'
       write(ADM_LOG_FID,*) &
            ' *** Not found namelist! STOP!!'
       call ADM_proc_stop
    else if(ierr > 0) then
       write(*,*) &
            'Msg : Sub[NUDGE_setup]/Mod[NUDGE]'
       write(*,*) &
            ' *** WARNING : Not appropriate names in namelist! STOP!!'
       call ADM_proc_stop
    end if
    !
    write(ADM_LOG_FID, nml=NUDGEPARAM)
    !
    kmax_in = ADM_kmax - ADM_kmin + 1
    !
    if(trim(VEL_TYPE) == 'UV') then
       varmax = 2
       nU = 1
       nV = 2
       do mt=1, ndg_months 
          basenames(mt,nU) = trim(inbasedir(mt)) // '/' // trim(U_basename)
          basenames(mt,nV) = trim(inbasedir(mt)) // '/' // trim(V_basename)
       end do
    else if(trim(VEL_TYPE) == 'VXVYVZ') then
       varmax = 3
       nVX = 1
       nVY = 2
       nVZ = 3
       do mt=1, ndg_months 
          basenames(mt,nVX) = trim(inbasedir(mt)) // '/' // trim(VX_basename)
          basenames(mt,nVY) = trim(inbasedir(mt)) // '/' // trim(VY_basename)
          basenames(mt,nVZ) = trim(inbasedir(mt)) // '/' // trim(VZ_basename)
       end do
    else
       write(ADM_LOG_FID,*) 'Msg : Sub[NUDGE_setup]/Mod[NUDGE]'
       write(ADM_LOG_FID,*) 'Invalid VEL_TYPE'
       write(ADM_LOG_FID,*) trim(VEL_TYPE)
       write(ADM_LOG_FID,*) 'STOP!!'
       call ADM_proc_stop
    end if
    !
    if(w_ndg) then
       varmax = varmax + 1
       nW = varmax
       do mt=1, ndg_months 
          basenames(mt,nW) = trim(inbasedir(mt)) // '/' // trim(PVEL_basename)
       end do
    end if
    !
    if(tem_ndg) then
       varmax = varmax + 1
       nT = varmax
       do mt=1, ndg_months 
          basenames(mt,nT) = trim(inbasedir(mt)) // '/' // trim(T_basename)
       end do
    end if
    !
    if(pres_ndg) then
       varmax = varmax + 1
       nP = varmax
       do mt=1, ndg_months 
          basenames(mt,nP) = trim(inbasedir(mt)) // '/' // trim(PRES_basename)
       end do
    end if
    !
    if(qv_ndg) then
       varmax = varmax + 1
       nQV = varmax
       do mt=1, ndg_months 
          basenames(mt,nQV) = trim(inbasedir(mt)) // '/' // trim(QV_basename)
       end do
    end if
    !
    !
    if(input_short_data) then
       allocate(scale_factor(ndg_months,varmax))
       allocate(offset(ndg_months,varmax))
    end if
    !
    !--- Check if input data files exist
    if(input_short_data) input_size = 2
    do n=1, varmax
       do mt=1, ndg_months
          do l=1, ADM_lall
             !
             rgnid = ADM_prc_tab(l, ADM_prc_me)
             call MISC_make_idstr(fname, trim(basenames(mt,n)), 'rgn', rgnid)
             fid = MISC_get_available_fid()
             open(fid, file=trim(fname), form='unformatted', &
                  access='direct', recl=input_size*ADM_gall*kmax_in, status='old', &
                  iostat=ierr)
             close(fid)
             !
             if(ierr /= 0) then
                write(ADM_LOG_FID,*) 'Msg : Sub[NUDGE_setup]/Mod[NUDGE]'
                write(ADM_LOG_FID,*) 'Cannot open the data file for nudging.'
                write(ADM_LOG_FID,*) trim(fname)
                write(ADM_LOG_FID,*) 'STOP!!'
                call ADM_proc_stop
             end if
             !
          end do
          !
          if(input_short_data) then
             fid = MISC_get_available_fid()
             fname = trim(basenames(mt,n)) // '.sinfo'
             open(fid,file=trim(fname),form='formatted',status='old', action='read', &
                  iostat=ierr)
             if(ierr /= 0) then
                write(ADM_LOG_FID,*) 'Msg : Sub[NUDGE_setup]/Mod[NUDGE]'
                write(ADM_LOG_FID,*) 'Cannot open the short data info file for nudging.'
                write(ADM_LOG_FID,*) trim(fname)
                write(ADM_LOG_FID,*) 'STOP!!'
                call ADM_proc_stop
             end if
             !
             read(fid,*) scale_factor(mt,n)
             read(fid,*) offset(mt,n)
             close(fid)
          end if
          !
       end do
    end do
    !
    !--- check vertical layer parameters
    ndg_kmin0 = max(1,min(ndg_kmin0,ADM_vlayer))
    ndg_kmin1 = max(1,min(ndg_kmin1,ADM_vlayer))
    ndg_kmax0 = max(1,min(ndg_kmax0,ADM_vlayer)) 
    ndg_kmax1 = max(1,min(ndg_kmax1,ADM_vlayer))
    !
    k0 = ndg_kmin0
    k1 = ndg_kmin1
    ndg_kmin0 = min(k0,k1)
    ndg_kmin1 = max(k0,k1)
    k0 = ndg_kmax0
    k1 = ndg_kmax1
    ndg_kmax0 = min(k0,k1)
    ndg_kmax1 = max(k0,k1)
    !
    if( ndg_kmin1 > ndg_kmax0 ) then
       write(ADM_LOG_FID,*) 'Msg : Sub[NUDGE_setup]/Mod[NUDGE]'
       write(ADM_LOG_FID,*) 'ERROR: Invalid vertical layers'
       write(*,*) 'stop!!'
       call AdM_proc_stop
    end if
    !
    ndg_kmin0 = ndg_kmin0 + ADM_kmin - 1
    ndg_kmin1 = ndg_kmin1 + ADM_kmin - 1
    ndg_kmax0 = ndg_kmax0 + ADM_kmin - 1
    ndg_kmax1 = ndg_kmax1 + ADM_kmin - 1
    !
    allocate(rtau_vel(ADM_gall,ADM_kall,ADM_lall))
    allocate(rtau_vel_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(rtau_tem(ADM_gall,ADM_kall,ADM_lall))
    allocate(rtau_tem_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(rtau_pres(ADM_gall,ADM_kall,ADM_lall))
    allocate(rtau_pres_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    allocate(rtau_qv(ADM_gall,ADM_kall,ADM_lall))
    allocate(rtau_qv_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL))
    !
    rtau_vel(:,:,:) = 1.0d0/tau_vel
    rtau_vel_pl(:,:,:) = 1.0d0/tau_vel
    rtau_tem(:,:,:) = 1.0d0/tau_tem
    rtau_tem_pl(:,:,:) = 1.0d0/tau_tem
    rtau_pres(:,:,:) = 1.0d0/tau_pres
    rtau_pres_pl(:,:,:) = 1.0d0/tau_pres
    rtau_qv(:,:,:)  = 1.0d0/tau_qv
    rtau_qv_pl(:,:,:)  = 1.0d0/tau_qv
    !
    do k=1, ADM_kall
       if( k < ndg_kmin0 ) then
          rtau_vel(:,k,:)    = 0.0d0
          rtau_vel_pl(:,k,:) = 0.0d0
          rtau_tem(:,k,:)    = 0.0d0
          rtau_tem_pl(:,k,:) = 0.0d0
          rtau_pres(:,k,:)   = 0.0d0
          rtau_pres_pl(:,k,:)= 0.0d0
          rtau_qv(:,k,:)     = 0.0d0
          rtau_qv_pl(:,k,:)  = 0.0d0
       else if ( ndg_kmin0 <= k .and. k < ndg_kmin1 ) then
          fact = 0.5d0*(1.0d0 + cos(CNST_PI* &
               (GRD_gz(k)-GRD_gz(ndg_kmin1))/(GRD_gz(ndg_kmin0)-GRD_gz(ndg_kmin1))))
          rtau_vel(:,k,:)    = rtau_vel(:,k,:)    *fact
          rtau_vel_pl(:,k,:) = rtau_vel_pl(:,k,:) *fact
          rtau_tem(:,k,:)    = rtau_tem(:,k,:)    *fact
          rtau_tem_pl(:,k,:) = rtau_tem_pl(:,k,:) *fact
          rtau_pres(:,k,:)   = rtau_pres(:,k,:)   *fact
          rtau_pres_pl(:,k,:)= rtau_pres_pl(:,k,:)*fact
          rtau_qv(:,k,:)     = rtau_qv(:,k,:)     *fact
          rtau_qv_pl(:,k,:)  = rtau_qv_pl(:,k,:)  *fact
       end if
       !
       if ( ndg_kmax0 <= k .and. k < ndg_kmax1 ) then
          fact = 0.5d0*(1.0d0 + cos(CNST_PI* &
               (GRD_gz(k)-GRD_gz(ndg_kmax0))/(GRD_gz(ndg_kmax1)-GRD_gz(ndg_kmax0))))
          rtau_vel(:,k,:)    = rtau_vel(:,k,:)    *fact
          rtau_vel_pl(:,k,:) = rtau_vel_pl(:,k,:) *fact
          rtau_tem(:,k,:)    = rtau_tem(:,k,:)    *fact
          rtau_tem_pl(:,k,:) = rtau_tem_pl(:,k,:) *fact
          rtau_pres(:,k,:)   = rtau_pres(:,k,:)   *fact
          rtau_pres_pl(:,k,:)= rtau_pres_pl(:,k,:)*fact
          rtau_qv(:,k,:)     = rtau_qv(:,k,:)     *fact
          rtau_qv_pl(:,k,:)  = rtau_qv_pl(:,k,:)  *fact
       else if( ndg_kmax1 <= k ) then
          rtau_vel(:,k,:)    = 0.0d0
          rtau_vel_pl(:,k,:) = 0.0d0
          rtau_tem(:,k,:)    = 0.0d0
          rtau_tem_pl(:,k,:) = 0.0d0
          rtau_pres(:,k,:)   = 0.0d0
          rtau_pres_pl(:,k,:)= 0.0d0
          rtau_qv(:,k,:)     = 0.0d0
          rtau_qv_pl(:,k,:)  = 0.0d0
       end if
    end do

    !
    ! 2008/09/10 [Add] M.Hara
    if (wt_ndg) then
       allocate(wt_ndg_weight(ADM_gall,ADM_lall))
       allocate(wt_ndg_weight_pl(ADM_GALL_PL,ADM_LALL_PL))
       call get_wt_ndg_weight(wt_ndg_weight(:,:),wt_ndg_weight_pl(:,:))
       do k=1, ADM_kall
          rtau_vel(:,k,:)    = rtau_vel(:,k,:)    * wt_ndg_weight(:,:)
          rtau_vel_pl(:,k,:) = rtau_vel_pl(:,k,:) * wt_ndg_weight_pl(:,:)
          rtau_tem(:,k,:)    = rtau_tem(:,k,:)    * wt_ndg_weight(:,:)
          rtau_tem_pl(:,k,:) = rtau_tem_pl(:,k,:) * wt_ndg_weight_pl(:,:)
          rtau_pres(:,k,:)   = rtau_pres(:,k,:)   * wt_ndg_weight(:,:)
          rtau_pres_pl(:,k,:)= rtau_pres_pl(:,k,:)* wt_ndg_weight_pl(:,:)
          rtau_qv(:,k,:)     = rtau_qv(:,k,:)     * wt_ndg_weight(:,:)
          rtau_qv_pl(:,k,:)  = rtau_qv_pl(:,k,:)  * wt_ndg_weight_pl(:,:)
       end do
       deallocate(wt_ndg_weight)
       deallocate(wt_ndg_weight_pl)
    end if
    !
    write(ADM_LOG_FID,*) '********* Nudging factors (1.0/tau) **************************'
    write(ADM_LOG_FID,*) '      z      velocity      temp.     pressure      qv'
    do k=1, ADM_kall
       write(ADM_LOG_FID,'(f8.2,4e14.6)') GRD_gz(k), rtau_vel(1,k,1), rtau_tem(1,k,1), &
            rtau_pres(1,k,1), rtau_qv(1,k,1)
    end do
    rtau_tem(:,:,:)    = rtau_tem(:,:,:) * CNST_CV
    rtau_tem_pl(:,:,:) = rtau_tem_pl(:,:,:) * CNST_CV
    rtau_pres(:,:,:)    = rtau_pres(:,:,:) * CNST_CV / CNST_RAIR * VMTR_GSGAM2(:,:,:)
    rtau_pres_pl(:,:,:) = rtau_pres_pl(:,:,:) * CNST_CV / CNST_RAIR * VMTR_GSGAM2_pl(:,:,:)
    !
    !--- allocate variables    ! ADM_kall => kmax_in
    allocate(var1(ADM_gall,kmax_in,ADM_lall,varmax))
    allocate(var2(ADM_gall,kmax_in,ADM_lall,varmax))
    allocate(var_bs(ADM_gall,ADM_kall,ADM_lall,varmax))
    allocate(var_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,varmax))
    var1 = 0.0 
    var2 = 0.0
    var_bs = 0.0d0
    var_bs_pl = 0.0d0
    !
    allocate(vel_bs(ADM_gall,ADM_kall,ADM_lall,nvel_max))
    allocate(vel_bs_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL,nvel_max))
    vel_bs = 0.0d0
    vel_bs_pl = 0.0d0
    !
    lsteps_dy = 86400/int(TIME_DTL)
    ndg_intv = lsteps_dy/ndg_times_dy
    !
    cmonth = ndg_start_mt
    cyear = ndg_start_yr
    !
    rl = 1
    mtnum = 1
    !
    call calendar_daymo(dyinmt, cyear, cmonth)
    rlinmt = dyinmt*ndg_times_dy
    ! [Add] 12/02/07 T.Seiki 
    call calendar_dayyr(dyinyr, cyear)
    rlinyr = dyinyr*ndg_times_dy
    !
    idate(1) = cyear
    idate(2) = cmonth
    idate(3) = ndg_start_dy 
    idate(4) = ndg_start_hr
    idate(5) = ndg_start_mn
    idate(6) = ndg_start_sc
    call calendar_yh2ss(ss, idate)
    if( ss < ctime ) then
       write(ADM_LOG_FID,*) 'Error about nudging start time! STOP!!' 
       call ADM_proc_stop
    end if
    ndg_start_st = int((ss - ctime)/dtime)
    !
    idate(1) = cyear
    ! [mod] 12/02/07 T.Seiki
!!$ idate(2) = cmonth
    if( opt_yearly_data )then
       idate(2) = 1
    else
       idate(2) = cmonth
    end if
    idate(3) = 1
    idate(4:6) = 0
    call calendar_yh2ss(ss0, idate)
    if(ss < ss0) then
       write(*,*) 'ss < ss0! abort'
       call ADM_proc_stop
    end if
    rl = int((ss - ss0)/dble(86400/ndg_times_dy)) + 1
    !
    ! for debug
    ! if(ADM_prc_me == ADM_prc_run_master) then 
    !    write(*,'(a,I8)') 'Nudging start step = ', ndg_start_st
    !    write(*,*) 'start record no. =', rl
    ! end if
    !
    !
    if(ndg_start_st == 0) then
       ndg_flag = .true.
       if(ADM_prc_me == ADM_prc_run_master) then 
          write(*,*) 'Nudging starts!'
       end if
    else
       ndg_flag = .false.
    end if
    !
    return
    !
  end subroutine ndg_setup
  !---------------------------------------------------------------------------------
  subroutine ndg_nudging_uvtp(  &
       rho,       rho_pl,       &
       vx,        vx_pl,        &
       vy,        vy_pl,        &
       vz,        vz_pl,        &
       w,         w_pl,         &
       tem,       tem_pl,       &
       pre,       pre_pl,       &
!       frhog,     frhog_pl,     &! S.Iga deleted 090908
       frhogvx,   frhogvx_pl,   &
       frhogvy,   frhogvy_pl,   &
       frhogvz,   frhogvz_pl,   &
       frhogw,    frhogw_pl,    &
       frhoge,    frhoge_pl,    &
       frhogetot, frhogetot_pl, &
       out_tendency             &
       )
    !
    use mod_adm, only :        &
         ADM_prc_me,           &
         ADM_prc_pl,           &
         ADM_prc_run_master,   &
         ADM_kmin,             &
         ADM_kmax,             &
         ADM_gall,             &
         ADM_kall,             &
         ADM_KNONE,            &
         ADM_lall,             &
         ADM_GALL_PL,          &
         ADM_LALL_PL
    use mod_vmtr, only : &
         VMTR_GSGAM2,    &
         VMTR_GSGAM2_pl, &
         VMTR_GSGAM2H,   &
         VMTR_GSGAM2H_pl
    use mod_history, only : &
         history_in
    use mod_cnst, only : &
         CNST_CV,   &
         CNST_EGRAV
    use mod_gmtr, only : &
         GMTR_P_var,  &
         GMTR_P_var_pl, &
         GMTR_P_IX,   &
         GMTR_P_IY,   &
         GMTR_P_IZ,   &
         GMTR_P_JX,   &
         GMTR_P_JY,   &
         GMTR_P_JZ
    use mod_oprt, only : &
         OPRT_horizontalize_vec
    use mod_grd, only : &
         GRD_afac, &
         GRD_bfac
    !
    implicit none
    !
    real(8), intent(in) :: rho(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: rho_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: vz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: vz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: w(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: w_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: tem(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: tem_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(in) :: pre(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(in) :: pre_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
!    real(8), intent(inout) :: frhog(ADM_gall,ADM_kall,ADM_lall)
!    real(8), intent(inout) :: frhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogvx(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogvx_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogvy(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogvy_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogvz(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogvz_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogw(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogw_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhoge(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhoge_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8), intent(inout) :: frhogetot(ADM_gall,ADM_kall,ADM_lall)
    real(8), intent(inout) :: frhogetot_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    logical, intent(in) :: out_tendency
    !
    real(8) :: frhogvx_old(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: frhogvy_old(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: frhogvz_old(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: frhogw_old(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: frhoge_old(ADM_gall,ADM_kall,ADM_lall)
    !
    real(8) :: u(ADM_gall,ADM_kall)
    real(8) :: v(ADM_gall,ADM_kall)
    !
    real(8) :: rhog(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rhog_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    real(8) :: rho_h(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: rho_h_pl(ADM_GALL_PL,ADM_kall,ADM_LALL_PL)
    !
    integer :: l, k
    !
    !
    !
    if(.not. ndg_flag) return
    !
    rhog(:,:,:) = rho(:,:,:) * VMTR_GSGAM2(:,:,:)
    rhog_pl(:,:,:) = rho_pl(:,:,:) * VMTR_GSGAM2_pl(:,:,:)
    !
    frhogvx_old = frhogvx
    frhogvy_old = frhogvy
    frhogvz_old = frhogvz
    frhogw_old  = frhogw
    frhoge_old  = frhoge
    !
    if(w_ndg) then
       do l=1, ADM_lall 
          do k =ADM_kmin,ADM_kmax+1
             rho_h(:,k,l) = 0.5D0 *         &
                  ( GRD_afac(k) * rho(:,k  ,l) &
                  + GRD_bfac(k) * rho(:,k-1,l) &
                  ) 
          end do
          rho_h(:,ADM_kmin-1,l) = rho_h(:,ADM_kmin,l)
          vel_bs(:,:,l,nvel_w) = - var_bs(:,:,l,nW) / (rho_h(:,:,l)*CNST_EGRAV)
          vel_bs(:,ADM_kmin-1:ADM_kmin,l,nvel_w) = 0.0d0
          vel_bs(:,ADM_kmax+1,l,nvel_w) = 0.0d0
       end do
       if(ADM_prc_me==ADM_prc_pl) then
          do l=1, ADM_lall_pl
             do k =ADM_kmin,ADM_kmax+1
                rho_h_pl(:,k,l) = 0.5D0 *         &
                     ( GRD_afac(k) * rho_pl(:,k  ,l) &
                     + GRD_bfac(k) * rho_pl(:,k-1,l) &
                  )
             end do
             rho_h_pl(:,ADM_kmin-1,l) = rho_h_pl(:,ADM_kmin,l)
             vel_bs_pl(:,:,l,nvel_w) = - var_bs_pl(:,:,l,nW) / (rho_h_pl(:,:,l)*CNST_EGRAV)
             vel_bs_pl(:,ADM_kmin-1:ADM_kmin,l,nvel_w) = 0.0d0
             vel_bs_pl(:,ADM_kmax+1,l,nvel_w) = 0.0d0
          end do
       end if
    end if
    !
    !
    !--- nudge winds
    !
    frhogvx(:,:,:) = frhogvx(:,:,:) &
         + rtau_vel(:,:,:)*( vel_bs(:,:,:,nvel_vx) - vx(:,:,:) ) * rhog(:,:,:)
    frhogvy(:,:,:) = frhogvy(:,:,:) &
         + rtau_vel(:,:,:)*( vel_bs(:,:,:,nvel_vy) - vy(:,:,:) ) * rhog(:,:,:)
    frhogvz(:,:,:) = frhogvz(:,:,:) &
         + rtau_vel(:,:,:)*( vel_bs(:,:,:,nvel_vz) - vz(:,:,:) ) * rhog(:,:,:)
    if(w_ndg) then
       frhogw(:,:,:) = frhogw(:,:,:) &
            + rtau_vel(:,:,:)*( vel_bs(:,:,:,nvel_w) - w(:,:,:) ) &
            * rho_h(:,:,:) * VMTR_GSGAM2H(:,:,:)
    end if
    !
    if(ADM_prc_me == ADM_prc_pl) then
       frhogvx_pl(:,:,:) = frhogvx_pl(:,:,:) &
            + rtau_vel_pl(:,:,:)*( vel_bs_pl(:,:,:,nvel_vx) - vx_pl(:,:,:) ) * rhog_pl(:,:,:)
       frhogvy_pl(:,:,:) = frhogvy_pl(:,:,:) &
            + rtau_vel_pl(:,:,:)*( vel_bs_pl(:,:,:,nvel_vy) - vy_pl(:,:,:) ) * rhog_pl(:,:,:)
       frhogvz_pl(:,:,:) = frhogvz_pl(:,:,:) &
            + rtau_vel_pl(:,:,:)*( vel_bs_pl(:,:,:,nvel_vz) - vz_pl(:,:,:) ) * rhog_pl(:,:,:)
       if(w_ndg) then
          frhogw_pl(:,:,:) = frhogw_pl(:,:,:) &
               + rtau_vel_pl(:,:,:)*( vel_bs_pl(:,:,:,nvel_w) - w_pl(:,:,:) ) &
               * rho_h_pl(:,:,:) * VMTR_GSGAM2H_pl(:,:,:)
       end if
    end if
    !
    call OPRT_horizontalize_vec(  &
         frhogvx, frhogvx_pl,     & !--- inout
         frhogvy, frhogvy_pl,     & !--- inout
         frhogvz, frhogvz_pl)       !--- inout
    !
    !--- nudge energy
    if(tem_ndg) then
       frhoge(:,:,:) = frhoge(:,:,:) &
            + rtau_tem(:,:,:)*( var_bs(:,:,:,nT) - tem(:,:,:) ) * rhog(:,:,:)
       frhogetot(:,:,:) = frhogetot(:,:,:) &
            + rtau_tem(:,:,:)*( var_bs(:,:,:,nT) - tem(:,:,:) ) * rhog(:,:,:)
       !
       if(ADM_prc_me == ADM_prc_pl) then
          frhoge_pl(:,:,:) = frhoge_pl(:,:,:) &
               + rtau_tem_pl(:,:,:)*( var_bs_pl(:,:,:,nT) - tem_pl(:,:,:) ) * rhog_pl(:,:,:)
          frhogetot_pl(:,:,:) = frhogetot_pl(:,:,:) &
               + rtau_tem_pl(:,:,:)*( var_bs_pl(:,:,:,nT) - tem_pl(:,:,:) ) * rhog_pl(:,:,:)
       end if
    end if
    !
    if(pres_ndg) then
       frhoge(:,:,:) = frhoge(:,:,:) &
            + rtau_pres(:,:,:)*( var_bs(:,:,:,nP) - pre(:,:,:) ) 
       frhogetot(:,:,:) = frhogetot(:,:,:) &
            + rtau_pres(:,:,:)*( var_bs(:,:,:,nP) - pre(:,:,:) ) 
       !
       if(ADM_prc_me == ADM_prc_pl) then
          frhoge_pl(:,:,:) = frhoge_pl(:,:,:) &
               + rtau_pres_pl(:,:,:)*( var_bs_pl(:,:,:,nP) - pre_pl(:,:,:) )
          frhogetot_pl(:,:,:) = frhogetot_pl(:,:,:) &
               + rtau_pres_pl(:,:,:)*( var_bs_pl(:,:,:,nP) - pre_pl(:,:,:) )
       end if
    end if
    !
    !
    !--- output tendency
    if(out_tendency) then
       rhog(:,:,:) = 1.0d0/rhog(:,:,:)
       frhogvx_old = (frhogvx - frhogvx_old)*rhog
       frhogvy_old = (frhogvy - frhogvy_old)*rhog
       frhogvz_old = (frhogvz - frhogvz_old)*rhog
       if(w_ndg) frhogw_old = (frhogw - frhogw_old)*rhog
       frhoge_old  = (frhoge  - frhoge_old)*rhog / CNST_CV
       !
       do l=1, ADM_lall
          do k=1,ADM_kall
             u(:,k) &
                  = (frhogvx_old(:,k,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IX) &
                  + frhogvy_old(:,k,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IY) &
                  + frhogvz_old(:,k,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_IZ))
             v(:,k) &
                  = (frhogvx_old(:,k,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JX) &
                  + frhogvy_old(:,k,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JY) &
                  + frhogvz_old(:,k,l)*GMTR_P_var(:,ADM_KNONE,l,GMTR_P_JZ))
          end do
          !
          call history_in('ndg_du',u)
          call history_in('ndg_dv',v)
          call history_in('ndg_dtem',frhoge_old(:,:,l))
          if(w_ndg) then
             call history_in('ndg_dw',frhogw_old(:,:,l))
             call history_in('ndg_w',vel_bs(:,:,l,nvel_w))
          end if
          !
       end do
    end if
    !
    !
    return
    !
  end subroutine ndg_nudging_uvtp
  !---------------------------------------------------------------------------------
  subroutine ndg_nudging_q( &
       rhog,    &
       rhogq,   &
       dt       &
       )
    !
    use mod_adm, only : &
         ADM_gall_in,  &
         ADM_kall, &
         ADM_lall
    use mod_runconf, only : &
         TRC_VMAX, &
         I_QV, I_QC
    use mod_gtl, only : &
         GTL_clip_region
    use mod_history, only : &
         history_in
    !
    implicit none                       
    !
    real(8), intent(inout) :: rhog(ADM_gall_in,ADM_kall,ADM_lall)
    real(8), intent(inout) :: rhogq(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)
    real(8), intent(in) :: dt
    !
    real(8) :: q(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)
    real(8) :: q_bs(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: drhogq(ADM_gall_in,ADM_kall,ADM_lall)
    real(8) :: rrhog(ADM_gall_in,ADM_kall,ADM_lall)
    integer :: l, n
    !
    rrhog(:,:,:) = 1.0d0/rhog(:,:,:)
    !
    do n=1, TRC_VMAX
       q(:,:,:,n) = rhogq(:,:,:,n)*rrhog(:,:,:)
    end do
    !
    rrhog(:,:,:) = rrhog(:,:,:)/dt
    !
    !--- nudge specific humidity        
    if(qv_ndg) then
       q_bs = 0.0d0
       call GTL_clip_region(var_bs(:,:,:,nQV),q_bs,1,ADM_kall)
       drhogq(:,:,:) = rtau_qv(:,:,:)*dt*( q_bs(:,:,:) - q(:,:,:,I_QV) ) * rhog(:,:,:)
       rhogq(:,:,:,I_QV) = rhogq(:,:,:,I_QV) + drhogq(:,:,:)
       drhogq(:,:,:) = drhogq(:,:,:) - min(0.0d0,rhogq(:,:,:,I_QV))
       rhogq(:,:,:,I_QV) = max(0.0d0,rhogq(:,:,:,I_QV))
       rhog(:,:,:) = rhog(:,:,:) + drhogq(:,:,:)
       !
       do l=1, ADM_lall
          call history_in('ndg_dqv',drhogq(:,:,l)*rrhog(:,:,l))
       end do
       !
    end if
    !
    return
    !
  end subroutine ndg_nudging_q
  !---------------------------------------------------------------------------------
  subroutine ndg_update_var
    !
    use mod_adm, only :     &
         ADM_gall,          & 
         ADM_kall,          &
         ADM_lall,          &
         ADM_kmin,          &
         ADM_kmax,          &
         ADM_KNONE,         &
         ADM_prc_me,        &
         ADM_prc_pl,        &
         ADM_GALL_PL,       &
         ADM_LALL_PL,       &
         ADM_prc_run_master
    use mod_comm, only : &
         comm_var
    use mod_oprt, only : &
         OPRT_horizontalize_vec
    use mod_gmtr, only : &
         GMTR_P_var,    &
         GMTR_P_var_pl, &
         GMTR_P_LAT,  &
         GMTR_P_IX,   &
         GMTR_P_IY,   &
         GMTR_P_IZ,   &
         GMTR_P_JX,   &
         GMTR_P_JY,   &
         GMTR_P_JZ
    use mod_time, only : &
         TIME_CSTEP
    use mod_history, only : &
         history_in
    !
    implicit none
    !
    real(8) :: f
    integer :: n, k, l
    integer :: cstep
    !
    cstep = TIME_CSTEP + 1
    if(.not.ndg_flag) then
       if(cstep == ndg_start_st) then
          ndg_flag = .true.
          if(ADM_prc_me == ADM_prc_run_master) then 
             write(*,*) 'Nudging starts!'
          end if
       end if
    end if
    !
    if(.not. ndg_flag) return
    !
    if(rd1st) then
       call set_rl_mtnum
       do n=1, varmax
          call read_var(mtnum,n,var1(:,:,:,n))
       end do
       rl = rl + 1
       !
       call set_rl_mtnum
       do n=1, varmax
          call read_var(mtnum,n,var2(:,:,:,n))
       end do
       rl = rl + 1
       !
       stepv1 = cstep
       stepv2 = cstep + ndg_intv
       !
       rd1st = .false.
    else
       call set_rl_mtnum
       if(mod((cstep-ndg_start_st), ndg_intv) == 0) then
          var1    = var2
          do n=1, varmax
             call read_var(mtnum,n,var2(:,:,:,n))
          end do 
          rl = rl + 1
          !
          stepv1 = stepv2
          stepv2 = stepv2 + ndg_intv
       end if
    end if
    !
    f = dble(cstep - stepv1)/dble(ndg_intv)
    !
    var_bs(:,ADM_kmin:ADM_kmax,:,:) =  &
         var1(:,1:kmax_in,:,:)*(1.0d0 - f) + var2(:,1:kmax_in,:,:)*f
    var_bs(:,ADM_kmin-1,:,:) = var_bs(:,ADM_kmin,:,:)
    var_bs(:,ADM_kmax+1,:,:) = var_bs(:,ADM_kmax,:,:) 
    !
    call comm_var(         &
         var_bs,var_bs_pl, &
         ADM_kall,         &
         varmax,           &
         comm_type = 1,    &
         NSval_fix=.true.  &
         )
    !
    if( trim(VEL_TYPE) == 'UV') then
       do l=1, ADM_lall
          call history_in('ndg_u',var_bs(:,:,l,nU))
          call history_in('ndg_v',var_bs(:,:,l,nV))
       end do
    end if
    !
    if(tem_ndg) then
       do l=1, ADM_lall
          call history_in('ndg_tem',var_bs(:,:,l,nT))
       end do
    end if
    !
    !
    if(trim(VEL_TYPE) == 'UV') then
       do k=1, ADM_kall
          vel_bs(:,k,:,nvel_vx) = var_bs(:,k,:,nU)*GMTR_P_var(:,ADM_KNONE,:,GMTR_P_IX) &
               + var_bs(:,k,:,nV)*GMTR_P_var(:,ADM_KNONE,:,GMTR_P_JX) 
          vel_bs(:,k,:,nvel_vy) = var_bs(:,k,:,nU)*GMTR_P_var(:,ADM_KNONE,:,GMTR_P_IY) &
               + var_bs(:,k,:,nV)*GMTR_P_var(:,ADM_KNONE,:,GMTR_P_JY) 
          vel_bs(:,k,:,nvel_vz) = var_bs(:,k,:,nU)*GMTR_P_var(:,ADM_KNONE,:,GMTR_P_IZ) &
               + var_bs(:,k,:,nV)*GMTR_P_var(:,ADM_KNONE,:,GMTR_P_JZ)
       end do
       if(ADM_prc_me == ADM_prc_pl) then
          do k=1, ADM_kall
             vel_bs_pl(:,k,:,nvel_vx) = var_bs_pl(:,k,:,nU)*GMTR_P_var_pl(:,ADM_KNONE,:,GMTR_P_IX) &
                  + var_bs_pl(:,k,:,nV)*GMTR_P_var_pl(:,ADM_KNONE,:,GMTR_P_JX) 
             vel_bs_pl(:,k,:,nvel_vy) = var_bs_pl(:,k,:,nU)*GMTR_P_var_pl(:,ADM_KNONE,:,GMTR_P_IY) &
                  + var_bs_pl(:,k,:,nV)*GMTR_P_var_pl(:,ADM_KNONE,:,GMTR_P_JY) 
             vel_bs_pl(:,k,:,nvel_vz) = var_bs_pl(:,k,:,nU)*GMTR_P_var_pl(:,ADM_KNONE,:,GMTR_P_IZ) &
                  + var_bs_pl(:,k,:,nV)*GMTR_P_var_pl(:,ADM_KNONE,:,GMTR_P_JZ)
          end do
       end if
    else if(trim(VEL_TYPE) == 'VXVYVZ') then
       vel_bs(:,:,:,nvel_vx)    = var_bs(:,:,:,nVX)
       vel_bs_pl(:,:,:,nvel_vx) = var_bs_pl(:,:,:,nVX)
       vel_bs(:,:,:,nvel_vy)    = var_bs(:,:,:,nVY)
       vel_bs_pl(:,:,:,nvel_vy) = var_bs_pl(:,:,:,nVY)
       vel_bs(:,:,:,nvel_vz)    = var_bs(:,:,:,nVZ)
       vel_bs_pl(:,:,:,nvel_vz) = var_bs_pl(:,:,:,nVZ)
    end if
    !
    call OPRT_horizontalize_vec(  &
         vel_bs(:,:,:,nvel_vx), vel_bs_pl(:,:,:,nvel_vx),     & !--- inout
         vel_bs(:,:,:,nvel_vy), vel_bs_pl(:,:,:,nvel_vy),     & !--- inout
         vel_bs(:,:,:,nvel_vz), vel_bs_pl(:,:,:,nvel_vz))       !--- inout
    !
    !
    return
    !
  end subroutine ndg_update_var
  !--------------------------------------------------------------------------------
  subroutine read_var(mo, nvar, var)
    !
    use mod_misc, only :          &
         MISC_get_available_fid,  &
         MISC_make_idstr
    use mod_adm, only :           &
         ADM_prc_me,              &
         ADM_prc_tab,             &
         ADM_gall,                &
         ADM_kall,                &
         ADM_kmin,                &
         ADM_kmax,                &
         ADM_lall
    !
    implicit none
    !
    integer, intent(in) :: mo
    integer, intent(in) :: nvar
    real(8), intent(out) :: var(ADM_gall,kmax_in,ADM_lall)
    !
    real(4) :: tmp(ADM_gall,kmax_in)
    integer(2) :: tmp2(ADM_gall,kmax_in)
    character(128) :: fname
    integer :: l
    integer :: rgnid
    integer :: fid
    !
    var = 0.0d0
    !
    if(.not. input_short_data) then
       do l=1, ADM_lall
          rgnid = ADM_prc_tab(l, ADM_prc_me)
          call MISC_make_idstr(fname, trim(basenames(mo,nvar)), 'rgn', rgnid)
          fid = MISC_get_available_fid()
          open(fid, file=trim(fname), form='unformatted', access='direct',  &
               recl=4*ADM_gall*kmax_in, status='old')
          !
          read(fid,rec=rl) tmp(:,:)
          var(:,:,l) = dble(tmp(:,:))
          !
          close(fid)
       end do
    else
       do l=1, ADM_lall
          rgnid = ADM_prc_tab(l, ADM_prc_me)
          call MISC_make_idstr(fname, trim(basenames(mo,nvar)), 'rgn', rgnid)
          fid = MISC_get_available_fid()
          open(fid, file=trim(fname), form='unformatted', access='direct',  &
               recl=2*ADM_gall*kmax_in, status='old')
          !
          read(fid,rec=rl) tmp2(:,:)
          var(:,:,l) = tmp2(:,:) * scale_factor(mo,nvar) + offset(mo,nvar)
          !
          close(fid)
       end do
    end if
    !
    !
    return
    !
  end subroutine read_var
  !--------------------------------------------------------------------------------
  subroutine set_rl_mtnum
    !
    use mod_calendar, only : &
         calendar_daymo,     &
         calendar_dayyr   ! [Add] 12/02/07 T.Seiki
    !
    implicit none
    !
    integer :: dyinmt
    integer :: dyinyr  ! [Add] 12/02/07 T.Seiki
    !
    ! [Mod] 12/02/07 T.Seiki
    if( .not. opt_yearly_data )then 
       if(rl > rlinmt) then
          rl = 1
          mtnum = mtnum + 1
          ! 
          cmonth = cmonth + 1
          if(cmonth > 12) then
             cmonth = 1
             if(.not. periodic) cyear = cyear + 1
          end if
          !
          call calendar_daymo(dyinmt, cyear, cmonth)
          !
          rlinmt = dyinmt*ndg_times_dy
          !
       end if
!!! ! [Add] 12/02/07 T.Seiki
    else
       if(rl > rlinyr) then
          rl = 1
          mtnum = mtnum + 1 ! here, mtnum is "yaer" number
          ! 
          cyear = cyear + 1
          !
          call calendar_dayyr(dyinyr, cyear)
          !
          rlinyr = dyinyr*ndg_times_dy
          !
       end if
    end if
    !
    return
    !
  end subroutine set_rl_mtnum
!-----------------------------------------------------------------------
   subroutine get_wt_ndg_weight(weight, weight_pl)
!
!  2008/09/10 [Add] M.Hara
!-----------------------------------------------------------------------
   use mod_adm, only : ADM_gall, ADM_lall, ADM_gall_pl, ADM_lall_pl
   use mod_cnst, only : CNST_ERADIUS
   use mod_gmtr, only : GMTR_lon, GMTR_lat, GMTR_lon_pl, GMTR_lat_pl
   implicit none

   real(8), intent(out) :: weight(ADM_gall,ADM_lall)
   real(8), intent(out) :: weight_pl(ADM_gall_pl,ADM_lall_pl)

   real(8) :: dist
   integer :: l, g

   do l =  1, ADM_lall
      do g = 1, ADM_gall
         call get_great_circle_dist(CNST_ERADIUS, &
&                                   wt_ndg_lon_c, wt_ndg_lat_c, &
&                                   GMTR_lon(g,l), GMTR_lat(g,l), &
&                                   dist)
         if (dist < wt_ndg_halo1) then
            weight(g,l) = wt_ndg_min
         else if (dist >= wt_ndg_halo1 .and. dist <= wt_ndg_halo2) then
            weight(g,l) = wt_ndg_min + &
&              ((dist-wt_ndg_halo1)*(wt_ndg_max-wt_ndg_min)/(wt_ndg_halo2-wt_ndg_halo1))
         else if (dist > wt_ndg_halo2) then
            weight(g,l) = wt_ndg_max
         end if
      end do
   end do

   do l =  1, ADM_lall_pl
      do g = 1, ADM_gall_pl
         call get_great_circle_dist(CNST_ERADIUS, &
&                                   wt_ndg_lon_c, wt_ndg_lat_c, &
&                                   GMTR_lon_pl(g,l), GMTR_lat_pl(g,l), &
&                                   dist)

         if (dist < wt_ndg_halo1) then
            weight_pl(g,l) = wt_ndg_min
         else if (dist >= wt_ndg_halo1 .and. dist <= wt_ndg_halo2) then
            weight_pl(g,l) = wt_ndg_min + &
&              ((dist-wt_ndg_halo1)*(wt_ndg_max-wt_ndg_min)/(wt_ndg_halo2-wt_ndg_halo1))
         else if (dist > wt_ndg_halo2) then
            weight_pl(g,l) = wt_ndg_max
         end if
      end do
   end do

   return
   end subroutine get_wt_ndg_weight
!-----------------------------------------------------------------------
   subroutine get_great_circle_dist(r, lon1, lat1, lon2, lat2, dist)
!  The formuation is Vincentry (1975)
!  http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
!
!  2008/09/10 [Add] M.Hara
!-----------------------------------------------------------------------
   USE mod_cnst, only : CNST_PI
   implicit none
   real(kind=8), intent(in)  :: r ! radius in meter
   real(kind=8), intent(in)  :: lon1, lat1 ! in degree
   real(kind=8), intent(in)  :: lon2, lat2 ! in degree
   real(kind=8), intent(out) :: dist ! distance of the two points in meter

   real(kind=8) :: rlon1, rlat1 ! in radian
   real(kind=8) :: rlon2, rlat2 ! in radian
   real(kind=8) :: drlon ! in radian

   rlon1 = lon1*CNST_PI/180.0d0
   rlat1 = lat1*CNST_PI/180.0d0
   rlon2 = lon2*CNST_PI/180.0d0
   rlat2 = lat2*CNST_PI/180.0d0
   drlon  = abs(rlon2-rlon1)

   dist = r*atan2( &
&     sqrt( &
&           ((cos(rlat2)*sin(drlon))**2)+ &
&           (((cos(rlat1)*sin(rlat2))-(sin(rlat1)*cos(rlat2)*cos(drlon)))**2) ), &
&           (sin(rlat1)*sin(rlat2))+(cos(rlat1)*cos(rlat2)*cos(drlon)))

   return
   end subroutine get_great_circle_dist
  !----------------------------------------------------------------------------------
end module mod_ndg
!------------------------------------------------------------------------------------
