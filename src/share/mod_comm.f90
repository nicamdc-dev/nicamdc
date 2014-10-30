!-------------------------------------------------------------------------------
!
!+  communication module
!
!-------------------------------------------------------------------------------
module mod_comm
  !-----------------------------------------------------------------------------
  !
  !++ description:
  !       this module is for the communication based on mpi library.
  !
  !++ Current Corresponding Author : K.Goto, H.Tomita
  !
  !++ History:
  !      Version    Date      Comment
  !      -----------------------------------------------------------------------
  !      0.00       04-02-17  Imported from igdc-4.33
  !                 06-09-??  K.Goto bug fix
  !                 06-10-08  S.Iga  add namelist (&COMMPARAM  max_varmax)
  !                 07-11-07  T.Mitsui add varmax check option(opt_check_varmax)
  !                 09-03-10  H.Tomita : Transplanting COMM_data_transfer2 from
  !                                      mod[mod_varcomm].
  !                 09-03-10  H.Tomita : rename COMM_data_transfer2 to COMM_var.
  !                 09-09-17  S.Iga : Add debug option and barrier option
  !                 10-06-07  S.Iga: new grid is implemented
  !                              (only the attribute of max_comm_xxx and
  !                              max_comm is changed from parameter to variable)
  !                 11-01-24  C.Kodama: Reduce memory usage in large rlevel.
  !                              (provided by Terai-san @ RIKEN)
  !                              Modified line: (20101207 teraim)
  !                 11-04-26  C.Kodama: default value of opt_check_varmax is changed to .true.
  !                                     and modify its behavior to abort when cmax exceeds max_maxvar*ADM_kall.
  !                 11-05-06  Y.Yamada: Merge tuning code with original code
  !                              (provided by Yamamoto-san @ NEC)
  !                              Modified line: !=org=
  !                 11-07-21  T.Ohno: A public variable 'comm_pl' is added.
  !                           If 'comm_pl' is false, pole data is not used in
  !                           COMM_data_transfer and COMM_var.
  !                 11-11-30  S.Iga (commit) : Modification around COMM_var,
  !                              suggested and modified by T.Inoue on 11-10-24
  !                 11-12-14  T.Seiki : allocatable variables are not permitted in type structure.
  !                              allocatable => pointer  (only @ SR16000 and ES)
  !                 12-03-26  T.Seiki : bug-fix if opt_comm_dbg=.true.
  !                 12-06-27  T.Ohno : bug fix for simulations at which
  !                           'comm_pl' is false
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_vlink_nmax
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ public procedure
  !
  public :: COMM_setup
  public :: COMM_data_transfer
  public :: COMM_data_transfer_nopl ! T.Ohno 110721
  public :: COMM_var
  public :: COMM_Stat_sum
  public :: COMM_Stat_sum_eachlayer
  public :: COMM_Stat_avg
  public :: COMM_Stat_max
  public :: COMM_Stat_min

  !-----------------------------------------------------------------------------


  !
  !++ private parameters & variables
  !
  integer,parameter,private::max_comm_r2r=9
!  integer,parameter,private::max_comm_r2p=ADM_vlink_nmax*2  !S.Iga100607 del
!  integer,parameter,private::max_comm_p2r=ADM_vlink_nmax*2  !S.Iga100607 del
!  integer,parameter,private::max_comm=max_comm_r2r+max_comm_r2p+max_comm_p2r !S.Iga100607 del
  integer,private,save::max_comm_r2p!S.Iga100607
  integer,private,save::max_comm_p2r!S.Iga100607
  integer,private,save::max_comm!S.Iga100607
!  integer,parameter,private::max_varmax=32
  integer,save,private::max_varmax=100 ! Iga(061008)
  logical,save,private::opt_check_varmax = .true. ! T.Mitsui 07/11/07  ! [mod] C.Kodama .false. -> .true.
  real(8),save,private::diag_varmax=0.d0           ! T.Mitsui 07/11/07
  !
  logical,save,private::opt_comm_dbg = .false.     ! S.Iga 09/09/XX
  real(8),save,private::dbg_sendbuf_init           ! S.Iga 09/09/XX
  real(8),save,private::dbg_recvbuf_init           ! S.Iga 09/09/XX
  logical,save,private::opt_comm_barrier = .false. ! S.Iga 09/09/XX
  integer,save,private,allocatable::dbg_areq_save(:,:) ! S.Iga 09/09/XX
  integer,save,private::dbg_tcount = 1 ! count comm_data_transfer is called  S.Iga 09/09/XX
  !
  !
  integer,parameter,private::ptr_prcid=1
  integer,parameter,private::ptr_lrgnid=2
  !
  integer,parameter,private::elemsize_comm=3
  integer,parameter,private::SIZE_COMM=1
  integer,parameter,private::LRGNID_COMM=2
  integer,parameter,private::BASE_COMM=3
  !
  integer,parameter,private::elemsize_copy=3
  integer,parameter,private::SIZE_COPY=1
  integer,parameter,private::LRGNID_COPY=2
  integer,parameter,private::SRC_LRGNID_COPY=3
  !
  !--------------------------------------------------
  integer,private,save::rank_me
  integer,private,save::max_comm_prc
!  integer,private,save::maxdatasize
  integer,private,save::maxdatasize_s
  integer,private,save::maxdatasize_r
  !
  integer,private,save::maxn
  integer,save,private::maxm
  integer,save,private::maxl
  !----
  integer,private,save::maxn_pl
  integer,save,private::maxm_pl
  integer,save,private::maxl_pl
  !----
  integer,private,save::maxn_r2r
  integer,save,private::maxm_r2r
  integer,save,private::maxl_r2r
  !----
  integer,private,save::maxn_r2p
  integer,save,private::maxm_r2p
  integer,save,private::maxl_r2p
  !----
  integer,private,save::maxn_p2r
  integer,save,private::maxm_p2r
  integer,save,private::maxl_p2r
  !----
  integer,private,save::maxn_sgp
  integer,save,private::maxm_sgp
  integer,save,private::maxl_sgp
  !
  integer,allocatable,private,save::prc_tab_rev(:,:)
  !
  integer,allocatable,private,save::clist(:)
  !
  !--------------------------------------------------
  !  for send
  !--------------------------------------------------
  integer,allocatable,private,save::nsmax(:,:)
  integer,allocatable,public,save::sendinfo(:,:,:,:)
  integer,allocatable,public,save::sendlist(:,:,:,:)
  integer,allocatable,private,save::nsmax_pl(:,:)
  integer,allocatable,private,save::sendinfo_pl(:,:,:,:)
  integer,allocatable,private,save::sendlist_pl(:,:,:,:)
  !
  !--------------------------------------------------
  !  for copy
  !--------------------------------------------------
  integer,allocatable,private,save::ncmax_r2r(:)
  integer,allocatable,private,save::copyinfo_r2r(:,:,:)
  integer,allocatable,private,save::recvlist_r2r(:,:,:)
  integer,allocatable,private,save::sendlist_r2r(:,:,:)
  !--------------------------------------------------
  integer,allocatable,private,save::ncmax_r2p(:)
  integer,allocatable,private,save::copyinfo_r2p(:,:,:)
  integer,allocatable,private,save::recvlist_r2p(:,:,:)
  integer,allocatable,private,save::sendlist_r2p(:,:,:)
  !--------------------------------------------------
  integer,allocatable,private,save::ncmax_p2r(:)
  integer,allocatable,private,save::copyinfo_p2r(:,:,:)
  integer,allocatable,private,save::recvlist_p2r(:,:,:)
  integer,allocatable,private,save::sendlist_p2r(:,:,:)
  !--------------------------------------------------
  integer,allocatable,private,save::ncmax_sgp(:)
  integer,allocatable,private,save::copyinfo_sgp(:,:,:)
  integer,allocatable,private,save::recvlist_sgp(:,:,:)
  integer,allocatable,private,save::sendlist_sgp(:,:,:)
  !--------------------------------------------------
  !
  !--------------------------------------------------
  !  for recv
  !--------------------------------------------------
  integer,allocatable,private,save::nrmax(:,:)
  integer,allocatable,public,save::recvinfo(:,:,:,:)
  integer,allocatable,public,save::recvlist(:,:,:,:)
  integer,allocatable,private,save::lrmax_pl(:)
  integer,allocatable,private,save::nrmax_pl(:,:)
  integer,allocatable,private,save::recvinfo_pl(:,:,:,:)
  integer,allocatable,private,save::recvlist_pl(:,:,:,:)
  !--------------------------------------------------
  !
  integer,allocatable,private::temp_sendorder(:,:)
  integer,allocatable,private::temp_recvorder(:,:)
  integer,allocatable,private::temp_dest_rgn(:,:,:)
  integer,allocatable,private::temp_src_rgn(:,:,:)
  integer,allocatable,private::temp_dest_rgn_pl(:,:,:)
  integer,allocatable,private::temp_src_rgn_pl(:,:,:)
  !integer,allocatable,private::temp_sb(:,:,:) !(20101207)removed by teraim
  integer,allocatable,private::tsb(:)
  !
  integer,allocatable,private,save::ssize(:,:)
  integer,allocatable,private,save::sendtag(:,:)
  integer,allocatable,private,save::somax(:)
  integer,allocatable,private,save::destrank(:,:)
  real(8),allocatable,public,save::sendbuf(:,:)
  integer,allocatable,private,save::rsize(:,:)
  integer,allocatable,private,save::recvtag(:,:)
  integer,allocatable,private,save::romax(:)
  integer,allocatable,private,save::sourcerank(:,:)
  real(8),allocatable,public,save::recvbuf(:,:)
  !
  integer,allocatable,private,save::n_nspl(:,:)
  !
  integer,allocatable,private,save::n_hemisphere_copy(:,:,:)
  integer,allocatable,private,save::s_hemisphere_copy(:,:,:)
  !
  !--------------------------------------------------
  integer,allocatable,private,save::tempbuf2D(:,:)
  integer,allocatable,private,save::tempbuf3D(:,:,:)
  integer,allocatable,private,save::tempbuf4D(:,:,:,:)
  integer,allocatable,private,save::tempbuf3D2(:,:,:)

  integer,allocatable,private,save::rsize_r2r(:,:,:)
  integer,allocatable,private,save::ssize_r2r(:,:,:)
  integer,allocatable,private,save::sourceid_r2r(:,:,:)
  integer,allocatable,private,save::destid_r2r(:,:,:)
  !integer,allocatable,private,save::mrecv_r2r(:,:,:) !(20101207)removed by teraim
  integer,allocatable,private,save::msend_r2r(:,:,:)
  integer,allocatable,private,save::maxcommrecv_r2r(:,:)
  integer,allocatable,private,save::maxcommsend_r2r(:,:)
  !integer,allocatable,private,save::recvtag_r2r(:,:,:) !(20101207)removed by teraim
  !integer,allocatable,private,save::sendtag_r2r(:,:,:) !(20101207)removed by teraim
  integer,allocatable,private,save::rlist_r2r(:,:,:,:)
  integer,allocatable,private,save::qlist_r2r(:,:,:,:)
  integer,allocatable,private,save::slist_r2r(:,:,:,:)
  integer,private,save::max_datasize_r2r
  real(8),allocatable,private::recvbuf_r2r(:,:,:)
  real(8),allocatable,private::sendbuf_r2r(:,:,:)
  !
  integer,allocatable,private,save::rsize_r2p(:,:,:)
  integer,allocatable,private,save::ssize_r2p(:,:,:)
  integer,allocatable,private,save::source_prc_r2p(:,:,:)
  integer,allocatable,private,save::source_rgn_r2p(:,:,:)
  integer,allocatable,private,save::dest_prc_r2p(:,:,:)
  integer,allocatable,private,save::maxcommrecv_r2p(:,:)
  integer,allocatable,private,save::maxcommsend_r2p(:,:)
  integer,allocatable,private,save::recvtag_r2p(:,:,:)
  integer,allocatable,private,save::sendtag_r2p(:,:,:)
  integer,allocatable,private,save::rlist_r2p(:,:,:,:)
  integer,allocatable,private,save::qlist_r2p(:,:,:,:)
  integer,allocatable,private,save::slist_r2p(:,:,:,:)
  integer,save,private :: max_datasize_r2p
  real(8),allocatable,private::recvbuf_r2p(:,:,:)
  real(8),allocatable,private::sendbuf_r2p(:,:,:)
  !
  integer,allocatable,private,save::recvtag_p2r(:,:)
  integer,allocatable,private,save::sendtag_p2r(:,:)
  real(8),allocatable,private::sendbuf_p2r(:,:)
  real(8),allocatable,private::recvbuf_p2r(:,:)
  !----------------------------------------------
  integer,allocatable,private,save::dest_rank_all(:,:,:)
  integer,allocatable,private,save::src_rank_all(:,:,:)
  !----------------------------------------------

!!  real(8),allocatable,public::  comm_dbg_recvbuf(:,:,:) !iga
!!  real(8),allocatable,public::  comm_dbg_sendbuf(:,:,:)  !iga



  !
  integer,allocatable,private,save::imin(:),imax(:) &
       ,jmin(:),jmax(:) &
       ,gmin(:),gmax(:) &
       ,gall(:)
  !
  integer,allocatable,private,save::nmin_nspl(:),nmax_nspl(:) &
       ,pmin_nspl(:),pmax_nspl(:) &
       ,lmin_nspl(:),lmax_nspl(:) &
       ,gmin_nspl(:),gmax_nspl(:) &
       ,gall_nspl(:)
  integer,allocatable,private,save::pl_index(:,:,:,:)
  !
  !----------------------------------------------------------------
  !++ public variablesizc
  !  integer,allocatable,public,save::izc(:,:,:)
  !  integer,allocatable,public,save::itc(:,:,:,:)
  !  integer,allocatable,public,save::itc2(:,:,:,:)
  !  integer,allocatable,public,save::ntr(:,:,:,:,:)
  !  integer,parameter,public::order=2
  !  integer,parameter,public::noupfi=(order+1)*(order+2)/2
  !----------------------------------------------------------------
  !
  integer,private,save::halomax
  integer,private,save::kmax

  integer,private,save :: comm_call_count=0
  real(8),private,save::time_total= 0.D0
  real(8),private,save::time_pre = 0.D0
  real(8),private,save::time_bar1= 0.D0
  real(8),private,save::time_bar2= 0.D0
  real(8),private,save::time_sbuf= 0.D0
  real(8),private,save::time_recv= 0.D0
  real(8),private,save::time_send= 0.D0
  real(8),private,save::time_copy= 0.D0
  real(8),private,save::time_wait= 0.D0
  real(8),private,save::time_rbuf= 0.D0
  real(8),private,save::time_copy_sgp= 0.D0
  real(8),private,save::size_total= 0.D0
  real(8),private,save::comm_count= 0.D0
  real(8),private::t(0:12)
  !
  !(20101207)added by teraim
  type type_tempsb
    integer::num
    ! 2011/12/14 [Mod] T.Seiki
!!$ integer,allocatable::col(:)
!!$ integer,allocatable::val(:)
    integer,pointer :: col(:)
    integer,pointer :: val(:)
  end type
  type(type_tempsb),allocatable::tempsb(:)
  integer,parameter::max_size=10 ! This is not optimal value.
  !
  logical, public, save :: comm_pl = .true. ! T.Ohno 110721

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !(20101207) added by teraim
  subroutine init_tempsb
    use mod_adm, only : ADM_rgn_nmax
    implicit none
    integer::i
    !
    allocate(tempsb(ADM_rgn_nmax+2))
    tempsb(:)%num=0
    !
    do i=1, ADM_rgn_nmax+2
      allocate(tempsb(i)%col(max_size))
      allocate(tempsb(i)%val(max_size))
      tempsb(i)%col(:)=-1
      tempsb(i)%val(:)=-1
    enddo
  end subroutine
  !(20101207) added by teraim
  subroutine finalize_tempsb
    use mod_adm, only : ADM_rgn_nmax
    implicit none
    integer::i
    !
    do i=1, ADM_rgn_nmax+2
      deallocate(tempsb(i)%col)
      deallocate(tempsb(i)%val)
    enddo
    deallocate(tempsb)
  end subroutine
  !(20101207) added by teraim
  subroutine add_tempsb(icol, irow, ival)
    implicit none
    integer,intent(in)::icol,irow,ival
    !
    if(ival > 0) then
      if(tempsb(irow)%num < max_size) then
        tempsb(irow)%num = tempsb(irow)%num + 1
        tempsb(irow)%col(tempsb(irow)%num)=icol
        tempsb(irow)%val(tempsb(irow)%num)=ival
      else
        write(*,*)"range of list is over."
        stop
      endif
    endif
  end subroutine
  !(20101207) added by teraim
  subroutine show_tempsb
    use mod_adm, only : ADM_rgn_nmax
    implicit none
    integer::i,j
    !
    write(*,*)"show"
    do i=1, ADM_rgn_nmax+2
      do j=1, max_size
        if(tempsb(i)%val(j) > 0) then
          write(*,*)"(i,j)=",i,j," : show_list=(",i,tempsb(i)%col(j),tempsb(i)%val(j),")"
        endif
      enddo
    enddo
  end subroutine
  !(20101207) added by teraim
  subroutine get_tempsb(icol, irow, ret)
    implicit none
    integer,intent(in)::icol,irow
    integer,intent(out)::ret
    integer::i
    !
    ret = 0
    do i=1, max_size
      if(tempsb(irow)%col(i) == icol) then
        ret = tempsb(irow)%val(i)
        exit
      endif
    enddo
  end subroutine
  !-----------------------------------------------------------------------------
  subroutine COMM_setup( &
       max_hallo_num,    & !--- IN : number of hallo regions
       debug             ) !--- IN : debug flag
    use mod_adm, only :    &
         !--- public parameters
         ADM_w,          &
         ADM_e,          &
         ADM_n,          &
         ADM_s,          &
         ADM_sw,         &
         ADM_nw,         &
         ADM_ne,         &
         ADM_se,         &
         ADM_rid,        &
         ADM_dir,        &
         ADM_vlink_nmax, &
         ADM_rgn_nmax_pl,&
         ADM_npl,        &
         ADM_spl,        &
         ADM_gslf_pl,    &
         ADM_prc_all,    &
         ADM_prc_rnum,   &
         ADM_prc_tab,    &
         ADM_prc_me,     &
         ADM_rgn_nmax,   &
         ADM_rgn_etab,   &
         ADM_rgn_vnum,   &
         ADM_rgn_vtab,   &
         ADM_rgn_vtab_pl,&
         ADM_gmin,       &
         ADM_gmax,       &
         ADM_gall_1d,    &
         ADM_lall,       &
         ADM_kall,       &
         ADM_prc_nspl,   &
         ADM_comm_world,&
         ADM_rgn2prc,& !(20101207)added by teraim
         !--- For namelist.
         ADM_CTL_FID,    &  !Iga(061008)
         ADM_LOG_FID        !Iga(061008)
    implicit none

    integer,intent(in),optional :: max_hallo_num
    logical,intent(in),optional :: debug

    integer :: i,j,l,n,m,p,q
    integer :: rgnid
    integer :: ierr
    !
    !integer :: nn,t,n1,n2,n3,n4,n5,n6
    !
    integer::lr,mr
    integer::ls,ms
    integer::nr,nc,ns
    integer::rs,cs,ss
    integer::nd,ld,pl,halo
    integer::in,jn
    integer::srgnid,rrgnid
    integer::ck
    integer::srank,drank
    integer::ro,so
    !
    integer::suf,g_1d
    suf(i,j,g_1d)=(g_1d)*((j)-1)+(i)
    !
    integer::rgnid1,rgnid2,ret !(20101207) added by teraim
    !
    ! Iga(061008) ==>
    namelist / COMMPARAM /   &
         max_varmax,         & ! max number of communication variables
         opt_check_varmax,   & ! check option of varmax [Add] T.Mitsui 07/11/07
         opt_comm_dbg,       & ! debug option of comm_data_transfer [Add] S.Iga 0909XX
         opt_comm_barrier      ! debug option of comm_data_transfer [Add] S.Iga 0909XX
    !
    max_comm_r2p=ADM_vlink_nmax*2!S.Iga100607
    max_comm_p2r=ADM_vlink_nmax*2!S.Iga100607
    max_comm=max_comm_r2r+max_comm_r2p+max_comm_p2r!S.Iga100607

    !--- < reading parameters > ---
    !
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=COMMPARAM,iostat=ierr)
    if(ierr<0) then
       write(ADM_LOG_FID,*) &
            'Msg : Sub[COMM_setup]/Mod[comm]'
       write(ADM_LOG_FID,*) &
            ' *** No namelist in paramter file.'
       write(ADM_LOG_FID,*) &
            ' *** Use default values.'
    else if(ierr>0) then
       write(*,*) &
            'Msg : Sub[COMM_setup]/Mod[comm]'
       write(*,*) &
            ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
    end if
    write(ADM_LOG_FID,COMMPARAM)
    ! <== Iga(061008)

    if (present(max_hallo_num)) then
       halomax=max_hallo_num
    else
       halomax=1
    endif
    !    if (present(max_k_num)) then
    !      kmax=max_k_num
    !    else
    !      kmax=ADM_kall
    !    endif
    kmax=ADM_kall
    !
    allocate(prc_tab_rev(ptr_prcid:ptr_lrgnid,ADM_rgn_nmax))
    !
    do p=1,ADM_prc_all
       do n=1,ADM_prc_rnum(p)
          prc_tab_rev(ptr_prcid,ADM_prc_tab(n,p))=p
          prc_tab_rev(ptr_lrgnid,ADM_prc_tab(n,p))=n
       end do
    end do

!    if (ADM_prc_me.eq.1)   write(*,*) 'ADM_prc_tab',ADM_prc_tab
!    if (ADM_prc_me.eq.1)   write(*,*) 'prc_tab_rev', prc_tab_rev

    if(ADM_prc_nspl(ADM_npl) < 0 .and. ADM_prc_nspl(ADM_spl) <0 ) comm_pl = .false. ! T.Ohno 110721

    !
    allocate(imin(halomax))
    allocate(imax(halomax))
    allocate(jmin(halomax))
    allocate(jmax(halomax))
    allocate(gmin(halomax))
    allocate(gmax(halomax))
    allocate(gall(halomax))
    !
    !(20101207) changed by teraim
    allocate(rsize_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(ssize_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(sourceid_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(destid_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    !allocate(mrecv_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    !allocate(msend_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    rgnid1=ADM_prc_tab(1,ADM_prc_me)
    rgnid2=ADM_prc_tab(ADM_prc_rnum(ADM_prc_me),ADM_prc_me)
    allocate(msend_r2r(ADM_rgn_nmax,halomax,rgnid1:rgnid2))
    allocate(maxcommrecv_r2r(halomax,ADM_rgn_nmax))
    allocate(maxcommsend_r2r(halomax,ADM_rgn_nmax))
    !allocate(recvtag_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    !allocate(sendtag_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    !
    imin(halomax)=(ADM_gmin-1)+halomax
    imax(halomax)=(ADM_gmax-1)+halomax
    jmin(halomax)=(ADM_gmin-1)+halomax
    jmax(halomax)=(ADM_gmax-1)+halomax
    gmin(halomax)=(ADM_gmin-1)+halomax
    gmax(halomax)=(ADM_gmax-1)+halomax
    gall(halomax)=(ADM_gall_1d-2)+2*halomax
    !
    max_datasize_r2r=(gmax(halomax)-gmin(halomax)+1)*halomax
    allocate(rlist_r2r(max_datasize_r2r,max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(qlist_r2r(max_datasize_r2r,max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(slist_r2r(max_datasize_r2r,max_comm_r2r,halomax,ADM_rgn_nmax))
    !
    allocate(recvbuf_r2r(max_datasize_r2r*kmax*max_varmax  &
         ,ADM_prc_rnum(ADM_prc_me),max_comm_r2r))
    allocate(sendbuf_r2r(max_datasize_r2r*kmax*max_varmax  &
         ,ADM_prc_rnum(ADM_prc_me),max_comm_r2r))

    allocate(n_hemisphere_copy(ADM_w:ADM_s,halomax,ADM_rgn_nmax))
    allocate(s_hemisphere_copy(ADM_w:ADM_s,halomax,ADM_rgn_nmax))

    allocate(tempbuf2D(halomax,ADM_prc_rnum(ADM_prc_me)))
    allocate(tempbuf3D(max_comm_r2r,halomax,ADM_prc_rnum(ADM_prc_me)))
    allocate(tempbuf4D(max_datasize_r2r,max_comm_r2r,halomax,ADM_prc_rnum(ADM_prc_me)))
    allocate(tempbuf3D2(ADM_w:ADM_s,halomax,ADM_prc_rnum(ADM_prc_me)))

    !
    rsize_r2r(:,:,:)=0
    ssize_r2r(:,:,:)=0
    sourceid_r2r(:,:,:)=-1
    destid_r2r(:,:,:)=-1
    !mrecv_r2r(:,:,:)=-1
    msend_r2r(:,:,:)=-1
    maxcommrecv_r2r(:,:)=max_comm_r2r
    maxcommsend_r2r(:,:)=max_comm_r2r
    !(20101207) removed by teraim
    !recvtag_r2r(:,:,:)=-1
    !sendtag_r2r(:,:,:)=-1
    !
    rlist_r2r(:,:,:,:)=-1
    qlist_r2r(:,:,:,:)=-1
    slist_r2r(:,:,:,:)=-1
    !
    n_hemisphere_copy(:,:,:)=0
    s_hemisphere_copy(:,:,:)=0
    !
    allocate(rsize_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(ssize_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(source_prc_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(source_rgn_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(dest_prc_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(maxcommrecv_r2p(ADM_npl:ADM_spl,halomax))
    allocate(maxcommsend_r2p(ADM_npl:ADM_spl,halomax))
    allocate(recvtag_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(sendtag_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    !
    max_datasize_r2p=halomax*(halomax+1)/2
    !
    allocate(rlist_r2p(max_datasize_r2p,max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(qlist_r2p(max_datasize_r2p,max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(slist_r2p(max_datasize_r2p,max_comm_r2p,ADM_npl:ADM_spl,halomax))
    !
    allocate(recvbuf_r2p(max_datasize_r2p*kmax*max_varmax &
         ,max_comm_r2p,ADM_npl:ADM_spl))
    allocate(sendbuf_r2p(max_datasize_r2p*kmax*max_varmax &
         ,max_comm_r2p,ADM_npl:ADM_spl))
    !
    rsize_r2p(:,:,:)=0
    ssize_r2p(:,:,:)=0
    source_prc_r2p(:,:,:)=-1
    source_rgn_r2p(:,:,:)=-1
    dest_prc_r2p(:,:,:)=-1
    maxcommrecv_r2p(:,:)=max_comm_r2p
    maxcommsend_r2p(:,:)=max_comm_r2p
    recvtag_r2p(:,:,:)=-1
    sendtag_r2p(:,:,:)=-1
    !
    rlist_r2p(:,:,:,:)=-1
    qlist_r2p(:,:,:,:)=-1
    slist_r2p(:,:,:,:)=-1
    !
!!!!!!!!!!!!!!!
    allocate(nmin_nspl(1:halomax))
    allocate(nmax_nspl(1:halomax))
    allocate(pmin_nspl(1:halomax))
    allocate(pmax_nspl(1:halomax))
    allocate(lmin_nspl(1:halomax))
    allocate(lmax_nspl(1:halomax))
    allocate(gmin_nspl(1:halomax))
    allocate(gmax_nspl(1:halomax))
    allocate(gall_nspl(1:halomax))
    nmin_nspl(halomax)=1
    nmax_nspl(halomax)=halomax+1
    pmin_nspl(halomax)=1
    pmax_nspl(halomax)=ADM_vlink_nmax
    lmin_nspl(halomax)=1
    lmax_nspl(halomax)=halomax
    gmin_nspl(halomax)=2
    gmax_nspl(halomax)=1+5*halomax*(halomax+1)/2
    gall_nspl(halomax)=1+5*halomax*(halomax+1)/2
    allocate(pl_index(nmin_nspl(halomax):nmax_nspl(halomax) &
         ,pmin_nspl(halomax):pmax_nspl(halomax) &
         ,lmin_nspl(halomax):lmax_nspl(halomax),halomax))
    !
    do halo=1,halomax
       !
       imin(halo)=(ADM_gmin-1)+halo
       imax(halo)=(ADM_gmax-1)+halo
       jmin(halo)=(ADM_gmin-1)+halo
       jmax(halo)=(ADM_gmax-1)+halo
       gmin(halo)=(ADM_gmin-1)+halo
       gmax(halo)=(ADM_gmax-1)+halo
       gall(halo)=(ADM_gall_1d-2)+2*halo
       !
       nmin_nspl(halo)=1
       nmax_nspl(halo)=halo+1
       pmin_nspl(halo)=1
       pmax_nspl(halo)=ADM_vlink_nmax
       lmin_nspl(halo)=1
       lmax_nspl(halo)=halo
       gmin_nspl(halo)=2
       gmax_nspl(halo)=1+5*halo*(halo+1)/2
       gall_nspl(halo)=1+5*halo*(halo+1)/2
       pl_index(:,:,:,halo)=-1
       do l=lmin_nspl(halo),lmax_nspl(halo)
          do p=pmin_nspl(halo),pmax_nspl(halo)
             do n=nmin_nspl(halo),l+1
                pl_index(n,p,l,halo)=n+(p-1)*l+(1+5*(l-1)*l/2)
             enddo
          enddo
          pl_index(l+1,pmax_nspl(halo),l,halo)=nmin_nspl(halo) &
               +(pmin_nspl(halo)-1)*l &
               +(1+5*(l-1)*l/2)
       enddo
       !
       ! --- r2p ----
       if(comm_pl)then
        do p=pmin_nspl(halo),pmax_nspl(halo)
           rsize_r2p(p,ADM_npl,halo)=halo*(halo+1)/2
           source_prc_r2p(p,ADM_npl,halo)=prc_tab_rev(ptr_prcid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_npl,p))
           source_rgn_r2p(p,ADM_npl,halo)=prc_tab_rev(ptr_lrgnid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_npl,p))
        enddo
        do p=pmin_nspl(halo),pmax_nspl(halo)
           q=0
           do ld=lmin_nspl(halo),lmax_nspl(halo)
              do nd=nmin_nspl(halo),ld
                 q=q+1
                 in=-nd+ld+nmin_nspl(halo)-lmin_nspl(halo)+imin(halo)
                 jn=-nd+nmin_nspl(halo)+(jmax(halo)-jmin(halo))+jmin(halo)
                 rlist_r2p(q,mod(p,ADM_vlink_nmax)+1,ADM_npl,halo) &
                      =pl_index(nd+1,p,ld,halo)
                 qlist_r2p(q,mod(p,ADM_vlink_nmax)+1,ADM_npl,halo)=suf(in,jn,gall(halo))
              enddo
           enddo
        enddo
        do p=pmin_nspl(halo),pmax_nspl(halo)
           rsize_r2p(p,ADM_spl,halo)=halo*(halo+1)/2
           source_prc_r2p(p,ADM_spl,halo)=prc_tab_rev(ptr_prcid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_spl,p))
           source_rgn_r2p(p,ADM_spl,halo)=prc_tab_rev(ptr_lrgnid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_spl,p))
        enddo
        do p=pmin_nspl(halo),pmax_nspl(halo)
           q=0
           do ld=lmin_nspl(halo),lmax_nspl(halo)
              do nd=nmin_nspl(halo),ld
                 q=q+1
                 in=nd-ld-nmin_nspl(halo) &
                      +lmin_nspl(halo)+(imax(halo)-imin(halo))+imin(halo)
                 jn=nd-nmin_nspl(halo)+jmin(halo)
                 rlist_r2p(q,p,ADM_spl,halo)=pl_index(nd,p,ld,halo)
                 qlist_r2p(q,p,ADM_spl,halo)=suf(in,jn,gall(halo))
              enddo
           enddo
        enddo
        maxcommrecv_r2p(ADM_npl,halo)=(pmax_nspl(halo)-pmin_nspl(halo)+1)
        maxcommrecv_r2p(ADM_spl,halo)=(pmax_nspl(halo)-pmin_nspl(halo)+1)
        !
        do pl=ADM_npl,ADM_spl
           do p=1,maxcommrecv_r2p(pl,halo)
              if (ADM_prc_me==source_prc_r2p(p,pl,halo)) then
                 dest_prc_r2p(p,pl,halo)=ADM_prc_nspl(pl)
                 ssize_r2p(p,pl,halo)=rsize_r2p(p,pl,halo)
                 do q=1,ssize_r2p(p,pl,halo)
                    slist_r2p(q,p,pl,halo)=qlist_r2p(q,p,pl,halo)
                 enddo
              endif
              sendtag_r2p(p,pl,halo)=pl+(ADM_spl-ADM_npl+1)*(p-1) &
                   +ADM_rgn_nmax**2+ADM_vlink_nmax*2
              recvtag_r2p(p,pl,halo)=sendtag_r2p(p,pl,halo)

!              write(*,*) 'sendtag_r2p',ADM_prc_me,p,pl,halo,sendtag_r2p(p,pl,halo)
           enddo
           maxcommsend_r2p(pl,halo)=(pmax_nspl(halo)-pmin_nspl(halo)+1)
        enddo
       endif
       !
       ! --- r2r ----
       do l=1,ADM_prc_rnum(ADM_prc_me)
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          m=0
          if (ADM_rgn_etab(ADM_dir,ADM_sw,rgnid)==ADM_ne) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_sw,rgnid)
                ! mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do j=jmin(halo)-halo,jmin(halo)-1
                   do i=imin(halo),imax(halo)
                      n=n+1
                      in=i
                      jn=j+jmax(halo)+1-jmin(halo)
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          if ((ADM_rgn_vnum(ADM_w,rgnid)==3)) then
             if (ADM_rgn_etab(ADM_dir,ADM_sw,rgnid)==ADM_se) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_sw,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)+1,imax(halo)+(j-(jmin(halo)-1))
                         n=n+1
                         in=j+jmax(halo)-2*jmin(halo)+imin(halo)+1
                         jn=-i+j+imin(halo)+jmax(halo)-jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                      do i=imin(halo)-halo,imin(halo)-1+(j-(jmin(halo)-1))
                         n=n+1
                         in=i+imax(halo)+1-imin(halo)
                         jn=j+jmax(halo)+1-jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
          endif
          !
          if (ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_ne) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_nw,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do i=imin(halo)-halo,imin(halo)-1
                   do j=jmin(halo)+(i-(imin(halo)-1)),jmax(halo)+(i-(imin(halo)-1))
                      n=n+1
                      in=i-j+imax(halo)-imin(halo)+jmin(halo)+1
                      jn=i+imax(halo)-2*imin(halo)+jmin(halo)+1
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          elseif (ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_se) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_nw,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do j=jmin(halo),jmax(halo)
                   do i=imin(halo)-halo,imin(halo)-1
                      n=n+1
                      in=i+imax(halo)-imin(halo)+1
                      jn=j
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          !
!!!!!
          if ((ADM_rgn_vnum(ADM_n,rgnid)==5)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2-1
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imin(halo)-halo,imin(halo)-1
                      do j=jmax(halo)+1+(i-(imin(halo)-1)) &
                           ,min(jmax(halo)+1,jmax(halo)+(halo-1)+(i-(imin(halo)-1)))
                         n=n+1
                         in=-j+jmax(halo)+imin(halo)+1
                         jn=i-j+imax(halo)-2*imin(halo)+jmax(halo)+jmin(halo)+2
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             if (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,3)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,3)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+2,jmax(halo)+halo
                      do i=imin(halo)+1,imin(halo)+1+(j-(jmax(halo)+2))
                         n=n+1
                         in=-i+j-2*jmax(halo)+jmin(halo)+imax(halo)+imin(halo)-1
                         jn=-i+imax(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
!!!!!
          !
          if (ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_nw) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_ne,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do j=jmax(halo)+1,jmax(halo)+halo
                   do i=imin(halo)+1+(j-(jmax(halo)+1)),imax(halo)+1+(j-(jmax(halo)+1))
                      n=n+1
                      in=j-jmax(halo)+imin(halo)-1
                      jn=-i+j+imax(halo)-jmax(halo)+jmin(halo)
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          elseif (ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_sw) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_ne,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                n=0
                do j=jmax(halo)+1,jmax(halo)+halo
                   do i=imin(halo),jmax(halo)
                      n=n+1
                      in=i
                      jn=j-jmax(halo)-1+jmin(halo)
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          !
          if (ADM_rgn_etab(ADM_dir,ADM_se,rgnid)==ADM_nw) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_se,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                n=0
                do j=jmin(halo),jmax(halo)
                   do i=imax(halo)+1,imax(halo)+halo
                      n=n+1
                      in=i-imax(halo)+imin(halo)-1
                      jn=j
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          if ((ADM_rgn_vnum(ADM_e,rgnid)==3)) then
             if (ADM_rgn_etab(ADM_dir,ADM_se,rgnid)==ADM_sw) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_se,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+halo
                      do j=jmin(halo)+1+(i-(imax(halo)+1)),jmax(halo)
                         n=n+1
                         in=i-j-imax(halo)+jmax(halo)+imin(halo)
                         jn=i-imax(halo)-1+jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                      do j=jmax(halo)+1+(i-(imax(halo)+1)),jmax(halo)+halo
                         n=n+1
                         in=i-imax(halo)-1+imin(halo)
                         jn=j-jmax(halo)-1+jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
          endif
          !
!!!!!
          if ((ADM_rgn_vnum(ADM_s,rgnid)==5)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+2,imax(halo)+halo
                      do j=jmin(halo)+1,jmin(halo)+1+(i-(imax(halo)+2))
                         n=n+1
                         in=-j+jmax(halo)+imin(halo)+1
                         jn=i-j-2*imax(halo)+imin(halo)+jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,3)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2-1
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,3)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imax(halo)+1+(j-(jmin(halo)-1)) &
                           ,min(imax(halo)+1,imax(halo)+(halo-1)+(j-(jmin(halo)-1)))
                         n=n+1
                         in=-i+j+imax(halo)+jmax(halo)+imin(halo)-2*jmin(halo)+2
                         jn=-i+imax(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !!
          !
          !!
          if ((ADM_rgn_vnum(ADM_w,rgnid)==4)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_w,rgnid,2)==ADM_n) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_w,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)-1+(j-(jmin(halo)-1)),imin(halo)-1
                         n=n+1
                         in=i-j+imax(halo)-imin(halo)-jmax(halo)+2*jmin(halo)
                         jn=i+imax(halo)-2*imin(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_w,rgnid,2)==ADM_e) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_w,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)-halo,imin(halo)-1
                         n=n+1
                         in=i+imax(halo)-imin(halo)+1
                         jn=j+jmax(halo)-jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_w,rgnid,2)==ADM_s) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_w,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)-halo,imin(halo)-halo+(j-(jmin(halo)-halo))
                         !do i=imin(halo)-halo,imin(halo)-1
                         !  do j=i+jmin(halo)-imin(halo),jmin(halo)-1
                         n=n+1
                         in=j+jmax(halo)+1-2*jmin(halo)+imin(halo)
                         jn=-i+j-imax(halo)+2*imin(halo)+jmax(halo)-jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             if (ADM_rgn_etab(ADM_dir,ADM_sw,rgnid)==ADM_se) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_sw,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)+(j-(jmin(halo)-1)),imax(halo)+(j-(jmin(halo)-1))
                         n=n+1
                         in=j+jmax(halo)-2*jmin(halo)+imin(halo)+1
                         jn=-i+j+imin(halo)+jmax(halo)-jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !!
          !
          !!
          if ((ADM_rgn_vnum(ADM_n,rgnid)==4)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_e) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo-1)
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imin(halo)-halo,imin(halo)-1
                      do j=jmax(halo)+1+(i-(imin(halo)-1)) &
                           ,jmax(halo)+(halo-1)+(i-(imin(halo)-1))
                         n=n+1
                         in=i-j+imax(halo)-imin(halo)+jmax(halo)+2
                         jn=i+imax(halo)-2*imin(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imin(halo)-(halo-1),imin(halo)-1
                      do j=jmax(halo)+1,jmax(halo)+1+(i-(imin(halo)-(halo-1)))
                         n=n+1
                         in=i+imax(halo)-imin(halo)+1
                         jn=j-jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_w) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+1,jmax(halo)+halo
                      do i=imin(halo)-(halo-1)+(j-(jmax(halo)+1)) &
                           ,imin(halo)+(j-(jmax(halo)+1))
                         n=n+1
                         in=j-jmax(halo)+imin(halo)-1
                         jn=-i+j+imin(halo)-jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !!
          !
          !!
          if ((ADM_rgn_vnum(ADM_e,rgnid)==4)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_e,rgnid,2)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_e,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+2,imax(halo)+halo
                      do j=jmax(halo)+1,jmax(halo)+1+(i-(imax(halo)+2))
                         n=n+1
                         in=j-jmax(halo)-1+imin(halo)
                         jn=-i+j+2*imax(halo)-imin(halo)-jmax(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_e,rgnid,2)==ADM_w) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_e,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+1,jmax(halo)+halo
                      do i=imax(halo)+1,imax(halo)+halo
                         n=n+1
                         in=i-imax(halo)+imin(halo)-1
                         jn=j-jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_e,rgnid,2)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_e,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+2,jmax(halo)+halo
                      do i=imax(halo)+1,imax(halo)+1+(j-(jmax(halo)+2))
                         n=n+1
                         in=i-j-imax(halo)+2*jmax(halo)+imin(halo)-jmin(halo)+1
                         jn=i-imax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             if (ADM_rgn_etab(ADM_dir,ADM_se,rgnid)==ADM_sw) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_se,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+halo
                      do j=jmin(halo)+1+(i-(imax(halo)+1)),jmax(halo)+1+(i-(imax(halo)+1))
                         n=n+1
                         in=i-j-imax(halo)+jmax(halo)+imin(halo)
                         jn=i-imax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
                !
             endif
          endif
          !!
          !
          !!
          if (ADM_rgn_vnum(ADM_s,rgnid)==4) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+(halo-1)
                      do j=jmin(halo)-(halo-1)+(i-(imax(halo)+1)),jmin(halo)-1
                         n=n+1
                         in=i-imax(halo)-1+imin(halo)
                         jn=j+jmax(halo)+1-jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_e) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo-1)
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imax(halo)+1+(j-(jmin(halo)-1)) &
                           ,imax(halo)+(halo-1)+(j-(jmin(halo)-1))
                         n=n+1
                         in=j+jmax(halo)-2*jmin(halo)+imin(halo)+1
                         jn=-i+j+imax(halo)+jmax(halo)-jmin(halo)+2
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_w) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+halo
                      do j=jmin(halo)-(halo-1)+(i-(imax(halo)+1)) &
                           ,jmin(halo)+(i-(imax(halo)+1))
                         n=n+1
                         in=i-j-imax(halo)+jmin(halo)+imin(halo)-1
                         jn=i-imax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !
          maxcommrecv_r2r(halo,rgnid)=m
          !
          if ((ADM_rgn_vnum(ADM_w,rgnid)==3)) then
             if ((ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_ne)) then
                n_hemisphere_copy(ADM_w,halo,rgnid)=1
             elseif ((ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_se)) then
                s_hemisphere_copy(ADM_w,halo,rgnid)=1
             endif
          endif
          if ((ADM_rgn_vnum(ADM_n,rgnid)==5)) then
             n_hemisphere_copy(ADM_n,halo,rgnid)=1
          endif
          if ((ADM_rgn_vnum(ADM_s,rgnid)==3)) then
             n_hemisphere_copy(ADM_s,halo,rgnid)=1
          endif
          if ((ADM_rgn_vnum(ADM_e,rgnid)==3)) then
             if ((ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_nw)) then
                n_hemisphere_copy(ADM_e,halo,rgnid)=1
             elseif ((ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_sw)) then
                s_hemisphere_copy(ADM_e,halo,rgnid)=1
             endif
          endif
          if ((ADM_rgn_vnum(ADM_s,rgnid)==5)) then
             s_hemisphere_copy(ADM_s,halo,rgnid)=1
          endif
          if ((ADM_rgn_vnum(ADM_n,rgnid)==3)) then
             s_hemisphere_copy(ADM_n,halo,rgnid)=1
          endif
          !
       enddo !loop l
       !
       !(20101207) removed by teraim
       !do rrgnid=1,ADM_rgn_nmax
       !   do srgnid=1,ADM_rgn_nmax
       !      sendtag_r2r(rrgnid,halo,srgnid)=rrgnid+ADM_rgn_nmax*(srgnid-1)
       !      recvtag_r2r(srgnid,halo,rrgnid)=sendtag_r2r(rrgnid,halo,srgnid)
!      !       write(*,*) 'sendtag_r2r',ADM_prc_me,rrgnid,srgnid,sendtag_r2r(rrgnid,halo,srgnid)
       !   enddo
       !enddo
       !
    enddo !loop halo



    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D(:,:,l) = rsize_r2r(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D,                                     &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        rsize_r2r,                                     &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        ADM_comm_world,                            &
                        ierr                                           )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D(:,:,l) = sourceid_r2r(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D,                                     &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        sourceid_r2r,                                  &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        ADM_comm_world,                            &
                        ierr                                           )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf2D(:,l) = maxcommrecv_r2r(:,rgnid)
    enddo

    call MPI_Allgather( tempbuf2D,                        &
                        halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                      &
                        maxcommrecv_r2r,                  &
                        halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                      &
                        ADM_comm_world,               &
                        ierr                              )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf4D(:,:,:,l) = rlist_r2r(:,:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf4D,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        rlist_r2r,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        ADM_comm_world,                                             &
                        ierr                                                            )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf4D(:,:,:,l) = qlist_r2r(:,:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf4D,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        qlist_r2r,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        ADM_comm_world,                                             &
                        ierr                                                            )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D2(:,:,l) = n_hemisphere_copy(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D2,                         &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        n_hemisphere_copy,                  &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        ADM_comm_world,                 &
                        ierr                                )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D2(:,:,l) = s_hemisphere_copy(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D2,                         &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        s_hemisphere_copy,                  &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        ADM_comm_world,                 &
                        ierr                                )

    call MPI_BARRIER(ADM_comm_world,ierr)

    do halo=1,halomax
       do ls=1,ADM_prc_rnum(ADM_prc_me)
          srgnid=ADM_prc_tab(ls,ADM_prc_me)
          ms=0
          do lr=1,ADM_rgn_nmax
             rrgnid=lr
             do mr=1,maxcommrecv_r2r(halo,rrgnid)
                if (srgnid==sourceid_r2r(mr,halo,rrgnid)) then
                   ms=ms+1
                   !
                   !(20101207)added by teraim
                   if(ADM_rgn2prc(srgnid)==ADM_prc_me) then
                     msend_r2r(rrgnid,halo,srgnid)=ms
                   else
                     write(*,*)"This process is abort because irregular access in msend_r2r."
                     exit
                   endif
                   !
                   destid_r2r(ms,halo,srgnid)=rrgnid
                   ssize_r2r(ms,halo,srgnid)=rsize_r2r(mr,halo,rrgnid)
                   do n=1,rsize_r2r(mr,halo,rrgnid)
                      slist_r2r(n,ms,halo,srgnid)=qlist_r2r(n,mr,halo,rrgnid)
                   enddo
                endif
             enddo
          enddo
          maxcommsend_r2r(halo,srgnid)=ms
       enddo
    enddo !loop halo
    !
    call mpi_barrier(ADM_comm_world,ierr)
    do l=1,ADM_rgn_nmax
       call mpi_bcast(                  &
            destid_r2r(1,1,l),            &
            max_comm_r2r*halomax,  &
            mpi_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_comm_world,             &
            ierr)
    end do
    do l=1,ADM_rgn_nmax
       call mpi_bcast(                  &
            ssize_r2r(1,1,l),            &
            max_comm_r2r*halomax,  &
            mpi_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_comm_world,             &
            ierr)
    end do
    !(20101207)removed by teraim
    !do l=1,ADM_rgn_nmax
    !   call mpi_bcast(                  &
    !        msend_r2r(1,1,l),            &
    !        ADM_rgn_nmax*halomax,  &
    !        mpi_integer,                &
    !        prc_tab_rev(ptr_prcid,l)-1, &
    !        ADM_comm_world,             &
    !        ierr)
    !end do
    do l=1,ADM_rgn_nmax
       call mpi_bcast(                  &
            slist_r2r(1,1,1,l),            &
            max_comm_r2r*max_datasize_r2r*halomax,  &
            mpi_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_comm_world,             &
            ierr)
    end do
    do l=1,ADM_rgn_nmax
       call mpi_bcast(                  &
            maxcommsend_r2r(1,l),            &
            1*halomax,  &
            mpi_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_comm_world,             &
            ierr)
    end do
    !
    allocate(sendbuf_p2r(kmax*max_varmax*2, &
         ADM_rgn_nmax_pl))
    allocate(recvbuf_p2r(kmax*max_varmax*2, &
         ADM_rgn_nmax_pl))
    allocate(recvtag_p2r(max_comm_p2r,ADM_npl:ADM_spl))
    allocate(sendtag_p2r(max_comm_p2r,ADM_npl:ADM_spl))
    do pl=ADM_npl,ADM_spl
       do p=1,ADM_vlink_nmax
          recvtag_p2r(p,pl)=ADM_rgn_nmax*ADM_rgn_nmax+p+ADM_vlink_nmax*(pl-1)
          sendtag_p2r(p,pl)=ADM_rgn_nmax*ADM_rgn_nmax+p+ADM_vlink_nmax*(pl-1)

!          write(*,*) 'sendtag_p2r',ADM_prc_me,p,pl,halo,sendtag_p2r(p,pl)

       enddo
    enddo
    !
    allocate(clist(max_varmax))
    !
    call mpi_barrier(ADM_comm_world,ierr)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  re-setup comm_table !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rank_me=ADM_prc_me-1
    max_comm_prc=min(ADM_prc_all,max_comm_r2r*ADM_lall+2*max_comm_r2p)
    !
    allocate(n_nspl(ADM_npl:ADM_spl,halomax))
    do halo=1,halomax
       n_nspl(ADM_npl,halo)=suf(imin(halo)+0,jmax(halo)+1,gall(halo))
       n_nspl(ADM_spl,halo)=suf(imax(halo)+1,jmin(halo)+0,gall(halo))
    enddo
    !
    allocate(temp_sendorder(0:ADM_prc_all-1,halomax))
    allocate(temp_recvorder(0:ADM_prc_all-1,halomax))
    !
    !--------------------------------------------------
    allocate(romax(halomax))
    allocate(somax(halomax))
    allocate(sourcerank(max_comm_prc,halomax))
    allocate(destrank(max_comm_prc,halomax))
    allocate(rsize(max_comm_prc,halomax))
    allocate(ssize(max_comm_prc,halomax))
    romax(:)=0
    somax(:)=0
    sourcerank(:,:)=-1
    destrank(:,:)=-1
    rsize(:,:)=0
    ssize(:,:)=0
    !--------------------------------------------------
    !
    maxn=((gmax(halomax)-gmin(halomax)+1)+2)*halomax
    maxm=max_comm_r2r+1
    maxl=ADM_lall+2
    !----
    maxn_pl=halomax*(halomax+1)/2
    maxm_pl=ADM_vlink_nmax
    maxl_pl=(ADM_spl-ADM_npl+1)
    !----
    maxn_r2r=(gmax(halomax)-gmin(halomax)+1)*halomax
    maxm_r2r=max_comm_r2r
    maxl_r2r=ADM_lall
    !----
    maxn_r2p=halomax*(halomax+1)/2
    maxm_r2p=ADM_vlink_nmax
    maxl_r2p=(ADM_spl-ADM_npl+1)
    !----
    maxn_p2r=1
    maxm_p2r=ADM_vlink_nmax
    maxl_p2r=(ADM_spl-ADM_npl+1)
    !----
    maxn_sgp=halomax
    maxm_sgp=4
    maxl_sgp=12
    !
    !--------------------------------------------------
    !  for send
    !--------------------------------------------------
    allocate(nsmax(max_comm_prc,halomax))
    allocate(sendinfo(elemsize_comm,maxm*maxl,max_comm_prc,halomax))
    allocate(sendlist(maxn,maxm*maxl,max_comm_prc,halomax))
    nsmax(:,:)=0
    sendinfo(:,:,:,:)=0
    sendlist(:,:,:,:)=0
    allocate(nsmax_pl(max_comm_prc,halomax))
    allocate(sendinfo_pl(elemsize_comm,maxm_pl*maxl_pl,max_comm_prc,halomax))
    allocate(sendlist_pl(maxn_pl,maxm_pl*maxl_pl,max_comm_prc,halomax))
    nsmax_pl(:,:)=0
    sendinfo_pl(:,:,:,:)=0
    sendlist_pl(:,:,:,:)=0
    !--------------------------------------------------
    !
    !--------------------------------------------------
    !  for copy
    !--------------------------------------------------
    allocate(ncmax_r2r(halomax))
    allocate(copyinfo_r2r(elemsize_copy,maxm_r2r*maxl_r2r,halomax))
    allocate(recvlist_r2r(maxn_r2r,maxm_r2r*maxl_r2r,halomax))
    allocate(sendlist_r2r(maxn_r2r,maxm_r2r*maxl_r2r,halomax))
    ncmax_r2r(:)=0
    copyinfo_r2r(:,:,:)=0
    recvlist_r2r(:,:,:)=0
    sendlist_r2r(:,:,:)=0
    !--------------------------------------------------
    allocate(ncmax_r2p(halomax))
    allocate(copyinfo_r2p(elemsize_copy,maxm_r2p*maxl_r2p,halomax))
    allocate(recvlist_r2p(maxn_r2p,maxm_r2p*maxl_r2p,halomax))
    allocate(sendlist_r2p(maxn_r2p,maxm_r2p*maxl_r2p,halomax))
    ncmax_r2p(:)=0
    copyinfo_r2p(:,:,:)=0
    recvlist_r2p(:,:,:)=0
    sendlist_r2p(:,:,:)=0
    !--------------------------------------------------
    allocate(ncmax_p2r(halomax))
    allocate(copyinfo_p2r(elemsize_copy,maxm_p2r*maxl_p2r,halomax))
    allocate(recvlist_p2r(maxn_p2r,maxm_p2r*maxl_p2r,halomax))
    allocate(sendlist_p2r(maxn_p2r,maxm_p2r*maxl_p2r,halomax))
    ncmax_p2r(:)=0
    copyinfo_p2r(:,:,:)=0
    recvlist_p2r(:,:,:)=0
    sendlist_p2r(:,:,:)=0
    !--------------------------------------------------
    !
    !--------------------------------------------------
    !  for recv
    !--------------------------------------------------
    allocate(nrmax(max_comm_prc,halomax))
    allocate(recvinfo(elemsize_comm,maxm*maxl,max_comm_prc,halomax))
    allocate(recvlist(maxn,maxm*maxl,max_comm_prc,halomax))
    nrmax(:,:)=0
    recvinfo(:,:,:,:)=0
    recvlist(:,:,:,:)=0
    allocate(nrmax_pl(max_comm_prc,halomax))
    allocate(recvinfo_pl(elemsize_comm,maxm_pl*maxl_pl,max_comm_prc,halomax))
    allocate(recvlist_pl(maxn_pl,maxm_pl*maxl_pl,max_comm_prc,halomax))
    nrmax_pl(:,:)=0
    recvinfo_pl(:,:,:,:)=0
    recvlist_pl(:,:,:,:)=0
    !--------------------------------------------------
    allocate(temp_dest_rgn(maxm*maxl,max_comm_prc,halomax))
    allocate(temp_src_rgn(maxm*maxl,max_comm_prc,halomax))
    allocate(temp_dest_rgn_pl(maxm_pl*maxl_pl,max_comm_prc,halomax))
    allocate(temp_src_rgn_pl(maxm_pl*maxl_pl,max_comm_prc,halomax))
    temp_dest_rgn(:,:,:)=0
    temp_dest_rgn_pl(:,:,:)=0
    temp_src_rgn(:,:,:)=0
    temp_src_rgn_pl(:,:,:)=0
    !--------------------------------------------------
    !
    do halo=1,halomax
       do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          do m=1,maxcommrecv_r2r(halo,rgnid)
             srank=prc_tab_rev(ptr_prcid,sourceid_r2r(m,halo,rgnid))-1
             if (srank/=rank_me) then
                ck=0
                loop_ro1:do ro=1,romax(halo)
                   if (srank==sourcerank(ro,halo)) exit loop_ro1
                   ck=ck+1
                enddo loop_ro1
                if (ck==romax(halo)) then
                   romax(halo)=romax(halo)+1
                   ro=romax(halo)
                   sourcerank(ro,halo)=srank
                   temp_recvorder(srank,halo)=ro
                endif
                ro=temp_recvorder(srank,halo)
                nrmax(ro,halo)=nrmax(ro,halo)+1
                nr=nrmax(ro,halo)
                recvinfo(SIZE_COMM,nr,ro,halo)=rsize_r2r(m,halo,rgnid)
                recvinfo(LRGNID_COMM,nr,ro,halo)=l
                temp_src_rgn(nr,ro,halo)=sourceid_r2r(m,halo,rgnid)
                rs=recvinfo(SIZE_COMM,nr,ro,halo)
                rsize(ro,halo)=rsize(ro,halo)+rs
                do n=1,rs
                   recvlist(n,nr,ro,halo)=rlist_r2r(n,m,halo,rgnid)
                enddo
             else
                ncmax_r2r(halo)=ncmax_r2r(halo)+1
                nc=ncmax_r2r(halo)
                copyinfo_r2r(SIZE_COPY,nc,halo)=rsize_r2r(m,halo,rgnid)
                copyinfo_r2r(LRGNID_COPY,nc,halo)=l
                copyinfo_r2r(SRC_LRGNID_COPY,nc,halo) &
                     =prc_tab_rev(ptr_lrgnid,sourceid_r2r(m,halo,rgnid))
                cs=copyinfo_r2r(SIZE_COPY,nc,halo)
                srgnid=sourceid_r2r(m,halo,rgnid)
                do n=1,cs
                   recvlist_r2r(n,nc,halo)=rlist_r2r(n,m,halo,rgnid)
                   !
                   !(20101207)added by teraim
                   if(ADM_rgn2prc(srgnid)==ADM_prc_me) then
                     sendlist_r2r(n,nc,halo)=slist_r2r(n,msend_r2r(rgnid,halo,srgnid),halo,srgnid)
                   else
                     write(*,*)"This process is abort because irregular access is msend_r2r."
                     exit
                   endif
                   !
                enddo
             endif
          enddo !loop m
          !enddo !loop l
          !!
          !do l=1,ADM_lall
          !  rgnid=ADM_prc_tab(l,ADM_prc_me)
          do m=1,maxcommsend_r2r(halo,rgnid)
             drank=prc_tab_rev(ptr_prcid,destid_r2r(m,halo,rgnid))-1
             if (drank/=rank_me) then
                ck=0
                loop_so1:do so=1,somax(halo)
                   if (drank==destrank(so,halo)) exit loop_so1
                   ck=ck+1
                enddo loop_so1
                if (ck==somax(halo)) then
                   somax(halo)=somax(halo)+1
                   so=somax(halo)
                   destrank(so,halo)=drank
                   temp_sendorder(drank,halo)=so
                endif
                so=temp_sendorder(drank,halo)
                nsmax(so,halo)=nsmax(so,halo)+1
                ns=nsmax(so,halo)
                sendinfo(SIZE_COMM,ns,so,halo)=ssize_r2r(m,halo,rgnid)
                sendinfo(LRGNID_COMM,ns,so,halo)=l
                temp_dest_rgn(ns,so,halo)=destid_r2r(m,halo,rgnid)
                ss=sendinfo(SIZE_COMM,ns,so,halo)
                ssize(so,halo)=ssize(so,halo)+ss
                do n=1,ss
                   sendlist(n,ns,so,halo)=slist_r2r(n,m,halo,rgnid)
                enddo
             endif
          enddo !loop m
       enddo !loop l
       !enddo !loop halo
       !do halo=1,halomax
       if(comm_pl) call re_setup_pl_comm_info ! T.Ohno 110721
    enddo !loop halo
    deallocate(temp_sendorder)
    deallocate(temp_recvorder)
    !
    !allocate(temp_sb(ADM_rgn_nmax+2,halomax,ADM_rgn_nmax+2)) !(20101207) removed by teraim
    allocate(tsb(somax(halomax)))
    !temp_sb(:,:,:)=0 !(20101207) removed by teraim

    call init_tempsb !(20101207) added by teraim

    do halo=1,halomax
       tsb(:)=0
       do so=1,somax(halo)
          do ns=1,nsmax(so,halo)
             ss=sendinfo(SIZE_COMM,ns,so,halo)
             srgnid=ADM_prc_tab(sendinfo(LRGNID_COMM,ns,so,halo),ADM_prc_me)
             rrgnid=temp_dest_rgn(ns,so,halo)
             sendinfo(BASE_COMM,ns,so,halo)=tsb(so)
             !temp_sb(rrgnid,halo,srgnid)=tsb(so) !(20101207)removed by teraim
             call add_tempsb(rrgnid, srgnid, tsb(so)) !(20101207)added by teraim
             tsb(so)=tsb(so)+ss
          enddo
          do ns=1,nsmax_pl(so,halo)
             ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
             pl=sendinfo_pl(LRGNID_COMM,ns,so,halo)
             srgnid=ADM_rgn_nmax+pl
             rrgnid=temp_dest_rgn_pl(ns,so,halo)
             sendinfo_pl(BASE_COMM,ns,so,halo)=tsb(so)
             !temp_sb(rrgnid,halo,srgnid)=tsb(so) !(20101207)removed by teraim
             call add_tempsb(rrgnid, srgnid, tsb(so)) !(20101207)added by teraim
             tsb(so)=tsb(so)+ss
          enddo
       enddo
    enddo
    deallocate(tsb)
    !
    !(20101207)removed by teraim
    !call mpi_barrier(ADM_comm_world,ierr)
    !do l=1,ADM_rgn_nmax
    !   call mpi_bcast(                  &
    !        temp_sb(1,1,l),        &
    !        (ADM_rgn_nmax+2)*halomax,      &
    !        mpi_integer,                &
    !        prc_tab_rev(ptr_prcid,l)-1, &
    !        ADM_comm_world,             &
    !        ierr)
    !enddo
    !do pl=ADM_npl,ADM_spl
    !   call mpi_bcast(                  &
    !        temp_sb(1,1,ADM_rgn_nmax+pl),       &
    !        (ADM_rgn_nmax+2)*halomax,      &
    !        mpi_integer,                &
    !        ADM_prc_nspl(pl)-1,         &
    !        ADM_comm_world,             &
    !        ierr)
    !enddo
    !call mpi_barrier(ADM_comm_world,ierr)
    !
    !(20101207)added by teraim
    call mpi_barrier(ADM_comm_world,ierr)
    do l=1,ADM_rgn_nmax
       call mpi_bcast(                  &
            tempsb(l)%col,              &
            max_size,                   &
            mpi_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_comm_world,         &
            ierr)
       call mpi_bcast(                  &
            tempsb(l)%val,              &
            max_size,                   &
            mpi_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_comm_world,         &
            ierr)
    enddo
    if(comm_pl) then ! T.Ohno 110721
      do pl=ADM_npl,ADM_spl
         call mpi_bcast(                    &
              tempsb(ADM_rgn_nmax+pl)%col,  &
              max_size,                     &
              mpi_integer,                  &
              ADM_prc_nspl(pl)-1,           &
              ADM_comm_world,           &
              ierr)
         call mpi_bcast(                    &
              tempsb(ADM_rgn_nmax+pl)%val,  &
              max_size,                     &
              mpi_integer,                  &
              ADM_prc_nspl(pl)-1,           &
              ADM_comm_world,           &
              ierr)
      enddo
    endif ! T.Ohno 110721
    call mpi_barrier(ADM_comm_world,ierr)
    !
    !call show_tempsb !(20101209) added by teraim
    !
    do halo=1,halomax
       do ro=1,romax(halo)
          do nr=1,nrmax(ro,halo)
             rrgnid=ADM_prc_tab(recvinfo(LRGNID_COMM,nr,ro,halo),ADM_prc_me)
             srgnid=temp_src_rgn(nr,ro,halo)
             !recvinfo(BASE_COMM,nr,ro,halo)=temp_sb(rrgnid,halo,srgnid) !(20101207)removed by teraim
             !(20101207) added by teraim
             call get_tempsb(rrgnid,srgnid,ret)
             recvinfo(BASE_COMM,nr,ro,halo)=ret
          enddo
          do nr=1,nrmax_pl(ro,halo)
             pl=recvinfo_pl(LRGNID_COMM,nr,ro,halo)
             rrgnid=pl+ADM_rgn_nmax
             srgnid=temp_src_rgn_pl(nr,ro,halo)
             !recvinfo_pl(BASE_COMM,nr,ro,halo)=temp_sb(rrgnid,halo,srgnid) !(20101207)removed by teraim
             !(20101207) added by teraim
             call get_tempsb(rrgnid,srgnid,ret)
             recvinfo_pl(BASE_COMM,nr,ro,halo)=ret
          enddo
       enddo !loop ro
    enddo !loop halo
    deallocate(temp_dest_rgn)
    deallocate(temp_dest_rgn_pl)
    deallocate(temp_src_rgn)
    deallocate(temp_src_rgn_pl)
    !deallocate(temp_sb) !(20101207)removed by teraim
    call finalize_tempsb !(20101207)added by teraim
    !
    allocate(recvtag(romax(halomax),halomax))
    allocate(sendtag(somax(halomax),halomax))
    recvtag(:,:)=-1
    sendtag(:,:)=-1
    do halo=1,halomax
       do ro=1,romax(halo)
          recvtag(ro,halo)=rank_me
       enddo
       do so=1,somax(halo)
          sendtag(so,halo)=destrank(so,halo)
       enddo
    enddo
!    maxdatasize=(max_comm_r2r*(gmax(halomax)-gmin(halomax)+1)+2*max_comm_r2p*(halomax+1)/2)*halomax*kmax*max_varmax
!    maxdatasize=(maxn_r2r*maxm_r2r*maxl_r2r+maxn_r2p*maxm_r2p*maxl_r2p+maxn_p2r*maxm_p2r*maxl_p2r)*kmax*max_varmax
    maxdatasize_s=0
    do so=1,somax(halomax)
      maxdatasize_s=maxdatasize_s+ssize(so,halomax)*kmax*max_varmax
    enddo
    maxdatasize_r=0
    do ro=1,romax(halomax)
      maxdatasize_r=maxdatasize_r+rsize(ro,halomax)*kmax*max_varmax
    enddo
    allocate(recvbuf(maxdatasize_r,romax(halomax)))
    allocate(sendbuf(maxdatasize_s,somax(halomax)))
    recvbuf(:,:)=0.D0
    sendbuf(:,:)=0.D0

!!    allocate(comm_dbg_recvbuf(maxdatasize_r,romax(halomax),2)) !iga
!!    allocate(comm_dbg_sendbuf(maxdatasize_s,somax(halomax),2)) !iga
!!    comm_dbg_recvbuf=CNST_UNDEF !iga
!!    comm_dbg_sendbuf=CNST_UNDEF !iga

    !
    allocate(ncmax_sgp(halomax))
    allocate(copyinfo_sgp(elemsize_copy,maxm_sgp*maxl_sgp,halomax))
    allocate(recvlist_sgp(maxn_sgp,maxm_sgp*maxl_sgp,halomax))
    allocate(sendlist_sgp(maxn_sgp,maxm_sgp*maxl_sgp,halomax))
    ncmax_sgp(:)=0
    copyinfo_sgp(:,:,:)=0
    recvlist_sgp(:,:,:)=0
    sendlist_sgp(:,:,:)=0
    do halo=1,halomax
       ncmax_sgp(halo)=0
       do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          if (n_hemisphere_copy(ADM_w,halo,rgnid)==1) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmin(halo)-n,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo),gmin(halo)-n,gall(halo))
             enddo
          endif
          if ((n_hemisphere_copy(ADM_n,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo),gmax(halo)+n+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmax(halo)+1,gall(halo))
             enddo
          endif
          if ((n_hemisphere_copy(ADM_e,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+n+1,gall(halo))
             enddo
          endif
          if ((n_hemisphere_copy(ADM_s,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+1,gmin(halo)-n,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmin(halo),gall(halo))
             enddo
          endif
          if (s_hemisphere_copy(ADM_w,halo,rgnid)==1) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo),gmin(halo)-n,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmin(halo)-n,gall(halo))
             enddo
          endif
          if ((s_hemisphere_copy(ADM_n,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmax(halo)+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo),gmax(halo)+n+1,gall(halo))
             enddo
          endif
          if ((s_hemisphere_copy(ADM_e,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+n+1,gall(halo))
             enddo
          endif
          if ((s_hemisphere_copy(ADM_s,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmin(halo),gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+1,gmin(halo)-n,gall(halo))
             enddo
          endif
       enddo !loop l
    enddo !loop halo
    !
    !-- for output_info  ---
    allocate( src_rank_all(max_comm_prc,halomax,ADM_prc_all))
    allocate(dest_rank_all(max_comm_prc,halomax,ADM_prc_all))
    src_rank_all(:,:,:)=-1
    dest_rank_all(:,:,:)=-1
    src_rank_all(:,:,ADM_prc_me)=sourcerank(:,:)
    dest_rank_all(:,:,ADM_prc_me)=destrank(:,:)
    call mpi_barrier(ADM_comm_world,ierr)
    do l=1,ADM_prc_all
       call mpi_bcast(                  &
            src_rank_all(1,1,l),        &
            max_comm_prc*halomax,       &
            mpi_integer,                &
            l-1,                        &
            ADM_comm_world,             &
            ierr)
       call mpi_bcast(                  &
            dest_rank_all(1,1,l),       &
            max_comm_prc*halomax,       &
            mpi_integer,                &
            l-1,                        &
            ADM_comm_world,             &
            ierr)
    enddo
    !
    call MPI_Barrier(ADM_COMM_WORLD,ierr)
    !
    !--- output for debug
    if(present(debug)) then
       if(debug) call output_info
    end if
    !
    ! <== iga for dbg 090917
    if (opt_comm_dbg) then
!       dbg_sendbuf_init = -1d66  * (ADM_prc_me+1000)
       dbg_sendbuf_init = -888d66
       dbg_recvbuf_init = -777d66
       allocate(dbg_areq_save(2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4),4))
       dbg_areq_save(:,:) = -999 ! [Add] 12/03/26 T.Seiki
    endif
    ! iga for dbg 090916 ==>
    contains
    !
    subroutine re_setup_pl_comm_info
       do pl=ADM_npl,ADM_spl
          if (ADM_prc_me==ADM_prc_nspl(pl)) then
             do p=1,maxcommrecv_r2p(pl,halo)
                srank=source_prc_r2p(p,pl,halo)-1
                if (srank/=rank_me) then
                   ck=0
                   loop_ro2:do ro=1,romax(halo)
                      if (srank==sourcerank(ro,halo)) exit loop_ro2
                      ck=ck+1
                   enddo loop_ro2
                   if (ck==romax(halo)) then
                      romax(halo)=romax(halo)+1
                      ro=romax(halo)
                      sourcerank(ro,halo)=srank
                      temp_recvorder(srank,halo)=ro
                   endif
                   ro=temp_recvorder(srank,halo)
                   nrmax_pl(ro,halo)=nrmax_pl(ro,halo)+1
                   nr=nrmax_pl(ro,halo)
                   recvinfo_pl(SIZE_COMM,nr,ro,halo)=rsize_r2p(p,pl,halo)
                   recvinfo_pl(LRGNID_COMM,nr,ro,halo)=pl
                   temp_src_rgn_pl(nr,ro,halo)=ADM_prc_tab(source_rgn_r2p(p,pl,halo),srank+1)
                   rs=recvinfo_pl(SIZE_COMM,nr,ro,halo)
                   rsize(ro,halo)=rsize(ro,halo)+rs
                   do n=1,rs
                      recvlist_pl(n,nr,ro,halo)=rlist_r2p(n,p,pl,halo)
                   enddo
                else
                   ncmax_r2p(halo)=ncmax_r2p(halo)+1
                   nc=ncmax_r2p(halo)
                   copyinfo_r2p(SIZE_COPY,nc,halo)=rsize_r2p(p,pl,halo)
                   copyinfo_r2p(LRGNID_COPY,nc,halo)=pl
                   copyinfo_r2p(SRC_LRGNID_COPY,nc,halo)=source_rgn_r2p(p,pl,halo)
                   cs=copyinfo_r2p(SIZE_COPY,nc,halo)
                   do n=1,cs
                      recvlist_r2p(n,nc,halo)=rlist_r2p(n,p,pl,halo)
                      sendlist_r2p(n,nc,halo)=slist_r2p(n,p,pl,halo)
                   enddo
                endif
             enddo !loop p
             !
             do p=1,ADM_vlink_nmax
                rgnid=ADM_rgn_vtab_pl(ADM_rid,pl,p)
                drank=prc_tab_rev(ptr_prcid,rgnid)-1
                if (drank/=rank_me) then
                   ck=0
                   loop_so2:do so=1,somax(halo)
                      if (drank==destrank(so,halo)) exit loop_so2
                      ck=ck+1
                   enddo loop_so2
                   if (ck==somax(halo)) then
                      somax(halo)=somax(halo)+1
                      so=somax(halo)
                      destrank(so,halo)=drank
                      temp_sendorder(drank,halo)=so
                   endif
                   so=temp_sendorder(drank,halo)
                   nsmax_pl(so,halo)=nsmax_pl(so,halo)+1
                   ns=nsmax_pl(so,halo)
                   sendinfo_pl(SIZE_COMM,ns,so,halo)=1
                   sendinfo_pl(LRGNID_COMM,ns,so,halo)=pl
                   temp_dest_rgn_pl(ns,so,halo)=rgnid
                   ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
                   ssize(so,halo)=ssize(so,halo)+ss
                   do n=1,ss
                      sendlist_pl(n,ns,so,halo)=ADM_gslf_pl
                   enddo
                endif
             enddo !loop p
          endif
          !
          do p=1,ADM_vlink_nmax
             rgnid=ADM_rgn_vtab_pl(ADM_rid,pl,p)
             drank=prc_tab_rev(ptr_prcid,rgnid)-1
             if (rank_me==drank) then
                srank=ADM_prc_nspl(pl)-1
                if (srank/=rank_me) then
                   ck=0
                   loop_ro3:do ro=1,romax(halo)
                      if (srank==sourcerank(ro,halo)) exit loop_ro3
                      ck=ck+1
                   enddo loop_ro3
                   if (ck==romax(halo)) then
                      romax(halo)=romax(halo)+1
                      ro=romax(halo)
                      sourcerank(ro,halo)=srank
                      temp_recvorder(srank,halo)=ro
                   endif
                   ro=temp_recvorder(srank,halo)
                   nrmax(ro,halo)=nrmax(ro,halo)+1
                   nr=nrmax(ro,halo)
                   recvinfo(SIZE_COMM,nr,ro,halo)=1
                   recvinfo(LRGNID_COMM,nr,ro,halo)=prc_tab_rev(ptr_lrgnid,rgnid)
                   temp_src_rgn(nr,ro,halo)=ADM_rgn_nmax+pl
                   rs=recvinfo(SIZE_COMM,nr,ro,halo)
                   rsize(ro,halo)=rsize(ro,halo)+rs
                   do n=1,rs
                      recvlist(n,nr,ro,halo)=n_nspl(pl,halo)
                   enddo
                else
                   ncmax_p2r(halo)=ncmax_p2r(halo)+1
                   nc=ncmax_p2r(halo)
                   copyinfo_p2r(SIZE_COPY,nc,halo)=1
                   copyinfo_p2r(LRGNID_COPY,nc,halo)=prc_tab_rev(ptr_lrgnid,rgnid)
                   copyinfo_p2r(SRC_LRGNID_COPY,nc,halo)=pl
                   cs=copyinfo_p2r(SIZE_COPY,nc,halo)
                   do n=1,cs
                      recvlist_p2r(n,nc,halo)=n_nspl(pl,halo)
                      sendlist_p2r(n,nc,halo)=ADM_gslf_pl
                   enddo
                endif
             endif
          enddo !loop p
          !
          do p=1,maxcommsend_r2p(pl,halo)
             srank=source_prc_r2p(p,pl,halo)-1
             if (rank_me==srank) then
                rgnid=ADM_rgn_vtab_pl(ADM_rid,pl,p)
                drank=ADM_prc_nspl(pl)-1
                if (drank/=rank_me) then
                   ck=0
                   loop_so3:do so=1,somax(halo)
                      if (drank==destrank(so,halo)) exit loop_so3
                      ck=ck+1
                   enddo loop_so3
                   if (ck==somax(halo)) then
                      somax(halo)=somax(halo)+1
                      so=somax(halo)
                      destrank(so,halo)=drank
                      temp_sendorder(drank,halo)=so
                   endif
                   so=temp_sendorder(drank,halo)
                   nsmax(so,halo)=nsmax(so,halo)+1
                   ns=nsmax(so,halo)
                   sendinfo(SIZE_COMM,ns,so,halo)=ssize_r2p(p,pl,halo)
                   sendinfo(LRGNID_COMM,ns,so,halo)=prc_tab_rev(ptr_lrgnid,rgnid)
                   temp_dest_rgn(ns,so,halo)=ADM_rgn_nmax+pl
                   ss=sendinfo(SIZE_COMM,ns,so,halo)
                   ssize(so,halo)=ssize(so,halo)+ss
                   do n=1,ss
                      sendlist(n,ns,so,halo)=slist_r2p(n,p,pl,halo)
                   enddo
                endif
             endif
          enddo !loop p
       enddo !loop pl
    end subroutine re_setup_pl_comm_info

  end subroutine COMM_setup
  !-----------------------------------------------------------------------------
  subroutine output_info
    use mod_adm, only :    &
         !--- public parameters
         ADM_log_fid,    &
         !--- public variables
         ADM_prc_all
    !
    implicit none
    !
    integer::halo
    !
    integer :: varmax
    integer :: cmax
    !
    !integer :: srgnid,rrgnid
    !
    integer :: ns,nr
    integer :: so,sl,sb,ss
    integer :: ro,rl,rb,rs
    integer :: l,n
    !
    write(ADM_log_fid,*)
    write(ADM_log_fid,*) &
         'msg : sub[output_info]/mod[comm]'
    write(ADM_log_fid,*) &
         'version : comm.f90.test5.2.1_wtime'
    write(ADM_log_fid,*) &
         '---------------------------------------&
         &       commnication table  start       &
         &---------------------------------------'
    !
    varmax=1
    cmax=kmax*varmax
    do halo=1,halomax
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &       halo region =',halo,'           &
            &---------------------------------------'
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                count                  &
            &---------------------------------------'
       write(ADM_log_fid,*) &
            'romax =',romax(halo) &
            ,'somax =',somax(halo)
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                send                   &
            &---------------------------------------'
       do so=1,somax(halo)
          write(ADM_log_fid,*) &
               'so =',so   &
               ,'mrank =',rank_me   &
               ,'drank =',destrank(so,halo)
       enddo
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                recv                   &
            &---------------------------------------'
       do ro=1,romax(halo)
          write(ADM_log_fid,*) &
               'ro =',ro   &
               ,'mrank =',rank_me   &
               ,'srank =',sourcerank(ro,halo)
       enddo
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                table                   &
            &---------------------------------------'
       do l=1,ADM_prc_all
          do n=1,max_comm_prc
             if (dest_rank_all(n,halo,l)==-1) cycle
             write(ADM_log_fid,*) &
                  'n =',n   &
                  ,'rank =',l-1   &
                  ,'dest_rank =',dest_rank_all(n,halo,l) &
                  ,' src_rank =', src_rank_all(n,halo,l)
          enddo
       enddo
       do so=1,somax(halo)
          do ns=1,nsmax(so,halo)
             ss=sendinfo(SIZE_COMM,ns,so,halo)
             sl=sendinfo(LRGNID_COMM,ns,so,halo)
             sb=sendinfo(BASE_COMM,ns,so,halo)*cmax
             write(ADM_log_fid,*) &
                  'so =',so   &
                  ,'rank =',rank_me   &
                  ,' dest_rank =', destrank(so,halo) &
                  ,' ss =', ss &
                  ,' sl =', sl &
                  ,' sb =', sb
          enddo
          do ns=1,nsmax_pl(so,halo)
             ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
             sl=sendinfo_pl(LRGNID_COMM,ns,so,halo)
             sb=sendinfo_pl(BASE_COMM,ns,so,halo)*cmax
             write(ADM_log_fid,*) &
                  'so =',so   &
                  ,'rank =',rank_me   &
                  ,' dest_rank =', destrank(so,halo) &
                  ,' ss =', ss &
                  ,' sl =', sl &
                  ,' sb =', sb
          enddo
       enddo !loop so
       do ro=1,romax(halo)
          do nr=1,nrmax(ro,halo)
             rs=recvinfo(SIZE_COMM,nr,ro,halo)
             rl=recvinfo(LRGNID_COMM,nr,ro,halo)
             rb=recvinfo(BASE_COMM,nr,ro,halo)*cmax
             write(ADM_log_fid,*) &
                  'ro =',ro   &
                  ,'rank =',rank_me   &
                  ,' src_rank =', sourcerank(ro,halo) &
                  ,' rs =', rs &
                  ,' rl =', rl &
                  ,' rb =', rb
          enddo
          do nr=1,nrmax_pl(ro,halo)
             rs=recvinfo_pl(SIZE_COMM,nr,ro,halo)
             rl=recvinfo_pl(LRGNID_COMM,nr,ro,halo)
             rb=recvinfo_pl(BASE_COMM,nr,ro,halo)*cmax
             write(ADM_log_fid,*) &
                  'ro =',ro   &
                  ,'rank =',rank_me   &
                  ,' src_rank =', sourcerank(ro,halo) &
                  ,' rs =', rs &
                  ,' rl =', rl &
                  ,' rb =', rb
          enddo
       enddo !loop ro
    enddo !loop halo
    !
    write(ADM_log_fid,*) &
         '---------------------------------------&
         &       commnication table  end         &
         &---------------------------------------'
    !
    !call ADM_proc_stop
    return
    !
  end subroutine output_info
  !-----------------------------------------------------------------------------
  subroutine output_time
    use mod_adm, only :    &
         !--- public parameters
         ADM_log_fid,    &
         !--- public variables
         ADM_prc_all,    &
         ADM_comm_world
    !
    implicit none
    !
    real(8)::           &
         ave_time_total,   &
         ave_time_pre ,    &
         ave_time_bar1,    &
         ave_time_bar2,    &
         ave_time_sbuf,    &
         ave_time_recv,    &
         ave_time_send,    &
         ave_time_copy,    &
         ave_time_wait,    &
         ave_time_rbuf,    &
         ave_time_copy_sgp
    real(8)::           &
         min_time_total,   &
         min_time_pre ,    &
         min_time_bar1,    &
         min_time_bar2,    &
         min_time_sbuf,    &
         min_time_recv,    &
         min_time_send,    &
         min_time_copy,    &
         min_time_wait,    &
         min_time_rbuf,    &
         min_time_copy_sgp
    real(8)::           &
         max_time_total,   &
         max_time_pre ,    &
         max_time_bar1,    &
         max_time_bar2,    &
         max_time_sbuf,    &
         max_time_recv,    &
         max_time_send,    &
         max_time_copy,    &
         max_time_wait,    &
         max_time_rbuf,    &
         max_time_copy_sgp
    real(8)::           &
         ave_size_total,   &
         min_size_total,   &
         max_size_total
    real(8)::           &
         ave_comm_count,   &
         min_comm_count,   &
         max_comm_count
    integer::ierr
    !
    write(ADM_log_fid,*)
    write(ADM_log_fid,*) &
         'msg : sub[output_time]/mod[comm]'
    write(ADM_log_fid,*) &
         'version : comm.f90.test5.2.1_wtime'
    write(ADM_log_fid,*) &
         '---------------------------------------&
         &       commnication time               &
         &---------------------------------------'
    !
    call mpi_reduce(time_total,ave_time_total,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_pre ,ave_time_pre ,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_bar1,ave_time_bar1,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_bar2,ave_time_bar2,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_sbuf,ave_time_sbuf,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_recv ,ave_time_recv ,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_send ,ave_time_send ,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_copy,ave_time_copy,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_wait,ave_time_wait,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_rbuf,ave_time_rbuf,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_copy_sgp,ave_time_copy_sgp,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(time_total ,min_time_total ,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_pre ,min_time_pre ,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_bar1,min_time_bar1,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_bar2,min_time_bar2,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_sbuf,min_time_sbuf,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_recv,min_time_recv ,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_send,min_time_send ,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_copy,min_time_copy,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_wait,min_time_wait,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_rbuf,min_time_rbuf,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_copy_sgp,min_time_copy_sgp,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(time_total ,max_time_total ,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_pre ,max_time_pre ,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_bar1,max_time_bar1,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_bar2,max_time_bar2,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_sbuf,max_time_sbuf,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_recv,max_time_recv ,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_send,max_time_send ,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_copy,max_time_copy,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_wait,max_time_wait,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_rbuf,max_time_rbuf,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(time_copy_sgp,max_time_copy_sgp,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    !-----
    call mpi_reduce(size_total,ave_size_total,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(size_total,min_size_total,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(size_total,max_size_total,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    call mpi_reduce(comm_count,ave_comm_count,1,mpi_double_precision,mpi_sum,0,ADM_comm_world,ierr)
    call mpi_reduce(comm_count,min_comm_count,1,mpi_double_precision,mpi_min,0,ADM_comm_world,ierr)
    call mpi_reduce(comm_count,max_comm_count,1,mpi_double_precision,mpi_max,0,ADM_comm_world,ierr)
    write(ADM_log_fid,*) "------- Using unit in the following tables is msec & KB -------"
    !-----
    write(ADM_log_fid,*) "---------  average by comm_call_count[",comm_call_count,"] --------"
    write(ADM_log_fid,'(a,1e11.4)') "time_pre      =",1000*time_pre /comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_bar1     =",1000*time_bar1/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_recv     =",1000*time_recv/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_sbuf     =",1000*time_sbuf/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_send     =",1000*time_send/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_copy     =",1000*time_copy/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_wait     =",1000*time_wait/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_rbuf     =",1000*time_rbuf/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_copy_sgp =",1000*time_copy_sgp/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_bar2     =",1000*time_bar2/comm_call_count
    write(ADM_log_fid,'(a,1e11.4)') "time_total    =",1000*time_total/comm_call_count
    write(ADM_log_fid,'(a,1f11.4)') "size_total    =",size_total*8/1024/comm_call_count
    write(ADM_log_fid,'(a,1f11.4)') "comm_count    =",comm_count/comm_call_count
    !-----
    if (rank_me==0) then
       write(ADM_log_fid,*) "---------  average  by ADM_prc_all & comm_call_count --------"
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_pre      =",1000*ave_time_pre /ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_bar1     =",1000*ave_time_bar1/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_recv     =",1000*ave_time_recv/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_sbuf     =",1000*ave_time_sbuf/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_send     =",1000*ave_time_send/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_copy     =",1000*ave_time_copy/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_wait     =",1000*ave_time_wait/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_rbuf     =",1000*ave_time_rbuf/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_copy_sgp =",1000*ave_time_copy_sgp/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_bar2     =",1000*ave_time_bar2/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "ave_time_total    =",1000*ave_time_total/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1f11.4)') "ave_size_total    =",ave_size_total*8/1024/ADM_prc_all/comm_call_count
       write(ADM_log_fid,'(a,1f11.4)') "ave_comm_count    =",ave_comm_count/ADM_prc_all/comm_call_count
       write(ADM_log_fid,*) "---------  minimum  --------"
       write(ADM_log_fid,'(a,1e11.4)') "min_time_pre      =",1000*min_time_pre/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_bar1     =",1000*min_time_bar1/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_recv     =",1000*min_time_recv/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_sbuf     =",1000*min_time_sbuf/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_send     =",1000*min_time_send/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_copy     =",1000*min_time_copy/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_wait     =",1000*min_time_wait/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_rbuf     =",1000*min_time_rbuf/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_copy_sgp =",1000*min_time_copy_sgp/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_bar2     =",1000*min_time_bar2/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "min_time_total    =",1000*min_time_total/comm_call_count
       write(ADM_log_fid,'(a,1f11.4)') "min_size_total    =",min_size_total*8/1024/comm_call_count
       write(ADM_log_fid,'(a,1f11.4)') "min_comm_count    =",min_comm_count/comm_call_count
       write(ADM_log_fid,*) "---------  maximum  --------"
       write(ADM_log_fid,'(a,1e11.4)') "max_time_pre      =",1000*max_time_pre/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_bar1     =",1000*max_time_bar1/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_recv     =",1000*max_time_recv/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_sbuf     =",1000*max_time_sbuf/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_send     =",1000*max_time_send/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_copy     =",1000*max_time_copy/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_wait     =",1000*max_time_wait/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_rbuf     =",1000*max_time_rbuf/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_copy_sgp =",1000*max_time_copy_sgp/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_bar2     =",1000*max_time_bar2/comm_call_count
       write(ADM_log_fid,'(a,1e11.4)') "max_time_total    =",1000*max_time_total/comm_call_count
       write(ADM_log_fid,'(a,1f11.4)') "max_size_total    =",max_size_total*8/1024/comm_call_count
       write(ADM_log_fid,'(a,1f11.4)') "max_comm_count    =",max_comm_count/comm_call_count
       write(ADM_log_fid,*) "----------------------------"
    endif
    !
    !call ADM_proc_stop
    return
    !
  end subroutine output_time
  !-----------------------------------------------------------------------------
  subroutine COMM_data_transfer(&
       var,var_pl,              & !--- INOUT : variables comminicated.
       trn,                     & !--- IN : commnucation flag for each variable
       hallo_num)                 !--- IN : number of hallo
    use mod_adm, only :      &
         ADM_proc_stop,      & ! [add] C.Kodama 2011.04.26
         ADM_vlink_nmax,     &
         ADM_lall,           &
         ADM_comm_world, &
         ADM_LOG_FID,        &
         ADM_kall
    use mod_cnst, only: &
       CNST_undef
    implicit none

    real(8),intent(inout)::var(:,:,:,:)
    real(8),intent(inout)::var_pl(:,:,:,:)
    !
    logical,intent(in),optional::trn(:)
    !
    integer,intent(in),optional::hallo_num
    !
    integer :: shp(4)
    integer :: varmax
    integer :: cmax
    integer::halo

    integer :: acount
    integer :: areq(2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4))
    integer :: stat(mpi_status_size,2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4))
    integer :: ierr
    !
    integer :: k,m,n
    integer :: nr,nc,ns
    integer :: sl,so,sb,ss
    integer :: rl,ro,rb,rs
    integer :: cl,scl,cs
!    integer :: i_dbg !iga
    logical :: dbg_ierr !iga
    !---------------------------------------------------------------------------

    if ( opt_comm_barrier ) then
       call DEBUG_rapstart('Node imbalance adjustment')
       call MPI_BARRIER( ADM_comm_world, ierr )
       call DEBUG_rapend  ('Node imbalance adjustment')
    endif

    call DEBUG_rapstart('COMM data_transfer')

    ! --- pre-process
    !
    comm_call_count=comm_call_count+1
    !
    t(0)=mpi_wtime()
    !
    shp=shape(var)
    kmax = shp(2)
    !
    if (present(trn)) then
       varmax=0
!cdir novector
       do n=1,shp(4)
          if(trn(n)) then
             varmax=varmax+1
             clist(varmax)=n
          endif
       enddo
    else
       varmax=shp(4)
!cdir novector
       do n=1,varmax
          clist(n)=n
       enddo
    endif
    !
    if (present(hallo_num)) then
       halo=hallo_num
    else
       halo=1
    endif

    if (halo.ne.1) then
       write(ADM_log_fid,*) 'halo=',halo
    endif

    !
    cmax=varmax*kmax
    ! [Add] 07/11/07 T.Mitsui for check of varmax
    if( opt_check_varmax ) then
       ! <- [rep] C.Kodama 2011.04.26
       if ( cmax > max_varmax * ADM_kall ) then
          write(ADM_LOG_FID,*)  'error: cmax >  max_varmax * ADM_kall, stop!'
          write(ADM_LOG_FID,*)  'cmax=', cmax, 'max_varmax*ADM_kall=', max_varmax*ADM_kall
          call ADM_proc_stop
       end if
!       equiv_varmax = real( varmax*kmax )/real( ADM_kall ) ! assuming variables are all 3-Dimension
!       diag_varmax  = max( equiv_varmax, diag_varmax )     ! diagnose max value
!       write(ADM_LOG_FID,'(a)'  )   ' *** max_varmax, varmax, kmax, equivalent varmax, diagnosed max_varmax '
!       write(ADM_LOG_FID,'(5x, 3i8, 2f18.1)') max_varmax, varmax, kmax, equiv_varmax, diag_varmax
       ! ->
    end if
    !
    t(1)=mpi_wtime()
    time_pre=time_pre+(t(1)-t(0))
    !call mpi_barrier(ADM_comm_world,ierr)
    t(2)=mpi_wtime()
    time_bar1=time_bar1+(t(2)-t(1))
    !-----------------------------------------
    ! call mpi_isend
    !-----------------------------------------

!!    comm_dbg_recvbuf(:,:,:)=CNST_UNDEF
!!    comm_dbg_sendbuf(:,:,:)=CNST_UNDEF
!!    comm_dbg_recvbuf(:,:,1)=recvbuf(:,:)
!!    comm_dbg_sendbuf(:,:,1)=sendbuf(:,:)

    if (opt_comm_dbg) then !iga
       sendbuf(:,:)= dbg_sendbuf_init
       recvbuf(:,:)= dbg_recvbuf_init
       if (somax(halo)>(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4)) then
          write(ADM_log_fid,*) 'somax inconsistency'
          write(ADM_log_fid,*) somax(halo),(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4)
       endif
    endif


    do ro=1,romax(halo)
       !
       if (opt_comm_dbg) then
          if (rsize(ro,halo)*cmax > size(recvbuf(:,ro)) ) then
             write(ADM_log_fid,*) 'ro',ro,'buf-dim=',shape(recvbuf),'rsize=',rsize(ro,halo),'cmax=',cmax
          endif
       endif
       !
       call mpi_irecv(recvbuf(1,ro)              &
            ,rsize(ro,halo)*cmax        &
            ,mpi_double_precision       &
            ,sourcerank(ro,halo)        &
            ,recvtag(ro,halo)           &
            ,ADM_comm_world         &
            ,areq(ro)                   &
            ,ierr)
       !
       if (opt_comm_dbg) then
          if (ierr.ne.0) then
             write(ADM_log_fid,*) 'mpi_irecv info start=='
             write(ADM_log_fid,*) 'ierr=',ierr
             write(ADM_log_fid,*) 'ro=',ro
             write(ADM_log_fid,*) 'rsize=',rsize(ro,halo)
             write(ADM_log_fid,*) 'cmax=',cmax
             write(ADM_log_fid,*) 'areq=',areq(ro)
             write(ADM_log_fid,*) 'mpi_irecv info end=='
          endif
          ! [fix] 12/03/26 T.Seiki
!!$       dbg_areq_save(:,1)=areq(:)
          dbg_areq_save(ro,1)=areq(ro)
       endif
    enddo

    t(3)=mpi_wtime()
    time_recv=time_recv+(t(3)-t(2))
    !-----------------------------------------
    !  var -> sendbuf
    !-----------------------------------------
    do so=1,somax(halo)



       t(4)=mpi_wtime()
       do ns=1,nsmax(so,halo)
          ss=sendinfo(SIZE_COMM,ns,so,halo)
          sl=sendinfo(LRGNID_COMM,ns,so,halo)
          sb=sendinfo(BASE_COMM,ns,so,halo)*cmax
!=org=!cdir novector
!=org=          do n=1,ss
!=org=!cdir unroll=3
!=org=             do m=1,varmax
!=org=!!cdir shortloop
!=org=                do k=1,kmax
!=org=                   sendbuf(k+(m-1)*kmax+(n-1)*cmax+sb,so) &
!=org=                        =var(sendlist(n,ns,so,halo),k,sl,clist(m))
!=org=                enddo
!=org=             enddo
!=org=!             if (ADM_prc_me==1) write(ADM_log_fid,*) 'n',n,sendlist(n,ns,so,halo)
!=org=          enddo

          do m=1,varmax
!cdir outerunroll=8
             do k=1,kmax
                do n=1,ss
                   sendbuf(n+(k-1)*ss+(m-1)*ss*kmax+sb,so) &
                        =var(sendlist(n,ns,so,halo),k,sl,clist(m))
                enddo
             enddo
          enddo


       enddo
       !-----------------------------------------
       !
       !-----------------------------------------
       !  var_pl -> sendbuf
       !-----------------------------------------
       do ns=1,nsmax_pl(so,halo)
          ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
          sl=sendinfo_pl(LRGNID_COMM,ns,so,halo)
          sb=sendinfo_pl(BASE_COMM,ns,so,halo)*cmax
!=org=!cdir novector
!=org=          do n=1,ss
!=org=!cdir unroll=3
!=org=             do m=1,varmax
!=org=!!cdir shortloop
!=org=                do k=1,kmax
!=org=                   sendbuf(k+(m-1)*kmax+(n-1)*cmax+sb,so) &
!=org=                        =var_pl(sendlist_pl(n,ns,so,halo),k,sl,clist(m))
!=org=                enddo
!=org=             enddo
!=org=          enddo
          do m=1,varmax
!cdir outerunroll=8
             do k=1,kmax
                do n=1,ss
                   sendbuf(n+(k-1)*ss+(m-1)*ss*kmax+sb,so) &
                        =var_pl(sendlist_pl(n,ns,so,halo),k,sl,clist(m))
                enddo
             enddo
          enddo

       enddo
       t(5)=mpi_wtime()
       time_sbuf=time_sbuf+(t(5)-t(4))


!       write(ADM_log_fid,*) 'send count=',i_dbg,'prc=',adm_prc_me
       !-----------------------------------------
       !
       !-----------------------------------------
       ! call mpi_isend
       !-----------------------------------------

       if (opt_comm_dbg) then
          if (ssize(so,halo)*cmax>size(sendbuf(:,so))) then
             write(ADM_log_fid,*) 'so',so,'buf-dim=',shape(sendbuf),'rsize=',ssize(so,halo),'cmax=',cmax
          endif
       endif

!       write(*,*) 'me=',ADM_prc_me,'sendtag',sendtag(:,:)
!       write(*,*) 'sendbuf',sendbuf(:,:)
!       write(*,*) 'me=',ADM_prc_me,'destrank',destrank(:,:)


       call mpi_isend(sendbuf(1,so)              &
            ,ssize(so,halo)*cmax        &
            ,mpi_double_precision       &
            ,destrank(so,halo)          &
            ,sendtag(so,halo)           &
            ,ADM_comm_world             &
            ,areq(so+romax(halo))       &
            ,ierr)
       if (opt_comm_dbg) then
          if (ierr.ne.0) then
             write(ADM_log_fid,*) 'mpi_isend info start=='
             write(ADM_log_fid,*) 'ierr=',ierr
             write(ADM_log_fid,*) 'so=',so
             write(ADM_log_fid,*) 'ssize=',ssize(so,halo)
             write(ADM_log_fid,*) 'cmax=',cmax
             write(ADM_log_fid,*) 'areq=',areq(so+romax(halo))
             write(ADM_log_fid,*) 'mpi_isend info end=='
          endif
          ! [fix] 12/03/26 T.Seiki
!!$       dbg_areq_save(:,2)=areq(:)
          dbg_areq_save(so+romax(halo),2)=areq(so+romax(halo))
       endif


       t(6)=mpi_wtime()
       time_send=time_send+(t(6)-t(5))
       size_total=size_total+ssize(so,halo)*cmax
       comm_count=comm_count+1
    enddo !loop so
    !-----------------------------------------
    !
    t(7)=mpi_wtime()
    !---------------------------------------------------
    !  var -> var (region to region copy in same rank)
    !---------------------------------------------------
    do nc=1,ncmax_r2r(halo)
       cs=copyinfo_r2r(SIZE_COPY,nc,halo)
       cl=copyinfo_r2r(LRGNID_COPY,nc,halo)
       scl=copyinfo_r2r(SRC_LRGNID_COPY,nc,halo)
!=org=!cdir novector
!=org=       do n=1,cs
!=org=!cdir unroll=3
!=org=          do m=1,varmax
!=org=!!cdir shortloop
!=org=             do k=1,kmax
!=org=                var(recvlist_r2r(n,nc,halo),k,cl ,clist(m)) &
!=org=                     =var(sendlist_r2r(n,nc,halo),k,scl,clist(m))
!=org=             enddo
!=org=          enddo
!=org=       enddo
       do m=1,varmax
!cdir outerunroll=8
          do k=1,kmax
             do n=1,cs
                var(recvlist_r2r(n,nc,halo),k,cl ,clist(m)) &
                     =var(sendlist_r2r(n,nc,halo),k,scl,clist(m))
             enddo
          enddo
       enddo

    enddo
    !------------------------------------------
    !
    !------------------------------------------
    !  var -> var_pl ( data copy in same rank)
    !------------------------------------------
    do nc=1,ncmax_r2p(halo)
       cs=copyinfo_r2p(SIZE_COPY,nc,halo)
       cl=copyinfo_r2p(LRGNID_COPY,nc,halo)
       scl=copyinfo_r2p(SRC_LRGNID_COPY,nc,halo)
!=org=!cdir novector
!=org=       do n=1,cs
!=org=!cdir unroll=3
!=org=          do m=1,varmax
!=org=!!cdir shortloop
!=org=             do k=1,kmax
!=org=                var_pl(recvlist_r2p(n,nc,halo),k,cl,clist(m)) &
!=org=                     =var(sendlist_r2p(n,nc,halo),k,scl,clist(m))
!=org=             enddo
!=org=          enddo
!=org=       enddo
       do m=1,varmax
!cdir outerunroll=8
          do k=1,kmax
             do n=1,cs
                var_pl(recvlist_r2p(n,nc,halo),k,cl,clist(m)) &
                     =var(sendlist_r2p(n,nc,halo),k,scl,clist(m))
             enddo
          enddo
       enddo

    enddo
    !-----------------------------------------
    !
    !-----------------------------------------
    !  var_pl -> var (data copy in same rank)
    !-----------------------------------------
    do nc=1,ncmax_p2r(halo)
       cs=copyinfo_p2r(SIZE_COPY,nc,halo)
       cl=copyinfo_p2r(LRGNID_COPY,nc,halo)
       scl=copyinfo_p2r(SRC_LRGNID_COPY,nc,halo)
!=org=!cdir novector
!=org=       do n=1,cs
!=org=!cdir unroll=3
!=org=          do m=1,varmax
!=org=!!cdir shortloop
!=org=             do k=1,kmax
!=org=                var(recvlist_p2r(n,nc,halo),k,cl,clist(m)) &
!=org=                     =var_pl(sendlist_p2r(n,nc,halo),k,scl,clist(m))
!=org=             enddo
!=org=          enddo
!=org=       enddo
       do m=1,varmax
!cdir outerunroll=8
          do k=1,kmax
             do n=1,cs
                var(recvlist_p2r(n,nc,halo),k,cl,clist(m)) &
                     =var_pl(sendlist_p2r(n,nc,halo),k,scl,clist(m))
             enddo
          enddo
       enddo

    enddo
    !-----------------------------------------
    !
    !-----------------------------------------
    t(8)=mpi_wtime()
    time_copy=time_copy+(t(8)-t(7))
    acount=romax(halo)+somax(halo)
    call mpi_waitall(acount,areq,stat,ierr)

    t(9)=mpi_wtime()
    time_wait=time_wait+(t(9)-t(8))

    if (opt_comm_dbg) then  !================== dbg start
       if (ierr.ne.0) then
          write(ADM_log_fid,*) 'mpi_wait info start=='
          write(ADM_log_fid,*) 'ierr=',ierr
          write(ADM_log_fid,*) 'acount=',acount
          write(ADM_log_fid,*) 'stat=',stat(:,1:acount)
          write(ADM_log_fid,*) 'areq=',areq(1:acount)
          write(ADM_log_fid,*) 'mpi_wait info end=='
       endif
       if (acount > size(areq)) then
          write(ADM_log_fid,*) 'acount, size(areq)',acount,size(areq)
       endif
       ! [fix] 12/03/26 T.Seiki
!!$    dbg_areq_save(:,3)=areq(:)
       dbg_areq_save(1:acount,3)=areq(1:acount)
    endif  !=================================== dbg end

    if (opt_comm_dbg) then  !================== dbg start
       dbg_ierr=.false.
       do ro=1,romax(halo)
          do n=1,rsize(ro,halo)*cmax
             if (recvbuf(n,ro) == dbg_recvbuf_init) then
                dbg_ierr=.true.
             endif
             if (recvbuf(n,ro) == dbg_sendbuf_init) then
                dbg_ierr=.true.
             endif
          enddo
       enddo
       if (dbg_ierr) then
          write(ADM_log_fid,*) 'communication is not completed!'
          write(ADM_log_fid,*) 'dbg_tcount=',dbg_tcount
          write(ADM_log_fid,*) 'romax=',romax,'rsize=',rsize(ro,halo),'cmax=',cmax
          do ro=1,romax(halo)
             do n=1,rsize(ro,halo)*cmax
                if (recvbuf(n,ro) == dbg_recvbuf_init) then
                   write(ADM_log_fid,*) 'n=',n,'ro=',ro,'recvbuf=',recvbuf(n,ro)
                endif
             enddo
          enddo
          write(ADM_log_fid,*) 'areq after irecv:',dbg_areq_save(1:acount,1)
          write(ADM_log_fid,*) 'areq after isend',dbg_areq_save(1:acount,2)
          write(ADM_log_fid,*) 'areq after waitall',dbg_areq_save(1:acount,3)
!          write(ADM_log_fid,*) 'areq after barrier',dbg_areq_save(1:acount,4)

          write(ADM_log_fid,*) 'ierr of mpi_waitall=',ierr
          write(ADM_log_fid,*) 'acount of mpi_waitall=',acount
          write(ADM_log_fid,*) 'stat of mpi_waitall=',stat(:,1:acount)

       endif
       dbg_tcount=dbg_tcount+1
    endif  !=================================== dbg end


    if  (opt_comm_barrier) then
       call mpi_barrier(ADM_comm_world,ierr)
       if (ierr.ne.0) then
          write(ADM_log_fid,*) 'mpi_barriert info start=='
          write(ADM_log_fid,*) 'ierr=',ierr
          write(ADM_log_fid,*) 'mpi_barrier info end=='
       endif
!       if (opt_comm_dbg) then
!          dbg_areq_save(:,4)=areq(:)
!       endif
    endif


    !-----------------------------------------
    !
!    i_dbg=0 !iga
    do ro=1,romax(halo)
       !-----------------------------------------
       !  recvbuf -> var ( recieve in region )
       !-----------------------------------------
       do nr=1,nrmax(ro,halo)
          rs=recvinfo(SIZE_COMM,nr,ro,halo)
          rl=recvinfo(LRGNID_COMM,nr,ro,halo)
          rb=recvinfo(BASE_COMM,nr,ro,halo)*cmax
!=org=!cdir novector
!=org=          do n=1,rs
!=org=!cdir unroll=3
!=org=             do m=1,varmax
!=org=!!cdir shortloop
!=org=                do k=1,kmax
!=org=                   var(recvlist(n,nr,ro,halo),k,rl,clist(m)) &
!=org=                        =recvbuf(k+(m-1)*kmax+(n-1)*cmax+rb,ro)
!=org=                enddo
!=org=             enddo
!=org=!             if (ADM_prc_me==2) write(ADM_log_fid,*) 'ro,n,nr,',ro,n,nr,recvlist(n,nr,ro,halo)
!=org=          enddo
          do m=1,varmax
!cdir outerunroll=8
             do k=1,kmax
                do n=1,rs
                   var(recvlist(n,nr,ro,halo),k,rl,clist(m)) &
                        =recvbuf(n+(k-1)*rs+(m-1)*rs*kmax+rb,ro)
                enddo
             enddo
!             if (ADM_prc_me==2) write(ADM_log_fid,*) 'ro,n,nr,',ro,n,nr,recvlist(n,nr,ro,halo)
          enddo
!          i_dbg=i_dbg+max(rs,0)
       enddo

       !-----------------------------------------
       !
       !-----------------------------------------
       !  recvbuf -> var_pl ( recieve in pole )
       !-----------------------------------------
       do nr=1,nrmax_pl(ro,halo)
          rs=recvinfo_pl(SIZE_COMM,nr,ro,halo)
          rl=recvinfo_pl(LRGNID_COMM,nr,ro,halo)
          rb=recvinfo_pl(BASE_COMM,nr,ro,halo)*cmax
!=org=!cdir novector
!=org=          do n=1,rs
!=org=!cdir unroll=3
!=org=             do m=1,varmax
!=org=!!cdir shortloop
!=org=                do k=1,kmax
!=org=                   var_pl(recvlist_pl(n,nr,ro,halo),k,rl,clist(m)) &
!=org=                        =recvbuf(k+(m-1)*kmax+(n-1)*cmax+rb,ro)
!=org=                enddo
!=org=             enddo
!=org=          enddo
          do m=1,varmax
!cdir outerunroll=8
             do k=1,kmax
                do n=1,rs
                   var_pl(recvlist_pl(n,nr,ro,halo),k,rl,clist(m)) &
                        =recvbuf(n+(k-1)*rs+(m-1)*rs*kmax+rb,ro)
                enddo
             enddo
          enddo

       enddo
    enddo !loop ro
    t(10)=mpi_wtime()
    time_rbuf=time_rbuf+(t(10)-t(9))

!    write(ADM_log_fid,*) 'recv count=',i_dbg,'prc=',adm_prc_me

    !-----------------------------------------
    !
    !-----------------------------------------
    !  copy data around singular point
    !-----------------------------------------
    do nc=1,ncmax_sgp(halo)
       cs=copyinfo_sgp(SIZE_COPY,nc,halo)
       cl=copyinfo_sgp(LRGNID_COPY,nc,halo)
       scl=copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)
!=org=!cdir novector
!=org=       do n=1,cs
!=org=!cdir unroll=3
!=org=          do m=1,varmax
!=org=!!cdir shortloop
!=org=             do k=1,kmax
!=org=                var(recvlist_sgp(n,nc,halo),k,cl ,clist(m)) &
!=org=                     =var(sendlist_sgp(n,nc,halo),k,scl,clist(m))
!=org=             enddo
!=org=          enddo
!=org=       enddo

       do m=1,varmax
!cdir outerunroll=8
          do k=1,kmax
             do n=1,cs
                var(recvlist_sgp(n,nc,halo),k,cl ,clist(m)) &
                     =var(sendlist_sgp(n,nc,halo),k,scl,clist(m))
             enddo
          enddo
       enddo

    enddo
    t(11)=mpi_wtime()
    time_copy_sgp=time_copy_sgp+(t(11)-t(10))
    !!
    !call mpi_barrier(ADM_comm_world,ierr)
    t(12)=mpi_wtime()
    time_bar2=time_bar2+(t(12)-t(11))
    time_total=time_total+(t(12)-t(0))

    call DEBUG_rapend('COMM data_transfer')

    return
  end subroutine COMM_data_transfer
  !-----------------------------------------------------------------------------
  subroutine COMM_var( &
       var,            & !--- INOUT : variables
       var_pl,         & !--- INOUT : variables at poles
       KNUM,           & !--- IN : number of layers
       NNUM            ) !--- IN : number of variables
    use mod_adm, only :     &
         ADM_gall_pl,       &
         ADM_lall_pl,       &
         ADM_prc_me,        &
         ADM_prc_pl,        &
         ADM_KNONE,         &
         ADM_gall,          &
         ADM_gmin,          &
         ADM_gmax,          &
         ADM_lall,          &
         ADM_gall_1d,       &
         ADM_rgnid_npl_mng, &
         ADM_rgnid_spl_mng, &
         ADM_prc_tab,       &
         ADM_prc_npl,       &
         ADM_prc_spl,       &
         ADM_GSLF_PL,       &
         ADM_NPL,           &
         ADM_SPL,           &
         ADM_COMM_WORLD,&
         ADM_rgn2prc
    !
    implicit none

    integer, intent(in) :: KNUM
    integer, intent(in) :: NNUM
    real(8), intent(inout) :: var(ADM_gall,KNUM,ADM_lall,NNUM)
    real(8), intent(inout) :: var_pl(ADM_gall_pl,KNUM,ADM_lall_pl,NNUM)

    integer :: ireq(4),istat(MPI_STATUS_SIZE),ierr
    integer :: l
    !
    real(8) :: v_npl(KNUM,NNUM)
    real(8) :: v_spl(KNUM,NNUM)
!!! FIX 20111024:s
    real(8) :: v_npl_recv(KNUM,NNUM)
    real(8) :: v_spl_recv(KNUM,NNUM)
!!! FIX 20111024:e
    !
    integer :: rgnid
    !
    integer :: i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( opt_comm_barrier ) then
       call DEBUG_rapstart('Node imbalance adjustment')
       call MPI_BARRIER( ADM_comm_world, ierr )
       call DEBUG_rapend  ('Node imbalance adjustment')
    endif

    call DEBUG_rapstart('COMM var')

    if ( comm_pl ) then ! T.Ohno 110721

       !--- recv pole value
       if(ADM_prc_me==ADM_prc_npl) then
          call MPI_IRECV(            &
               v_npl_recv,                & !--- starting address
               KNUM * NNUM,          & !--- number of array
               MPI_DOUBLE_PRECISION, & !--- type
               ADM_rgn2prc(ADM_rgnid_npl_mng)-1, & !--- source rank
               ADM_NPL,              & !--- tag
               ADM_COMM_WORLD,       & !--- world
               ireq(3),                 & !--- request id
               ierr)                   !--- error id
       end if
       !
       if(ADM_prc_me==ADM_prc_spl) then
          call MPI_IRECV(            &
               v_spl_recv,                & !--- starting address
               KNUM * NNUM,          & !--- number of array
               MPI_DOUBLE_PRECISION, & !--- type
               ADM_rgn2prc(ADM_rgnid_spl_mng)-1, & !--- srouce rank
               ADM_SPL,              & !--- tag
               ADM_COMM_WORLD,       & !--- world
               ireq(4),                 & !--- request id
               ierr)                   !--- error id
       end if

       !--- send pole value
       do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          !
          !--- north pole
          if(rgnid==ADM_rgnid_npl_mng) then
             v_npl(1:KNUM,1:NNUM)&
                  = var(suf(ADM_gmin,ADM_gmax+1),1:KNUM,l,1:NNUM)
             call MPI_ISEND(            &
                  v_npl,                & !--- starting address
                  KNUM * NNUM,          & !--- number of array
                  MPI_DOUBLE_PRECISION, & !--- type
                  ADM_prc_npl-1,        & !--- dest rank
                  ADM_NPL,              & !--- tag
                  ADM_COMM_WORLD,       & !--- world
                  ireq(1),                 & !--- request id
                  ierr)                   !--- error id
          end if
          !
          !--- south pole
          if(rgnid==ADM_rgnid_spl_mng) then
             v_spl(1:KNUM,1:NNUM)&
                  = var(suf(ADM_gmax+1,ADM_gmin),1:KNUM,l,1:NNUM)
             call MPI_ISEND(            &
                  v_spl,                & !--- starting address
                  KNUM * NNUM,          & !--- number of array
                  MPI_DOUBLE_PRECISION, & !--- type
                  ADM_prc_spl-1,        & !--- dest rank
                  ADM_SPL,              & !--- tag
                  ADM_COMM_WORLD,       & !--- world
                  ireq(2),                 & !--- request id
                  ierr)                   !--- error id
          end if
       end do

       do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          !
          !--- north pole
          if(rgnid==ADM_rgnid_npl_mng) then
             call MPI_WAIT(ireq(1),istat,ierr)
          end if
          !
          !--- south pole
          if(rgnid==ADM_rgnid_spl_mng) then
             call MPI_WAIT(ireq(2),istat,ierr)
          end if
       end do

       if(ADM_prc_me==ADM_prc_npl) then
          call MPI_WAIT(ireq(3),istat,ierr)
          var_pl(ADM_GSLF_PL,1:KNUM,ADM_NPL,1:NNUM) &
               = v_npl_recv(1:KNUM,1:NNUM)
       end if
       !
       if(ADM_prc_me==ADM_prc_spl) then
          call MPI_WAIT(ireq(4),istat,ierr)
          var_pl(ADM_GSLF_PL,:,ADM_SPL,:) = v_spl_recv(:,:)
       end if

    end if

    call COMM_data_transfer(var,var_pl)

    var(suf(ADM_gall_1d,1),:,:,:) = var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
    var(suf(1,ADM_gall_1d),:,:,:) = var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    call DEBUG_rapend('COMM var')

    return
  end subroutine COMM_var

  !-----------------------------------------------------------------------------
  ! [add] T.Ohno 110721
  subroutine COMM_data_transfer_nopl(&
       var,                     & !--- INOUT : variables comminicated.
       trn,                     & !--- IN : commnucation flag for each variable
       hallo_num)                 !--- IN : number of hallo
    use mod_adm, only :      &
         ADM_proc_stop,      & ! [add] C.Kodama 2011.04.26
         ADM_vlink_nmax,     &
         ADM_lall,           &
         ADM_comm_world, &
         ADM_kall
    use mod_cnst, only : &
       CNST_undef
    implicit none
    !
    real(8),intent(inout)::var(:,:,:,:)
    !
    logical,intent(in),optional::trn(:)
    !
    integer,intent(in),optional::hallo_num
    !
    integer :: shp(4)
    integer :: varmax
    integer :: cmax
    integer::halo

    integer :: acount
    integer :: areq(2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4))
    integer :: stat(mpi_status_size,2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4))
    integer :: ierr
    !
    integer :: k,m,n
    integer :: nr,nc,ns
    integer :: sl,so,sb,ss
    integer :: rl,ro,rb,rs
    integer :: cl,scl,cs
!    integer :: i_dbg !iga
    logical :: dbg_ierr !iga
  !=============================================================================
    !
    ! --- pre-process
    !
    comm_call_count=comm_call_count+1
    !
    t(0)=mpi_wtime()
    !
    shp=shape(var)
    kmax = shp(2)
    !
    if (present(trn)) then
       varmax=0
!cdir novector
       do n=1,shp(4)
          if(trn(n)) then
             varmax=varmax+1
             clist(varmax)=n
          endif
       enddo
    else
       varmax=shp(4)
!cdir novector
       do n=1,varmax
          clist(n)=n
       enddo
    endif
    !
    if (present(hallo_num)) then
       halo=hallo_num
    else
       halo=1
    endif

    if (halo.ne.1) then
       write(ADM_log_fid,*) 'halo=',halo
    endif

    !
    cmax=varmax*kmax
    ! [Add] 07/11/07 T.Mitsui for check of varmax
    if( opt_check_varmax ) then
       ! <- [rep] C.Kodama 2011.04.26
       if ( cmax > max_varmax * ADM_kall ) then
          write(ADM_LOG_FID,*)  'error: cmax >  max_varmax * ADM_kall, stop!'
          write(ADM_LOG_FID,*)  'cmax=', cmax, 'max_varmax*ADM_kall=', max_varmax*ADM_kall
          call ADM_proc_stop
       end if
       ! ->
    end if
    !
    t(1)=mpi_wtime()
    time_pre=time_pre+(t(1)-t(0))
    t(2)=mpi_wtime()
    time_bar1=time_bar1+(t(2)-t(1))
    !-----------------------------------------
    ! call mpi_isend
    !-----------------------------------------
    if (opt_comm_dbg) then !iga
       sendbuf(:,:)= dbg_sendbuf_init
       recvbuf(:,:)= dbg_recvbuf_init
       if (somax(halo)>(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4)) then
          write(ADM_log_fid,*) 'somax inconsistency'
          write(ADM_log_fid,*) somax(halo),(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4)
       endif
    endif


    do ro=1,romax(halo)
       !
       if (opt_comm_dbg) then
          if (rsize(ro,halo)*cmax > size(recvbuf(:,ro)) ) then
             write(ADM_log_fid,*) 'ro',ro,'buf-dim=',shape(recvbuf),'rsize=',rsize(ro,halo),'cmax=',cmax
          endif
       endif
       !
       call mpi_irecv(recvbuf(1,ro)              &
            ,rsize(ro,halo)*cmax        &
            ,mpi_double_precision       &
            ,sourcerank(ro,halo)        &
            ,recvtag(ro,halo)           &
            ,ADM_comm_world         &
            ,areq(ro)                   &
            ,ierr)
       !
       if (opt_comm_dbg) then
          if (ierr.ne.0) then
             write(ADM_log_fid,*) 'mpi_irecv info start=='
             write(ADM_log_fid,*) 'ierr=',ierr
             write(ADM_log_fid,*) 'ro=',ro
             write(ADM_log_fid,*) 'rsize=',rsize(ro,halo)
             write(ADM_log_fid,*) 'cmax=',cmax
             write(ADM_log_fid,*) 'areq=',areq(ro)
             write(ADM_log_fid,*) 'mpi_irecv info end=='
          endif
          ! [fix] 12/03/26 T.Seiki
!!$       dbg_areq_save(:,1)=areq(:)
          dbg_areq_save(ro,1)=areq(ro)
       endif
    enddo

    t(3)=mpi_wtime()
    time_recv=time_recv+(t(3)-t(2))
    !-----------------------------------------
    !  var -> sendbuf
    !-----------------------------------------
    do so=1,somax(halo)



       t(4)=mpi_wtime()
       do ns=1,nsmax(so,halo)
          ss=sendinfo(SIZE_COMM,ns,so,halo)
          sl=sendinfo(LRGNID_COMM,ns,so,halo)
          sb=sendinfo(BASE_COMM,ns,so,halo)*cmax
          do m=1,varmax
!cdir outerunroll=8
             do k=1,kmax
                do n=1,ss
                   sendbuf(n+(k-1)*ss+(m-1)*ss*kmax+sb,so) &
                        =var(sendlist(n,ns,so,halo),k,sl,clist(m))
                enddo
             enddo
          enddo


       enddo
       !-----------------------------------------
       t(5)=mpi_wtime()
       time_sbuf=time_sbuf+(t(5)-t(4))


       !-----------------------------------------
       !
       !-----------------------------------------
       ! call mpi_isend
       !-----------------------------------------

       if (opt_comm_dbg) then
          if (ssize(so,halo)*cmax>size(sendbuf(:,so))) then
             write(ADM_log_fid,*) 'so',so,'buf-dim=',shape(sendbuf),'rsize=',ssize(so,halo),'cmax=',cmax
          endif
       endif

       call mpi_isend(sendbuf(1,so)              &
            ,ssize(so,halo)*cmax        &
            ,mpi_double_precision       &
            ,destrank(so,halo)          &
            ,sendtag(so,halo)           &
            ,ADM_comm_world             &
            ,areq(so+romax(halo))       &
            ,ierr)
       if (opt_comm_dbg) then
          if (ierr.ne.0) then
             write(ADM_log_fid,*) 'mpi_isend info start=='
             write(ADM_log_fid,*) 'ierr=',ierr
             write(ADM_log_fid,*) 'so=',so
             write(ADM_log_fid,*) 'ssize=',ssize(so,halo)
             write(ADM_log_fid,*) 'cmax=',cmax
             write(ADM_log_fid,*) 'areq=',areq(so+romax(halo))
             write(ADM_log_fid,*) 'mpi_isend info end=='
          endif
          ! [fix] 12/03/26 T.Seiki
!!$       dbg_areq_save(:,2)=areq(:)
          dbg_areq_save(so+romax(halo),2)=areq(so+romax(halo))
       endif


       t(6)=mpi_wtime()
       time_send=time_send+(t(6)-t(5))
       size_total=size_total+ssize(so,halo)*cmax
       comm_count=comm_count+1
    enddo !loop so
    !-----------------------------------------
    !
    t(7)=mpi_wtime()
    !---------------------------------------------------
    !  var -> var (region to region copy in same rank)
    !---------------------------------------------------
    do nc=1,ncmax_r2r(halo)
       cs=copyinfo_r2r(SIZE_COPY,nc,halo)
       cl=copyinfo_r2r(LRGNID_COPY,nc,halo)
       scl=copyinfo_r2r(SRC_LRGNID_COPY,nc,halo)
       do m=1,varmax
!cdir outerunroll=8
          do k=1,kmax
             do n=1,cs
                var(recvlist_r2r(n,nc,halo),k,cl ,clist(m)) &
                     =var(sendlist_r2r(n,nc,halo),k,scl,clist(m))
             enddo
          enddo
       enddo

    enddo
    !-----------------------------------------
    !
    !-----------------------------------------
    t(8)=mpi_wtime()
    time_copy=time_copy+(t(8)-t(7))
    acount=romax(halo)+somax(halo)
    call mpi_waitall(acount,areq,stat,ierr)

    t(9)=mpi_wtime()
    time_wait=time_wait+(t(9)-t(8))

    if (opt_comm_dbg) then  !================== dbg start
       if (ierr.ne.0) then
          write(ADM_log_fid,*) 'mpi_wait info start=='
          write(ADM_log_fid,*) 'ierr=',ierr
          write(ADM_log_fid,*) 'acount=',acount
          write(ADM_log_fid,*) 'stat=',stat(:,1:acount)
          write(ADM_log_fid,*) 'areq=',areq(1:acount)
          write(ADM_log_fid,*) 'mpi_wait info end=='
       endif
       if (acount > size(areq)) then
          write(ADM_log_fid,*) 'acount, size(areq)',acount,size(areq)
       endif
       ! [fix] 12/03/26 T.Seiki
!!$    dbg_areq_save(:,3)=areq(:)
       dbg_areq_save(1:acount,3)=areq(1:acount)
    endif  !=================================== dbg end

    if (opt_comm_dbg) then  !================== dbg start
       dbg_ierr=.false.
       do ro=1,romax(halo)
          do n=1,rsize(ro,halo)*cmax
             if (recvbuf(n,ro) == dbg_recvbuf_init) then
                dbg_ierr=.true.
             endif
             if (recvbuf(n,ro) == dbg_sendbuf_init) then
                dbg_ierr=.true.
             endif
          enddo
       enddo
       if (dbg_ierr) then
          write(ADM_log_fid,*) 'communication is not completed!'
          write(ADM_log_fid,*) 'dbg_tcount=',dbg_tcount
          write(ADM_log_fid,*) 'romax=',romax,'rsize=',rsize(ro,halo),'cmax=',cmax
          do ro=1,romax(halo)
             do n=1,rsize(ro,halo)*cmax
                if (recvbuf(n,ro) == dbg_recvbuf_init) then
                   write(ADM_log_fid,*) 'n=',n,'ro=',ro,'recvbuf=',recvbuf(n,ro)
                endif
             enddo
          enddo
          write(ADM_log_fid,*) 'areq after irecv:',dbg_areq_save(1:acount,1)
          write(ADM_log_fid,*) 'areq after isend',dbg_areq_save(1:acount,2)
          write(ADM_log_fid,*) 'areq after waitall',dbg_areq_save(1:acount,3)

          write(ADM_log_fid,*) 'ierr of mpi_waitall=',ierr
          write(ADM_log_fid,*) 'acount of mpi_waitall=',acount
          write(ADM_log_fid,*) 'stat of mpi_waitall=',stat(:,1:acount)

       endif
       dbg_tcount=dbg_tcount+1
    endif  !=================================== dbg end


    if  (opt_comm_barrier) then
       call mpi_barrier(ADM_comm_world,ierr)
       if (ierr.ne.0) then
          write(ADM_log_fid,*) 'mpi_barriert info start=='
          write(ADM_log_fid,*) 'ierr=',ierr
          write(ADM_log_fid,*) 'mpi_barrier info end=='
       endif
    endif


    !-----------------------------------------
    !
    do ro=1,romax(halo)
       !-----------------------------------------
       !  recvbuf -> var ( recieve in region )
       !-----------------------------------------
       do nr=1,nrmax(ro,halo)
          rs=recvinfo(SIZE_COMM,nr,ro,halo)
          rl=recvinfo(LRGNID_COMM,nr,ro,halo)
          rb=recvinfo(BASE_COMM,nr,ro,halo)*cmax
          do m=1,varmax
!cdir outerunroll=8
             do k=1,kmax
                do n=1,rs
                   var(recvlist(n,nr,ro,halo),k,rl,clist(m)) &
                        =recvbuf(n+(k-1)*rs+(m-1)*rs*kmax+rb,ro)
                enddo
             enddo
          enddo
       enddo
    enddo !loop ro
    t(10)=mpi_wtime()
    time_rbuf=time_rbuf+(t(10)-t(9))

    !-----------------------------------------
    !
    !-----------------------------------------
    !  copy data around singular point
    !-----------------------------------------
    do nc=1,ncmax_sgp(halo)
       cs=copyinfo_sgp(SIZE_COPY,nc,halo)
       cl=copyinfo_sgp(LRGNID_COPY,nc,halo)
       scl=copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)
       do m=1,varmax
!cdir outerunroll=8
          do k=1,kmax
             do n=1,cs
                var(recvlist_sgp(n,nc,halo),k,cl ,clist(m)) &
                     =var(sendlist_sgp(n,nc,halo),k,scl,clist(m))
             enddo
          enddo
       enddo

    enddo
    t(11)=mpi_wtime()
    time_copy_sgp=time_copy_sgp+(t(11)-t(10))
    !!
    !call mpi_barrier(ADM_comm_world,ierr)
    t(12)=mpi_wtime()
    time_bar2=time_bar2+(t(12)-t(11))
    time_total=time_total+(t(12)-t(0))
    !
  end subroutine COMM_data_transfer_nopl

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum( localsum, globalsum )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all,        &
       ADM_prc_me
    implicit none

    real(8), intent(in)  :: localsum
    real(8), intent(out) :: globalsum

    real(8) :: sendbuf(1)
    real(8) :: recvbuf(ADM_prc_all)

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localsum

       call MPI_Allgather( sendbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           recvbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           ADM_COMM_WORLD,       &
                           ierr                  )

       globalsum = sum( recvbuf(:) )
    else
       globalsum = localsum
    endif

    return
  end subroutine COMM_Stat_sum

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum_eachlayer( kall, localsum, globalsum )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all,        &
       ADM_prc_me
    implicit none

    integer, intent(in)  :: kall
    real(8), intent(in)  :: localsum (kall)
    real(8), intent(out) :: globalsum(kall)

    real(8) :: sendbuf(kall)
    integer :: displs (ADM_prc_all)
    integer :: counts (ADM_prc_all)
    real(8) :: recvbuf(kall,ADM_prc_all)

    integer :: ierr
    integer :: k, p
    !---------------------------------------------------------------------------

    do p = 1, ADM_prc_all
       displs(p) = (p-1) * kall
       counts(p) = kall
    enddo

    if ( COMM_pl ) then
       sendbuf(:) = localsum(:)

       call MPI_Allgatherv( sendbuf,              &
                            kall,                 &
                            MPI_DOUBLE_PRECISION, &
                            recvbuf,              &
                            counts,               &
                            displs,               &
                            MPI_DOUBLE_PRECISION, &
                            ADM_COMM_WORLD,       &
                            ierr                  )

       do k = 1, kall
          globalsum(k) = sum( recvbuf(k,:) )
       enddo
    else
       do k = 1, kall
          globalsum(k) = localsum(k)
       enddo
    endif

    return
  end subroutine COMM_Stat_sum_eachlayer

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_avg( localavg, globalavg )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all,        &
       ADM_prc_me
    implicit none

    real(8), intent(in)  :: localavg
    real(8), intent(out) :: globalavg

    real(8) :: sendbuf(1)
    real(8) :: recvbuf(ADM_prc_all)

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localavg

       call MPI_Allgather( sendbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           recvbuf,              &
                           1,                    &
                           MPI_DOUBLE_PRECISION, &
                           ADM_COMM_WORLD,       &
                           ierr                  )

       globalavg = sum( recvbuf(:) ) / real(ADM_prc_all,kind=8)
    else
       globalavg = localavg
    endif

    !write(ADM_LOG_FID,*) 'COMM_Stat_avg', sendbuf(1), recvbuf(:)

    return
  end subroutine COMM_Stat_avg

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_max( localmax, globalmax )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all,        &
       ADM_prc_me
    implicit none

    real(8), intent(in)  :: localmax
    real(8), intent(out) :: globalmax

    real(8) :: sendbuf(1)
    real(8) :: recvbuf(ADM_prc_all)

    integer :: ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmax

    call MPI_Allgather( sendbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        recvbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        ADM_COMM_WORLD,       &
                        ierr                  )

    globalmax = maxval( recvbuf(:) )

    !write(ADM_LOG_FID,*) 'COMM_Stat_max', sendbuf(1), recvbuf(:)

    return
  end subroutine COMM_Stat_max

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_min( localmin, globalmin )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all,        &
       ADM_prc_me
    implicit none

    real(8), intent(in)  :: localmin
    real(8), intent(out) :: globalmin

    real(8) :: sendbuf(1)
    real(8) :: recvbuf(ADM_prc_all)

    integer :: ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmin

    call MPI_Allgather( sendbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        recvbuf,              &
                        1,                    &
                        MPI_DOUBLE_PRECISION, &
                        ADM_COMM_WORLD,       &
                        ierr                  )

    globalmin = minval( recvbuf(:) )

    !write(ADM_LOG_FID,*) 'COMM_Stat_min', sendbuf(1), recvbuf(:)

    return
  end subroutine COMM_Stat_min

end module mod_comm
!-------------------------------------------------------------------------------
