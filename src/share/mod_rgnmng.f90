!-------------------------------------------------------------------------------
!>
!! Region management module
!!
!! @par Description
!!         This module is for region-process topology management
!!
!! @author NICAM developers
!!
!! @par History
!! @li      2004-02-17 (H.Tomita  ) [NEW]
!! @li      2007-10-22 (T.Mitsui  ) change value of rgn_nlim
!! @li      2010-06-07 (S.Iga     ) new grid systems (XTMS)
!! @li      2011-07-21 (T.Ohno    ) new grid systems (PERIODIC-1DMD)
!! @li      2011-07-21 (M.Hara    ) new grid systems (1DMD-ON-SPHERE)
!! @li      2013-01-16 (H.Yashiro ) NICAM-DC
!!
!<
module mod_rgnmng
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: RGNMNG_setup
  public :: RGNMNG_input
  public :: RGNMNG_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !====== Information for processes-region relationship ======

  integer, public, allocatable :: ADM_prc_rnum(:)  !------ Number of regions mangeged by each process
  integer, public, allocatable :: ADM_prc_tab(:,:)  !------ Table of regions managed by each process

  !> Table of edge link information
  !>
  !> ADM_rgn_etab( RID:DIR,    &
  !>               SW:SE,      &
  !>               ADM_rgn_all )
  !>
  integer, public, allocatable, save :: ADM_rgn_etab(:,:,:)

  !> Table of process ID from region ID
  !>
  !> ADM_rgn2prc(ADM_rgn_all)
  !>
  integer, public, allocatable, save :: ADM_rgn2prc(:)

  !------ Table of n-vertex-link(?) at the region vertex
  !<-----
  !<----- ADM_rgn_vnum( ADM_W:ADM_S, &
  !<-----               ADM_rgn_nmax )
  !<-----
  integer, public, allocatable, save :: ADM_rgn_vnum(:,:)

  !------ Table of vertex link information
  !<-----
  !<----- ADM_rgn_vtab( ADM_RID:ADM_DIR, &
  !<-----               ADM_W:ADM_S,     &
  !<-----               ADM_rgn_nmax,    &
  !<-----               ADM_vlink_nmax   )
  !<-----
  integer, public, allocatable, save :: ADM_rgn_vtab(:,:,:,:)

  !------ Table of vertex link information for poles
  !<-----
  !<----- ADM_rgn_vtab_pl( ADM_RID:ADM_DIR, &
  !<-----                  ADM_rgn_nmax_pl, &
  !<-----                  ADM_vlink_nmax   )
  !<-----
  integer, public, allocatable, save :: ADM_rgn_vtab_pl(:,:,:)

  !------ Region ID (reguler) of north/south pole management
  integer, public, save :: ADM_rgn_npl
  integer, public, save :: ADM_rgn_spl

  integer, public :: ADM_prc_npl        !< process ID which have the north pole
  integer, public :: ADM_prc_spl        !< process ID which have the south pole

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: RGNMNG_generate_ico

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: debug = .false. !< debug option

  character(len=IO_FILECHR), private :: RGNMNG_in_fname  = '' !< input  file name for region management file
  character(len=IO_FILECHR), private :: RGNMNG_out_fname = '' !< output file name for region management file

  integer, private :: RGNMNG_nprc !< # of processor for generate region management file

  !------ Maximum number of regions managed by 1 process.
  integer, private, parameter :: ADM_rgn_lim = 2560

  character(len=ADM_NSYS), private :: mapping_type = '' ! mapping method [add] C.Kodama 2011/12/14
                                                 ! = ''       : standard
                                                 ! = 'K-TERAI': terai mapping for K computer

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine RGNMNG_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_rlevel,    &
       ADM_rgn_all,   &
       ADM_prc_all,   &
       ADM_XTMS_K,    &
       ADM_MLCP_S
    implicit none

    namelist / PARAM_RGNMNG / &
        RGNMNG_in_fname,  &
        RGNMNG_out_fname, &
        RGNMNG_nprc,      &
        debug

    character(len=ADM_MAXFNAME), intent(in) :: fname

    integer, intent(in)  :: rlevel
    integer, intent(in)  :: rgn_nmax
    integer, intent(in)  :: prc_all
    integer, intent(in)  :: XTMS_K
    integer, intent(in)  :: XTMS_MLCP_S
    integer, intent(out) :: rgn_etab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,rgn_nmax)
    integer, intent(out) :: prc_rnum(prc_all)
    integer, intent(out) :: prc_tab (rgn_nlim,prc_all)
    integer, intent(out) :: rgn2prc (rgn_nmax)
    !---------------------------------------------------------------------------

    ! fill default
    RGNMNG_nprc = ADM_prc_all

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[rgnmng]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RGNMNG,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Use default.'
    elseif ( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_RGNMNG. STOP.'
       call ADM_proc_stop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_RGNMNG)

    allocate( RGN_nrgn(ADM_prc_all) )

    allocate( RGN_prc_tab(ADM_llim,ADM_prc_all) )

    allocate( RGN_rgn2prc(ADM_rgn_all) )

    allocate( RGN_edge_tab( RGN_RID:RGN_DIR, &
                            RGN_SW:RGN_SE,   &
                            ADM_rgn_all      ) )

    allocate( ADM_vert_nrgn( RGN_W:RGN_S,   &
                             ADM_rgn_all  ) )

    allocate( ADM_vert_tab( RGN_RID:RGN_DIR, &
                            RGN_W:RGN_S,     &
                            ADM_rgn_all,     &
                            ADM_vlink_nmax   ) )

    allocate( ADM_vert_tab_pl( RGN_RID:RGN_DIR, &
                               ADM_rgn_all_pl,  &
                               ADM_vlink_nmax   ) )


    if ( RGNMNG_in_fname /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** input file is not specified.'

       if ( topology_type == I_default ) then
          call RGNMNG_generate_ico( ADM_rlevel,   & ! [IN]
                                    ADM_rgn_all, & ! [IN]
                                    RGNMNG_nprc,  & ! [IN]
                                    RGNMNG_rgn_etab, & ! [OUT]
                                    RGNMNG_prc_rnum, & ! [OUT]
                                    RGNMNG_prc_tab,  & ! [OUT]
                                    RGNMNG_rgn2prc   ) ! [OUT]
       endif
    else
       call RGNMNG_input( fname,    & ! [IN]
                          rgn_etab, & ! [OUT]
                          prc_rnum, & ! [OUT]
                          prc_tab,  & ! [OUT]
                          rgn2prc   ) ! [OUT]
    endif


!    if ( ADM_debug ) then
!       do l = 1, ADM_lall
!          rgnid = ADM_prc_tab(l,ADM_prc_me)
!          if( IO_L ) write(IO_FID_LOG,*) ' --- Link information for region', rgnid
!
!          if( IO_L ) write(IO_FID_LOG,*) '     < edge link >   --- ( rgnid , edgid )'
!          do d = ADM_SW, ADM_SE
!             if( IO_L ) write(IO_FID_LOG,*) '     (',rgnid,',',d,') -> ',         &
!                                  '(', ADM_rgn_etab(ADM_RID,d,rgnid),   &
!                                  ',', ADM_rgn_etab(ADM_DIR,d,rgnid), ')'
!          enddo
!
!          if( IO_L ) write(IO_FID_LOG,*) '     < vertex link > --- ( rgnid , edgid )'
!          do d = ADM_W, ADM_S
!             if( IO_L ) write(IO_FID_LOG,*) '     (',rgnid,',',d,') : ', ADM_rgn_vnum(d,rgnid), 'point link'
!             do v = 1, ADM_rgn_vnum(d,rgnid)
!                if( IO_L ) write(IO_FID_LOG,*) '                -> ',                  &
!                                     '(', ADM_rgn_vtab(ADM_RID,d,rgnid,v),   &
!                                     ',', ADM_rgn_vtab(ADM_DIR,d,rgnid,v), ')'
!             enddo
!          enddo
!
!       enddo
!
!       if( IO_L ) write(IO_FID_LOG,*) ' --- Table of corresponding between region ID and process ID'
!       if( IO_L ) write(IO_FID_LOG,*) '    region ID :  process ID'
!       do l = 1, ADM_lall
!          if( IO_L ) write(IO_FID_LOG,'(I13,I14)') l, ADM_rgn2prc(l)
!       enddo
!    endif

    return
  end subroutine RGNMNG_setup

  !-----------------------------------------------------------------------------
  subroutine RGNMNG_input( &
      fname,    &
      rgn_etab, &
      prc_rnum, &
      prc_tab,  &
      rgn2prc   )
    use mod_misc, only: &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_prc_stop, &
       rgn_nlim, &
       ADM_RID, &
       ADM_DIR, &
       ADM_SW, &
       ADM_NW, &
       ADM_NE, &
       ADM_SE, &
       rgn_nmax, &
       prc_all
    implicit none

    character(len=ADM_MAXFNAME), intent(in) :: fname

    integer, intent(out) :: rgn_etab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,rgn_nmax)
    integer, intent(out) :: prc_rnum(prc_all)
    integer, intent(out) :: prc_tab (rgn_nlim,prc_all)
    integer, intent(out) :: rgn2prc (rgn_nmax)

    integer :: num_of_rgn !--- number of region

    namelist / rgn_info / &
         num_of_rgn

    integer :: rgnid                    !--- region ID
    integer :: sw(ADM_RID:ADM_DIR) = -1 !--- south-west region info
    integer :: nw(ADM_RID:ADM_DIR) = -1 !--- nouth-west region info
    integer :: ne(ADM_RID:ADM_DIR) = -1 !--- nouth-east region info
    integer :: se(ADM_RID:ADM_DIR) = -1 !--- south-east region info

    namelist / rgn_link_info / &
         rgnid, &
         sw,    &
         nw,    &
         ne,    &
         se                     

    integer :: num_of_proc !--- number of run-processes

    namelist /proc_info/ &
         num_of_proc

    integer :: peid                         !--- process ID
    integer :: num_of_mng                   !--- number of regions be managed
    integer :: mng_rgnid(rgn_nlim) = -1 !--- managed region ID

    namelist /rgn_mng_info/ &
         peid,       &
         num_of_mng, &
         mng_rgnid            

    integer :: fid, ierr
    integer :: l, p
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** input mnginfo:',trim(fname)

    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       ! ERROR if filename are not defined
       if ( ierr /= 0 ) then
          write(ADM_LOG_FID,*) 'xxx mnginfo file is not found! STOP. ', trim(fname)
          call ADM_proc_stop
       endif

       read(fid,nml=rgn_info)
       if ( num_of_rgn /= rgn_nmax ) then
          write(ADM_LOG_FID,*) 'xxx No match for region number! STOP.'
          write(ADM_LOG_FID,*) 'xxx rgn_nmax= ',rgn_nmax,' num_of_rgn=',num_of_rgn
          call ADM_proc_stop
       endif

       do l = 1, rgn_nmax
          read(fid,nml=rgn_link_info)

          rgn_etab(ADM_RID:ADM_DIR,ADM_SW,rgnid) = sw(ADM_RID:ADM_DIR)
          rgn_etab(ADM_RID:ADM_DIR,ADM_NW,rgnid) = nw(ADM_RID:ADM_DIR)
          rgn_etab(ADM_RID:ADM_DIR,ADM_NE,rgnid) = ne(ADM_RID:ADM_DIR)
          rgn_etab(ADM_RID:ADM_DIR,ADM_SE,rgnid) = se(ADM_RID:ADM_DIR)
       enddo

       read(fid,nml=proc_info)
       if ( num_of_proc /= prc_all ) then
          write(ADM_LOG_FID,*) ' xxx No match for  process number! STOP.'
          write(ADM_LOG_FID,*) ' xxx prc_all= ',prc_all,' num_of_proc=',num_of_proc
          call ADM_proc_stop
       endif

       prc_tab(:,:) = -1 ! [Fix] 11/06/30  T.Seiki, fill undefined value
       do p = 1, prc_all
          read(fid,nml=rgn_mng_info)

          prc_rnum(peid)   = num_of_mng
          prc_tab (:,peid) = mng_rgnid(:)
          do l = 1, prc_rnum(peid)
             rgn2prc(prc_tab(l,peid)) = peid
          enddo
       enddo

    close(fid)

    return
  end subroutine RGNMNG_input

  !-----------------------------------------------------------------------------
  subroutine RGNMNG_output( &
      fname,    &
      rgn_etab, &
      prc_rnum, &
      prc_tab   )
    use mod_misc, only: &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_prc_stop, &
       rgn_nlim, &
       ADM_RID, &
       ADM_DIR, &
       ADM_SW, &
       ADM_NW, &
       ADM_NE, &
       ADM_SE, &
       rgn_nmax, &
       prc_all
    implicit none

    character(len=ADM_MAXFNAME), intent(in) :: fname

    integer, intent(in) :: rgn_etab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,rgn_nmax)
    integer, intent(in) :: prc_rnum(prc_all)
    integer, intent(in) :: prc_tab (rgn_nlim,prc_all)

    integer :: num_of_rgn !--- number of region

    namelist / rgn_info / &
         num_of_rgn

    integer :: rgnid                    !--- region ID
    integer :: sw(ADM_RID:ADM_DIR) = -1 !--- south-west region info
    integer :: nw(ADM_RID:ADM_DIR) = -1 !--- nouth-west region info
    integer :: ne(ADM_RID:ADM_DIR) = -1 !--- nouth-east region info
    integer :: se(ADM_RID:ADM_DIR) = -1 !--- south-east region info

    namelist / rgn_link_info / &
         rgnid, &
         sw,    &
         nw,    &
         ne,    &
         se                     

    integer :: num_of_proc !--- number of run-processes

    namelist /proc_info/ &
         num_of_proc

    integer :: peid                         !--- process ID
    integer :: num_of_mng                   !--- number of regions be managed
    integer :: mng_rgnid(rgn_nlim) = -1 !--- managed region ID

    namelist /rgn_mng_info/ &
         peid,       &
         num_of_mng, &
         mng_rgnid            

    integer :: fid, ierr
    integer :: l, p
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '*** output mnginfo:',trim(fname)

    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted'  )

       num_of_rgn = rgn_nmax
       write(fid,nml=rgn_info)

       do l = 1, rgn_nmax
          rgnid = l
          sw(ADM_RID:ADM_DIR) = rgn_etab(ADM_RID:ADM_DIR,ADM_SW,rgnid)
          nw(ADM_RID:ADM_DIR) = rgn_etab(ADM_RID:ADM_DIR,ADM_NW,rgnid)
          ne(ADM_RID:ADM_DIR) = rgn_etab(ADM_RID:ADM_DIR,ADM_NE,rgnid)
          se(ADM_RID:ADM_DIR) = rgn_etab(ADM_RID:ADM_DIR,ADM_SE,rgnid)

          write(fid,nml=rgn_link_info)
       enddo

       num_of_proc = prc_all
       read(fid,nml=proc_info)

       do p = 1, prc_all
          peid = p
          num_of_mng   = prc_rnum(peid)
          mng_rgnid(:) = prc_tab (:,peid)

          write(fid,nml=rgn_mng_info)
       enddo

    close(fid)

    return
  end subroutine RGNMNG_output

  !-----------------------------------------------------------------------------
  Subroutine RGNMNG_generate_ico( &
      rlevel,   &
      rgn_nmax, &
      prc_all,  &
      rgn_etab, &
      prc_rnum, &
      prc_tab,  &
      rgn2prc   )
    implicit none

    integer, intent(in)  :: rlevel
    integer, intent(in)  :: rgn_nmax
    integer, intent(in)  :: prc_all
    integer, intent(out) :: rgn_etab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,rgn_nmax)
    integer, intent(out) :: prc_rnum(prc_all)
    integer, intent(out) :: prc_tab (rgn_nlim,prc_all)
    integer, intent(out) :: rgn2prc (rgn_nmax)

    integer :: lall_1d, lall

    integer, parameter :: nmax_dmd = 10
    integer :: dmd_data(ADM_SW:ADM_SE,nmax_dmd)

    integer :: d_nb, i_nb, j_nb, rgnid_nb, direction
    integer :: d, i, j, rgnid
    integer :: l, p, tmp
    integer :: fid
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** generate region management table'
    write(ADM_LOG_FID,*) '*** Topology: default icosahedral'

    dmd_data(:, 1) = (/  6, 5, 2,10 /)
    dmd_data(:, 2) = (/ 10, 1, 3, 9 /)
    dmd_data(:, 3) = (/  9, 2, 4, 8 /)
    dmd_data(:, 4) = (/  8, 3, 5, 7 /)
    dmd_data(:, 5) = (/  7, 4, 1, 6 /)
    dmd_data(:, 6) = (/  7, 5, 1,10 /)
    dmd_data(:, 7) = (/  8, 4, 5, 6 /)
    dmd_data(:, 8) = (/  9, 3, 4, 7 /)
    dmd_data(:, 9) = (/ 10, 2, 3, 8 /)
    dmd_data(:,10) = (/  6, 1, 2, 9 /)

    lall_1d  = 2**rlevel
    lall     = lall_1d * lall_1d

    !--- make region link table
    do d = 1, nmax_dmd
    do i = 1, lall_1d
    do j = 1, lall_1d
       rgnid = (d-1)*lall + (j-1)*lall_1d + i

       !--- ADM_SW
       if ( j == 1 ) then
          if ( d <= 5 ) then
             i_nb = i
             j_nb = lall_1d
             d_nb = dmd_data(ADM_SW,d)
             direction = ADM_NE
          else
             i_nb = lall_1d
             j_nb = lall_1d+1-i
             d_nb = dmd_data(ADM_SW,d)
             direction = ADM_SE
          endif
       else
          i_nb = i
          j_nb = j-1
          d_nb = d
          direction = ADM_NE
       endif
       rgnid_nb = (d_nb-1)*lall + (j_nb-1)*lall_1d + i_nb

       rgn_etab(ADM_RID,ADM_SW,rgnid) = rgnid_nb
       rgn_etab(ADM_DIR,ADM_SW,rgnid) = direction

       !--- ADM_NW
       if ( i == 1 ) then
          if ( d <= 5 ) then
             i_nb = lall_1d+1-j
             j_nb = lall_1d
             d_nb = dmd_data(ADM_NW,d)
             direction = ADM_NE
          else
             i_nb = lall_1d
             j_nb = j
             d_nb = dmd_data(ADM_NW,d)
             direction = ADM_SE
          endif
       else
          i_nb = i-1
          j_nb = j
          d_nb = d
          direction = ADM_SE
       endif
       rgnid_nb = (d_nb-1)*lall + (j_nb-1)*lall_1d + i_nb

       rgn_etab(ADM_RID,ADM_NW,rgnid) = rgnid_nb
       rgn_etab(ADM_DIR,ADM_NW,rgnid) = direction

       !--- ADM_NE
       if ( j == lall_1d ) then
          if ( d <= 5 ) then
             i_nb = 1
             j_nb = lall_1d+1-i
             d_nb = dmd_data(ADM_NE,d)
             direction = ADM_NW
          else
             i_nb = i
             j_nb = 1
             d_nb = dmd_data(ADM_NE,d)
             direction = ADM_SW
          endif
       else
          i_nb = i
          j_nb = j+1
          d_nb = d
          direction = ADM_SW
       endif
       rgnid_nb = (d_nb-1)*lall + (j_nb-1)*lall_1d + i_nb

       rgn_etab(ADM_RID,ADM_NE,rgnid) = rgnid_nb
       rgn_etab(ADM_DIR,ADM_NE,rgnid) = direction

       !--- ADM_SE
       if ( i == lall_1d ) then
          if ( d <= 5 ) then
             i_nb = 1
             j_nb = j
             d_nb = dmd_data(ADM_SE,d)
             direction = ADM_NW
          else
             i_nb = lall_1d+1-j
             j_nb = 1
             d_nb = dmd_data(ADM_SE,d)
             direction = ADM_SW
          endif
       else
          i_nb = i+1
          j_nb = j
          d_nb = d
          direction = ADM_NW
       endif
       rgnid_nb = (d_nb-1)*lall + (j_nb-1)*lall_1d + i_nb

       rgn_etab(ADM_RID,ADM_SE,rgnid) = rgnid_nb
       rgn_etab(ADM_DIR,ADM_SE,rgnid) = direction

    enddo
    enddo
    enddo

    if ( mod(rgn_nmax,prc_all) /= 0 ) then
       write(*,*) 'invalid number of process!', rgn_nmax, prc_all
       stop
    else
       prc_rnum(:) = rgn_nmax / prc_all
    endif

    !--- make region-pe relationship
    prc_tab(:,:) = -1
    if ( RGNMNG_PRCMAP_TYPE == 'K-TERAI' ) then ! [Add] C.Kodama 2011/12/14
       do p = 1, prc_all
          if ( prc_rnum(p) > 2 ) then
             write(*,*) 'RGNMNG_PRCMAP_TYPE = K-TERAI'
             write(*,*) 'More than two regions is not allowed. stop'
             stop
          endif

          do l = 1, prc_rnum(p)
             ! if 1reg-1prc, rgnid = p
             ! if 2reg-1prc, rgnid = 2p or 2p-1
             rgnid = (p-1)*prc_rnum(p) + l
          enddo

          if ( mod(rgnid-1,2*lall) < lall ) then
             tmp =             ( rgnid-1 - mod(rgnid-1,lall) ) / (2*lall)
          else
             tmp = 10 - ( lall + rgnid-1 - mod(rgnid-1,lall) ) / (2*lall)
          endif

          prc_tab(l,p) = tmp * lall + mod(rgnid-1,lall) + 1
       enddo
    else
       do p = 1, prc_all
       do l = 1, prc_rnum(p)
          prc_tab(l,p) = (p-1)*prc_rnum(p) + l
       enddo
       enddo
    endif

    do p = 1, prc_all
       write(ADM_LOG_FID,*) 'peid=', p, 'regid=', prc_tab(l,p)

       do l = 1, prc_rnum(p)
          rgn2prc(prc_tab(l,p)) = p
       enddo
    enddo

    return
  end Subroutine RGNMNG_generate_ico

end module mod_rgnmng
!-------------------------------------------------------------------------------
