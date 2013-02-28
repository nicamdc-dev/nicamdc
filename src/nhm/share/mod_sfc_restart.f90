!-------------------------------------------------------------------------------
! 
!+  restart sfc I/O module
!
!-------------------------------------------------------------------------------
module mod_sfc_restart
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the sfc restart file I/O.
  !       
  ! 
  !++ Current Corresponding Author : T.Mitsui
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      06-04-18   Add this module
  !                08-03-10   T.Mitsui: add output of intermediate restart file
  !                11-09-03   H.Yashiro : New I/O
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only : &
    ADM_MAXFNAME
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: sfc_restart_setup
  public :: sfc_restart_output_all
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(LEN=ADM_MAXFNAME), private, save :: input_io_mode  = 'LEGACY' ! [add] H.Yashiro 20110819
  character(LEN=ADM_MAXFNAME), private, save :: output_io_mode = 'LEGACY' ! [add] H.Yashiro 20110819
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine sfc_restart_setup ( &
       ctime                   )
    !
    use mod_adm, only :     &
         ADM_proc_stop, &
         ADM_MAXFNAME,      &
         ADM_LOG_FID,       &
         ADM_GALL_PL,       &
         ADM_LALL_PL,       &
         ADM_gall,          &
         ADM_lall,          &
         ADM_KNONE,         &
         ADM_VMISS,         &
         adm_prc_me,        &
         adm_prc_pl
    use mod_runconf, only : &
         NRDIR,             &
         NRBND
    use mod_fio, only : & ! [add] H.Yashiro 20110826
         FIO_input
    use mod_sfcvar, only : &
         sfcvar_set2,      &
         I_ALBEDO_SFC,     &
         sfcvar_comm
    !
    implicit none
    !
    real(8), intent(in) :: ctime
    !--- albedo
    real(8) :: albedo_sfc   (ADM_gall   ,ADM_KNONE,ADM_lall   ,NRDIR,NRBND)
    real(8) :: albedo_sfc_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,NRDIR,NRBND)
    real(8), allocatable :: tmp(:,:,:)
    real(8), allocatable :: tmp_pl(:,:,:)
    !
    real(8) :: ALB_INIT = 0.2d0
    namelist /nm_sfc_init/ ALB_INIT
    !
    integer :: ij,k,k1,k2,l
    integer :: knum
    integer :: level
    !
    allocate(tmp(ADM_gall,NRDIR*NRBND,ADM_lall) )
    allocate(tmp_pl(ADM_gall_pl,NRDIR*NRBND,ADM_lall_pl) )    
    tmp(:,:,:) = ALB_INIT
    tmp_pl(:,:,:) = ALB_INIT
    !
    knum=NRDIR*NRBND
    call sfc_restart_input( &
         'ALB_SFC',          & !--- in
         tmp(:,:,:),         & !--- inout
         tmp_pl(:,:,:),      & !--- inout
         knum,               & !--- in  !-- default layer num for the data
         level              )  !--- out !-- modified layer num of the data
    !
    if(level /= 0) then
       write(ADM_LOG_FID,*) ' Msg : Sub[sfc_restart_setup]/Mod[sfc_restart]'
       write(ADM_LOG_FID,*) ' Msg : layer num differ from default: ALB_SFC ',knum,' to ',level
    end if
    do l=1, adm_lall
       k=0
       do k2=1,NRBND
          do k1=1,NRDIR
             k = k + 1
             if((level /= 0).and. k>level) then
                k = level
             end if
             do ij=1,ADM_gall
                albedo_sfc(ij,ADM_KNONE,l,k1,k2)=tmp(ij,k,l)
             end do
          end do
       end do
    end do
    if(adm_prc_me == adm_prc_pl) then
       do l=1, adm_lall_pl
          k=0
          do k2=1,NRBND
             do k1=1,NRDIR
                k = k + 1
             if((level /= 0).and. k>level) then
                k = level
             end if
                do ij=1,ADM_gall_pl
                   albedo_sfc_pl(ij,ADM_KNONE,l,k1,k2)=tmp_pl(ij,k,l)
                end do
             end do
          end do
       end do
    else
       albedo_sfc_pl(:,:,:,:,:)= ADM_VMISS
    end if
    deallocate(tmp)
    deallocate(tmp_pl)

!!! debug
!!$    if(adm_prc_me == 1) then
!!$       write(*,*) ' INPUT albedo (100,1,2,:,:) :' , albedo_sfc(100,1,2,:,:) 
!!$    end if
!!!
    !
    call sfcvar_set2( albedo_sfc, albedo_sfc_pl, I_ALBEDO_SFC, NRDIR, NRBND )
    call sfcvar_comm( comm_type=2 )
    !
    write(ADM_LOG_FID,*) ' Msg : Sub[sfc_restart_setup]/Mod[sfc_restart]'
    write(ADM_LOG_FID,*) &
         '   *** input sfc restart data at the time : ', ctime

    return
  end subroutine sfc_restart_setup
  !-----------------------------------------------------------------------------
  subroutine  sfc_restart_output_all( &
      ctime, cdate,     & !--- IN
      albsfc_out_atland ) !--- IN [add] H.Yashiro 20120512
    use mod_adm, only : &
         ADM_LOG_FID,   &
         ADM_gall,      &
         ADM_gall_pl,   &
         ADM_KNONE,     &
         ADM_lall,      &
         ADM_lall_pl,   &
         ADM_VMISS,     &
         adm_prc_me,    &
         adm_prc_pl
    use mod_runconf, only : &
         NRDIR,             &
         NRBND
    use mod_sfcvar, only : &
         sfcvar_get2,      &
         sfcvar_comm,      &
         I_ALBEDO_SFC
    !
    implicit none
    real(8), intent(in) :: ctime
    character(len=14), intent(in) :: cdate ! 08/03/10 [Add] T.Mitsui
    logical, intent(in) :: albsfc_out_atland    ! [add] H.Yashiro 20120512
    !--- albedo
    real(8) :: albedo_sfc   (ADM_gall   ,ADM_KNONE,ADM_lall   ,NRDIR,NRBND)
    real(8) :: albedo_sfc_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,NRDIR,NRBND)
    real(8), allocatable :: tmp(:,:,:)
    real(8), allocatable :: tmp_pl(:,:,:)
    !
    integer :: ij,k,k1,k2,l
    !---------------------------------------------------------------------------

    ! [add] H.Yashiro 20120512
    ! Skip if albedo_sfc is treated in land_restart
    if( albsfc_out_atland ) return

    call sfcvar_comm( comm_type=2 )
    call sfcvar_get2( albedo_sfc, albedo_sfc_pl, I_ALBEDO_SFC, NRDIR, NRBND )
    !
    allocate(tmp   (ADM_gall,   NRDIR*NRBND,ADM_lall   ) )
    allocate(tmp_pl(ADM_gall_pl,NRDIR*NRBND,ADM_lall_pl) )

    do l = 1, adm_lall
       k = 0
       do k2 = 1, NRBND
       do k1 = 1, NRDIR
          k = k + 1
          do ij = 1, ADM_gall
             tmp(ij,k,l) = albedo_sfc(ij,ADM_KNONE,l,k1,k2)
          enddo
       enddo
       enddo
    enddo

    if(adm_prc_me == adm_prc_pl) then
       do l = 1, adm_lall_pl
          k = 0
          do k2 = 1, NRBND
          do k1 = 1, NRDIR
             k = k + 1
             do ij = 1, ADM_gall_pl
                tmp_pl(ij,k,l) = albedo_sfc_pl(ij,ADM_KNONE,l,k1,k2)
             enddo
          enddo
          enddo
       enddo
    else
       tmp_pl(:,:,:)=ADM_VMISS
    end if

    !--- output
    call sfc_restart_output( &
         'ALB_SFC',           & !--- in
         ctime,               & !--- in [add] H.Yashiro 20110819
         cdate,               & !--- in 08/03/10 [Add] T.Mitsui
         tmp(:,:,:),          & !--- in
         tmp_pl(:,:,:),       & !--- in
         NRDIR*NRBND        )   !--- in

    deallocate(tmp )
    deallocate(tmp_pl )
    !
!!! debug
!!$    if(adm_prc_me == 1) then
!!$       write(*,*) ' OUTPUT albedo (100,1,2,:,:) :' , albedo_sfc(100,1,2,:,:) 
!!$    end if
!!!
    write(ADM_LOG_FID,*) ' Msg : Sub[sfc_restart_output_all]/Mod[sfc_restart]'
    write(ADM_LOG_FID,*) &
         '   *** output sfc restart data at the time : ', ctime
    return
  end subroutine sfc_restart_output_all
  !-----------------------------------------------------------------------------
  subroutine sfc_restart_input( &
       DNAME,                  & !--- in
       gdata,                  & !--- inout
       gdata_pl,               & !--- inout
       knum,                   & !--- in
       klev                    ) !--- out
    use mod_misc, only :        &
         MISC_make_idstr,       &
         MISC_get_available_fid,&
         MISC_msg_nmerror
    use mod_adm, only :         &
         ADM_proc_stop, &
         ADM_MAXFNAME,          &
         ADM_LOG_FID,           &
         ADM_CTL_FID,           &
         ADM_GALL_PL,           &
         ADM_LALL_PL,           &
         ADM_prc_me,            &
         ADM_prc_pl,            &
         ADM_gall,              &
         ADM_lall,              &
         ADM_prc_tab
    use mod_gtl, only :         &
         GTL_input_var2,        &
         GTL_input_var2_da
    use mod_comm, only :     &
         comm_var
    use mod_fio, only : & ! [add] H.Yashiro 20110826
         FIO_input
    !
    implicit none
    !
    character(LEN=*), intent(in) :: DNAME     ! data name
    integer, intent(in) :: knum
    real(8), intent(inout) :: gdata(ADM_gall,knum,ADM_lall)          ! data
    real(8), intent(inout) :: gdata_pl(ADM_gall_pl,knum,ADM_lall_pl) ! data (pole)
    integer, optional, intent(out) :: klev                           ! data level(if unchanged level=0)
    !
    character(ADM_MAXFNAME) :: dataname = ''
    character(ADM_MAXFNAME) :: filename = ''
    !
    character(ADM_MAXFNAME) :: fname
    integer :: rgnid
    integer :: fid
    !
    logical :: direct_access = .false.
    !
    logical :: otag
    integer :: ierr
    integer :: l
    !
    integer :: level = 0
    integer :: knum_in
    !
    real(8), allocatable :: tmp(:,:,:,:)
    real(8), allocatable :: tmp_pl(:,:,:,:)
    !
    namelist /nm_sfc_restart_input/ &
         dataname, &
         filename, &
         level,    &
         direct_access, &
         input_io_mode   ! [add] H.Yashiro 20110826
    !
    rewind(ADM_CTL_FID)
    otag = .false.
    findtag: do
       read (ADM_CTL_FID, nml= nm_sfc_restart_input, iostat= ierr, end= 1001)
       if ( DNAME == dataname ) then
          write(ADM_LOG_FID,nml= nm_sfc_restart_input)
          otag = .true.
          exit findtag
       end if
    end do findtag
1001 continue
    !
    call MISC_msg_nmerror( &
         ierr,             & !--- in
         ADM_LOG_FID,      & !--- in
         'nm_sfc_restart_input',   & !--- in
         'sfc_restart_input', & !--- in
         'sfc_restart'    & !--- in
         )
    !
    knum_in = knum
    if (present(klev))then
       klev    = level
       if(level /= 0) then
          knum_in = level
       end if
    end if
    !
    if ( .not. otag ) then
       write(ADM_LOG_FID,*) &
            '###FILEREAD_DATA_SURFACE: warning no tag for:', DNAME
       ! not set gdata
       return
    end if
    !
    write(ADM_LOG_FID,*) &
         'Msg : Sub[sfc_restart_input]/Mod[sfc_restart]'
    write(ADM_LOG_FID,*) '### data=', DNAME, ': file name=', trim(filename)

    ! -> [add] H.Yashiro 20110826
    if ( input_io_mode == 'ADVANCED' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,input : ',trim(input_io_mode)
       direct_access = .true.
    elseif( input_io_mode == 'LEGACY' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,input : ',trim(input_io_mode)
    else
       write(ADM_LOG_FID,*) 'xxx Invalid input_io_mode!',trim(input_io_mode)
       call ADM_proc_stop
    endif
    ! <- [add] H.Yashiro 20110826

    if(direct_access) then
       allocate(tmp   (ADM_gall,   knum_in,ADM_lall,   1))
       allocate(tmp_pl(ADM_gall_pl,knum_in,ADM_lall_pl,1))

       ! -> [add] H.Yashiro 20110826
       if ( input_io_mode == 'ADVANCED' ) then

          call FIO_input( tmp(:,:,:,1),filename,'albedo_sfc','GSALB',1,1,1 )

       elseif( input_io_mode == 'LEGACY' ) then
       ! <- [add] H.Yashiro 20110826

          call GTL_input_var2_da(   &
               trim(filename),      &
               tmp(:,:,:,1), 1, knum_in,     &
               recnum=1, input_size=8 )

       ! -> [add] H.Yashiro 20110826
       endif
       ! <- [add] H.Yashiro 20110826

       call comm_var(   &
            tmp,tmp_pl, &
            knum_in,    &
            1,          &
            comm_type=2,&
            NSval_fix=.true.&
            )
       !
       gdata(:,1:knum_in,:) = tmp(:,1:knum_in,:,1)
       If(ADM_prc_me==ADM_prc_pl) Then
          gdata_pl(:,1:knum_in,:) = tmp_pl(:,1:knum_in,:,1)
       end If
       !
       deallocate(tmp)
       deallocate(tmp_pl)
    else
       !--- input
       Do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          Call MISC_make_idstr(fname,Trim(filename),'rgn',rgnid)
          fid = MISC_get_available_fid()
          Open(fid,file=Trim(fname),form='unformatted',status='old')
          read(fid) gdata(:,1:knum_in,l)
          Close(fid)
          !
       End Do
       !
       If(ADM_prc_me==ADM_prc_pl) Then
          fname=Trim(filename)//'.pl'
          fid = MISC_get_available_fid()
          Open(fid,file=Trim(fname),form='unformatted',status='old')
          read(fid) gdata_pl(:,1:knum_in,:)
          Close(fid)
          !
       end If
    end if
    !
    return
  end subroutine sfc_restart_input
  !-----------------------------------------------------------------------------
  subroutine sfc_restart_output( &
       DNAME,                  & !--- in
       ctime,                  & !--- in [add] H.Yashiro 20110819
       cdate,                  & !--- in 08/03/10 [Add] T.Mitsui
       gdata,                  & !--- in
       gdata_pl,               & !--- in
       knum                  )   !--- in
    
    use mod_misc, only :        &
         MISC_make_idstr,       &
         MISC_get_available_fid,&
         MISC_msg_nmerror
    use mod_adm, only :         &
         ADM_proc_stop,         &
         ADM_MAXFNAME,          &
         ADM_LOG_FID,           &
         ADM_CTL_FID,           &
         ADM_GALL_PL,           &
         ADM_LALL_PL,           &
         ADM_prc_me,            &
         ADM_prc_pl,            &
         ADM_gall,              &
         ADM_lall,              &
         ADM_prc_tab
    use mod_gtl, only :         &
         GTL_output_var2,       &
         GTL_output_var2_da  
    use mod_fio, only : & ! [add] H.Yashiro 20110819
         FIO_output, &
         FIO_HMID,   &
         FIO_REAL8
    implicit none
    !
    character(LEN=*), intent(in) :: DNAME     ! data name
    real(8),          intent(in) :: ctime     ! time [add] H.Yashiro 20110819
    character(len=14),intent(in) :: cdate     ! 08/03/10 [Add] T.Mitsui
    integer, intent(in) :: knum
    real(8), intent(in) :: gdata(ADM_gall,knum,ADM_lall)          ! data
    real(8), intent(in) :: gdata_pl(ADM_gall_pl,knum,ADM_lall_pl) ! data (pole)
    !    
    character(ADM_MAXFNAME) :: dataname = ''
    character(ADM_MAXFNAME) :: filename = ''
    !
    character(ADM_MAXFNAME) :: fname
    integer :: rgnid
    integer :: fid
    !   
    logical :: otag
    integer :: ierr
    integer :: l
    !
    logical :: direct_access = .false.
    !
    namelist /nm_sfc_restart_output/ &
         dataname,      &
         filename,      &
         direct_access, &
         output_io_mode   ! [add] H.Yashiro 20110826

    ! -> [add] H.Yashiro 20110819
    character(LEN=FIO_HMID) :: basename
    character(LEN=FIO_HMID) :: desc = 'INITIAL/RESTART DATA of ALBEDO'
    ! <- [add] H.Yashiro 20110819
    !---------------------------------------------------------------------------

    rewind(ADM_CTL_FID)
    otag = .false.
    findtag: do
       read (ADM_CTL_FID, nml=nm_sfc_restart_output, iostat=ierr, end=1001)
       if ( DNAME == trim(dataname) ) then
          write(ADM_LOG_FID, nml=nm_sfc_restart_output)
          otag = .true.
          exit findtag
       end if
    end do findtag
1001 continue
    !
    call MISC_msg_nmerror( &
         ierr,             & !--- in
         ADM_LOG_FID,      & !--- in
         'nm_sfc_restart_output',   & !--- in
         'sfc_restart_output', & !--- in
         'sfc_restart'    & !--- in
         )
    !
    if ( .not. otag ) then
       write(ADM_LOG_FID,*) &
            'Msg : Sub[sfc_restart_output]/Mod[sfc_restart]'
       write(ADM_LOG_FID,*) &
            '### warning no tag for:', DNAME
       write(ADM_LOG_FID,*) &
            '### use DNAME for restart filename'
       filename = 'restart_'//trim(DNAME)
    end if

    ! -> [add] H.Yashiro 20110819
    if ( output_io_mode == 'ADVANCED' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,output: ',trim(output_io_mode)
    elseif( output_io_mode == 'LEGACY' ) then
       write(ADM_LOG_FID,*) '*** io_mode for restart,output: ',trim(output_io_mode)
    else
       write(ADM_LOG_FID,*) 'xxx Invalid output_io_mode!',trim(output_io_mode)
       call ADM_proc_stop
    endif
    ! <- [add] H.Yashiro 20110819

    ! 08/03/10 [Mod] T.Mitsui
!!$ write(ADM_LOG_FID,*) '### data=', DNAME, ': file name=', trim(filename)
    write(ADM_LOG_FID,*) '### data=', DNAME, ': file name=', trim(filename)//trim(cdate)

    ! -> [add] H.Yashiro 20110819
    if ( output_io_mode == 'ADVANCED' ) then
       basename=trim(filename)//trim(cdate)

       call FIO_output( gdata(:,:,:), basename, desc, '',            &
                       'albedo_sfc', 'Surface Albedo', '', '0-1',    &
                        FIO_REAL8, 'GSALB', 1, knum, 1, ctime, ctime )

    elseif( output_io_mode == 'LEGACY' ) then
    ! <- [add] H.Yashiro 20110819

    if(direct_access) then
       call GTL_output_var2_da(&
            ! 08/03/10 [Mod] T.Mitsui
!!$         trim(filename),    &
            trim(filename)//trim(cdate),    &
            gdata(:,:,:),      &
            1,knum,            &
            recnum=1,          &
            output_size=8      &
            )
    else
       !--- output
       Do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          ! 08/03/10 [Mod] T.Mitsui
!!$       Call MISC_make_idstr(fname,Trim(filename),'rgn',rgnid)
          Call MISC_make_idstr(fname,Trim(filename)//trim(cdate),'rgn',rgnid)
          fid = MISC_get_available_fid()
          Open(fid,file=Trim(fname),form='unformatted',status='replace')
          Write(fid) gdata(:,:,l)
          Close(fid)
       End Do
       If(ADM_prc_me==ADM_prc_pl) Then
          ! 08/03/10 [Mod] T.Mitsui
!!$       fname=Trim(filename)//'.pl'
          fname=Trim(filename)//trim(cdate)//'.pl'
          fid = MISC_get_available_fid()
          Open(fid,file=Trim(fname),form='unformatted',status='replace')
          Write(fid) gdata_pl(:,:,:)
          Close(fid)
       End If
       !
    end if

    write(ADM_LOG_FID,*) &
         'Msg : Sub[sfc_restart_output]/Mod[sfc_sfc_restart]', trim(DNAME)

    ! -> [add] H.Yashiro 20110819
    endif  !--- io_mode
    ! <- [add] H.Yashiro 20110819

    return
  end subroutine sfc_restart_output
  !-----------------------------------------------------------------------------
end module mod_sfc_restart
