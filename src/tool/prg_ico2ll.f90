!-------------------------------------------------------------------------------
!
!+  program ico2ll
!
!-------------------------------------------------------------------------------
program ico2ll
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This program converts the icosahedral grid data
  !       to latitude-longitude grid.
  !
  !
  !++ Current Corresponding Author : M.Satoh, S.Iga, Y.Niwa, H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      03-02-17   Imported from igdc-4.34
  !                05-12-02   M.Satoh: modified according to mod_history
  !                05-12-21   M.Satoh: introduce input_size option
  !                                    to read double precision data
  !                05-12-26   S.Iga  : add options to output separately in time
  !                06-01-06   S.Iga  : Bug fix for separate_time = 'yes'
  !                06-03-22   S.Iga  : change 'JLY' to 'JUL'
  !                06-05-08   S.Iga  : Bug fix. (Before the modification,
  !                                    option 'no' and 'template' didn't
  !                                    work on the gbnode ).
  !                06-05-19   Y.Niwa : bug fix
  !                06-08-04   H.Tomita : add 'header_strings'
  !                06-08-21   H.Tomita : del duplicated sentense: latlon_data
  !                07-03-29   T.Mitsui : saving memory latlon_data_swap(ij,k=1)
  !                07-10-24   W.Yanase: debug in the calculation of 8-byte integer
  !                08-10-08   W.Yanase: use %h2 instead of %f2 in GrADS template
  !                09-02-04   M.Satoh:  add option access_icodata
  !                09-03-16   H.Taniguchi: add option access_icodata
  !                                        for sequential output of icodata
  !                09-09-28   Y.Yamada: add option compress and negative
  !                                        for compressing a floating data
  !                                            to an integer data and
  !                                        for negative value
  !                10-08-03   T.Mitsui: bug(?) fix for gfortran
  !                13-06-13   T.Seiki : application to z* coordinate,
  !                                     add option to use the nearest neighbor method
  !                14-05-09   C.Kodama : character(128) -> character(256)
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules (shared)
  !
  use mpi
  use mod_precision
  use mod_stdio
  use mod_cnst, only : &
       CNST_UNDEF,     & ! 05/12/21 M.Satoh
       CNST_UNDEF4
  use mod_calendar, only : &
       calendar_ym2dd
  ! Y.Yamada 09-09-09 ->
!  use mod_comp, only :&
!       comp_output,   &
!       comp_ctlgen
  !,   &
       !       comp_8to2
  ! Y.Yamada 09-09-09 <-
  implicit none

  character(256) :: title = 'NICAM_MODEL_OUTPUT' ! 06/08/04 H.Tomita
  character(256) :: llmap_dir = './'
  character(256) :: llmap_base = 'llmap'
  character(256) :: info_fname = 'output.info'
  logical :: lon_swap = .false.
  !
  real(8) :: z(1000)
  integer :: rlevel=0
  integer :: glevel=5
  integer :: kall = 1
  integer :: tstart = 1
  integer :: tend = 1
  integer :: kdum = 0
  character(256) :: input_dir = './'
  character(256) :: output_dir = './'
  integer, parameter :: allowable_vmax = 100
  character(256) :: base_name(allowable_vmax)
  character(256) :: header_strings ='' ! H.Tomita 06/08/04 add this option
  integer :: vmax = 1
  integer :: rec_counter
  integer :: input_size = 4  ! 05/12/21 M.Satoh
  ! S.Iga (051226)=>
  character(40) :: separate_time = 'no' ! no      : like *******.grd
  ! yes     : like *******.time01234.grd
  ! template: like *******.2005-12-01-23h50m.grd (not yet)
  character(256) :: timename
  integer :: init_year=0,init_month=1,init_day=1,init_hour=0,init_min=0   ! for the case template
  integer :: init_timenumber=1   ! for the case yes
  integer(8) :: initsecond, absosecond! for the case template
  integer :: absoday! for the case template
  character(20) :: sec2initplate,sec2template,timeincrement
  character(20) :: dc1,dc2
  ! S.Iga (051226)<=
  logical :: sepdir = .false. !S.Iga060508
  character(20) :: access_icodata = 'direct' ! 'direct' or 'sequential' ! 09/01/31 M.Satoh
  !                                          ! or 'sequential-time'     ! 09/03/16 H.Taniguchi

  logical :: opt_zstar2z=.false.                       ! 13/06/13 T.Seiki
  logical :: opt_lagintrpl=.true.                      ! 13/06/13 T.Seiki
  character(len=256) :: topo_fname="topog"             ! 13/06/13 T.Seiki
  character(len=256) :: vgrid_fname="vgrid_used.dat"   ! 13/06/13 T.Seiki
  logical :: opt_nearest_neighbor = .false.            ! 13/06/13 T.Seiki
  !
  namelist / ico2ll_param /&
       glevel,             &
       rlevel,             &
       input_dir,          &
       output_dir,         &
       info_fname,         &
       llmap_dir,          &
       llmap_base,         &
       lon_swap,           &
       tstart,             &  ! t  of icofile
       input_size,         &  ! 05/12/21 M.Satoh
       ! S.Iga 051226=>
       separate_time,      &
       init_year,          &  ! for the case template (it is used only for output filename and ctl file)
       init_month,         &  ! for the case template (it is used only for output filename and ctl file)
       init_day,           &  ! for the case template (it is used only for output filename and ctl file)
       init_hour,          &  ! for the case template (it is used only for output filename and ctl file)
       init_min,           &  ! for the case template (it is used only for output filename and ctl file)
       init_timenumber,    &  ! for the case yes (it is used only for output filename)
       sepdir,             &  ! whether the input directory is separated.
       ! S.Iga 051226<=
       header_strings,     &  ! header of file names
       access_icodata,     &  ! access type 09/01/31 M.Satoh
       opt_nearest_neighbor,& ! [Add] 13/06/13 T.Seiki
       opt_zstar2z,        &  ! [Add] 13/06/13 T.Seiki
       opt_lagintrpl,      &  ! [Add] 13/06/13 T.Seiki
       topo_fname,         &  ! [Add] 13/06/13 T.Seiki
       vgrid_fname            ! [Add] 13/06/13 T.Seiki
  !
  integer :: llmap_fid
  integer :: ierr
  !
  integer :: i,j
  integer :: imax,jmax
  real(8),allocatable :: lon(:),lat(:), lon2(:)
  real(8),allocatable :: topog(:,:)        ! [Add] 13/06/13 T.Seiki
  real(8),allocatable :: gz(:), gzh(:)     ! [Add] 13/06/13 T.Seiki
  real(8),allocatable :: vz(:,:), vzh(:,:) ! [Add] 13/06/13 T.Seiki
  real(8)             :: ztop=CNST_UNDEF4  ! [Add] 13/06/13 T.Seiki
  integer             :: kvgrid            ! [Add] 13/06/13 T.Seiki
  integer             :: kmin, kmax, kdim
  integer, parameter  :: khalo=1
  !
  character(256) :: fname

  integer :: t,l,k,v
  integer :: lall,gall
  integer :: it ! 09/03/16 H.Taniguchi

  integer :: n

  integer,allocatable :: lon_index(:),lat_index(:)
  integer,allocatable :: n1_index(:),n2_index(:),n3_index(:)
  real(8),allocatable :: w1(:),w2(:),w3(:)

  integer,allocatable :: max_num_latlon(:)
  integer,allocatable :: nstart(:), nend(:)

  real(4),allocatable :: ico_data(:,:)
  real(8),allocatable :: ico_data8(:,:)   ! 05/12/21 M.Satoh
  real(4),allocatable :: latlon_data(:,:,:)
  real(4),allocatable :: latlon_data_swap(:,:,:)
  !
  integer :: fid,ctl_fid,ofid,ifid
  !  character(256) :: ofname, ifname
  character(256) :: ifname
  !
  integer :: info_fid
  integer :: tall
  real(8) :: tintv
  !
  real(8) :: pi
  integer :: fnum
  integer :: ij ! 05/12/21 M.Satoh

  call MPI_Init(ierr)

  ctl_fid = IO_get_available_fid()
  open(CTL_FID,             &
       file='ico2ll.cnf',   &
       form='formatted',    &
       status='old',        &
       iostat=ierr)
  if(ierr/=0) then
     write(*,*) 'Cannot open PARAMETER file!'
     stop
  endif
  rewind(ctl_fid)
  read(ctl_fid,nml=ico2ll_param)
  close(ctl_fid)
  !
  info_fid = IO_get_available_fid()
  open(info_fid,file=trim(info_fname),form='formatted',status='old' ,iostat=ierr)
  if(ierr/=0) then
     write(*,*) 'nothing info file!'
     stop
  endif

  ! 05/12/02 M.Satoh: move from below
  lall = (2**rlevel)**2 * 10
  allocate(max_num_latlon(lall))
  allocate(nstart(lall))
  allocate(nend(lall))
  gall = (2**(glevel-rlevel)+2)*(2**(glevel-rlevel)+2)

  llmap_fid = IO_get_available_fid()
  open(llmap_fid,file=trim(trim(llmap_dir)//'/'//trim(llmap_base))//'.info',&
       form='unformatted',status='old' ,iostat=ierr)
  if(ierr/=0) then
     write(*,*) 'Cannot open llmap info file!'
     stop
  endif
  read(llmap_fid) imax
  allocate(lon(imax))
  read(llmap_fid) lon(:)
  read(llmap_fid) jmax
  allocate(lat(jmax))
  read(llmap_fid) lat(:)
  close(llmap_fid)

  allocate(lon_index(imax*jmax))
  allocate(lat_index(imax*jmax))
  allocate(n1_index(imax*jmax))
  allocate(n2_index(imax*jmax))
  allocate(n3_index(imax*jmax))
  allocate(w1(imax*jmax))
  allocate(w2(imax*jmax))
  allocate(w3(imax*jmax))
  !
  allocate(lon2(imax))
  allocate(topog(gall,lall)) ! [Add] 13/06/13 T.Seiki
  !
  do l= 1,lall
     call IO_make_idstr(fname,trim(trim(llmap_dir)//'/'//trim(llmap_base)),'rgn',l)
     ! write(*,*) trim(fname)
     fid = IO_get_available_fid()
     open(fid,file=trim(fname),form='unformatted',status='old')
     read(fid) max_num_latlon(l)
     nend(l) = sum(max_num_latlon(1:l))
     nstart(l) = nend(l) - max_num_latlon(l)+1
     if(max_num_latlon(l)/=0) then
        read(fid) lon_index(nstart(l):nend(l))
        read(fid) lat_index(nstart(l):nend(l))
        read(fid) n1_index(nstart(l):nend(l))
        read(fid) n2_index(nstart(l):nend(l))
        read(fid) n3_index(nstart(l):nend(l))
        read(fid) w1(nstart(l):nend(l))
        read(fid) w2(nstart(l):nend(l))
        read(fid) w3(nstart(l):nend(l))
     endif
     close(fid)
     ! [Add] 13/06/13 T.Seiki
     if ( opt_zstar2z )then
        call IO_make_idstr(fname,trim(topo_fname),'rgn',l)
        fid = IO_get_available_fid()
        open(fid,file=trim(fname),form='unformatted',status='old')
        read(fid) topog(:,l)
        close(fid)
     endif
     !
  enddo
  ! [Add] 13/06/13 T.Seiki
  if ( opt_zstar2z )then

     open(fid,file=trim(vgrid_fname), form="unformatted",status="old" )
     read(fid) kvgrid
     kmin=khalo+1
     kmax=khalo+kvgrid
     kdim=kmax+1
     allocate(gz (kdim))
     allocate(gzh(kdim))
     gz(:) =0.d0
     gzh(:)=0.d0
     read(fid) gz(1:kdim)
     read(fid) gzh(1:kdim)
     close(fid)
     ztop=gzh(kmax+1)-gzh(kmin)
     allocate(vz (gall,kdim))
     allocate(vzh(gall,kdim))
  endif

  !
  pi=4.0d0*atan(1.0d0)

  ! 05/12/02 M.Satoh
  fnum = 0
  do

     read(info_fid,*,end=1000) tall, tintv
     read(info_fid,*) kall
     do k=1,kall
        read(info_fid,*) z(k)
     enddo
     read(info_fid,*) vmax
     do v=1,vmax
        fnum = fnum + 1
        read(info_fid,*) base_name(v)
        write(*,*) 'file number=', fnum, base_name(v)
     enddo
     !
     allocate(ico_data(gall,kall))
     if ( input_size == 8 ) then ! 05/12/21 M.Satoh
        allocate(ico_data8(gall,kall))
     endif
     allocate(latlon_data(imax,jmax,kall))
     if(lon_swap) then
        ! T.Mitsui 07.03.29
!!$        allocate(latlon_data_swap(imax,jmax,kall))
        allocate(latlon_data_swap(imax,jmax,1))
     endif

     !S.Iga051226 =>
     call calendar_ym2dd(absoday,init_year,init_month,init_day)
! W.Yanase  2007/10/24   debug in calculation of 8-byte integer
     initsecond= absoday*int(86400,8) + init_hour*3600 + init_min*60
     !S.Iga051226 <=
     !
     tend = tstart + tall - 1
     !
     do v=1,vmax
        open(unit=fid,file=trim(output_dir)//'/'//trim(header_strings)//trim(base_name(v))//'.ctl',&
             form='formatted', status='replace')

        if (trim(separate_time) == 'no') then ! S.Iga051226
           write(unit=fid,fmt='(2a)')   'DSET ','^'//trim(header_strings)//trim(base_name(v))//'.grd'
           ! S.Iga051226=>
           elseif (trim(separate_time) == 'yes') then
           write(unit=fid,fmt='(2a)')   'DSET ','^'//trim(header_strings)//trim(base_name(v))//'.time0%y4'//'.grd'
           write(unit=fid,fmt='(1a)')   'OPTIONS TEMPLATE '
           elseif (trim(separate_time) == 'template') then
!   W. Yanase 081008   use %h2 instead of %f2 for GrADS template
           write(unit=fid,fmt='(2a)')   'DSET ','^'//trim(header_strings)//trim(base_name(v))//'.%y4-%m2-%d2-%h2h%n2m'//'.grd'
           write(unit=fid,fmt='(1a)')   'OPTIONS TEMPLATE '
        endif
        ! S.Iga051226<=

        write(unit=fid,fmt='(2a)')    'TITLE ',trim(title)
        write(unit=fid,fmt='(1a)')    'OPTIONS BIG_ENDIAN '
        write(unit=fid,fmt='(a,e12.5)')  'UNDEF ', real(-99.9E+33,4)
        write(unit=fid,fmt='(a,i5,a)')   'XDEF ',imax, ' LEVELS'
        if(lon_swap) then
           lon2(1:imax/2) = lon(imax/2+1:imax)
           lon2(imax/2+1:imax) = lon(1:imax/2)+2*pi
           write(unit=fid,fmt='(5x,5f10.4)')    (lon2(i)*180.0D0/pi,i=1,imax)
        else
           write(unit=fid,fmt='(5x,5f10.4)')    (lon(i)*180.0D0/pi,i=1,imax)
        endif
        write(unit=fid,fmt='(a,i5,a)')   'YDEF ',jmax, ' LEVELS'
        write(unit=fid,fmt='(5x,5f10.4)')    (lat(j)*180.0D0/pi,j=1,jmax)
        if(kall == 1 ) then
           write(unit=fid,fmt='(a,i5,a,2i5)')   'ZDEF ',kall, ' LINEAR',1,1
        else
           write(unit=fid,fmt='(a,i5,a)')   'ZDEF ',kall-2*kdum, ' LEVELS'
           write(unit=fid,fmt='(5f12.3)')  ((z(k)),k=kdum+1,kall-kdum)  ! 14-05-09 10->12
        endif
        if (trim(separate_time) == 'no') then ! S.Iga051226
           dc1 = trim(sec2initplate(initsecond)) ! S.Iga060508
           dc2 = trim(timeincrement(int(tintv))) ! S.Iga060508
           write(unit=fid,fmt='(a,i5,2a,1x,a)')  'TDEF ',tend,' LINEAR ',&
                trim(dc1),trim(dc2)
           ! S.Iga051226=>
        elseif (trim(separate_time) == 'yes') then
           !        write(unit=fid,fmt='(a,i5,2a,1x,a)')  'TDEF ',tend,' LINEAR ','00:00Z01JAN0000','1yr'

           !        S.Iga060106
           call make_idstr(timename,'','',init_timenumber)

           !        write(*,*) trim(timename)
           !        stop
           write(unit=fid,fmt='(a,i5,2a,1x,a)')  'TDEF ',tend,' LINEAR ','00:00Z01JAN'//timename(3:6),'1yr'
        elseif (trim(separate_time) == 'template') then
           dc1 = trim(sec2initplate(initsecond)) ! S.Iga060508
           dc2 = trim(timeincrement(int(tintv))) ! S.Iga060508
           write(unit=fid,fmt='(a,i5,2a,1x,a)')  'TDEF ',tend,' LINEAR ',&
                trim(dc1),trim(dc2)
        endif
        ! S.Iga051226<=
        write(unit=fid,fmt='(a,i5)')     'VARS ',1
        if(kall == 1) then
           write(unit=fid,fmt='(a,2i5,1x,a)')    trim(base_name(v)),0, 99, 'NONE'
        else
           write(unit=fid,fmt='(a,2i5,1x,a)')    trim(base_name(v)),kall-2*kdum, 99, 'NONE'
        endif
        write(unit=fid,fmt='(a)')        'ENDVARS '
        close(fid)
        !
        ofid = IO_get_available_fid()
        if (trim(separate_time) == 'no') then !S.Iga051226
           !
           open(ofid,file=trim(trim(output_dir)//'/'//trim(header_strings)//trim(base_name(v))&
                //'.grd'),  form='unformatted',access='direct', &
                recl=imax*jmax*4,status='unknown')
           rec_counter = (kall-2*kdum)*(tstart-1)+1
           !     write(*,*) '#### ',trim(base_name(v))
           ! 09-09-09 Y.Yamada ->
!           if( compress )then
!              ofid2b=IO_get_available_fid()
!              !              write(*,*) '#### ',trim(trim(output_dir)//'/'//trim(base_name(v))//'_2byte.grd'),ofid2b
!              open(ofid2b,file=trim(trim(output_dir)//'/'//trim(header_strings)//trim(base_name(v))&
!                   //'_2byte.grd'),form='unformatted')
!              write(ofid2b)dummy
!              call comp_ctlgen(   &
!                   output_dir,     &
!                   header_strings, &
!                   base_name(v),   &
!                   title,          &
!                   dc1,            &
!                   dc2,            &
!                   imax,           &
!                   jmax,           &
!                   kall,           &
!                   kdum,           &
!                   tend,           &
!                   z,              &
!                   lon,            &
!                   lat,            &
!                   pi,             &
!                   lon_swap        &
!                   )
!           endif
           ! 09-09-09 Y.Yamada <-
           !
        endif
        do t=tstart, tend
           !S.Iga051226=>
           if (trim(separate_time) == 'yes') then
              !        S.Iga060106
              call make_idstr(timename,&
                   trim(output_dir)//'/'//trim(header_strings)//trim(base_name(v)),'time',t-tstart+init_timenumber)
              open(ofid,file=trim(timename)//'.grd',&
                   form='unformatted',access='direct', &
                   recl=imax*jmax*4,status='unknown')
              rec_counter = 1
              !     write(*,*) '#### ',trim(base_name(v))
              elseif (trim(separate_time) == 'template') then
              absosecond = initsecond + (t-tstart) * tintv
              open(ofid,file=trim(output_dir)//'/'//trim(header_strings)//trim(base_name(v))//'.'&
                   //trim(sec2template(absosecond))//'.grd',&
                   form='unformatted',access='direct', &
                   recl=imax*jmax*4,status='unknown')
              rec_counter = 1
           endif
           !S.Iga051226<=

           do l=1, lall

!              write(*,*) 'l=', l, ifname

              !
              if (sepdir) then
                 call IO_make_idstr(ifname,trim(trim(input_dir)//'/'//&
                      trim(base_name(v))//'/'//trim(base_name(v))),'rgn',l)
              else
                 call IO_make_idstr(ifname,trim(trim(input_dir)//'/'//trim(base_name(v))),'rgn',l)
              endif
              !write(*,*) trim(ifname)
              ifid = IO_get_available_fid()

              ! 05/12/21 M.Satoh =>
              !           open(ifid,file=trim(ifname),form='unformatted',access='direct',recl=kall*gall*4,status='old')
              !              read(ifid,rec=t-tstart+1) ico_data(:,:)

              if ( access_icodata == 'sequential' ) then ! 09/01/31 M.Satoh
                 open(ifid,file=trim(ifname),form='unformatted', &
                      access='sequential', status='old')
              else if ( access_icodata == 'sequential-time' ) then ! 09/03/16 H.Taniguchi
                 open(ifid,file=trim(ifname),form='unformatted', &
                      access='sequential', status='old')
              else
                 open(ifid,file=trim(ifname),form='unformatted', &
                      access='direct',recl=kall*gall*input_size,status='old')
              endif

              if ( input_size == 4 ) then
                 if ( access_icodata == 'sequential' ) then ! 09/01/31 M.Satoh
                    read(ifid) ico_data(:,:)
                 else if ( access_icodata == 'sequential-time' ) then ! 09/03/16 H.Taniguchi
                    do it=1,t-tstart+1
                       read(ifid) ico_data(:,:)
                    enddo
                 else
                    read(ifid,rec=t-tstart+1) ico_data(:,:)
                 endif
              else if ( input_size == 8 ) then
                 if ( access_icodata == 'sequential' ) then ! 09/01/31 M.Satoh
                    read(ifid) ico_data8(:,:)
                 else if ( access_icodata == 'sequential-time' ) then ! 09/03/16 H.Taniguchi
                    do it=1,t-tstart+1
                       read(ifid) ico_data8(:,:)
                    enddo
                 else
                    read(ifid,rec=t-tstart+1) ico_data8(:,:)
                 endif

                 do k = 1, kall
                    do ij = 1, gall
                       if ( ico_data8(ij,k) == CNST_UNDEF ) then
                          ico_data(ij,k) = CNST_UNDEF4
                       else
                          ico_data(ij,k) = ico_data8(ij,k)
                       endif
                    enddo
                 enddo
              endif
              ! [Add] 13/06/13 T.Seiki
              ! vertical interpolation (z* => z coordinate)
              if( opt_zstar2z .and. kall == kvgrid )then
                 do k=1, kdim
                    do ij = 1, gall
                       vzh(ij,k) = topog(ij,l) + (ztop-topog(ij,l))/ztop*gzh(k)
                       vz(ij,k)  = topog(ij,l) + (ztop-topog(ij,l))/ztop*gz(k)
                    enddo
                 enddo
                 if( opt_lagintrpl )then ! lagrangian interpolation
                    call VINTRPL_z_level ( &
                         gall, kdim,       & ! in
                         kmin, kmax,       & ! in
                         gz, gzh,          & ! in zstar coordinate
                         vz, vzh,          & ! in z     coordinate
                         ico_data(:,1:kall)) ! inout
                 else                    ! linear interpolation
                    call VINTRPL_z_level2( &
                         gall, kdim,       & ! in
                         kmin, kmax,       & ! in
                         gz, gzh,          & ! in zstar coordinate
                         vz, vzh,          & ! in z     coordinate
                         ico_data(:,1:kall)) ! inout
                 endif
              endif
              ! 05/12/21 M.Satoh <=

!              write(*,'(A2,I5,A30,2E30.20)') 'l=', l, ifname, maxval(ico_data(1:gall,1:kall)), &
!                   minval(ico_data(1:gall,1:kall))

              !
              if(max_num_latlon(l)/=0) then
                 ! [Add] 13/06/13 T.Seiki
                 if( .not. opt_nearest_neighbor )then
                    do k=kdum+1,kall-kdum
                       do n = nstart(l),nend(l)
                          if(ico_data(n1_index(n),k)==CNST_UNDEF4) then
                             latlon_data(lon_index(n),lat_index(n),k) = CNST_UNDEF4
                          else if(ico_data(n2_index(n),k)==CNST_UNDEF4) then
                             latlon_data(lon_index(n),lat_index(n),k) = CNST_UNDEF4
                          else if(ico_data(n3_index(n),k)==CNST_UNDEF4) then
                             latlon_data(lon_index(n),lat_index(n),k) = CNST_UNDEF4
                          else
! del 06/08/21 M.Satoh
!                          latlon_data(lon_index(n),lat_index(n),k) &
!                               = w1(n)*ico_data(n1_index(n),k)     &
!                               + w2(n)*ico_data(n2_index(n),k)     &
!                               + w3(n)*ico_data(n3_index(n),k)
                             latlon_data(lon_index(n),lat_index(n),k) &
                                  = w1(n)*ico_data(n1_index(n),k)     &
                                  + w2(n)*ico_data(n2_index(n),k)     &
                                  + w3(n)*ico_data(n3_index(n),k)
                          endif
                       enddo
                    enddo
                 else
                    ! [Add] 13/06/13 T.Seiki, finding the value at the nearest neighbor grid
                    do k=kdum+1,kall-kdum
                       do n = nstart(l),nend(l)
                          if     (max(w1(n),w2(n),w3(n))==w1(n) .and. ico_data(n1_index(n),k) /= CNST_UNDEF4)then
                             latlon_data(lon_index(n),lat_index(n),k) = ico_data(n1_index(n),k)
                          else if(max(w2(n),w3(n)      )==w2(n) .and. ico_data(n2_index(n),k) /= CNST_UNDEF4)then
                             latlon_data(lon_index(n),lat_index(n),k) = ico_data(n2_index(n),k)
                          else if(                                    ico_data(n3_index(n),k) /= CNST_UNDEF4)then
                             latlon_data(lon_index(n),lat_index(n),k) = ico_data(n3_index(n),k)
                          else
                             latlon_data(lon_index(n),lat_index(n),k) = CNST_UNDEF4
                          endif
                       enddo
                    enddo
                 endif
              endif
              close(ifid)
              !
           enddo
           if(lon_swap) then
              ! Mod 07.03.29 T.Mitsui saving memory
              do k=1,kall
                 latlon_data_swap(1:imax/2,      1:jmax, 1) = latlon_data(imax/2+1:imax, 1:jmax, k)
                 latlon_data_swap(imax/2+1:imax, 1:jmax, 1) = latlon_data(1:imax/2     , 1:jmax, k)
                 latlon_data(:,:,k) = latlon_data_swap(:,:,1)
              enddo
!!$              latlon_data_swap(1:imax/2,1:jmax,1:kall) = latlon_data(imax/2+1:imax,1:jmax,1:kall)
!!$              latlon_data_swap(imax/2+1:imax,1:jmax,1:kall) = latlon_data(1:imax/2,1:jmax,1:kall)
!!$              latlon_data(:,:,:) = latlon_data_swap(:,:,:)
           endif
           do k=kdum+1,kall-kdum
              write(ofid,rec=rec_counter) latlon_data(:,:,k)
              !
!              if( compress )then
!                 call COMP_output(        &
!                      imax,               &
!                      jmax,               &
!                      rec_counter,        &
!                      ofid2b,             &
!                      latlon_data(:,:,k), &
!                      negative,           &
!                      k,&!XS
!                      base_name(v)        &
!                      )
!              endif
              rec_counter = rec_counter+1
           enddo
!           write(*,*) 't=',t,' : done', maxval(latlon_data), minval(latlon_data)
           write(*,*) 't=',t,' : done'
           if ((trim(separate_time) == 'yes').or.(trim(separate_time) == 'template')) then
              close(ofid)
           endif
        enddo
        if (trim(separate_time) == 'no') then
           close(ofid)
           ! 09-09-19 Y.Yamada ->
!           if(compress)then
!              close(ofid2b)
!           endif
        else
!           if(compress)then
!              close(ofid2b)
!           endif
           ! <-
        endif
     enddo

     ! 05/12/02 M.Satoh
     deallocate(ico_data)
     if ( input_size == 8 ) then ! 05/12/21 M.Satoh
        deallocate(ico_data8)
     endif
     deallocate(latlon_data)
     if(lon_swap) then
        deallocate(latlon_data_swap)
     endif

  enddo ! 05/12/02 M.Satoh
1000 continue
  close(info_fid)

  ! 05/12/02 M.Satoh
  deallocate(max_num_latlon)
  deallocate(nstart)
  deallocate(nend)
  deallocate(lon)
  deallocate(lat)
  deallocate(lon_index)
  deallocate(lat_index)
  deallocate(n1_index)
  deallocate(n2_index)
  deallocate(n3_index)
  deallocate(w1)
  deallocate(w2)
  deallocate(w3)
  deallocate(lon2)

  call MPI_Finalize(ierr)

  stop
contains
  ! [Add] 13/06/13  T.Seiki imported from mod_vintrpl.f90
  subroutine VINTRPL_z_level( ijdim, kdim, kmin, kmax, gz, gzh, vz, vzh, v )
    implicit none
    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim    ! corresponds to ADM_kall
    integer, intent(in)    :: kmin
    integer, intent(in)    :: kmax
    !
    real(8), intent(in)    :: gz(kdim)
    real(8), intent(in)    :: gzh(kdim)
    real(8), intent(in)    :: vz(ijdim,kdim)
    real(8), intent(in)    :: vzh(ijdim,kdim)
    !
    real(4), intent(inout) :: v(ijdim,kmin:kmax)
    real(8) :: tmp(ijdim,1:kdim)
    integer :: k,n,kk

    integer :: kp
    real(8) :: lag_intpl_quadra, lag_intpl_linear
    real(8) :: z,z1,p1,z2,p2,z3,p3
    lag_intpl_quadra(z,z1,p1,z2,p2,z3,p3)       &
         = ((z-z2)*(z-z3))/((z1-z2)*(z1-z3))*p1 &
         + ((z-z1)*(z-z3))/((z2-z1)*(z2-z3))*p2 &
         + ((z-z1)*(z-z2))/((z3-z1)*(z3-z2))*p3
    lag_intpl_linear(z,z1,p1,z2,p2) &
         = (z-z2)/(z1-z2)*p1 &
         + (z1-z)/(z1-z2)*p2
    !
    tmp(:,kmin:kmax) = v(:,kmin:kmax)
    tmp(:,1)         = v(:,kmin)
    tmp(:,kdim)      = v(:,kmax)
    !
    do k=kmin, kmax
       do n=1, ijdim
          if(      (gz(k)<vzh(n,kmin))&
               .or.(gz(k)>vzh(n,kmax+1)) ) then
             v(n,k) = CNST_UNDEF
          else
             do kk = kmin-1,kmax+1
                kp=kk
                if(gz(k)<vz(n,kk)) exit
             enddo
             if (kp<=kmin) then
                kp=kmin+1
             elseif(kp>=kmax) then
                kp=kmax-1
             endif
             !
             if ( tmp(n,kp+1) == CNST_UNDEF ) then
                if ( tmp(n,kp) == CNST_UNDEF ) then
                   v(n,k) = tmp(n,kp-1)
                else if ( tmp(n,kp-1) == CNST_UNDEF ) then
                   v(n,k) = tmp(n,kp)
                else
                   v(n,k) &
                        = lag_intpl_linear(gz(k), &
                        vz(n,kp  ),tmp(n,kp  ),   &
                        vz(n,kp-1),tmp(n,kp-1)    )
                endif
             else if ( tmp(n,kp) == CNST_UNDEF ) then
                v(n,k) = tmp(n,kp-1)
             else if ( tmp(n,kp-1) == CNST_UNDEF ) then
                v(n,k) = tmp(n,kp)
             else
                v(n,k)                        &
                     = lag_intpl_quadra(gz(k),&
                     vz(n,kp+1),tmp(n,kp+1),  &
                     vz(n,kp  ),tmp(n,kp  ),  &
                     vz(n,kp-1),tmp(n,kp-1) )
             endif
          endif
       enddo
    enddo
    return
  end subroutine VINTRPL_z_level
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_z_level2( ijdim, kdim, kmin, kmax, gz, gzh, vz, vzh, v )
    implicit none
    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: kmin
    integer, intent(in)    :: kmax
    !
    real(8), intent(in)    :: gz(kdim)
    real(8), intent(in)    :: gzh(kdim)
    real(8), intent(in)    :: vz(ijdim,kdim)
    real(8), intent(in)    :: vzh(ijdim,kdim)
    !
    real(4), intent(inout) :: v(ijdim,kmin:kmax)
    !
    real(8) :: tmp(ijdim,kdim)
    integer :: k,n,kk
    integer :: kp
    real(8) :: lag_intpl
    real(8) :: z,z1,p1,z2,p2
    lag_intpl(z,z1,p1,z2,p2)  &
         = (z-z2)/(z1-z2)*p1  &
         + (z1-z)/(z1-z2)*p2
    tmp(:,kmin:kmax) = v(:,kmin:kmax)
    tmp(:,1)         = v(:,kmin)
    tmp(:,kdim)      = v(:,kmax)
    !
    do k=kmin, kmax
       do n=1, ijdim
          if(      (gz(k)<vzh(n,kmin)  ) &
               .or.(gz(k)>vzh(n,kmax+1)) ) then
             v(n,k) = CNST_UNDEF
          else
             do kk = kmin-1, kmax+1
                kp=kk
                if(gz(k)<vz(n,kk)) exit
             enddo
             if    (kp<=kmin  ) then
                v(n,k)=tmp(n,kmin)
             elseif(kp>=kmax+1) then
                v(n,k)=tmp(n,kmax)
             else
                !
                if ( tmp(n,kp) == CNST_UNDEF ) then
                   v(n,k) = tmp(n,kp-1)
                else if ( tmp(n,kp-1) == CNST_UNDEF ) then
                   v(n,k) = tmp(n,kp)
                else
                   v(n,k)                        &
                        = lag_intpl(gz(k),       &
                        vz(n,kp  ),tmp(n,kp  ),  &
                        vz(n,kp-1),tmp(n,kp-1) )
                endif
             endif
          endif
       enddo
    enddo
  end subroutine VINTRPL_z_level2

  !
end program ico2ll

!S.Iga051226 =>
!-- output grads-like template part  like 01JAN0000
function sec2initplate(absosecond) & !in
     result(plate)
  !
  use mod_calendar, only : &
       calendar_dd2ym
  !
  implicit none
  integer(8) :: absosecond
  integer :: iday,day,year,month,min,hour
  ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: plate
  character(len=20):: plate
  character(4):: nyear
  character(2):: nday, nhour,nmin
  character(3):: nmonth(12)
  data nmonth /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/
  !
  iday= absosecond/86400
  hour=mod(absosecond,int(86400,8))/3600 !! Y.Niwa 86400 -> int(86400,8)
  min=mod(absosecond,int(3600,8))/60 !! Y.Niwa 3600 -> int(3600,8)
  call calendar_dd2ym(year,month,day,iday)
  write(nyear,'(I4.4)') year
  write(nday,'(I2.2)') day
  write(nmin,'(I2.2)') min
  write(nhour,'(I2.2)') hour
  plate=nhour//':'//nmin//'Z'//nday//nmonth(month)//nyear
  !
end function sec2initplate


!-- output grads-like template part  like 2005-12-01-23h50m
function sec2template(absosecond) & !in
     result(plate)
  !
  use mod_calendar, only : &
       calendar_dd2ym
  !
  implicit none
  integer(8) :: absosecond
  integer :: iday,day,year,month,min,hour
  ! [Mod] 10/08/03 T.Mitsui
!!$  character(*):: plate
  character(len=20):: plate
  character(4):: nyear
  character(2):: nday, nhour,nmin,nmonth
  !
  iday= absosecond/86400
  hour=mod(absosecond,int(86400,8))/3600 !! Y.Niwa 86400 -> int(86400,8)
  min=mod(absosecond,int(3600,8))/60 !! Y.Niwa 3600 -> int(3600,8)
  call calendar_dd2ym(year,month,day,iday)
  write(nyear,'(I4.4)') year
  write(nday,'(I2.2)') day
  write(nmin,'(I2.2)') min
  write(nhour,'(I2.2)') hour
  write(nmonth,'(I2.2)') month
  plate=nyear//'-'//nmonth//'-'//nday//'-'//nhour//'h'//nmin//'m'
  !
end function sec2template

function timeincrement(isec)&
     result(plate)
  !
  implicit none
  integer :: isec
  character(20):: plate
  character(18):: tmp
  write(tmp,*) max(isec/60, 1)
  plate=trim(tmp)//'mn'
  !
end function timeincrement

!S.Iga051226 <=

!S.Iga060106 =>
subroutine make_idstr( &
     str,                   & !--- INOUT : file name
     headstr,               & !--- IN : header name
     ext,                   & !--- IN : extention string
     rank )                   !--- IN : ID number
  !
  implicit none
  !
  character(*), intent(inout) :: str   !--- strings
  character(*), intent(in) :: headstr  !--- header strings
  character(*), intent(in) :: ext      !--- extention( eg. ***.dat )
  integer, intent(in) :: rank          !--- number(+1)
  !
  character(10) :: cnum
  character(5) cnum1
  !
  write(cnum,'(I10.10)') rank
  cnum1(1:5) = cnum(6:10)
  str=headstr//'.'//trim(ext)//cnum1
  !
end subroutine make_idstr

