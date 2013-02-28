!-------------------------------------------------------------------------------
!
!+  Program ll2ico
!
!-------------------------------------------------------------------------------
program ll2ico
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This program is the conversion program from lat-lon data
  !       to icosahedral data.
  ! 
  !++ Current Code Owner: H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      09-??-??   H.Tomita  : Add this program.
  !                12-06-04   R.Yoshida : add direct access output, and add "kmax=k+2"
  !                12-06-07   R.Yoshida : allow large number of record length
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules (shared)
  !
  use mod_adm, only : &
       ADM_NSYS,      &
       ADM_GALL_PL,   &
       ADM_LALL_PL
  use mod_grd, only : &
       GRD_XDIR,      &
       GRD_YDIR,      &
       GRD_ZDIR
  use mod_misc, only : &
       misc_get_available_fid,&
       misc_make_idstr, &
       MISC_get_latlon

  implicit none

  character(128) :: lldata_fname = 'restart_LLGDATA.dat'
  character(128) :: lldata_dim_fname = 'LLGDATA_dimension.info'
  !
  character(128) :: hgrid_fname = 'grid'
  character(128) :: output_data = 'rst'
  integer :: glevel = 5 
  integer :: rlevel = 1 
  integer :: num_of_rgn, num_of_gp


  integer :: ctl_fid, fid, fid_grd,fid_odat
  integer :: i,j,k,l,n,t
  character(ADM_NSYS) :: z_coordinate_type
  real(8), allocatable :: lon(:), lat(:), zh(:), z(:), dz(:)
  real(8), allocatable :: data(:,:,:)  
  integer :: imax,jmax,kall, tmax, time_increment, kmax
  character(128) :: fname
  !
  integer :: ierr
  real(8),allocatable :: x_ico(:,:),lat_ico(:), lon_ico(:), data_ico(:,:)
  real(8),allocatable :: x_ico_pl(:,:,:),lat_ico_pl(:,:), lon_ico_pl(:,:), data_ico_pl(:,:,:)
  integer :: gall_1d
  real(8) :: dlon, dlat
  real(8) :: a,b,c,d
  !
  integer(8) :: irecl                        ! for long record length ! 12/06/07 R.Yoshida
  character(20) :: access_outdata = 'direct' ! 'direct' or 'sequential' ! 12/06/01 R.Yoshida
  logical :: read_conf = .false.

  namelist / ll2icoparam / &
       lldata_fname,       &
       lldata_dim_fname,   &
       hgrid_fname,        &
       output_data,        &
       access_outdata,     &
       glevel,             &
       rlevel,             &
       read_conf
  !
  !--- open and read configuration file
  ctl_fid = 10
  open(CTL_FID,               &
       file='ll2ico.cnf',     &
       form='formatted',      &
       status='old'           &
       )
  read(CTL_FID,nml=ll2icoparam,iostat=ierr)
  write(*,nml=ll2icoparam) 
  if(ierr<0) then
     write(*,*) &
          'Msg : Prg[ll2ico]'
     write(*,*) &
          ' *** Not found namelist.'
     write(*,*) &
          ' *** Use default values.'
  else if(ierr>0) then
     write(*,*) &
          'Msg : Prg[ll2ico]'
     write(*,*) &
          ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
  end if
  close(ctl_fid)
  !
  !---- read lat-lon-z-t dimension data
  fid = MISC_get_available_fid()
  open(fid,file=lldata_dim_fname,form='unformatted',status='old')
  !
  !------ longitude
  read(fid) imax
  allocate(lon(-1:imax+2))
  read(fid) lon(1:imax)
  !
  !------ latitude
  read(fid) jmax
  allocate(lat(jmax))
  read(fid) lat(:) 
  !
  !------ z
  read(fid) z_coordinate_type
  read(fid) kall
  kmax = kall + 2       ! 12/06/04 R.Yoshida
  !allocate(zh(kall+1))  ! 12/06/04 R.Yoshida
  allocate(z(kmax))
  allocate(zh(kmax))
  read(fid) z(:)
  read(fid) zh(:)
  !
  allocate(dz(kall+1))
  !allocate(z(kall))  ! 12/06/04 R.Yoshida
  do k=1,kall+1
     !z(k) = 0.5D0*(zh(k+1)+zh(k))  ! 12/06/04 R.Yoshida
     dz(k) = abs(zh(k+1)-zh(k))
  end do
  !
  !------ time
  read(fid) tmax
  read(fid) time_increment
  !
  close(fid)
  !
  !--- debug info 
  write(*,*) 'data structure :',imax,jmax,kall,kmax
  write(*,*) 'number of data :',tmax
  !
  !--- set incremental value and longitudinal hallo(2)
  dlon=(lon(imax)-lon(1))/(dble(imax)-1.0D0)
  dlat=(lat(jmax)-lat(1))/(dble(jmax)-1.0D0)
  lon(0) = lon(1)-dlon
  lon(-1) = lon(0)-dlon
  lon(imax+1) = lon(imax)+dlon
  lon(imax+2) = lon(imax+1)+dlon
  !
  !--- memory allocation for latlon data
  !allocate(data(-1:imax+2,jmax,kall))
  allocate(data(-1:imax+2,jmax,kmax))    ! 12/06/04 R.Yoshida
  !
  !--- setup icosahedral info and memory allocation
  num_of_rgn=(2**rlevel)*(2**rlevel)*10
  num_of_gp=(2**(glevel-rlevel)+2)*(2**(glevel-rlevel)+2)
  allocate(x_ico(num_of_gp,GRD_XDIR:GRD_ZDIR))
  allocate(lat_ico(num_of_gp))
  allocate(lon_ico(num_of_gp))
!  allocate(data_ico(num_of_gp,kall))
  allocate(data_ico(num_of_gp,kmax))    ! 12/06/04 R.Yoshida

  allocate(x_ico_pl(ADM_GALL_PL,ADM_LALL_PL,GRD_XDIR:GRD_ZDIR))
  allocate(lat_ico_pl(ADM_GALL_PL,ADM_LALL_PL))
  allocate(lon_ico_pl(ADM_GALL_PL,ADM_LALL_PL))
!  allocate(data_ico_pl(ADM_GALL_PL,kall,ADM_LALL_PL))
  allocate(data_ico_pl(ADM_GALL_PL,kmax,ADM_LALL_PL))    ! 12/06/04 R.Yoshida
  !
  !--- open latlon data
  fid = MISC_get_available_fid()
  irecl = int(imax,kind=8)*int(jmax,kind=8)*int(kmax,kind=8)*8_8    ! 12/06/07 R.Yoshida
!  open(fid,file=lldata_fname,form='unformatted',status='old')
!  open(fid,file=lldata_fname,form='unformatted',status='old',access='direct',recl=8*imax*jmax*kall)
  open(fid,file=lldata_fname,form='unformatted',status='old',access='direct',recl=irecl)    ! 12/06/04 R.Yoshida
  !
  write (*,'(A)') "---------------------------"  ! 12/06/01 R.Yoshida
  !
  !--- main loop
  do t=1,tmax
     !
     write (*,'(" Time Loop:",I4)') t   ! 12/06/07 R.Yoshida
     !
     !--- read data of 1 layer and fix to cyclic condition (longitudinal)
     read(fid,rec=t) data(1:imax,:,:)                                 ! 12/06/04 R.Yoshida
     if( read_conf ) write (*,'(" Read Conf. (t=",I4,", k=1)  Max:",F16.8," Min:",F16.8)') t, &  ! 12/06/01 R.Yoshida
     maxval(data(1:imax,:,:)), minval(data(1:imax,:,:))    ! 12/06/04 R.Yoshida
!     do k=1,kall
!        read(fid) data(1:imax,1:jmax,k)
!        if(k == 1)then  ! 12/06/01 R.Yoshida
!           write (*,'(" Read Conf. (t=",I3,", k=1)  Max:",F10.4," Min:",F10.4)') t, maxval(data), minval(data)
!        endif
!     end do
     data(0,:,:) = data(imax,:,:)
     data(-1,:,:) = data(imax-1,:,:)
     data(imax+1,:,:) = data(1,:,:)
     data(imax+2,:,:) = data(2,:,:)
     !
     !--- regular region
     do l = 1, num_of_rgn
        call MISC_make_idstr(fname,trim(hgrid_fname),'rgn',l)
        fid_grd = MISC_get_available_fid()
        open(fid_grd,file=trim(fname),status='old',form='unformatted',iostat=ierr)
        if(ierr/=0) then
           write(*,*) &
                '   *** No grid file.',trim(fname)
           stop
        end if
        read(fid_grd) gall_1d
        if(gall_1d**2/=num_of_gp) then
           write(*,*) 'inconsistent in number of gridpoint!'
           stop
        end if
        !
        read(fid_grd) x_ico(:,GRD_XDIR)
        read(fid_grd) x_ico(:,GRD_YDIR)
        read(fid_grd) x_ico(:,GRD_ZDIR)
        !
        !
        close(fid_grd)
        !
        !------ setup lat-lon on icosahedral gridpoints   
        do n=1,num_of_gp
           call MISC_get_latlon(     &
                lat_ico(n), lon_ico(n),    &
                x_ico(n,GRD_XDIR), x_ico(n,GRD_YDIR), x_ico(n,GRD_ZDIR) )
        end do
        !
        !------ bi-linear interpolation
        do n=1, num_of_gp
           i = floor((lon_ico(n)-lon(1))/dlon)+1
           j = min(max(floor((lat_ico(n)-lat(1))/dlat)+1,1),jmax-1)
           a = (lon(i+1)-lon_ico(n))/dlon
           b = (lon_ico(n)-lon(i))/dlon
           c = (lat(j+1)-lat_ico(n))/dlat
           d = (lat_ico(n)-lat(j))/dlat
           !do k=1,kall
           do k=1,kmax       ! 12/06/04 R.Yoshida
              data_ico(n,k) &
                   = c*(a*data(i,j  ,k)+b*data(i+1,j  ,k)) &
                   + d*(a*data(i,j+1,k)+b*data(i+1,j+1,k))
           end do
        end do
        !
        !------ output 
        call MISC_make_idstr(fname,trim(output_data),'rgn',l)
        fid_odat = MISC_get_available_fid()
!        open(fid_odat,file=trim(fname),access='append',form='unformatted')
        if ( access_outdata == 'sequential' ) then ! 12/06/01 R.Yoshida
           open(fid_odat,file=trim(fname),access='sequential',position='append', form='unformatted')
           write(fid_odat) data_ico(:,:)
        else
           !open(fid_odat,file=trim(fname),access='direct',recl=8*num_of_gp*kall,form='unformatted')
           open(fid_odat,file=trim(fname),access='direct',recl=8*num_of_gp*kmax,form='unformatted')   ! 12/06/04 R.Yoshida
           write(fid_odat, rec=t) data_ico(:,:)
        endif
        !
        close(fid_odat)
     end do
     !
     !--- pole region
     fname=trim(hgrid_fname)//'.pl'
     fid_grd = MISC_get_available_fid()
     open(fid_grd,file=trim(fname),status='old',form='unformatted',iostat=ierr)
     if(ierr/=0) then
        write(*,*) &
             '   *** No grid file.',trim(fname)
        stop
     end if
     !
     read(fid_grd) x_ico_pl(:,:,:)
     !
     close(fid_grd)
     !
     !------ setup lat-lon on icosahedral gridpoints   
     do l=1,ADM_LALL_PL
        do n=1,ADM_GALL_PL
           call MISC_get_latlon(                    &
                lat_ico_pl(n,l), lon_ico_pl(n,l),    &
                x_ico_pl(n,l,GRD_XDIR), x_ico_pl(n,l,GRD_YDIR), x_ico_pl(n,l,GRD_ZDIR) )
        end do
     end do
     !
     !------ bi-linear interpolation
     do l=1,ADM_LALL_PL
        do n=1,ADM_GALL_PL
           i = floor((lon_ico_pl(n,l)-lon(1))/dlon)+1
           j = min(max(floor((lat_ico_pl(n,l)-lat(1))/dlat)+1,1),jmax-1)
           a = (lon(i+1)-lon_ico_pl(n,l))/dlon
           b = (lon_ico_pl(n,l)-lon(i))/dlon
           c = (lat(j+1)-lat_ico_pl(n,l))/dlat
           d = (lat_ico_pl(n,l)-lat(j))/dlat
           !do k=1,kall
           do k=1,kmax      ! 12/06/04 R.Yoshida
              data_ico_pl(n,k,l) &
                   = c*(a*data(i,j  ,k)+b*data(i+1,j  ,k)) &
                   + d*(a*data(i,j+1,k)+b*data(i+1,j+1,k))
           end do
        end do
     end do
     !
     !------ output 
     fname=trim(output_data)//'.pl'
     fid_odat = MISC_get_available_fid()
!     open(fid_odat,file=trim(fname),access='append',form='unformatted')
     open(fid_odat,file=trim(fname),access='sequential',position='append',form='unformatted')
     write(fid_odat) data_ico_pl(:,:,:)
     close(fid_odat)
  end do
  !
  close(fid)



end program ll2ico
