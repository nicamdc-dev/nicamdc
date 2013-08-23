!-------------------------------------------------------------------------------
!
!+  program ico2ll (NEW I/O)
!
!-------------------------------------------------------------------------------
program fio_ico2ll
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This program converts from data on dataicosahedral grid (new I/O format)
  !       to that on latitude-longitude grid.
  !       (some part of source code is imported from ico2ll.f90)
  !
  !++ Current Corresponding Author : H.Yashiro
  !
  !++ Contributer of ico2ll.f90 : M.Satoh, S.Iga, Y.Niwa, H.Tomita, T.Mitsui,
  !                               W.Yanase,  H.Taniguchi, Y.Yamada
  !
  !++ History: 
  !      Version   Date      Comment 
  !      -----------------------------------------------------------------------
  !      0.90      11-09-07  H.Yashiro : [NEW] partially imported from ico2ll.f90
  !      0.95      12-04-19  H.Yashiro : [mod] deal large record length
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_misc, only : &
    MISC_get_available_fid, &
    MISC_make_idstr
  use mod_cnst, only : &
    CNST_UNDEF, &
    CNST_UNDEF4
  use mod_calendar, only : &
    calendar_ss2yh
  use mod_fio, only : &
    FIO_HSHORT,       &
    FIO_HMID,         &
    FIO_HLONG,        &
    FIO_REAL4,        &
    FIO_REAL8,        &
    FIO_BIG_ENDIAN,   &
    FIO_ICOSAHEDRON,  &
    FIO_IGA_LCP,      &
    FIO_IGA_MLCP,     &
    FIO_INTEG_FILE,   &
    FIO_SPLIT_FILE,   &
    FIO_FREAD,        &
    headerinfo,       &
    datainfo
  use mod_mnginfo_light, only : &
    MNG_mnginfo_input,   &
    MNG_mnginfo_noinput, &
    MNG_PALL,            &
    MNG_prc_rnum,        &
    MNG_prc_tab
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  integer, parameter :: max_nvar   = 500
  integer, parameter :: max_nstep  = 1500
  integer, parameter :: max_nlayer = 200

  integer, parameter :: flim = 1
  integer,      save :: fmax

  !--- NAMELIST
  integer                   :: glevel              = -1
  integer                   :: rlevel              = -1
  character(LEN=FIO_HSHORT) :: grid_topology       = 'ICOSAHEDRON'
                                                   ! 'LCP'
                                                   ! 'MLCP'
  logical                   :: complete            = .false.
  character(LEN=FIO_HLONG)  :: mnginfo             = ''
  character(LEN=FIO_HLONG)  :: layerfile_dir       = ''
  character(LEN=FIO_HLONG)  :: llmap_base          = ''
  character(LEN=FIO_HLONG)  :: infile(flim)        = ''
  integer                   :: step_str            = 1
  integer                   :: step_end            = max_nstep
  character(LEN=FIO_HLONG)  :: outfile_dir         = '.'
  character(LEN=FIO_HSHORT) :: outfile_prefix      = ''
  integer                   :: outfile_rec         = 1
  logical                   :: lon_swap            = .false.
  logical                   :: devide_template     = .false.
  logical                   :: output_grads        = .true.
  logical                   :: output_gtool        = .false.
  logical                   :: datainfo_nodep_pe   = .false.   !   <- can be .true. if data header do not depend on pe.
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''
  character(LEN=FIO_HSHORT) :: large_memory_var(max_nvar) = '' ! [add] 13-04-18

  logical                   :: help = .false.

  namelist /OPTION/ glevel,            &
                    rlevel,            &
                    grid_topology,     &
                    complete,          &
                    mnginfo,           &
                    layerfile_dir,     &
                    llmap_base,        &
                    infile,            &
                    step_str,          &
                    step_end,          &
                    outfile_dir,       &
                    outfile_prefix,    &
                    outfile_rec,       &
                    lon_swap,          &
                    devide_template,   &
                    output_grads,      &
                    output_gtool,      &
                    datainfo_nodep_pe, &  ! [add] 13-04-18
                    selectvar,         &
                    large_memory_var,  &  ! [add] 13-04-18
                    help
  !-----------------------------------------------------------------------------
  character(LEN=FIO_HLONG) :: infname   = ""
  character(LEN=FIO_HLONG) :: outbase   = ""
  character(LEN=FIO_HLONG) :: layerfile = ""
  integer                  :: fmode
  integer                  :: gtopology
  logical                  :: allvar = .true.

  ! ll grid coordinate
  integer              :: imax, jmax
  real(8), allocatable :: lon(:), lat(:)
  real(8), allocatable :: lon_tmp(:) ! [add] 13-04-18

  ! ico2ll weight mapping
  integer              :: LALL
  integer, allocatable :: num_llgrid(:)
  integer, allocatable :: llstr(:), llend(:)
  integer, allocatable :: lon_idx(:), lat_idx(:)
  integer, allocatable :: n1(:), n2(:), n3(:)
  real(8), allocatable :: w1(:), w2(:), w3(:)

  ! ico data information
  integer, allocatable :: ifid(:)
  integer, allocatable :: prc_tab_C(:)
  type(headerinfo) hinfo 
  type(datainfo)   dinfo 

  integer                                :: num_of_data
  integer                                :: nvar
  character(LEN=FIO_HSHORT), allocatable :: var_name(:)
  character(LEN=FIO_HMID),   allocatable :: var_desc(:)
  character(LEN=FIO_HSHORT), allocatable :: var_unit(:)
  character(LEN=FIO_HSHORT), allocatable :: var_layername(:)
  integer,                   allocatable :: var_datatype(:)
  integer,                   allocatable :: var_nlayer(:)
  integer,                   allocatable :: var_nstep(:)
  integer(8),                allocatable :: var_time_str(:)
  integer(8),                allocatable :: var_dt(:)
  real(8),                   allocatable :: var_zgrid(:,:)
  ! header
  character(LEN=16),         allocatable :: var_gthead(:,:)

  ! ico data
  integer              :: GALL 
  real(4), allocatable :: data4allrgn(:)
  real(8), allocatable :: data8allrgn(:)
  real(4), allocatable :: icodata4(:,:,:)

  ! ll data
  real(4), allocatable :: lldata(:,:,:)
  real(4), allocatable :: temp(:,:)


  character(LEN=FIO_HLONG) :: fname
  character(LEN=20)        :: tmpl
  character(LEN=16)        :: gthead(64)
  integer(8)               :: nowsec
  integer(8)               :: recsize ! [mod] 12-04-19 H.Yashiro
  integer                  :: kmax, num_of_step, step, date_str(6)

  logical :: addvar
  integer :: fid, did, ofid, ierr, irec
  integer :: v, t, p, l, k, n, i, j
  real(8) :: pi
  !=============================================================================

  pi = 4.D0 * atan( 1.D0 ) ! [add] 13-04-18

  !--- read option and preprocess
  call readoption !! set fmax, infile

  if ( glevel==-1 .or. rlevel==-1 ) then
     write(*,*) "xxx Set glevel, rlevel. STOP"
     stop
  endif
  if ( step_str < 1 .or. step_end < 1 ) then
     write(*,*) "xxx step must be >= 1. STOP"
     stop
  elseif( step_str > step_end ) then
     write(*,*) "xxx step_str must be < step_end. STOP"
     stop
  endif

  if ( grid_topology=="ICOSAHEDRON" ) then
     gtopology = FIO_ICOSAHEDRON
  elseif( grid_topology=="LCP" ) then
     gtopology = FIO_IGA_LCP
  elseif( grid_topology=="MLCP" ) then
     gtopology = FIO_IGA_MLCP
  else
     write(*,*) "Unknown type of Grid toporogy:",grid_topology
     stop
  endif

  if (output_gtool) then
     output_grads    = .false.
     devide_template = .false.
     outfile_rec = 1
  endif

  if ( trim(selectvar(1)) /= '' ) then
     allvar = .false.
  endif

  !--- prepare region infomation
  if (complete) then ! all region
    fmode = FIO_INTEG_FILE
    call MNG_mnginfo_noinput( rlevel )
  else               ! region specified by mnginfo
    fmode = FIO_SPLIT_FILE
    call MNG_mnginfo_input( rlevel, trim(mnginfo) )
  endif

  !--- setup
  call fio_syscheck()

  !#########################################################
  !--- Read lat-lon grid information
  fid = MISC_get_available_fid()
  open(fid, file=trim(llmap_base)//'.info',form='unformatted',status='old',iostat=ierr)
     if (ierr/=0) then
        write(*,*) 'Cannot open llmap info file!',trim(llmap_base)//'.info'
        stop
     endif

     read(fid) imax
     allocate(lon(imax))
     read(fid) lon(:)
     read(fid) jmax
     allocate(lat(jmax))
     read(fid) lat(:)
  close(fid)

  !--- Read lat-lon weight map
  !write(*,*) '*** io_mode for llmap : LEGACY'

  GALL = ( (2**(glevel-rlevel))+2 ) &
       * ( (2**(glevel-rlevel))+2 )

  LALL = 10 * (2**rlevel) * (2**rlevel)

  allocate( num_llgrid(LALL) )
  allocate( llstr     (LALL) )
  allocate( llend     (LALL) )

  allocate( lon_idx(imax*jmax) )
  allocate( lat_idx(imax*jmax) )
  allocate( n1     (imax*jmax) )
  allocate( n2     (imax*jmax) )
  allocate( n3     (imax*jmax) )
  allocate( w1     (imax*jmax) )
  allocate( w2     (imax*jmax) )
  allocate( w3     (imax*jmax) )

  do l = 1, LALL
     call MISC_make_idstr(fname,trim(llmap_base),'rgn',l)
     fid = MISC_get_available_fid()
     open(fid,file=trim(fname),form='unformatted',status='old',iostat=ierr)
        if (ierr/=0) then
           write(*,*) 'Cannot open llmap file!',trim(fname)
           stop
        endif

        read(fid) num_llgrid(l)
        llend(l) = sum( num_llgrid(1:l) )
        llstr(l) = llend(l) - num_llgrid(l) + 1

        if ( num_llgrid(l)/=0 ) then
           read(fid) lon_idx( llstr(l):llend(l) )
           read(fid) lat_idx( llstr(l):llend(l) )
           read(fid) n1     ( llstr(l):llend(l) )
           read(fid) n2     ( llstr(l):llend(l) )
           read(fid) n3     ( llstr(l):llend(l) )
           read(fid) w1     ( llstr(l):llend(l) )
           read(fid) w2     ( llstr(l):llend(l) )
           read(fid) w3     ( llstr(l):llend(l) )
        endif
     close(fid)
  enddo

  ! Read icodata information
  !write(*,*) '*** io_mode for data  : ADVANCED'

  allocate( ifid(MNG_PALL) )
  allocate( var_nstep    (max_nvar) )
  allocate( var_name     (max_nvar) )
  allocate( var_desc     (max_nvar) )
  allocate( var_unit     (max_nvar) )
  allocate( var_layername(max_nvar) )
  allocate( var_datatype (max_nvar) )
  allocate( var_nlayer   (max_nvar) )
  allocate( var_time_str (max_nvar) )
  allocate( var_dt       (max_nvar) )
  allocate( var_zgrid    (max_nlayer, max_nvar) )
  allocate( var_gthead   (64, max_nvar) )


  do v = 1, max_nvar
     if ( trim(large_memory_var(v)) == '' ) exit
     !write(*,*) v, trim(large_memory_var(v)), len_trim(large_memory_var(v))
     call fio_register_vname_tmpdata( trim(large_memory_var(v)), &
          len_trim(large_memory_var(v)) )
  end do

  do p = 1, MNG_PALL
     write(*,*) 'p=', p

     if (complete) then ! all region
        infname = trim(infile(1))//'.rgnall'
     else
        call fio_mk_fname(infname,trim(infile(1)),'pe',p-1,6)
     endif
     allocate( prc_tab_C(MNG_prc_rnum(p))   )
     prc_tab_C(:) = MNG_prc_tab(:,p)-1

     call fio_put_commoninfo( fmode,           &
                              FIO_BIG_ENDIAN,  &
                              gtopology,       &
                              glevel,          &
                              rlevel,          &
                              MNG_prc_rnum(p), &
                              prc_tab_C        )

     call fio_register_file(ifid(p),trim(infname))
     call fio_fopen(ifid(p),FIO_FREAD)
     
     ! <-- [add] C.Kodama 13.04.18
     if( datainfo_nodep_pe .and. p > 1 ) then
        ! assume that datainfo do not depend on pe.
        call fio_read_pkginfo( ifid(p) )
        call fio_valid_pkginfo( ifid(p) )
        call fio_copy_datainfo( ifid(p), ifid(1) )
     else if ( .not. datainfo_nodep_pe  &
               .and. trim(large_memory_var(1)) /= '' ) then
        ! also read field data and store them in temporary buffer
        call fio_read_allinfo_tmpdata( ifid(p) )
     else
        ! normal way to read pkginfo and datainfo
        call fio_read_allinfo( ifid(p) )
     end if
     ! -->
     
     if ( p == 1 ) then ! only once
        allocate( hinfo%rgnid(MNG_prc_rnum(p)) )

        call fio_get_pkginfo(ifid(p),hinfo)

        num_of_data = hinfo%num_of_data
        write(*,*) '*** get variable informations'
        write(*,*) 'num_of_data    : ', num_of_data

        nvar = 0
        do did = 0, num_of_data-1
           call fio_get_datainfo(ifid(p),did,dinfo)

           if (allvar) then ! output all variables
              addvar = .true.
           else             ! select valiables to output
              addvar = .false.

              do v = 1, max_nvar
                 if ( trim(selectvar(v)) == trim(dinfo%varname) ) then
                    addvar = .true.
                    exit
                 elseif( trim(selectvar(v)) == '' ) then
                    exit
                 endif
              enddo
           endif

           do v = 1, nvar
              if ( trim(var_name(v)) == trim(dinfo%varname) ) then
                 var_nstep(v) = var_nstep(v) + 1

                 if( var_nstep(v) == 2 ) var_dt(v) = dinfo%time_start - var_time_str(v)

                 if( var_nstep(v) == step_str ) var_time_str(v) = dinfo%time_start ! [mod] H.Yashiro 20111003

                 addvar = .false.
                 exit
              endif
           enddo

           if (addvar) then
              nvar = nvar + 1
              var_nstep    (nvar) = 1
              var_name     (nvar) = dinfo%varname
              var_desc     (nvar) = dinfo%description
              var_unit     (nvar) = dinfo%unit
              var_layername(nvar) = dinfo%layername
              var_datatype (nvar) = dinfo%datatype
              var_nlayer   (nvar) = dinfo%num_of_layer
              var_time_str (nvar) = dinfo%time_start
              var_dt       (nvar) = dinfo%time_end - dinfo%time_start

              layerfile = trim(layerfile_dir)//'/'//trim(dinfo%layername)//'.txt'

              fid = MISC_get_available_fid()
              open(fid,file=trim(layerfile),form='formatted',status='old',iostat=ierr)
                 if ( ierr /= 0 ) then
                    write(*,*) 'xxx layerfile doesnt exist!', trim(layerfile)
                    stop
                 endif

                 read(fid,*) kmax
                 do k = 1, kmax
                    read(fid,'(F16.4)') var_zgrid(k,nvar)
                 enddo
              close(fid)
           endif

        enddo !--- did LOOP
     endif !--- PE=000000

     deallocate( prc_tab_C )
     call fio_fclose(ifid(p)) ! [add] 13-04-18
  enddo !--- PE LOOP

  if ( nvar == 0 ) then
     write(*,*) 'No variables to convert. Finish.'
     stop
  endif

  write(*,*) '########## Variable List ########## '
  write(*,*) 'ID |NAME            |STEPS|Layername       |START FROM         |DT [sec]'
  do v = 1, nvar
     call calendar_ss2yh( date_str(:), real(var_time_str(v),kind=8) )
     write(tmpl,'(I4.4,"/",I2.2,"/",I2.2,1x,I2.2,":",I2.2,":",I2.2)') date_str(:)
     write(*,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A19,A1,I8)') &
              v,'|',var_name(v),'|',var_nstep(v),'|',var_layername(v),'|', tmpl,'|', var_dt(v)
  enddo

  write(*,*) '*** convert start : PaNDa format to lat-lon data'

  !#########################################################
  !--- start weighting summation
  do v = 1, nvar

     kmax = var_nlayer(v)
     recsize = int(imax,kind=8)*int(jmax,kind=8)*int(kmax,kind=8)*4_8 ! [mod] 12-04-19 H.Yashiro

     !--- open output file
     outbase = trim(outfile_dir)//'/'//trim(outfile_prefix)//trim(var_name(v))
     ofid = MISC_get_available_fid()
 
     num_of_step = min(step_end,var_nstep(v)) - step_str + 1  ! [mov] 13-04-18

     if (.not. devide_template) then
        if (output_grads) then

           write(*,*) 'Output: ', trim(outbase)//'.grd', recsize, imax, jmax, kmax
           open( unit   = ofid,                  &
                 file   = trim(outbase)//'.grd', &
                 form   = 'unformatted',         &
                 access = 'direct',              &
                 recl   = recsize,               &
                 status = 'unknown'              )
           irec = 1

           if ( outfile_rec > 1 ) then
              write(*,*) 'Change output record position : start from step ', outfile_rec
              irec = outfile_rec
           endif

        elseif(output_gtool) then

           write(*,*) 'Output: ', trim(outbase)//'.gt3'
           open( unit   = ofid,                  &
                 file   = trim(outbase)//'.gt3', &
                 form   = 'unformatted',         &
                 access = 'sequential',          &
                 status = 'unknown'              )

           ! [mod] H.Yashiro 20111003
           call makegtoolheader( var_gthead(:,v),     &
                                 outfile_dir,         &
                                 var_name(v),         &
                                 var_desc(v),         &
                                 var_unit(v),         &
                                 var_layername(v),    &
                                 imax,                &
                                 jmax,                &
                                 var_nlayer(v),       &
                                 lon,                 &
                                 lat,                 &
                                 var_zgrid(1:kmax,v), &
                                 var_dt(v),           &
                                 lon_swap             )

        endif
     endif


     do t = 1, num_of_step

        allocate( lldata(imax,jmax,kmax) )
        lldata(:,:,:) = CNST_UNDEF4

        nowsec = var_time_str(v) + (t-1)*var_dt(v)

        !--- open output file (every timestep)
        if (devide_template) then
           tmpl   = sec2template(nowsec)
           write(*,*)
           write(*,*) 'Output: ', trim(outbase)//'.'//trim(tmpl)//'.grd', recsize

           open( unit   = ofid,             &
                 file   = trim(outbase)//'.'//trim(tmpl)//'.grd', &
                 form   = 'unformatted',    &
                 access = 'direct',         &
                 recl   = recsize,          &
                 status = 'unknown'         )
           irec = 1 
        endif

        step = t-1 + step_str

        do p = 1, MNG_PALL
           call fio_fopen(ifid(p),FIO_FREAD)  ! [add] 13-04-18
           if ( t==1 ) write(*,'(A10)',advance='no') ' ->region:'

           allocate( data4allrgn(GALL*kmax*MNG_prc_rnum(p)) )
           allocate( data8allrgn(GALL*kmax*MNG_prc_rnum(p)) )
           allocate( icodata4   (GALL,kmax,MNG_prc_rnum(p)) )
           data4allrgn(:)  = CNST_UNDEF4
           data8allrgn(:)  = CNST_UNDEF
           icodata4(:,:,:) = CNST_UNDEF4

           !--- seek data ID and get information
           call fio_seek_datainfo(did,ifid(p),var_name(v),step)
           call fio_get_datainfo(ifid(p),did,dinfo)

           !--- verify
           if ( did == -1 ) then
              write(*,*) 'xxx data not found! varname:',trim(var_name(v)),", step : ",step
              stop
           endif

           !--- read from pe000xx file


           if ( dinfo%datatype == FIO_REAL4 ) then
              if ( trim(large_memory_var(1)) /= '' ) then
                 call fio_read_data_tmpdata(ifid(p),did,data4allrgn(:))
              else
                 call fio_read_data(ifid(p),did,data4allrgn(:))
              end if

           elseif( dinfo%datatype == FIO_REAL8 ) then
              if ( trim(large_memory_var(1)) /= '' ) then
                 call fio_read_data_tmpdata(ifid(p),did,data8allrgn(:))
              else
                 call fio_read_data(ifid(p),did,data8allrgn(:))
              end if

              data4allrgn(:) = real(data8allrgn(:),kind=4)
              where( data8allrgn(:) == CNST_UNDEF )
                 data4allrgn(:) = CNST_UNDEF4
              endwhere

           endif
           icodata4(:,:,:) = reshape( data4allrgn(:), shape(icodata4) )

           do l = 1, MNG_prc_rnum(p)
              if ( t==1 ) write(*,'(1x,I5.5)',advance='no') MNG_prc_tab(l,p)

              !--- ico -> lat-lon
              if ( num_llgrid(MNG_prc_tab(l,p)) /= 0 ) then
                 do k = 1, kmax
                 do n = llstr(MNG_prc_tab(l,p)), llend(MNG_prc_tab(l,p))
                    if (      icodata4(n1(n),k,l) == CNST_UNDEF4 &
                         .OR. icodata4(n2(n),k,l) == CNST_UNDEF4 &
                         .OR. icodata4(n3(n),k,l) == CNST_UNDEF4 ) then

                       lldata(lon_idx(n),lat_idx(n),k) = CNST_UNDEF4
                    else
                       lldata(lon_idx(n),lat_idx(n),k) = real(w1(n),kind=4) * icodata4(n1(n),k,l) &
                                                       + real(w2(n),kind=4) * icodata4(n2(n),k,l) &
                                                       + real(w3(n),kind=4) * icodata4(n3(n),k,l)
                    endif
                 enddo
                 enddo
              endif

           enddo ! region LOOP

           if ( t==1 ) write(*,*)
           deallocate( data4allrgn )
           deallocate( data8allrgn )
           deallocate( icodata4    )
           call fio_fclose(ifid(p)) ! [add] 13-04-18
        enddo ! PE LOOP

        !--- swap longitude
        if (lon_swap) then
           allocate( temp(imax,jmax) )
           do k = 1, kmax ! Mod 07.03.29 T.Mitsui saving memory 
              temp(1:imax/2,     :) = lldata(imax/2+1:imax,:,k)
              temp(imax/2+1:imax,:) = lldata(1:imax/2     ,:,k)
              lldata(:,:,k)         = temp(:,:)
           enddo
           deallocate( temp )
        endif

        !--- output lat-lon data file
        if (output_grads) then
           write(ofid,rec=irec) lldata(:,:,:)
           irec = irec + 1
        elseif(output_gtool) then
           write(var_gthead(25,v),'(I16)') int( nowsec/3600,kind=4 )
           write(var_gthead(27,v),'(A16)') calendar_ss2cc_gtool(nowsec)
           gthead(:) = var_gthead(:,v)

           write(ofid) gthead(:)
           write(ofid) lldata(:,:,:)
        endif

        !--- close output file
        if (devide_template) then
           close(ofid)
        endif

        write(*,*) ' +append step:', step
        deallocate( lldata )
     enddo ! step LOOP

     !--- close output file
     if (.not. devide_template) then
        close(ofid)
     endif

     if (output_grads) then
        call makegradsctl( outfile_dir,         &
                           outfile_prefix,      &
                           var_name(v),         &
                           imax,                &
                           jmax,                &
                           kmax,                &
                           lon,                 &
                           lat,                 &
                           var_zgrid(1:kmax,v), &
                           var_nstep(v),        &
                           var_time_str(v),     &
                           var_dt(v),           &
                           lon_swap,            &
                           devide_template      )
     endif

  enddo ! variable LOOP

! [del] 13-04-18
!  do p = 1, MNG_PALL
!     call fio_fclose(ifid(p))
!  enddo

contains
  !-----------------------------------------------------------------------------
  !> read option
  !-----------------------------------------------------------------------------
  subroutine readoption
    use mod_misc, only : &
      MISC_get_available_fid
    use mod_tool_option, only: &
      OPT_convert, &
      OPT_fid
    implicit none

    integer :: io
    !---------------------------------------------------------------------------

    ! --- Set option
    OPT_fid = MISC_get_available_fid()
    open(OPT_fid,status='SCRATCH')

      call OPT_convert( fmax )

      read(OPT_fid,nml=OPTION,iostat=io)

    close(OPT_fid)

    if (      io /= 0     &
         .OR. fmax == 0   &
         .OR. fmax > flim &
         .OR. help        ) call helpoption

  end subroutine readoption

  !-----------------------------------------------------------------------------
  !> display help for option and abort
  !-----------------------------------------------------------------------------
  subroutine helpoption
    implicit none
    !---------------------------------------------------------------------------

    write(*,OPTION)

    stop
  end subroutine helpoption

  !-----------------------------------------------------------------------------
  subroutine makegradsctl( &
      outfile_dir,    &
      outfile_prefix, &
      varname,        &
      imax,           &
      jmax,           &
      kmax,           &
      lon,            &
      lat,            &
      alt,            &
      nstep,          &
      time_str,       &
      dt,             &
      lon_swap,       &
      devide_template )
    implicit none

    character(LEN=128), intent(in) :: outfile_dir
    character(LEN=16),  intent(in) :: outfile_prefix
    character(LEN=16),  intent(in) :: varname
    integer,            intent(in) :: imax
    integer,            intent(in) :: jmax
    integer,            intent(in) :: kmax
    real(8),            intent(in) :: lon(imax)
    real(8),            intent(in) :: lat(jmax)
    real(8),            intent(in) :: alt(kmax)
    integer,            intent(in) :: nstep
    integer(8),         intent(in) :: time_str
    integer(8),         intent(in) :: dt
    logical,            intent(in) :: lon_swap
    logical,            intent(in) :: devide_template
 
    real(8) :: pi
    real(8) :: temp(imax)

    character(LEN=32)  :: outfile
    integer            :: fid
    character(LEN=20)  :: s1, s2
    integer            :: i, j, k
    !---------------------------------------------------------------------------
    pi = 4.D0 * atan( 1.D0 )

    outfile = trim(outfile_prefix)//trim(var_name(v))

    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(outfile_dir)//'/'//trim(outfile)//'.ctl', &
          form   = 'formatted', &
          status = 'replace'    )

       ! S.Iga051226=>
       if ( devide_template ) then
          ! W. Yanase 081008   use %h2 instead of %f2 for GrADS template
          write(fid,'(A)') 'DSET ^'//trim(outfile)//'.%y4-%m2-%d2-%h2h%n2m'//'.grd'
          write(fid,'(A)') 'OPTIONS TEMPLATE '
       else
          write(fid,'(A)') 'DSET ^'//trim(outfile)//'.grd'
       endif
       ! S.Iga051226<=

       write(fid,'(A)')      'TITLE NICAM data output'
       write(fid,'(A)')      'OPTIONS BIG_ENDIAN '
       write(fid,'(A,E12.5)') 'UNDEF ', real( -99.9E+33, kind=4 )

       write(fid,'(A,I5,A)') 'XDEF ', imax, ' LEVELS'
       if (lon_swap) then
          temp(1:imax/2)      = lon(imax/2+1:imax)
          temp(imax/2+1:imax) = lon(1:imax/2) + 2*pi
          write(fid,'(10(1x,F9.4))') (temp(i)*180.D0/pi,i=1,imax)
       else
          write(fid,'(10(1x,F9.4))') (lon(i) *180.D0/pi,i=1,imax)
       endif

       write(fid,'(A,I5,A)')    'YDEF ',jmax, ' LEVELS'
       write(fid,'(10(1x,F9.4))') (lat(j)*180.D0/pi,j=1,jmax)

       if ( kmax == 1 ) then
          write(fid,'(A,I5,A,2I5)') 'ZDEF ', kmax, ' LINEAR', 1, 1
       else
          write(fid,'(A,I5,A)') 'ZDEF ', kmax, ' LEVELS'
          write(fid,'(10(1x,F9.2))') (alt(k),k=1,kmax)
       endif

       s1 = trim( sec2initplate(time_str) ) ! S.Iga060508
       s2 = trim( timeincrement(int(dt)) )  ! S.Iga060508
       write(fid,'(A,I5,2A,1x,A)') 'TDEF ',nstep, ' LINEAR ', trim(s1), trim(s2)

       write(fid,'(a,i5)') 'VARS ', 1
       if ( kmax == 1 ) then
          write(fid,'(a,2i5,1x,a)') trim(varname), 0, 99, 'NONE'
       else
          write(fid,'(a,2i5,1x,a)') trim(varname), kmax, 99, 'NONE'
       endif
       write(fid,'(a)') 'ENDVARS '
    close(fid)

    write(*,'(A,A)') 'Generate ',trim(outfile)//'.ctl'

  end subroutine makegradsctl

  !-----------------------------------------------------------------------------
  subroutine makegtoolheader( &
      gthead,      & !--- OUT
      outfile_dir, &
      varname,     &
      description, &
      unit,        &
      layername,   &
      imax,        &
      jmax,        &
      kmax,        &
      lon,         &
      lat,         &
      alt,         &
      dt,          &
      lon_swap     )
    implicit none

    character(LEN=16),         intent(out) :: gthead(64)
    character(LEN=FIO_HLONG),  intent( in) :: outfile_dir
    character(LEN=FIO_HSHORT), intent( in) :: varname
    character(LEN=FIO_HMID),   intent( in) :: description
    character(LEN=FIO_HSHORT), intent( in) :: unit
    character(LEN=FIO_HSHORT), intent( in) :: layername
    integer,                   intent( in) :: imax
    integer,                   intent( in) :: jmax
    integer,                   intent( in) :: kmax
    real(8),                   intent( in) :: lon(imax)
    real(8),                   intent( in) :: lat(jmax)
    real(8),                   intent( in) :: alt(kmax)
    integer(8),                intent( in) :: dt
    logical,                   intent( in) :: lon_swap
 
    character(LEN=16) :: axhead(64)
    character(LEN=16) :: hitem
    character(LEN=32) :: htitle
    character(LEN=16) :: gt_axisx
    character(LEN=16) :: gt_axisy
    character(LEN=16) :: kdate

    integer           :: ndttm(8)
    character(LEN=10) :: ndate, ntime, nzone

    real(8) :: pi
    real(8) :: temp(imax)
    real(8) :: dx, lonp1(imax+1)
    integer :: i
    !---------------------------------------------------------------------------
    pi = 4.D0 * atan( 1.D0 )

    hitem  = trim(varname)
    do i=1,16
       if( hitem(i:i)=='_' ) hitem(i:i)  = '-' ! escape underbar
    enddo
    htitle(1:32) = description(1:32) ! trim to 32char
    do i=1,32
       if( htitle(i:i)=='_' ) htitle(i:i) = '-' ! escape underbar
    enddo

    write(gt_axisx,'(A,I4.4)') 'LON', imax
    write(gt_axisy,'(A,I4.4)') 'LAT', jmax

    call date_and_time(ndate, ntime, nzone, ndttm)
    write(kdate,'(I4.4,I2.2,I2.2,1x,I2.2,I2.2,I2.2,1x)') ndttm(1),ndttm(2),ndttm(3),ndttm(5),ndttm(6),ndttm(7)

    gthead(:) = ' '
    write(gthead( 1),'(I16)'  ) 9010
    write(gthead( 2),'(A16)'  ) 'NICAM'
    write(gthead( 3),'(A16)'  ) hitem
    write(gthead(12),'(I16)'  ) 1
    write(gthead(13),'(I16)'  ) 1
    write(gthead(14),'(A16)'  ) htitle(1:16)
    write(gthead(15),'(A16)'  ) htitle(17:32)
    write(gthead(16),'(A16)'  ) unit

    write(gthead(26),'(A16)'  ) 'HOUR            '
    write(gthead(28),'(I16)'  ) int(dt/3600,kind=4)
    write(gthead(29),'(A16)'  ) gt_axisx ! from info file
    write(gthead(30),'(I16)'  ) 1
    write(gthead(31),'(I16)'  ) imax
    write(gthead(32),'(A16)'  ) gt_axisy ! from info file
    write(gthead(33),'(I16)'  ) 1
    write(gthead(34),'(I16)'  ) jmax
    write(gthead(35),'(A16)'  ) layername
    write(gthead(36),'(I16)'  ) 1
    write(gthead(37),'(I16)'  ) kmax
    write(gthead(38),'(A16)'  ) 'UR4'
    write(gthead(39),'(E16.7)') real(-99.9E+33,4)
    write(gthead(40),'(E16.7)') real(-99.9E+33,4)
    write(gthead(41),'(E16.7)') real(-99.9E+33,4)
    write(gthead(42),'(E16.7)') real(-99.9E+33,4)
    write(gthead(43),'(E16.7)') real(-99.9E+33,4)
    write(gthead(44),'(I16)'  ) 1
    write(gthead(46),'(I16)'  ) 0
    write(gthead(47),'(E16.7)') 0. 
    write(gthead(48),'(I16)'  ) 0
    write(gthead(60),'(A16)'  ) kdate
    write(gthead(62),'(A16)'  ) kdate
    write(gthead(61),'(A16)'  ) 'NICAM'
    write(gthead(63),'(A16)'  ) 'NICAM'
    write(gthead(64),'(I16)'  ) imax*jmax*kmax

    !--- Generate axis file
    axhead(:) = ' '
    write(axhead( 1),'(I16)'  ) 9010
    write(axhead( 2),'(A16)'  ) 'AXLOC'
    write(axhead(12),'(I16)'  ) 1
    write(axhead(13),'(I16)'  ) 1
    write(axhead(25),'(I16)'  ) 0
    write(axhead(26),'(A16)'  ) 'SEC'
    write(axhead(27),'(A16)'  ) kdate
    write(axhead(28),'(I16)'  ) 1
    write(axhead(30),'(I16)'  ) 1
    write(axhead(33),'(I16)'  ) 1
    write(axhead(34),'(I16)'  ) 1
    write(axhead(36),'(I16)'  ) 1
    write(axhead(37),'(I16)'  ) 1
    write(axhead(38),'(A16)'  ) 'UR4'
    write(axhead(39),'(E16.7)') -999.0
    write(axhead(44),'(I16)'  ) 1
    write(axhead(46),'(I16)'  ) 0
    write(axhead(47),'(E16.7)') 0. 
    write(axhead(48),'(I16)'  ) 0
    write(axhead(60),'(A16)'  ) kdate
    write(axhead(62),'(A16)'  ) kdate
    write(axhead(61),'(A16)'  ) 'NICAM'
    write(axhead(63),'(A16)'  ) 'NICAM'

    fid = MISC_get_available_fid()
    open( unit   = fid,           &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(gt_axisx)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then

       dx = lon(2)-lon(1)

       if (lon_swap) then
          temp(1:imax/2)      = lon(imax/2+1:imax)
          temp(imax/2+1:imax) = lon(1:imax/2) + 2*pi

          lonp1(1) = temp(1) - dx/2
       else
          lonp1(1) = lon(1)  - dx/2
       endif
       if ( abs(lonp1(1)) < 1.D-10 ) lonp1(1) = 0.D0

       do i = 2, imax+1
          lonp1(i) = lonp1(i-1) + dx
       enddo

       write(axhead( 3),'(A16)'  ) trim(gt_axisx)
       write(axhead(29),'(A16)'  ) trim(gt_axisx)
       write(axhead(31),'(I16)'  ) imax+1
       write(axhead(40),'(E16.7)')   0.E0
       write(axhead(41),'(E16.7)') 360.E0
       write(axhead(42),'(E16.7)')  10.E0
       write(axhead(43),'(E16.7)')  30.E0
       write(axhead(64),'(I16)'  ) imax+1

       write(fid) axhead
       write(fid) real( lonp1(1:imax+1)/pi*180.D0, kind=4 )

       close(fid)
    endif

    open( unit   = fid,          &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(gt_axisy)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then
       write(axhead( 3),'(A16)'  ) trim(gt_axisy)
       write(axhead(29),'(A16)'  ) trim(gt_axisy)
       write(axhead(31),'(I16)'  ) jmax
       write(axhead(40),'(E16.7)') -90.E0
       write(axhead(41),'(E16.7)')  90.E0
       write(axhead(42),'(E16.7)')  10.E0
       write(axhead(43),'(E16.7)')  30.E0
       write(axhead(64),'(I16)'  ) jmax

       write(fid) axhead
       write(fid) real( lat(1:jmax)/pi*180.D0, kind=4 )

       close(fid)
    endif

    open( unit   = fid,          &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(layername)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then
       write(axhead( 3),'(A16)'  ) trim(layername)
       write(axhead(29),'(A16)'  ) trim(layername)
       write(axhead(31),'(I16)'  ) kmax
       write(axhead(40),'(E16.7)')     0.E0
       write(axhead(41),'(E16.7)') real(maxval(alt),kind=4)
       write(axhead(42),'(E16.7)')  1000.E0
       write(axhead(43),'(E16.7)') 10000.E0
       write(axhead(64),'(I16)'  ) kmax

       write(fid) axhead
       write(fid) real( alt(1:kmax), kind=4 )

       close(fid)
    endif


  end subroutine makegtoolheader

  !S.Iga051226 =>
  !-----------------------------------------------------------------------------
  function sec2initplate(datesec) result(template)
    !-- output grads-like template part  like 01JAN0000
    implicit none

    integer(8)        :: datesec
    ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: template
    character(LEN=20) :: template

    integer :: d(6)

    character(LEN=3) :: nmonth(12)
    data nmonth / 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC' /
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call calendar_ss2yh( d(:), real(datesec,kind=8) )

    write(template,'(I2.2,A1,I2.2,A1,I2.2,A3,I4.4)') &
                              d(4), ':', d(5), 'Z', d(3), nmonth(d(2)), d(1)

  end function sec2initplate

  !-----------------------------------------------------------------------------
  function sec2template(datesec) result(template)
    !-- output grads-like template part  like 2005-12-01-23h50m
    implicit none

    integer(8)        :: datesec
    ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: template
    character(LEN=20) :: template

    integer :: d(6)
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call calendar_ss2yh( d(:), real(datesec,kind=8) )

    write(template,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1)') &
                          d(1), '-', d(2), '-', d(3), '-', d(4), 'h', d(5), 'm'

  end function sec2template

  !-----------------------------------------------------------------------------
  function timeincrement(isec) result(template)
    implicit none

    integer       :: isec
    character(20) :: template

    character(18):: tmp
    !---------------------------------------------------------------------------

    write(tmp,*) max(isec/60, 1)

    template = trim(tmp)//'mn'

  end function timeincrement
  !S.Iga051226 <=

  !-----------------------------------------------------------------------------
  function calendar_ss2cc_gtool(datesec) result(template)
    !--- calendar, sec. -> character (YYYYMMDD HHMMSS)
    implicit none

    integer(8)        :: datesec
    character(LEN=16) :: template

    integer :: d(6), i
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call calendar_ss2yh( d(:), real(datesec,kind=8) )

    write (template,'(i4.4,i2.2,i2.2,1x,i2.2,i2.2,i2.2,1x)') (d(i),i=1,6)

  end function calendar_ss2cc_gtool

end program fio_ico2ll
!-------------------------------------------------------------------------------

