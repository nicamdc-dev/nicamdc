!-------------------------------------------------------------------------------
!> Program FIO pe2pe
!!
!! @par Description
!!          combine pe0000x format data
!!          ( packaged NICAM data format : PaNDa )
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
program fio_pe2pe
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use iso_c_binding
  use mod_fio_common
  use mod_stdio
  use mod_mnginfo_light, only: &
     MNG_mnginfo_input, &
     MNG_PALL,          &
     MNG_prc_rnum,      &
     MNG_prc_tab,       &
     MNG_rgn2prc
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ C interface
  !
  include 'mod_fio_panda.inc'

  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  integer, parameter :: flim = 1
  integer,      save :: fmax

  !--- NAMELIST
  integer               :: glevel       = -1
  integer               :: rlevel_in    = -1
  character(len=H_LONG) :: mnginfo_in   = ''
  character(len=H_LONG) :: infile(flim) = ''
  integer               :: rlevel_out   = -1
  character(len=H_LONG) :: mnginfo_out  = ''
  character(len=H_LONG) :: outfile      = ''
  logical               :: use_mpi      = .true.
  logical               :: help         = .false.

  namelist /OPTION/ glevel,      &
                    rlevel_in,   &
                    mnginfo_in,  &
                    infile,      &
                    rlevel_out,  &
                    mnginfo_out, &
                    outfile,     &
                    use_mpi,     &
                    help

  !-----------------------------------------------------------------------------
  character(len=H_LONG)  :: fname      = ''
  character(len=H_LONG)  :: fname_dump = ''
  integer                :: oprec

  integer                :: num_of_rgn_in
  integer                :: num_of_rgn_out

  integer                :: MNG_PALL_in
  integer, allocatable   :: MNG_prc_rnum_in (:)
  integer, allocatable   :: MNG_prc_tab_in  (:,:)
  integer, allocatable   :: rgnid_in        (:)
  integer                :: MNG_PALL_out
  integer, allocatable   :: MNG_prc_rnum_out(:)
  integer, allocatable   :: MNG_prc_tab_out (:,:)
  integer, allocatable   :: rgnid_out       (:)

  integer                :: ndivide
  integer, allocatable   :: rgnid_div(:)

  type(headerinfo_panda)            :: hinfo
  type(datainfo_panda), allocatable :: dinfo(:)

  character(len=H_SHORT) :: var_name_file
  character(len=H_MID)   :: pkg_desc
  character(len=H_LONG)  :: pkg_note
  integer                :: nmax_data

  integer                :: GALL_in
  integer                :: GALL_out
  integer                :: KALL
  integer                :: LALL_in
  integer                :: LALL_out
  real(SP), target, allocatable :: data4_1D(:)
  real(DP), target, allocatable :: data8_1D(:)
  real(SP), allocatable  :: data4_3D_in (:,:,:)
  real(SP), allocatable  :: data4_3D_out(:,:,:)
  real(DP), allocatable  :: data8_3D_in (:,:,:)
  real(DP), allocatable  :: data8_3D_out(:,:,:)

  ! for MPI
  integer                :: prc_nall, prc_nlocal
  integer                :: prc_myrank, pstr, pend
  integer                :: fid_log
  character(len=6)       :: rankstr

  integer  :: p, l, n
  integer  :: fid, did, odid, ierr, fid_dump
  !=====================================================================

  !--- read option and preprocess
  call readoption !! set fmax, infile

  !--- prepare region infomation
  num_of_rgn_in = 10 * 4**rlevel_in
  GALL_in       = ( (2**(glevel-rlevel_in))+2 ) &
                * ( (2**(glevel-rlevel_in))+2 )

  call MNG_mnginfo_input( rlevel_in, trim(mnginfo_in) )
  MNG_PALL_in = MNG_PALL
  LALL_in     = MNG_prc_rnum(1)

  allocate( MNG_prc_rnum_in(MNG_PALL_in) )
  allocate( MNG_prc_tab_in (num_of_rgn_in,MNG_PALL_in) )
  MNG_prc_rnum_in(:)   = MNG_prc_rnum(:)
  MNG_prc_tab_in (:,:) = MNG_prc_tab (:,:)

  deallocate( MNG_prc_rnum )
  deallocate( MNG_prc_tab  )
  deallocate( MNG_rgn2prc  )

  num_of_rgn_out = 10 * 4**rlevel_out
  GALL_out       = ( (2**(glevel-rlevel_out))+2 ) &
                 * ( (2**(glevel-rlevel_out))+2 )

  call MNG_mnginfo_input( rlevel_out, trim(mnginfo_out) )
  MNG_PALL_out = MNG_PALL
  LALL_out     = MNG_prc_rnum(1)

  allocate( MNG_prc_rnum_out(MNG_PALL_out) )
  allocate( MNG_prc_tab_out (num_of_rgn_out,MNG_PALL_out) )
  MNG_prc_rnum_out(:)   = MNG_prc_rnum(:)
  MNG_prc_tab_out (:,:) = MNG_prc_tab (:,:)

  MNG_PALL = -1
  deallocate( MNG_prc_rnum )
  deallocate( MNG_prc_tab  )
  deallocate( MNG_rgn2prc  )

  allocate( rgnid_in (LALL_in)  )
  allocate( rgnid_out(LALL_out) )

  ndivide = 2**(rlevel_out-rlevel_in)
  allocate( rgnid_div(ndivide*ndivide) )

  fid_log = IO_get_available_fid()
  if ( use_mpi ) then
     !--- Parallel Excution, No communication
     call MPI_Init(ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, prc_nall,   ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD, prc_myrank, ierr)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     write(rankstr,'(I6.6)') prc_myrank
     open(fid_log, file='msg_p2p.pe'//trim(rankstr) )
     write(fid_log,*) '+++ Parallel Execution, Use MPI'

     if ( mod( MNG_PALL_in, prc_nall) /= 0 ) then
        write(fid_log,*) '*** Invalid processor number, STOP:', MNG_PALL_in, prc_nall
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
     elseif( mod( MNG_PALL_out, prc_nall) /= 0 ) then
        write(fid_log,*) '*** Invalid processor number, STOP:', MNG_PALL_out, prc_nall
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
     endif
  else
     open(fid_log, file='msg.serial' )
     write(fid_log,*) '+++ Serial Execution'
     prc_nall   = 1
     prc_myrank = 0
  endif

  !--- setup
  ierr = fio_syscheck()

  write(fid_log,*) '*** pe2pe start : PaNDa format to PaNDa format data'

  prc_nlocal = MNG_PALL_in / prc_nall
  pstr       = prc_myrank*prc_nlocal + 1
  pend       = prc_myrank*prc_nlocal + prc_nlocal
  write(fid_log,*)
  write(fid_log,*) '*** For Input *** '
  write(fid_log,*) '*** Number of Total .pexxxxxx files : ', MNG_PALL_in
  write(fid_log,*) '*** Number of PE to packing precess : ', prc_nall
  write(fid_log,*) '*** The rank of this process        : ', prc_myrank
  write(fid_log,*) '*** Number of files for this rank   : ', prc_nlocal
  write(fid_log,*) '*** file ID to pack                 : ', pstr-1, ' - ', pend-1
  write(fid_log,*)

  do p = pstr, pend
     write(fid_log,'(1x,A,I6.6)') '+pe= ', p-1

     call fio_mk_fname(fname,cstr(infile(1)),cstr('pe'),p-1,6)

     fid  = fio_register_file(fname)
     ierr = fio_fopen(fid,FIO_FREAD)
     ! put information from 1st input file
     ierr = fio_put_commoninfo_fromfile(fid,FIO_BIG_ENDIAN)

     ierr = fio_read_allinfo(fid)
     ierr = fio_get_pkginfo(hinfo,fid)

     if ( p == pstr ) then
        call fstr(pkg_desc, hinfo%description)
        call fstr(pkg_note, hinfo%note       )
        nmax_data = hinfo%num_of_data
        allocate( dinfo(0:nmax_data-1) )
     endif
     rgnid_in(:) = MNG_prc_tab_out(1:LALL_in,p)-1

     call fstr(fname)
     write(fid_log,*) '+input: ', trim(fname), ' (n=', nmax_data, ')'

     do did = 0, nmax_data-1
        ! get datainfo from input file
        ierr = fio_get_datainfo(dinfo(did),fid,did)
        call fstr(var_name_file, dinfo(did)%varname )

        write(fid_log,'(1x,A,I4.4,2A)') '+variable id= ', did, ', name: ', trim(var_name_file)

        KALL  = dinfo(did)%num_of_layer
        oprec = FIO_preclist(dinfo(did)%datatype)

        ! read PaNDa->write dump
        if ( dinfo(did)%datatype == FIO_REAL4 ) then

           allocate( data4_1D    (GALL_in *KALL*LALL_in ) )
           allocate( data4_3D_in (GALL_in ,KALL,LALL_in ) )
           allocate( data4_3D_out(GALL_out,KALL,ndivide*ndivide) )

           ierr = fio_read_data(fid,did,c_loc(data4_1D(:)))

           data4_3D_in = reshape( data4_1D,shape(data4_3D_in) )

           do l = 1, LALL_in
              write(fid_log,'(1x,A,I6.6,2ES11.3)') ' +in (rgnid,min,max)= ',   &
                                                   rgnid_in(l),                &
                                                   minval(data4_3D_in(:,:,l)), &
                                                   maxval(data4_3D_in(:,:,l))

              call devideregion_SP( glevel,              & ! [IN]
                                    rlevel_in,           & ! [IN]
                                    rlevel_out,          & ! [IN]
                                    ndivide,             & ! [IN]
                                    GALL_in,             & ! [IN]
                                    GALL_out,            & ! [IN]
                                    KALL,                & ! [IN]
                                    rgnid_in    (l),     & ! [IN]
                                    data4_3D_in (:,:,l), & ! [IN]
                                    rgnid_div   (:),     & ! [OUT]
                                    data4_3D_out(:,:,:)  ) ! [OUT]

              do n = 1, ndivide*ndivide
                 write(fid_log,'(1x,A,I6.6,2ES11.3)') ' +out(rgnid,min,max)= ',    &
                                                      rgnid_div(n),                &
                                                      minval(data4_3D_out(:,:,n)), &
                                                      maxval(data4_3D_out(:,:,n))

                 write(fname_dump,'(A4,I4.4,I6.6)') 'dump', did, rgnid_div(n)

                 fid_dump = IO_get_available_fid()
                 open( unit   = fid_dump,            &
                       file   = trim(fname_dump),    &
                       form   = 'unformatted',       &
                       access = 'direct',            &
                       recl   = GALL_out*KALL*oprec, &
                       status = 'new'                )
                     write(fid_dump,rec=1) data4_3D_out(:,:,n)
                 close(fid_dump)
              enddo
           enddo

           deallocate( data4_1D     )
           deallocate( data4_3D_in  )
           deallocate( data4_3D_out )

        elseif( dinfo(did)%datatype == FIO_REAL8 ) then

           allocate( data8_1D    (GALL_in *KALL*LALL_in ) )
           allocate( data8_3D_in (GALL_in ,KALL,LALL_in ) )
           allocate( data8_3D_out(GALL_out,KALL,ndivide*ndivide) )

           ierr = fio_read_data(fid,did,c_loc(data8_1D(:)))

           data8_3D_in = reshape( data8_1D,shape(data8_3D_in) )

           do l = 1, LALL_in
              write(fid_log,'(1x,A,I6.6,2ES11.3)') ' +in (rgnid,min,max)= ',   &
                                                   rgnid_in(l),                &
                                                   minval(data8_3D_in(:,:,l)), &
                                                   maxval(data8_3D_in(:,:,l))

              call devideregion_DP( glevel,              & ! [IN]
                                    rlevel_in,           & ! [IN]
                                    rlevel_out,          & ! [IN]
                                    ndivide,             & ! [IN]
                                    GALL_in,             & ! [IN]
                                    GALL_out,            & ! [IN]
                                    KALL,                & ! [IN]
                                    rgnid_in    (l),     & ! [IN]
                                    data8_3D_in (:,:,l), & ! [IN]
                                    rgnid_div   (:),     & ! [OUT]
                                    data8_3D_out(:,:,:)  ) ! [OUT]

              do n = 1, ndivide*ndivide
                 write(fid_log,'(1x,A,I6.6,2ES11.3)') ' +out(rgnid,min,max)= ',    &
                                                      rgnid_div(n),                &
                                                      minval(data8_3D_out(:,:,n)), &
                                                      maxval(data8_3D_out(:,:,n))

                 write(fname_dump,'(A4,I4.4,I6.6)') 'dump', did, rgnid_div(n)

                 fid_dump = IO_get_available_fid()
                 open( unit   = fid_dump,            &
                       file   = trim(fname_dump),    &
                       form   = 'unformatted',       &
                       access = 'direct',            &
                       recl   = GALL_out*KALL*oprec, &
                       status = 'new'                )
                     write(fid_dump,rec=1) data8_3D_out(:,:,n)
                 close(fid_dump)
              enddo
           enddo

           deallocate( data8_1D     )
           deallocate( data8_3D_in  )
           deallocate( data8_3D_out )

        endif
     enddo

     ierr = fio_fclose(fid)
  enddo

  prc_nlocal = MNG_PALL_out / prc_nall
  pstr       = prc_myrank*prc_nlocal + 1
  pend       = prc_myrank*prc_nlocal + prc_nlocal
  write(fid_log,*)
  write(fid_log,*) '*** For Output *** '
  write(fid_log,*) '*** Number of Total .pexxxxxx files : ', MNG_PALL_out
  write(fid_log,*) '*** Number of PE to packing precess : ', prc_nall
  write(fid_log,*) '*** The rank of this process        : ', prc_myrank
  write(fid_log,*) '*** Number of files for this rank   : ', prc_nlocal
  write(fid_log,*) '*** file ID to pack                 : ', pstr-1, ' - ', pend-1
  write(fid_log,*)

  do p = pstr, pend
!     write(fid_log,'(1x,A,I6.6)') '+pe= ', p-1

     rgnid_out(1:LALL_out) = MNG_prc_tab_out(1:LALL_out,p)-1

     ierr = fio_put_commoninfo( FIO_SPLIT_FILE,  &
                                FIO_BIG_ENDIAN,  &
                                FIO_ICOSAHEDRON, &
                                glevel,          &
                                rlevel_out,      &
                                LALL_out,        &
                                rgnid_out        )

     call fio_mk_fname(fname, trim(outfile),'pe',p-1,6)

     fid = fio_register_file(fname)
     ierr = fio_fopen(fid,FIO_FREAD)
     ierr = fio_put_write_pkginfo(fid,cstr(pkg_desc),cstr(pkg_note))

     write(fid_log,*) '+output : ', trim(fname), ' (n=', nmax_data, ')'

     do did = 0, nmax_data-1
        write(fid_log,'(1x,A,I4.4,2A)') '+variable id= ', did, ', name: ', trim(var_name_file)

        KALL  = dinfo(did)%num_of_layer
        oprec = FIO_preclist(dinfo(did)%datatype)

        if ( p == pstr ) then
           dinfo(did)%datasize = GALL_out * LALL_out * KALL * oprec
        endif

        odid = fio_put_write_datainfo(fid,dinfo(did))

        do l = 1, LALL_out
           write(fname_dump,'(A4,2I4.4,I6.6)') 'dump', did, odid, rgnid_out(l)

           ! read dump->write PaNDa
           if ( dinfo(did)%datatype == FIO_REAL4 ) then

              allocate( data4_1D(GALL_out*KALL) )

              open( unit   = fid_dump,            &
                    file   = trim(fname_dump),    &
                    form   = 'unformatted',       &
                    access = 'direct',            &
                    recl   = GALL_out*KALL*oprec, &
                    status = 'old'                )
                  read(fid_dump,rec=1) data4_1D
              close(fid_dump)

              ierr = fio_write_data_1rgn(fid,odid,c_loc(data4_1D))

              deallocate( data4_1D )

           elseif( dinfo(did)%datatype == FIO_REAL8 ) then

              allocate( data8_1D(GALL_out*KALL) )

              open( unit   = fid_dump,            &
                    file   = trim(fname_dump),    &
                    form   = 'unformatted',       &
                    access = 'direct',            &
                    recl   = GALL_out*KALL*oprec, &
                    status = 'old'                )
                  read(fid_dump,rec=1) data8_1D
              close(fid_dump)

              ierr = fio_write_data_1rgn(fid,odid,c_loc(data8_1D))

              deallocate( data8_1D )
           endif
        enddo ! l-loop

     enddo ! data loop

     ierr = fio_fclose(fid)
  enddo ! PE loop

  write(fid_log,*) 'convert finished.'

  if ( use_mpi ) then
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
  endif

  close(fid_log)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> read option
  subroutine readoption
    use mod_tool_option, only: &
      OPT_convert, &
      OPT_fid
    implicit none

    integer  :: io
    !---------------------------------------------------------------------------

    ! --- Set option
    OPT_fid = IO_get_available_fid()
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
  subroutine helpoption
    implicit none
    !---------------------------------------------------------------------------

    write(*,OPTION)

    stop
  end subroutine helpoption

  !-----------------------------------------------------------------------------
  !> region devider
  subroutine devideregion_SP( &
       glevel,     &
       rlevel_in,  &
       rlevel_out, &
       ndivide,    &
       GALL_in,    &
       GALL_out,   &
       KALL,       &
       rgnid_in,   &
       data_in,    &
       rgnid_div,  &
       data_out    )
    implicit none

    integer,  intent(in)  :: glevel
    integer,  intent(in)  :: rlevel_in
    integer,  intent(in)  :: rlevel_out
    integer,  intent(in)  :: ndivide
    integer,  intent(in)  :: GALL_in
    integer,  intent(in)  :: GALL_out
    integer,  intent(in)  :: KALL
    integer,  intent(in)  :: rgnid_in
    real(SP), intent(in)  :: data_in  (GALL_in,KALL)
    integer,  intent(out) :: rgnid_div(ndivide*ndivide)
    real(SP), intent(out) :: data_out (GALL_out,KALL,ndivide*ndivide)

    integer :: grid_1d_in, grid_1d_out
    integer :: rgn_1d_in, rgn_1d_out
    integer :: dmdid
    integer :: ri_in, rj_in, rij_in
    integer :: ri_out, rj_out, rij_out
    integer :: gi_in, gj_in, gij_in
    integer :: gi_out, gj_out, gij_out

    integer :: i, j, n
    !---------------------------------------------------------------------------

    grid_1d_in  = 2**(glevel-rlevel_in ) + 2
    grid_1d_out = 2**(glevel-rlevel_out) + 2

    rgn_1d_in   = 2**rlevel_in
    rgn_1d_out  = 2**rlevel_out

    dmdid  =     rgnid_in/rgn_1d_in**2
    rij_in = mod(rgnid_in,rgn_1d_in**2)
    rj_in  =     rij_in/rgn_1d_in
    ri_in  = mod(rij_in,rgn_1d_in)

    do j = 1, ndivide
    do i = 1, ndivide
       n = (j-1)*ndivide + i

       rj_out  = rj_in*ndivide
       ri_out  = ri_in*ndivide
       rij_out = rj_out*rgn_1d_out + ri_out &
               + (j-1) *rgn_1d_out + (i-1)

       rgnid_div(n) = dmdid*rgn_1d_out*rgn_1d_out + rij_out

       do gj_out = 1, grid_1d_out
       do gi_out = 1, grid_1d_out
          gj_in = (j-1)*(grid_1d_out-2) + gj_out
          gi_in = (i-1)*(grid_1d_out-2) + gi_out

          gij_in  = (gj_in -1)*grid_1d_in  + gi_in
          gij_out = (gj_out-1)*grid_1d_out + gi_out

          data_out(gij_out,:,n) = data_in(gij_in,:)
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine devideregion_SP

  !-----------------------------------------------------------------------------
  !> region devider
  subroutine devideregion_DP( &
       glevel,     &
       rlevel_in,  &
       rlevel_out, &
       ndivide,    &
       GALL_in,    &
       GALL_out,   &
       KALL,       &
       rgnid_in,   &
       data_in,    &
       rgnid_div,  &
       data_out    )
    implicit none

    integer,  intent(in)  :: glevel
    integer,  intent(in)  :: rlevel_in
    integer,  intent(in)  :: rlevel_out
    integer,  intent(in)  :: ndivide
    integer,  intent(in)  :: GALL_in
    integer,  intent(in)  :: GALL_out
    integer,  intent(in)  :: KALL
    integer,  intent(in)  :: rgnid_in
    real(DP), intent(in)  :: data_in  (GALL_in,KALL)
    integer,  intent(out) :: rgnid_div(ndivide*ndivide)
    real(DP), intent(out) :: data_out (GALL_out,KALL,ndivide*ndivide)

    integer :: grid_1d_in, grid_1d_out
    integer :: rgn_1d_in, rgn_1d_out
    integer :: dmdid
    integer :: ri_in, rj_in, rij_in
    integer :: ri_out, rj_out, rij_out
    integer :: gi_in, gj_in, gij_in
    integer :: gi_out, gj_out, gij_out

    integer :: i, j, n
    !---------------------------------------------------------------------------

    grid_1d_in  = 2**(glevel-rlevel_in ) + 2
    grid_1d_out = 2**(glevel-rlevel_out) + 2

    rgn_1d_in   = 2**rlevel_in
    rgn_1d_out  = 2**rlevel_out

    dmdid  =     rgnid_in/rgn_1d_in**2
    rij_in = mod(rgnid_in,rgn_1d_in**2)
    rj_in  =     rij_in/rgn_1d_in
    ri_in  = mod(rij_in,rgn_1d_in)

    do j = 1, ndivide
    do i = 1, ndivide
       n = (j-1)*ndivide + i

       rj_out  = rj_in*ndivide
       ri_out  = ri_in*ndivide
       rij_out = rj_out*rgn_1d_out + ri_out &
               + (j-1) *rgn_1d_out + (i-1)

       rgnid_div(n) = dmdid*rgn_1d_out*rgn_1d_out + rij_out

       do gj_out = 1, grid_1d_out
       do gi_out = 1, grid_1d_out
          gj_in = (j-1)*(grid_1d_out-2) + gj_out
          gi_in = (i-1)*(grid_1d_out-2) + gi_out

          gij_in  = (gj_in -1)*grid_1d_in  + gi_in
          gij_out = (gj_out-1)*grid_1d_out + gi_out

          data_out(gij_out,:,n) = data_in(gij_in,:)
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine devideregion_DP

end program fio_pe2pe
