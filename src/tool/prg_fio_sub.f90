!-------------------------------------------------------------------------------
!> Program FIO sub
!!
!! @par Description
!!          substitute pe0000x format data
!!          ( packaged NICAM data format : PaNDa )
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
program fio_sub
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
     MNG_mnginfo_input,   &
     MNG_PALL,            &
     MNG_prc_rnum
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
  integer, parameter :: max_nvar = 1000

  integer, parameter :: flim = 2
  integer,      save :: fmax

  !--- NAMELIST
  integer               :: glevel       = -1
  integer               :: rlevel       = -1
  character(len=H_LONG) :: mnginfo      = ''
  character(len=H_LONG) :: infile(flim) = ''
  character(len=H_LONG) :: outfile      = ''
  logical               :: use_mpi      = .true.
  integer               :: pe_str       =  0
  integer               :: pe_end       = -1
  logical               :: help         = .false.

  namelist /OPTION/ glevel,  &
                    rlevel,  &
                    mnginfo, &
                    infile,  &
                    outfile, &
                    use_mpi, &
                    pe_str,  &
                    pe_end,  &
                    help

  !-----------------------------------------------------------------------------
  character(len=H_LONG)  :: infname1 = ''
  character(len=H_LONG)  :: infname2 = ''
  character(len=H_LONG)  :: outfname = ''

  type(headerinfo_panda) :: hinfo
  type(datainfo_panda)   :: dinfo1
  type(datainfo_panda)   :: dinfo2

  character(len=H_MID)   :: pkg_desc
  character(len=H_LONG)  :: pkg_note
  integer                :: nmax_data

  integer                :: nvar
  character(len=H_SHORT) :: var_name (max_nvar)
  character(len=H_SHORT) :: var_name_file
  integer                :: var_nstep(max_nvar)

  integer                :: GALL
  integer                :: KALL
  integer                :: LALL
  real(SP), target, allocatable :: data4_1D1(:)
  real(DP), target, allocatable :: data8_1D1(:)
  real(SP), target, allocatable :: data4_1D2(:)
  real(DP), target, allocatable :: data8_1D2(:)

  ! for MPI
  integer                :: pe_all
  integer                :: prc_nall, prc_nlocal
  integer                :: prc_myrank, pstr, pend
  integer                :: fid_log
  character(len=6)       :: rankstr

  logical  :: addvar
  integer  :: p, v, vid
  integer  :: ifid1, ifid2, idid, ofid, odid, ierr
  !=====================================================================

  !--- read option and preprocess
  call readoption !! set fmax, infile

  !--- prepare region infomation
  call MNG_mnginfo_input( rlevel, trim(mnginfo) )

  GALL = ( (2**(glevel-rlevel))+2 ) &
       * ( (2**(glevel-rlevel))+2 )

  fid_log = IO_get_available_fid()
  if ( use_mpi ) then
     !--- Parallel Excution, No communication
     call MPI_Init(ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, prc_nall,   ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD, prc_myrank, ierr)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     write(rankstr,'(I6.6)') prc_myrank
     open(fid_log, file='msg_sub.pe'//trim(rankstr) )
     write(fid_log,*) '+++ Parallel Execution, Use MPI'

     if( mod( MNG_PALL, prc_nall) /= 0)then
        write(fid_log,*) '*** Invalid processor number, STOP:', MNG_PALL, prc_nall
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

  if ( pe_end >= 0 ) then
     pe_all = pe_end - pe_str + 1
     write(fid_log,*) '*** pe range is specified. '
     write(fid_log,*) '*** pe(all,start,end)=',pe_all,pe_str,pe_end
  else
     pe_all = MNG_PALL
  endif

  if ( mod( pe_all, prc_nall) /= 0 ) then
     write(fid_log,*) '*** Invalid processor number, STOP:', pe_all, prc_nall
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
     stop
  endif

  prc_nlocal = pe_all / prc_nall
  pstr       = prc_myrank*prc_nlocal + pe_str + 1
  pend       = prc_myrank*prc_nlocal + pe_str + prc_nlocal
  write(fid_log,*) '*** Number of Total .pexxxxxx files: ', MNG_PALL
  write(fid_log,*) '*** Number of PE to packing precess: ', prc_nall
  write(fid_log,*) '*** The rank of this process       : ', prc_myrank
  write(fid_log,*) '*** Number of files for this rank  : ', prc_nlocal
  write(fid_log,*) '*** file ID to pack                : ', pstr-1, ' - ', pend-1

  !--- setup
  ierr = fio_syscheck()

  write(fid_log,*) '*** combine start : PaNDa format to PaNDa format data'

  do p = pstr, pend
     write(fid_log,*) '+pe:', p-1
     LALL = MNG_prc_rnum(p)

     call fio_mk_fname(infname1,cstr(infile(1)),cstr('pe'),p-1,6)
     call fio_mk_fname(infname2,cstr(infile(2)),cstr('pe'),p-1,6)
     call fio_mk_fname(outfname,cstr(outfile)  ,cstr('pe'),p-1,6)
     call fstr(outfname)
     write(fid_log,*) '++output : ', trim(outfname)

     ifid1 = fio_register_file(infname1)
     ierr  = fio_fopen(ifid1,FIO_FREAD)
     ! put information from 1st input file
     ierr  = fio_put_commoninfo_fromfile(ifid1,FIO_BIG_ENDIAN)

     ierr  = fio_read_allinfo(ifid1)
     ierr  = fio_get_pkginfo(hinfo,ifid1)
     call fstr(pkg_desc, hinfo%description)
     call fstr(pkg_note, hinfo%note       )
     nmax_data = hinfo%num_of_data
     call fstr(infname1)
     write(fid_log,*) '++input', 1, ' : ', trim(infname1), '(n=', hinfo%num_of_data, ')'

     ifid2 = fio_register_file(infname2)
     ierr  = fio_fopen(ifid2,FIO_FREAD)
     ierr  = fio_read_allinfo(ifid2)
     ierr  = fio_get_pkginfo(hinfo,ifid2)
     call fstr(infname2)
     write(fid_log,*) '++input', 2, ' : ', trim(infname2), '(n=', hinfo%num_of_data, ')'

     if ( hinfo%num_of_data /= nmax_data ) then
        write(fid_log,*) '*** Mismatch number of data, STOP:', hinfo%num_of_data, nmax_data
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
     endif

     ofid = fio_register_file(cstr(outfname))
     ierr = fio_fopen(ofid,FIO_FWRITE)
     ierr = fio_put_write_pkginfo(ofid,cstr(pkg_desc),cstr(pkg_note))

     nvar = 0
     do idid = 0, nmax_data-1
        ! get datainfo from input file
        ierr = fio_get_datainfo(dinfo1,ifid1,idid)
        ierr = fio_get_datainfo(dinfo2,ifid2,idid)
        KALL = dinfo1%num_of_layer
        call fstr(var_name_file, dinfo1%varname)

        if ( dinfo2%num_of_layer /= KALL ) then
           write(fid_log,*) '*** Mismatch number of layer, STOP:', dinfo2%num_of_layer, KALL
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           call MPI_FINALIZE(ierr)
           stop
        endif

        addvar = .true.
        do v = 1, nvar
           if ( var_name(v) == var_name_file ) then
              var_nstep(v) = var_nstep(v) + 1
              vid = v
              addvar = .false.
              exit
           endif
        enddo

        if (addvar) then
           nvar = nvar + 1
           var_nstep(nvar) = 1
           var_name (nvar) = var_name_file
           vid = nvar
        endif

        dinfo1%step = var_nstep(vid)

        ! read->write data
        if ( dinfo1%datatype == FIO_REAL4 ) then
           allocate( data4_1D1(GALL*KALL*LALL) )
           ierr = fio_read_data(ifid1,idid,c_loc(data4_1D1(:)))
           if ( dinfo2%datatype == FIO_REAL4 ) then
              allocate( data4_1D2(GALL*KALL*LALL) )
              ierr = fio_read_data(ifid2,idid,c_loc(data4_1D2(:)))
              data4_1D1(:) = data4_1D1(:) - data4_1D2(:)
              deallocate( data4_1D2 )
           elseif( dinfo2%datatype == FIO_REAL8 ) then
              allocate( data8_1D2(GALL*KALL*LALL) )
              ierr = fio_read_data(ifid2,idid,c_loc(data8_1D2(:)))
              data4_1D1(:) = data4_1D1(:) - real(data8_1D2(:),kind=4)
              deallocate( data8_1D2 )
           endif
           odid = fio_put_write_datainfo_data(ofid,dinfo1,c_loc(data4_1D1(:)))
           deallocate( data4_1D1 )
        elseif( dinfo1%datatype == FIO_REAL8 ) then
           allocate( data8_1D1(GALL*KALL*LALL) )
           ierr = fio_read_data(ifid1,idid,c_loc(data8_1D1(:)))
           if ( dinfo2%datatype == FIO_REAL4 ) then
              allocate( data4_1D2(GALL*KALL*LALL) )
              ierr = fio_read_data(ifid2,idid,c_loc(data4_1D2(:)))
              data8_1D1(:) = data8_1D1(:) - real(data4_1D2(:),kind=8)
              deallocate( data4_1D2 )
           elseif( dinfo2%datatype == FIO_REAL8 ) then
              allocate( data8_1D2(GALL*KALL*LALL) )
              ierr = fio_read_data(ifid2,idid,c_loc(data8_1D2(:)))
              data8_1D1(:) = data8_1D1(:) - data8_1D2(:)
              deallocate( data8_1D2 )
           endif
           odid = fio_put_write_datainfo_data(ofid,dinfo1,c_loc(data8_1D1(:)))
           deallocate( data8_1D1 )
        endif

     enddo

     ierr = fio_fclose(ifid1)
     ierr = fio_fclose(ifid2)
     ierr = fio_fclose(ofid)
  enddo ! PE loop

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

end program fio_sub
