!-------------------------------------------------------------------------------
!
! interface description between C and Fortran
!
!-------------------------------------------------------------------------------

!> struct for data infomation
type, BIND(C) :: datainfo_panda
   character(c_char) :: varname    (FIO_HSHORT) !< variable name
   character(c_char) :: description(FIO_HMID)   !< variable description
   character(c_char) :: unit       (FIO_HSHORT) !< unit of variable
   character(c_char) :: layername  (FIO_HSHORT) !< layer name
   character(c_char) :: note       (FIO_HLONG)  !< longer note of variable
   integer(8)        :: datasize                !< data size
   integer           :: datatype                !< data type(0,1,2,3)
   integer           :: num_of_layer            !< number of layer
   integer           :: step
   integer(8)        :: time_start
   integer(8)        :: time_end
endtype datainfo_panda

!> struct for package infomation
type, BIND(C) :: headerinfo_panda
   character(c_char) :: fname      (FIO_HLONG) !< file name
   character(c_char) :: description(FIO_HMID)  !< variable description
   character(c_char) :: note       (FIO_HLONG) !< longer note of variable
   integer           :: num_of_data            !< number of data
   integer           :: fmode                  !< file mode(0,1,2)
   integer           :: endiantype             !< endian type(0,1,2)
   integer           :: grid_topology          !< grid topology(0,1,2)
   integer           :: glevel                 !< glevel
   integer           :: rlevel                 !< rlevel
   integer           :: num_of_rgn             !< number of region
   type(c_ptr)       :: rgnid                  !< array of region id
endtype headerinfo_panda

interface

   subroutine fio_mk_fname(fname, base, ext, i, y) BIND(C)
     import c_int, c_char
     character(c_char), intent(in)        :: fname(*)
     character(c_char), intent(in)        :: base(*)
     character(c_char), intent(in)        :: ext(*)
     integer(c_int),    intent(in), value :: i
     integer(c_int),    intent(in), value :: y
   end subroutine fio_mk_fname

   function fio_syscheck() result(ierr) BIND(C)
     import c_int
     integer(c_int) :: ierr
   end function fio_syscheck

   function fio_put_commoninfo(fmode, endiantype, grid_topology, glevel, rlevel, num_of_rgn, rgnid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fmode
     integer(c_int), intent(in), value :: endiantype
     integer(c_int), intent(in), value :: grid_topology
     integer(c_int), intent(in), value :: glevel
     integer(c_int), intent(in), value :: rlevel
     integer(c_int), intent(in), value :: num_of_rgn
     integer(c_int), intent(in)        :: rgnid(*)
     integer(c_int) :: ierr
   end function fio_put_commoninfo

   function fio_put_commoninfo_fromfile(fid, endiantype) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: endiantype
     integer(c_int) :: ierr
   end function fio_put_commoninfo_fromfile

   function fio_get_pkginfo(hinfo, fid) result(ierr) BIND(C)
     import c_int
     import headerinfo_panda
     type(headerinfo_panda), intent(out)     :: hinfo
     integer(c_int),         intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_get_pkginfo

   function fio_get_datainfo(ditem, fid, did) result(ierr) BIND(C)
     import c_int
     import datainfo_panda
     type(datainfo_panda), intent(out)       :: ditem
     integer(c_int),       intent(in), value :: fid
     integer(c_int),       intent(in), value :: did
     integer(c_int) :: ierr
   end function fio_get_datainfo

   function fio_seek_datainfo(fid, varname, step) result(did) BIND(C)
     import c_int, c_char
     integer(c_int),    intent(in), value :: fid
     character(c_char), intent(in)        :: varname(*)
     integer(c_int),    intent(in), value :: step
     integer(c_int) :: did
   end function fio_seek_datainfo

   function fio_fopen(fid, mode) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: mode
     integer(c_int) :: ierr
   end function fio_fopen

   function fio_fclose(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_fclose

   function fio_write_pkginfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_write_pkginfo

   function fio_read_pkginfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_read_pkginfo

   function fio_write_data_1rgn(fid, did, data) result(ierr) BIND(C)
     import c_int, c_ptr
     import datainfo_panda
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: did
     type(c_ptr),    intent(in), value :: data
     integer(c_int) :: ierr
   end function fio_write_data_1rgn

   function fio_read_data(fid, did, data) result(ierr) BIND(C)
     import c_int, c_ptr
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: did
     type(c_ptr),    intent(in), value :: data
     integer(c_int) :: ierr
   end function fio_read_data

   function fio_read_data_tmpdata(fid, did, data) result(ierr) BIND(C)
     import c_int, c_ptr
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: did
     type(c_ptr),    intent(in), value :: data
     integer(c_int) :: ierr
   end function fio_read_data_tmpdata

   function fio_read_data_1layer(fid, did, k, kmax, l, lmax, data) result(ierr) BIND(C)
     import c_int, c_ptr
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: did
     integer(c_int), intent(in), value :: k
     integer(c_int), intent(in), value :: kmax
     integer(c_int), intent(in), value :: l
     integer(c_int), intent(in), value :: lmax
     type(c_ptr),    intent(in), value :: data
     integer(c_int) :: ierr
   end function fio_read_data_1layer

   function fio_register_file(fname) result(fid) BIND(C)
     import c_int, c_char
     character(c_char), intent(in) :: fname(*)
     integer(c_int) :: fid
   end function fio_register_file

   function fio_put_write_pkginfo(fid, description, note) result(ierr) BIND(C)
     import c_int, c_char
     integer(c_int),    intent(in), value :: fid
     character(c_char), intent(in)        :: description(*)
     character(c_char), intent(in)        :: note(*)
     integer(c_int) :: ierr
   end function fio_put_write_pkginfo

   function fio_valid_pkginfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_valid_pkginfo

   function fio_valid_pkginfo_validrgn(fid, rgnid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in)        :: rgnid(*)
     integer(c_int) :: ierr
   end function fio_valid_pkginfo_validrgn

   function fio_put_write_datainfo_data(fid, ditem, data) result(did) BIND(C)
     import c_int, c_ptr
     import datainfo_panda
     integer(c_int),       intent(in), value :: fid
     type(datainfo_panda), intent(in)        :: ditem
     type(c_ptr),          intent(in), value :: data
     integer(c_int) :: did
   end function fio_put_write_datainfo_data

   function fio_put_write_datainfo(fid, ditem) result(did) BIND(C)
     import c_int, c_ptr
     import datainfo_panda
     integer(c_int),       intent(in), value :: fid
     type(datainfo_panda), intent(in)        :: ditem
     integer(c_int) :: did
   end function fio_put_write_datainfo

   function fio_read_allinfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_read_allinfo

   function fio_read_allinfo_validrgn(fid, rgnid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in)        :: rgnid(*)
     integer(c_int) :: ierr
   end function fio_read_allinfo_validrgn

   function fio_read_allinfo_tmpdata(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int) :: ierr
   end function fio_read_allinfo_tmpdata

   function fio_register_vname_tmpdata(varname) result(ierr) BIND(C)
     import c_int, c_char
     character(c_char), intent(in) :: varname(*)
     integer(c_int) :: ierr
   end function fio_register_vname_tmpdata

   function fio_copy_datainfo(fid, fid_org) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: fid_org
     integer(c_int) :: ierr
   end function fio_copy_datainfo

   function fio_dump_finfo(fid, endiantype, dumptype) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: fid
     integer(c_int), intent(in), value :: endiantype
     integer(c_int), intent(in), value :: dumptype
     integer(c_int) :: ierr
   end function fio_dump_finfo

   subroutine free(ptr) BIND(C)
     import c_ptr
     type(c_ptr), intent(in), value :: ptr
   end subroutine free

   function orgico_setup(gl, rl) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), value :: gl
     integer(c_int), intent(in), value :: rl
     integer(c_int) :: ierr
   end function orgico_setup

   subroutine orgico_mk_fname(fname, base, i, y) BIND(C)
     import c_int, c_char
     character(c_char), intent(in)        :: fname(*)
     character(c_char), intent(in)        :: base(*)
     integer(c_int),    intent(in), value :: i
     integer(c_int),    intent(in), value :: y
   end subroutine orgico_mk_fname

   function orgico_setbasename(basename) result(ierr) BIND(C)
     import c_int, c_char
     character(c_char), intent(in) :: basename(*)
     integer(c_int) :: ierr
   end function orgico_setbasename

   function orgico_readgriddata(fname, data, i) result(ierr) BIND(C)
     import c_int, c_char, c_ptr
     character(c_char), intent(in)        :: fname(*)
     type(c_ptr),       intent(in), value :: data
     integer(c_int),    intent(in), value :: i
     integer(c_int) :: ierr
   end function orgico_readgriddata

   function orgico_readdata_seq(fname, did, nol, data, esize) result(ierr) BIND(C)
     import c_int, c_char, c_ptr
     character(c_char), intent(in)        :: fname(*)
     integer(c_int),    intent(in), value :: did
     type(c_ptr),       intent(in), value :: data
     integer(c_int),    intent(in), value :: nol
     integer(c_int),    intent(in), value :: esize
     integer(c_int) :: ierr
   end function orgico_readdata_seq

   function orgico_readdata_dir(fname, did, nol, data, esize) result(ierr) BIND(C)
     import c_int, c_char, c_ptr
     character(c_char), intent(in)        :: fname(*)
     integer(c_int),    intent(in), value :: did
     type(c_ptr),       intent(in), value :: data
     integer(c_int),    intent(in), value :: nol
     integer(c_int),    intent(in), value :: esize
     integer(c_int) :: ierr
   end function orgico_readdata_dir

end interface
