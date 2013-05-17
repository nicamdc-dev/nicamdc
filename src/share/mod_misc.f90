!-------------------------------------------------------------------------------
!>
!! miscellaneous module
!!
!! @par Description
!!         This module contains miscellaneous subroutines.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2005-11-15 (M.Satoh)   [mod] MISC_make_idstr
!! @li      2006-02-12 (S.Iga)     add critical value in MISC_triangle_area for the case angle=0
!!                                 integer(4) -> integer  for 'second'
!! @li      2006-02-25 (S.Iga)     bugfix on 06-02-12
!! @li      2006-08-10 (W.Yanase)  bug(?)fix on 06-08-10
!! @li      2006-08-22 (Y.Niwa)    bug fix
!! @li      2007-01-26 (H.Tomita)  Add func[MISC_gammafunc].
!! @li      2008-04-12 (T.Mitsui)  Add func[MISC_igammafunc].
!!                                 Add func[MISC_iigammafunc].
!! @li      2009-07-17 (Y.Yamada)  Add func[MISC_triangle_area_q].
!! @li      2011-11-11 (H.Yashiro) [Add] vector calculation suite
!!
!<
module mod_misc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters, variables & subroutines
  !
  public :: MISC_make_idstr        !--- make file name with a number
  public :: MISC_in_region         !--- judge whether inner region or not
  public :: MISC_frndm             !--- random generator (import from dcl)
  public :: MISC_get_available_fid !--- get an available file ID
  public :: MISC_get_fid           !--- get an information of open or not
  public :: MISC_get_latlon        !--- calculate (lat,lon) from (x,y,z)
  public :: MISC_triangle_area     !--- calculate triangle area.
  public :: MISC_triangle_area_q   !--- calculate triangle area at quad precision.  ![Add] Y.Yamada 09/07/17
  public :: MISC_mk_gmtrvec        !--- calculate normal and tangential vecs.
  public :: MISC_sec2ymdhms        !--- convert second to (yr,mo,dy,hr,mn,sc).
  public :: MISC_ymdhms2sec        !--- convert (yr,mo,dy,hr,mn,sc) to second.
  public :: MISC_msg_nmerror       !--- output error message
  public :: MISC_gammafunc         !--- Gamma function
  public :: MISC_igammafunc        !--- Incomplete Gamma function       ![Add] T.Mitsui 08/04/12
  public :: MISC_iigammafunc       !--- inverse function of igammafunc  ![Add] T.Mitsui 08/04/12
  ![Add] H.Yashiro 11/11/11
  public :: MISC_3dvec_cross         !--- calc exterior product for 3D vector
  public :: MISC_3dvec_dot           !--- calc interior product for 3D vector
  public :: MISC_3dvec_abs           !--- calc absolute value for 3D vector
  public :: MISC_3dvec_angle         !--- calc angle for 3D vector
  public :: MISC_3dvec_intersec      !--- calc intersection of two 3D vectors
  public :: MISC_3dvec_anticlockwise !--- sort 3D vectors anticlockwise
  public :: MISC_3Dvec_triangle      !--- calc triangle area on sphere, more precise
  public :: MISC_get_cartesian       !--- calc (x,y,z) from (lat,lon)
  public :: MISC_get_distance        !--- calc horizontal distance on the sphere

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters, variables & subroutines
  !
  !--- minimum available fid
  integer, parameter, private :: min_fid = 7
  !
  !--- maximum available fid.
  integer, parameter, private :: max_fid = 99

  !--- number of separated file starts from 0 ?
  logical, private, parameter :: NSTR_ZERO_START = .true.

  !--- digit of separated file
  integer, private, save      :: NSTR_MAX_DIGIT  = 5
  !
  !<----  if running on ES, 
  !<----      NSTR_ZERO_START = .TRUE.
  !<----      NSTR_MAX_DIGIT  = 5
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_make_idstr
  !>
  subroutine MISC_make_idstr( &
       str,     & !--- [OUT]
       prefix,  & !--- [IN]
       ext,     & !--- [IN]
       numID,   & !--- [IN]
       digit    ) !--- [IN]
    implicit none

    character(len=*), intent(out) :: str    ! combined string (file name)
    character(len=*), intent(in)  :: prefix ! prefix
    character(len=*), intent(in)  :: ext    ! extention( e.g. .rgn ) 
    integer,          intent(in)  :: numID  ! number

    integer, optional, intent(in) :: digit  ! digit

    character(len=128) :: rankstr
    integer            :: setdigit
    !---------------------------------------------------------------------------

    if ( NSTR_ZERO_START ) then
       write(rankstr,'(I128.128)') numID-1
    else
       write(rankstr,'(I128.128)') numID
    endif

    if ( present(digit) ) then
       setdigit = digit
    else
       setdigit = NSTR_MAX_DIGIT
    endif

    rankstr(1:setdigit) = rankstr(128-(setdigit-1):128)
    rankstr(setdigit+1:128) = ' '

    str = trim(prefix)//'.'//trim(ext)//trim(rankstr) ! -> prefix.ext00000

    return
  end subroutine MISC_make_idstr

  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return f
  !>
  function MISC_in_region( &
       clat,  &
       clon,  &
       alpha, &
       tlat,  &
       tlon   ) &
       result(f)
    implicit none

    real(8), intent(in) :: clat  ! center of circle ( latitude )
    real(8), intent(in) :: clon  ! center of circle ( longitude )
    real(8), intent(in) :: alpha ! search radius
    real(8), intent(in) :: tlat  ! target point ( latitude )
    real(8), intent(in) :: tlon  ! target point ( longitude )
    logical :: f                 ! OUT : judgement

    real(8) :: cp(3)
    real(8) :: tp(3)
    real(8) :: cdott
    !---------------------------------------------------------------------------

    cp(1) = cos(clat) * cos(clon)
    cp(2) = cos(clat) * sin(clon)
    cp(3) = sin(clat)

    tp(1) = cos(tlat) * cos(tlon)
    tp(2) = cos(tlat) * sin(tlon)
    tp(3) = sin(tlat)

    cdott = cp(1)*tp(1) + cp(2)*tp(2) + cp(3)*tp(3) ! inner production
    cdott = min( max( cdott, -1.0D0 ), 1.D0 )

    f = .false.
    if( acos(cdott) <= alpha ) f = .true.

    return
  end function MISC_in_region

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_frndm
  !>
  subroutine MISC_frndm( &
       ndim,   &
       iseed0, &
       random  )
    implicit none

    integer, intent(in)    :: ndim         ! number of array (1D)
    integer, intent(in)    :: iseed0       ! initial seed for randomization
    real(8), intent(inout) :: random(ndim) ! random values

    integer :: iseed, i
    !---------------------------------------------------------------------------

    iseed = iseed0
    do i = 1, ndim
       random(i) = rngu3(iseed) - 0.5D0
    enddo

    return
  contains
    !---------------------------------------------------------------------------
    !>
    !> RANDOM NUMBER GENERATOR (? METHOD BY KUNUTH) from dcl4.2
    !> @return rngu3
    !>
    function rngu3(iseed)
      integer :: iseed
      real(8) :: rngu3

      integer, parameter :: mbig  = 100000000
      integer, parameter :: mseed = 161803398
      integer, parameter :: mz    = 0
      real(8), parameter :: fac   = 1.D-9

      integer, save      :: ma(55)
      integer, save      :: inext, inextp

      integer :: mj, mk, i, ii, k
      !-------------------------------------------------------------------------

      if ( iseed /= 0 ) then
         mj = mseed - iabs(iseed)
         mj = mod(mj,mbig)
         ma(55) = mj
         mk = 1

         do i = 1,54
            ii = mod(21*i,55)
            ma(ii) = mk
            mk = mj-mk
            if(mk < mz) mk = mk+mbig
            mj = ma(ii)
         enddo

         do k = 1, 4
         do i = 1, 55
            ma(i) = ma(i) - ma(1+mod(i+30,55))
            if(ma(i) < mz)ma(i)=ma(i)+mbig
         enddo
         enddo

         inext  = 0
         inextp = 31
         iseed  = 0
      endif

      inext = inext + 1
      if(inext == 56) inext = 1

      inextp = inextp + 1
      if(inextp == 56) inextp = 1

      mj = ma(inext) - ma(inextp)
      if(mj < mz) mj = mj + mbig

      ma(inext) = mj

      rngu3 = mj * fac

      return
    end function rngu3

  end subroutine MISC_frndm

  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return
  !>
  function MISC_get_available_fid()  &
       result(fid)                     !--- file id
    !
    implicit none
    !
    integer :: fid
    logical :: i_opened
    !
    do fid = min_fid,max_fid
       INQUIRE (fid, OPENED=I_OPENED)
       if(.not.I_OPENED) return
    enddo
    !
  end function MISC_get_available_fid

  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return
  !>
  function MISC_get_fid( fname )  &
       result(fid)                   !--- file ID
    !
    implicit none
    !
    character(*), intent(in) :: fname
    !
    integer :: fid
    logical :: i_opened
    !
    INQUIRE (FILE=trim(fname), OPENED=i_opened, NUMBER=fid)
    if(.not.i_opened) fid = -1
    !
  end function MISC_get_fid

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_get_latlon
  !>
  subroutine MISC_get_latlon( &
       lat, lon,              & !--- INOUT : latitude and longitude
       x, y, z )                !--- IN : Cartesian coordinate ( on the sphere )
    !
    implicit none
    !
    real(8),intent(inout) :: lat, lon
    real(8),intent(in) :: x,y,z
    !
    real(8), parameter :: epsilon = 1.0D-99
    real(8) :: leng,leng_xy
    !
    leng=sqrt(x*x+y*y+z*z)
    !
    ! --- vector length equals to zero.
    if(leng<epsilon) then
       lat=0.0D0
       lon=0.0D0
       return
    endif
    ! --- vector is parallele to z axis.
    if(z/leng>=1.0D0) then
       lat=asin(1.0D0)
       lon=0.0D0
       return
    elseif(z/leng<=-1.0D0) then
       lat=asin(-1.0D0)
       lon=0.0D0
       return
    endif
    ! --- not parallele to z axis
    lat=asin(z/leng)
    !
    leng_xy=sqrt(x*x+y*y)
    if(leng_xy<epsilon) then
       lon=0.0D0
       return
    endif
    if(x/leng_xy>=1.0D0) then
       lon=acos(1.0D0)
       if(y<0.0D0) lon=-lon
       return
    elseif(x/leng_xy<=-1.0D0) then
       lon=acos(-1.0D0)
       if(y<0.0D0) lon=-lon
       return
    endif
    lon=acos(x/leng_xy)
    if(y<0.0D0) lon=-lon
    return
  end subroutine MISC_get_latlon
  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return
  !>
  function MISC_triangle_area( &
       a, b, c,                & !--- IN : three points vectors on a sphere.
       polygon_type,           & !--- IN : sphere triangle or plane one?
       radius,                 & !--- IN : radius
       critical)               & !--- IN : critical value to handle the case
                                 !         of ang=0 (optional) S.Iga060209
       result( area )            !--- OUT : triangle area
    !
    implicit none
    !
    integer, parameter :: ix = 1
    integer, parameter :: iy = 2
    integer, parameter :: iz = 3
    real(8), parameter :: pi  = 3.14159265358979323846D0
    !
    real(8) :: area
    real(8),intent(in) :: a(ix:iz),b(ix:iz),c(ix:iz)
    character(len=*), intent(in) :: polygon_type
    real(8) :: radius
    !
    !
    real(8) :: v01(ix:iz)
    real(8) :: v02(ix:iz)
    real(8) :: v03(ix:iz)
    !
    real(8) :: v11(ix:iz)
    real(8) :: v12(ix:iz)
    real(8) :: v13(ix:iz)
    !
    real(8) :: v21(ix:iz)
    real(8) :: v22(ix:iz)
    real(8) :: v23(ix:iz)
    real(8) :: w21(ix:iz)
    real(8) :: w22(ix:iz)
    real(8) :: w23(ix:iz)
    real(8) :: w11(ix:iz)
    real(8) :: w12(ix:iz)
    real(8) :: w13(ix:iz)

    real(8) :: v1(ix:iz)
    real(8) :: v2(ix:iz)
    real(8) :: w(ix:iz)
    !
    real(8) :: fac11,fac12,fac13
    real(8) :: fac21,fac22,fac23
    !
    real(8) :: r_v01_x_v01,r_v02_x_v02,r_v03_x_v03
    !
    real(8) :: ang(3)
    real(8) :: len
    !
    ! S.Iga060209=>
    real(8), optional:: critical
    ! S.Iga060209>=
    real(8):: epsi


    ! S.Iga060209=>
    if (.not.present(critical)) then
       !       critical = 0 !1d-10e
       epsi = 0.D0 !1d-10e
    else
       epsi=critical  !060224
    endif
    ! S.Iga060209>=


    !
    if(trim(polygon_type)=='ON_PLANE') then
       !
       !---- Note : On a plane,
       !----        area = | ourter product of two vectors |.
       v1(ix:iz)=b(ix:iz)-a(ix:iz)
       v2(ix:iz)=c(ix:iz)-a(ix:iz)
       !----
       w(ix) = v1(iy)*v2(iz)-v2(iy)*v1(iz)
       w(iy) = v1(iz)*v2(ix)-v2(iz)*v1(ix)
       w(iz) = v1(ix)*v2(iy)-v2(ix)*v1(iy)
       !
       area=0.5D0*sqrt(w(ix)*w(ix)+w(iy)*w(iy)+w(iz)*w(iz))
       !
       len = sqrt(a(ix)*a(ix)+a(iy)*a(iy)+a(iz)*a(iz))
       area=area*(radius/len)*(radius/len)
       !
       return
    elseif(trim(polygon_type)=='ON_SPHERE') then
       !
       !---- NOTE : On a unit sphere,
       !----        area = sum of angles - pi.
       v01(ix:iz)=a(ix:iz)
       v11(ix:iz)=b(ix:iz)-a(ix:iz)
       v21(ix:iz)=c(ix:iz)-a(ix:iz)
       !----
       v02(ix:iz)=b(ix:iz)
       v12(ix:iz)=a(ix:iz)-b(ix:iz)
       v22(ix:iz)=c(ix:iz)-b(ix:iz)
       !----
       v03(ix:iz)=c(ix:iz)
       v13(ix:iz)=a(ix:iz)-c(ix:iz)
       v23(ix:iz)=b(ix:iz)-c(ix:iz)
       !----
       r_v01_x_v01&
            =1.0D0/(v01(ix)*v01(ix)+v01(iy)*v01(iy)+v01(iz)*v01(iz))
       fac11=(v01(ix)*v11(ix)+v01(iy)*v11(iy)+v01(iz)*v11(iz))&
            *r_v01_x_v01
       fac21=(v01(ix)*v21(ix)+v01(iy)*v21(iy)+v01(iz)*v21(iz))&
            *r_v01_x_v01
       !---- Escape for the case arg=0 (S.Iga060209)
       area = 0.D0
       if ((v12(ix)**2+v12(iy)**2+v12(iz)**2) * r_v01_x_v01 <= epsi**2 ) return
       if ((v13(ix)**2+v13(iy)**2+v13(iz)**2) * r_v01_x_v01 <= epsi**2 ) return
       if ((v23(ix)**2+v23(iy)**2+v23(iz)**2) * r_v01_x_v01 <= epsi**2 ) return

       !----       !
       w11(ix)=v11(ix)-fac11*v01(ix)
       w11(iy)=v11(iy)-fac11*v01(iy)
       w11(iz)=v11(iz)-fac11*v01(iz)
       !
       w21(ix)=v21(ix)-fac21*v01(ix)
       w21(iy)=v21(iy)-fac21*v01(iy)
       w21(iz)=v21(iz)-fac21*v01(iz)
       !
       ang(1)=(w11(ix)*w21(ix)+w11(iy)*w21(iy)+w11(iz)*w21(iz))&
            /( sqrt(w11(ix)*w11(ix)+w11(iy)*w11(iy)+w11(iz)*w11(iz))&
            * sqrt(w21(ix)*w21(ix)+w21(iy)*w21(iy)+w21(iz)*w21(iz)) )
       if(ang(1)>1.0D0) ang(1) = 1.0D0
       if(ang(1)<-1.0D0) ang(1) = -1.0D0
       ang(1)=acos(ang(1))
       !
       r_v02_x_v02&
            =1.0D0/(v02(ix)*v02(ix)+v02(iy)*v02(iy)+v02(iz)*v02(iz))
       fac12=(v02(ix)*v12(ix)+v02(iy)*v12(iy)+v02(iz)*v12(iz))&
            *r_v02_x_v02
       fac22=(v02(ix)*v22(ix)+v02(iy)*v22(iy)+v02(iz)*v22(iz))&
            *r_v02_x_v02
       !
       w12(ix)=v12(ix)-fac12*v02(ix)
       w12(iy)=v12(iy)-fac12*v02(iy)
       w12(iz)=v12(iz)-fac12*v02(iz)
       !
       w22(ix)=v22(ix)-fac22*v02(ix)
       w22(iy)=v22(iy)-fac22*v02(iy)
       w22(iz)=v22(iz)-fac22*v02(iz)
       !
       ang(2)=(w12(ix)*w22(ix)+w12(iy)*w22(iy)+w12(iz)*w22(iz))&
            /( sqrt(w12(ix)*w12(ix)+w12(iy)*w12(iy)+w12(iz)*w12(iz))&
            *sqrt(w22(ix)*w22(ix)+w22(iy)*w22(iy)+w22(iz)*w22(iz)) )
       if(ang(2)>1.0D0) ang(2) = 1.0D0
       if(ang(2)<-1.0D0) ang(2) = -1.0D0
       ang(2)=acos(ang(2))
       !
       r_v03_x_v03&
            =1.0D0/(v03(ix)*v03(ix)+v03(iy)*v03(iy)+v03(iz)*v03(iz))
       fac13=(v03(ix)*v13(ix)+v03(iy)*v13(iy)+v03(iz)*v13(iz))&
            *r_v03_x_v03
       fac23=(v03(ix)*v23(ix)+v03(iy)*v23(iy)+v03(iz)*v23(iz))&
            *r_v03_x_v03
       !
       w13(ix)=v13(ix)-fac13*v03(ix)
       w13(iy)=v13(iy)-fac13*v03(iy)
       w13(iz)=v13(iz)-fac13*v03(iz)
       !
       w23(ix)=v23(ix)-fac23*v03(ix)
       w23(iy)=v23(iy)-fac23*v03(iy)
       w23(iz)=v23(iz)-fac23*v03(iz)
       !
       ang(3)=(w13(ix)*w23(ix)+w13(iy)*w23(iy)+w13(iz)*w23(iz))&
            /( sqrt(w13(ix)*w13(ix)+w13(iy)*w13(iy)+w13(iz)*w13(iz))&
            *sqrt(w23(ix)*w23(ix)+w23(iy)*w23(iy)+w23(iz)*w23(iz)) )
       if(ang(3)>1.0D0) ang(3) = 1.0D0
       if(ang(3)<-1.0D0) ang(3) = -1.0D0
       ang(3)=acos(ang(3))
       !----
       area=(ang(1)+ang(2)+ang(3)-pi)*radius*radius
       !
       return
       !
    endif
    !
  end function MISC_triangle_area
  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return
  !>
  function MISC_triangle_area_q( &
       a, b, c,                & !--- IN : three points vectors on a sphere.
       polygon_type,           & !--- IN : sphere triangle or plane one?
       radius,                 & !--- IN : radius
       critical)               & !--- IN : critical value to handle the case
                                 !         of ang=0 (optional) S.Iga060209
       result( area )            !--- OUT : triangle area
    !
    implicit none
    !
    integer, parameter :: ix = 1
    integer, parameter :: iy = 2
    integer, parameter :: iz = 3
    real(8), parameter :: pi  = 3.14159265358979323846E0_8
    !
    real(8) :: area
    real(8),intent(in) :: a(ix:iz),b(ix:iz),c(ix:iz)
    character(len=*), intent(in) :: polygon_type
    real(8) :: radius
    !
    !
    real(8) :: v01(ix:iz)
    real(8) :: v02(ix:iz)
    real(8) :: v03(ix:iz)
    !
    real(8) :: v11(ix:iz)
    real(8) :: v12(ix:iz)
    real(8) :: v13(ix:iz)
    !
    real(8) :: v21(ix:iz)
    real(8) :: v22(ix:iz)
    real(8) :: v23(ix:iz)
    real(8) :: w21(ix:iz)
    real(8) :: w22(ix:iz)
    real(8) :: w23(ix:iz)
    real(8) :: w11(ix:iz)
    real(8) :: w12(ix:iz)
    real(8) :: w13(ix:iz)

    real(8) :: v1(ix:iz)
    real(8) :: v2(ix:iz)
    real(8) :: w(ix:iz)
    !
    real(8) :: fac11,fac12,fac13
    real(8) :: fac21,fac22,fac23
    !
    real(8) :: r_v01_x_v01,r_v02_x_v02,r_v03_x_v03
    !
    real(8) :: ang(3)
    real(8) :: len
    real(8) :: a16(3)
    real(8) :: area16

    real(8), optional:: critical
    real(8):: epsi


    if ( .not. present(critical)) then
       epsi = 0.E0_8
    else
       epsi = real(critical,kind=8)
    endif

    if(trim(polygon_type)=='ON_PLANE') then
       !---- Note : On a plane,
       !----        area = | ourter product of two vectors |.
       v1(ix:iz)= real(b(ix:iz),kind=8) - real(a(ix:iz),kind=8)
       v2(ix:iz)= real(c(ix:iz),kind=8) - real(a(ix:iz),kind=8)

       w(ix) = v1(iy)*v2(iz)-v2(iy)*v1(iz)
       w(iy) = v1(iz)*v2(ix)-v2(iz)*v1(ix)
       w(iz) = v1(ix)*v2(iy)-v2(ix)*v1(iy)

       area16 = 0.5E0_8 * sqrt(w(ix)*w(ix)+w(iy)*w(iy)+w(iz)*w(iz))

       a16(ix:iz) = real(a(ix:iz),kind=8)
       len = sqrt( a16(ix)*a16(ix)+a16(iy)*a16(iy)+a16(iz)*a16(iz) )
       area16 = area16 * (radius/len)*(radius/len)

       area = real(area16,kind=8)

       return
    elseif(trim(polygon_type)=='ON_SPHERE') then
       !
       !---- NOTE : On a unit sphere,
       !----        area = sum of angles - pi.
       v01(ix:iz)=real(a(ix:iz),kind=8)
       v11(ix:iz)=real(b(ix:iz),kind=8)-real(a(ix:iz),kind=8)
       v21(ix:iz)=real(c(ix:iz),kind=8)-real(a(ix:iz),kind=8)
       !----
       v02(ix:iz)=real(b(ix:iz),kind=8)
       v12(ix:iz)=real(a(ix:iz),kind=8)-real(b(ix:iz),kind=8)
       v22(ix:iz)=real(c(ix:iz),kind=8)-real(b(ix:iz),kind=8)
       !----
       v03(ix:iz)=real(c(ix:iz),kind=8)
       v13(ix:iz)=real(a(ix:iz),kind=8)-real(c(ix:iz),kind=8)
       v23(ix:iz)=real(b(ix:iz),kind=8)-real(c(ix:iz),kind=8)
       !----
       r_v01_x_v01&
            =1.0D0/(v01(ix)*v01(ix)+v01(iy)*v01(iy)+v01(iz)*v01(iz))
       fac11=(v01(ix)*v11(ix)+v01(iy)*v11(iy)+v01(iz)*v11(iz))&
            *r_v01_x_v01
       fac21=(v01(ix)*v21(ix)+v01(iy)*v21(iy)+v01(iz)*v21(iz))&
            *r_v01_x_v01
       !---- Escape for the case arg=0 (S.Iga060209)
       area=0E0_8

       if ((v12(ix)**2+v12(iy)**2+v12(iz)**2) * r_v01_x_v01 <= epsi**2 ) return
       if ((v13(ix)**2+v13(iy)**2+v13(iz)**2) * r_v01_x_v01 <= epsi**2 ) return
       if ((v23(ix)**2+v23(iy)**2+v23(iz)**2) * r_v01_x_v01 <= epsi**2 ) return

       !----       !
       w11(ix)=v11(ix)-fac11*v01(ix)
       w11(iy)=v11(iy)-fac11*v01(iy)
       w11(iz)=v11(iz)-fac11*v01(iz)
       !
       w21(ix)=v21(ix)-fac21*v01(ix)
       w21(iy)=v21(iy)-fac21*v01(iy)
       w21(iz)=v21(iz)-fac21*v01(iz)
       !
       ang(1)=(w11(ix)*w21(ix)+w11(iy)*w21(iy)+w11(iz)*w21(iz))&
            /( sqrt(w11(ix)*w11(ix)+w11(iy)*w11(iy)+w11(iz)*w11(iz))&
            * sqrt(w21(ix)*w21(ix)+w21(iy)*w21(iy)+w21(iz)*w21(iz)) )
       if( ang(1) >  1.0E0_8 ) ang(1) =  1.E0_8
       if( ang(1) <- 1.0E0_8 ) ang(1) = -1.E0_8
       ang(1)=acos(ang(1))
       !
       r_v02_x_v02&
            =1.0D0/(v02(ix)*v02(ix)+v02(iy)*v02(iy)+v02(iz)*v02(iz))
       fac12=(v02(ix)*v12(ix)+v02(iy)*v12(iy)+v02(iz)*v12(iz))&
            *r_v02_x_v02
       fac22=(v02(ix)*v22(ix)+v02(iy)*v22(iy)+v02(iz)*v22(iz))&
            *r_v02_x_v02
       !
       w12(ix)=v12(ix)-fac12*v02(ix)
       w12(iy)=v12(iy)-fac12*v02(iy)
       w12(iz)=v12(iz)-fac12*v02(iz)
       !
       w22(ix)=v22(ix)-fac22*v02(ix)
       w22(iy)=v22(iy)-fac22*v02(iy)
       w22(iz)=v22(iz)-fac22*v02(iz)
       !
       ang(2)=(w12(ix)*w22(ix)+w12(iy)*w22(iy)+w12(iz)*w22(iz))&
            /( sqrt(w12(ix)*w12(ix)+w12(iy)*w12(iy)+w12(iz)*w12(iz))&
            *sqrt(w22(ix)*w22(ix)+w22(iy)*w22(iy)+w22(iz)*w22(iz)) )

       if( ang(2) >  1.0E0_8 ) ang(2) =  1.E0_8
       if( ang(2) <- 1.0E0_8 ) ang(2) = -1.E0_8

       ang(2)=acos(ang(2))

       r_v03_x_v03&
            =1.0D0/(v03(ix)*v03(ix)+v03(iy)*v03(iy)+v03(iz)*v03(iz))
       fac13=(v03(ix)*v13(ix)+v03(iy)*v13(iy)+v03(iz)*v13(iz))&
            *r_v03_x_v03
       fac23=(v03(ix)*v23(ix)+v03(iy)*v23(iy)+v03(iz)*v23(iz))&
            *r_v03_x_v03
       !
       w13(ix)=v13(ix)-fac13*v03(ix)
       w13(iy)=v13(iy)-fac13*v03(iy)
       w13(iz)=v13(iz)-fac13*v03(iz)
       !
       w23(ix)=v23(ix)-fac23*v03(ix)
       w23(iy)=v23(iy)-fac23*v03(iy)
       w23(iz)=v23(iz)-fac23*v03(iz)
       !
       ang(3)=(w13(ix)*w23(ix)+w13(iy)*w23(iy)+w13(iz)*w23(iz))&
            /( sqrt(w13(ix)*w13(ix)+w13(iy)*w13(iy)+w13(iz)*w13(iz))&
            *sqrt(w23(ix)*w23(ix)+w23(iy)*w23(iy)+w23(iz)*w23(iz)) )

       if( ang(3) >  1.0E0_8 ) ang(3) =  1.E0_8
       if( ang(3) <- 1.0E0_8 ) ang(3) = -1.E0_8

       ang(3)=acos(ang(3))
       !----
       area16=(ang(1)+ang(2)+ang(3)-pi)*radius*radius

       area = real(area16,kind=8)

       return
    endif

  end function MISC_triangle_area_q

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_mk_gmtrvec
  !>
  subroutine MISC_mk_gmtrvec( &
       vs, ve,                & !--- IN : vectors at the start and the end
       tv,                    & !--- INOUT : tangential vector
       nv,                    & !--- INOUT : normal vector
       polygon_type,          & !--- IN : polygon type
       radius )                 !--- IN : radius
    !
    implicit none
    !
    integer, parameter :: ix = 1
    integer, parameter :: iy = 2
    integer, parameter :: iz = 3
    !
    real(8),intent(in) :: vs(ix:iz)
    real(8),intent(in) :: ve(ix:iz)
    real(8),intent(inout) :: tv(ix:iz)
    real(8),intent(inout) :: nv(ix:iz)
    character(len=*), intent(in) :: polygon_type
    real(8), intent(in) :: radius
    real(8) :: vec_len
    !
    real(8) :: len
    real(8) :: fact_nv,fact_tv
    !
    if(trim(polygon_type)=='ON_SPHERE') then
       !
       !--- NOTE : Length of a geodesic line is caluclatd
       !---        by (angle X radius).
       vec_len = sqrt(vs(ix)*vs(ix)+vs(iy)*vs(iy)+vs(iz)*vs(iz))
       if(vec_len/=0.0D0) then
          len=acos((vs(ix)*ve(ix)+vs(iy)*ve(iy)+vs(iz)*ve(iz))/vec_len/vec_len)&
               *radius
       else
          len = 0.0D0
       endif
       !
    elseif(trim(polygon_type)=='ON_PLANE') then
       !
       !--- NOTE : Length of a line
       len=sqrt(&
            +(vs(ix)-ve(ix))*(vs(ix)-ve(ix))&
            +(vs(iy)-ve(iy))*(vs(iy)-ve(iy))&
            +(vs(iz)-ve(iz))*(vs(iz)-ve(iz))&
            )
    endif
    !
    !
    !--- calculate normal and tangential vecors
    nv(ix)=vs(iy)*ve(iz)-vs(iz)*ve(iy)
    nv(iy)=vs(iz)*ve(ix)-vs(ix)*ve(iz)
    nv(iz)=vs(ix)*ve(iy)-vs(iy)*ve(ix)
    tv(ix)=ve(ix)-vs(ix)
    tv(iy)=ve(iy)-vs(iy)
    tv(iz)=ve(iz)-vs(iz)
    !
    !--- scaling to radius ( normal vector )
    fact_nv=len/sqrt(nv(ix)*nv(ix)+nv(iy)*nv(iy)+nv(iz)*nv(iz))
    nv(ix)=nv(ix)*fact_nv
    nv(iy)=nv(iy)*fact_nv
    nv(iz)=nv(iz)*fact_nv
    !
    !--- scaling to radius ( tangential vector )
    fact_tv=len/sqrt(tv(ix)*tv(ix)+tv(iy)*tv(iy)+tv(iz)*tv(iz))
    tv(ix)=tv(ix)*fact_tv
    tv(iy)=tv(iy)*fact_tv
    tv(iz)=tv(iz)*fact_tv
    !
    return
    !
  end subroutine MISC_mk_gmtrvec
  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_sec2ymdhms
  !>
  subroutine MISC_sec2ymdhms(           &
       ctime_second,                    & !--- IN : target time [sec]
       year, month, day, hour, min, sec & !--- OUT : year, month day hour minite, second
       )
    !
    implicit none
    !
    real(8), intent(in)  :: ctime_second
    integer, intent(out) :: sec
    integer, intent(out) :: min
    integer, intent(out) :: hour
    integer, intent(out) :: day
    integer, intent(out) :: month
    integer, intent(out) :: year
    !
    !    integer(4) :: second
    integer :: second      ! S.Iga 060212
    !
    second = idnint(ctime_second)
    year = second/(60*60*24*30*12)
    !
    second = mod(second,60*60*24*30*12)
    month = second/(60*60*24*30)
    !
    second = mod(second,60*60*24*30)
    day = second/(60*60*24)
    !
    second = mod(second,60*60*24)
    hour = second/(60*60)
    !
    second = mod(second,60*60)
    min = second/(60)
    !
    second = mod(second,60)
    sec = second

  end subroutine MISC_sec2ymdhms
  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_ymdhms2sec
  !>
  subroutine MISC_ymdhms2sec(            &
       ctime_second,                     & !--- OUT : target time [sec]
       year, month, day, hour, min, sec  & !--- IN : year, month day hour minite, second
       )

    implicit none

    real(8), intent(out) :: ctime_second
    integer, intent(in)  :: sec
    integer, intent(in)  :: min
    integer, intent(in)  :: hour
    integer, intent(in)  :: day
    integer, intent(in)  :: month
    integer, intent(in)  :: year

    ctime_second = year  * 60.D0 * 60.D0 * 24.D0 * 30.D0 * 12.D0 &
                 + month * 60.D0 * 60.D0 * 24.D0 * 30.D0         &
                 + day   * 60.D0 * 60.D0 * 24.D0                 &
                 + hour  * 60.D0 * 60.D0                         &
                 + min   * 60.D0                                 &
                 + sec

    return
  end subroutine MISC_ymdhms2sec

  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine MISC_msg_nmerror
  !>
  subroutine MISC_msg_nmerror( &
       ierr,                 & !--- IN : error id
       fid,                  & !--- IN : file id
       namelist_name,        & !--- IN : namelist name
       sub_name,             & !--- IN : subroutine name
       mod_name              & !--- IN : module name
       )
    implicit none
    integer, intent(in) ::  ierr
    integer, intent(in) ::  fid
    character(len=*) :: namelist_name
    character(len=*) :: sub_name
    character(len=*) :: mod_name
    if(ierr<0) then
       write(fid,*) &
            'Msg : Sub[',trim(sub_name),']/Mod[',trim(mod_name),']'
       write(fid,*) &
            ' *** Not found namelist. ',trim(namelist_name)
       write(fid,*) &
            ' *** Use default values.'
    elseif(ierr>0) then !--- fatal error
       write(*,*) &
            'Msg : Sub[',trim(sub_name),']/Mod[',trim(mod_name),']'
       write(*,*) &
            ' *** WARNING : Not appropriate names in namelist!! ',&
            trim(namelist_name),' CHECK!!'
    endif
  end subroutine MISC_msg_nmerror
  !-----------------------------------------------------------------------------
  !>
  !> Description of the function %NAME
  !> @return
  !>
  function MISC_gammafunc( xx ) result(f)
    implicit none
    real(8), intent(in) :: xx
    real(8) :: f
    real(8) :: coef(6)=(/&
         +76.18009172947146D0,&
         -86.50532032941677D0,&
         +24.01409824083091D0,&
         -1.231739572450155D0,&
         +0.1208650973866179D-2,&
         -0.5395239384953D-5&
         /)
    integer :: j
    real(8) :: x,y,tmp,ser

    x=xx
    y=x
    tmp=x+5.5D0
    tmp = tmp - (x+0.5)*log(tmp)
    ser=1.000000000190015D0
    do j=1,6
       y=y+1
       ser = ser+coef(j)/y
    enddo
    f = exp(-tmp+log(2.5066282746310005D0*ser/x))
  end function MISC_gammafunc
  !-----------------------------------------------------------------------------
  ![Add] T.Mitsui 08/04/12
  !>
  !> Section 6.2 in "Numerical Recipes in C"
  !> function is represented by(6.2.1)
  !> g(x,alpha)=1/gamma(alpha)*integral_0^x { exp(-t) * t^(alpha-1) }dt
  !> g(0)=0 and g(+infinity)=1
  !> @return g
  !>
  function MISC_igammafunc( x, alpha, nm ) result(g)
    implicit none

    real(8), intent(in) :: x
    real(8), intent(in) :: alpha

    real(8), intent(in), optional :: nm ! nmax
    real(8) :: g

    real(8), parameter :: eps=1.D-30
    real(8) :: a,b,c,d,h,an,del
    real(8) :: gm, lgm

    integer :: nmax, n
    !---------------------------------------------------------------------------

    gm  = MISC_gammafunc(alpha)
    lgm = log(gm)

    if ( x < alpha + 1.D0 )then ! series expansion (6.2.5)

       if ( present(nm) ) then
          nmax = int(nm)
       else
          ! nmax=10 is enough to represent Droplet Size Distribution
          ! relative error is below 1% arround local maximum(x=alpha-1)
          nmax = 10
       endif

       a = 1.D0 / alpha
       g = a

       do n = 1, nmax
          a = a*x/(alpha+n)
          g = g+a
       enddo
       g = g*exp( -x + alpha*log(x) - lgm )

    else ! continued fraction expansion (6.2.6)

       if ( present(nm) ) then
          nmax = int(nm)
       else
          ! nmax=2 is enough to represent Droplet Size Distribution
          ! relative error is below 1% arround local maximum(x=alpha-1)
          nmax = 2
       endif

       b = x+1.d0-alpha
       c = 1.d0/eps
       d = 1.d0/b
       h=d
       do n=1, nmax
          an = -n*(n-alpha)
          b  = b + 2.d0
          d  = an*d+b
          if(abs(d)<eps) d=eps
          c  = b+an/c
          if(abs(c)<eps) c=eps
          d  = 1.d0/d
          del= d*c
          h  = h*del
       enddo
       g = 1.d0 - exp( -x+alpha*log(x) - lgm )*h
    endif

    return
  end function MISC_igammafunc

  !-----------------------------------------------------------------------------
  ![Add] T.Mitsui 08/04/12
  !>
  !> Description of the function %NAME
  !> @return
  !>
  function MISC_iigammafunc( p, alpha ) result(x)
    ! p = MISC_igammafunc( x, alpha, nm )
    ! Didonato and Morris, 1986,
    ! Computation of the Incomplete Gamma Function Ratios and their Inverse,
    ! ACM Transactions on Mathmatical Software, Vol.12, No.4,
    !
    implicit none
    !
    real(8), intent(in) :: alpha  !
    real(8), intent(in) :: p      !
    real(8) :: x
    !
    real(8) :: s,s2,s3,s4,s5
    real(8) :: t,t2,t3,t4
    real(8) :: tau
    real(8) :: m
    real(8) :: a_n, a_d
    real(8) :: sa ! sqrt(alpha)
    !
    real(8), parameter :: a0=3.31125922108741d0
    real(8), parameter :: a1=11.6616720288968d0
    real(8), parameter :: a2=4.28342155967104d0
    real(8), parameter :: a3=0.213623493715853d0
    real(8), parameter :: b1=6.61053765625462d0
    real(8), parameter :: b2=6.40691597760039d0
    real(8), parameter :: b3=1.27364489782223d0
    real(8), parameter :: b4=0.361170810188420d-1
    !
    real(8), parameter :: eps = 1.d-10
    !
    if(alpha < 1.d0)then
       write(*,*) "Invalid parameter, alpha=",alpha
       write(*,*) "This routine does not support."
       return
    endif
    if( (alpha > 1.d0-eps) .and. (alpha < 1.d0+eps) )then ! alpha = 1
       x = -log(1.d0 - p)
       return
    endif
    !
    if(p>=0.5d0)then
       m  = 1.d0
       tau= 1.d0-p
    else
       m  =-1.d0
       tau= p
    endif
    t    = sqrt(-2.d0*log(tau))
    t2   = t*t
    t3   = t2*t
    t4   = t2*t2
    a_n  = a0   + a1*t + a2*t2 + a3*t3
    a_d  = 1.d0 + b1*t + b2*t2 + b3*t3 + b4*t4
    s    = m*(t-a_n/a_d) !(32)
    !
    s2   = s*s
    s3   = s*s2
    s4   = s2*s2
    s5   = s3*s2
    sa   = sqrt(alpha)
    !
    x    = alpha + s*sa + (s2-1.d0)/3.d0 + (s3 - 7.d0*s)/(36.d0*sa) &
         - (3.d0*s4 + 7.d0*s2 - 16.d0)/(810.d0*alpha)               &
         + (9.d0*s5 + 256.d0*s3 - 433.d0*s)/(38880.d0*alpha*sa)
    !
    return
  end function MISC_iigammafunc

  !-----------------------------------------------------------------------------
  subroutine MISC_3dvec_cross( nv, a, b, c, d )
    ! exterior product of vector a->b and c->d
    implicit none

    real(8), intent(out) :: nv(3)                  ! normal vector
    real(8), intent(in ) :: a(3), b(3), c(3), d(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------

    nv(1) = ( b(2)-a(2) ) * ( d(3)-c(3) ) &
          - ( b(3)-a(3) ) * ( d(2)-c(2) )
    nv(2) = ( b(3)-a(3) ) * ( d(1)-c(1) ) &
          - ( b(1)-a(1) ) * ( d(3)-c(3) )
    nv(3) = ( b(1)-a(1) ) * ( d(2)-c(2) ) &
          - ( b(2)-a(2) ) * ( d(1)-c(1) )

    return
  end subroutine MISC_3dvec_cross

  !-----------------------------------------------------------------------------
  subroutine MISC_3dvec_dot( l, a, b, c, d )
    ! interior product of vector a->b and c->d
    implicit none

    real(8), intent(out) :: l
    real(8), intent(in ) :: a(3), b(3), c(3), d(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------
    ! if a=c=zero-vector and b=d, result is abs|a|^2

    l = ( b(1)-a(1) ) * ( d(1)-c(1) ) &
      + ( b(2)-a(2) ) * ( d(2)-c(2) ) &
      + ( b(3)-a(3) ) * ( d(3)-c(3) ) 

    return
  end subroutine MISC_3dvec_dot

  !-----------------------------------------------------------------------------
  subroutine MISC_3dvec_abs( l, a )
    ! length of vector o->a
    implicit none

    real(8), intent(out) :: l
    real(8), intent(in ) :: a(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------

    l = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    l = sqrt(l)

    return
  end subroutine MISC_3dvec_abs

  !---------------------------------------------------------------------
  subroutine MISC_3dvec_angle( angle, a, b, c )
    ! calc angle between two vector(b->a,b->c)
    implicit none

    real(8), intent(out) :: angle
    real(8), intent(in ) :: a(3), b(3), c(3)

    real(8) :: nv(3), nvlenS, nvlenC
    !---------------------------------------------------------------------

    call MISC_3dvec_dot  ( nvlenC, b, a, b, c )
    call MISC_3dvec_cross( nv(:),  b, a, b, c )
    call MISC_3dvec_abs  ( nvlenS, nv(:) )
    angle = atan2( nvlenS, nvlenC )

    return
  end subroutine MISC_3dvec_angle

  !-----------------------------------------------------------------------------
  function MISC_3Dvec_triangle( &
       a, b, c,      & !--- IN : three point vectors on a sphere.
       polygon_type, & !--- IN : sphere triangle or plane one?
       radius      ) & !--- IN : radius
       result(area)    !--- OUT : triangle area
    implicit none

    real(8) :: area
    real(8),          intent( in) :: a(3), b(3), c(3)
    character(len=*), intent( in) :: polygon_type
    real(8),          intent( in) :: radius

    real(8), parameter :: o(3) = 0.D0

    ! ON_PLANE
    real(8) :: abc(3)
    real(8) :: prd, r

    ! ON_SPHERE
    real(8) :: angle(3)
    real(8) :: oaob(3), oaoc(3)
    real(8) :: oboc(3), oboa(3)
    real(8) :: ocoa(3), ocob(3)
    real(8) :: abab, acac
    real(8) :: bcbc, baba
    real(8) :: caca, cbcb

    real(8), parameter :: eps = 1.D-16
    real(8) :: pi
    !---------------------------------------------------------------------------

    pi = atan(1.D0) * 4.D0

    area = 0.D0

    if(trim(polygon_type)=='ON_PLANE') then
       !
       !---- Note : On a plane,
       !----        area = | ourter product of two vectors |.
       !
       call MISC_3dvec_cross( abc(:), a(:), b(:), a(:), c(:) )
       call MISC_3dvec_abs( prd, abc(:) )
       call MISC_3dvec_abs( r  , a(:)   )

       prd = 0.5D0 * prd !! triangle area
       if ( r < eps * radius ) then
          print *, "zero length?", a(:)
       else
          r = 1.D0 / r   !! 1 / length
       endif

       area = prd * r*r * radius*radius

    elseif(trim(polygon_type)=='ON_SPHERE') then
       !
       !---- NOTE : On a unit sphere,
       !----        area = sum of angles - pi.
       !

       ! angle 1
       call MISC_3dvec_cross( oaob(:), o(:), a(:), o(:), b(:) )
       call MISC_3dvec_cross( oaoc(:), o(:), a(:), o(:), c(:) )
       call MISC_3dvec_abs( abab, oaob(:) )
       call MISC_3dvec_abs( acac, oaoc(:) )

       if ( abab < eps * radius .OR. acac < eps * radius ) then
          write(*,'(A,3(E20.10))') "zero length abab or acac:", abab/radius, acac/radius
          return
       endif

       call MISC_3dvec_angle( angle(1), oaob(:), o(:), oaoc(:) )

       ! angle 2
       call MISC_3dvec_cross( oboc(:), o(:), b(:), o(:), c(:) )
       oboa(:) = -oaob(:)
       call MISC_3dvec_abs( bcbc, oboc(:) )
       baba = abab

       if ( bcbc < eps * radius .OR. baba < eps * radius ) then
          write(*,'(A,3(E20.10))') "zero length bcbc or baba:", bcbc/radius, baba/radius
          return
       endif

       call MISC_3dvec_angle( angle(2), oboc(:), o(:), oboa(:) )

       ! angle 3
       ocoa(:) = -oaoc(:)
       ocob(:) = -oboc(:)
       caca = acac
       cbcb = bcbc

       if ( caca < eps * radius .OR. cbcb < eps * radius ) then
          write(*,'(A,3(E20.10))') "zero length caca or cbcb:", caca/radius, cbcb/radius
          return
       endif

       call MISC_3dvec_angle( angle(3), ocoa(:), o(:), ocob(:) )

       ! calc area
       area = ( angle(1)+angle(2)+angle(3)-pi ) * radius*radius

    endif

    return
  end function MISC_3Dvec_triangle

  !-----------------------------------------------------------------------------
  subroutine MISC_3dvec_intersec( ifcross, p, a, b, c, d )
    ! judge intersection of two vector
    implicit none

    logical, intent(out) :: ifcross
    ! .true. : line a->b and c->d intersect
    ! .false.: line a->b and c->d do not intersect and p = (0,0)
    real(8), intent(out) :: p(3) ! intersection point
    real(8), intent(in ) :: a(3), b(3), c(3), d(3)

    real(8), parameter :: o(3) = 0.D0

    real(8)            :: oaob(3), ocod(3), cdab(3)
    real(8)            :: ip, length
    real(8)            :: angle_aop, angle_pob, angle_aob
    real(8)            :: angle_cop, angle_pod, angle_cod

    real(8), parameter :: eps = 1.D-12
    !---------------------------------------------------------------------

    call MISC_3dvec_cross( oaob, o, a, o, b )
    call MISC_3dvec_cross( ocod, o, c, o, d )
    call MISC_3dvec_cross( cdab, o, ocod, o, oaob )

    call MISC_3dvec_abs  ( length, cdab )
    call MISC_3dvec_dot  ( ip, o, cdab, o, a )

    p(:) = cdab(:) / sign(length,ip)
!    write(ADM_LOG_FID,*), "p:", p(:)

    call MISC_3dvec_angle( angle_aop, a, o, p )
    call MISC_3dvec_angle( angle_pob, p, o, b )
    call MISC_3dvec_angle( angle_aob, a, o, b )
!    write(ADM_LOG_FID,*), "angle a-p-b:", angle_aop, angle_pob, angle_aob

    call MISC_3dvec_angle( angle_cop, c, o, p )
    call MISC_3dvec_angle( angle_pod, p, o, d )
    call MISC_3dvec_angle( angle_cod, c, o, d )
!    write(ADM_LOG_FID,*), "angle c-p-d:", angle_cop, angle_pod, angle_cod

!    write(ADM_LOG_FID,*), "judge:", angle_aob-(angle_aop+angle_pob), angle_cod-(angle_cop+angle_pod)

    ! --- judge intersection
    if (       abs(angle_aob-(angle_aop+angle_pob)) < eps &
         .AND. abs(angle_cod-(angle_cop+angle_pod)) < eps &
         .AND. abs(angle_aop) > eps                       &
         .AND. abs(angle_pob) > eps                       &
         .AND. abs(angle_cop) > eps                       &
         .AND. abs(angle_pod) > eps                       ) then
       ifcross = .true.
    else
       ifcross = .false.
       p(:) = 0.D0
    endif

    return
  end subroutine MISC_3dvec_intersec

  !---------------------------------------------------------------------
  subroutine MISC_3dvec_anticlockwise( vertex, nvert )
    ! bubble sort anticlockwise by angle
    implicit none

    integer, intent(in)    :: nvert
    real(8), intent(inout) :: vertex(nvert,3)

    real(8), parameter :: o(3) = 0.D0
    real(8)            :: v1(3), v2(3), v3(3)
    real(8)            :: xp(3), ip
    real(8)            :: angle1, angle2

    real(8), parameter :: eps = 1.D-12

    integer :: i, j
    !---------------------------------------------------------------------

    do j = 2  , nvert-1
    do i = j+1, nvert
       v1(:) = vertex(1,:)
       v2(:) = vertex(j,:)
       v3(:) = vertex(i,:)

       call MISC_3dvec_cross( xp(:), v1(:), v2(:), v1(:), v3(:) )
       call MISC_3dvec_dot  ( ip, o(:), v1(:), o(:), xp(:) )

       if ( ip < -eps ) then ! right hand : exchange
!          write(ADM_LOG_FID,*) 'exchange by ip', i, '<->',j
          vertex(i,:) = v2(:)
          vertex(j,:) = v3(:)
       endif

    enddo
    enddo

    v1(:) = vertex(1,:)
    v2(:) = vertex(2,:)
    v3(:) = vertex(3,:)
    ! if 1->2->3 is on the line
    call MISC_3dvec_cross( xp(:), v1(:), v2(:), v1(:), v3(:) )
    call MISC_3dvec_dot  ( ip, o(:), v1(:), o(:), xp(:) )
    call MISC_3dvec_angle( angle1, v1(:), o, v2(:) )
    call MISC_3dvec_angle( angle2, v1(:), o, v3(:) )
!    write(ADM_LOG_FID,*) ip, angle1, angle2, abs(angle1)-abs(angle2)

    if (       abs(ip)                 < eps  &      ! on the same line
         .AND. abs(angle2)-abs(angle1) < 0.D0 ) then ! which is far?
!       write(ADM_LOG_FID,*) 'exchange by angle', 2, '<->', 3
       vertex(2,:) = v3(:)
       vertex(3,:) = v2(:)
    endif

    v2(:) = vertex(nvert  ,:)
    v3(:) = vertex(nvert-1,:)
    ! if 1->nvert->nvert-1 is on the line
    call MISC_3dvec_cross( xp(:), v1(:), v2(:), v1(:), v3(:) )
    call MISC_3dvec_dot  ( ip, o(:), v1(:), o(:), xp(:) )
    call MISC_3dvec_angle( angle1, v1(:), o, v2(:) )
    call MISC_3dvec_angle( angle2, v1(:), o, v3(:) )
!    write(ADM_LOG_FID,*) ip, angle1, angle2, abs(angle1)-abs(angle2)

    if (       abs(ip)                 < eps  &      ! on the same line
         .AND. abs(angle2)-abs(angle1) < 0.D0 ) then ! which is far?
!       write(ADM_LOG_FID,*) 'exchange by angle', nvert, '<->', nvert-1
       vertex(nvert,  :) = v3(:)
       vertex(nvert-1,:) = v2(:)
    endif

    return
  end subroutine MISC_3dvec_anticlockwise

  ![Add] H.Yashiro 11/06/01
  !-----------------------------------------------------------------------------
  subroutine MISC_get_cartesian( &
       x, y, z,  & !--- INOUT : Cartesian coordinate ( on the sphere )
       lat, lon, & !--- IN    : latitude and longitude, radian
       radius    ) !--- IN    : radius
    implicit none

    real(8),intent(inout) :: x, y, z
    real(8),intent(in)    :: lat, lon
    real(8),intent(in)    :: radius
    !---------------------------------------------------------------------------

    x = radius * cos(lat) * cos(lon)
    y = radius * cos(lat) * sin(lon) 
    z = radius * sin(lat)

    return
  end subroutine MISC_get_cartesian

  !-----------------------------------------------------------------------
  ! Get horizontal distance on the sphere
  !  2008/09/10 [Add] M.Hara
  !  The formuation is Vincentry (1975)
  !  http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
  !  2012/11/07 [Mod] H.Yashiro
  subroutine MISC_get_distance( &
       r,    &
       lon1, &
       lat1, &
       lon2, &
       lat2, &
       dist  )
    implicit none

    real(8), intent(in)  :: r          ! radius in meter
    real(8), intent(in)  :: lon1, lat1 ! in radian
    real(8), intent(in)  :: lon2, lat2 ! in radian
    real(8), intent(out) :: dist       ! distance of the two points in meter

    real(8) :: gmm, gno_x, gno_y
    !-----------------------------------------------------------------------

    gmm = sin(lat1) * sin(lat2) &
        + cos(lat1) * cos(lat2) * cos(lon2-lon1)

    gno_x = gmm * ( cos(lat2) * sin(lon2-lon1) )
    gno_y = gmm * ( cos(lat1) * sin(lat2) &
                  - sin(lat1) * cos(lat2) * cos(lon2-lon1) )

    dist = r * atan2( sqrt(gno_x*gno_x+gno_y*gno_y), gmm )

    return
  end subroutine MISC_get_distance

  !-----------------------------------------------------------------------------
end module mod_misc
!-------------------------------------------------------------------------------


