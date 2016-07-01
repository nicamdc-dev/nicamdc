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
!! @li      2009-07-17 (Y.Yamada)  Add func[MISC_triangle_area_q].
!! @li      2011-11-11 (H.Yashiro) [Add] vector calculation suite
!!
!<
module mod_misc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters, variables & subroutines
  !
  public :: MISC_make_idstr        !--- make file name with a number
  public :: MISC_get_available_fid !--- get an available file ID
  public :: MISC_get_fid           !--- get an information of open or no
  public :: MISC_get_latlon        !--- calculate (lat,lon) from (x,y,z)
  public :: MISC_msg_nmerror       !--- output error message
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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> make extention with process number
  subroutine MISC_make_idstr( &
       str,    &
       prefix, &
       ext,    &
       numID,  &
       digit   )
    implicit none

    character(len=*),  intent(out) :: str    !< combined extention string
    character(len=*),  intent(in)  :: prefix !< prefix
    character(len=*),  intent(in)  :: ext    !< extention ( e.g. .rgn )
    integer,           intent(in)  :: numID  !< number
    integer, optional, intent(in)  :: digit  !< digit

    logical, parameter            :: NSTR_ZERO_START = .true. ! number of separated file starts from 0 ?
    integer, parameter            :: NSTR_MAX_DIGIT  = 5      ! digit of separated file

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
  !> Search and get available machine id
  !> @return fid
  function MISC_get_available_fid() result(fid)
    implicit none

    integer :: fid

    integer, parameter :: min_fid =  7 !< minimum available fid
    integer, parameter :: max_fid = 99 !< maximum available fid

    logical :: i_opened
    !---------------------------------------------------------------------------

    do fid = min_fid, max_fid
       inquire(fid,opened=i_opened)
       if( .NOT. i_opened ) return
    enddo

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
    real(RP),intent(inout) :: lat, lon
    real(RP),intent(in) :: x,y,z
    !
    real(RP), parameter :: epsilon = 1.E-30_RP
    real(RP) :: leng,leng_xy
    !
    leng=sqrt(x*x+y*y+z*z)
    !
    ! --- vector length equals to zero.
    if(leng<epsilon) then
       lat=0.0_RP
       lon=0.0_RP
       return
    endif
    ! --- vector is parallele to z axis.
    if(z/leng>=1.0_RP) then
       lat=asin(1.0_RP)
       lon=0.0_RP
       return
    elseif(z/leng<=-1.0_RP) then
       lat=asin(-1.0_RP)
       lon=0.0_RP
       return
    endif
    ! --- not parallele to z axis
    lat=asin(z/leng)
    !
    leng_xy=sqrt(x*x+y*y)
    if(leng_xy<epsilon) then
       lon=0.0_RP
       return
    endif
    if(x/leng_xy>=1.0_RP) then
       lon=acos(1.0_RP)
       if(y<0.0_RP) lon=-lon
       return
    elseif(x/leng_xy<=-1.0_RP) then
       lon=acos(-1.0_RP)
       if(y<0.0_RP) lon=-lon
       return
    endif
    lon=acos(x/leng_xy)
    if(y<0.0_RP) lon=-lon
    return
  end subroutine MISC_get_latlon

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
  subroutine MISC_3dvec_cross( nv, a, b, c, d )
    ! exterior product of vector a->b and c->d
    implicit none

    real(RP), intent(out) :: nv(3)                  ! normal vector
    real(RP), intent(in ) :: a(3), b(3), c(3), d(3) ! x,y,z(cartesian)
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

    real(RP), intent(out) :: l
    real(RP), intent(in ) :: a(3), b(3), c(3), d(3) ! x,y,z(cartesian)
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

    real(RP), intent(out) :: l
    real(RP), intent(in ) :: a(3) ! x,y,z(cartesian)
    !---------------------------------------------------------------------------

    l = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    l = sqrt(l)

    return
  end subroutine MISC_3dvec_abs

  !---------------------------------------------------------------------
  subroutine MISC_3dvec_angle( angle, a, b, c )
    ! calc angle between two vector(b->a,b->c)
    implicit none

    real(RP), intent(out) :: angle
    real(RP), intent(in ) :: a(3), b(3), c(3)

    real(RP) :: nv(3), nvlenS, nvlenC
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

    real(RP) :: area
    real(RP),          intent(in)  :: a(3), b(3), c(3)
    character(len=*), intent(in)  :: polygon_type
    real(RP),          intent(in)  :: radius

    real(RP), parameter :: o(3) = 0.0_RP

    ! ON_PLANE
    real(RP) :: abc(3)
    real(RP) :: prd, r

    ! ON_SPHERE
    real(RP) :: angle(3)
    real(RP) :: oaob(3), oaoc(3)
    real(RP) :: oboc(3), oboa(3)
    real(RP) :: ocoa(3), ocob(3)
    real(RP) :: abab, acac
    real(RP) :: bcbc, baba
    real(RP) :: caca, cbcb

    real(RP), parameter :: eps = 1.E-16_RP
    real(RP) :: pi
    !---------------------------------------------------------------------------

    pi = atan(1.0_RP) * 4.0_RP

    area = 0.0_RP

    if(trim(polygon_type)=='ON_PLANE') then
       !
       !---- Note : On a plane,
       !----        area = | ourter product of two vectors |.
       !
       call MISC_3dvec_cross( abc(:), a(:), b(:), a(:), c(:) )
       call MISC_3dvec_abs( prd, abc(:) )
       call MISC_3dvec_abs( r  , a(:)   )

       prd = 0.5_RP * prd !! triangle area
       if ( r < eps * radius ) then
          print *, "zero length?", a(:)
       else
          r = 1.0_RP / r   !! 1 / length
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
    real(RP), intent(out) :: p(3) ! intersection point
    real(RP), intent(in ) :: a(3), b(3), c(3), d(3)

    real(RP), parameter :: o(3) = 0.0_RP

    real(RP)            :: oaob(3), ocod(3), cdab(3)
    real(RP)            :: ip, length
    real(RP)            :: angle_aop, angle_pob, angle_aob
    real(RP)            :: angle_cop, angle_pod, angle_cod

    real(RP), parameter :: eps = 1.E-12_RP
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
       p(:) = 0.0_RP
    endif

    return
  end subroutine MISC_3dvec_intersec

  !---------------------------------------------------------------------
  subroutine MISC_3dvec_anticlockwise( vertex, nvert )
    ! bubble sort anticlockwise by angle
    implicit none

    integer, intent(in)    :: nvert
    real(RP), intent(inout) :: vertex(nvert,3)

    real(RP), parameter :: o(3) = 0.0_RP
    real(RP)            :: v1(3), v2(3), v3(3)
    real(RP)            :: xp(3), ip
    real(RP)            :: angle1, angle2

    real(RP), parameter :: eps = 1.E-12_RP

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

    if (       abs(ip)                 < eps    &      ! on the same line
         .AND. abs(angle2)-abs(angle1) < 0.0_RP ) then ! which is far?
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

    if (       abs(ip)                 < eps    &      ! on the same line
         .AND. abs(angle2)-abs(angle1) < 0.0_RP ) then ! which is far?
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

    real(RP),intent(inout) :: x, y, z
    real(RP),intent(in)    :: lat, lon
    real(RP),intent(in)    :: radius
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

    real(RP), intent(in)  :: r          ! radius in meter
    real(RP), intent(in)  :: lon1, lat1 ! in radian
    real(RP), intent(in)  :: lon2, lat2 ! in radian
    real(RP), intent(out) :: dist       ! distance of the two points in meter

    real(RP) :: gmm, gno_x, gno_y
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


