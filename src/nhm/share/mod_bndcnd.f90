!-------------------------------------------------------------------------------
!>
!! Boundary conditions module
!!
!! @par Description
!!         This module provides the subroutines for boundary conditions.
!!
!! @author H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita) Imported from igdc-4.33
!! @li      2011-07-22 (T.Ohno)   Add subroutines for plane hgrid systems.
!!
!<
module mod_bndcnd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: BNDCND_setup

  public :: BNDCND_all
  public :: BNDCND_thermo
  public :: BNDCND_rhovxvyvz
  public :: BNDCND_rhow

  public  :: BNDCND_all_plane   ! [add] T.Ohno 110722
  public  :: BNDCND_rhow_plane  ! [add] T.Ohno 110722
  public  :: BNDCND_rhov2_plane ! [add] T.Ohno 110722
  private :: BNDCND_w_plane     ! [add] T.Ohno 110722

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  !--- Vertical boundary condition for momentum at the top
  character(len=ADM_NSYS), private, save :: BND_TYPE_M_TOP    != 'RIGID' : rigid surface
                                                              != 'FREE'  : free surface

  !--- Vertical boundary condition for momentum at the ground
  character(len=ADM_NSYS), private, save :: BND_TYPE_M_BOTTOM != 'RIGID' : rigid surface
                                                              != 'FREE'  : free surface

  !--- Vertical boundary condition for temperature at the top
  character(len=ADM_NSYS), private, save :: BND_TYPE_T_TOP    != 'TEM' : tem(kmax+1) = tem(kmax)
                                                              != 'EPL' : lagrange extrapolation

  !--- Vertical boundary condition for temperature at the ground
  character(len=ADM_NSYS), private, save :: BND_TYPE_T_BOTTOM != 'FIX' : tems fix
                                                              != 'TEM' : tem(kmin-1) = tem(kmin)
                                                              != 'EPL' : lagrange extrapolation

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !>
  !> Description of the subroutine CNST_setup
  !>
  subroutine BNDCND_setup
    use mod_adm, only: &
       ADM_CTL_FID, &
       ADM_proc_stop
    use mod_cnst, only: &
       CNST_TEMS0
    implicit none

    namelist / BNDCNDPARAM / &
         BND_TYPE_M_TOP,    &
         BND_TYPE_M_BOTTOM, &
         BND_TYPE_T_TOP,    &
         BND_TYPE_T_BOTTOM

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- set default
    BND_TYPE_M_TOP    = 'FREE'
    BND_TYPE_M_BOTTOM = 'RIGID'
    BND_TYPE_T_TOP    = 'TEM'
    BND_TYPE_T_BOTTOM = 'TEM'

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[bndcnd]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=BNDCNDPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** BNDCNDPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,BNDCNDPARAM)

    if    ( BND_TYPE_M_TOP == 'RIGID' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    top   ) : rigid'
    elseif( BND_TYPE_M_TOP == 'FREE' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    top   ) : free'
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_M_TOP. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_M_BOTTOM == 'RIGID' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    bottom) : rigid'
    elseif( BND_TYPE_M_BOTTOM == 'FREE' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    bottom) : free'
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_M_BOTTOM. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_T_TOP == 'TEM' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, top   ) : equal to uppermost atmosphere'
    elseif( BND_TYPE_T_TOP == 'EPL' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, top   ) : lagrange extrapolation'
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_T_TOP. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_T_BOTTOM == 'FIX' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, bottom) : fixed'
       write(ADM_LOG_FID,*) '***           boundary temperature (CNST_TEMS0) : ', CNST_TEMS0
    elseif( BND_TYPE_T_BOTTOM == 'TEM' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, bottom) : equal to lowermost atmosphere'
    elseif( BND_TYPE_T_BOTTOM == 'EPL' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, bottom) : lagrange extrapolation'
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_T_BOTTOM. STOP.'
       call ADM_proc_stop
    endif

    return
  end subroutine BNDCND_setup

  !-----------------------------------------------------------------------------
  !------
  !------ Boundary condition setting for thermodynamical variables
  !------    1. calculation region : (:,[ADM_kmin-1,kmax+1],:)
  !------    2. boundary types are controled in the sub[BNDCND_setup].
  !------
  subroutine BNDCND_thermo( &
       ijdim, &
       tem,   &
       rho,   &
       pre,   &
       phi    )                   
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_RAIR,  &
       CNST_EGRAV, &
       CNST_TEMS0
    implicit none

    integer, intent(in)    :: ijdim           ! number of horizontal grid
    real(8), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(8), intent(inout) :: rho(ijdim,kdim) ! density  
    real(8), intent(inout) :: pre(ijdim,kdim) ! pressure
    real(8), intent(in)    :: phi(ijdim,kdim) ! geopotential

    integer :: ij

    real(8) :: z,z1,z2,z3,p1,p2,p3
    real(8) :: lag_intpl
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3))/((z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3))/((z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2))/((z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------

    !--- set the TOP boundary of temperature
    select case( trim(BND_TYPE_T_TOP) )
    case('TEM') 
       tem(:,kmax+1) = tem(:,kmax) ! dT/dz = 0
    case('EPL') 
       do ij = 1, ijdim
          z  = phi(ij,kmax+1) / CNST_EGRAV
          z1 = phi(ij,kmax  ) / CNST_EGRAV
          z2 = phi(ij,kmax-1) / CNST_EGRAV
          z3 = phi(ij,kmax-2) / CNST_EGRAV

          tem(ij,kmax+1) = lag_intpl( z ,                 &
                                      z1, tem(ij,kmax  ), &
                                      z2, tem(ij,kmax-1), &
                                      z3, tem(ij,kmax-2)  )
       enddo
    endselect

    !--- set the BOTTOM boundary of temperature
    select case( trim(BND_TYPE_T_BOTTOM) )
    case('FIX')
       tem(:,kmin-1) = CNST_TEMS0
    case('TEM')
       tem(:,kmin-1) = tem(:,kmin) ! dT/dz = 0
    case('EPL') 
       do ij = 1, ijdim
          z1 = phi(ij,kmin+2) / CNST_EGRAV
          z2 = phi(ij,kmin+1) / CNST_EGRAV
          z3 = phi(ij,kmin  ) / CNST_EGRAV
          z  = phi(ij,kmin-1) / CNST_EGRAV

          tem(ij,kmin-1) = lag_intpl( z,                  &
                                      z1, tem(ij,kmin+2), &
                                      z2, tem(ij,kmin+1), &
                                      z3, tem(ij,kmin  )  )
       enddo
    endselect

    do ij = 1, ijdim

       !--- set the boundary of pressure ( hydrostatic balance )
       pre(ij,kmax+1) = pre(ij,kmax-1) - rho(ij,kmax) * ( phi(ij,kmax+1) - phi(ij,kmax-1) )
       pre(ij,kmin-1) = pre(ij,kmin+1) - rho(ij,kmin) * ( phi(ij,kmin-1) - phi(ij,kmin+1) )

       !--- set the boundary of density
       rho(ij,kmax+1) = pre(ij,kmax+1) / ( CNST_RAIR * tem(ij,kmax+1) )
       rho(ij,kmin-1) = pre(ij,kmin-1) / ( CNST_RAIR * tem(ij,kmin-1) )

    enddo

    return
  end subroutine BNDCND_thermo

  !-----------------------------------------------------------------------------
  !------ Boundary condition setting for rhogvx
  subroutine BNDCND_rhovxvyvz( &
       ijdim,  &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    implicit none

    integer, intent(in)    :: ijdim
    real(8), intent(in)    :: rhog  (ijdim,kdim)
    real(8), intent(inout) :: rhogvx(ijdim,kdim)
    real(8), intent(inout) :: rhogvy(ijdim,kdim)
    real(8), intent(inout) :: rhogvz(ijdim,kdim)

    integer :: ij
    !---------------------------------------------------------------------------

    !------ top ( rhogvx, rhogvy, rhogvz )
    select case( trim(BND_TYPE_M_TOP) )
    case('RIGID')
       do ij = 1, ijdim
          rhogvx(ij,kmax+1) = -rhogvx(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvy(ij,kmax+1) = -rhogvy(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvz(ij,kmax+1) = -rhogvz(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
       enddo
    case('FREE')
       do ij = 1, ijdim
          rhogvx(ij,kmax+1) =  rhogvx(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvy(ij,kmax+1) =  rhogvy(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvz(ij,kmax+1) =  rhogvz(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
       enddo
    endselect

    !------ bottom ( rhogvx, rhogvy, rhogvz )
    select case( trim(BND_TYPE_M_BOTTOM) )
    case('RIGID')
       do ij = 1, ijdim
          rhogvx(ij,kmin-1) = -rhogvx(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvy(ij,kmin-1) = -rhogvy(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvz(ij,kmin-1) = -rhogvz(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
       enddo
    case('FREE')
       do ij = 1, ijdim
          rhogvx(ij,kmin-1) =  rhogvx(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvy(ij,kmin-1) =  rhogvy(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvz(ij,kmin-1) =  rhogvz(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
       enddo
    endselect

    return
  end subroutine BNDCND_rhovxvyvz

  !-----------------------------------------------------------------------------
  !------
  !------ Boundary condition setting for rhogw
  !------
  subroutine BNDCND_rhow( &
       ijdim,   &
       rhogvx,  &
       rhogvy,  &
       rhogvz,  &
       rhogw,   &
       c2wfact  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax, &
       ADM_VMISS
    implicit none

    integer, intent(in)    :: ijdim
    real(8), intent(in)    :: rhogvx (ijdim,kdim)
    real(8), intent(in)    :: rhogvy (ijdim,kdim)
    real(8), intent(in)    :: rhogvz (ijdim,kdim)
    real(8), intent(inout) :: rhogw  (ijdim,kdim)
    real(8), intent(in)    :: c2wfact(6,ijdim,kdim)

    integer :: ij, k
    !---------------------------------------------------------------------------

    select case( trim(BND_TYPE_M_TOP) )
    case('RIGID') ! rhow / G^{1/2} = 0.D0
       rhogw(:,kmax+1) = 0.D0
    case('FREE')
       k = kmax+1
       do ij = 1, ijdim
          rhogw(ij,k) = - ( c2wfact(1,ij,k) * rhogvx(ij,k  ) &
                          + c2wfact(2,ij,k) * rhogvx(ij,k-1) &
                          + c2wfact(3,ij,k) * rhogvy(ij,k  ) &
                          + c2wfact(4,ij,k) * rhogvy(ij,k-1) &
                          + c2wfact(5,ij,k) * rhogvz(ij,k  ) &
                          + c2wfact(6,ij,k) * rhogvz(ij,k-1) )
       enddo
    endselect

    select case( trim(BND_TYPE_M_BOTTOM) )
    case('RIGID') ! rhow / G^{1/2} = 0.D0
       rhogw(:,kmin  ) = 0.D0
    case('FREE')
       k = kmin
       do ij = 1, ijdim
          rhogw(ij,k) = - ( c2wfact(1,ij,k) * rhogvx(ij,k  ) &
                          + c2wfact(2,ij,k) * rhogvx(ij,k-1) &
                          + c2wfact(3,ij,k) * rhogvy(ij,k  ) &
                          + c2wfact(4,ij,k) * rhogvy(ij,k-1) &
                          + c2wfact(5,ij,k) * rhogvz(ij,k  ) &
                          + c2wfact(6,ij,k) * rhogvz(ij,k-1) )
       enddo
    endselect

    rhogw(:,1:kmin-1) = ADM_VMISS

    return
  end subroutine BNDCND_rhow

  !-----------------------------------------------------------------------------
  !------
  !------ Boundary condition setting for all variables.
  !------    1. calculation region (:,[kmin,kmax+1],:) 
  !------       for  rhogw & w.
  !------    2. calculation region (:,[kmin-1,kmax+1],:) 
  !------       for  the other variables.
  !------
  subroutine BNDCND_all( &
       ijdim,      &
       rho,        &
       vx,         &
       vy,         &
       vz,         &
       w,          &
       ein,        &
       tem,        &
       pre,        &
       rhog,       &
       rhogvx,     &
       rhogvy,     &
       rhogvz,     &
       rhogw,      &
       rhoge,      &
       gsqrtgam2,  &
       gsqrtgam2h, &
       phi,        &
       c2wfact     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax, &
       ADM_VMISS
    use mod_cnst, only: &
       CNST_CV
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    implicit none

    integer, intent(in)    :: ijdim
    real(8), intent(inout) :: rho(ijdim,kdim)
    real(8), intent(inout) :: vx (ijdim,kdim)
    real(8), intent(inout) :: vy (ijdim,kdim)
    real(8), intent(inout) :: vz (ijdim,kdim)
    real(8), intent(inout) :: w  (ijdim,kdim)
    real(8), intent(inout) :: ein(ijdim,kdim)
    real(8), intent(inout) :: tem(ijdim,kdim)
    real(8), intent(inout) :: pre(ijdim,kdim)

    real(8), intent(inout) :: rhog  (ijdim,kdim)
    real(8), intent(inout) :: rhogvx(ijdim,kdim)
    real(8), intent(inout) :: rhogvy(ijdim,kdim)
    real(8), intent(inout) :: rhogvz(ijdim,kdim)
    real(8), intent(inout) :: rhogw (ijdim,kdim)
    real(8), intent(inout) :: rhoge (ijdim,kdim)

    real(8), intent(in)    :: gsqrtgam2 (ijdim,kdim)
    real(8), intent(in)    :: gsqrtgam2h(ijdim,kdim)
    real(8), intent(in)    :: phi       (ijdim,kdim)
    real(8), intent(in)    :: c2wfact   (ijdim,kdim)

    integer :: ij, k
    !---------------------------------------------------------------------------

    !
    !--- Thermodynamical variables ( tem, th, rho, pre, rhoge )
    !
    call BNDCND_thermo( ijdim, & !--- [IN]
                        tem,   & !--- [INOUT]
                        rho,   & !--- [INOUT]
                        pre,   & !--- [INOUT]
                        phi    ) !--- [IN]

    rhog(:,kmax+1) = rho(:,kmax+1) * gsqrtgam2(:,kmax+1)
    rhog(:,kmin-1) = rho(:,kmin-1) * gsqrtgam2(:,kmin-1)

    !
    !--- internal energy ( ein, rhoge ), q = 0 at boundary
    !
    ein  (:,kmax+1) = CNST_CV * tem(:,kmax+1)
    ein  (:,kmin-1) = CNST_CV * tem(:,kmin-1)

    rhoge(:,kmax+1) = rhog(:,kmax+1) * ein(:,kmax+1)
    rhoge(:,kmin-1) = rhog(:,kmin-1) * ein(:,kmin-1)

    !
    !--- Momentum ( rhogvx, rhogvy, rhogvz, vx, vy, vz )
    !
    call BNDCND_rhovxvyvz( ijdim,  & !--- [IN]
                           rhog,   & !--- [IN]
                           rhogvx, & !--- [INOUT]
                           rhogvy, & !--- [INOUT]
                           rhogvz  ) !--- [INOUT]

    vx(:,kmax+1) = rhogvx(:,kmax+1) / rhog(:,kmax+1)
    vy(:,kmax+1) = rhogvy(:,kmax+1) / rhog(:,kmax+1)
    vz(:,kmax+1) = rhogvz(:,kmax+1) / rhog(:,kmax+1)

    vx(:,kmin-1) = rhogvx(:,kmin-1) / rhog(:,kmin-1)
    vy(:,kmin-1) = rhogvy(:,kmin-1) / rhog(:,kmin-1)
    vz(:,kmin-1) = rhogvz(:,kmin-1) / rhog(:,kmin-1)

    !
    !--- Momentum ( rhogw, w )
    !
    call BNDCND_rhow( ijdim,   & !--- [IN]
                      rhogvx,  & !--- [IN]
                      rhogvy,  & !--- [IN]
                      rhogvz,  & !--- [IN]
                      rhogw,   & !--- [INOUT]
                      c2wfact  ) !--- [IN]

    k = kmax+1
    do ij = 1, ijdim
       w(ij,k) = rhogw(ij,k) / ( gsqrtgam2h(ij,k) * 0.5D0 * ( GRD_afac(k) * rho(ij,k  ) &
                                                            + GRD_bfac(k) * rho(ij,k-1) ) )
    enddo

    k = kmin
    do ij = 1, ijdim
       w(ij,k) = rhogw(ij,k) / ( gsqrtgam2h(ij,k) * 0.5D0 * ( GRD_afac(k) * rho(ij,k  ) &
                                                            + GRD_bfac(k) * rho(ij,k-1) ) )
    enddo

    w(:,1:kmin-1) = ADM_VMISS

    return
  end subroutine BNDCND_all

  !-----------------------------------------------------------------------------
  subroutine BNDCND_rhow_plane(&
       ijdim,            & !--- IN : number of horizontal grid
       rhogvx,           & !--- IN : rho*Vx   ( gam2 X G^{1/2} )
       rhogvy,           & !--- IN : rho*Vy   ( gam2 X G^{1/2} )
       rhogw,            & !--- INOUT : rho*w ( gam2 X G^{1/2} )
       gsqrtgam2,        & !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,       & !--- IN : G^{1/2} at the cell wall
       g3xh,             & !--- IN : G3X at the cell wall
       g3yh,             & !--- IN : G3Y at the cell wall
       g3zh              & !--- IN : G3Z at the cell wall
       )
    ! [add] T.Ohno 110722
    !------
    !------ Boundary condition setting for rhow only.
    !------    1. calculation region (:,[kmin,kmax+1],:)
    !------
    !
    use mod_adm, only :  &
         kdim => ADM_kall,       &
         kmin => ADM_kmin,       &
         kmax => ADM_kmax,       &
         ADM_VMISS
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac
    implicit none

    integer, intent(in)    :: ijdim
    real(8), intent(in) :: rhogvx(ijdim,kdim)
    real(8), intent(in) :: rhogvy(ijdim,kdim)
    real(8), intent(in) :: gsqrtgam2(ijdim,kdim)
    real(8), intent(in) :: gsqrtgam2h(ijdim,kdim)
    real(8), intent(in) :: g3xh(ijdim,kdim)
    real(8), intent(in) :: g3yh(ijdim,kdim)
    real(8), intent(in) :: g3zh(ijdim,kdim)
    !
    real(8), intent(inout) :: rhogw(ijdim,kdim)
    !
    integer :: k
    !
    select case(trim(BND_TYPE_M_BOTTOM))
    case('FREE')
       !
       !--- rhow/G^{1/2} + G3X*rhovx + G3Y*rhovy + G3Z*rhovz = 0
       rhogw(:,1:kmin-1) = ADM_VMISS
       do k = kmin, kmax+1, (kmax-kmin+1)
          rhogw(:,k) = -(                                       &
               +( ( GRD_afac(k)/gsqrtgam2(:,k  )*rhogvx(:,k  )  &
               +GRD_bfac(k)/gsqrtgam2(:,k-1)*rhogvx(:,k-1) )    &
               * 0.5D0*gsqrtgam2h(:,k)*g3xh(:,k)                &
               +( GRD_afac(k)/gsqrtgam2(:,k  )*rhogvy(:,k  )    &
               +GRD_bfac(k)/gsqrtgam2(:,k-1)*rhogvy(:,k-1) )    &
               * 0.5D0**gsqrtgam2h(:,k)*g3yh(:,k) ))            & 
               * gsqrtgam2h(:,k) 
       end do
       !
    case('RIGID')
       !
       !--- rhow/G^{1/2} =0.0D0
       rhogw(:,1:kmin-1) = ADM_VMISS
       do k = kmin, kmax+1, (kmax-kmin+1)
          rhogw(:,k) = 0.0D0
       end do
       !
    end select
    !
  end subroutine BNDCND_rhow_plane
  !-----------------------------------------------------------------------------
  subroutine BNDCND_rhov2_plane( &
       ijdim,              &  !--- IN : number of horizontal grid
       rhog,               &  !--- IN :    rho     ( gam2 X G^{1/2} )
       rhogvx,             &  !--- INOUT : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy              &  !--- INOUT : rho*Vy  ( gam2 X G^{1/2} )
       )
    ! [add] T.Ohno 110722
    !------
    !------ Boundary condition setting for rhogvx
    !------    * calculation region (:,[kmin-1,kmax+1],:) 
    !------
    !
    use mod_adm, only :  &
         ADM_LOG_FID,    &
         kdim => ADM_kall,       &
         kmin => ADM_kmin,       &
         kmax => ADM_kmax
    implicit none
    !
    integer, intent(in)    :: ijdim
    real(8), intent(in) :: rhog(ijdim,kdim)
    real(8), intent(inout) :: rhogvx(ijdim,kdim)
    real(8), intent(inout) :: rhogvy(ijdim,kdim)
    !
    !
    !------ bottom ( rhogvx, rhogvy, rhogvz )
    select case(trim(BND_TYPE_M_BOTTOM))
    case('FREE')
       rhogvx(:,kmin-1)                    &
            = rhogvx(:,kmin)               &
            / rhog(:,kmin)                 &
            * rhog(:,kmin-1)
       rhogvy(:,kmin-1)                    &
            = rhogvy(:,kmin)               &
            / rhog(:,kmin)                 &
            * rhog(:,kmin-1)
    case('RIGID')
       rhogvx(:,kmin-1)                    &
            = -rhogvx(:,kmin)              &
            / rhog(:,kmin)                 &
            * rhog(:,kmin-1)
       rhogvy(:,kmin-1)                    &
            = -rhogvy(:,kmin)              &
            / rhog(:,kmin)                 &
            * rhog(:,kmin-1)
       !
    case default
       write(ADM_LOG_FID,*) &
            'Msg : Sub[BNDCND_all]/Mod[bndcnd]'
       write(ADM_LOG_FID,*) &
            ' **** Warning : invalid t_top_type',BND_TYPE_M_BOTTOM
    end select
    !
    !------ top ( rhogvx, rhogvy, rhogvz ) : stress free
    rhogvx(:,kmax+1)                    &
         = rhogvx(:,kmax)               &
         / rhog(:,kmax)                 &
         * rhog(:,kmax+1)
    rhogvy(:,kmax+1)                    &
         = rhogvy(:,kmax)               &
         / rhog(:,kmax)                 &
         * rhog(:,kmax+1)
    !
    return
    !
  end subroutine BNDCND_rhov2_plane
  !-----------------------------------------------------------------------------
  subroutine BNDCND_all_plane(&
       ijdim,           &  !--- IN : number of horizontal grid
       vx,              &  !--- INOUT : Vx 
       vy,              &  !--- INOUT : Vy 
       w,               &  !--- INOUT : w 
       tem,             &  !--- INOUT : temp. 
       rho,             &  !--- INOUT : density
       pre,             &  !--- INOUT : pressure 
       ein,             &  !--- INOUT : internal energy
       rhog,            &  !--- INOUT : rho  ( gam2 X G^{1/2} )
       rhogvx,          &  !--- INOUT : rho*Vx  ( gam2 X G^{1/2} )
       rhogvy,          &  !--- INOUT : rho*Vy  ( gam2 X G^{1/2} )
       rhogw,           &  !--- INOUT : rho*w   ( gam2 X G^{1/2} )
       rhoge,           &  !--- INOUT : rho*ein ( gam2 X G^{1/2} )
       phi,             &  !--- IN : geopotential
       gsqrtgam2,       &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,      &  !--- IN : G^{1/2} at the cell wall
       gzxh,            &  !--- IN : GZX at the cell wall
       gzyh             &  !--- IN : GZY at the cell wall
       )
    ! [add] T.Ohno 110722
    !------
    !------ Boundary condition setting for all variables.
    !------    1. calculation region (:,[kmin,kmax+1],:) 
    !------       for  rhogw & w.
    !------    2. calculation region (:,[kmin-1,kmax+1],:) 
    !------       for  the other variables.
    !------
    !
    use mod_adm, only :  &
         ADM_LOG_FID,    &
         kdim => ADM_kall,       &
         kmin => ADM_kmin,       &
         kmax => ADM_kmax
    use mod_cnst, only : &
         CNST_CV
    !
    implicit none
    integer, intent(in)    :: ijdim
    real(8), intent(inout) :: vx(ijdim,kdim)
    real(8), intent(inout) :: vy(ijdim,kdim)
    real(8), intent(inout) :: w(ijdim,kdim)
    !
    real(8), intent(inout) :: tem(ijdim,kdim)
    real(8), intent(inout) :: rho(ijdim,kdim)
    real(8), intent(inout) :: pre(ijdim,kdim)
    real(8), intent(inout) :: ein(ijdim,kdim)
    !
    real(8), intent(inout) :: rhog(ijdim,kdim)
    real(8), intent(inout) :: rhogvx(ijdim,kdim)
    real(8), intent(inout) :: rhogvy(ijdim,kdim)
    real(8), intent(inout) :: rhogw(ijdim,kdim)
    real(8), intent(inout) :: rhoge(ijdim,kdim)
    !
    real(8), intent(in) :: phi(ijdim,kdim)
    real(8), intent(in) :: gsqrtgam2(ijdim,kdim)
    real(8), intent(in) :: gsqrtgam2h(ijdim,kdim)
    real(8), intent(in) :: gzxh(ijdim,kdim)
    real(8), intent(in) :: gzyh(ijdim,kdim)
    !
    !
    !--- Thermodynamical variables ( tem, th, rho, pre )
    call BNDCND_thermo( &
         ijdim,         & !-- in
         tem,           & !-- inout
         rho,           & !-- inout
         pre,           & !-- inout
         phi )            !-- in
    rhog(:,kmin-1) = rho(:,kmin-1)*gsqrtgam2(:,kmin-1)
    rhog(:,kmax+1) = rho(:,kmax+1)*gsqrtgam2(:,kmax+1)
    !
    !--- Momentum ( vx, vy, vz, rhogvx, rhogvy, rhogvz, w, rhogw )
    !
    !------ bottom ( vx, vy, vz )
    select case(trim(BND_TYPE_M_BOTTOM))
    case('FREE')
       vx(:,kmin-1)    = vx(:,kmin)
       vy(:,kmin-1)    = vy(:,kmin)
    case('RIGID')
       vx(:,kmin-1)    = -vx(:,kmin)
       vy(:,kmin-1)    = -vy(:,kmin)
    case default
       write(ADM_LOG_FID,*) &
            'Msg : Sub[BNDCND_all]/Mod[bndcnd]'
       write(ADM_LOG_FID,*) &
            ' **** Warning : invalid t_top_type',BND_TYPE_M_BOTTOM
    end select
    !
    !------ top ( vx, vy, vz ) : stress free
    vx(:,kmax+1)    =  vx(:,kmax)
    vy(:,kmax+1)    =  vy(:,kmax)
    !
    !------ top & bottom ( rhogvx, rhogvy, rhogvz )
    rhogvx(:,kmin-1)                  &
         = rho(:,kmin-1)*vx(:,kmin-1) &
         * gsqrtgam2(:,kmin-1)
    rhogvy(:,kmin-1)                  &
         = rho(:,kmin-1)*vy(:,kmin-1) &
         * gsqrtgam2(:,kmin-1)
    rhogvx(:,kmax+1)                  &
         = rho(:,kmax+1)*vx(:,kmax+1) &
         * gsqrtgam2(:,kmax+1)
    rhogvy(:,kmax+1)                  &
         = rho(:,kmax+1)*vy(:,kmax+1) &
         * gsqrtgam2(:,kmax+1)
    !
    !------- top & bottom ( w, rhogw )
    call BNDCND_w_plane(  &
         ijdim,     & !--- IN
         rho,       & !--- IN
         rhogvx,    & !--- IN
         rhogvy,    & !--- IN
         rhogw,     & !--- INOUT
         w,         & !--- INOUT
         gsqrtgam2, & !--- IN
         gsqrtgam2h,& !--- IN
         gzxh,      & !--- IN
         gzyh       & !--- IN
         )
    !
    !--- internal energy ( ein, rhoge )
    ein(:,kmin-1)=CNST_CV*tem(:,kmin-1)
    ein(:,kmax+1)=CNST_CV*tem(:,kmax+1)
    rhoge(:,kmin-1)                      &
         = rho(:,kmin-1) * ein(:,kmin-1) &
         * gsqrtgam2(:,kmin-1)
    rhoge(:,kmax+1)                      &
         = rho(:,kmax+1) * ein(:,kmax+1) &
         * gsqrtgam2(:,kmax+1)
    !
    return
    !
  end subroutine BNDCND_all_plane
  !-----------------------------------------------------------------------------
  subroutine BNDCND_w_plane(  &
       ijdim,           &  !--- IN : number of horizontal grid
       rho,             & !--- IN : rho      ( physical )
       rhogvx,          & !--- IN : rho*Vx   ( gam2 X G^{1/2} )
       rhogvy,          & !--- IN : rho*Vy   ( gam2 X G^{1/2} )
       rhogw,           & !--- INOUT : rho*w ( gam2 X G^{1/2} )
       w,               & !--- INOUT : w     ( physical )
       gsqrtgam2,       & !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,      & !--- IN : G^{1/2} at the cell wall
       g3xh,            & !--- IN : G3X at the cell wall
       g3yh             & !--- IN : G3Y at the cell wall
       )
    ! [add] T.Ohno 110722
    !------
    !------ Boundary condition setting for rhogw and w.
    !------    1. calculation region (:,[kmin,kmax+1],:)
    !------
    !
    use mod_adm, only :  &
         kdim => ADM_kall,       &
         kmin => ADM_kmin,       &
         kmax => ADM_kmax,       &
         ADM_VMISS
    use mod_grd, only :  &
         GRD_afac,       &
         GRD_bfac
    !
    implicit none
    !
    integer, intent(in)    :: ijdim
    real(8), intent(in) :: rho(ijdim,kdim)
    real(8), intent(in) :: rhogvx(ijdim,kdim)
    real(8), intent(in) :: rhogvy(ijdim,kdim)
    real(8), intent(in) :: gsqrtgam2(ijdim,kdim)
    real(8), intent(in) :: gsqrtgam2h(ijdim,kdim)
    real(8), intent(in) :: g3xh(ijdim,kdim)
    real(8), intent(in) :: g3yh(ijdim,kdim)
    !
    real(8), intent(inout) :: rhogw(ijdim,kdim)
    real(8), intent(inout) :: w(ijdim,kdim)
    !
    integer :: k
    !
    select case(trim(BND_TYPE_M_BOTTOM))
    case('FREE')
       !
       !--- rhow/G^{1/2} + G3X*rhovx + G3Y*rhovy + G3Z*rhovz = 0
       rhogw(:,1:kmin-1) = ADM_VMISS
       w(:,1:kmin-1) = ADM_VMISS
       do k = kmin, kmax+1, (kmax-kmin+1)
          rhogw(:,k) = -(                                        &
               +( ( GRD_afac(k)/gsqrtgam2(:,k  )*rhogvx(:,k  )   &
               +GRD_bfac(k)/gsqrtgam2(:,k-1)*rhogvx(:,k-1) )     &
               * 0.5D0*gsqrtgam2h(:,k)*g3xh(:,k)                 &
               +( GRD_afac(k)/gsqrtgam2(:,k  )*rhogvy(:,k  )     &
               +GRD_bfac(k)/gsqrtgam2(:,k-1)*rhogvy(:,k-1) )     &
               * 0.5D0**gsqrtgam2h(:,k)*g3yh(:,k) ))             & 
               * gsqrtgam2h(:,k) 
          w(:,k)                                                 &
               = rhogw(:,k)                                      &
               / gsqrtgam2h(:,k)                                 &
               /( ( GRD_afac(k) * rho(:,k)                       &
               + GRD_bfac(k) * rho(:,k-1) ) *0.5D0 )
       end do
       !
    case('RIGID')
       !
       !--- rhow/G^{1/2} =0.0D0
       rhogw(:,1:kmin-1) = ADM_VMISS
       w(:,1:kmin-1) = ADM_VMISS
       do k = kmin, kmax+1, (kmax-kmin+1)
          rhogw(:,k) = 0.0D0
          w(:,k)    = 0.0D0
       end do
       !
    end select
    !
  end subroutine BNDCND_w_plane
  !-----------------------------------------------------------------------------
end module mod_bndcnd
!-------------------------------------------------------------------------------
