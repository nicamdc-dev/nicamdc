!-------------------------------------------------------------------------------
!> Module boundary conditions
!!
!! @par Description
!!          This module provides the subroutines for boundary conditions
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_bndcnd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: BNDCND_setup
  public :: BNDCND_all
  public :: BNDCND_thermo
  public :: BNDCND_rhovxvyvz
  public :: BNDCND_rhow

  interface BNDCND_thermo
     module procedure BNDCND_thermo_SP
     module procedure BNDCND_thermo_DP
  end interface BNDCND_thermo

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
  !--- Vertical boundary condition for temperature at the top
  character(len=H_SHORT), private :: BND_TYPE_T_TOP    != 'TEM' : tem(kmax+1) = tem(kmax) (default)
                                                       != 'EPL' : lagrange extrapolation

  !--- Vertical boundary condition for temperature at the ground
  character(len=H_SHORT), private :: BND_TYPE_T_BOTTOM != 'TEM' : tem(kmin-1) = tem(kmin) (default)
                                                       != 'EPL' : lagrange extrapolation

  !--- Vertical boundary condition for momentum at the top
  character(len=H_SHORT), private :: BND_TYPE_M_TOP    != 'RIGID' : rigid surface
                                                       != 'FREE'  : free surface (default)

  !--- Vertical boundary condition for momentum at the ground
  character(len=H_SHORT), private :: BND_TYPE_M_BOTTOM != 'RIGID' : rigid surface (default)
                                                       != 'FREE'  : free surface

  logical, private :: is_top_tem   = .false.
  logical, private :: is_top_epl   = .false.
  logical, private :: is_btm_tem   = .false.
  logical, private :: is_btm_epl   = .false.
  logical, private :: is_top_rigid = .false.
  logical, private :: is_top_free  = .false.
  logical, private :: is_btm_rigid = .false.
  logical, private :: is_btm_free  = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine BNDCND_setup
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / BNDCNDPARAM / &
       BND_TYPE_T_TOP,    &
       BND_TYPE_T_BOTTOM, &
       BND_TYPE_M_TOP,    &
       BND_TYPE_M_BOTTOM

    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- set default
    BND_TYPE_T_TOP    = 'TEM'
    BND_TYPE_T_BOTTOM = 'TEM'
    BND_TYPE_M_TOP    = 'FREE'
    BND_TYPE_M_BOTTOM = 'RIGID'

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[bndcnd]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=BNDCNDPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** BNDCNDPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=BNDCNDPARAM)

    if    ( BND_TYPE_T_TOP == 'TEM' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (temperature, top   ) : equal to uppermost atmosphere'
       is_top_tem = .true.
    elseif( BND_TYPE_T_TOP == 'EPL' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (temperature, top   ) : lagrange extrapolation'
       is_top_epl = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_T_TOP. STOP.'
       call PRC_MPIstop
    endif

    if    ( BND_TYPE_T_BOTTOM == 'TEM' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (temperature, bottom) : equal to lowermost atmosphere'
       is_btm_tem = .true.
    elseif( BND_TYPE_T_BOTTOM == 'EPL' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (temperature, bottom) : lagrange extrapolation'
       is_btm_epl = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_T_BOTTOM. STOP.'
       call PRC_MPIstop
    endif

    if    ( BND_TYPE_M_TOP == 'RIGID' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (momentum,    top   ) : rigid'
       is_top_rigid = .true.
    elseif( BND_TYPE_M_TOP == 'FREE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (momentum,    top   ) : free'
       is_top_free  = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_M_TOP. STOP.'
       call PRC_MPIstop
    endif

    if    ( BND_TYPE_M_BOTTOM == 'RIGID' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (momentum,    bottom) : rigid'
       is_btm_rigid = .true.
    elseif( BND_TYPE_M_BOTTOM == 'FREE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Boundary setting type (momentum,    bottom) : free'
       is_btm_free  = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_M_BOTTOM. STOP.'
       call PRC_MPIstop
    endif

    return
  end subroutine BNDCND_setup

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for all prognostic variables.
  subroutine BNDCND_all( &
       ijdim,      &
       kdim,       &
       ldim,       &
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
       phi,        &
       c2wfact,    &
       c2wfact_Gz  )
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    use mod_const, only: &
       CONST_CVdry
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: ldim
    real(RP), intent(inout) :: rho       (ijdim,kdim,ldim) ! density
    real(RP), intent(inout) :: vx        (ijdim,kdim,ldim) ! horizontal wind (x)
    real(RP), intent(inout) :: vy        (ijdim,kdim,ldim) ! horizontal wind (y)
    real(RP), intent(inout) :: vz        (ijdim,kdim,ldim) ! horizontal wind (z)
    real(RP), intent(inout) :: w         (ijdim,kdim,ldim) ! vertical   wind
    real(RP), intent(inout) :: ein       (ijdim,kdim,ldim) ! internal energy
    real(RP), intent(inout) :: tem       (ijdim,kdim,ldim) ! temperature
    real(RP), intent(inout) :: pre       (ijdim,kdim,ldim) ! pressure
    real(RP), intent(inout) :: rhog      (ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogvx    (ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogvy    (ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogvz    (ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogw     (ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhoge     (ijdim,kdim,ldim)
    real(RP), intent(in)    :: gsqrtgam2 (ijdim,kdim,ldim)
    real(RP), intent(in)    :: phi       (ijdim,kdim,ldim)   ! geopotential
    real(RP), intent(in)    :: c2wfact   (ijdim,kdim,2,ldim)
    real(RP), intent(in)    :: c2wfact_Gz(ijdim,kdim,6,ldim)

    integer  :: kmin, kmax
    real(RP) :: CVdry

    integer  :: ij, l
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcopy(rho,vx,vy,vz,w,ein,tem,pre) &
    !$acc pcopy(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge) &
    !$acc pcopyin(gsqrtgam2,phi,c2wfact,c2wfact_Gz)

    kmin = ADM_kmin
    kmax = ADM_kmax

    CVdry = CONST_CVdry

    !--- Thermodynamical variables ( rho, ein, tem, pre, rhog, rhoge ), q = 0 at boundary

    do l = 1, ldim
       call BNDCND_thermo( ijdim,      & ! [IN]
                           tem(:,:,l), & ! [INOUT]
                           rho(:,:,l), & ! [INOUT]
                           pre(:,:,l), & ! [INOUT]
                           phi(:,:,l)  ) ! [IN]
    enddo

    !$acc kernels pcopy(rhog,ein,rhoge) pcopyin(rho,tem,gsqrtgam2)
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,rhog,rhoge,ein,rho,tem,gsqrtgam2,CVdry)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          rhog (ij,kmax+1,l) = rho(ij,kmax+1,l) * gsqrtgam2(ij,kmax+1,l)
          rhog (ij,kmin-1,l) = rho(ij,kmin-1,l) * gsqrtgam2(ij,kmin-1,l)

          ein  (ij,kmax+1,l) = CVdry * tem(ij,kmax+1,l)
          ein  (ij,kmin-1,l) = CVdry * tem(ij,kmin-1,l)

          rhoge(ij,kmax+1,l) = rhog(ij,kmax+1,l) * ein(ij,kmax+1,l)
          rhoge(ij,kmin-1,l) = rhog(ij,kmin-1,l) * ein(ij,kmin-1,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !--- Momentum ( rhogvx, rhogvy, rhogvz, vx, vy, vz )

    call BNDCND_rhovxvyvz( ijdim,         & ! [IN]
                           kdim,          & ! [IN]
                           ldim,          & ! [IN]
                           rhog  (:,:,:), & ! [IN]
                           rhogvx(:,:,:), & ! [INOUT]
                           rhogvy(:,:,:), & ! [INOUT]
                           rhogvz(:,:,:)  ) ! [INOUT]

    !$acc kernels pcopy(vx,vy,vz) pcopyin(rhogvx,rhogvy,rhogvz,rhog)
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,vx,vy,vz,rhog,rhogvx,rhogvy,rhogvz)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          vx(ij,kmax+1,l) = rhogvx(ij,kmax+1,l) / rhog(ij,kmax+1,l)
          vy(ij,kmax+1,l) = rhogvy(ij,kmax+1,l) / rhog(ij,kmax+1,l)
          vz(ij,kmax+1,l) = rhogvz(ij,kmax+1,l) / rhog(ij,kmax+1,l)

          vx(ij,kmin-1,l) = rhogvx(ij,kmin-1,l) / rhog(ij,kmin-1,l)
          vy(ij,kmin-1,l) = rhogvy(ij,kmin-1,l) / rhog(ij,kmin-1,l)
          vz(ij,kmin-1,l) = rhogvz(ij,kmin-1,l) / rhog(ij,kmin-1,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !--- Momentum ( rhogw, w )
    do l = 1, ldim
       call BNDCND_rhow( ijdim,               & ! [IN]
                         rhogvx    (:,:,l),   & ! [IN]
                         rhogvy    (:,:,l),   & ! [IN]
                         rhogvz    (:,:,l),   & ! [IN]
                         rhogw     (:,:,l),   & ! [INOUT]
                         c2wfact_Gz(:,:,:,l)  ) ! [IN]
    enddo

    !$acc kernels pcopy(w) pcopyin(rhog,rhogw,c2wfact)
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,w,rhogw,rhog,c2wfact)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          w(ij,kmax+1,l) = rhogw(ij,kmax+1,l) / ( c2wfact(ij,kmax+1,1,l) * rhog(ij,kmax+1,l) &
                                                + c2wfact(ij,kmax+1,2,l) * rhog(ij,kmax  ,l) )
          w(ij,kmin  ,l) = rhogw(ij,kmin  ,l) / ( c2wfact(ij,kmin  ,1,l) * rhog(ij,kmin  ,l) &
                                                + c2wfact(ij,kmin  ,2,l) * rhog(ij,kmin-1,l) )
          w(ij,kmin-1,l) = 0.0_RP
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_all

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for thermodynamical variables
  subroutine BNDCND_thermo_SP( &
       ijdim, &
       tem,   &
       rho,   &
       pre,   &
       phi    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    use mod_const, only: &
       CONST_GRAV, &
       CONST_Rdry
    implicit none

    integer,  intent(in)    :: ijdim           ! number of horizontal grid
    real(SP), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(SP), intent(inout) :: rho(ijdim,kdim) ! density
    real(SP), intent(inout) :: pre(ijdim,kdim) ! pressure
    real(SP), intent(in)    :: phi(ijdim,kdim) ! geopotential

    integer  :: kmin, kmax
    real(SP) :: GRAV, Rdry

    integer  :: ij

    real(SP) :: z,z1,z2,z3,p1,p2,p3
    real(SP) :: lag_intpl
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3))/((z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3))/((z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2))/((z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcopy(tem,rho,pre) &
    !$acc pcopyin(phi)

    kmin = ADM_kmin
    kmax = ADM_kmax

    GRAV = CONST_GRAV
    Rdry = CONST_Rdry

    !$acc kernels pcopy(tem,pre,rho) pcopyin(phi)
    !$omp parallel do default(none),private(ij,z,z1,z2,z3), &
    !$omp shared(ijdim,kmax,kmin,is_top_tem,is_top_epl,is_btm_tem,is_btm_epl, &
    !$omp        rho,pre,tem,phi,Rdry,GRAV)
    do ij = 1, ijdim

       if    ( is_top_tem ) then

          tem(ij,kmax+1) = tem(ij,kmax) ! dT/dz = 0

       elseif( is_top_epl ) then
          z  = phi(ij,kmax+1) / GRAV
          z1 = phi(ij,kmax  ) / GRAV
          z2 = phi(ij,kmax-1) / GRAV
          z3 = phi(ij,kmax-2) / GRAV

          tem(ij,kmax+1) = lag_intpl( z ,                 &
                                      z1, tem(ij,kmax  ), &
                                      z2, tem(ij,kmax-1), &
                                      z3, tem(ij,kmax-2)  )
       endif

       if   ( is_btm_tem ) then

          tem(ij,kmin-1) = tem(ij,kmin) ! dT/dz = 0

       elseif( is_btm_epl ) then
          z1 = phi(ij,kmin+2) / GRAV
          z2 = phi(ij,kmin+1) / GRAV
          z3 = phi(ij,kmin  ) / GRAV
          z  = phi(ij,kmin-1) / GRAV

          tem(ij,kmin-1) = lag_intpl( z,                  &
                                      z1, tem(ij,kmin+2), &
                                      z2, tem(ij,kmin+1), &
                                      z3, tem(ij,kmin  )  )
       endif

       !--- set the boundary of pressure ( hydrostatic balance )
       pre(ij,kmax+1) = pre(ij,kmax-1) - rho(ij,kmax) * ( phi(ij,kmax+1) - phi(ij,kmax-1) )
       pre(ij,kmin-1) = pre(ij,kmin+1) - rho(ij,kmin) * ( phi(ij,kmin-1) - phi(ij,kmin+1) )

       !--- set the boundary of density ( equation of state )
       rho(ij,kmax+1) = pre(ij,kmax+1) / ( Rdry * tem(ij,kmax+1) )
       rho(ij,kmin-1) = pre(ij,kmin-1) / ( Rdry * tem(ij,kmin-1) )

    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_thermo_SP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for thermodynamical variables
  subroutine BNDCND_thermo_DP( &
       ijdim, &
       tem,   &
       rho,   &
       pre,   &
       phi    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    use mod_const, only: &
       CONST_GRAV, &
       CONST_Rdry
    implicit none

    integer,  intent(in)    :: ijdim           ! number of horizontal grid
    real(DP), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(DP), intent(inout) :: rho(ijdim,kdim) ! density
    real(DP), intent(inout) :: pre(ijdim,kdim) ! pressure
    real(DP), intent(in)    :: phi(ijdim,kdim) ! geopotential

    integer  :: kmin, kmax
    real(DP) :: GRAV, Rdry

    integer  :: ij

    real(DP) :: z,z1,z2,z3,p1,p2,p3
    real(DP) :: lag_intpl
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3))/((z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3))/((z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2))/((z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcopy(tem,rho,pre) &
    !$acc pcopyin(phi)

    kmin = ADM_kmin
    kmax = ADM_kmax

    GRAV = CONST_GRAV
    Rdry = CONST_Rdry

    !$acc kernels pcopy(tem,pre,rho) pcopyin(phi)
    !$omp parallel do default(none),private(ij,z,z1,z2,z3), &
    !$omp shared(ijdim,kmax,kmin,is_top_tem,is_top_epl,is_btm_tem,is_btm_epl, &
    !$omp        rho,pre,tem,phi,Rdry,GRAV)
    do ij = 1, ijdim

       if    ( is_top_tem ) then

          tem(ij,kmax+1) = tem(ij,kmax) ! dT/dz = 0

       elseif( is_top_epl ) then
          z  = phi(ij,kmax+1) / GRAV
          z1 = phi(ij,kmax  ) / GRAV
          z2 = phi(ij,kmax-1) / GRAV
          z3 = phi(ij,kmax-2) / GRAV

          tem(ij,kmax+1) = lag_intpl( z ,                 &
                                      z1, tem(ij,kmax  ), &
                                      z2, tem(ij,kmax-1), &
                                      z3, tem(ij,kmax-2)  )
       endif

       if   ( is_btm_tem ) then

          tem(ij,kmin-1) = tem(ij,kmin) ! dT/dz = 0

       elseif( is_btm_epl ) then
          z1 = phi(ij,kmin+2) / GRAV
          z2 = phi(ij,kmin+1) / GRAV
          z3 = phi(ij,kmin  ) / GRAV
          z  = phi(ij,kmin-1) / GRAV

          tem(ij,kmin-1) = lag_intpl( z,                  &
                                      z1, tem(ij,kmin+2), &
                                      z2, tem(ij,kmin+1), &
                                      z3, tem(ij,kmin  )  )
       endif

       !--- set the boundary of pressure ( hydrostatic balance )
       pre(ij,kmax+1) = pre(ij,kmax-1) - rho(ij,kmax) * ( phi(ij,kmax+1) - phi(ij,kmax-1) )
       pre(ij,kmin-1) = pre(ij,kmin+1) - rho(ij,kmin) * ( phi(ij,kmin-1) - phi(ij,kmin+1) )

       !--- set the boundary of density ( equation of state )
       rho(ij,kmax+1) = pre(ij,kmax+1) / ( Rdry * tem(ij,kmax+1) )
       rho(ij,kmin-1) = pre(ij,kmin-1) / ( Rdry * tem(ij,kmin-1) )

    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_thermo_DP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for horizontal momentum
  subroutine BNDCND_rhovxvyvz( &
       ijdim,  &
       kdim,   &
       ldim,   &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz  )
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: ldim
    real(RP), intent(in)    :: rhog  (ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogvx(ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogvy(ijdim,kdim,ldim)
    real(RP), intent(inout) :: rhogvz(ijdim,kdim,ldim)

    integer  :: kmin, kmax

    integer  :: ij, l
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcopy(rhogvx,rhogvy,rhogvz) &
    !$acc pcopyin(rhog)

    kmin = ADM_kmin
    kmax = ADM_kmax

    !$acc kernels pcopy(rhogvx,rhogvy,rhogvz) pcopyin(rhog)
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,is_top_rigid,is_top_free,is_btm_rigid,is_btm_free, &
    !$omp        rhog,rhogvx,rhogvy,rhogvz)
    do l  = 1, ldim
       !$omp do
       do ij = 1, ijdim

          if    ( is_top_rigid ) then
             rhogvx(ij,kmax+1,l) = -rhogvx(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvy(ij,kmax+1,l) = -rhogvy(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvz(ij,kmax+1,l) = -rhogvz(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
          elseif( is_top_free  ) then
             rhogvx(ij,kmax+1,l) =  rhogvx(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvy(ij,kmax+1,l) =  rhogvy(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvz(ij,kmax+1,l) =  rhogvz(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
          endif

          if    ( is_btm_rigid ) then
             rhogvx(ij,kmin-1,l) = -rhogvx(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvy(ij,kmin-1,l) = -rhogvy(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvz(ij,kmin-1,l) = -rhogvz(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
          elseif( is_btm_free  ) then
             rhogvx(ij,kmin-1,l) =  rhogvx(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvy(ij,kmin-1,l) =  rhogvy(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvz(ij,kmin-1,l) =  rhogvz(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
          endif

       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_rhovxvyvz

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for vertical momentum
  subroutine BNDCND_rhow( &
       ijdim,  &
       rhogvx, &
       rhogvy, &
       rhogvz, &
       rhogw,  &
       c2wfact )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(in)    :: rhogvx (ijdim,kdim)
    real(RP), intent(in)    :: rhogvy (ijdim,kdim)
    real(RP), intent(in)    :: rhogvz (ijdim,kdim)
    real(RP), intent(inout) :: rhogw  (ijdim,kdim)
    real(RP), intent(in)    :: c2wfact(ijdim,kdim,6)

    integer  :: kmin, kmax

    integer  :: ij
    !---------------------------------------------------------------------------
    !$acc  data &
    !$acc& pcopy(rhogw) &
    !$acc& pcopyin(rhogvx,rhogvy,rhogvz,c2wfact)

    kmin = ADM_kmin
    kmax = ADM_kmax

    !$acc kernels pcopy(rhogw) pcopyin(rhogvx,rhogvy,rhogvz,c2wfact)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmax,kmin,is_top_rigid,is_top_free,is_btm_rigid,is_btm_free, &
    !$omp        rhogw,c2wfact,rhogvx,rhogvy,rhogvz)
    do ij = 1, ijdim

       if    ( is_top_rigid ) then
          rhogw(ij,kmax+1) = 0.0_RP
       elseif( is_top_free  ) then
          rhogw(ij,kmax+1) = - ( c2wfact(ij,kmax+1,1) * rhogvx(ij,kmax+1) &
                               + c2wfact(ij,kmax+1,2) * rhogvx(ij,kmax  ) &
                               + c2wfact(ij,kmax+1,3) * rhogvy(ij,kmax+1) &
                               + c2wfact(ij,kmax+1,4) * rhogvy(ij,kmax  ) &
                               + c2wfact(ij,kmax+1,5) * rhogvz(ij,kmax+1) &
                               + c2wfact(ij,kmax+1,6) * rhogvz(ij,kmax  ) )
       endif

       if    ( is_btm_rigid ) then
          rhogw(ij,kmin) = 0.0_RP
       elseif( is_btm_free  ) then
          rhogw(ij,kmin) = - ( c2wfact(ij,kmin,1) * rhogvx(ij,kmin  ) &
                             + c2wfact(ij,kmin,2) * rhogvx(ij,kmin-1) &
                             + c2wfact(ij,kmin,3) * rhogvy(ij,kmin  ) &
                             + c2wfact(ij,kmin,4) * rhogvy(ij,kmin-1) &
                             + c2wfact(ij,kmin,5) * rhogvz(ij,kmin  ) &
                             + c2wfact(ij,kmin,6) * rhogvz(ij,kmin-1) )
       endif

       rhogw(ij,kmin-1) = 0.0_RP

    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_rhow

end module mod_bndcnd
