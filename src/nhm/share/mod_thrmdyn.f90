!-------------------------------------------------------------------------------
!
!+  Thermodynamics variables module
!
!-------------------------------------------------------------------------------
module mod_thrmdyn
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines in which the thermodyanics variables
  !       are calculated.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                07-01-26   H.Tomita: Support all type of EIN_TYPE
  !                08-06-13   T.Mitsui: mod thrmdyn_ent as strict difinition
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only :       &
       kdim => ADM_kall,    &
       kmin => ADM_kmin,    &
       kmax => ADM_kmax
  use mod_cnst,  only :&
       CNST_RAIR,      &
       CNST_RVAP,      &
       CNST_CV,        &
       CNST_CP,        &
       CNST_PRE00,     &
       CNST_KAPPA,     &
       CNST_EPSV,      &
       CNST_TEM00,     &
       CNST_PSAT0,     &
       CNST_LHF0,      &
       CNST_LH0,       & 
       CNST_LHF00,     &
       CNST_LH00
  use mod_runconf, only : &
       nqmax => TRC_VMAX ,&
       NQW_STR,NQW_END,   &
       I_QV,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       I_QH,              &
       CVW,CPW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: thrmdyn_rho
  public :: thrmdyn_ein
  public :: thrmdyn_eth
  public :: thrmdyn_tem
  public :: thrmdyn_pre
  public :: thrmdyn_th
  public :: thrmdyn_tempreth
  public :: thrmdyn_tempre
  public :: thrmdyn_cv
  public :: thrmdyn_cp
  public :: thrmdyn_qd
  public :: thrmdyn_ent

  public :: THRMDYN_qd_ijkl
  public :: THRMDYN_rho_ijkl
  public :: THRMDYN_ein_ijkl
  public :: THRMDYN_tempre_ijkl
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_rho( &
       ijdim,             & !--- IN : number of horizontal grid
       rho,               & !--- OUT : density     
       pre,               & !--- IN  : pressure    
       tem,               & !--- IN  : temperature 
       qd,                & !--- IN  : dry concentration 
       q )                  !--- IN  : water concentration 
    !------ 
    !------ Density calculation in all regions
    !------    1. calculation region of rho
    !------                   : (:,:,:)
    !------ 
    implicit none
    !
    integer, intent(in) :: ijdim
    real(8), intent(out) :: rho(1:ijdim,1:kdim)
    real(8), intent(in)  :: pre(1:ijdim,1:kdim)
    real(8), intent(in)  :: tem(1:ijdim,1:kdim)
    real(8), intent(in)  :: qd(1:ijdim,1:kdim)
    real(8), intent(in)  :: q(1:ijdim,1:kdim,1:nqmax)
    !
    rho(:,:) = pre(:,:) / tem(:,:) &
         / ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )
    !
  end subroutine thrmdyn_rho
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_ein( &
       ijdim,             & !--- IN : number of horizontal grid
       ein,               &  !--- OUT : internal energy
       tem,               &  !--- IN  : temperature
       qd,                &  !--- IN  : dry concentration 
       q )                   !--- IN  : water concentration 
    !------ 
    !------ Internal energy calculation in all regions
    !------    1. calculation region of ein
    !------                   : (:,:)
    !------ 
    !------ CAUTION : ein = CV*T*qd + CVV*T*qv + CPL*T*qc
    !------           ein_moist = CV*T*qd + (CVV*T+LH00)*qv + CPL*T*qc
    !              
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: ein(1:ijdim,1:kdim)
    real(8), intent(in)  :: tem(1:ijdim,1:kdim)
    real(8), intent(in)  :: qd(1:ijdim,1:kdim)
    real(8), intent(in)  :: q(1:ijdim,1:kdim,1:nqmax)
    !
    real(8) :: cv(1:ijdim,1:kdim)
    !
    call thrmdyn_cv( &
         ijdim,      & !--- in
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    ein(:,:) = cv(:,:) * tem(:,:)
    !
  end subroutine thrmdyn_ein
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_eth( &
       ijdim,             &  !--- IN : number of horizontal grid
       eth,               &  !--- OUT : enthalpy
       ein,               &  !--- IN  : internal energy
       pre,               &  !--- IN  : pressure
       rho )                 !--- IN  : density
    !------ 
    !------ Enthalpy calculation in all regions
    !------    1. calculation region of eth
    !------                   : (:,:)
    !------ 
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: eth(1:ijdim,1:kdim)
    real(8), intent(in)  :: ein(1:ijdim,1:kdim)
    real(8), intent(in)  :: pre(1:ijdim,1:kdim)
    real(8), intent(in)  :: rho(1:ijdim,1:kdim)
    !
    eth(:,:) = ein(:,:) + pre(:,:) / rho(:,:)
    !
  end subroutine thrmdyn_eth
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tem( &
       ijdim,             &  !--- IN : number of horizontal grid
       tem,               &  !--- OUT  : temperature       
       ein,               &  !--- IN : internal energy
       qd,                &  !--- IN  : dry concentration 
       q )                   !--- IN  : water concentration 
    !------ 
    !------ temperature calculation in all regions
    !------    1. calculation region of tem
    !------                   : (:,:)
    !------ 
    !------ CAUTION : ein = CV*T*qd + CVV*T*qv + CPL*T*qc
    !------           ein_moist = CV*T*qd + (CVV*T+LH00)*qv + CPL*T*qc
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out)  :: tem(1:ijdim,1:kdim)
    real(8), intent(in)   :: ein(1:ijdim,1:kdim)
    real(8), intent(in)  :: qd(1:ijdim,1:kdim)
    real(8), intent(in)  :: q(1:ijdim,1:kdim,1:nqmax)
    !
    real(8)  :: cv(1:ijdim,1:kdim)
    !
    call thrmdyn_cv( &
         ijdim,      & !--- in
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    tem(:,:) = ein(:,:) / cv(:,:)
    !
  end subroutine thrmdyn_tem
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_pre( &
       ijdim,             &  !--- IN : number of horizontal grid
       pre,               &  !--- OUT : pressure
       tem,               &  !--- IN  : temperature       
       rho,               &  !--- IN  : density
       qd,                &  !--- IN  : dry concentration 
       q )                   !--- IN  : water concentration 

    !------ 
    !------ Pressure calculation in all regions
    !------    1. calculation region of pre
    !------                   : (:,:)
    !------ 
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: pre(1:ijdim,1:kdim)
    real(8), intent(in)  :: tem(1:ijdim,1:kdim)
    real(8), intent(in)  :: rho(1:ijdim,1:kdim)
    real(8), intent(in)  :: qd(1:ijdim,1:kdim)
    real(8), intent(in)  :: q(1:ijdim,1:kdim,1:nqmax)
    !
    pre(:,:) = rho(:,:) * tem(:,:) &
         * ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )
    !
  end subroutine thrmdyn_pre
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_th( &
       ijdim,            &  !--- IN : number of horizontal grid
       th,               &  !--- OUT : potential temperature
       tem,              &  !--- IN  : temperature       
       pre )                !--- IN  : pressure
    !------ 
    !------ Potential temperature calculation in all regions
    !------    1. calculation region of th
    !------                   : (:,:)
    !------ 
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: th(1:ijdim,1:kdim)
    real(8), intent(in)  :: tem(1:ijdim,1:kdim)
    real(8), intent(in)  :: pre(1:ijdim,1:kdim)
    !
    real(8) :: p0k
    !
    p0k=CNST_PRE00**CNST_KAPPA
    !
    th(:,:) =tem(:,:)*(abs(pre(:,:))**(-CNST_KAPPA))*p0k
    !
  end subroutine thrmdyn_th
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempreth( &
       ijdim,                  &  !--- IN : number of horizontal grid
       tem,                    &  !--- OUT : temperature              
       pre,                    &  !--- OUT : pressure
       th,                     &  !--- OUT : potential temperature
       ein,                    &  !--- IN  : internal energy
       rho,                    &  !--- IN  : density
       qd,                     &  !--- IN  : dry concentration 
       q )                        !--- IN  : water concentration 
    !------ 
    !------  Calculation of temperature, pressure, 
    !------  potential temperature in all regions.
    !------    1. calculation region  : (:,:)
    !------ 
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: tem(1:ijdim,1:kdim)
    real(8), intent(out) :: pre(1:ijdim,1:kdim)
    real(8), intent(out) :: th(1:ijdim,1:kdim)
    real(8), intent(in)  :: ein(1:ijdim,1:kdim)
    real(8), intent(in)  :: rho(1:ijdim,1:kdim)
    real(8), intent(in)  :: qd(1:ijdim,1:kdim)
    real(8), intent(in)  :: q(1:ijdim,1:kdim,1:nqmax)
    !
    real(8) :: cv(1:ijdim,1:kdim)
    !
    real(8) :: p0k
    !
    p0k=CNST_PRE00**CNST_KAPPA
    !
    call thrmdyn_cv( &
         ijdim,      & !--- in
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    !
    tem(:,:) = ein(:,:) / cv(:,:)
    !
    pre(:,:) = rho(:,:) * tem(:,:) &
         * ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )
    th(:,:) =tem(:,:)*(abs(pre(:,:))**(-CNST_KAPPA))*p0k
    !
  end subroutine thrmdyn_tempreth
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempre( &
       ijdim,                &  !--- IN : number of horizontal grid
       tem,                  &  !--- OUT : temperature              
       pre,                  &  !--- OUT : pressure
       ein,                  &  !--- IN  : internal energy
       rho,                  &  !--- IN  : density
       qd,                   &  !--- IN  : dry concentration 
       q )                      !--- IN  : water concentration 
    !------ 
    !------  Calculation of temperature, pressure, 
    !------  potential temperature in all regions.
    !------    1. calculation region  : (:,:)
    !------ 
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: tem(ijdim,kdim)
    real(8), intent(out) :: pre(ijdim,kdim)
    real(8), intent(in)  :: ein(ijdim,kdim)
    real(8), intent(in)  :: rho(ijdim,kdim)
    real(8), intent(in)  :: qd (ijdim,kdim)
    real(8), intent(in)  :: q  (ijdim,kdim,nqmax)

    real(8) :: cv(1:ijdim,1:kdim)

    call thrmdyn_cv( &
         ijdim,      & !--- in
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in

    tem(:,:) = ein(:,:) / cv(:,:)

    pre(:,:) = rho(:,:) * tem(:,:) &
         * ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )

    return
  end subroutine thrmdyn_tempre

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cv( &
       ijdim,            & !--- IN : number of horizontal grid
       cva,              & !--- OUT : specific heat
       q,                & !--- IN  : mass concentration
       qd)                 !--- IN  : dry mass concentration
    !
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: cva(1:ijdim,1:kdim)
    real(8), intent(in) :: q(1:ijdim,1:kdim,1:nqmax)
    real(8), intent(in) :: qd(1:ijdim,1:kdim)
    !
    integer :: nq
    !
    cva(:,:) = qd(:,:) * CNST_CV
    do nq = NQW_STR,NQW_END
       cva(:,:) = cva(:,:) + q(:,:,nq) * CVW(nq)
    end do
    !
  end subroutine thrmdyn_cv
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cp( &
       ijdim,            & !--- IN : number of horizontal grid
       kdim_local,       & !--- IN :
       cpa,              & !--- OUT : specific heat
       q,                & !--- IN  : mass concentration
       qd )                !--- IN  : dry mass concentration
    !
    implicit none
    !
    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim_local
    real(8), intent(out) :: cpa(1:ijdim,1:kdim_local)
    real(8), intent(in) :: q(1:ijdim,1:kdim_local,1:nqmax)
    real(8), intent(in) :: qd(1:ijdim,1:kdim_local)
    !
    integer :: nq
    !
    cpa(:,:) = qd(:,:) * CNST_CP
    do nq = NQW_STR,NQW_END
       cpa(:,:) = cpa(:,:) + q(:,:,nq) * CPW(nq)
    end do
    !
  end subroutine thrmdyn_cp
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_qd( &
       ijdim,            & !--- IN : number of horizontal grid
       qd,               & !--- OUT  : dry mass concentration
       q )                 !--- IN  : mass concentration
    !
    implicit none
    !
    integer, intent(in)  :: ijdim
    real(8), intent(out) :: qd(1:ijdim,1:kdim)
    real(8), intent(in) :: q(1:ijdim,1:kdim,1:nqmax)
    !
    integer :: nq
    !
    qd(:,:) = 1.0D0
    do nq = NQW_STR,NQW_END
       qd(:,:) = qd(:,:) - q(:,:,nq)
    end do
    !
  end subroutine thrmdyn_qd
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_ent( &
       ijdim,             & !--- IN : 
       kdim_local,        & !--- IN :
       ent,               & !--- OUT :
       tem,               & !--- IN
       pre,               & !--- IN
       q,                 & !--- IN
       qd                 & !--- IN
       )
    use mod_cnst, only : &
         CNST_RAIR,      &
         CNST_RVAP,      &
         CNST_CI,        &
         CNST_CL,        &
         CNST_CP,        &
         CNST_CPV,       &
         CNST_PRE00,     &
         CNST_EPSV,      &
         CNST_TEM00,     &
         CNST_PSAT0,     &
         CNST_LH0,       &
         CNST_LHF0
    use mod_runconf, only:  &
         nqmax => TRC_VMAX, &
         I_QC, I_QR, I_QI, I_QS, I_QG, I_QH
    !    
    implicit none
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim_local
    real(8), intent(out) :: ent(1:ijdim,1:kdim_local)
    real(8), intent(in) :: tem(1:ijdim,1:kdim_local)
    real(8), intent(in) :: pre(1:ijdim,1:kdim_local)
    real(8), intent(in) :: q(1:ijdim,1:kdim_local,1:nqmax)
    real(8), intent(in) :: qd(1:ijdim,1:kdim_local)
    !
    real(8) :: pred(1:ijdim,1:kdim_local)
    real(8) :: prev(1:ijdim,1:kdim_local)

    real(8), parameter :: PREMIN = 1.d-10
    !
    integer :: nq
    !
    ! entropy
!!$    ! ent = qd * sd + qv * sv + qc * sc
    ! ent = qd * sd + qv * sv + qc * sc + qr * sr ...
    pred(:,:) = pre(:,:) * CNST_EPSV * qd(:,:) &
         / ( CNST_EPSV * qd(:,:) + q(:,:,I_QV) )
    prev(:,:) = pre(:,:) * q(:,:,I_QV) &
         / ( CNST_EPSV * qd(:,:) + q(:,:,I_QV) )
    pred(:,:) = max( pred(:,:), PREMIN )
    prev(:,:) = max( prev(:,:), PREMIN )
    
    !     call thrmdyn_cp( &
    !          ijdim,      & !--- IN
    !          kdim_local, & !--- IN
    !          cpa,        & !--- OUT
    !          q,          & !--- IN
    !          qd )          !--- IN
    ! 08/06/13 [fix] T.Mitsui, strict definition of entropy is given by 
    ! (8),(11),(13) Satoh (2003), Mon.Wea.Rev.
!!$     cpa(:,:) = CNST_CP
!!$     !
!!$     ent(:,:) = cpa(:,:) * log( tem(:,:) / CNST_TEM00 ) &
!!$          - qd(:,:) * CNST_RAIR * log ( pred(:,:) / CNST_PRE00 ) &
!!$          + q(:,:,I_QV) * ( - CNST_RVAP * log ( ( prev(:,:) / CNST_PSAT0 ) ) &
!!$          + CNST_LH0 / CNST_TEM00 )
    !
    ! write(*,*) cpa(1,1), ent(1,1)
    !
!!$     do nq = NQW_STR, NQW_END
!!$        if ( nq == I_QI .or. nq == I_QS &
!!$             .or. nq == I_QG .or. nq == I_QH ) then
!!$           ent(:,:) = ent(:,:) - q(:,:,nq) * CNST_LHF0 / CNST_TEM00
!!$        end if
!!$     end do
    !
    ent(:,:) = qd(:,:) * CNST_CP   * log ( tem(:,:)  / CNST_TEM00 )   &
         -     qd(:,:) * CNST_RAIR * log ( pred(:,:) / CNST_PRE00 )   &
         + q(:,:,I_QV) * CNST_CPV  * log ( tem(:,:)  / CNST_TEM00 )   &
         - q(:,:,I_QV) * CNST_RVAP * log ( prev(:,:) / CNST_PSAT0 )   &
         + q(:,:,I_QV) * CNST_LH0 / CNST_TEM00 
    do nq = NQW_STR, NQW_END
       if      ( (nq==I_QC) .or. (nq==I_QR) )then
          ent(:,:) = ent(:,:) + q(:,:,nq) * CNST_CL * log( tem(:,:) / CNST_TEM00 )
       else if ( (nq==I_QI) .or. (nq==I_QS) .or. (nq==I_QG) .or. (nq==I_QH) ) then
          ent(:,:) = ent(:,:) + q(:,:,nq) * CNST_CI * log( tem(:,:) / CNST_TEM00 ) &
               - q(:,:,nq) * CNST_LHF0 / CNST_TEM00
       end if
    end do
    !
    return
  end subroutine thrmdyn_ent

  !-----------------------------------------------------------------------------
  subroutine THRMDYN_qd_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       vmax,  &
       wstr,  &
       wend,  &
       qd,    &
       q      )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    integer, intent(in)  :: vmax                     ! total number of tracers
    integer, intent(in)  :: wstr                     ! start index of water mass tracers
    integer, intent(in)  :: wend                     ! end   index of water mass tracers
    real(8), intent(out) :: qd(ijdim,kdim,ldim)      ! dry air mass concentration [kg/kg]
    real(8), intent(in)  :: q (ijdim,kdim,ldim,vmax) ! tracer  mass concentration [kg/kg]

    integer :: iw
    !---------------------------------------------------------------------------

    qd(:,:,:) = 1.D0

    do iw = wstr, wend
       qd(:,:,:) = qd(:,:,:) - q(:,:,:,iw)
    enddo

    return
  end subroutine THRMDYN_qd_ijkl

  !-----------------------------------------------------------------------------
  subroutine THRMDYN_rho_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       rho,   &
       pre,   &
       tem,   &
       qd,    &
       qv     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(8), intent(out) :: rho(ijdim,kdim,ldim) ! density     [kg/m3]
    real(8), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure    [Pa]
    real(8), intent(in)  :: tem(ijdim,kdim,ldim) ! temperature [K]
    real(8), intent(in)  :: qd (ijdim,kdim,ldim) ! dry air     mass concentration [kg/kg]
    real(8), intent(in)  :: qv (ijdim,kdim,ldim) ! water vapor mass concentration [kg/kg]
    !---------------------------------------------------------------------------

    rho(:,:,:) = pre(:,:,:) / ( ( qd(:,:,:)*CNST_RAIR + qv(:,:,:)*CNST_RVAP ) * tem(:,:,:) )

    return
  end subroutine THRMDYN_rho_ijkl

  !-----------------------------------------------------------------------------
  subroutine THRMDYN_ein_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       vmax,  &
       wstr,  &
       wend,  &
       ein,   &
       tem,   &
       qd,    &
       q      )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    integer, intent(in)  :: vmax                      ! total number of tracers
    integer, intent(in)  :: wstr                      ! start index of water mass tracers
    integer, intent(in)  :: wend                      ! end   index of water mass tracers
    real(8), intent(out) :: ein(ijdim,kdim,ldim)      ! internal energy [J]
    real(8), intent(in)  :: tem(ijdim,kdim,ldim)      ! temperature     [K]
    real(8), intent(in)  :: qd (ijdim,kdim,ldim)      ! dry air mass concentration [kg/kg]
    real(8), intent(in)  :: q  (ijdim,kdim,ldim,vmax) ! tracer  mass concentration [kg/kg]

    real(8) :: cv(ijdim,kdim,ldim)

    integer :: iw
    !---------------------------------------------------------------------------

    cv(:,:,:) = qd(:,:,:) * CNST_CV

    do iw = wstr, wend
       cv(:,:,:) = cv(:,:,:) + q(:,:,:,iw) * CVW(iw)
    enddo

    ein(:,:,:) = tem(:,:,:) * cv(:,:,:)

    return
  end subroutine THRMDYN_ein_ijkl

  !-----------------------------------------------------------------------------
  subroutine THRMDYN_tempre_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       vmax,  &
       wstr,  &
       wend,  &
       tem,   &
       pre,   &
       ein,   &
       rho,   &
       qd,    &
       q      )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    integer, intent(in)  :: vmax                      ! total number of tracers
    integer, intent(in)  :: wstr                      ! start index of water mass tracers
    integer, intent(in)  :: wend                      ! end   index of water mass tracers
    real(8), intent(out) :: tem(ijdim,kdim,ldim)      ! temperature     [K]
    real(8), intent(out) :: pre(ijdim,kdim,ldim)      ! pressure    [Pa]
    real(8), intent(in)  :: ein(ijdim,kdim,ldim)      ! internal energy [J]
    real(8), intent(in)  :: rho(ijdim,kdim,ldim)      ! density     [kg/m3]
    real(8), intent(in)  :: qd (ijdim,kdim,ldim)      ! dry air mass concentration [kg/kg]
    real(8), intent(in)  :: q  (ijdim,kdim,ldim,vmax) ! tracer  mass concentration [kg/kg]

    real(8) :: cv(ijdim,kdim,ldim)

    integer :: iw
    !---------------------------------------------------------------------------

    cv(:,:,:) = qd(:,:,:) * CNST_CV

    do iw = wstr, wend
       cv(:,:,:) = cv(:,:,:) + q(:,:,:,iw) * CVW(iw)
    enddo

    tem(:,:,:) = ein(:,:,:) / cv(:,:,:)

    pre(:,:,:) = rho(:,:,:) * tem(:,:,:) * ( qd(:,:,:)*CNST_RAIR + q(:,:,:,wstr)*CNST_RVAP )

    return
  end subroutine THRMDYN_tempre_ijkl

end module mod_thrmdyn
!-------------------------------------------------------------------------------
