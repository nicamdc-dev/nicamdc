!-------------------------------------------------------------------------------
!
!+  experimental forcing driver module
!
!-------------------------------------------------------------------------------
module mod_af_driver
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module sets up and load the additional forcing process.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Add this module.
  !                06-04-19   Support Held-Suarez forcing.
  !                10-10-16   A.Noda : Support TWP-ICE exp.
  !                11-02-18   A.Noda : update for TWP-ICE exp.
  !                12-03-09   S.Iga: tuned (phase4-1)
  !                12-10-11   R.Yoshida : gcss and twpice are removed for DC test.
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  character(ADM_NSYS), public, save :: AF_TYPE = 'NONE' 
  !                                             = 'GCSS_CASE1'
  !                                             = 'HELD-SUAREZ'
  !                                             = 'TWP-ICE'

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: af_init
  public :: af_driver

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_init( AF_TYPE_in )
    use mod_af_heldsuarez, only: &
       af_heldsuarez_init
    implicit none

    character(len=*), intent(in) :: AF_TYPE_in
    !---------------------------------------------------------------------------

    AF_TYPE = trim(AF_TYPE_in)

    select case(AF_TYPE)
    case('NONE')
       !--- do nothing
    case('HELD-SUAREZ')
       call af_heldsuarez_init
    case default
       write(*,*) 'Msg : Sub[af_init]/Mod[af_driver]'
       write(*,*) ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
    end select

    return
  end subroutine af_init

  !-----------------------------------------------------------------------------
  subroutine af_driver( &
       ijdim,     &
       rho,       &
       pre,       &
       tem,       &
       vx,        &
       vy,        &
       vz,        &
       w,         &
       lat,       &
       z,         &
       frhog,     &
       frhogvx,   &
       frhogvy,   &
       frhogvz,   &
       frhogw,    &
       frhoge,    &
       frhogetot, &
       frhogq,    &
       gsgam2,    &
       gsgam2h    )
    use mod_adm, only: &
       kdim => ADM_kall
    use mod_cnst, only: &
       CNST_EGRAV
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_af_heldsuarez, only: &
       af_heldsuarez
    implicit none

    integer, intent(in)  :: ijdim
    real(8), intent(in)  :: rho      (ijdim,kdim)
    real(8), intent(in)  :: pre      (ijdim,kdim)
    real(8), intent(in)  :: tem      (ijdim,kdim)
    real(8), intent(in)  :: vx       (ijdim,kdim)
    real(8), intent(in)  :: vy       (ijdim,kdim)
    real(8), intent(in)  :: vz       (ijdim,kdim)
    real(8), intent(in)  :: w        (ijdim,kdim)
    real(8), intent(in)  :: lat      (ijdim)
    real(8), intent(in)  :: z        (ijdim,kdim)
    real(8), intent(out) :: frhog    (ijdim,kdim)
    real(8), intent(out) :: frhogvx  (ijdim,kdim)
    real(8), intent(out) :: frhogvy  (ijdim,kdim)
    real(8), intent(out) :: frhogvz  (ijdim,kdim)
    real(8), intent(out) :: frhogw   (ijdim,kdim)
    real(8), intent(out) :: frhoge   (ijdim,kdim)
    real(8), intent(out) :: frhogetot(ijdim,kdim)
    real(8), intent(out) :: frhogq   (ijdim,kdim,nqmax)
    real(8), intent(in)  :: gsgam2   (ijdim,kdim)
    real(8), intent(in)  :: gsgam2h  (ijdim,kdim)

    real(8) :: fvx(ijdim,kdim)
    real(8) :: fvy(ijdim,kdim)
    real(8) :: fvz(ijdim,kdim)
    real(8) :: fw (ijdim,kdim)
    real(8) :: fe (ijdim,kdim)
    real(8) :: fq (ijdim,kdim,nqmax)
    !
    integer :: nq
    !---------------------------------------------------------------------------

    select case(trim(AF_TYPE))
    case('HELD-SUAREZ')

       call af_HeldSuarez( ijdim, &
                           lat,   & ! [IN] : latitude  
                           rho,   & ! [IN] : density
                           vx,    & ! [IN] : Vx
                           vy,    & ! [IN] : Vy
                           vz,    & ! [IN] : Vz
                           w,     & ! [IN] : w
                           pre,   & ! [IN] : pressure
                           tem,   & ! [IN] : temperature
                           fvx,   & ! [OUT]: tend. of rhovx
                           fvy,   & ! [OUT]: tend. of rhovy
                           fvz,   & ! [OUT]: tend. of rhovz
                           fw,    & ! [OUT]: tend. of rhow
                           fe     ) ! [OUT]: tend. of rhoe

       fq  = 0.D0

    case default

       fvx = 0.D0
       fvy = 0.D0
       fvz = 0.D0
       fw  = 0.D0
       fe  = 0.D0
       fq  = 0.D0

    end select

    frhogvx  (:,:) = fvx(:,:) * rho(:,:) * gsgam2 (:,:)
    frhogvy  (:,:) = fvy(:,:) * rho(:,:) * gsgam2 (:,:)
    frhogvz  (:,:) = fvz(:,:) * rho(:,:) * gsgam2 (:,:)
    frhogw   (:,:) = fw (:,:) * rho(:,:) * gsgam2h(:,:)
    frhoge   (:,:) = fe (:,:) * rho(:,:) * gsgam2 (:,:)
    frhogetot(:,:) = fe (:,:) * rho(:,:) * gsgam2 (:,:)

    frhog(:,:) = 0.D0
    do nq = 1, nqmax
       frhogq   (:,:,nq) = fq(:,:,nq) * rho(:,:) * gsgam2(:,:)

       frhog    (:,:) = frhog    (:,:) + frhogq(:,:,nq)
       frhogetot(:,:) = frhogetot(:,:) + frhogq(:,:,nq) * CNST_EGRAV * z(:,:)
    enddo

    return
  end subroutine af_driver

end module mod_af_driver
!-------------------------------------------------------------------------------
