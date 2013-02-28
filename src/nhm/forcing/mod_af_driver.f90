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
  !
  !++ Used modules
  !
  use mod_adm, only :       &
       ADM_NSYS
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
  public :: af_init
  public :: af_driver
  !-----------------------------------------------------------------------------
  !
  !
  character(ADM_NSYS), public, save :: AF_TYPE = 'NONE' 
  !                                             = 'GCSS_CASE1'
  !                                             = 'HELD-SUAREZ'
  !                                             = 'TWP-ICE'
  !++ Private variables
!  character(ADM_NSYS), private, save :: AF_TYPE = 'NONE' 
  !
contains
  !-----------------------------------------------------------------------------
  subroutine af_init ( AF_TYPE_in )

    use mod_af_heldsuarez, only : &
         af_heldsuarez_init

    implicit none
    character(len=*), intent(in) :: AF_TYPE_in
    !
    AF_TYPE = trim(AF_TYPE_in)
    !
    select case(trim(AF_TYPE))
    case('NONE')
       !--- nothing
    case('HELD-SUAREZ')
       call af_heldsuarez_init
    case default
       write(*,*) &
            'Msg : Sub[af_init]/Mod[af_driver]'
       write(*,*) &
            ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
    end select
  end subroutine af_init
  !-----------------------------------------------------------------------------
  !------------- a part of porting vars is deleted: 12/10/12 R.Yoshida
  subroutine af_driver(&
       ijdim,            &  !--- IN : number of horizontal grid
       rho,              &  !--- IN : density
       pre,              &  !--- IN : pressure
       tem,              &  !--- IN : temperature
       vx,               &  !--- IN : Vx
       vy,               &  !--- IN : Vy
       vz,               &  !--- IN : Vz
       w,                &  !--- IN : Vz
       lat,              &  !--- IN : latitude
       z,                &  !--- IN : height
       frhog,            &  !--- INOUT : forcing for rhog   ( gam2 X G^{1/2} )
       frhogvx,          &  !--- INOUT : forcing for rhogvx ( gam2 X G^{1/2} )
       frhogvy,          &  !--- INOUT : forcing for rhogvy ( gam2 X G^{1/2} )
       frhogvz,          &  !--- INOUT : forcing for rhogvz ( gam2 X G^{1/2} )
       frhogw,           &  !--- INOUT : forcing for rhogw  ( gam2 X G^{1/2} )
       frhoge,           &  !--- INOUT : forcing for rhoge  ( gam2 X G^{1/2} )
       frhogetot,        &  !--- INOUT : forcing for rhogetot ( gam2 X G^{1/2} )
       frhogq,           &  !--- INOUT : forcing for rhogq  ( gam2 X G^{1/2} )
       gsgam2,           &  !--- IN : gsgam2
       gsgam2h,          &  !--- IN : gsgam2
       add_type          &  !--- IN : append or renew
       )
    !
    use mod_adm, only :       &
       ADM_LOG_FID,         &
       kdim => ADM_kall

    use mod_runconf, only :     &
       nqmax => TRC_VMAX
    use mod_cnst, only : &
       CNST_EGRAV,     &
       CNST_PRE00,     &
       CNST_KAPPA
    use mod_af_heldsuarez, only : &
         af_heldsuarez
    !
    implicit none
    !
    integer, intent(in) :: ijdim
    real(8), intent(in) :: rho(1:ijdim,1:kdim)
    real(8), intent(in) :: pre(1:ijdim,1:kdim)
    real(8), intent(in) :: tem(1:ijdim,1:kdim)
    real(8), intent(in) :: vx(1:ijdim,1:kdim)
    real(8), intent(in) :: vy(1:ijdim,1:kdim)
    real(8), intent(in) :: vz(1:ijdim,1:kdim)
    real(8), intent(in) :: w(1:ijdim,1:kdim)
    real(8), intent(in) :: lat(1:ijdim)
    real(8), intent(in) :: z(1:ijdim,1:kdim)
    !
    real(8), intent(inout) :: frhog(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhogvx(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhogvy(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhogvz(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhogw(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhoge(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhogetot(1:ijdim,1:kdim)
    real(8), intent(inout) :: frhogq(1:ijdim,1:kdim,1:nqmax)
    !
    real(8), intent(in) :: gsgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsgam2h(1:ijdim,1:kdim)
    !
    character(*), intent(in) :: add_type
    !
    real(8) :: fvx(1:ijdim,1:kdim)
    real(8) :: fvy(1:ijdim,1:kdim)
    real(8) :: fvz(1:ijdim,1:kdim)
    real(8) :: fw(1:ijdim,1:kdim)
    real(8) :: fe(1:ijdim,1:kdim)
    real(8) :: fq(1:ijdim,1:kdim,1:nqmax)
    !
    integer :: nq
    !
    select case(trim(AF_TYPE))
    case('HELD-SUAREZ')
       call af_HeldSuarez(     &
            ijdim,             &
            lat,               &  !--- IN : latitude  
            rho,               &  !--- IN : density
            vx,                &  !--- IN : Vx
            vy,                &  !--- IN : Vy
            vz,                &  !--- IN : Vz
            w,                 &  !--- IN : w
            pre,               &  !--- IN : pressure
            tem,               &  !--- IN : temperature
            fvx,               &  !--- INOUT : tend. of rhovx
            fvy,               &  !--- INOUT : tend. of rhovy
            fvz,               &  !--- INOUT : tend. of rhovz
            fw,                &  !--- INOUT : tend. of rhow
            fe                 &  !--- INOUT : tend. of rhoe
            )
       fq=0.0D0
    case default
       !write(ADM_LOG_FID,*) &
       !     'Msg : Sub[af_driver]/Mod[af_driver]'
       !write(ADM_LOG_FID,*) &
       !     ' **** NO FORCING'
       fvx= 0.0D0
       fvy= 0.0D0
       fvz= 0.0D0
       fw = 0.0D0
       fe = 0.0D0
       fq = 0.0D0
    end select
    !
    if(trim(add_type)=='APPEND') then
       frhogvx(:,:) = frhogvx(:,:) + fvx(:,:)*rho(:,:)*gsgam2(:,:)
       frhogvy(:,:) = frhogvy(:,:) + fvy(:,:)*rho(:,:)*gsgam2(:,:)
       frhogvz(:,:) = frhogvz(:,:) + fvz(:,:)*rho(:,:)*gsgam2(:,:)
       frhogw (:,:) = frhogw (:,:) + fw(:,:)*rho(:,:)*gsgam2h(:,:)
       frhoge (:,:) = frhoge (:,:) + fe(:,:)*rho(:,:)*gsgam2(:,:)
       frhogetot (:,:) = frhogetot (:,:) + fe(:,:)*rho(:,:)*gsgam2(:,:)
       do nq=1,nqmax
          frhog(:,:) = frhog(:,:) + fq(:,:,nq)*rho(:,:)*gsgam2(:,:)
          frhogq (:,:,nq) = frhogq (:,:,nq) + fq(:,:,nq)*rho(:,:)*gsgam2(:,:)
          !--- add potential energy
          frhogetot (:,:) = frhogetot (:,:) + fq(:,:,nq)*rho(:,:)*gsgam2(:,:)*CNST_EGRAV*z(:,:)
       end do
    elseif(trim(add_type)=='RENEW') then
       frhogvx(:,:) = fvx(:,:)*rho(:,:)*gsgam2(:,:)
       frhogvy(:,:) = fvy(:,:)*rho(:,:)*gsgam2(:,:)
       frhogvz(:,:) = fvz(:,:)*rho(:,:)*gsgam2(:,:)
       frhogw (:,:) = fw(:,:)*rho(:,:)*gsgam2h(:,:)
       frhoge (:,:) = fe(:,:)*rho(:,:)*gsgam2(:,:)
       frhogetot (:,:) = fe(:,:)*rho(:,:)*gsgam2(:,:)
       frhog(:,:) = 0.0D0
       do nq=1,nqmax
          frhog(:,:) =  + fq(:,:,nq)*rho(:,:)*gsgam2(:,:)
          frhogq (:,:,nq) = fq(:,:,nq)*rho(:,:)*gsgam2(:,:)
          !--- add potential energy
          frhogetot (:,:) = frhogetot (:,:) + fq(:,:,nq)*rho(:,:)*gsgam2(:,:)*CNST_EGRAV*z(:,:)
       end do
    else
       write(*,*) 'ERROR!', add_type
    end if

  end subroutine af_driver
  !-----------------------------------------------------------------------------
end module mod_af_driver
!-------------------------------------------------------------------------------
