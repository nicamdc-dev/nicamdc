!-------------------------------------------------------------------------------
!
!+  Artificial Forcing initialization module
!
!-------------------------------------------------------------------------------
module mod_forcing_init
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module holds the routines which are used to perform 
  !       model start-up operations for the individual domains.  
  !       This is the stage after inputting wrfinput and before
  !       calling 'integrate'.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      12-10-11  R.Yoshida: convert from physics_init for dry dyn-core experiments
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only :  &
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
  public :: forcing_init
  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------  
  ! [Mod] 12/10/12 R.Yoshida
!!$  subroutine artforce_init
  subroutine forcing_init( ctime )

    !
    !++ Used modules
    !
    use mod_adm, only :  &
         ADM_NSYS,       &
         ADM_CTL_FID,    &
         ADM_LOG_FID,    &
         ADM_kmin,       &
         ADM_kmax,       &
         ADM_kall
    use mod_runconf, only : &
         AF_TYPE
    use mod_af_driver, only : &
         af_init
    !
    implicit none
    !
    ! 09/04/14 [Add] T.Mitsui
    real(8), intent(in) :: ctime
    !
    !--- additional forcing initialization
    call af_init ( AF_TYPE )
    !
    return
    !
  end subroutine forcing_init
  !-----------------------------------------------------------------------------
end module mod_forcing_init
!-------------------------------------------------------------------------------
