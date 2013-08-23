!-------------------------------------------------------------------------------
!>
!! Tracer variable module
!!
!! @par Description
!!         This module contains the chemical or general-perpose tracer variables
!!
!! @author H.Yashiro
!!
!! @par History
!! @li      2012-11-07 (H.Yashiro)  [NEW]
!!
!<
module mod_chemvar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: chemvar_getid

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                 public, parameter :: CHEM_TRC_vmax = 4
  character(len=16),       public, save      :: CHEM_TRC_name(CHEM_TRC_vmax) ! short name  of tracer
  character(len=ADM_NSYS), public, save      :: CHEM_TRC_desc(CHEM_TRC_vmax) ! description of tracer

  data CHEM_TRC_name / 'passive1', &
                       'passive2', &
                       'passive3', &
                       'passive4'  /

  data CHEM_TRC_desc / 'passive_tracer_no1', &
                       'passive_tracer_no2', &
                       'passive_tracer_no3', &
                       'passive_tracer_no4'  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  function chemvar_getid( tracername )
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    character(len=*), intent(in) :: tracername
    integer                      :: chemvar_getid

    character(len=16) :: tname
    integer           :: itrc
    !---------------------------------------------------------------------------

    tname = trim(tracername)

    chemvar_getid = -1

    do itrc = 1, CHEM_TRC_vmax
       if ( tname == CHEM_TRC_name(itrc) ) then
          chemvar_getid = itrc
          return
       endif
    enddo

    if ( chemvar_getid <= 0 ) then
       write(ADM_LOG_FID,*) 'xxx INDEX does not exist =>', tname
       call ADM_proc_stop
    endif

  end function chemvar_getid

end module mod_chemvar
!-------------------------------------------------------------------------------
