!-------------------------------------------------------------------------------
!
!+ Thuburn limiter module
!
!-------------------------------------------------------------------------------
module mod_advlim_thuburn
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for limitting values of q at cell walls in advection
  !       calculation.
  !       (Reference: Thuburn, 1996) 
  !        These are imported from sub[src_update_tracer/mod_src] and 
  !        sub[OPRT_divergence2/mod_oprt].
  !
  !++ Current COrresponding Author: Y.Niwa
  !
  !++ History:
  !      Version   Date      Comment
  !      -----------------------------------------------------------------------
  !      0.00      08-01-24  Imported from mod_src and mod_oprt
  !                11-09-27  T.Seiki: merge optimization for K by M.Terai and RIST
  !                11-11-28  Y.Yamada: Merge Terai-san timer code
  !                                           into the original code.
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: advlim_thuburn_v

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
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine advlim_thuburn_v( &
       q_h, q_h_pl, &
       q,   q_pl,   &
       ck,  ck_pl,  &
       d,   d_pl    )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    implicit none

    real(8), intent(inout) :: q_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: q     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(8), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(8), intent(in)    :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: qin_min    (ADM_gall,   ADM_kall,2)
    real(8) :: qin_min_pl (ADM_gall_pl,ADM_kall,2)
    real(8) :: qin_max    (ADM_gall,   ADM_kall,2)
    real(8) :: qin_max_pl (ADM_gall_pl,ADM_kall,2)

    real(8) :: qout_min   (ADM_gall,   ADM_kall)
    real(8) :: qout_min_pl(ADM_gall_pl,ADM_kall)
    real(8) :: qout_max   (ADM_gall,   ADM_kall)
    real(8) :: qout_max_pl(ADM_gall_pl,ADM_kall)

    real(8) :: q_mp1_min
    real(8) :: q_mp1_min_pl
    real(8) :: q_mp1_max
    real(8) :: q_mp1_max_pl

    real(8) :: c_in
    real(8) :: c_in_pl
    real(8) :: c_out
    real(8) :: c_out_pl

    real(8) :: c_qin_min
    real(8) :: c_qin_min_pl
    real(8) :: c_qin_max
    real(8) :: c_qin_max_pl

    real(8) :: mask, mask1, mask2
    real(8) :: tmp

    integer :: n, k, l
    !---------------------------------------------------------------------------

    !--- inflow & outflow limiter
    do l = 1, ADM_lall

       do k = 1, ADM_kall
          do n = 1, ADM_gall
             qin_min(n,k,1) = CNST_MAX_REAL
             qin_max(n,k,1) =-CNST_MAX_REAL
          enddo
          do n = 1, ADM_gall
             qin_min(n,k,2) = CNST_MAX_REAL
             qin_max(n,k,2) =-CNST_MAX_REAL
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax+1
          do n = 1,ADM_gall
             if ( ck(n,k,l,1) > 0.D0 ) then
                qin_min(n,k-1,2) = min( q(n,k,l), q(n,k-1,l) )
                qin_max(n,k-1,2) = max( q(n,k,l), q(n,k-1,l) )
             else
                qin_min(n,k,  1) = min( q(n,k,l), q(n,k-1,l) )
                qin_max(n,k,  1) = max( q(n,k,l), q(n,k-1,l) )
             endif
          enddo
       enddo

       do k = 1, ADM_kall
          do n = 1, ADM_gall

             q_mp1_min = min( qin_min(n,k,1), qin_min(n,k,2) )
             if( q_mp1_min ==  CNST_MAX_REAL ) q_mp1_min = q(n,k,l)
             q_mp1_min = max( 0.D0, q_mp1_min )

             q_mp1_max = max( qin_max(n,k,1), qin_max(n,k,2) )
             if( q_mp1_max == -CNST_MAX_REAL ) q_mp1_max = q(n,k,l)

             mask1 = 0.5D0 - sign( 0.5D0, ck(n,k,l,1) )
             mask2 = 0.5D0 - sign( 0.5D0, ck(n,k,l,2) )

             c_in  = (        mask1 ) * ck(n,k,l,1) &
                   + (        mask2 ) * ck(n,k,l,2)
             c_out = ( 1.D0 - mask1 ) * ck(n,k,l,1) &
                   + ( 1.D0 - mask2 ) * ck(n,k,l,2)

             c_qin_max = ( mask1 ) * ( ck(n,k,l,1) * qin_max(n,k,1) ) &
                       + ( mask2 ) * ( ck(n,k,l,2) * qin_max(n,k,2) )
             c_qin_min = ( mask1 ) * ( ck(n,k,l,1) * qin_min(n,k,1) ) &
                       + ( mask2 ) * ( ck(n,k,l,2) * qin_min(n,k,2) )

             if ( abs(c_out) <= CNST_EPS_ZERO ) then
                qout_min(n,k) = q(n,k,l)
                qout_max(n,k) = q(n,k,l)
             else
                qout_min(n,k) = ( q(n,k,l)                                  &
                                 - c_qin_max                                &
                                 - q_mp1_max * ( 1.D0-c_in-c_out+d(n,k,l) ) &
                                 ) / c_out
                qout_max(n,k) = ( q(n,k,l)                                  &
                                 - c_qin_min                                &
                                 - q_mp1_min * ( 1.D0-c_in-c_out+d(n,k,l) ) &
                                 ) / c_out
             endif

          enddo
       enddo

       !--- apply limiter
       do k = ADM_kmin, ADM_kmax+1
          do n = 1, ADM_gall
             mask = 0.5D0 - sign( 0.5D0, ck(n,k,l,1) )

             tmp = (        mask ) * min( max( q_h(n,k,l), qin_min(n,k  ,1) ), qin_max(n,k  ,1) ) &
                 + ( 1.D0 - mask ) * min( max( q_h(n,k,l), qin_min(n,k-1,2) ), qin_max(n,k-1,2) )

             q_h(n,k,l) = (        mask ) * max( min( tmp, qout_max(n,k-1) ), qout_min(n,k-1) ) &
                        + ( 1.D0 - mask ) * max( min( tmp, qout_max(n,k  ) ), qout_min(n,k  ) )
          enddo
       enddo

    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl

          do k = 1, ADM_kall
             do n = 1, ADM_gall_pl
                qin_min_pl(n,k,1) =  CNST_MAX_REAL
                qin_max_pl(n,k,1) = -CNST_MAX_REAL
             enddo
             do n = 1, ADM_gall_pl
                qin_min_pl(n,k,2) =  CNST_MAX_REAL
                qin_max_pl(n,k,2) = -CNST_MAX_REAL
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do n = 1, ADM_gall_pl
                if ( ck_pl(n,k,l,1) > 0.D0 ) then
                   qin_min_pl(n,k-1,2) = min( q_pl(n,k,l), q_pl(n,k-1,l) )
                   qin_max_pl(n,k-1,2) = max( q_pl(n,k,l), q_pl(n,k-1,l) )
                else
                   qin_min_pl(n,k,  1) = min( q_pl(n,k,l), q_pl(n,k-1,l) )
                   qin_max_pl(n,k,  1) = max( q_pl(n,k,l), q_pl(n,k-1,l) )
                endif
             enddo
          enddo

          do k = 1, ADM_kall
             do n = 1, ADM_gall_pl

                q_mp1_min_pl = min( qin_min_pl(n,k,1), qin_min_pl(n,k,2) )
                if( q_mp1_min_pl ==  CNST_MAX_REAL ) q_mp1_min_pl = q_pl(n,k,l)
                q_mp1_min_pl = max( 0.D0, q_mp1_min_pl)

                q_mp1_max_pl = max( qin_max_pl(n,k,1), qin_max_pl(n,k,2) )
                if( q_mp1_max_pl == -CNST_MAX_REAL ) q_mp1_max_pl = q_pl(n,k,l)

                mask1 = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,1) )
                mask2 = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,2) )

                c_in_pl  = (        mask1 ) * ck_pl(n,k,l,1) &
                         + (        mask2 ) * ck_pl(n,k,l,2)
                c_out_pl = ( 1.D0 - mask1 ) * ck_pl(n,k,l,1) &
                         + ( 1.D0 - mask2 ) * ck_pl(n,k,l,2)

                c_qin_max_pl = ( mask1 ) * ( ck_pl(n,k,l,1) * qin_max_pl(n,k,1) ) &
                             + ( mask2 ) * ( ck_pl(n,k,l,2) * qin_max_pl(n,k,2) )
                c_qin_min_pl = ( mask1 ) * ( ck_pl(n,k,l,1) * qin_min_pl(n,k,1) ) &
                             + ( mask2 ) * ( ck_pl(n,k,l,2) * qin_min_pl(n,k,2) )

                if ( abs(c_out_pl) < CNST_EPS_ZERO ) then
                   qout_min_pl(n,k) = q_pl(n,k,l)
                   qout_max_pl(n,k) = q_pl(n,k,l)
                else
                   qout_min_pl(n,k) = ( q_pl(n,k,l)                                          &
                                      - c_qin_max_pl                                         &
                                      - q_mp1_max_pl * ( 1.D0-c_in_pl-c_out_pl+d_pl(n,k,l) ) &
                                      ) / c_out_pl
                   qout_max_pl(n,k) = ( q_pl(n,k,l)                                          &
                                      - c_qin_min_pl                                         &
                                      - q_mp1_min_pl * ( 1.D0-c_in_pl-c_out_pl+d_pl(n,k,l) ) &
                                      ) / c_out_pl
                endif

             enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
             do n = 1, ADM_gall_pl
                mask = 0.5D0 - sign( 0.5D0, ck_pl(n,k,l,1) )

                tmp = (        mask ) * min( max( q_h_pl(n,k,l), qin_min_pl(n,k  ,1) ), qin_max_pl(n,k  ,1) ) &
                    + ( 1.D0 - mask ) * min( max( q_h_pl(n,k,l), qin_min_pl(n,k-1,2) ), qin_max_pl(n,k-1,2) )

                q_h_pl(n,k,l) = (        mask ) * max( min( tmp, qout_max_pl(n,k-1) ), qout_min_pl(n,k-1) ) &
                              + ( 1.D0 - mask ) * max( min( tmp, qout_max_pl(n,k  ) ), qout_min_pl(n,k  ) )
             enddo
          enddo

       enddo
    endif

    return
  end subroutine advlim_thuburn_v

end module mod_advlim_thuburn
!-------------------------------------------------------------------------------------
