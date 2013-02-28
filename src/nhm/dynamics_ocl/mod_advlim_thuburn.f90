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
  !                                   (see !OCL directive )
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
  public :: advlim_thuburn_h

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

#ifdef _FJTIMER_
call timer_sta(2700)
#endif

    !--- inflow & outflow limiter
!OCL SERIAL
    do l = 1, ADM_lall

!OCL PARALLEL
       do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do n = 1, ADM_gall
             qin_min(n,k,1) = CNST_MAX_REAL
             qin_max(n,k,1) =-CNST_MAX_REAL
          enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION
          do n = 1, ADM_gall
             qin_min(n,k,2) = CNST_MAX_REAL
             qin_max(n,k,2) =-CNST_MAX_REAL
          enddo
       enddo

       do k = ADM_kmin, ADM_kmax+1
!OCL PARALLEL
!OCL SIMD
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

!OCL PARALLEL
       do k = 1, ADM_kall
!OCL SIMD
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
!OCL SERIAL
       do k = ADM_kmin, ADM_kmax+1
!OCL PARALLEL
!OCL SIMD
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
!OCL SERIAL
       do l = 1, ADM_lall_pl

!OCL PARALLEL
          do k = 1, ADM_kall
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do n = 1, ADM_gall_pl
                qin_min_pl(n,k,1) =  CNST_MAX_REAL
                qin_max_pl(n,k,1) = -CNST_MAX_REAL
             enddo
!OCL XFILL, PREFETCH_STRONG, LOOP_NOFUSION, UNROLL('full')
             do n = 1, ADM_gall_pl
                qin_min_pl(n,k,2) =  CNST_MAX_REAL
                qin_max_pl(n,k,2) = -CNST_MAX_REAL
             enddo
          enddo

!OCL SERIAL
          do k = ADM_kmin, ADM_kmax+1
!OCL PARALLEL
!OCL SIMD
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

!OCL PARALLEL
          do k = 1, ADM_kall
!OCL SIMD
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
!OCL PARALLEL
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

#ifdef _FJTIMER_
call timer_end(2700)
#endif

    return
  end subroutine advlim_thuburn_v

  !-----------------------------------------------------------------------------------
  subroutine advlim_thuburn_h( &
       qa,   qa_pl,   &
       q,    q_pl,    &
       c,    c_pl,    &
       d,    d_pl     &
       )
    !
    use mod_adm, only :   &
         !--- public parameters
         ADM_W,           &
         ADM_prc_pl,      &
         ADM_AI,          &
         ADM_AIJ,         &
         ADM_AJ,          &
         ADM_lall_pl,     &
         ADM_GMIN_PL,     &
         ADM_GSLF_PL,     &
         ADM_gall_pl,     &
         ADM_GMAX_PL,     &
         !--- public variables
         ADM_prc_me,      &
         ADM_prc_tab,     &
         ADM_rgn_vnum,    &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall,        &
         ADM_comm_run_world
    use mod_comm, only : &
         COMM_data_transfer
    use mod_cnst, only : &
         CNST_MAX_REAL,  &
         CNST_EPS_ZERO
    !
    implicit none
    !
    real(8), intent(inout) :: qa(ADM_gall,   ADM_kall,ADM_lall   ,ADM_AI:ADM_AJ)
    real(8), intent(inout) :: qa_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: q(ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in) :: q_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: c(ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(8), intent(in) :: c_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: d(ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in) :: d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !    logical, intent(in), optional :: NON_NEG=.false.
    !
    integer, parameter :: q_out_k_min=1
    integer, parameter :: q_out_k_max=2
    real(8)  :: wrk(ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(8)  :: wrk_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    !
    real(8)  :: q_in_min(ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(8)  :: q_in_min_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(8)  :: q_in_max(ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(8)  :: q_in_max_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    !
    real(8)  :: q_m1_k_min(ADM_gall)
    real(8)  :: q_m1_k_min_pl
    real(8)  :: q_m1_k_max(ADM_gall)
    real(8)  :: q_m1_k_max_pl
    real(8)  :: c_in_sum(ADM_gall)
    real(8)  :: c_in_sum_pl
    real(8)  :: c_out_sum(ADM_gall)
    real(8)  :: c_out_sum_pl
    real(8)  :: c_qin_sum_max(ADM_gall)
    real(8)  :: c_qin_sum_max_pl
    real(8)  :: c_qin_sum_min(ADM_gall)
    real(8)  :: c_qin_sum_min_pl
    !
    integer :: l,n,k
    integer :: rgnid
    !
    integer :: np1(ADM_gall_pl)
    integer :: nm1(ADM_gall_pl)
    !
    integer :: nstart,nend

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !
    !
    nm1(ADM_GMIN_PL) = ADM_GMAX_PL
    do n=ADM_GMIN_PL+1,ADM_GMAX_PL
       nm1(n) = n-1
    enddo
    !
    do n=ADM_GMIN_PL,ADM_GMAX_PL-1
       np1(n) = n+1
    enddo
    np1(ADM_GMAX_PL) = ADM_GMIN_PL
    !
    !--- calculation of inflow limiter
!OCL SERIAL
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
!OCL PARALLEL
       do k = 1, ADM_kall
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          q_in_min(:,k,l,:) = CNST_MAX_REAL
          q_in_max(:,k,l,:) =-CNST_MAX_REAL
          !
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax,  ADM_gmax  )
!OCL SIMD
          do n=nstart,nend
             if(c(n  ,k,l,1)<=0.0D0) then
                q_in_min(n,k,l,1) = min(q(n,k,l),q(n+1,k,l),q(n+1+ADM_gall_1d,k,l),q(n-ADM_gall_1d,k,l))
                q_in_max(n,k,l,1) = max(q(n,k,l),q(n+1,k,l),q(n+1+ADM_gall_1d,k,l),q(n-ADM_gall_1d,k,l))
             else
                q_in_min(n+1,k,l,4) = min(q(n,k,l),q(n+1,k,l),q(n+1+ADM_gall_1d,k,l),q(n-ADM_gall_1d,k,l))
                q_in_max(n+1,k,l,4) = max(q(n,k,l),q(n+1,k,l),q(n+1+ADM_gall_1d,k,l),q(n-ADM_gall_1d,k,l))
             endif
          enddo
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax,  ADM_gmax  )
!OCL SIMD
          do n=nstart,nend
             !
             if(c(n  ,k,l,2)<=0.0D0) then
                q_in_min(n,k,l,2) = min(q(n,k,l),q(n+1+ADM_gall_1d,k,l),q(n+1,k,l),q(n+ADM_gall_1d,k,l))
                q_in_max(n,k,l,2) = max(q(n,k,l),q(n+1+ADM_gall_1d,k,l),q(n+1,k,l),q(n+ADM_gall_1d,k,l))
             else
                q_in_min(n+1+ADM_gall_1d,k,l,5) = min(q(n,k,l),q(n+1+ADM_gall_1d,k,l),q(n+1,k,l),q(n+ADM_gall_1d,k,l))
                q_in_max(n+1+ADM_gall_1d,k,l,5) = max(q(n,k,l),q(n+1+ADM_gall_1d,k,l),q(n+1,k,l),q(n+ADM_gall_1d,k,l))
             endif
             !
          enddo
          !
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
!OCL SIMD
          do n=nstart,nend
             if(c(n  ,k,l,3)<=0.0D0) then
                q_in_min(n,k,l,3) = min(q(n,k,l),q(n+ADM_gall_1d,k,l),q(n+1+ADM_gall_1d,k,l),q(n-1,k,l))
                q_in_max(n,k,l,3) = max(q(n,k,l),q(n+ADM_gall_1d,k,l),q(n+1+ADM_gall_1d,k,l),q(n-1,k,l))
             else
                q_in_min(n+ADM_gall_1d,k,l,6) = min(q(n,k,l),q(n+ADM_gall_1d,k,l),q(n+1+ADM_gall_1d,k,l),q(n-1,k,l))
                q_in_max(n+ADM_gall_1d,k,l,6) = max(q(n,k,l),q(n+ADM_gall_1d,k,l),q(n+1+ADM_gall_1d,k,l),q(n-1,k,l))
             endif
             !
          enddo
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             n = suf(ADM_gmin-1,ADM_gmin-1)
             if(c(n  ,k,l,2)<=0.0D0) then
                q_in_min(n,k,l,2) = min(q(n,k,l),q(n+1+ADM_gall_1d,k,l),q(n+1+1+ADM_gall_1d,k,l),q(n+ADM_gall_1d,k,l))
                q_in_max(n,k,l,2) = max(q(n,k,l),q(n+1+ADM_gall_1d,k,l),q(n+1+1+ADM_gall_1d,k,l),q(n+ADM_gall_1d,k,l))
             else
                q_in_min(n+1+ADM_gall_1d,k,l,5) = min(q(n,k,l),q(n+1+ADM_gall_1d,k,l),&
                     q(n+1+1+ADM_gall_1d,k,l),q(n+ADM_gall_1d,k,l))
                q_in_max(n+1+ADM_gall_1d,k,l,5) = max(q(n,k,l),q(n+1+ADM_gall_1d,k,l),&
                     q(n+1+1+ADM_gall_1d,k,l),q(n+ADM_gall_1d,k,l))
             endif
          endif
          !
       enddo
       !
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1,ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
             q_in_min_pl(:,k,l,:) = CNST_MAX_REAL
             q_in_max_pl(:,k,l,:) =-CNST_MAX_REAL
!OCL SIMD
             do n=ADM_GMIN_PL,ADM_GMAX_PL
                if(c_pl(n,k,l)<=0.0D0) then
                   q_in_min_pl(n,k,l,1) = min(q_pl(ADM_GSLF_PL,k,l),q_pl(n,k,l),q_pl(nm1(n),k,l),q_pl(np1(n),k,l))
                   q_in_max_pl(n,k,l,1) = max(q_pl(ADM_GSLF_PL,k,l),q_pl(n,k,l),q_pl(nm1(n),k,l),q_pl(np1(n),k,l))
                else
                   q_in_min_pl(n,k,l,2) = min(q_pl(ADM_GSLF_PL,k,l),q_pl(n,k,l),q_pl(nm1(n),k,l),q_pl(np1(n),k,l))
                   q_in_max_pl(n,k,l,2) = max(q_pl(ADM_GSLF_PL,k,l),q_pl(n,k,l),q_pl(nm1(n),k,l),q_pl(np1(n),k,l))
                endif
             enddo

          enddo
       enddo
    endif
    !
    !--- calcluation outflow limiter
!OCL SERIAL
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
!OCL PARALLEL
       do k = 1, ADM_kall
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax+1,ADM_gmax+1)
!OCL SIMD
          do n = nstart,nend
             q_m1_k_min(n) = min(q_in_min(n,k,l,1),q_in_min(n,k,l,2),q_in_min(n,k,l,3),&
                  q_in_min(n,k,l,4),q_in_min(n,k,l,5),q_in_min(n,k,l,6))
             if(q_m1_k_min(n)==CNST_MAX_REAL) q_m1_k_min(n) = q(n,k,l)
             !        q_m1_k_min(n) = max(SMALL_ZERO,q_m1_k_min(n))
             q_m1_k_max(n) = max(q_in_max(n,k,l,1),q_in_max(n,k,l,2),q_in_max(n,k,l,3),&
                  q_in_max(n,k,l,4),q_in_max(n,k,l,5),q_in_max(n,k,l,6))
             if(q_m1_k_max(n)==-CNST_MAX_REAL) q_m1_k_max(n) = q(n,k,l)
             !
          enddo
          !
          nstart = suf(ADM_gmin,  ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
!OCL SIMD
          do n = nstart,nend
             c_in_sum(n) &
                  = (0.5D0-sign(0.5D0,c(n,k,l,1)))*c(n,k,l,1)&
                  + (0.5D0-sign(0.5D0,c(n,k,l,2)))*c(n,k,l,2)&
                  + (0.5D0-sign(0.5D0,c(n,k,l,3)))*c(n,k,l,3)&
                  + (0.5D0-sign(0.5D0,c(n,k,l,4)))*c(n,k,l,4)&
                  + (0.5D0-sign(0.5D0,c(n,k,l,5)))*c(n,k,l,5)&
                  + (0.5D0-sign(0.5D0,c(n,k,l,6)))*c(n,k,l,6)
             c_out_sum(n) &
                  = (0.5D0+sign(0.5D0,c(n,k,l,1)))*c(n,k,l,1)&
                  + (0.5D0+sign(0.5D0,c(n,k,l,2)))*c(n,k,l,2)&
                  + (0.5D0+sign(0.5D0,c(n,k,l,3)))*c(n,k,l,3)&
                  + (0.5D0+sign(0.5D0,c(n,k,l,4)))*c(n,k,l,4)&
                  + (0.5D0+sign(0.5D0,c(n,k,l,5)))*c(n,k,l,5)&
                  + (0.5D0+sign(0.5D0,c(n,k,l,6)))*c(n,k,l,6)
             !
             c_qin_sum_max(n) &
                  = (0.5D0-sign(0.5D0,c(n,k,l,1)))*(c(n,k,l,1)*q_in_max(n,k,l,1))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,2)))*(c(n,k,l,2)*q_in_max(n,k,l,2))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,3)))*(c(n,k,l,3)*q_in_max(n,k,l,3))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,4)))*(c(n,k,l,4)*q_in_max(n,k,l,4))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,5)))*(c(n,k,l,5)*q_in_max(n,k,l,5))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,6)))*(c(n,k,l,6)*q_in_max(n,k,l,6))

             c_qin_sum_min(n) &
                  = (0.5D0-sign(0.5D0,c(n,k,l,1)))*(c(n,k,l,1)*q_in_min(n,k,l,1))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,2)))*(c(n,k,l,2)*q_in_min(n,k,l,2))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,3)))*(c(n,k,l,3)*q_in_min(n,k,l,3))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,4)))*(c(n,k,l,4)*q_in_min(n,k,l,4))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,5)))*(c(n,k,l,5)*q_in_min(n,k,l,5))&
                  + (0.5D0-sign(0.5D0,c(n,k,l,6)))*(c(n,k,l,6)*q_in_min(n,k,l,6))
             if(abs(c_out_sum(n))<CNST_EPS_ZERO) then
                wrk(n,k,l,q_out_k_min) = q(n,k,l)
                wrk(n,k,l,q_out_k_max) = q(n,k,l)
             else
                wrk(n,k,l,q_out_k_min) = ( &
                     q(n,k,l)-c_qin_sum_max(n)&
                     -q_m1_k_max(n)*(1.0D0-c_in_sum(n)-c_out_sum(n)+d(n,k,l)) &
                     )/c_out_sum(n)
                wrk(n,k,l,q_out_k_max) = ( &
                     q(n,k,l)-c_qin_sum_min(n)&
                     -q_m1_k_min(n)*(1.0D0-c_in_sum(n)-c_out_sum(n)+d(n,k,l)) &
                     )/c_out_sum(n)
             endif
          enddo
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1,ADM_lall_pl
!OCL PARALLEL,SIMD
          do k = 1, ADM_kall
             q_m1_k_min_pl &
                  = min(q_in_min_pl(ADM_GMIN_PL,k,l,1),q_in_min_pl(ADM_GMIN_PL+1,k,l,1),&
                  q_in_min_pl(ADM_GMIN_PL+2,k,l,1),q_in_min_pl(ADM_GMIN_PL+3,k,l,1),q_in_min_pl(ADM_GMIN_PL+4,k,l,1))
             if(q_m1_k_min_pl== CNST_MAX_REAL) q_m1_k_min_pl=q_pl(ADM_GSLF_PL,k,l)
             !          q_m1_k_min_pl = max(SMALL_ZERO,q_m1_k_min_pl)
             q_m1_k_max_pl &
                  = max(q_in_min_pl(ADM_GMIN_PL,k,l,1),q_in_min_pl(ADM_GMIN_PL+1,k,l,1),&
                  q_in_min_pl(ADM_GMIN_PL+2,k,l,1),q_in_min_pl(ADM_GMIN_PL+3,k,l,1),q_in_min_pl(ADM_GMIN_PL+4,k,l,1))
             if(q_m1_k_max_pl==-CNST_MAX_REAL) q_m1_k_max_pl=q_pl(ADM_GSLF_PL,k,l)
             !
             c_in_sum_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL  ,k,l)))*c_pl(ADM_GMIN_PL  ,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+1,k,l)))*c_pl(ADM_GMIN_PL+1,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+2,k,l)))*c_pl(ADM_GMIN_PL+2,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+3,k,l)))*c_pl(ADM_GMIN_PL+3,k,l)&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+4,k,l)))*c_pl(ADM_GMIN_PL+4,k,l)
             c_out_sum_pl &
                  = (0.5D0+sign(0.5D0,c_pl(ADM_GMIN_PL  ,k,l)))*c_pl(ADM_GMIN_PL  ,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_GMIN_PL+1,k,l)))*c_pl(ADM_GMIN_PL+1,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_GMIN_PL+2,k,l)))*c_pl(ADM_GMIN_PL+2,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_GMIN_PL+3,k,l)))*c_pl(ADM_GMIN_PL+3,k,l)&
                  + (0.5D0+sign(0.5D0,c_pl(ADM_GMIN_PL+4,k,l)))*c_pl(ADM_GMIN_PL+4,k,l)
             c_qin_sum_max_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL  ,k,l)))*(c_pl(ADM_GMIN_PL  ,k,l)*q_in_max_pl(ADM_GMIN_PL  ,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+1,k,l)))*(c_pl(ADM_GMIN_PL+1,k,l)*q_in_max_pl(ADM_GMIN_PL+1,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+2,k,l)))*(c_pl(ADM_GMIN_PL+2,k,l)*q_in_max_pl(ADM_GMIN_PL+2,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+3,k,l)))*(c_pl(ADM_GMIN_PL+3,k,l)*q_in_max_pl(ADM_GMIN_PL+3,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+4,k,l)))*(c_pl(ADM_GMIN_PL+4,k,l)*q_in_max_pl(ADM_GMIN_PL+4,k,l,1))
             c_qin_sum_min_pl &
                  = (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL  ,k,l)))*(c_pl(ADM_GMIN_PL  ,k,l)*q_in_min_pl(ADM_GMIN_PL  ,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+1,k,l)))*(c_pl(ADM_GMIN_PL+1,k,l)*q_in_min_pl(ADM_GMIN_PL+1,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+2,k,l)))*(c_pl(ADM_GMIN_PL+2,k,l)*q_in_min_pl(ADM_GMIN_PL+2,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+3,k,l)))*(c_pl(ADM_GMIN_PL+3,k,l)*q_in_min_pl(ADM_GMIN_PL+3,k,l,1))&
                  + (0.5D0-sign(0.5D0,c_pl(ADM_GMIN_PL+4,k,l)))*(c_pl(ADM_GMIN_PL+4,k,l)*q_in_min_pl(ADM_GMIN_PL+4,k,l,1))
             !
             if(abs(c_out_sum_pl)<CNST_EPS_ZERO) then
                wrk_pl(ADM_GSLF_PL,k,l,q_out_k_min) = q_pl(ADM_GSLF_PL,k,l)
                wrk_pl(ADM_GSLF_PL,k,l,q_out_k_max) = q_pl(ADM_GSLF_PL,k,l)
             else
                wrk_pl(ADM_GSLF_PL,k,l,q_out_k_min) = ( &
                     q_pl(ADM_GSLF_PL,k,l)-c_qin_sum_max_pl&
                     -q_m1_k_max_pl*(1.0D0-c_in_sum_pl-c_out_sum_pl+d_pl(ADM_GSLF_PL,k,l)) &
                     )/c_out_sum_pl
                wrk_pl(ADM_GSLF_PL,k,l,q_out_k_max) = ( &
                     q_pl(ADM_GSLF_PL,k,l)-c_qin_sum_min_pl&
                     -q_m1_k_min_pl*(1.0D0-c_in_sum_pl-c_out_sum_pl+d_pl(ADM_GSLF_PL,k,l)) &
                     )/c_out_sum_pl
             endif
          enddo
       enddo
    endif
    !

    call COMM_data_transfer(wrk,wrk_pl)

    !
    !---- apply inflow limiter
!OCL SERIAL
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
!OCL PARALLEL
       do k = 1, ADM_kall
          !
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n=nstart,nend
             qa(n,k,l,ADM_AI) &
                  =(0.5D0-sign(0.5D0,c(n,k,l,1)))&
                  *min(max(qa(n,k,l,ADM_AI),q_in_min(n,k,l,1)),q_in_max(n,k,l,1))&
                  +(0.5D0+sign(0.5D0,c(n,k,l,1)))&
                  *min(max(qa(n,k,l,ADM_AI),q_in_min(n+1,k,l,4)),q_in_max(n+1,k,l,4))
	  enddo
	  !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          do n=nstart,nend
             qa(n,k,l,ADM_AIJ) &
                  =(0.5D0-sign(0.5D0,c(n,k,l,2)))&
                  *min(max(qa(n,k,l,ADM_AIJ),q_in_min(n,k,l,2)),q_in_max(n,k,l,2))&
                  +(0.5D0+sign(0.5D0,c(n,k,l,2)))&
                  *min(max(qa(n,k,l,ADM_AIJ),q_in_min(n+1+ADM_gall_1d,k,l,5)),q_in_max(n+1+ADM_gall_1d,k,l,5))
          enddo
	  !
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          do n=nstart,nend
             qa(n,k,l,ADM_AJ) &
                  =(0.5D0-sign(0.5D0,c(n,k,l,3)))&
                  *min(max(qa(n,k,l,ADM_AJ),q_in_min(n,k,l,3)),q_in_max(n,k,l,3))&
                  +(0.5D0+sign(0.5D0,c(n,k,l,3)))&
                  *min(max(qa(n,k,l,ADM_AJ),q_in_min(n+ADM_gall_1d,k,l,6)),q_in_max(n+ADM_gall_1d,k,l,6))
          enddo
       enddo
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1,ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
             do n=ADM_GMIN_PL,ADM_GMAX_PL
                qa_pl(n,k,l) &
                  =(0.5D0-sign(0.5D0,c_pl(n,k,l)))&
                  *min(max(qa_pl(n,k,l),q_in_min_pl(n,k,l,1)),q_in_max_pl(n,k,l,1))&
                  +(0.5D0+sign(0.5D0,c_pl(n,k,l)))&
                  *min(max(qa_pl(n,k,l),q_in_min_pl(n,k,l,2)),q_in_max_pl(n,k,l,2))
             enddo
          enddo
       enddo
    endif
    !
    !
    !---- apply outflow limitter
!OCL SERIAL
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
!OCL PARALLEL
       do k = 1, ADM_kall
          !
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             qa(n,k,l,ADM_AI) &
                  =(0.5D0-sign(0.5D0,c(n,k,l,1)))&
                  *max(min(qa(n,k,l,ADM_AI),wrk(n+1,k,l,q_out_k_max)),wrk(n+1,k,l,q_out_k_min))&
                  +(0.5D0+sign(0.5D0,c(n,k,l,1)))&
                  *max(min(qa(n,k,l,ADM_AI),wrk(n,k,l,q_out_k_max)),wrk(n,k,l,q_out_k_min))
          enddo
	  !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          do n=nstart,nend
             qa(n,k,l,ADM_AIJ) &
                  =(0.5D0-sign(0.5D0,c(n,k,l,2)))&
                  *max(min(qa(n,k,l,ADM_AIJ),wrk(n+1+ADM_gall_1d,k,l,q_out_k_max)),wrk(n+1+ADM_gall_1d,k,l,q_out_k_min))&
                  +(0.5D0+sign(0.5D0,c(n,k,l,2)))&
                  *max(min(qa(n,k,l,ADM_AIJ),wrk(n,k,l,q_out_k_max)),wrk(n,k,l,q_out_k_min))
	  enddo
	  !
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          do n=nstart,nend
             qa(n,k,l,ADM_AJ) &
                  =(0.5D0-sign(0.5D0,c(n,k,l,3)))&
                  *max(min(qa(n,k,l,ADM_AJ),wrk(n+ADM_gall_1d,k,l,q_out_k_max)),wrk(n+ADM_gall_1d,k,l,q_out_k_min))&
                  +(0.5D0+sign(0.5D0,c(n,k,l,3)))&
                  *max(min(qa(n,k,l,ADM_AJ),wrk(n,k,l,q_out_k_max)),wrk(n,k,l,q_out_k_min))
          enddo
       enddo
    enddo
    if(ADM_prc_me==ADM_prc_pl) then
!OCL SERIAL
       do l=1,ADM_lall_pl
!OCL PARALLEL
          do k = 1, ADM_kall
             do n = ADM_GMIN_PL,ADM_GMAX_PL
                qa_pl(n,k,l)&
                     =(0.5D0-sign(0.5D0,c_pl(n,k,l)))&
                     *max(min(qa_pl(n,k,l),wrk_pl(n,k,l,q_out_k_max)),wrk_pl(n,k,l,q_out_k_min))&
                     +(0.5D0+sign(0.5D0,c_pl(n,k,l)))&
                     *max(min(qa_pl(n,k,l),wrk_pl(ADM_GSLF_PL,k,l,q_out_k_max)),wrk_pl(ADM_GSLF_PL,k,l,q_out_k_min))
             enddo
          enddo
       enddo
    endif
    !
    return
    !
  end subroutine advlim_thuburn_h
  !-----------------------------------------------------------------------------------
end module mod_advlim_thuburn
!-------------------------------------------------------------------------------------
