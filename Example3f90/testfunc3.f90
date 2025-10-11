!test AD module: 2D Twiss parameters (alpha,beta) tracking through a FODO lattice
       program testfunc
       use TPSAad
       implicit none
       !# of control variables
       integer, parameter :: noi = 7
       !7 control variables
       type (dtpsad) :: x1,x2,x3,x4,x5,x6,x7
       !2D Twiss parameters (beta_x,alpha_x,beta_y,alpha_y)
       type (dtpsad), dimension(4) :: twiss 
       integer :: i, nperd
       real*8 :: xk

       !initialize the AD module
       call dtpsad_Initialize(noi)

       !initialize each lattice variable
       !drift length
       call assign_tpsad(x1,0.2d0,1)
       !quad length
       call assign_tpsad(x2,0.1d0,2)
       !quad strengh (k)
       call assign_tpsad(x3,29.6d0,3)
       !drift length
       call assign_tpsad(x4,0.4d0,4)
       !quad length
       call assign_tpsad(x5,0.1d0,5)
       !quad strengh (k)
       call assign_tpsad(x6,-29.6d0,6)
       !drift length
       call assign_tpsad(x7,0.2d0,7)
       
       !initial Twiss parameters
       do i = 1, 4
         twiss(1) = 1.d0
         twiss(2) = 0.5d0
         twiss(3) = 1.d0
         twiss(4) = -0.5d0
       enddo

       !loop through nperd of the FODO lattice
       nperd = 1
       do i = 1, nperd
         !drift
         xk = 0.0d0
         call quadmap(twiss,x1,x2,xk)
         !quad
         xk = x3%map(1)
         call quadmap(twiss,x2,x3,xk)
         !drift
         xk = 0.0d0
         call quadmap(twiss,x4,x2,xk)
         !quad
         xk = x6%map(1)
         call quadmap(twiss,x5,x6,xk)
         !drift
         xk = 0.0d0
         call quadmap(twiss,x7,x2,xk)
       enddo

       print*,"Beta_x and its derivatives w.r.t. 7 variables:"
       print*,twiss(1)%map(1:noi+1)
       print*,"Alpha_x and its derivatives w.r.t. 7 variables:"
       print*,twiss(2)%map(1:noi+1)
       print*,"Beta_y and its derivatives w.r.t. 7 variables:"
       print*,twiss(3)%map(1:noi+1)
       print*,"Alpha_y and its derivatives w.r.t. 7 variables:"
       print*,twiss(4)%map(1:noi+1)

       end program

       subroutine quadmap(tpsaPtc,tau,outin,xk)
        use TPSAad
        implicit none
        type (dtpsad), intent(in) :: tau,outin
        type (dtpsad), dimension(4) , intent(inout) ::tpsaPtc
        real*8, intent(in) :: xk
        type (dtpsad), dimension(4) ::tpsaPtcTemp
        type (dtpsad) :: co,si,ch,sh,sqk,gammax,gammay

        gammax = (1.0d0+pow(tpsaPtc(2),2))/tpsaPtc(1)
        gammay = (1.0d0+pow(tpsaPtc(4),2))/tpsaPtc(3)
        if(xk.eq.0.0d0)then
          tpsaPtc(1)=tpsaPtc(1) - 2*tpsaPtc(2)*tau + tau*tau*gammax
          tpsaPtc(3)=tpsaPtc(3) - 2*tpsaPtc(4)*tau + tau*tau*gammay
        elseif(xk.gt.0.0d0)then
          sqk=sqrt(outin)
          co=cos(sqk*tau)
          si=sin(sqk*tau)
          ch=cosh(sqk*tau)
          sh=sinh(sqk*tau)
          tpsaPtcTemp(1)=tpsaPtc(1)*co*co - &
                        2*tpsaPtc(2)*co*si/sqk + pow(si/sqk,2)*gammax
          tpsaPtcTemp(2)=tpsaPtc(1)*co*si*sqk + &
                         tpsaPtc(2)*(co*co-si*si) - co*si/sqk*gammax
          tpsaPtcTemp(3)=tpsaPtc(3)*ch*ch - &
                         2*tpsaPtc(4)*ch*sh/sqk + pow(sh/sqk,2)*gammay
          tpsaPtcTemp(4)=-tpsaPtc(3)*ch*sh*sqk + &
                          tpsaPtc(4)*(ch*ch+sh*sh) - ch*sh/sqk*gammay
          tpsaPtc(1) = tpsaPtcTemp(1)
          tpsaPtc(2) = tpsaPtcTemp(2)
          tpsaPtc(3) = tpsaPtcTemp(3)
          tpsaPtc(4) = tpsaPtcTemp(4)
        elseif(xk.lt.0.0d0)then
          sqk=sqrt(-outin)
          co=cos(sqk*tau)
          si=sin(sqk*tau)
          ch=cosh(sqk*tau)
          sh=sinh(sqk*tau)
          tpsaPtcTemp(1)=tpsaPtc(1)*ch*ch - &
                         2*tpsaPtc(2)*ch*sh/sqk + pow(sh/sqk,2)*gammax
          tpsaPtcTemp(2)=-tpsaPtc(1)*ch*sh*sqk + &
                          tpsaPtc(2)*(ch*ch+sh*sh) - ch*sh/sqk*gammax
          tpsaPtcTemp(3)=tpsaPtc(3)*co*co - &
                        2*tpsaPtc(4)*co*si/sqk + pow(si/sqk,2)*gammay
          tpsaPtcTemp(4)=tpsaPtc(3)*co*si*sqk + &
                         tpsaPtc(4)*(co*co-si*si) - co*si/sqk*gammay
          tpsaPtc(1) = tpsaPtcTemp(1)
          tpsaPtc(2) = tpsaPtcTemp(2)
          tpsaPtc(3) = tpsaPtcTemp(3)
          tpsaPtc(4) = tpsaPtcTemp(4)
        endif

       end subroutine quadmap
