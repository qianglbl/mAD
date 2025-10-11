!test AD module a 2D particle tracking through a FODO lattice
       program testfunc
       use TPSAad
       implicit none
       !# of control variables
       integer, parameter :: noi = 7
       !7 control variables
       type (dtpsad) :: x1,x2,x3,x4,x5,x6,x7
       !particle coordinates (x,px,y,py)
       type (dtpsad), dimension(4) :: pt 
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
       
       !initial particle coordinates
       do i = 1, 4
         pt(1) = 1.d-3
         pt(2) = 1.d-3
         pt(3) = -1.d-3
         pt(4) = -1.d-3
       enddo

       !loop through nperd of the FODO lattice
       nperd = 1
       do i = 1, nperd
         !drift
         xk = 0.0d0
         call quadmap(pt,x1,x2,xk)
         !quad
         xk = x3%map(1)
         call quadmap(pt,x2,x3,xk)
         !drift
         xk = 0.0d0
         call quadmap(pt,x4,x2,xk)
         !quad
         xk = x6%map(1)
         call quadmap(pt,x5,x6,xk)
         !drift
         xk = 0.0d0
         call quadmap(pt,x7,x2,xk)
       enddo

       print*,"X coordinate and its derivatives w.r.t. 7 variables:"
       print*,pt(1)%map(1:noi+1)
       print*,"Px coordinate and its derivatives w.r.t. 7 variables:"
       print*,pt(2)%map(1:noi+1)
       print*,"Y coordinate and its derivatives w.r.t. 7 variables:"
       print*,pt(3)%map(1:noi+1)
       print*,"Py coordinate and its derivatives w.r.t. 7 variables:"
       print*,pt(4)%map(1:noi+1)

       end program

       subroutine quadmap(tpsaPtc,tau,outin,xk)
        use TPSAad
        implicit none
        type (dtpsad), intent(in) :: tau,outin
        type (dtpsad), dimension(4) , intent(inout) ::tpsaPtc
        real*8, intent(in) :: xk
        type (dtpsad), dimension(4) ::tpsaPtcTemp
        type (dtpsad) :: co,si,ch,sh,sqk

        if(xk.eq.0.0d0)then
          tpsaPtc(1)=tpsaPtc(1)+ tpsaPtc(2)*tau
          tpsaPtc(3)=tpsaPtc(3)+ tpsaPtc(4)*tau
        elseif(xk.gt.0.0d0)then
          sqk=sqrt(outin)
          co=cos(sqk*tau)
          si=sin(sqk*tau)
          ch=cosh(sqk*tau)
          sh=sinh(sqk*tau)
          tpsaPtcTemp(1)=tpsaPtc(1)*co+ &
                           tpsaPtc(2)*si/sqk
          tpsaPtcTemp(2)=-tpsaPtc(1)*si*sqk+ &
                            tpsaPtc(2)*co
          tpsaPtcTemp(3)=tpsaPtc(3)*ch+ &
                           tpsaPtc(4)*sh/sqk
          tpsaPtcTemp(4)=tpsaPtc(3)*sh*sqk+ &
                           tpsaPtc(4)*ch
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
          tpsaPtcTemp(1)=tpsaPtc(1)*ch+ &
                         tpsaPtc(2)*sh/sqk
          tpsaPtcTemp(2)=tpsaPtc(1)*sh*sqk+ &
                         tpsaPtc(2)*ch
          tpsaPtcTemp(3)=tpsaPtc(3)*co+ &
                         tpsaPtc(4)*si/sqk
          tpsaPtcTemp(4)=-tpsaPtc(3)*si*sqk+ &
                         tpsaPtc(4)*co
          tpsaPtc(1) = tpsaPtcTemp(1)
          tpsaPtc(2) = tpsaPtcTemp(2)
          tpsaPtc(3) = tpsaPtcTemp(3)
          tpsaPtc(4) = tpsaPtcTemp(4)
        endif

       end subroutine quadmap
