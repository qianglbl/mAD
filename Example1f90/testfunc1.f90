!test AD module with a single function of 5 variables
       program testfunc
       use TPSAad
       implicit none
       integer, parameter :: noi = 5

       type (dtpsad) :: x1,x2,x3,x4,x5,func

       !initialize the AD module
       call dtpsad_Initialize(noi)

       !initialize each variable
       call assign_tpsad(x1,3.0d0,1)
       call assign_tpsad(x2,0.1d0,2)
       call assign_tpsad(x3,0.1d0,3)
       call assign_tpsad(x4,1.0d0,4)
       call assign_tpsad(x5,1.0d0,5)

       !test function
       func = (2*cos(x1/x2)+x1/x2+exp(x2))*x3+2*x4+sinh(x5)

       print*,"func value and its derivatives w.r.t. 5 variables:"
       print*,func%map(1:noi+1)

       end program
