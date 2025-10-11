!****************************
!
!*** Copyright Notice ***
!
!A multi-language auto differentiation package (mAD) Copyright (c) 2025, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
!
!If you have questions about your rights to use or distribute this software,
!please contact Berkeley Lab's Intellectual Property Office at
!IPO@lbl.gov.
!
!NOTICE.  This Software was developed under funding from the U.S. Department
!of Energy and the U.S. Government consequently retains certain rights.  As
!such, the U.S. Government has been granted for itself and others acting on
!its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
!Software to reproduce, distribute copies to the public, prepare derivative
!works, and perform publicly and display publicly, and to permit others to do so.
!****************************
! Fortran module for forward auto-differentiation
! Ji Qiang - LBNL, jqiang@lbl.gov
!

Module TPSAad
  
  implicit none

  !dimmax: maximum # of variables to be differentiable
  !For better performance, this number should be the # differentiable variables.
  integer, parameter :: dimmax = 5
  integer :: TPS_Dim

  type dtpsad
    integer(kind=4), private  :: terms = dimmax+1
    double precision, dimension(dimmax+1) :: map=0.0d0
  end type dtpsad
  
  interface assignment(=)
    module procedure TpsaEq, TpsaEqN
  end interface

  interface assign_tpsad
    module procedure assignV1,assignV2
  end interface
  
  interface operator (+)
    module procedure TpsaAdd,TpsaAddN1,TpsaAddN2,TpsaAddN3,TpsaAddN4
  end interface
  
  interface operator (-)
    module procedure TpsaDec,TpsaDecN1,TpsaDecN2,TpsaDecN3,TpsaDecN4,TpsaDecN5
  end interface
  
  interface operator (*)
    module procedure TpsaMul,TpsaMulN1,TpsaMulN2,TpsaMulN3,TpsaMulN4
  end interface
  
  interface operator (/)
    module procedure TpsaDiv,TpsaDivN1,TpsaDivN2,TpsaDivN3,TpsaDivN4
  end interface
  
  interface exp
    module procedure TpsaExp
  end interface
  
  interface log
    module procedure TpsaLog
  end interface
  
  interface sqrt
    module procedure TpsaSqrt
  end interface
  
  interface pow
    module procedure TpsaPowd,TpsaPowi
  end interface
  
  interface sin
    module procedure TpsaSin
  end interface
  
  interface cos
    module procedure TpsaCos
  end interface
  
  interface tan
    module procedure TpsaTan
  end interface

  interface asin
    module procedure TpsaAsin
  end interface
  
  interface acos
    module procedure TpsaAcos
  end interface

  interface atan
    module procedure TpsaAtan
  end interface

  interface sinh
    module procedure TpsaSinh
  end interface
  
  interface cosh
    module procedure TpsaCosh
  end interface
  
  interface tanh
    module procedure TpsaTanh
  end interface
  
  interface asinh
    module procedure TpsaAsinh
  end interface

  interface acosh
    module procedure TpsaAcosh
  end interface

  interface atanh
    module procedure TpsaAtanh
  end interface

contains

    subroutine dtpsad_Initialize(nvar)
    
        implicit none
        integer, intent(in) :: nvar
        
        TPS_Dim   = nvar
    end subroutine


    subroutine init_tpsad(this,nvar)
        implicit none
        type(dtpsad), intent(out) :: this
        integer :: nvar

        this%terms = nvar+1
        this%map(1:this%terms) = 0.0d0
    end subroutine
    
    subroutine assignV1(this,a)
        implicit none
        type(dtpsad) :: this
        double precision :: a
        
        this%terms= TPS_Dim + 1
        this%map(1) = a
        this%map(2:this%terms) = 0.0d0

    end subroutine assignV1
    
    subroutine assignV2(this,a,i_var)
        implicit none
        type(dtpsad) :: this
        integer          :: i_var,i
        double precision :: a
        
          this%terms=TPS_Dim+1;
          this%map(1) = a
          this%map(2:this%terms)  = 0.0d0
          this%map(i_var+1) = 1.0d0
    end subroutine assignV2
    
!    Type(dtpsad) Function TpsaAdd(M,N)
!      implicit none
!      type(dtpsad) ,Intent(In) :: M,N
!      integer :: i
!      
!      TpsaAdd%terms = M%terms
!      i = M%terms
!      TpsaAdd%map(1:i) = M%map(1:i) + N%map(1:i)
!
!    End Function TpsaAdd

    Type(dtpsad) Function TpsaAdd(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M,N

      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms) +  N%map(1:M%terms)
    End Function TpsaAdd
    
    Type(dtpsad) Function TpsaAddN1(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      
      this%terms = M%terms
      this%map(1) = M%map(1) + N
      this%map(2:M%terms) = M%map(2:M%terms)
    End Function TpsaAddN1
    
    Type(dtpsad) Function TpsaAddN2(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N

      this%terms = M%terms
      this%map(1) = M%map(1) + N
      this%map(2:M%terms) = M%map(2:M%terms)
      
    End Function TpsaAddN2

    Type(dtpsad) Function TpsaAddN3(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1) = M%map(1) + N
      this%map(2:M%terms) = M%map(2:M%terms)
    End Function TpsaAddN3

    Type(dtpsad) Function TpsaAddN4(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1) = M%map(1) + N
      this%map(2:M%terms) = M%map(2:M%terms)

    End Function TpsaAddN4

    Type(dtpsad) Function TpsaDec(M,N)
      implicit none
      type(dtpsad) ,Intent(In) :: M,N
      integer :: i
      
        TpsaDec%terms = M%terms
        i = M%terms
        TpsaDec%map(1:i) = M%map(1:i) - N%map(1:i)
    End Function TpsaDec
    
    Type(dtpsad) Function TpsaDecN1(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      
      this%terms = M%terms
      this%map(1) = M%map(1) - N
      this%map(2:M%terms) = M%map(2:M%terms) 

    End Function TpsaDecN1
    
    Type(dtpsad) Function TpsaDecN2(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      
      this%terms = M%terms
      this%map(1) = N- M%map(1) 
      this%map(2:M%terms) = -M%map(2:M%terms) 
    End Function TpsaDecN2

    Type(dtpsad) Function TpsaDecN3(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1) = M%map(1) - N
      this%map(2:M%terms) = M%map(2:M%terms) 

    End Function TpsaDecN3

    Type(dtpsad) Function TpsaDecN4(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1) = N- M%map(1) 
      this%map(2:M%terms) = -M%map(2:M%terms) 
    End Function TpsaDecN4

    Type(dtpsad) Function TpsaDecN5(M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M

      this%terms = M%terms
      this%map(1:M%terms) = -M%map(1:M%terms) 
    End Function TpsaDecN5


    Type(dtpsad) Function TpsaMul(M,N) result(this)
      implicit none
      type(dtpsad) ,Intent(In) :: M,N
      integer :: i
      
      this%terms = max(M%terms,N%terms)
      this%map(1) = M%map(1)*N%map(1)
      this%map(2:M%terms) = M%map(1)*N%map(2:N%terms)+&
                               M%map(2:M%terms)*N%map(1)
    End Function TpsaMul
    
    Type(dtpsad) Function TpsaMulN1(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      
      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms)*N
    End Function TpsaMulN1
    
    Type(dtpsad) Function TpsaMulN2(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      
      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms)*N
    End Function TpsaMulN2

    Type(dtpsad) Function TpsaMulN3(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms)*N
    End Function TpsaMulN3

    Type(dtpsad) Function TpsaMulN4(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms)*N
    End Function TpsaMulN4

    Type(dtpsad) Function TpsaDiv(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M,N
      integer ::i
      
      if(abs(N%map(1)) <1.0d-14) then
        write(*,*) "Error: Divide by zero0, in CTPS"
        stop
        return
      endif
      
      this%terms = M%terms
!      i = this%terms
      i = M%terms
      this%map(1) = M%map(1) / N%map(1)
      this%map(2:i) = (M%map(2:i)*N%map(1)-N%map(2:i)*M%map(1))/(N%map(1)**2)
      
    End Function TpsaDiv
    
    Type(dtpsad) Function TpsaDivN1(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      
      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms)/N

    End Function TpsaDivN1

    Type(dtpsad) Function TpsaDivN3(M,N) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N

      this%terms = M%terms
      this%map(1:M%terms) = M%map(1:M%terms)/N

    End Function TpsaDivN3
    
    Type(dtpsad) Function TpsaDivN2(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      double precision, Intent(in) :: N
      integer :: i
      
      if(abs(M%map(1)) <1.0d-14) then
        write(*,*) "Error: Divide by zero2, in CTPS"
        return
      endif

      this%terms = M%terms
      i = M%terms
      this%map(1) = N/M%map(1)
      this%map(2:i) = -M%map(2:i)*N/(M%map(1)**2)

    End Function TpsaDivN2

    Type(dtpsad) Function TpsaDivN4(N,M) result(this)
      implicit none
      type(dtpsad)    , Intent(in) :: M
      integer, Intent(in) :: N
      integer :: i

      if(abs(M%map(1)) <1.0d-14) then
        write(*,*) "Error: Divide by zero4, in CTPS"
        return
      endif

      this%terms = M%terms
      i = M%terms
      this%map(1) = N/M%map(1)
      this%map(2:i) = -M%map(2:i)*N/(M%map(1)**2)

    End Function TpsaDivN4
    
    Subroutine TpsaEq(M,N)
      implicit none
      type(dtpsad) ,Intent(out) :: M
      type(dtpsad) ,Intent(In)  :: N
      
      M%terms  = N%terms
      M%map(1:M%terms) = N%map(1:N%terms)
    End Subroutine TpsaEq
    
    Subroutine TpsaEqN(M,N)
      implicit none
      type(dtpsad) ,Intent(out) :: M
      double precision ,Intent(In)  :: N
      
      M%terms= TPS_Dim + 1
      M%map(1) = N
      M%map(2:M%terms) = 0.0d0

    End Subroutine TpsaEqN
    
    Logical Function TpsaIsEq(M,N) result(res)
      implicit none
      type(dtpsad) ,Intent(in) :: M,N
      integer :: i,si
      
      si = M%terms
      if(si .ne. size(N%map)) then 
        res = .false.
        return
      endif
      
      do i=1,si
        if(abs(M%map(i)) > 1.0d-30 .or. abs(N%map(i)) > 1.0d-30) then
          if(abs(M%map(i) / N%map(i)- 1.0d0) > 1.0d-15) then
            res = .false.
            print*, "BBB",M%map(i),N%map(i)
            return
          endif
        endif
      enddo
      
      res = .true.
    End Function TpsaIsEq

    Logical Function TpsaIsEq2(M,N) result(res)
      implicit none
      type(dtpsad) ,Intent(in) :: M,N
      integer :: i,si

      si = M%terms
      if(si .ne. size(N%map)) then
        res = .false.
        return
      endif

      do i=1,si
        if(abs(M%map(i)) > 1.0d-30 .or. abs(N%map(i)) > 1.0d-30) then
          if(abs(M%map(i) / N%map(i)- 1.0d0) > 1.0d-15) then
            res = .false.
            print*, "BBB",M%map(i),N%map(i)
            return
          endif
        endif
      enddo

      res = .true.
    End Function TpsaIsEq2
    
!    Logical Function TpsaNotEq(M,N) result(res)
!      type(dtpsad) ,Intent(in) :: M,N
!      
!      res = .not. M.eq.N
!    End Function TpsaNotEq
    
    Type(dtpsad) Function TpsaInv(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      integer      :: i
      
      this%terms = M%terms
      i = M%terms
      if (abs(M%map(1)) < 1d-14) then 
        write(*,*) "Error: Divide by zeroInv, in CTPS"
      endif

      this%map(1) = 1.0d0/M%map(1)
      this%map(2:i) = -M%map(2:i)/(M%map(1)**2)

    End Function TpsaInv
    
    Type(dtpsad) Function TpsaExp(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      integer      :: i

      this%terms = M%terms
      this%map(1) = exp(M%map(1))
      this%map(2:M%terms) = exp(M%map(1))*M%map(2:M%terms)

    End Function TpsaExp
    
    
    Type(dtpsad) Function TpsaLog(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M

      this%terms = M%terms
      if(abs(M%map(1)).lt.1.0d-15) then
        print*,"zero in log function:"
        stop
      endif

      this%map(1) = log(M%map(1))
      this%map(2:M%terms) = M%map(2:M%terms)/M%map(1)

    End Function TpsaLog
    
    Type(dtpsad) Function TpsaSqrt(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M

      this%terms = M%terms
      this%map(1) = sqrt(M%map(1))
      this%map(2:M%terms) = 0.5d0*M%map(2:M%terms)/sqrt(M%map(1))

    End Function TpsaSqrt
    
    Type(dtpsad) Function TpsaPowd(M,b) result(this)
      implicit none
      type(dtpsad),     Intent(in) :: M
      double precision, Intent(in) :: b
      
      this%terms = M%terms
      if(abs(b-1.0d0)<1d-14) then
        this = M
        return
      else if (abs(b)<1e-14) then
        call assign_tpsad(this,1.0d0)
        return
      else
        this%map(1) = (M%map(1))**b
        this%map(2:M%terms) = b*(M%map(1))**(b-1)*M%map(2:M%terms)
      endif
    End Function TpsaPowd

    Type(dtpsad) Function TpsaPowi(M,b) result(this)
      implicit none
      type(dtpsad),     Intent(in) :: M
      integer, Intent(in) :: b

      this%terms = M%terms
      if(b == 1) then
        this = M
        return
      else if (b == 0) then
        call assign_tpsad(this,1.0d0)
        return
      else
        this%map(1) = (M%map(1))**b
        this%map(2:M%terms) = b*(M%map(1))**(b-1)*M%map(2:M%terms)
      endif
    End Function TpsaPowi
    
    Type(dtpsad) Function TpsaSin(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0, sin_a0, cos_a0

      this%terms = M%terms
      a0     = M%map(1)
      sin_a0 = dsin(a0)
      cos_a0 = dcos(a0)
      this%map(1) = sin_a0
      this%map(2:M%terms) = cos_a0*M%map(2:M%terms)

    End Function TpsaSin
    
    Type(dtpsad) Function TpsaCos(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0, sin_a0, cos_a0

      this%terms = M%terms
      a0     = M%map(1)
      sin_a0 = dsin(a0)
      cos_a0 = dcos(a0)
      this%map(1) = cos_a0
      this%map(2:M%terms) = -sin_a0*M%map(2:M%terms)

    End Function TpsaCos
    
    Type(dtpsad) Function TpsaTan(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M

      this%terms = M%terms
      this%map(1) = dtan(M%map(1))
      this%map(2:M%terms) = (1.0d0+dtan(M%map(1))**2)*M%map(2:M%terms)
    
    End Function TpsaTan

    Type(dtpsad) Function TpsaAsin(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0

      this%terms = M%terms
      a0 = M%map(1)
      if(abs(a0)>=1) then
        print*,"wrong argument in asin function!"
        stop
      endif
      this%map(1) = asin(a0)
      this%map(2:M%terms) = M%map(2:M%terms)/sqrt(1-a0**2)

    End Function TpsaAsin

    Type(dtpsad) Function TpsaAcos(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0

      this%terms = M%terms
      a0 = M%map(1)
      if(abs(a0)>=1) then
        print*,"wrong argument in acos function!"
        stop
      endif
      this%map(1) = acos(a0)
      this%map(2:M%terms) = -M%map(2:M%terms)/sqrt(1-a0**2)

    End Function TpsaAcos

    Type(dtpsad) Function TpsaAtan(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0

      this%terms = M%terms
      a0 = M%map(1)
      this%map(1) = atan(a0)
      this%map(2:M%terms) = M%map(2:M%terms)/(1+a0**2)

    End Function TpsaAtan
    
    Type(dtpsad) Function TpsaSinh(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M

      this%terms = M%terms
      this%map(1) = sinh(M%map(1))
      this%map(2:M%terms) = cosh(M%map(1))*M%map(2:M%terms)

    End Function TpsaSinh
    
    Type(dtpsad) Function TpsaCosh(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M

      this%terms = M%terms
      this%map(1) = cosh(M%map(1))
      this%map(2:M%terms) = sinh(M%map(1))*M%map(2:M%terms)

    End Function TpsaCosh

    Type(dtpsad) Function TpsaTanh(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M

      this%terms = M%terms
      this%map(1) = dtanh(M%map(1))
      this%map(2:M%terms) = (1.0d0-dtanh(M%map(1))**2)*M%map(2:M%terms)

    End Function TpsaTanh

    Type(dtpsad) Function TpsaAsinh(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0

      this%terms = M%terms
      a0 = M%map(1)
      this%map(1) = asinh(a0)
      this%map(2:M%terms) = M%map(2:M%terms)/sqrt(1+a0**2)

    End Function TpsaAsinh

    Type(dtpsad) Function TpsaAcosh(M) result(this)
      implicit none
      type(dtpsad), Intent(in) :: M
      double precision :: a0

      this%terms = M%terms
      a0 = M%map(1)
      if(abs(a0)<=1) then
        print*,"wrong argument in acosh function!"
        stop
      endif
      this%map(1) = acosh(a0)
      this%map(2:M%terms) = M%map(2:M%terms)/sqrt(a0**2-1)

    End Function TpsaAcosh

    Type(dtpsad) Function TpsaAtanh(M) result(this)
      implicit none
      !type(dtpsad), Intent(in) :: M
      type(dtpsad), Intent(in) :: M
      double precision :: a0

      this%terms = M%terms
      a0 = M%map(1)
      if(abs(a0)>=1) then
        print*,"wrong argument in atanh function!"
        stop
      endif
      this%map(1) = atanh(a0)
      this%map(2:M%terms) = M%map(2:M%terms)/(1-a0**2)

    End Function TpsaAtanh

End Module TPSAad
