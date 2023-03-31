      real*8 function dslnff(x)
      implicit none
      real*8 dsff,x
      dslnff=dsff(exp(x))*exp(x)
      return
      end
