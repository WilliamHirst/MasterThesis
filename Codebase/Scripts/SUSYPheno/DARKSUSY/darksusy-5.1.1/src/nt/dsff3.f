      real*8 function dsff3(x)
      implicit none
      real*8 dsff2,x
      dsff3=dsff2(exp(x))*exp(x)*1.0d15
      return
      end
