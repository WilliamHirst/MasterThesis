      real*8 function dsfff3(x)
      implicit none
      real*8 dsfff2,x
      dsfff3=dsfff2(exp(x))*exp(x)*1.0d15
      return
      end
