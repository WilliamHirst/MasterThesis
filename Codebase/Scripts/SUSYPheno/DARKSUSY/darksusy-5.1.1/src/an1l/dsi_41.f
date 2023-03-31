c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsi_41(a,b,c)
      implicit none
      real*8 a,b,c,dssubka,tdelta1,delta2,sum1,sum2,
     &       r1,r2,delta,res1,res2,sum3,sum4
      sum1=0.d0
      sum2=0.d0
      res1=0.d0
      res2=0.d0
      tdelta1=(a-c/2.d0+b-1.d0)**2+4.d0*(a-c/2.d0)
      delta2=(b-a-1.d0)**2-4.d0*a
c first two terms
      r1=1.d0+a-b-c/2.d0
      r2=1.d0-a-b+c/2.d0
      delta=tdelta1
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(r1,delta)-dssubka(r2,delta)
      else
        call dssubkb(r1,r2,delta,res1,res2)
        sum1=sum1+res1
        sum2=sum2+res2
      endif
c second two terms
      r1=1.d0-a-b
      r2=1.d0+a-b
      delta=delta2
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(r1,delta)-dssubka(r2,delta)
      else
        call dssubkb(r1,r2,delta,res1,res2)
        sum1=sum1+res1
        sum2=sum2+res2
      endif
c term in c
      delta=c*(c-4.d0*b)
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(c,delta)
      else
        sum3=sqrt(delta)/2.d0*
     &    log(c/4.d0*(1.d0+sqrt(1.d0-4.d0*b/c))**2)
        sum4=-sqrt(delta)/2.d0
        sum1=sum1+sum3
        sum2=sum2+sum4
      endif
c term in log(b)
      sum2=sum2+(1.d0-b+c/4.d0)
      dsi_41=(sum1+sum2*log(b))/(a-c/4.d0)
      return
      end


