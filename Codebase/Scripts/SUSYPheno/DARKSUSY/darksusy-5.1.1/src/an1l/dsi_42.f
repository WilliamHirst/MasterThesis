c====================================================================
c
c   auxiliary function used in:  
c   dsanzgpar.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsi_42(a,b,d,c)
      implicit none
      real*8 a,b,d,c,dssubka,
     &    tdelta1,hdelta1,delta2,hdelta2,delta3,r1,r2,delta,
     &    res1,res2,sum1,sum2,sum3,sum4
      sum1=0.d0
c term in log(b)
      sum2=0.d0
c term in log(d)
      sum3=0.d0
c term in log(c)
      sum4=0.d0
      res1=0.d0
      res2=0.d0
      tdelta1=(a-c/2.d0+b-1.d0)**2+4.d0*(a-c/2.d0)
      hdelta1=(a-c/2.d0+b-d)**2+4.d0*(a-c/2.d0)*d
      delta2=(b-a-1.d0)**2-4.d0*a
      hdelta2=(b-a-d)**2-4.d0*a*d
      delta3=(1.d0+d-c)**2-4.d0*d
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
c third two terms
      r1=d+a-b-c/2.d0
      r2=d-a-b+c/2.d0
      delta=hdelta1
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(r1,delta)-dssubka(r2,delta)
      else
        call dssubkb(r1,r2,delta,res1,res2)
        sum1=sum1+res1
        sum2=sum2+res2
        sum3=sum3-res2
      endif
c fourth two terms
      r1=d-a-b
      r2=d+a-b
      delta=hdelta2
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(r1,delta)-dssubka(r2,delta)
      else
        call dssubkb(r1,r2,delta,res1,res2)
        sum1=sum1+res1
        sum2=sum2+res2
        sum3=sum3-res2
      endif
c fifth term
      r1=c+1.d0-d
      delta=delta3
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(r1,delta)
      else
        call dssubkc(r1,delta,res1,res2)
        sum1=sum1+res1
        sum4=sum4+res2
      endif
c sixth term
      r1=c-1.d0+d
      delta=delta3
      if(delta.lt.0.d0) then
        sum1=sum1+dssubka(r1,delta)
      else
        call dssubkc(r1,delta,res1,res2)
        sum1=sum1+res1
        sum3=sum3+res2
        sum4=sum4+res2
      endif
c term in log(b)
      sum2=sum2+(1.d0+d-2.d0*b-c/2.d0)
c term in log(d)
      sum3=sum3-(1.d0+d-2.d0*b-c/2.d0)/2.d0
      dsi_42=(sum1+sum2*log(b)+sum3*log(d)
     &      +sum4*log(c))/(2.d0*(a-c/4.d0))
      return
      end


