      function dsmqpole4loop(k,mqmq)
      implicit none
      include 'dsmssm.h'
      real*8 dsmqpole4loop,mqmq
      integer k
      real*8 dsralph34loop,dsrmq4loop,sum,x,a,fac,ans
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
c     no running except for quarks
c
      if (k.ne.ku.and.k.ne.kd.and.k.ne.kc.and.
     &     k.ne.ks.and.k.ne.kt.and.k.ne.kb) then
         dsmqpole4loop=mqmq
         return
      endif
c
      if (k.eq.kt) then
         ans=mass(kt)
      elseif (k.eq.kb) then
         a=dsralph34loop(mqmq)/pi
         sum=0.d0
         x=dsrmq4loop(mqmq,ku)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         x=dsrmq4loop(mqmq,kd)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         x=dsrmq4loop(mqmq,ks)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         x=dsrmq4loop(mqmq,kc)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         fac=1.d0+a*(4.d0/3.d0+
     &        a*(13.44339625693176d0-1.041366911171631d0*5.d0+
     &        4.d0/3.d0*sum+a*75.d0))
         ans=mqmq*fac
      elseif (k.eq.kc) then
         a=dsralph34loop(mqmq)/pi
         sum=0.d0
         x=dsrmq4loop(mqmq,ku)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         x=dsrmq4loop(mqmq,kd)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         x=dsrmq4loop(mqmq,ks)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         fac=1.d0+a*(4.d0/3.d0+
     &        a*(13.44339625693176d0-1.041366911171631d0*4.d0+
     &        4.d0/3.d0*sum+a*96.d0))
         ans=mqmq*fac
      elseif (k.eq.ks) then
         a=dsralph34loop(mqmq)/pi
         sum=0.d0
         x=dsrmq4loop(mqmq,ku)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         x=dsrmq4loop(mqmq,kd)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         fac=1.d0+a*(4.d0/3.d0+
     &        a*(13.44339625693176d0-1.041366911171631d0*3.d0+
     &        4.d0/3.d0*sum+a*119.d0))
         ans=mqmq*fac
      elseif (k.eq.kd) then
         a=dsralph34loop(mqmq)/pi
         sum=0.d0
         x=dsrmq4loop(mqmq,ku)/mqmq
         sum=sum+x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
         fac=1.d0+a*(4.d0/3.d0+
     &        a*(13.44339625693176d0-1.041366911171631d0*2.d0+
     &        4.d0/3.d0*sum+a*143.d0))
         ans=mqmq*fac
      elseif (k.eq.ku) then
         a=dsralph34loop(mqmq)/pi
         sum=0.d0
         fac=1.d0+a*(4.d0/3.d0+
     &        a*(13.44339625693176d0-1.041366911171631d0*1.d0+
     &        4.d0/3.d0*sum+a*168.d0))
         ans=mqmq*fac
      endif
      dsmqpole4loop=ans
      return
      end
