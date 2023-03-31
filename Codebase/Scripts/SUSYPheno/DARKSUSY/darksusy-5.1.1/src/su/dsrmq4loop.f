      function dsrmq4loop(mscale,k)
c     
c     run the mass of particle k from the 
c     msbar scale mz to the msbar scale mscale
c     using 4 loop formulas in  Chetyrkin, Kuehn, and 
c     Steinhauser, hep-ph/0004189 (Thanks to Andre
c     Hoang for providing the reference)
c     
c     Paolo Gondolo (paolo@physics.utah.edu) 2008-06-30
c     
      implicit none
      include 'dsmssm.h'
      integer k
      real*8 dsrmq4loop,mscale
      real*8 ans,x,cc,ccold,zeta,a,mu1,mu2,mk,dsralph34loop
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      integer n1,n2,n,s,p
      real*8 beta0(3:6),b1(3:6),b2(3:6),b3(3:6),ell(3:6),muh(4:6),
     &     c0(3:6),c1(3:6),c2(3:6),c3(3:6)
      data beta0/2.25d0,2.083333333333333d0,
     &     1.916666666666667d0,1.75d0/
      data b1/1.777777777777778d0,1.54d0,
     &     1.260869565217391d0,0.9285714285714286d0/
      data b2/4.471064814814815d0,3.047638888888889d0,
     &     1.474788647342995d0,0.2901785714285714d0/
      data b3/20.99023981042311d0,15.06597453710647d0,
     &     9.835916430959708d0,5.518490496829723d0/
      data c0/0.4444444444444444d0,0.48d0,
     &     0.5217391304347826d0,0.5714285714285714d0/
      data c1/1.685185185185185d0,1.753333333333333d0,
     &     1.833333333333333d0,1.928571428571429d0/
      data c2/5.520091095254772d0,4.774579325315020d0,
     &     3.871232980248869d0,2.764956467163064d0/
      data c3/19.67240376254819d0,13.10540814442490,
     &     5.757037169774640,2.576096786623935d0/
c
      mu2=mscale
c
c     no running except for quarks
c
      if (k.ne.ku.and.k.ne.kd.and.k.ne.kc.and.
     &     k.ne.ks.and.k.ne.kt.and.k.ne.kb) then
         dsrmq4loop=mk
         return
      endif
c
c     no running below mqmq for c,b,t
c     no running below 2 GeV for u,d,s
c
      if (k.eq.kt.and.mu2.lt.mtmt) then
         dsrmq4loop=mtmt
         return
      endif
      if (k.eq.kb.and.mu2.lt.mbmb) then
         dsrmq4loop=mbmb
         return
      endif
      if (k.eq.kc.and.mu2.lt.mcmc) then
         dsrmq4loop=mcmc
         return
      endif
      if (k.eq.ks.and.mu2.lt.2.d0) then
         dsrmq4loop=ms2gev
         return
      endif
      if (k.eq.kd.and.mu2.lt.2.d0) then
         dsrmq4loop=ms2gev
         return
      endif
      if (k.eq.ku.and.mu2.lt.2.d0) then
         dsrmq4loop=ms2gev
         return
      endif
c
c     start from mqmq for c,b,t
c     start from 2 GeV for u,d,s
c
      if (k.eq.kt) then
         mk=mtmt
         mu1=mtmt
         n1=6
      else if (k.eq.kb) then
         mk=mbmb
         mu1=mbmb
         n1=5
      else if (k.eq.kc) then
         mk=mcmc
         mu1=mcmc
         n1=4
      else if (k.eq.ks) then
         mk=ms2gev
         mu1=2.d0
         n1=3
      else if (k.eq.ku) then
         mk=mu2gev
         mu1=2.d0
         n1=3
      else if (k.eq.kd) then
         mk=md2gev
         mu1=2.d0
         n1=3
      endif
c
      muh(4)=mcmc
      muh(5)=mbmb
      muh(6)=mtmt
c
      if (mu2.ge.mtmt) then
         n2=6
      else if (mu2.ge.mbmb) then
         n2=5
      else if (mu2.ge.mcmc) then
         n2=4
      else
         n2=3
      endif
c
      x=dsralph34loop(mu1)/pi
      cc=x**c0(n1)*(1.d0+x*((c1(n1)-b1(n1)*c0(n1))+
     &     x*(0.5d0*((c1(n1)-b1(n1)*c0(n1))**2+c2(n1)-b1(n1)*c1(n1)+
     &        b1(n1)**2*c0(n1)-b2(n1)*c0(n1))+
     &     x*(1.d0/3.d0*(c3(n1)-b1(n1)*c2(n1)+b1(n1)**2*c1(n1)
     &        -b2(n1)*c1(n1)-b1(n1)**3*c0(n1)
     &        -2.d0*b1(n1)*b2(n1)*c0(n1)-b3(n1)*c0(n1))))))
      ans=mk/cc
c
      if (n1.ne.n2) then
         if (n2.gt.n1) then
            s=+1
         else
            s=-1
         endif
         do n=n1,n2-s,s
            p=n+s
            x=dsralph34loop(muh(max(n,p)))/pi
            ans=ans*(1.d0-s*x**2*(89.d0/432.d0)
     &           +x**3*(-s*1.847626742224786d0-s*
     &           min(n,n+s)*0.02472760936815077))
            ccold=x**c0(n)*(1.d0+x*((c1(n)-b1(n)*c0(n))+
     &           x*(0.5d0*((c1(n)-b1(n)*c0(n))**2+c2(n)-b1(n)*c1(n)+
     &           b1(n)**2*c0(n)-b2(n)*c0(n))+
     &           x*(1.d0/3.d0*(c3(n)-b1(n)*c2(n)+b1(n)**2*c1(n)
     &           -b2(n)*c1(n)-b1(n)**3*c0(n)
     &           -2.d0*b1(n)*b2(n)*c0(n)-b3(n)*c0(n))))))
            cc=x**c0(p)*(1.d0+x*((c1(p)-b1(p)*c0(p))+
     &           x*(0.5d0*((c1(p)-b1(p)*c0(p))**2+c2(p)-b1(p)*c1(p)+
     &           b1(p)**2*c0(p)-b2(p)*c0(p))+
     &           x*(1.d0/3.d0*(c3(p)-b1(p)*c2(p)+b1(p)**2*c1(p)
     &           -b2(p)*c1(p)-b1(p)**3*c0(p)
     &           -2.d0*b1(p)*b2(p)*c0(p)-b3(p)*c0(p))))))
c            downward: cdown/cup = cold/cnew
c            upward: cup/cdown = cold/cnew
            ans=ans*ccold/cc
         enddo
      endif
c
      x=dsralph34loop(mu2)/pi
      cc=x**c0(n2)*(1.d0+x*((c1(n2)-b1(n2)*c0(n2))+
     &     x*(0.5d0*((c1(n2)-b1(n2)*c0(n2))**2+c2(n2)-b1(n2)*c1(n2)+
     &        b1(n2)**2*c0(n2)-b2(n2)*c0(n2))+
     &     x*(1.d0/3.d0*(c3(n2)-b1(n2)*c2(n2)+b1(n2)**2*c1(n2)
     &        -b2(n2)*c1(n2)-b1(n2)**3*c0(n2)
     &        -2.d0*b1(n2)*b2(n2)*c0(n2)-b3(n2)*c0(n2))))))
      ans=ans*cc
c
      dsrmq4loop=ans
      return
      end
