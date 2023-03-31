      function dsmsbaralph3(alph3_1,mu_1,mu_2)
c
c     computes the running msbar alpha_s using
c     formulas collected in Chetyrkin, Kuehn, and 
c     Steinhauser, hep-ph/0004189 (Thanks to Andre
c     Hoang for providing the reference)
c
c     Paolo Gondolo (paolo@physics.utah.edu) 2008-06-30
c
      implicit none
      include 'dsmssm.h'
      real*8 dsmsbaralph3,alph3_1,mu_1,mu_2
      real*8 beta0(3:6),b1(3:6),b2(3:6),b3(3:6),ell(3:6),muh(4:6),
     &     c0(3:6),c1(3:6),c2(3:6),c3(3:6)
      real*8 lh,a,ll
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      integer nf_1,nf_2,n,p,s
      logical first
      data first/.true./
      save first
c     coefficients in the beta and gamma functions
c     b1(N), c1(N), etc. corresponds to N flavors
c     only N=3,..,6 are needed
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
      save ell
c
c     ell is log(Lambda^2/mu_1^2)
c
      muh(4)=mcmc
      muh(5)=mbmb
      muh(6)=mtmt
c
      if (mu_1.ge.muh(6)) then
         nf_1=6
      else if (mu_1.ge.muh(5)) then
         nf_1=5
      else if (mu_1.ge.muh(4)) then
         nf_1=4
      else
         nf_1=3
      endif
      a=alph3_1/pi
      ell(nf_1)=-1.d0/beta0(nf_1)*(1.d0/a+b1(nf_1)*dlog(a)+
     &     (b2(nf_1)-b1(nf_1)**2)*a+
     &     (b3(nf_1)*0.5d0-b1(nf_1)*b2(nf_1)+b1(nf_1)**2*0.5d0)*a**2)-
     &     (b1(nf_1)/beta0(nf_1))*dlog(beta0(nf_1))
c
      if (mu_2.ge.muh(6)) then
         nf_2=6
      else if (mu_2.ge.muh(5)) then
         nf_2=5
      else if (mu_2.ge.muh(4)) then
         nf_2=4
      else
         nf_2=3
      endif
      if (nf_2.ne.nf_1) then
         if (nf_2.lt.nf_1) then
            s=-1
         else
            s=+1
         endif
         do p=nf_1+s,nf_2,s
            n=p-s
            lh=2.d0*dlog(muh(max(n,p))/mu_1)-ell(n)
            ell(p)=ell(n)+1.d0/beta0(p)*(
     &           (beta0(p)-beta0(n))*lh+
     &           (b1(p)-b1(n))*dlog(lh)-b1(p)*dlog(beta0(p)/beta0(n))+
     &           1.d0/(beta0(n)*lh)*(b1(n)*(b1(p)-b1(n))*dlog(lh)+
     &           b1(p)**2-b1(n)**2-b2(p)+b2(n)-s*11.d0/72.d0)+
     &           1.d0/(beta0(n)*lh)**2*(-b1(n)**2*0.5d0*(b1(p)-b1(n))*
     &           (dlog(lh))**2+b1(n)*(-b1(p)*(b1(p)-b1(n))+b2(p)-b2(n)
     &           -11.d0/72.d0)*dlog(lh)+0.5d0*(-b1(p)**3-b1(n)**3-b3(p)+
     &           b3(n))+b1(p)*(b1(n)**2+b2(p)-b2(n)+s*11.d0/72.d0)-
     &           s*0.9720566866267066d0+s*0.08465149176954733*min(p,n)))
         enddo
      endif
c
      ll=2.d0*dlog(mu_2/mu_1)-ell(nf_2)
      a=1.d0/(beta0(nf_2)*ll)
     &     -b1(nf_2)*dlog(ll)/(beta0(nf_2)*ll)**2
     &     +1.d0/(beta0(nf_2)*ll)**3*(b1(nf_2)**2*(
     &       (dlog(ll))**2-dlog(ll)-1.d0)+b2(nf_2))
     &     +1.d0/(beta0(nf_2)*ll)**4*(b1(nf_2)**3*
     &       (-(dlog(ll))**3+2.5d0*(dlog(ll))**2+
     &        2.d0*dlog(ll)-0.5d0)-
     &        3.d0*b1(nf_2)*b2(nf_2)*dlog(ll)+b3(nf_2)*0.5d0)
      dsmsbaralph3=a*pi
      return
      end
