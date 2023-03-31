      subroutine dsddlsff(q,a,z,ff)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Calculates F^2(q) 
c     
c     F(q) is Helm type form factor computed using
c     R1 from Lewin and Smith with s = 0.9 fm and with a and c 
c     using the least squares fit from Table IIIa in Lewin and Smith.
c     
c     Parameters and Equations found in:
c     Lewin J D and Smith P F 1996 Astropart. Phys. 6 87-112.
c     Fricke G et al. 1995 Atomic Data and Nuclear Data Tables 60 177-285.
c     Duda, Kemper, and Gondolo JCAP 0704:012,2007. 
c     
c     
c     Authors: Gintaras Duda gkduda@creighton.edu (2007-06-27)
c              (w/ students George Reifenberger and Katherine Garrett)
c
c     Input: a: Mass number of Element
c            z: Atomic Number of Element
c            q: Momentum transfer in GeV
c
c     Output: ff: Form Factor squared (normalized to F^2(0) = 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-------------------------------------------------------

      implicit none
      include 'dsmpconst.h'
      real*8 q,ff
      integer a,z
      real*8 f,s,ha,cpa,hc
      real*8 r1,x,x2,y,rsq,r1tmp

      if (a.eq.1) then
         ff = 1.d0
      else
         s = 0.9d0              ! fm
         ha = 0.52d0            ! fm
         hc = 1.23d0*exp(log(dble(a))/3.d0)-0.6d0
         cpa = (7.d0/3.d0)*(pi*ha)**2
         rsq = hc*hc+cpa
         r1tmp = rsq-5.d0*s*s
         if (r1tmp.gt.0.d0) then
            r1 = sqrt(r1tmp)
         else
            r1 = sqrt(rsq)
         endif 
         x = dabs(q*r1*fermiGeV)
         y = q*s*fermiGeV
         if (x.gt.5.d-8) then
            f = 3.d0*(sin(x)-x*cos(x))/x**3
         else
            x2 = x*x
            f = 1.d0+x2*(-0.1d0+x2*(1.d0/280.d0+x2*(-1.d0/15120.d0+
     &           x2/1330560.d0)))
         endif
         ff = f**2*exp(-y**2)
      endif
      
      return
      end
