      subroutine dsddfermi(q,a,z,ff)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Calculates F^2(q) 
c     
c     F(q) is Helm type form factor computed using
c     R1 from Lewin and Smith with s = 0.9 fm and with a and c 
c     directly from the Two-Parameter Fermi distribution in table VIII 
c     in Fricke et al. (instead of using the least squares fit from
c     Table IIIa in Lewin and Smith).
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
      include 'dsio.h'
      include 'dsmpconst.h'
      real*8 q,ff
      integer a,z
      real*8 f,s,ha,cpa,hc
      real*8 r1,x,x2,y,rsq,r1tmp

ccccc start testing for nuclei

      if (a.eq.1) then
         ff=1.d0
         return

      else if ((a.eq.12).and.(z.eq.6)) then
        hc = 2.0005d0
        ha = .5807d0

      else if ((a.eq.13).and.(z.eq.6)) then
        hc = 1.9958d0
        ha = .5807d0

      else if ((a.eq.14).and.(z.eq.6)) then
        hc = 2.0445d0
        ha = .5807d0

      else if ((a.eq.16).and.(z.eq.8)) then
        hc = 2.4130d0
        ha = .5807d0

      else if ((a.eq.18).and.(z.eq.8)) then
        hc = 2.5540d0
        ha = .5807d0

      else if ((a.eq.23).and.(z.eq.11)) then
        hc = 2.9393d0
        ha = .523d0

      else if ((a.eq.27).and.(z.eq.13)) then
        hc = 3.0554d0
        ha = .5807d0

      else if ((a.eq.28).and.(z.eq.14)) then
        hc = 3.1544d0
        ha = .5807d0

      else if ((a.eq.29).and.(z.eq.14)) then
        hc = 3.1482d0
        ha = .5807d0

      else if ((a.eq.30).and.(z.eq.14)) then
        hc = 3.1720d0
        ha = .5807d0

      else if ((a.eq.32).and.(z.eq.16)) then
        hc = 3.3816d0
        ha = .5807d0

      else if ((a.eq.34).and.(z.eq.16)) then
        hc = 3.4175d0
        ha = .5807d0
 
      else if ((a.eq.36).and.(z.eq.16)) then
        hc = 3.4411d0
        ha = .5807d0

      else if ((a.eq.36).and.(z.eq.18)) then
        hc = 3.5845d0
        ha = .5807d0

      else if ((a.eq.38).and.(z.eq.18)) then
        hc = 3.6025d0
        ha = .5807d0

      else if ((a.eq.40).and.(z.eq.18)) then
        hc = 3.6416d0
        ha = .5807d0

      else if ((a.eq.40).and.(z.eq.20)) then
        hc = 3.7221d0
        ha = .5807d0

      else if ((a.eq.42).and.(z.eq.20)) then
        hc = 3.7690d0
        ha = .5807d0

      else if ((a.eq.43).and.(z.eq.20)) then
        hc = 3.7477d0
        ha = .5807d0

      else if ((a.eq.44).and.(z.eq.20)) then
        hc = 3.7843d0
        ha = .5807d0

      else if ((a.eq.46).and.(z.eq.20)) then
        hc = 3.7537d0
        ha = .5807d0

      else if ((a.eq.48).and.(z.eq.20)) then
        hc = 3.7231d0
        ha = .5807d0

      else if ((a.eq.52).and.(z.eq.24)) then
        hc = 3.9742d0
        ha = .5807d0

      else if ((a.eq.56).and.(z.eq.26)) then
        hc = 4.1198d0
        ha = .5807d0

      else if ((a.eq.58).and.(z.eq.28)) then
        hc = 4.1772d0
        ha = .5807d0

      else if ((a.eq.64).and.(z.eq.30)) then
        hc = 4.3041d0
        ha = .5807d0

      else if ((a.eq.70).and.(z.eq.32)) then
        hc = 4.5687d0
        ha = .5807d0

      else if ((a.eq.72).and.(z.eq.32)) then
        hc = 4.5926d0
        ha = .5807d0

      else if ((a.eq.73).and.(z.eq.32)) then
        hc = 4.6015d0
        ha = .5807d0

      else if ((a.eq.74).and.(z.eq.32)) then
        hc = 4.6185d0
        ha = .5807d0

      else if ((a.eq.76).and.(z.eq.32)) then
        hc = 4.6294d0
        ha = .5807d0

      else if ((a.eq.72).and.(z.eq.32)) then
        hc = 4.5926d0
        ha = .5807d0

      else if ((a.eq.74).and.(z.eq.32)) then
        hc = 4.6185d0
        ha = .5807d0

      else if ((a.eq.76).and.(z.eq.32)) then
        hc = 4.6294d0
        ha = .5807d0

      else if ((a.eq.88).and.(z.eq.38)) then
        hc = 4.8399d0
        ha = .5807d0

      else if ((a.eq.90).and.(z.eq.40)) then
        hc = 4.9075d0
        ha = .5807d0

      else if ((a.eq.96).and.(z.eq.42)) then
        hc = 5.0711d0
        ha = .5807d0

      else if ((a.eq.98).and.(z.eq.42)) then
        hc = 5.1054d0
        ha = .5807d0

      else if ((a.eq.106).and.(z.eq.46)) then
        hc = 5.2847d0
        ha = .5807d0

      else if ((a.eq.108).and.(z.eq.46)) then
        hc = 5.3184d0
        ha = .5807d0

      else if ((a.eq.127).and.(z.eq.53)) then
        hc = 5.5931d0
        ha = .523d0

      else if ((a.eq.129).and.(z.eq.54)) then
        hc = 5.6315d0
        ha = .523d0

      else if ((a.eq.131).and.(z.eq.54)) then
        hc = 5.6384d0
        ha = .523d0

      else if ((a.eq.132).and.(z.eq.54)) then
        hc = 5.6460d0
        ha = .523d0

      else if ((a.eq.134).and.(z.eq.54)) then
        hc = 5.6539d0
        ha = .523d0

      else if ((a.eq.184).and.(z.eq.74)) then
        hc = 6.517d0
        ha = .53536d0

      else if ((a.eq.186).and.(z.eq.74)) then
        hc = 6.583d0
        ha = .48023d0

      else if ((a.eq.196).and.(z.eq.78)) then
        hc = 6.5597d0
        ha = .5807d0

      else if ((a.eq.204).and.(z.eq.82)) then
        hc = 6.6169d0
        ha = .5807d0

      else if ((a.eq.206).and.(z.eq.82)) then
        hc = 6.6311d0
        ha = .5807d0

      else if ((a.eq.207).and.(z.eq.82)) then
        hc = 6.6367d0
        ha = .5807d0

      else if ((a.eq.208).and.(z.eq.82)) then
        hc = 6.6468d0
        ha = .5807d0
 
      else
         if (prtlevel.gt.0)
     &        write (*,*) 'dsddfermi: Fermi parameters unavailable ',
     &        'for A=',a,' Z=',z
         ff = -1.d0
         return

      endif

ccccc lewin-smith calculations, hc formula removed

      if (a.eq.1) then
         ff = 1.d0
      else
         s = 0.9d0              ! fm
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
