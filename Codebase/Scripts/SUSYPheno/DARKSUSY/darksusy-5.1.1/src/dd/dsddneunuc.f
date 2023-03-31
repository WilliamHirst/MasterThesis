      subroutine dsddneunuc(sigsip,sigsin,sigsdp,sigsdn)
c_______________________________________________________________________
c  neutralino nucleon cross section.
c  common:
c    'dssusy.h' - file with susy common blocks
c  output:
c    sigsip, sigsin : proton and neutron spin-independent cross sections
c    sigsdp, sigsdn : proton and neutron spin-dependent cross sections
c    units: cm^2
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995,2002
c     13-sep-94 pg no drees-nojiri twist-2 terms
c     22-apr-95 pg important bug corrected [ft -> ft mp/mq]
c     06-apr-02 pg drees-nojiri treatment added
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'
      include 'dsddcom.h'
      include 'dshacom.h'
      include 'dsmpconst.h'

      real*8 sigsip,sigsdp,sigsin,sigsdn
      integer kx
      real*8 fkinp,fkinn,gps,gns,gpa,gna

      kx = kn(kln)
      mx = mass(kx)
      fkinp = 4./pi*(m_p*mx/(m_p+mx))**2
      fkinn = 4./pi*(m_n*mx/(m_n+mx))**2

      call dsddgpgn(gps,gns,gpa,gna)

      sigsip = fkinp*(gps/2.d0)**2   ! in gev^-2
      sigsin = fkinn*(gns/2.d0)**2   ! in gev^-2

      sigsdp = fkinp*(3.d0/4.d0)*(gpa)**2 ! in gev^-2
      sigsdn = fkinn*(3.d0/4.d0)*(gna)**2 ! in gev^-2

      sigsip = gev2cm2*sigsip
      sigsin = gev2cm2*sigsin
      sigsdp = gev2cm2*sigsdp
      sigsdn = gev2cm2*sigsdn

      return
      end
