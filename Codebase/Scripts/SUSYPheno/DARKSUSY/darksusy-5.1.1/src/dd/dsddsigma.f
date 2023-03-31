      subroutine dsddsigma(n,a,z,si,sd)
c_______________________________________________________________________
c  neutralino nucleus cross section.
c  input:
c    n : number of nuclear species
c    a : n-dim integer array with mass numbers
c    z : n-dim integer array with atomic numbers
c  output:
c    si : n-dim real*8 array with spin-independent cross sections at q=0
c    sd : n-dim real*8 array with spin-dependent cross sections at q=0
c    units: cm^2
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c  modified paolo gondolo 2008-02-18 add spin-dependent part
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'
      include 'dsddcom.h'
      include 'dshacom.h'
      include 'dsnuclides.h'

      integer n,a(n),z(n)
      real*8 si(n),sd(n)
      integer kx,i,ii,dsnucldindx
      real*8 gev2cm2,mni,fkin,gps,gns,gpa,gna,muxi
      real*8 l2jjpp,l2jjnn,l2jjpn
      parameter (gev2cm2 = (197.327053d-16)**2)
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      
      call dsddgpgn(gps,gns,gpa,gna)

      kx = kn(kln)
      mx = mass(kx)

      do i=1,n
         ii=dsnucldindx(a(i),z(i))
         if (ii.ne.0) then
            mni = nucldm(ii)
            muxi = mni*mx/(mni+mx)
            fkin = gev2cm2*muxi**2/pi
            si(i) = fkin*(z(i)*gps+(a(i)-z(i))*gns)**2
            call dsddffsd(0.d0,a(i),z(i),l2jjpp,l2jjnn,l2jjpn)
            sd(i) = fkin*4.d0*
     &          (gpa**2*l2jjpp+gna**2*l2jjnn+gpa*gna*l2jjpn)
         endif
      enddo

      return
      end
