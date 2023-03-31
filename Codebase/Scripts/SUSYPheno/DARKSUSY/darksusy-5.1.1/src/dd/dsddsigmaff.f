      subroutine dsddsigmaff(e,n,a,z,siff,sdff)
c_______________________________________________________________________
c  neutralino nucleus cross sections times form factors.
c  NOTE: the spin-dependent cross section is available only for a 
c        limited number of nuclei
c  input:
c    e : real*8              : nuclear recoil energy in keV
c    n : integer             : number of nuclear species
c    a : n-dim integer array : mass numbers
c    z : n-dim integer array : atomic numbers
c  output:
c    siff : n-dim real*8 array : spin-independent cross
c           section times form factor
c    sdff : n-dim real*8 array : spin-dependent cross
c           section times form factor
c    units: cm^2
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c  modified: pg 040605 added missing factor of 4 in spin-dependent
c  modified: pg 080217 changed to \lambda^2 J (J+1)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'
      include 'dsddcom.h'
      include 'dshacom.h'
      include 'dsnuclides.h'

      integer n,a(n),z(n)
      real*8 e,siff(n),sdff(n)
      integer kx,i,aa,zz,ii,dsnucldindx
      real*8 fermiGeV,gev2cm2,q,ff,
     &     mni,muxi,gps,gns,gpa,gna
      real*8 l2jjpp,l2jjnn,l2jjpn
      parameter (fermiGeV = 1.d0/0.1973269602d0)
      parameter (gev2cm2 = (197.327053d-16)**2)
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      
      call dsddgpgn(gps,gns,gpa,gna)

      kx = kn(kln)
      mx = mass(kx)

      do i=1,n
         aa=a(i)
         zz=z(i)
         ii=dsnucldindx(aa,zz)
         if (ii.ne.0) then
           mni = nucldm(ii)
           muxi = mni*mx/(mni+mx)
           q = dsqrt(2.d0*mni*e*1.d-6) ! Bug fixed 060415, 1.d0->2.d0


c  SPIN-INDEPENDENT

           call dsddffsi(q,aa,zz,ff)
           siff(i) = gev2cm2*muxi**2/pi*ff*(zz*gps+(aa-zz)*gns)**2
           
c  SPIN-DEPENDENT
c

           call dsddffsd(q,aa,zz,l2jjpp,l2jjnn,l2jjpn)
           sdff(i) = gev2cm2*muxi**2*4.d0/pi*
     &          (gpa**2*l2jjpp+gna**2*l2jjnn+gpa*gna*l2jjpn)
         
        endif
      enddo

      return
      end
