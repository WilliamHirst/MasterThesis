      subroutine dsdddrde(t,e,n,a,z,stoich,rsi,rsd,modulation)
c_______________________________________________________________________
c  differential recoil rate
c  common:
c    'dssusy.h' - file with susy common blocks
c  input:
c    e : real*8              : nuclear recoil energy in keV
c    n : integer             : number of nuclear species
c    a : n-dim integer array : mass numbers
c    z : n-dim integer array : atomic numbers
c    stoich : n-dim integer array : stoichiometric coefficients
c    t : real*8 : time in days from 12:00UT Dec 31, 1999
c    modulation : integer : 0=no modulation 1=annual modulation
c      with no modulation (ie without earth velocity), t is irrelevant
c  output:
c    rsi : real*8 : spin-independent differential rate
c    rsd : real*8 : spin-dependent differential rate
c  units: counts/kg-day-keV
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c
c  2010/09/08 [C Savage]: fixed error in rate formula
c                         (missing factor of rhox)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'
      include 'dsddcom.h'
      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsnuclides.h'
      include 'dsmpconst.h'

      integer n,a(n),z(n),stoich(n)
      real*8 e,t,rsi,rsd,siff(10),sdff(10),cw(10),cwtmp
      integer kx,i,aa,zz,ii,dsnucldindx
      real*8 mni,muxi,norm,vmin,eta
      ! norm = (kg day keV)/(cm GeV^2 km/s) 
      !      = 86400/10 * [kg/(GeV/c^2)] * [c/(km/s)]^2
      parameter (norm=4.355983308d41)
      integer modulation,modulatio
      common /ddmodul/ modulatio

      cwtmp = 0.d0
      do i=1,n
         cw(i) = a(i)*stoich(i)
         cwtmp = cwtmp+cw(i)
      enddo
      do i=1,n
         cw(i) = cw(i)/cwtmp
      enddo
      
      kx = kn(kln)
      mx = mass(kx)
      call dsddsigmaff(e,n,a,z,siff,sdff)
      
      rsi = 0.d0
      rsd = 0.d0
      do i=1,n
         aa=a(i)
         zz=z(i)
         ii=dsnucldindx(aa,zz)
         if (ii.ne.0) then
            mni=nucldm(ii)
            muxi = mni*mx/(mni+mx)
            vmin = sqrt(mni*e*1.d-6/2.d0/muxi**2)*c_light
            call dsddeta(vmin,t,eta)
            rsi = rsi + cw(i)*siff(i)/2.d0/mx/muxi**2*eta
            rsd = rsd + cw(i)*sdff(i)/2.d0/mx/muxi**2*eta
         endif
      enddo
      rsi = norm*rhox*rsi
      rsd = norm*rhox*rsd

      return
      end
