***********************************************************************
*** The routine dsmhUED sets up masses for the mUED model; so far, only 
*** leptons and the LKP are implemented. The LKP is assumed to be the 
*** KK photon [roughly equal to the B^(1)] and assigned the same particle 
*** code as the lightest neutralino; correspondingly, KK leptons are 
*** assigned the same particle codes as sleptons.
*** 
***  input:  Rinv  -  inverse compactification radius [GeV]
***          LR     -  cutoff scale (given in units of Rinv)
***
*** author: torsten bringmann (troms@physto.se), 2009-05-10
***********************************************************************

      subroutine dsmhUED(Rinv,LR)
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'

      real*8 Rinv,LR

      real*8 mBW2(2,2),msl
      integer i
      real*8 pi
      parameter (pi=3.141592653589793238d0)
 
 
c      write(*,*) 'calculating mass splittings according to '//
c     &           'Cheng, Matchev & Schmaltz, PRD 66, 036005 (2002)...'

      call dssuconst_couplings  
      
c... set up B-W mass matrix and diagonalize it to get KK photon mass
      mBW2(1,1)=Rinv**2*(1.-39/2.*gyweak**2*zeta(3)/(16.*pi**4)
     &                 -1/6.*gyweak**2/(16.*pi**2)*log(LR**2))
     &                 +mass(kw)**2*(gyweak/g2weak)**2
      mBW2(2,2)=Rinv**2*(1.-5/2.*g2weak**2*zeta(3)/(16.*pi**4)
     &                 +15/2.*g2weak**2/(16.*pi**2)*log(LR**2))
     &                 +mass(kw)**2
      mBW2(1,2)=mass(kw)**2*(gyweak/g2weak)
      mBW2(1,2)=mBW2(2,1)
      mass(kn(1))=sqrt((mBW2(1,1)+mBW2(2,2)
     &            -sqrt((mBW2(1,1)-mBW2(2,2))**2+4*mBW2(1,2)))/2.)

c... set KK lepton doublet masses
      do 10 i=1,3
        mass(ksnu(i))=Rinv*(1.+(27*g2weak**2+9*gyweak**2)/16.
     &                         /(16.*pi**2)*log(LR**2))
        width(ksnu(i))=(gyweak/2.)**2*mass(kn(1))/(8.*pi)   ! (valid for dM<<1)
     &                 *((mass(ksnu(i))-mass(kn(1)))/mass(kn(1)))**2
 10   continue
      mass(kse(1))=mass(ksnu(1))+mass(ke)
      mass(ksmu(1))=mass(ksnu(1))+mass(kmu)
      mass(kstau(1))=mass(ksnu(1))+mass(ktau)
      width(kse(1))=(gyweak/2.)**2*mass(kn(1))/(8.*pi)      ! (valid for dM<<1)
     &              *((mass(kse(1))-mass(kn(1)))/mass(kn(1)))**2
      width(ksmu(1))=(gyweak/2.)**2*mass(kn(1))/(8.*pi)     ! (valid for dM<<1)
     &              *((mass(ksmu(1))-mass(kn(1)))/mass(kn(1)))**2
      width(kstau(1))=(gyweak/2.)**2*mass(kn(1))/(8.*pi)    ! (valid for dM<<1)
     &              *((mass(kstau(1))-mass(kn(1)))/mass(kn(1)))**2
 
c... set KK lepton singlet masses
      msl=Rinv*(1.+9/4.*gyweak**2/(16.*pi**2)*log(LR**2))
      mass(kse(2))=msl+mass(ke)
      mass(ksmu(2))=msl+mass(kmu)
      mass(kstau(2))=msl+mass(ktau)
      width(kse(2))=gyweak**2*mass(kn(1))/(8.*pi)           ! (valid for dM<<1)
     &                 *((mass(kse(2))-mass(kn(1)))/mass(kn(1)))**2
      width(ksmu(2))=gyweak**2*mass(kn(1))/(8.*pi)          ! (valid for dM<<1)
     &                 *((mass(ksmu(2))-mass(kn(1)))/mass(kn(1)))**2
      width(kstau(2))=gyweak**2*mass(kn(1))/(8.*pi)         ! (valid for dM<<1)
     &                 *((mass(kstau(2))-mass(kn(1)))/mass(kn(1)))**2
 

      return
      end
