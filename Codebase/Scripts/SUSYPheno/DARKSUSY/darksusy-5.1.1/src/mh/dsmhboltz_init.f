***********************************************************************
*** Prepares integration of Boltzmann equation and has to be called 
*** before dsmhboltz (called by dsmhtkd)
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      subroutine dsmhboltz_init
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'
      include 'dsidtag.h'

      real*8  ktmp,m0
      integer i,j, n, t


c... set up potential resonances and order them
      nres=0
      n=0
      if (mhtype.eq.1.or.mhtype.eq.2) then   ! SUSY or UED
        m0=mass(kn(1))
 10     n=n+1
        if (n.ge.10) goto 40 ! no more resonances
        if (n.eq.1) ktmp=(mass(ksnu(1))**2-m0**2)/(2.*m0)
        if (n.eq.2) ktmp=(mass(ksnu(2))**2-m0**2)/(2.*m0)
        if (n.eq.3) ktmp=(mass(ksnu(3))**2-m0**2)/(2.*m0)
        if (n.eq.4) then 
           ktmp=(mass(kse(1))**2-m0**2-mass(ke)**2)/(2.*m0)
           if (ktmp.lt.mass(ke)/10.) goto 10
        endif
        if (n.eq.5) then
           ktmp=(mass(kse(2))**2-m0**2-mass(ke)**2)/(2.*m0)
           if (ktmp.lt.mass(ke)/10.) goto 10
        endif
        if (n.eq.6) then
           ktmp=(mass(ksmu(1))**2-m0**2-mass(kmu)**2)/(2.*m0)
           if (ktmp.lt.mass(kmu)/10.) goto 10
        endif
        if (n.eq.7) then
           ktmp=(mass(ksmu(2))**2-m0**2-mass(kmu)**2)/(2.*m0)
           if (ktmp.lt.mass(kmu)/10.) goto 10
        endif
        if (n.eq.8) then
           ktmp=(mass(kstau(1))**2-m0**2-mass(ktau)**2)/(2.*m0)
           if (ktmp.lt.mass(ktau)/10.) goto 10
        endif
        if (n.eq.9) then
           ktmp=(mass(kstau(2))**2-m0**2-mass(ktau)**2)/(2.*m0)
           if (ktmp.lt.mass(ktau)/10.) goto 10
        endif

        if (ktmp.gt.m0/10.) goto 10    ! only keep resonances for relevant 
                                       ! (small) SM particle energies
        nres=nres+1
        j=1
 20     if (j.lt.nres) then
          if (resk(j).le.ktmp) then
            j=j+1
            goto 20
          endif
        endif 
        do 30 i=0,nres-j-1
          resk(nres-i)=resk(nres-i-1)
 30     continue
        resk(j)=ktmp
        goto 10 ! check next resonance
      endif

 40   return

      end





        

