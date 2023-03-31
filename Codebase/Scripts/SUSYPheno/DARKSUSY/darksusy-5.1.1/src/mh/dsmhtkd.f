***********************************************************************
*** The function dsmhtkd returns the kinetic coupling temperature [in MeV]
***
***  input: m0  - DM mass (in GeV)
***         how - integer flag
***
***  possible values for 'how':
***  1 - full calculation 
***      [recommended: following Bringmann, New J. Phys. 11, 105027 (2009)]
***  2 - fast calculation 
***      [neglecting SM masses and resonances in the scattering amplitude]
***  3 - "old" order-of-magnitude estimate 
***      [use only for comparison!]
***
*** author: torsten bringmann (troms@physto.se), 2010-01-10
***********************************************************************

      real*8 function dsmhtkd(m0,how)
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'

      integer how,i
      real*8  m0

      integer n,lspecies,conv
      real*8  tmpres,tmpres2,a,cn,cnsum(4),nj,sqrtgeff,sqrtg,sqrtgt
      real*8  dsmhfac,dsmhgamma
      real*8  xi,xf,yi,yf,rhs
      real*8  stepguess,stepmin
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c... these parameters are needed by dqagse
      real*8 aint,bint,epsabs,intres,abserr,dsmhcint2
      integer neval,ier, limit
      parameter (limit=20)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last

      external dsmhboltz,dsmhcint2


      tmpres=0d0
      dsmhtkd=0d0

      if (how.le.3) then !  use fast calculation to get a first estimate

        lspecies=4       ! start by assuming scattering with all leptons
        do 10 i=1,lspecies
          cnsum(i)=0d0
 10     continue
        do 50 i=1,6      ! add up contributing leptons
          call dsmhm2simp(m0,i,cn,n)
          if (n.le.-2) then
            write(*,*) 'badly defined parameter in dsmhm2simp: n=',n
            write(*,*) 'leaving dsmhtkd.f...'
            return
          endif
          if (i.le.3) cnsum(1)=cnsum(1)+cn ! neutrino contributions         
          if (i.ge.4) cnsum(i-2)=cnsum(i-3)+cn
  50    continue
        nj=(1.-1/2.**(n+3))*dsmhfac(n+4)
        if (n.le.4) nj=nj*zeta(n+4)


 100    if (lspecies.eq.1) sqrtgeff=1.83d0  ! only photons and neutrinos,
                                            ! reheating taken into account
        if (lspecies.eq.2) sqrtgeff=3.28d0  ! rel. e+,e-,nu and photons
        if (lspecies.eq.3) sqrtgeff=3.77d0  !  "   " and muons
        if (lspecies.eq.4) sqrtgeff=9.13d0  ! all leptons, geff at 1 GeV

        a = sqrt(10/2.**9/pi**9)/sqrtgeff 
     -      *cnsum(lspecies)*nj*mpl/m0
        tmpres = m0 / dsmhgamma((n+1d0)/(n+2d0))
     -             / ( (a/(n+2.))**(1/(n+2.)) )

c... check whether resulting Tkd is consistent with assumed number
c... of leptonic scattering partners
        if (lspecies.eq.4
     -       .and.tmpres.lt.(1.78)) then
           lspecies=3
           goto 100
        endif
        if (lspecies.eq.3
     -       .and.tmpres.lt.(1.2*105.6d-3)) then
           lspecies=2
           goto 100
        endif
        if (lspecies.eq.2
     -       .and.tmpres.lt.(1.2*0.511d-3)) then
           lspecies=1
           goto 100
        endif

c... Slightly improve quick estimate for g
  110   tmpres2=tmpres
        call dsmhdof(1d3*tmpres2,sqrtg,sqrtgt)
        tmpres=tmpres2*(1.+(sqrtg/sqrtgeff)**(1./(2.+n)))/2.
        sqrtgeff=sqrtg
        if (abs(tmpres-tmpres2)/tmpres.gt.mheps/2.) goto 110
      endif

c        if (how.eq.2.and.tmpres.gt.2*tqcdmin) write(*,*)
c     -    'WARNING: probably too high Tkd=',tmpres*1d3,
c     -    ' -- try full calculation'//
c     -    '(how=1 in dsmhtkd) instead!'


      if (how.eq.1) then   ! full calculation  
        call dsmhboltz_init    ! prepare integration of Boltzmann equation

c... to obtain an initial value for its integration, follow the Boltzmann
c... equation into the regime where small deviations from thermal eq. become possible
        xi=m0/tmpres/10.
  150   xf=1.5*xi
        call dsmhboltz(xf,0d0,rhs)
        if (rhs.gt.1/mheps/5.) then
          xi=xf
          goto 150
        endif
        call dsmhdof(1d3*m0/xi,sqrtg,sqrtgt)
        yi=xi/sqrtg

c... the integration of the Boltzmann Equation will proceed stepwise until
c... the required precision in Tkd is reached
  200   xf=1.2*xi
        yf=yi    ! starting with value yi, yf will contain
                 ! result of integration on return
        stepguess=xi*mheps/10.
        stepmin=xi*1D-20
        call dsmhinty(yf,xi,xf,mheps/10.,stepguess,stepmin,dsmhboltz,
     &                ier)
        if (ier.ne.0) then
           write(*,*) 'dsmhtkd: Cannot integrate Boltzmann Eq.!'
           write(*,*) 'Ti,Tistart,m0/yf,m0/yi : ',
     &                1d3*m0/xi,5d3*tmpres,m0/yf,m0/yi
c          if (xf.lt.1.05*xi) then
c            write (*,*) 'dsmhtkd: Cannot integrate Boltzmann Eq.!'
c            stop
c          endif
c          xf=xf/1.7
c          goto 200
        endif
        if (abs(yf-yi)/yf.ge.mheps/2.) then  ! not yet reached Tkd
          xi=xf
          yi=yf
          goto 200
        endif
        tmpres=m0/yf/sqrtg                   ! first estimate for Tkd

c... In case of large dg/dT, the above estimate for Tkd needs to be improved:
        conv=0
  300   tmpres2=tmpres
        call dsmhdof(1d3*tmpres2,sqrtg,sqrtgt)
        tmpres=(2*tmpres2+m0/yf/sqrtg)/3.
        conv=conv+1
        if (conv.gt.10) write(*,*) 'Slow convergence in dsmhTkd ',
     &                             '(how=2):Tkd = ',tmpres2,tmpres
        if (conv.gt.20) stop
        if (abs(tmpres-tmpres2)/tmpres.gt.mheps/2.) goto 300

      endif   ! how=1 (full calculation)


c... "old" estimate used in the literature (use only for comparison!)
      if (how.eq.3) then
        conv=0
  400   Tint=tmpres
        call dsmhdof(1d3*Tint,sqrtg,sqrtgt)  
        tmpres2=tmpres
        aint=1d-1
        bint=3d1
        epsabs=mheps*(bint-aint)*abs(dsmhcint2(bint)-dsmhcint2(aint))/2.
        call dqagse(dsmhcint2,aint,bint,epsabs,1d-3,20,intres,
     &         abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        if (ier.ne.0) write(*,*) 'DSmhTkd: integration problem for ',
     &     'how=3 at T = ',1d3*Tint,' MeV. Try different method!'
        tmpres=(2*tmpres2+m0*(1.089d3*m0/mpl*sqrtg/intres)**(1./(2.+n)))/3.
        conv=conv+1 
        if (conv.gt.10) write(*,*) 'Slow convergence in dsmhTkd ',
     &                             '(how=3):Tkd = ',tmpres2,tmpres
        if (conv.gt.20) then
           write(*,*) 'exiting dsmhTkd...'
           return
        endif
        if (abs(tmpres-tmpres2)/tmpres.gt.mheps) goto 400

      endif     ! how=3 (order-of-magnitude estimate from older literature)


      if (how.lt.1.or.how.gt.3) then
         write(*,*) 'ERROR in dsmhtkd -- unknown option: how=',how
         write(*,*) 'Stopping...'
         stop
      endif

      dsmhtkd=1d3*tmpres   ! give result in MeV

      return

      end





        

