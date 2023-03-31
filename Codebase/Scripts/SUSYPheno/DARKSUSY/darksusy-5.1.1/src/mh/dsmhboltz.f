***********************************************************************
*** dsmhboltz implements the Boltzmann equation for the evolution of
*** the WIMP temperature as dy/dx = rhs, expressed in the variables
***
***    y = m0*T_DM / (g**0.5 * T**2)
***    x = m0 / T
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      subroutine dsmhboltz(x,y,rhs)
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'
      include 'dsidtag.h'

      real*8 x,y,rhs

      integer n, lspecies, nsub, i
      real*8  a,m0,c,cn,nj,sqrtg,sqrtgt,subint(20),ctmp(20)
      real*8 dsmhfac
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c... these parameters are needed by dqagse
      real*8 aint,bint,epsabs,result,abserr
      integer neval,ier, limit
      parameter (limit=40)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last

      real*8 dsmhcint
      external dsmhcint

      rhs=0d0
      m0=mass(kn(1))

c... integrate collision term from T/10 to 30 T, which is a *very* good
c... approximation to the full integration [0,infinity],  taking account
c... of resonances by dividing the interval into suitable subintervals

      nsub=1
      Tint=m0/x
      aint=Tint/10.
      bint=30*Tint
c      write(*,*)
c      write(*,*) 'Total range : ',aint,bint
c      write(*,*) 'Resonances : ',nres,resk
c      write(*,*) '---------------'
      subint(1)=aint
      i=0        ! i indexes the resonance closest to start of present interval
 100  i=i+1
c      write(*,*) 'nsub,i,subint(i): ',nsub,i,subint(i)
      if (i.gt.nres.or.(resk(i)/2d0).ge.bint) goto 200 
      if (subint(nsub).ge.2*resk(i)) goto 100
      nsub=nsub+1
      if (subint(nsub-1).le.(resk(i)/2d0)) then
        subint(nsub)=resk(i)/1.9999
        i=i-1    ! the intervall covering this resonance will 
                 ! only be closed in next step
        goto 100
      endif
      if (i.eq.nres.or.(resk(i+1)/2d0).gt.2*resk(i)) then
        subint(nsub)=2.0001*resk(i)
      else
        subint(nsub)=(resk(i)+resk(i+1))/2.
      endif
      if (nsub.ge.limit-1.or.subint(nsub).ge.bint) then
        nsub=nsub-1
        goto 200
      else
        goto 100
      endif
c... no more resonances of interest found -- define last subinterval
  200 subint(nsub+1)=bint

      c=0
      do 300 i=1,nsub
        aint=subint(i)
        bint=subint(i+1)
c        write(*,*) 'Now integrating between ',aint,bint
        epsabs=5d-2*mheps*(bint-aint)*abs(dsmhcint(bint)-dsmhcint(aint))
        epsabs=epsabs*1.d20 ! integrand is multiplied by 1d20 for better integration properties
        call dqagse(dsmhcint,aint,bint,epsabs,5d-2*mheps,limit,result,
     &         abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        result=result/1.d20 ! correct for rescaling factor in integrand
        if (ier.ne.0) write(*,*) 'dsmhboltz: problem in integration of ',
     &     'subinterval ',i, '(of ',nsub,') at T = ',1d3*m0/x,' MeV',
     &    ' (ier = ',ier,' )'
        if (ier.eq.0) c=c+result
  300 continue


c... this implements, finally, the Boltzmann equation
      call dsmhdof(1d3*m0/x,sqrtg,sqrtgt)
      a=2*c*sqrt(45/4d0/pi**3)*x**2*mpl/m0
      rhs = a*(1.-sqrtg*y/x)/(sqrtg*sqrtgt)

      return

      end





        

