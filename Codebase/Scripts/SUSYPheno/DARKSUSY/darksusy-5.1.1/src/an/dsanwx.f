      real*8 function dsanwx(p)
c_______________________________________________________________________
c  Neutralino self-annihilation invariant rate.
c  Input:
c    p - initial cm momentum (real) for lsp annihilations
c  Output
c  BeginTex
c    \begin{displaymath}
c    W_{\rm{eff}} = \sum_{ij}\frac{p_{ij}}{p_{11}}
c    \frac{g_ig_j}{g_1^2} W_{ij} = 
c    \sum_{ij} \sqrt{\frac{[s-(m_{i}-m_{j})^2][s-(m_{i}+m_{j})^2]}
c    {s(s-4m_1^2)}} \frac{g_ig_j}{g_1^2} W_{ij}.
c    \end{displaymath}
c    where the $p$'s are the momenta, the $g$'s are the internal
c    degrees of freedom, the $m$'s are the masses and $W_{ij}$ is
c    the invariant annihilation rate for the included subprocess.
c  EndTex
c  uses dsabsq.
c  called by dsrdfunc, wirate, dsrdwintrp.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c  modified: joakim edsjo (edsjo@physto.se) 97-09-09
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      include 'dsandwcom.h'
      include 'dsprep.h'
      include 'dsidtag.h'
      real*8 p,dsanwxint,wtmp,dsandwdcosopt,dsandwdcos
      real*8 epsabs,epsrel,result,abserr,alist(5000),blist(5000),
     &  elist(5000),rlist(5000)
      real*8 al,be
      integer ier,iord(5000),last,limit,neval

      real*8 dsandwdcoss,eps,sum,y1,y2,dsandwdcosy
      external dsandwdcoss,dsandwdcosopt,dsandwdcos,dsandwdcosy

      real*8 dsandwdcosd
      external dsandwdcosd
c-----------------------------------------------------------------------
      real*8 pd
      common /gadint/ pd

      real*8 alph,bet
      common /yint/ alph,bet
      save /yint/
c-----------------------------------------------------------------------

      if (newmodelanwx) then
        newmodelanwx=.false.
        call dsanalbe(al,be)
        alph=al
        bet=be
      endif

      goto 90

      epsabs=1.d-2     !numerical accuracy
      epsrel=1.d-2
      limit=5000
      pd = p
      call dqagse(dsandwdcosd,-1.0d0,1.0d0,epsabs,epsrel,limit,result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      dsanwx=result/1.0d15
      goto 900


 90   eps=0.001d0
      sum=0.0d0
      pd = p

c      pd=16.064d0
c      call gewspecs(dsandwdcoss,-1.0d0,1.0d0,100,
c     & 'p1cth.dat                      ',33)
c      y1=(mco(1)**2+bet*pd**2)**(alph)
c      y2=(mco(1)**2)**(alph)
c      call gewspecs(dsandwdcosy,y1,y2,100,
c     & 'p1y.dat                        ',33)

c      pd=325.0d0
c      call gewspecs(dsandwdcoss,-1.0d0,1.0d0,100,
c     & 'p2cth.dat                      ',33)
c      y1=(mco(1)**2+bet*pd**2)**(alph)
c      y2=(mco(1)**2)**(alph)
c      call gewspecs(dsandwdcosy,y1,y2,100,
c     & 'p2y.dat                         ',33)

c      pd=1000.0d0
c      call gewspecs(dsandwdcoss,-1.0d0,1.0d0,100,
c     & 'p3cth.dat                      ',33)
c      y1=(mco(1)**2+bet*pd**2)**(alph)
c      y2=(mco(1)**2)**(alph)
c      call gewspecs(dsandwdcosy,y1,y2,100,
c     & 'p3y.dat                        ',33)
c      stop

      if (p.ge.2.5d0*mco(1)) then
c        call dsanalbe2(p,alph,be)
c        bet=be
        y1=(mco(1)**2+bet*pd**2)**(alph)
        y2=(mco(1)**2)**(alph)
        call dgadap(y1,y2,dsandwdcosy,eps,sum)
        dsanwx=sum/1.0d15
      else
        call dgadap(-1.0d0,1.0d0,dsandwdcoss,eps,sum)
        dsanwx=sum/1.0d15
      endif
c      write(*,*) rdtag,nr,p,dsanwx

c      call gadap(-0.99,0.99,dsandwdcoss,eps,sum)
c      dsanwx=dsanwx+dble(sum)/1.0d15
c      call gadap(0.99,1.0,dsandwdcoss,eps,sum)
c      dsanwx=dsanwx+dble(sum)/1.0d15
      goto 900

 100  if (p.lt.0.5d0*mco(1)) then
        dsanwx=dsanwxint(p,-1.0d0,1.0d0)
      else if ((p.ge.0.5d0*mco(1)).and.(p.lt.1.0d0*mco(1))) then
        dsanwx=dsanwxint(p,-1.0d0,-0.9d0)
        dsanwx=dsanwx+dsanwxint(p,-0.9d0,0.9d0)
        dsanwx=dsanwx+dsanwxint(p,0.9d0,1.0d0)
      else
        dsanwx=dsanwxint(p,-1.0d0,-0.999d0)
        dsanwx=dsanwx+dsanwxint(p,-0.999d0,-0.99d0)
        dsanwx=dsanwx+dsanwxint(p,-0.99d0,-0.97d0)
        dsanwx=dsanwx+dsanwxint(p,-0.97d0,-0.9d0)
        dsanwx=dsanwx+dsanwxint(p,-0.9d0,-0.5d0)
        dsanwx=dsanwx+dsanwxint(p,-0.5d0,0.5d0)
        dsanwx=dsanwx+dsanwxint(p,0.5d0,0.9d0)
        dsanwx=dsanwx+dsanwxint(p,0.9d0,0.97d0)
        dsanwx=dsanwx+dsanwxint(p,0.97d0,0.99d0)
        dsanwx=dsanwx+dsanwxint(p,0.99d0,0.999d0)
        dsanwx=dsanwx+dsanwxint(p,0.999d0,1.0d0)
      endif

 900  continue

      return

      end
