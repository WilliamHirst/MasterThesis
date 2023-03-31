**********************************************************************
*** dsgf2s2thw takes the Fermi constant (as measured from mu decay)
*** and calculates the sin^2(theta_W) at the MZ scale.
*** Formulas from PDG 2007, chapter 10.
*** Input: opt = 1, calculate the on-shell value (use generally)
***        opt = 2, calculate s_Z
***        opt = 3, calculate the MS-bar value at MZ (use for GUT relations,
***                 and other places that require this value)
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-07-01<
**********************************************************************

      real*8 function dsgf2s2thw(gf,alphmz,mz,mt,opt)
      implicit none

      real*8 aux,gf,alphmz,mz,mt
      real*8 dr0,rhot,A0,a,b,dr1
      real*8 alphlow,alph,cw2,s2w,c,ccorr
      integer opt
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      alphlow=1.d0/137.03599911d0


      rhot=3.d0*gf*mt**2/(8*sqrt(2.d0)*pi**2)
      A0=sqrt(pi*alphlow/(sqrt(2.d0)*gf))
      dr0=0.06654 ! PDG
c...bosonic and Higgs contributions, extracted from PDG by iteration
c      dr1=0.00273d0 ! (not OK)
      dr1=-0.0005d0
      a=(1.d0-dr0)/(rhot-(1-dr0-dr1))
      b=-(A0/mz)**2/(rhot-(1-dr0-dr1))
      cw2=-0.5d0*a+sqrt(0.25d0*a**2-b)
      s2w=1.d0-cw2

      if (opt.eq.2) then
c...Below is the s^2_MZ value. This should be close enough to what we want
        alph=alphmz-(1/127.918)+(1/128.91) ! convert from MS-bar to on-shell
        aux=pi*alph/(sqrt(2.d0)*gf*mz**2)
        s2w=0.5d0-sqrt(0.25d0-aux)
      
      elseif (opt.eq.3) then
         ccorr=0.0034d0 ! additional corrections
         c=1.d0+rhot/s2w*(1.d0-s2w) +ccorr ! s2w is on-shell here
         s2w=s2w*c
      endif
         
      dsgf2s2thw=s2w

      return
      end

