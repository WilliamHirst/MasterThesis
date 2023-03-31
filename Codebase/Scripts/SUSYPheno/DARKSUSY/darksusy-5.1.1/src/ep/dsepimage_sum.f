************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98, dec-04-02 (eb)
************************************************************************


************************************************************************
      real*8 function dsepimage_sum(deltav)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 deltav,nn,root,sum,term,charge,erf
      external erf
      if (deltav.le.0.0d0) then
         dsepimage_sum=0.0d0
         return
      endif
      sum=0.0d0
      charge=-1.0d0
      root=l_h/sqrt(4.0d0*k0tau*deltav)
      term=erf(root)
      sum=sum+term
      nn=1.0d0

 10   continue
c...The next line was the original image sum
c      term=erf((charge+2.0d0*nn)*root)+erf((charge-2.0d0*nn)*root)
c...Let's now use erfc instead of erf since this is numerically
c...more stable (only affects positrons that have lost a very large fraction
c...of their energy). Fix 2002-12-04 by Ted Baltz
      term=-erfc((charge+2.0d0*nn)*root)+erfc((-charge+2.0d0*nn)*root)
      sum=sum+term
      charge=-1.0d0*charge
      nn=nn+1.0d0
      if (abs(term).gt.sumtol*abs(sum)) goto 10
      dsepimage_sum=sum
      return
      end
