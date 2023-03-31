

**********************************************************************
*** dsntdkfbigu is the velocity distribution in terms of big u=(u/v_e)^2.
c damour's velocity distribution
c l. bergstrom 99-03-12
c u = (u/v_e)^2
c lambda = j_z^{now}/j_z^0
**********************************************************************
      real*8 function dsntdkfbigu(bigu)
      implicit none
      include 'dsntdkcom.h'
      real*8 bigu,epsilon,term1,term2,umax
      umax=2.6154d0 ! cut-off value
      epsilon=0.18377*dklambda
c      write(*,*) 'epsilon=',epsilon,'  bigu=',bigu
      if (bigu.ge.umax) then
      dsntdkfbigu=0.d0
      goto 999
      endif
      if (bigu.ge.(1.d0-epsilon)) then
      term1=sqrt(max(bigu-1.0d0+epsilon,0.0d0))
      else
      term1=0.d0
      endif
      if (bigu.ge.(1.d0+epsilon)) then
      term2=sqrt(max(bigu-1.0d0-epsilon,0.0d0))
      else
      term2=0.d0
      endif
      dsntdkfbigu=1./epsilon/(3.d0-bigu)**0.6*(term1-term2)
 999  continue
c      if (bigu.gt.0.8d0) then
c        write(*,*) 'bigu=',bigu,'  umax=',umax,dsntdkfbigu
c      endif
      return
      end
