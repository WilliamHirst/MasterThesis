      real*8 function dspbsigmavpbar(en)
c total inelastic cross section pbar + h
c tan and ng, j.phys.g 9 (1983) 227. formula 3.7
      implicit none
      include 'dsmpconst.h'
      real*8 en,tp,vp,sigma,cc
      parameter (cc=c_light/1.d5)
      tp=en-m_p
      vp=dsqrt(2*m_p*tp+tp**2)/en
      sigma=24.7d0*(1.d0+0.584d0*tp**(-0.115d0)
     &    +0.856d0*tp**(-0.566d0))
c sigma in mb
      dspbsigmavpbar=sigma*(vp*cc)
      return
      end
