      function dsdddn3(ms,mq,mx)
      implicit none
c
c     auxiliary function replacing the propagator for heavy squarks in
c     the drees-nojiri treatment of neutralino-nucleon scattering
c     twist-2 a^2-b^2 terms
c     dsdddn3 = - mq^2 mx^2 I3
c
      real*8 dsdddn3,ms,mq,mx
      real*8 delta,ll,ms2,mq2,mx2,mx4,epsilon,sqrtdelta,scale
      epsilon=1.d-8
      scale=ms*ms
      ms2=ms*ms/scale
      mq2=mq*mq/scale
      mx2=mx*mx/scale
      delta=4.d0*mq2*ms2-(mq2+ms2-mx2)**2
      mx4=mx2*mx2
      if (abs(delta).lt.epsilon*mx4) then
         delta=epsilon*mx4
         mx2=mq2+ms2-dsqrt(4.d0*mq2*ms2-delta)
      endif
      sqrtdelta=dsqrt(dabs(delta))
      if (delta.gt.0.d0) then
         ll=2.d0/sqrtdelta*atan(sqrtdelta/(mq2+ms2-mx2))
      else
         ll=1.d0/sqrtdelta*log((mq2+ms2-mx2+sqrtdelta)/
     &        (mq2+ms2-mx2-sqrtdelta))
      endif
      dsdddn3=-mq2*mx2*(3.d0*(mx2-mq2-ms2)/delta**2+
     &     ll/delta*(-1.d0+6.d0*mq2*ms2/delta))/scale
      return
      end
