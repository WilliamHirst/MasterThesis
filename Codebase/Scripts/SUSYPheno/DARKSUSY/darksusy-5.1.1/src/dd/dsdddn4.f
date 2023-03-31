      function dsdddn4(ms,mq,mx)
      implicit none
c
c     auxiliary function replacing the propagator for heavy squarks in
c     the drees-nojiri treatment of neutralino-nucleon scattering
c     twist-2 a^2+b^2 terms
c     dsdddn4 = I5 + 2 mx^2 I4
c
      real*8 dsdddn4,ms,mq,mx
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
      dsdddn4=(3.d0/2.d0/mx4*(log(ms2/mq2)-(ms2-mq2-mx2)*ll)+
     &     (4.d0+(mq2+ms2)/mx2-2.d0/ms2/mx2*(mx2-mq2)**2)/delta +
     &     6.d0*mq2/delta**2*(mx2+ms2-mq2) +
     &     ll/delta*(-2.d0*(ms2-mx2)-mq2-3.d0*mq2/mx2*(ms2-mq2)) +
     &     ll/delta**2*6.d0*mq2*(mq2*(mq2-ms2)-mx2*(2.d0*mq2+ms2-mx2))
     &     -2.d0/mx2/ms2)/scale**2
      return
      end
