

      function dsbsggxy(x,y)

***********************************************************************
* Function G'(x,y) in app. B p. 18 of                                 *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-22                   *
***********************************************************************

      implicit none
      real*8 x,y
      real*8 dsbsggxy


      dsbsggxy=x*(1.d0-(1.d0/(x-1.d0)+y/(x-y))*log(x))
     &         /((x-y)*(x-1.d0))

c...Now add the remaining part
c...when y=1, the expression is numerical unstable, so we add
c...limes of the expression directly


      if (abs(y-1.0d0).gt.1d-10) then
        dsbsggxy=dsbsggxy+
     &     y**2*log(y)/((x-y)**2*(y-1.d0))
      else
        dsbsggxy=dsbsggxy+1.d0/(x-1.d0)**2
      endif


      return
      end


