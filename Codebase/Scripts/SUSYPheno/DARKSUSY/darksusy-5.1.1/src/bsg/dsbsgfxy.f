
      function dsbsgfxy(x,y)

***********************************************************************
* Function F'(x,y) in app. B p. 18 of                                 *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-22                   *
***********************************************************************

      implicit none
      real*8 x,y
      real*8 dsbsgfxy


      dsbsgfxy=(x**2-y)*log(x)/((x-y)**2*(x-1.d0)**2)
     &     -y*log(y)/((x-y)**2*(y-1.d0))
     &     -1.d0/((x-y)*(x-1.d0))


      return
      end


