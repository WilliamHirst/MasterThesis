      function dsbsgh2xy(x,y)

***********************************************************************
* Function H_2(x,y) in (12) of Degrassi et al.,                       *
* hep-ph/0009337                                                      *
* x and y must be  positive numbers                                   *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      real*8 x,y
      real*8 dsbsgh2xy

c     The expression is numerical unstable for x=y
c     for which we therefore introduce directly the limes expression 

      if(x.eq.y) then
       dsbsgh2xy=(1.d0-y+log(y))/(-1.d0+y)**2
       return
      endif

      dsbsgh2xy=x*log(x)/((1.d0-x)*(x-y))
     &         +y*log(y)/((1.d0-y)*(y-x))

      return
      end


