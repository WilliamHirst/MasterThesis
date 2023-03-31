      function dsbsgwxy(x,y)

***********************************************************************
* Function W[x,y] in app. A p. 15 of                                  *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-08                   *
***********************************************************************

      implicit none
      real*8 x,y
      real*8 dsbsgwxy

c     This function vanish for x=y, but it is not 
c     obvious if the numerical calculation give that answer

      if(x.eq.y) then
       dsbsgwxy=0.d0
c       write(*,*) 
c     &  'dsbsgwxy called with argument 1 equal to argument 2'
       return
      endif      


      dsbsgwxy=(x+y-2.d0*x*y)/((x-1.d0)*(y-1.d0))
     &     +(x**3-2.d0*x*y+x**2*y)*log(x)
     &      /((x-1.d0)**2*(x-y))
     &     +(2.d0*x*y-x*y**2-y**3)*log(y)
     &      /((y-1.d0)**2*(x-y))


      return
      end


