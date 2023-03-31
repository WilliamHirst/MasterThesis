      real*8 function dsansumaa()
c_______________________________________________________________________
c  sum the amplitude matrix
c  author: joakim edsjo (edsjo@physto.se) 96-02-02
c          paolo gondolo 99-1-15 factor of 3 faster
c  called by: dwdcos
c=======================================================================
      implicit none
      include 'dsandiacom.h'
      real*8 sumf,aaaa(108)
      integer i
      equivalence (aa,aaaa)

      sumf=0.0d0
      do i=1,108
         sumf = sumf + aaaa(i)**2
      enddo
      dsansumaa = sumf

      end
