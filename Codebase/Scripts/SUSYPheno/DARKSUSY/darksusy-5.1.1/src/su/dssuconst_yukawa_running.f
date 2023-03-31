      subroutine dssuconst_yukawa_running
c_______________________________________________________________________
c  Calculate Yukawas with running masses
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa fixed (pg)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      real*8 s12,s23,s13,c12,c13,c23,d13,aux,mscale
      complex*16 ed13
      real*8 dsrmq,dsralph3
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      integer i

      call dssuconst_yukawa ! non-running default yukawas

      mscale=2.d0*mass(lsp)
      alph3=dsralph3(mscale)
      g3stro=sqrt(4.d0*pi*alph3)
c      write(*,*) 'alpha3=',alph3,' alph3mz=',alph3mz
      aux = g2weak/dsqrt(2.d0)/mass(kw)
      yukawa(ktau)= aux*dsrmq(mscale,ktau)/cosbe
      yukawa(kqu(2))= aux*dsrmq(mscale,kc)/sinbe
      yukawa(kqu(3))= aux*dsrmq(mscale,kt)/sinbe
      yukawa(kqd(3))= aux*dsrmq(mscale,kb)/cosbe

      end
