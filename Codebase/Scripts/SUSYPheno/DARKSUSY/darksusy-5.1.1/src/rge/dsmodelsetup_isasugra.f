      subroutine dsmodelsetup_isasugra
c=======================================================================
c  replacement for dsmodelsetup for using ISASUGRA
c  author: E.A.Baltz, 2001 eabaltz@alum.mit.edu
c=======================================================================
      implicit none
      integer i
      real*8 aux,mscale,dsralph3,dsrmq

      include 'dsmssm.h'
      include 'dsprep.h'

      mass(0)=1.d10

c------------------------------------------------------ global constants
c      call dssuconst    ! this call is done already in dsisasugra_darksusy

      mx=mass(kn(1))
      
c------------- reset running things

      call dssuconst_yukawa_running

c      mscale=2.d0*mass(lsp)
c      alph3=dsralph3(mscale)
c      g3stro=sqrt(4.d0*pi*alph3)
c      aux = g2weak/dsqrt(2.d0)/mass(kw)
c      yukawa(ktau)= aux*dsrmq(mscale,ktau)/cosbe
c      yukawa(kqu(2))= aux*dsrmq(mscale,kc)/sinbe
c      yukawa(kqu(3))= aux*dsrmq(mscale,kt)/sinbe
c      yukawa(kqd(3))= aux*dsrmq(mscale,kb)/cosbe
  
c-------------------------------------------------- some useful vertices

      call dsvertx

c--------------------------------------------------- and particle widths

c...Add Higgs widhts and QCD correction to them and to vertices
      call dshigwid

cc...Widths changed 020613 by Joakim Edsjo
c      width(kcha(1)) = 1.0d0
c      width(kcha(2)) = 1.0d0
c      width(kn(1)) = 1.0d0
c      width(kn(2)) = 1.0d0
c      width(kn(3)) = 1.0d0
c      width(kn(4)) = 1.0d0
c
c      do i=1,6
c        width(ksl(i))=0.2d0
c        width(ksnu(i))=0.2d0
c        width(ksqu(i))=0.2d0
c        width(ksqd(i))=0.2d0
c      enddo

c...Widths changed 050907 by Joakim Edsjo
c      width(kcha(1)) = mass(kcha(1))*0.01d0
c      width(kcha(2)) = mass(kcha(2))*0.01d0
c      width(kn(1)) = mass(kn(1))*0.01d0
c      width(kn(2)) = mass(kn(2))*0.01d0
c      width(kn(3)) = mass(kn(3))*0.01d0
c      width(kn(4)) = mass(kn(4))*0.01d0

c      do i=1,6
c        width(ksl(i))=mass(ksl(i))*0.005d0
c        width(ksnu(i))=mass(ksnu(i))*0.005d0
c        width(ksqu(i))=mass(ksqu(i))*0.005d0
c        width(ksqd(i))=mass(ksqd(i))*0.005d0
c      enddo

      call dsspwid

      end
