      subroutine  dsddsogff(q,a,z,ff) 
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calculates form factor from a charge density given as a sum of gaussians 
c
c     Parameters found in:
c     H.DeVries, C.W.DeJager and C.DeVries, At. Data and Nucl.
c         Data Tables.  36 (1987) 495.
c     See also: Duda, Kemper, and Gondolo JCAP 0704:012,2007. 
c    
c     Authors: Gintaras Duda gkduda@creighton.edu (2007-06-27)
c              (w/ students Ann Kemper and George Reifenberger)
c
c     Input: a: Mass number of Element
c            z: Atomic Number of Element
c            q: Momentum transfer in GeV
c
c     Output: ff: Form Factor squared (normalized to F^2(0) = 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'dsio.h'
      include 'dsmpconst.h'

      integer a,z
      real*8 q,ff
      real*8 capq(12),capr(12),x,gamma,f,rp,qr,b,
     &     qfermi,xg,sinxx
      integer i,p
      
      if (a.eq.1) then
         ff=1.d0
         return

      else if ((a.eq.12).and.(z.eq.6)) then
         capr(1) = 0.0d0
         capr(2) = 0.0d0
         capr(3) = 0.4d0
         capr(4) = 1.0d0
         capr(5) = 1.3d0
         capr(6) = 1.7d0
         capr(7) = 2.3d0
         capr(8) = 2.7d0
         capr(9) = 3.5d0
         capr(10) = 4.3d0
         capr(11) = 5.4d0
         capr(12) = 6.7d0
         
         capq(1) = 0.0d0
         capq(2) = 0.016690d0
         capq(3) = 0.050325d0
         capq(4) = 0.128621d0
         capq(5) = 0.180515d0
         capq(6) = 0.219097d0
         capq(7) = 0.278416d0
         capq(8) = 0.058779d0
         capq(9) = 0.057817d0
         capq(10) = 0.007739d0
         capq(11) = 0.002001d0
         capq(12) = 0.000007d0
         rp = 1.20d0
         p = 11
         
      else if ((a.eq.16).and.(z.eq.8)) then
         capr(1) = 0.4d0
         capr(2) = 1.1d0
         capr(3) = 1.9d0
         capr(4) = 2.2d0
         capr(5) = 2.7d0
         capr(6) = 3.3d0
         capr(7) = 4.1d0
         capr(8) = 4.6d0
         capr(9) = 5.3d0
         capr(10) = 5.6d0
         capr(11) = 5.9d0
         capr(12) = 6.4d0
         
         capq(1) = 0.057056d0
         capq(2) = 0.195701d0
         capq(3) = 0.311188d0
         capq(4) = 0.224321d0
         capq(5) = 0.059946d0
         capq(6) = 0.135714d0
         capq(7) = 0.000024d0
         capq(8) = 0.013961d0
         capq(9) = 0.000007d0
         capq(10) = 0.000002d0
         capq(11) = 0.002096d0
         capq(12) = 0.000002d0
         rp = 1.30d0
         p = 12
         
      else if ((a.eq.28).and.(z.eq.14)) then
         capr(1) = 0.4d0
         capr(2) = 1.0d0
         capr(3) = 1.9d0
         capr(4) = 2.4d0
         capr(5) = 3.2d0
         capr(6) = 3.6d0
         capr(7) = 4.1d0
         capr(8) = 4.6d0
         capr(9) = 5.1d0
         capr(10) = 5.5d0
         capr(11) = 6.0d0
         capr(12) = 6.9d0
         
         capq(1) = 0.033149d0
         capq(2) = 0.106452d0
         capq(3) = 0.206866d0
         capq(4) = 0.286391d0
         capq(5) = 0.250448d0
         capq(6) = 0.056944d0
         capq(7) = 0.016829d0
         capq(8) = 0.039630d0
         capq(9) = 0.000002d0
         capq(10) = 0.000938d0
         capq(11) = 0.000002d0
         capq(12) = 0.002366d0
         rp = 1.30d0
         p = 12
         
      else if ((a.eq.32).and.(z.eq.16))then
         capr(1) = 0.4d0
         capr(2) = 1.1d0
         capr(3) = 1.7d0
         capr(4) = 2.5d0
         capr(5) = 3.2d0
         capr(6) = 4.0d0
         capr(7) = 4.6d0
         capr(8) = 5.0d0
         capr(9) = 5.5d0
         capr(10) = 6.3d0
         capr(11) = 7.3d0
         capr(12) = 7.7d0
         
         capq(1) = 0.045356d0
         capq(2) = 0.067478d0
         capq(3) = 0.172560d0
         capq(4) = 0.324870d0
         capq(5) = 0.254889d0
         capq(6) = 0.101799d0
         capq(7) = 0.022166d0
         capq(8) = 0.002081d0
         capq(9) = 0.005616d0
         capq(10) = 0.000020d0
         capq(11) = 0.000020d0
         capq(12) = 0.003219d0
         rp = 1.35d0
         p = 12
         
      else if ((a.eq.40).and.(z.eq.20)) then
         capr(1) = 0.4d0
         capr(2) = 1.2d0
         capr(3) = 1.8d0
         capr(4) = 2.7d0
         capr(5) = 3.2d0
         capr(6) = 3.6d0
         capr(7) = 4.3d0
         capr(8) = 4.6d0
         capr(9) = 5.4d0
         capr(10) = 6.3d0
         capr(11) = 6.6d0
         capr(12) = 8.1d0
         
         capq(1) = 0.042870d0
         capq(2) = 0.056020d0
         capq(3) = 0.167853d0
         capq(4) = 0.317962d0
         capq(5) = 0.155450d0
         capq(6) = 0.161897d0
         capq(7) = 0.053763d0
         capq(8) = 0.032612d0
         capq(9) = 0.004803d0
         capq(10) = 0.004541d0
         capq(11) = 0.000015d0
         capq(12) = 0.002218d0
         rp = 1.45d0
         p = 12
         
      else if ((a.eq.48).and.(z.eq.20)) then
         capr(1) = 0.6d0
         capr(2) = 1.1d0
         capr(3) = 1.7d0
         capr(4) = 2.1d0
         capr(5) = 2.9d0
         capr(6) = 3.4d0
         capr(7) = 4.3d0
         capr(8) = 5.2d0
         capr(9) = 5.7d0
         capr(10) = 6.2d0
         capr(11) = 6.5d0
         capr(12) = 7.4d0
         
         capq(1) = 0.063035d0
         capq(2) = 0.011672d0
         capq(3) = 0.064201d0
         capq(4) = 0.203813d0
         capq(5) = 0.259070d0
         capq(6) = 0.307899d0
         capq(7) = 0.080585d0
         capq(8) = 0.008498d0
         capq(9) = 0.000025d0
         capq(10) = 0.000005d0
         capq(11) = 0.000004d0
         capq(12) = 0.001210d0
         rp = 1.45d0
         p = 12
         
      else if ((a.eq.208).and.(z.eq.82)) then
         capr(1) = 0.1d0
         capr(2) = 0.7d0
         capr(3) = 1.6d0
         capr(4) = 2.1d0
         capr(5) = 2.7d0
         capr(6) = 3.5d0
         capr(7) = 4.2d0
         capr(8) = 5.1d0
         capr(9) = 6.0d0
         capr(10) = 6.6d0
         capr(11) = 7.6d0
         capr(12) = 8.7d0
         
         capq(1) = 0.003845d0
         capq(2) = 0.009724d0
         capq(3) = 0.033093d0
         capq(4) = 0.000120d0
         capq(5) = 0.083107d0
         capq(6) = 0.080869d0
         capq(7) = 0.139957d0
         capq(8) = 0.260892d0
         capq(9) = 0.336013d0
         capq(10) = 0.033637d0
         capq(11) = 0.018729d0
         capq(12) = 0.000020d0
         rp = 1.70d0
         p = 12
         
      else
         if (prtlevel.gt.0) 
     &        write (*,*) 'dsddsogff: SOG parameters unavailable ',
     &        'for A=',a,' Z=',z
         ff = -1.d0
         return
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
      gamma = rp/(sqrt(3.0d0/2.0d0))
      f = 0.d0
      qfermi = q*fermiGeV
      
      do i = 1,p
         qr = qfermi*capr(i)
         xg = 2.d0*(capr(i)/gamma)**2
         if (qr.lt.0.000395809d0) then ! qr^16/362880 < 10^-60
            x=qr*qr
            sinxx = 1.d0+x*(-1.d0/6.d0+x*(1.d0/120.d0+x*
     &           (-1.d0/5040.d0)))
         else
            sinxx = dsin(qr)/qr
         endif
         b = capq(i)/(1.d0+xg)*(dcos(qr)+xg*sinxx)
         f = f+b
      enddo
      
      ff = f**2*exp(-0.5d0*(qfermi*gamma)**2)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end
