      program caprates

      implicit none
      include 'dsntcom.h'

      real*8 cap_si,cap_sd,cap_sinew,cap_sdnew
      real*8 dsntcapsunnum,dsntcapsuntab
      real*8 mx,sigsip,sigsdp

      integer i

      call dsinit
c      call dshmset('burksm')
      call dsntset('tab')

      sigsip=1.d-36*1.d-7
      sigsdp=1.d-36*1.d-3
      open(unit=67,file='caprates-comb.dat',
     &  form='formatted',status='unknown')
      do i=0,50
         mx=10**(dble(i)/10.d0+1)
         call dsntset('tab')
c         cap_si=dsntcapsuntab(mx,sigsip,0.d0)
c         cap_sd=dsntcapsuntab(mx,0.d0,sigsdp)
         cap_si=dsntcapsunnum(mx,sigsip,0.d0)
         cap_sd=dsntcapsunnum(mx,0.d0,sigsdp)
         call dsntset('tabcut')
c         cap_sinew=dsntcapsuntab(mx,sigsip,0.d0)
c         cap_sdnew=dsntcapsuntab(mx,0.d0,sigsdp)
         cap_sinew=dsntcapsunnum(mx,sigsip,0.d0)
         cap_sdnew=dsntcapsunnum(mx,0.d0,sigsdp)
         write(*,100) i,mx,cap_si,cap_sd,cap_sinew,cap_sdnew
         write(67,100) i,mx,cap_si,cap_sd,cap_sinew,cap_sdnew
      enddo
      
 100  format(1x,I2,1x,5(E12.6,1x))
      end
