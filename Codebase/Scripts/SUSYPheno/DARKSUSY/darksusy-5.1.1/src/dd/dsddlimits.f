**********************************************************************
*** function dsddlimits gives the limits on f*sigma as a function
*** of neutralino mass. f is the halo fraction of dm (f=1 for 0.3
*** gev/cm^3) and sigma is the cross section in pb.
*** input: mx (wimp mass in gev)
***        type (1=spin-independent, 2=spin-dependent)
*** output: limit on f*sigma (pb)
*** based upon r. bernabei et al, plb 389 (1996) 757.
*** author: j. edsjo
*** date: 98-03-19
**********************************************************************

      real*8 function dsddlimits(mx,type)
      implicit none

      real*8 mx,z,zpl
      integer type,zi

      real*8 sigma(0:20,2)

      data sigma/1.0d10,1.0d10,1.0d10,1d-2,8.0d-4,
     &  3.0d-4,4.3d-5,2.0d-5,7.0d-6,5.0d-6,
     &  5.5d-6,7.5d-6,1.2d-5,1.6d-5,2.5d-5,
     &  4.0d-5,6.5d-5,1.0d-4,1.6d-4,2.3d-4,
     &  4.0d-4,
     &  21*1.0d10/

      z=log10(mx/1.0d0)*5.0
      zi=int(z)
      zpl=z-dble(zi)

      dsddlimits=10**(
     &  (1.0d0-zpl)*log10(sigma(zi,type))+
     &  zpl*log10(sigma(zi+1,type)))

      return
      end





