      real*8 function dshrgacont(egath,jpsi,istat)
**********************************************************************
***   function dshrgacont gives the flux of gamma-rays with continuum
***   energy spectrum above the threshold egath (gev) from neutralino
***   annihilation in the halo.
***   
***   jpsi is the value of the integral of rho^2 along the line of
***   sight, and can be obtained with a call to dshmj.  jpsi can also be
***   the averaged value of j over a solid angle delta, obtained with a
***   call to dshmjave.
***   
***   dshrgacont uses the rescaled local density, while j uses the
***   unrescaled local density
***   
***   dshrgacont in units of cm^-2 s^-1 sr^-1
***   
***   in case of a clumpy halo the factor fdelta has to be added
***   
***   author: paolo gondolo (gondolo@mppmu.mpg.de) 00-07-19
**********************************************************************
      implicit none
      include 'dshmcom.h'
      real*8 egath,jpsi,dshrgacsusy
      integer istat
      dshrgacont=dshrgacsusy(egath,istat)*jpsi*(rhox/rho0)**2
      end


