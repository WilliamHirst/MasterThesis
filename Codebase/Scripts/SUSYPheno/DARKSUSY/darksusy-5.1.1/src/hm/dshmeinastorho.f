****************************************************************
*** General Einasto dark matter halo density profile,        *** 
*** intended mostly for dwarf galaxies.                      ***
*** This is a generalised (earlier) version of the n03 halo  ***
*** model.                                                   *** 
***                                                          ***  
*** radialdist = distance from centre of halo in kpc         ***
*** rs_e  = halo scale radius in kpc                         ***
*** rhos_e = halo scale density in gev/cm**3                  ***
*** n_e   = Einasto index                                    ***
***                                                          ***
*** Output: density in gev/cm**3                             ***
***                                                          ***
*** Author: Pat Scott (pat@fysik.su.se)                      ***      
*** Date: 2009-05-07                                         ***
****************************************************************

      real*8 function dshmeinastorho(radialdist)
      implicit none
      real*8 radialdist,x
      include 'dshmcom.h'
      x=radialdist/rs_e
      dshmeinastorho=rhos_e*dexp(-2.d0*n_e*(x**(1.d0/n_e)-1.d0))
      return
      end

