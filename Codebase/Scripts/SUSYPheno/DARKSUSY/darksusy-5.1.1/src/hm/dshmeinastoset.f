****************************************************************
*** Parameter resetting routine for the Einasto profile      *** 
***                                                          ***  
*** rs_e   = halo scale radius in kpc                        ***
*** rhos_e = halo scale density in gev/cm**3                 ***
*** n_e    = Einasto index                                   ***
*** r_0    = galactocentric distance                         ***
***                                                          ***
*** If any of the input arguments are negative, the          ***
*** respective parameter is left unchanged.                  ***
***                                                          ***
*** Author: Pat Scott (pat@fysik.su.se)                      ***      
*** Date: 2009-05-07                                         ***
****************************************************************

      subroutine dshmeinastoset(rs_in, rho_in, n_in, r_0_in)
      implicit none
      real*8 rs_in, rho_in, n_in, r_0_in
      include 'dshmcom.h'
      if (rs_in .gt. 0.d0) rs_e = rs_in       
      if (rho_in .gt. 0.d0) rhos_e = rho_in       
      if (n_in .gt. 0.d0) n_e = n_in       
      if (r_0_in .gt. 0.d0) r_0 = r_0_in       
      return
      end

