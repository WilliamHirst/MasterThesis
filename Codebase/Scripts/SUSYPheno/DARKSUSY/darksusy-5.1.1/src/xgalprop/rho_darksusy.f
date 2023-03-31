**********************************************************************
*** subroutine rho_darksusy
*** input: x, y, z - galactic coordinates in kpc
*** output: rho - halo density in GeV cm^-3
*** author: e.a. baltz 2/21/2006
**********************************************************************

      subroutine rho_darksusy(x,y,z,rho)
      implicit none

      real*8 x,y,z,rho,r,dshmrho
      r=sqrt(x**2+y**2+z**2)

      rho=dshmrho(r)
      
      return
      end
