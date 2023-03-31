***********************************************************************
*** dsntearthne gives the number density of electrons as a function
*** of the Earth's radius.
*** Input: Earth radius [m]
*** Output: n_e [cm^-3]
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2006-04-21
***********************************************************************

      real*8 function dsntearthne(r)
      implicit none

      include 'dsmpconst.h'

      real*8 dsntearthdens,r

      dsntearthne=dsntearthdens(r)*n_avogadro/2.0d0  ! number of electrons / cm^3

      return

      end
