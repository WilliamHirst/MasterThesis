      function dsbsgg8chi2(x)

***********************************************************************
* Function G^{chi,2}_8(x) in eq. (18) of                              *
* Ciuchini et al., hep-ph/9806308                                     *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 x,x1s,x2b
      real*8 ddilog
      real*8 dsbsgf83,dsbsgd1td,dsbsgd2d
      real*8 dsbsgg8chi2

c     Mass-squared ratios for the functions \Delta^(1)_{t,d}
c     and \Delta^(1)_d 
c     Note that we use the same fortran code for these 
c     \Delta's, but different input parameters 

      if (mass(kss1).gt.mass(kss2)) then
        x1s=mass(kss1)**2/mass(kgluin)**2
      else
        x1s=mass(kss2)**2/mass(kgluin)**2
      endif

      if (mass(ksb2).lt.mass(ksb1)) then
        x2b=mass(ksb2)**2/mass(kgluin)**2
      else
        x2b=mass(ksb1)**2/mass(kgluin)**2
      endif
      
      dsbsgg8chi2=-(4.d0/3.d0)*dsbsgf83(x)
     &    *(dsbsgd1td(x1s)+dsbsgd1td(x2b)
     &      +dsbsgd2d()-2.d0)
     &   -4.d0*x*(3.d0+4.d0*x)*ddilog(1.d0-1.d0/x)
     &    /(3.d0*(x-1.d0)**3)
     &   -(61.d0-39.d0*x)/(12.d0*(x-1.d0)**2)
     &   -(7.d0-60.d0*x-14.d0*x**2)
     &    *log(x)/(6.d0*(x-1.d0)**3)
     &   -14.d0*x*(log(x))**2/(3.d0*(x-1.d0)**3)


      return
      end


