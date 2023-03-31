      function dsbsgg8chij1(x,j)

***********************************************************************
* Function G^{chi,1}_8(x) in eq. (17) of                              *
* Ciuchini et al., hep-ph/9806308                                     *
* The expression has been extended to large tanbe, by replacing       *
* ln(m^2(kgluin)/m^2(\chi_j)) by ln((mu_w)^2)/m^2(\chi_j))            *
* as explained in Degrassi et al., hep-ph/0009337 p.11                *
* The input, j, is the nb of the chargino, it must be 1 or 2          *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-04-07                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      integer j
      real*8 x,x1s,x1b,muw,mt
      real*8 ddilog
      real*8 dsbsgf81,dsbsgd1td
      real*8 dsbsgg8chij1

c     The mass scale mu_w, here set to mt 

      mt=mtmt  ! mt at scale mt, value from  DarkSUSY
      muw=mt


c     Mass-squared ratios for the function \Delta^(1)_{t,d}


      if (mass(kss1).gt.mass(kss2)) then
        x1s=mass(kss1)**2/mass(kgluin)**2
      else
        x1s=mass(kss2)**2/mass(kgluin)**2
      endif

      if (mass(ksb2).lt.mass(ksb1)) then
        x1b=mass(ksb1)**2/mass(kgluin)**2
      else
        x1b=mass(ksb2)**2/mass(kgluin)**2
      endif

      dsbsgg8chij1=-(8.d0/9.d0)*dsbsgf81(x)
     &    *(dsbsgd1td(x1s)+dsbsgd1td(x1b)-1.d0
     &      +3.d0*log(muw**2/mass(kcha(j))**2))
     &   -x*(1210.d0-437.d0*x-1427.d0*x**2)/(648.d0*(x-1.d0)**3)
     &   -x**2*(49.d0+46.d0*x+9.d0*x**2)
     &    *ddilog(1.d0-1.d0/x)/(12.d0*(x-1.d0)**4)
     &   -x*(85.d0-603.d0*x-387.d0*x**2+78*x**3)
     &    *log(x)/(108.d0*(x-1.d0)**4)
     &   -13.d0*x**2*(log(x))**2/(3.d0*(x-1.d0)**4)


      return
      end


