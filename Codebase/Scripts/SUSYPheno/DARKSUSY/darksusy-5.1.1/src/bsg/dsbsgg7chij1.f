      function dsbsgg7chij1(x,j)

***********************************************************************
* Function G^{chi,1}_7(x) in eq. (15) of                              *
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
      real*8 dsbsgf71,dsbsgd1td
      real*8 dsbsgg7chij1



      mt=mtmt  ! mt at scale mt, value from DarkSUSY
      muw=mtmt


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

      dsbsgg7chij1=-(8.d0/9.d0)*dsbsgf71(x)
     &    *(dsbsgd1td(x1s)+dsbsgd1td(x1b)-1.d0
     &      +3.d0*log(muw**2/mass(kcha(j))**2))
     &   +x*(85.d0-347.d0*x+526.d0*x**2)/(243.d0*(1.d0-x)**3)
     &   +4.d0*x**2*(-8.d0+13.d0*x+6.d0*x**2)
     &    *ddilog(1.d0-1.d0/x)/(9.d0*(x-1.d0)**4)
     &  -4.d0*x*(20.d0-126.d0*x+144.d0*x**2+39*x**3)
     &    *log(x)/(81.d0*(x-1.d0)**4)
     &  +2.d0*x**2*(21.d0*x-10.d0)
     &    *(log(x))**2/(9.d0*(x-1.d0)**4)


      return
      end


