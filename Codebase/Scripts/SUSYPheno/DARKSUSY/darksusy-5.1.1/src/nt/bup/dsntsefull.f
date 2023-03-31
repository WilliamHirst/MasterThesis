

        real*8 function dsntsefull(mx,m_a,v_star,v_bar)
c--------------------------------------------------------------------------
c the gould function for capture in the earth, as given by the expression
c (a10) in gould ap.j. 321 (1987) 571
c mx: neutralino mass in gev
c m_a: nuclear mass
c v_star: velocity of earth with respect to wimps
c v_bar: velocity dispersion of wimps
c output in natural units (power of gev)
c l. bergstrom 1998-09-21
c--------------------------------------------------------------------------
        implicit none
        real*8 mx,v_star,v_bar,dsntmoderf,apil
        real*8 eta,etahat,ahatpl,ahatmin
        real*8 etapil,apilpl,apilmin,a0,b0,m_a,betaplus,v_02,q_a,
     1       r_a2,betaminus,v_02c,
     2       aa,ahat
        real*8 sla1
        include 'dshmcom.h'
       eta=sqrt(3./2.)*v_star/v_bar
       betaplus=4.d0*mx*m_a/(mx+m_a)**2   !
       betaminus=4.d0*mx*m_a/(mx-m_a)**2  !
       v_02=2./3.*v_bar**2   ! v_bar \sim 270 km/s gives v_0 \sim 220 km/s
       v_02c=v_02/(300000.)**2 ! use units where c=1
       r_a2=25.7d0*(0.3d0+0.91*m_a**0.3333d0)**2 !(dk 2.26) in gev**-2
       q_a=3./2./m_a/r_a2                        !(dk 2.26)
       a0=mx*v_02c/2./q_a                        !(jkg 9.16)
       b0=betaplus*a0
c       write(*,*) 'a0,b0: ',a0,b0
       etahat=eta/sqrt(1.+a0)
       etapil=eta/sqrt(1.+b0)
       aa=sqrt(betaminus*(13.2/v_bar)**2)!effective v_esc, cf after jkg (9.22)
c
       ahat=aa*sqrt(1.+a0)
       apil=aa*sqrt(1.+b0)
       ahatpl=ahat+etahat
       ahatmin=ahat-etahat
       apilpl=apil+etapil
       apilmin=apil-etapil
       sla1=exp(-a0*etahat**2)/sqrt(1.+a0)*(2.*dsntmoderf(etahat)+
     1      dsntmoderf(ahatmin)-dsntmoderf(ahatpl))-
     2      exp(-b0*etapil**2)/sqrt(1+b0)*exp(-(a0-b0)*aa**2)*
     3       (2.*dsntmoderf(etapil)+
     4      dsntmoderf(apilmin)-dsntmoderf(apilpl))
       sla1=sla1/2./eta/b0*sqrt(8./3./3.141592)
       dsntsefull=sla1
c       write(*,*) 'dsntsefull: ',mx/m_a,v_star,v_bar,dsntsefull
        return
        end
