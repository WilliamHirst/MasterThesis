       real*8 function dsntsefull(mx,m_a,v_star,v_bar,phi)
c--------------------------------------------------------------------------
c the gould function for capture in the earth, as given by the expression
c (a10) in gould ap.j. 321 (1987) 571
c mx: neutralino mass in gev
c m_a: nuclear mass
c v_star: velocity of earth with respect to wimps
c v_bar: velocity dispersion of wimps
c output in natural units (power of gev)
c l. bergstrom 1998-09-21
c Modified by J. Edsjo, 2003-11-22
c References:
c    dk: Damour and Krauss, Phys. Rev. D59 (1999) 063509.
c   jkg: Jungman, Kamionkowski and Griest, Phys. Rep. 267 (1996) 195.
c gould: Gould, ApJ 321 (1987) 571.
c--------------------------------------------------------------------------
       implicit none
       real*8 mx,v_star,v_bar,dsntmoderf,apil
       real*8 eta,etahat,ahatpl,ahatmin
       real*8 etapil,apilpl,apilmin,a0,b0,m_a,betaplus,v_02,q_a,
     1      r_a2,betaminus,v_02c,
     2      aa,ahat,c,vesc,phi
       real*8 sla1
       include 'dshmcom.h'
       parameter(c=299792.458d0) ! km/s

       eta=sqrt(3./2.)*v_star/v_bar
       betaplus=4.d0*mx*m_a/(mx+m_a)**2   !
       betaminus=4.d0*mx*m_a/(mx-m_a)**2  !
c...Nuclear radius (squared), gould (a8) or dk (2.28) in gev**-2
       r_a2=25.7d0*(0.3d0+0.91*m_a**0.3333d0)**2
c...Characteristic coherence energy, in gould it is eq. (A4) and it is
c...there called E0. It is equivalent to Q_A in dk (2.28). Note however that
c...q_a = 2*Q0 in jkg (9.12)
       q_a=3./2./m_a/r_a2
       a0=mx*(v_bar/c)**2/(3.0d0*q_a)    ! gould (A7)
       b0=betaplus*a0            ! gould (A7), jkg (9.16)
       etahat=eta/sqrt(1.+a0)    ! gould (A11), jkg (9.16)
       etapil=eta/sqrt(1.+b0)    ! gould (A11), jkg (9.17)
c...Escape velocity dependent on typical potential
c...This in effect means that we approximate all atoms of each element as
c...being located at the typical potential phi
       vesc=11.2d0*sqrt(phi)

       aa=sqrt(3.0d0/2.0d0*betaminus*(vesc/v_bar)**2)! gould (2.22)
c
       ahat=aa*sqrt(1.+a0)   ! gould (A12)
       apil=aa*sqrt(1.+b0)   ! gould (A12)
       ahatpl=ahat+etahat    ! gould (A13)
       ahatmin=ahat-etahat   ! gould (A13)
       apilpl=apil+etapil    ! gould (A13)
       apilmin=apil-etapil   ! gould (A13)
       sla1=exp(-a0*etahat**2)/sqrt(1.+a0)*(2.*dsntmoderf(etahat)+
     1      dsntmoderf(ahatmin)-dsntmoderf(ahatpl))-
     2      exp(-b0*etapil**2)/sqrt(1+b0)
     3        *exp(-(a0-b0)*aa**2)*
     4       (2.*dsntmoderf(etapil)+
     5      dsntmoderf(apilmin)-dsntmoderf(apilpl))
       sla1=sla1/2./eta/b0*sqrt(8./3./3.141592)
       dsntsefull=sla1  ! gould (A10), apart from sigma n n_w v_bar
c       write(*,*) 'dsntsefull: ',mx/m_a,v_star,v_bar,dsntsefull
        return
        end
