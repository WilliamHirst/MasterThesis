*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                        dsmpconst.h                               ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c mathematical and physical constants, p. gondolo 2011-11-03

      real*8 m_p,m_n,m_d,n_avogadro,c_light,pi,
     &     gev2cm3s,fermiGeV,gev2cm2,atomicmassunit,m_p_amu,m_n_amu
      parameter(
c     &     pi=3.1415926535897932384626433832795029D0, ! pi
     &     pi=4.0d0*datan(1.d0), ! pi
     &     m_p_amu=1.00727646688d0, ! proton mass [amu]
     &     m_n_amu=1.0086649156d0, ! neutron mass [amu]
     &     atomicmassunit=0.931494028d0, ! atomic mass unit [GeV/c^2]
c     &     m_p=0.938271998d0, ! proton mass [GeV/c^2]
c     &     m_n=0.9396d0,   ! neutron mass [GeV/c^2]
     &     m_p=m_p_amu*atomicmassunit, ! proton mass [GeV/c^2]
     &     m_n=m_n_amu*atomicmassunit,   ! neutron mass [GeV/c^2]
     &     m_d=1.875612762d0, ! deuteron mass [GeV/c^2]
     &     n_avogadro=6.022d23,  ! Avogradros number [number / mol ]
     &     c_light=299792.458d0, ! speed of light [km/s]
     &     gev2cm3s=0.38937966d-27*3.d10, ! conversion [GeV^2 cm^3/s]
     &     fermiGeV=1.d0/0.1973269602d0, ! conversion [GeV fm]
     &     gev2cm2=(197.327053d-16)**2) ! conversion [GeV^2 cm^2]

***                                                                 ***
*********************** end of dsmpconst.h ****************************
