***********************************************************************
*** note. this routine assumes a maxwell-boltzmann velocity
*** distribution and uses approximations in the jkg review,
*** Jungman, Kamionkowski and Griest, Phys. Rep. 267 (1996) 195.
*** In particular, it is assumed that the Sun's velocity is
*** sqrt(2/3)*vd_3d.
*** for more accurate results, use dsntcapearthfull instead.
*** for an arbitrary velocity distribution, use dsntcapearthnumi instead.
***********************************************************************

       real*8 function dsntcapearth(mx,sigsi)
c----------------------------------------------------------------------
c         capture rate in the earth
c         based on jungman, kamionkowski, griest review
c       mx: neutralino mass
c       sigsi: spin independent cross section in units of cm^2
c       vobs: average halo velocity
c       lars bergstrom 1995-12-14
c----------------------------------------------------------------------
       implicit none
       include 'dsmpconst.h'
       real*8 mx,sigsi,dsntlitlf_e
       include 'dshmcom.h'

c...Eq. (9-27) in jkg with the assumption m_chi >> m_proton
       dsntcapearth=sigsi*1.0d40*2.57d-13*3.1415/4./m_p**2 ! from sigsi to f_p
       dsntcapearth=dsntcapearth*dsntlitlf_e(mx,vd_3d)*2.4d28*(rhox/0.3)
       return
       end
