       real*8 function dsntcapsun(mx,sigsi,sigsd)
c----------------------------------------------------------------------
c         capture rate in the sun
c         based on jungman, kamionkowski, griest review
c       mx: neutralino mass
c       sigsi: spin independent cross section in units of cm^2
c       sigsd: spin dependent cross section in units of cm^2
c       vobs: average halo velocity
c       output:
c         capture rate in s^-1
c       lars bergstrom 1995-12-12
c----------------------------------------------------------------------
       implicit none
       include 'dsmpconst.h'
       real*8 mx,sigsi,sigsd,dsntlitlf_s,dsntss
       include 'dshmcom.h'
       dsntcapsun=sigsi*1.0d40*2.57d-13*3.1415/4./m_p**2 ! from sigsi to f_p
       dsntcapsun=dsntcapsun*dsntlitlf_s(mx,vd_3d)*2.4d37
c   add spin dependent piece:
       dsntcapsun=dsntcapsun+
     &   1.3d25*sigsd*1.0d40*dsntss(mx/m_p,vd_3d)/(mx*vd_3d/270.d0)
       dsntcapsun=dsntcapsun*(rhox/0.3)
       return
       end








