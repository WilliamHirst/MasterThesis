      real*8 function dsntdkcapea(mx,sigsi,sigsd)
c----------------------------------------------------------------------
c         capture rate in the earth
c         low-velocity population described by damour and krauss (1998)
c *** simple version: only change vobs -> v_dk \sim 3*v_esc_earth
c *** in a factor (9.22) in jkg
c         based on jungman, kamionkowski, griest review
c       mx: neutralino mass
c       sigsi: spin independent cross section in units of cm^2
c       vobs: average halo velocity
c       lars bergstrom 1998-09-15
c----------------------------------------------------------------------
      implicit none
      real*8 mx,sigsi,dsntlitlf_e,v_dk,v_0,delta_e,dsntdkgtot10
      real*8 sigsd
      include 'dshmcom.h'
      include 'dsmpconst.h'
      v_dk=3.*11.3d0 ! from dk after eq (6.11)
      v_0=sqrt(2./3.)*vobs
      delta_e=0.212/(v_0/220.d0)*dsntdkgtot10(mx,sigsi,sigsd) ! d-k eq (5.16)
c	   write(*,*) 'delta_e: ',delta_e
      dsntdkcapea=sigsi*1.0d40*2.57d-13*3.1415/4./m_p**2 ! from sigsi to f_p
      dsntdkcapea=dsntdkcapea*dsntlitlf_e(mx,v_dk)*
     &  2.4d28*(delta_e*rhox/0.3)
c  this is the staandard expression (vobs instead of v_dk):
c       dsntdkcapea=dsntdkcapea*dsntlitlf_e(mx,vobs)*2.4d28*(delta_e*rhox/0.3)
      return
      end
