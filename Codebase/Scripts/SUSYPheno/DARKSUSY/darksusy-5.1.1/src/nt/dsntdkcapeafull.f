       real*8 function dsntdkcapeafull(mx,sigsi,sigsd)
c----------------------------------------------------------------------
c         capture rate in the earth
c         low-velocity population described by damour and krauss (1998)
c *** full: use formulas by gould as reported in jkg
c
c       mx: neutralino mass
c       sigsi: spin independent cross section in units of cm^2
c       vobs: average halo velocity
c       lars bergstrom 1998-09-15
c----------------------------------------------------------------------
       implicit none
       real*8 mx,sigsi,v_dk,v_0,delta_e,dsntcapearthfull
       real*8 sigsd,rhodk,dsntdkgtot10,v_star,v_bar
       include 'dshmcom.h'
       v_dk=3.*11.3d0 ! from dk after eq (6.11)
       v_0=sqrt(2./3.)*vobs
       delta_e=0.212/(v_0/220.d0)*dsntdkgtot10(mx,sigsi,sigsd) ! d-k eq (5.16)
c	   write(*,*) 'delta_e: ',delta_e
c       delta_e=1.           !to test against old routines
       rhodk=rhox*delta_e
c       v_star=v_dk   ! dk 1
c       v_bar=v_dk ! velocity dispersion dk 1
       v_star=39.0d0    ! dk2
       v_bar=9.2d0      ! dk2
c       v_star=sqrt(2./3.)*vobs  !to test against old routines
c       v_bar=vobs               !to test against old routines
       dsntdkcapeafull=dsntcapearthfull(mx,sigsi,v_star,v_bar,rhodk)
       return
       end
