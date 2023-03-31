***********************************************************************
*** note. this routine uses the full expressions for the capture
*** rate in the earth from gould, apj 521 (1987) 571.
*** this routine replaces dsntcapearth which use the approximations
*** given in the jkg review.
***********************************************************************


       real*8 function dsntcapearth2(mx,sigsi)
c----------------------------------------------------------------------
c         capture rate in the earth
c       uses the full routines instead of jkg (as in dsntcapearth).
c *** full: use formulas by gould as reported in jkg
c
c       mx: neutralino mass
c       sigsi: spin independent cross section in units of cm^2
c       vbar: 3D WIMP velocity dispersion in the halo
c       vstar: Sun's velocity through the halo
c       lars bergstrom 1998-09-15
c----------------------------------------------------------------------
       implicit none
       real*8 mx,sigsi,dsntcapearthfull
       real*8 v_star,v_bar
       include 'dshmcom.h'
       v_star=v_sun
       v_bar=vd_3d
       dsntcapearth2=dsntcapearthfull(mx,sigsi,v_star,v_bar,rhox)
       return
       end





