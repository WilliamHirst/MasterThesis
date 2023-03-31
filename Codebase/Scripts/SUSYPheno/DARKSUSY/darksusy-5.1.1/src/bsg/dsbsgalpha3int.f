      function dsbsgalpha3int(al,mstart,m,nf)

***********************************************************************
* The coupling constant alpha_3 evaluated at the scale m              *
* given the value at a given scale mstart. nf effective quark flavours*
* are used in the running (if nf=7, 6 quark flavours and one squark   *
* flavor are used)                                                    *
* Uses eq. (42) of Ciuchini et al. hep-ph/9710335                     *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      integer nf
      real*8 m,mstart,al
      real*8 b0,b1,v 
      real*8 dsbsgalpha3int
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     define the numbers beta_0 and beta_1 from eq. (33)
c     depend on the input number n_f of eff. quark flavours

      if (nf.le.6) then

        b0=11.d0-2.d0*dble(nf)/3.d0
     
        b1=102.d0-38.d0*dble(nf)/3.d0

      else

        b0=41.0d0/6.0d0                ! UPDATE what should this be?
        b1=102.d0-38.d0*dble(6)/3.d0   ! UPDATE what should this be?

      endif


c     define the function v of eq. (42)
c     depends on the input scale m
c     and on the input n_f through beta_0      

      v=1.d0+b0*(al/(2.d0*pi))*log(m/mstart)

      dsbsgalpha3int=(al/v)*(1.d0-
     &       b1*al*log(v)/(b0*4.d0*pi*v))

      return
      end


