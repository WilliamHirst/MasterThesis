      function dsrdfuncs(u)
c_______________________________________________________________________
c  10^15 * dsrdfunc.
c  input:
c    u - integration variable
c  uses dsrdfunc
c  used for gaussian integration with gadap.f
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-01-17
c=======================================================================
      implicit none
      real*8 dsrdfuncs,u
      real*8 dsrdfunc,dsrdwintp
      external dsrdfunc,dsrdwintp

c-----------------------------------------------------------------------
      real*8 rdx
      common /gadint2/ rdx
c-----------------------------------------------------------------------

      dsrdfuncs=dsrdfunc(u,rdx,dsrdwintp)*1.0d15

      end








































