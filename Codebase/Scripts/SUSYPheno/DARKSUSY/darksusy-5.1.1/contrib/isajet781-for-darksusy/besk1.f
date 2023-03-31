
C-------------------------------------------------------
C#######################################################
C-------------------------------------------------------
C-------------------------------------------------------
      Function BesK1(z)
C-----------------------------------------------------------------------
C     The modified Bessel function of the second kind K(n, z)
C     for large arguments z.
C     These functions have a precision better then 1(.1) % for z > 1.5(3).
C     For more details see my notes on thermal average. CsB 12/2001
C-----------------------------------------------------------------------
      Implicit none
      Real*8 BesK1,BesK2,z,Pi
      Integer n
      Data Pi / 3.14159265358979323846D0 /

      BesK1 = (Sqrt(Pi/2.)*(1./z)**3.5*
     -    (105. - 120.*z + 384.*z**2 + 1024.*z**3))/(1024.)

      ENTRY BesK2(z)
      BesK2 = (Sqrt(Pi/2.)*(1./z)**3.5*
     -    (-315. + 840.*z + 1920.*z**2 + 1024.*z**3))/(1024.)

      If (z.LT.1.5) Print*, '!!! WARNING: Precision < 1% in BesK !!!'
      If (z.LT.3.) Print*,  '!!! WARNING: Precision < .1% in BesK !!!'
     
      Return
      End
