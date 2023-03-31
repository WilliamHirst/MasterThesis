
C-------------------------------------------------------
C#######################################################
C-------------------------------------------------------
      Function FA12_integr(a,b1,b2,xF)
C-----------------------------------------------------------------------
C     Function to calculate F_{12} for the velocity averaged
C     co-annihilation cross section.
C       x  = T/m is already integrated over analytically.
C       a  = Sqrt(s)/m.
C       bi = m_i/m.
C       xF is the upper limit of the x-integral.
C       m is a convenient scale (can be choosen e.g. as (m_1+m_2)/2).
C-----------------------------------------------------------------------
      Implicit none
      Real*8 FA12_integr,a,b1,b2,xF, a12,A1, Pi, dErfc
      Data Pi / 3.14159265358979323846D0 /

      a12 = a - b1 - b2
      A1  = dErfc(Sqrt(a12/xF))

      FA12_integr = A1*Sqrt(Pi/a12) + 
     -  2.*(3./(8.*a) - 15./(8.*b1) - 15./(8.*b2))*
     -   (-A1*Sqrt(a12*Pi) + Sqrt(xF)/Exp(a12/xF)) + 
     -  (2.*(-15./(128.*a**2) + 345./(128.*b1**2) - 45./(64.*a*b1) + 
     -       345./(128.*b2**2) - 45./(64.*a*b2) + 225./(64.*b1*b2))*
     -     (2.*A1*a12**1.5*Sqrt(Pi) + 
     -       (-2.*a12*Sqrt(xF) + xF**1.5)/Exp(a12/xF)))/3.
     

      Return
      End   
