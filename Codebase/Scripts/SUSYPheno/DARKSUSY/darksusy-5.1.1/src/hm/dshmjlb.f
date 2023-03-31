      function dshmjlb(l,b)
c
c     Astrophysical factor J(l,b) as a function of galactic longitude l and latitude b
c
      implicit none
      real*8 dshmjlb,l,b,dshmj,cospsi
      cospsi=cos(l)*cos(b)
      dshmjlb=dshmj(cospsi)
      return
      end
