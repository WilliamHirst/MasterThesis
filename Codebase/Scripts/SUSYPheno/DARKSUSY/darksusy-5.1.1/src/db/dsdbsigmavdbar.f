      real*8 function dsdbsigmavdbar(en)
c total inelastic cross section dbar + h
c Yad.Fiz.14:134-136,1971
c check Review of Particle Properties of 1992, rev [9] in DFS
      implicit none
      include 'dsmpconst.h'
      real*8 en,vp,sigma,cc
      parameter (cc=c_light/1.d5)
c v= c**2 * p/en
      vp=dsqrt(en**2-m_d**2)/en
      sigma=109.9d0
c... this is the only data point at P=13.3 GeV, retrieved from:
c... http://durpdg.dur.ac.uk/cgi-hepdata/hepreac/3410480
c sigma in mb
      dsdbsigmavdbar=sigma*(vp*cc)
      return
      end
