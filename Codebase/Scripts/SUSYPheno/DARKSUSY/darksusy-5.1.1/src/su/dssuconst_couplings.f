      subroutine dssuconst_couplings
c_______________________________________________________________________
c  useful constants
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa fixed (pg)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      real*8 s12,s23,s13,c12,c13,c23,d13,aux
      real*8 dsgf2s2thw
      complex*16 ed13
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      integer i

c...Determine sin^2 theta_W from GFermi, alphem at MZ (MS-bar) and MZ
      s2thw=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),1)

c...Now calculate sinthw and costhw from inputs
      sinthw=sqrt(s2thw)
      costhw=sqrt(1.0d0-s2thw)

c...Also calculate expressions at mZ
      s2wmz=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),3)
      swmz=sqrt(s2wmz)
      cwmz=sqrt(1.0d0-s2wmz)

      g2weak=sqrt(4.0d0*pi*alphem/s2thw)
      gyweak=g2weak*sinthw/costhw
      g3stro=sqrt(4.d0*pi*alph3mz)

c...And weak couplings at MZ
      g2wmz=sqrt(4.0d0*pi*alphem/s2wmz)
      gywmz=g2wmz*swmz/cwmz
      end
