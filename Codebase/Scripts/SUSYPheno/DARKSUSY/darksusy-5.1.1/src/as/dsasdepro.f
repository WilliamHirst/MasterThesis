**************************************************************
*** FUNCTION dsasdepro                                     ***
*** computes the denominator of a propagator               ***
***                                                        ***
*** input: mom2 is S,T,U;                                  ***
*** kkpp is the number of the particle in the propagator   ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************

      complex*16 function dsasdepro(mom2,kkpp)
      implicit none
      include 'dsmssm.h'
      real*8 mom2,intmass,intwid,diff,gammam,res1,res2
      integer kkpp
      intmass=mass(kkpp)
      intwid=width(kkpp)
      diff=mom2-intmass**2
      gammam=intmass*intwid
      if(dabs(diff).gt.1.d-3*gammam) then
        res1=diff+gammam**2/diff
        res1=1.d0/res1
        res2=diff**2+gammam**2
        res2=gammam/res2
      elseif(dabs(diff).lt.1.d3*gammam) then
        res1=diff**2+gammam**2
        res1=diff/res1
        res2=diff**2/gammam+gammam
        res2=1.d0/res2
      else
        res1=diff**2+gammam**2
        res2=res1
        res1=diff/res1
        res2=gammam/res2
      endif
      dsasdepro=res1-res2*dcmplx(0.d0,1.d0)
      return
      end








