*****************************************************************************
*** integrates a function f(z,y) over the kinematic variables
***
*** z = E_p/mx
*** y = (p+k)^2/(4 mx^2)
***
*** that appear for 3-body final states. (k is the photon momentum and p
*** that of one of the other final sates). The result is the
*** total p yield per annihilation
***
*** called by dsIByieldone
***
*** author: Torsten Bringmann (troms@physto.se)
*** date  : 2007-10-20
*****************************************************************************
      real*8 function dsIBf_intdxdy2(IBch,z,mwimp,mp1,mp2)
      implicit none
      include 'dsprep.h'  
      include 'dsibcom.h'

      real*8 a,b,z,mwimp,mp1,mp2
      real*8 dsf_int, dsIBintsel2
      external dsIBintsel2
      integer IBch

      dsIBf_intdxdy2=0D0

c... set parameters needed by dsIBintsel
      if (IBch.eq.4) then
          intch = -4
      else
          return
      endif
      ibcom_mx=mwimp
      ibcom_mp1=mp1
      ibcom_mp2=mp2

     
c...determine integration limits
      a=z
      b=0.999d0   ! be aware that very high positron energies are enhanced by IR 
                  ! photons a full treatment would have to include virtual photons
      if (a.ge.b) then
         return
      endif

      dsIBf_intdxdy2=dsf_int(dsIBintsel2,a,b,IBacc)

      end
