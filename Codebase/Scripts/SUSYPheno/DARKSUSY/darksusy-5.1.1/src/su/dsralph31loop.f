      real*8 function dsralph31loop(mscale)
      implicit none
      include 'dsmssm.h'
      real*8 mscale,alph3mb,alph3mt
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      logical first
      data first/.true./
      save first
c
      if(mscale.le.mtmt.and.mscale.ge.mbmb) then
        dsralph31loop=alph3mz/(1.d0+(11.d0-10.d0/3.d0)
     &    *alph3mz/(4.d0*pi)*2.d0*dlog(mscale/mass(kz)))
      elseif(mscale.lt.mbmb.and.mscale.ge.mcmc) then
        alph3mb=alph3mz/(1.d0+(11.d0-10.d0/3.d0)
     &    *alph3mz/(4.d0*pi)*2.d0*dlog(mbmb/mass(kz)))
        dsralph31loop=alph3mb/(1.d0+(11.d0-8.d0/3.d0)
     &    *alph3mb/(4.d0*pi)*2.d0*dlog(mscale/mbmb))
      elseif(mscale.gt.mtmt) then
        alph3mt=alph3mz/(1.d0+(11.d0-10.d0/3.d0)
     &    *alph3mz/(4.d0*pi)*2.d0*dlog(mtmt/mass(kz)))
        dsralph31loop=alph3mt/(1.d0+(11.d0-12.d0/3.d0)
     &    *alph3mt/(4.d0*pi)*2.d0*dlog(mscale/mtmt))
      else
        if (first) then
          write(*,*) ' ' 
          write(*,*) 'WARNING: ',
     &      'dsralph31loop called with mscale out of scale'
          write(*,*) 'mscale = ',mscale
          write(*,*) 'Using scale ',mcmc,' instead.'
          write(*,*) 'This warning will only be printed once.'
          first=.false.
        endif
        alph3mb=alph3mz/(1.d0+(11.d0-10.d0/3.d0)
     &    *alph3mz/(4.d0*pi)*2.d0*dlog(mbmb/mass(kz)))
        dsralph31loop=alph3mb/(1.d0+(11.d0-8.d0/3.d0)
     &    *alph3mb/(4.d0*pi)*2.d0*dlog(mcmc/mbmb))
      endif
      return
      end
