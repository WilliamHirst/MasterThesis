      real*8 function dsralph3(mscale)
      implicit none
      include 'dsmssm.h'
      include 'dsisasugra.h'

      real*8 mscale,dsralph31loop,dsralph34loop
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      if(roption.eq.'norun') then  
        dsralph3 = alph3mz
      elseif(roption.eq.'1loop') then
        dsralph3 = dsralph31loop(mscale) 
      elseif(roption.eq.'4loop') then
        dsralph3 = dsralph34loop(mscale) 
      elseif(roption.eq.'isasu') then
        dsralph3 = gss(3)**2/4.0d0/pi
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsralph3'
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end
