      real*8 function dsrmq(mscale,kpart)
      implicit none
      include 'dsmssm.h'
      include 'dsisasugra.h'

      real*8 mscale,dsrmq1loop,dsrmq4loop,cosbeta
      integer kpart

      if(roption.eq.'norun') then  
        if(kpart.eq.kc) then
          dsrmq=mass(kc)
        elseif(kpart.eq.kb) then
          dsrmq=mass(kb)
        elseif(kpart.eq.kt) then
          dsrmq=mass(kt)
        elseif(kpart.eq.ktau) then
          dsrmq=mass(ktau)
        else  
          write(*,*) 'dsrmq called for wrong particle'
          write(*,*) 'particle =  ',pname(kpart)
          write(*,*) 'rather than tau or c, b, t quark'
          write(*,*) 'program stopped'
          stop
        endif  
      elseif(roption.eq.'1loop') then
        if(kpart.eq.kc.or.kpart.eq.kb.or.kpart.eq.kt) then
          dsrmq=dsrmq1loop(mscale,kpart) 
        elseif(kpart.eq.ktau) then
          dsrmq=mass(ktau)
        else  
          write(*,*) 'dsrmq called for wrong particle'
          write(*,*) 'particle =  ',pname(kpart)
          write(*,*) 'rather than tau or c, b, t quark'
          write(*,*) 'program stopped'
          stop
        endif
      elseif(roption.eq.'4loop') then
        if(kpart.eq.kc.or.kpart.eq.kb.or.kpart.eq.kt) then
          dsrmq=dsrmq4loop(mscale,kpart) 
        elseif(kpart.eq.ktau) then
          dsrmq=mass(ktau)
        else  
          write(*,*) 'dsrmq called for wrong particle'
          write(*,*) 'particle =  ',pname(kpart)
          write(*,*) 'rather than tau or c, b, t quark'
          write(*,*) 'program stopped'
          stop
        endif
      elseif(roption.eq.'isasu') then  
        cosbeta=1.0d0/sqrt(1.0d0+tanbe*tanbe)
        if(kpart.eq.kc) then
          dsrmq=mass(kc)
        elseif(kpart.eq.kb) then
          dsrmq=gss(5)*cosbeta*vev
        elseif(kpart.eq.kt) then
          dsrmq=gss(6)*tanbe*cosbeta*vev
        elseif(kpart.eq.ktau) then
          dsrmq=gss(4)*cosbeta*vev
        else  
          write(*,*) 'dsrmq called for wrong particle'
          write(*,*) 'particle =  ',pname(kpart)
          write(*,*) 'rather than tau or c, b, t quark'
          write(*,*) 'program stopped'
          stop
        endif  
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsrmq'
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end
