      subroutine dsfindmtmt
      implicit none
      include 'dsmssm.h'
c
      if(roption.eq.'norun') then  
        mtmt=mass(kt)
      elseif(roption.eq.'1loop') then
        call dsfindmtmt1loop
      elseif(roption.eq.'4loop') then
        call dsfindmtmt4loop
      elseif(roption.eq.'isasu') then
        call dsfindmtmt4loop ! we do not use it, so anything goes
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsralph3'
        write(*,*) 'program stopped'
        stop
      endif

      return
      end  
