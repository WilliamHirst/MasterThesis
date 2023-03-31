      subroutine dsrdinit
c...initialize relic density routines
c...for now, just read the dof table
c...author: paolo gondolo 2007-12-29
      implicit none
      integer dsi_trim
      character*200 filename
      include 'dsrdcom.h'
      include 'dsdirver.h'
      if (dofcode.eq.1) then
         filename=dsinstall(:dsi_trim(dsinstall))//
     &        'share/DarkSUSY/'//'dsdofGG_150.dat'
         call dsrdreaddof(filename)
         rdinit=1234
      else if (dofcode.eq.2) then 
         filename=dsinstall(:dsi_trim(dsinstall))
     &        //'share/DarkSUSY/'//'dsdofHP_A.dat'
         call dsrdreaddof(filename)
         rdinit=1234
      else if (dofcode.eq.3) then 
         filename=dsinstall(:dsi_trim(dsinstall))
     &        //'share/DarkSUSY/'//'dsdofHP_B.dat'
         call dsrdreaddof(filename)
         rdinit=1234
      else if (dofcode.eq.4) then 
         filename=dsinstall(:dsi_trim(dsinstall))
     &        //'share/DarkSUSY/'//'dsdofHP_B2.dat'
         call dsrdreaddof(filename)
         rdinit=1234
      else if (dofcode.eq.5) then 
         filename=dsinstall(:dsi_trim(dsinstall))
     &        //'share/DarkSUSY/'//'dsdofHP_B3.dat'
         call dsrdreaddof(filename)
         rdinit=1234
      else if (dofcode.eq.6) then 
         filename=dsinstall(:dsi_trim(dsinstall))
     &        //'share/DarkSUSY/'//'dsdofHP_C.dat'
         call dsrdreaddof(filename)
         rdinit=1234
      else 
         write (*,*) 'dsrdinit: invalid dofcode - ',
     &        'set it with call dsrdset(''dof'',dofcode)'
         stop
      endif
      return
      end
