      subroutine dswhwarn(unit,hwarning)
c_______________________________________________________________________
c  write reasons for unphys<>0 to specified unit.
c  input:
c    unit - logical unit to write on (integer)
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsmssm.h'
      integer unit,hwarning
      if (unit.le.0) return
      if (hwarning.eq.0) then
        write(unit,*) 'No warnings issued from the Higgs routines.'
        return
      endif

      if (higloop.eq.3) then
        write(unit,*) 'Warnings from Carena-Espinosa-Quiros-Wagner:'
        if (btest(hwarning,0))
     &    write(unit,*) '  |mst^2-mst1^2|/(mst1^2+mst2^2)>0.5'
        if (btest(hwarning,1))
     &    write(unit,*) '  |msb^2-msb1^2|/(msb1^2+msb2^2)>0.5'
      elseif (higloop.eq.6) then
        write(unit,*) 'Warnings from FeynHiggsFast:'
        if (btest(hwarning,0))
     &    write(unit,*) 
     &  '  potential numerical problems at 1-loop (not serious)'
        if (btest(hwarning,1))
     &    write(unit,*) 
     &  '  potential numerical problems at 2-loop (not serious)'
        if (btest(hwarning,2))
     &    write(unit,*) 
     &      '  error with a not used H2 expression at 1-loop',
     &      ' (not important)'
        if (btest(hwarning,3))
     &    write(unit,*) 
     &      '  error with a not used H2 expression at 2-loop',
     &      ' (not important)'
        if (btest(hwarning,4))
     &    write(unit,*) 
     &      '  error with a not used H2 expression at 2-loop',
     &      ' (not important)'
        if (btest(hwarning,5))
     &    write(unit,*) '  1-loop Higgs sector not OK'
        if (btest(hwarning,6))
     &    write(unit,*) '  2-loop Higgs sector not OK'
        if (btest(hwarning,7))
     &    write(unit,*) '  stop or sbottom masses not OK'
      elseif (higloop.eq.5) then
         if (hwarning.ne.0) then
c          write(unit,*) '  Error in line ',hwarning,' of FeynHiggs.'
         endif
      else
          write(unit,*) '  unknown higloop = ',higloop
      endif

1000  format (1x,a:1x,f8.3:1x,f8.3)
1001  format (1x,a,1x,i4)
1002  format (1x,a,a)
      end
