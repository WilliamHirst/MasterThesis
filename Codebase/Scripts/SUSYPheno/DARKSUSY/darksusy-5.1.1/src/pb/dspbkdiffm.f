**********************************************************************
*** diffusion constant in units of 10^27 cm^2 s^-1
*** n=1 value in the halo,  n=2 value in the gas disk
*** form needed for routine dspbtd15beum
***
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dspbkdiffm(beta,rig,n)
      implicit none
      integer n
      real*8 beta,rig
      include 'dspbcom.h'
      if(n.eq.1) then
c dspbkdiff equal to the diffusion constant in the halo
        if(rig.lt.pbrig0.and.pbc(1:6).ne.'salati') then ! M.Gustafsson 2006-01-18
          dspbkdiffm=pbkh
        else  
          dspbkdiffm=pbkh*(rig/pbrig0)**(pbdelta)
        endif  
      else if(n.eq.2) then
c dspbkdiff equal to the diffusion constant in the gas
        if(rig.lt.pbrig0.and.pbc(1:6).ne.'salati') then ! M.Gustafsson 2006-01-18
          dspbkdiffm=pbkg
        else 
          dspbkdiffm=pbkg*(rig/pbrig0)**(pbdelta)
        endif
      else
        write(*,*) 'dspbkdiffm called with wrong integer n=',n
      endif
      dspbkdiffm=dspbkdiffm*beta
      return
      end

