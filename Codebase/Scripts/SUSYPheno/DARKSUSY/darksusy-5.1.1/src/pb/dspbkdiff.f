**********************************************************************
*** diffusion constant in units of 10^27 cm^2 s^-1
*** n=1 value in the halo,  n=2 value in the gas disk
***
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dspbkdiff(rig,n)
      implicit none
      integer n
      real*8 rig
      include 'dspbcom.h'
      if(n.eq.1) then
c dspbkdiff equal to the diffusion constant in the halo
        dspbkdiff=pbkh*(1.d0+rig/pbrig0)**(pbdelta) ! JE change 060215
      else if(n.eq.2) then
c dspbkdiff equal to the diffusion constant in the gas
        dspbkdiff=pbkg*(1.d0+rig/pbrig0)**(pbdelta) ! JE change 060215
      else
        write(*,*) 'dspbkdiff called with wrong integer n=',n
      endif
      return
      end
