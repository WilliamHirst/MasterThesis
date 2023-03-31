****************************************************************
*** spherically symmetric dark matter halo density profile   ***
***                                                          ***
*** Input:  radialdist = radial coordinate in kpc            ***
***         in a spherically symmetric coordinate system     ***
***         centered in the Galactic Center                  ***
***         if radialdist lower than the cut radius rhcut,   ***
***         radialdist is shifted to rhcut                   ***
*** Output: dshmsphrho = density in GeV/cm^3                 ***
***   e.g.: local halo density = rho0 = dshmsphrho(r_0)      ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
*** Modified: 2009-05-07, 2009-08-08 Pat Scott               ***
****************************************************************

      real*8 function dshmsphrho(radialdist)
      implicit none
      include 'dshmcom.h'
      real*8 radialdist,rr,dshmabgrho,dshmburrho,dshmn03rho
     &  ,dshmnumrho,dshmeinastorho!,dshmPLUMrho
c
      if(radialdist.lt.rhcut) then !.and. (halotype .ne. 'PLUM' .or. rhcut .lt. Rh_PLUM) ) then
        rr=rhcut
      else
        rr=radialdist
      endif  
c
      if(halotype.eq.'albega') then
c alpha-beta-gamma parametrization
        dshmsphrho=dshmabgrho(rr)
      elseif(halotype.eq.'burkert') then
c burkert parametrization
        dshmsphrho=dshmburrho(rr) 
      elseif(halotype.eq.'n03') then
c navarro et al.(2003) parametrization
        dshmsphrho=dshmn03rho(rr) 
      elseif(halotype.eq.'einasto') then
c einasto parametrization
        dshmsphrho=dshmeinastorho(rr) 
      elseif(halotype.eq.'numerical') then
c profile read from disc
        dshmsphrho=dshmnumrho(rr) 
      !elseif(halotype.eq.'PLUM') then
c PLUM
      !  dshmsphrho=dshmPLUMrho(rr) 
      else
c no other option set yet  
        write(*,*) 'dshmsphrho called wrong option halotype'
        write(*,*) 'halotype = ',halotype
        write(*,*) 'program stopped'
        stop
      endif
      return
      end
