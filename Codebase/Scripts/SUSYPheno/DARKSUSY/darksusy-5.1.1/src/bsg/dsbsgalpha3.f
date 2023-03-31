      function dsbsgalpha3(m)

***********************************************************************
* The coupling constant alpha_3 evaluated at the scale m              *
* using nf effective quark flavours (usually taken to be nf=5)        *
* Uses eq. (42) of Ciuchini et al. hep-ph/9710335                     *
* for the calculation of b --> s gamma                                *
* Note: This routines is strictly speaking only valid for mass scales *
* between mb and mt where nf=5 should be used.                        *
* author:Mia Schelke, schelke@physto.se, 2003-04-03                   *
***********************************************************************

      implicit none
      include 'dsmssm.h'
      real*8 m,mt,altmp
      real*8 al3mz
      real*8 dsbsgalpha3,dsbsgalpha3int

c     set the value of alpha_3(M_Z) 
c     use the value from DarkSUSY

      al3mz=alph3mz
      
c...Pick t mass at scale mt from DarkSUSY
      mt=mtmt

      if (m.le.mtmt) then   ! 5 active quarks

        dsbsgalpha3=dsbsgalpha3int(al3mz,mass(kz),m,5)

      else            ! 6 active quarks
c      elseif (m.gt.mtmt.and.m.le.1.5d0*mtmt) then ! 6 active quarks 
c     The line above should be used if we want to use the alternative 
c     6 quarks and one light squark at higher energies

        altmp=dsbsgalpha3int(al3mz,mass(kz),mtmt,5)
        dsbsgalpha3=dsbsgalpha3int(altmp,mtmt,m,6)

c      elseif (m.gt.1.5d0*mtmt) then ! 6 active quarks and one squark

c        altmp=dsbsgalpha3int(al3mz,mass(kz),mtmt,5)
c        altmp=dsbsgalpha3int(altmp,mtmt,1.5d0*mtmt,6)
c        dsbsgalpha3=dsbsgalpha3int(altmp,1.5d0*mtmt,m,7)

c      else

c        write(*,*) 'ERROR in dsbsgalpha3: scale of out range: ',m
c        stop

      endif

      return
      end


