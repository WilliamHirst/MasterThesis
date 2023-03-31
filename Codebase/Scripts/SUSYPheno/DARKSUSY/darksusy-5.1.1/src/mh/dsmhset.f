***********************************************************************
*** The routine dsmhset has to be called once before any microhalo (mh)
*** (mh) calculations; it makes all necessary initializations and sets 
*** some default settings.
*** 
***  c - character string specifying choice to be made 
***  possible values for c:
***  'susy' or 'default' - neutralino DM
***  'mUED'              - (mUED) Kaluza-Klein DM 
***  'user'              - user-defined DM model
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      subroutine dsmhset(c)
      implicit none

      include 'dsmhcom.h'
      include 'dsdirver.h'

      character*(*) c
      character*200 doffile
      integer dsi_trim

      mheps=1D-2   ! this sets the accuracy for the numerical
                   ! routines that determine Tkd

      doffile=dsinstall(:dsi_trim(dsinstall))
     &        //'share/DarkSUSY/'//'dsmhdof.dat'

      call dsmhreaddof(doffile)  ! reads in tabulated expressions
                                       ! for eff. degrees of freedom

      if (c.eq.'susy'.or.c.eq.'default') then

         mhtype=1
c         write(*,*) 'SUSY'

      elseif (c.eq.'mUED') then

         mhtype=2
         write(*,*) 'please call dsmhUED to initialize Kaluza-Klein DM'

      elseif (c.eq.'user') then

         mhtype=3
         write(*,*) 'user defined model'

c -> extend files dsmhboltz_init.f, dsmhm2.f and dsmhm2simp.f accordingly
c -> add file dsmhXXX.f to initialize model (for an example, see dsmhUED)


      else
         write(*,*) 'ERROR in dsmhset -- unknown option: ',c
         write(*,*) 'Stopping...'
         stop
      endif

      return

      end





        

