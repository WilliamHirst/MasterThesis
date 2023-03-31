**********************************************************************
*** function dsepdiff calculates the differential flux of
*** positrons for the energy egev as a result of
*** neutralino annihilation in the halo.
*** input: egev - positron energy in gev
***        how = 1 - dsepsigvdnde is used directly
***              2 - dsepsigvdnde is tabulated on first call, and then
***                  interpolated. (default)
***       dhow = 2 - diffusion model is tabulated on first call, and then
***                  interpolated
***              3 - as 2, but also write the table to disk at the
***                  first call
***              4 - read table from disk on first call, and use that for
***                  subsequent calls. If the file does not exist, it will
***                  be created (as in 3). (default)
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
***         joakim edsjo, edsjo@physto.se
*** date: jun-02-98
*** modified: 99-07-02 paolo gondolo : order of calls to hrsetup
*** modified: 01-10-19 add moskalenko + strong option (eab)
*** modified: 04-01-27 J. Edsjo: added file handling (dhow=3,4)
**********************************************************************

      real*8 function dsepdiff(egev,how,dhow)
      implicit none

      include 'dsidtag.h'
      include 'dshacom.h'
      include 'dsepcom.h'
      include 'dshmcom.h'
      include 'dsprep.h'
      include 'dsgalpropcom.h'
      include 'dsdirver.h'

c----------------------------------------------------------------------

      real*8 egev,sigvline,dsepsigvdnde,dsepspec,dsepipol
      real*8 dsepeeuncut,dsepeecut,dsepvvuncut,dsepvvcut,
     &  dsepwuncut,dsepwcut,dsepktdiff,dsepmsdiff,dsepgalpropdiff
      integer how,dhow,i
      external dsepsigvdnde,dsepipol,dsepeeuncut,dsepeecut,
     &  dsepvvuncut,dsepvvcut,dsepwuncut,dsepwcut
      character*128 epfile,scr

c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dsepdiff: dshasetup must be called',
     &    ' before any halo rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif
 
      sigvline=habr(15)*hasv  ! cm^3 s^-1

      if (dsepdiffmethod.eq.2) then ! kamionkowski & turner propagation
        dsepdiff=dsepktdiff(egev)
        return
      endif
      
      if (dsepdiffmethod.eq.3) then ! moskalenko & strong propagation
         dsepdiff=dsepmsdiff(egev)
         return
      endif
      
      if (dsepdiffmethod.eq.4) then ! galprop propagation
         if (.not.gpgfread) then
            write (6,*)
     +           'You must call dsgalprop_gettable(0) to use GALPROP'
            stop
         endif
         dsepdiff=dsepgalpropdiff(egev)
         return
      endif
      
      if (newmodelep.and.how.eq.2) then ! new model
         call dseptab(0.1d0,200)
         newmodelep=.false.
      endif


c... Use Baltz and Edsjo propagation

c...Make tables according to dhow

      if (dhow.lt.2.or.dhow.gt.4) then
        write(*,*) 'ERROR in dsepdiff: Wrong value for dhow: ',dhow
        stop
      endif

      if (epload) then

c...Generate file name
        epfile=dsinstall
        call dscharadd(epfile,'share/DarkSUSY/eptab-')
        call dscharadd(epfile,epc)
        call dscharadd(epfile,'-')
        call dscharadd(epfile,haloid)
        call dscharadd(epfile,'.dat')

        if (dhow.eq.2) then

          write(*,*) 'dsepdiff: starting tabulation...'
          call dsepmake_tables
          call dsepmake_tables2
          write(*,*) '  ...done'

        elseif (dhow.eq.3) then

          write(*,*) 'dsepdiff: starting tabulation...'
          call dsepmake_tables
          call dsepmake_tables2
          write(*,*) '  ...done'
          write(*,*) 'dsepdiff: writing table to file ',epfile
          open(unit=13,file=epfile,status='unknown',
     &      form='formatted')
          write(13,1001) dsversion
          write(13,1002) 
 1001     format('# Made with DarkSUSY version ',A)
 1002     format('#','...table(i)...')
          do i=1,22001
            write(13,1003) table(i)
          enddo
 1003     format(1x,e14.8)
          close(13)
          write(*,*) '  ...done'

        else  ! dhow.eq.4

          write(*,*) 'dsepdiff: reading table from file ',epfile
          open(unit=13,file=epfile,status='unknown',
     &      form='formatted')
          read(13,'(a)',end=100) scr
          read(13,'(a)',end=100) scr
          do i=1,22001
            read(13,1003,end=100) table(i)
          enddo
          close(13)
          call dsepmake_tables2  ! make the small tables
          write(*,*) '  ...done'
          goto 110

 100      close(13)
          write(*,*) 'WARNING in dsepdiff: the file ',epfile
          write(*,*) 'does either not exist or is incomplete.'
          write(*,*) 'Will create a new data file instead.'

          write(*,*) 'dsepdiff: starting tabulation...'
          call dsepmake_tables
          call dsepmake_tables2
          write(*,*) '  ...done'
          write(*,*) 'dsepdiff: writing table to file ',epfile
          open(unit=13,file=epfile,status='unknown',
     &      form='formatted')
          write(13,1001) dsversion
          write(13,1002) 
          do i=1,22001
            write(13,1003) table(i)
          enddo
          close(13)
          write(*,*) '  ...done'

 110      continue
        endif

        epload=.false.

      endif

c...Now  calculate the fluxes

      if (how.eq.1.and.dsepdiffmethod.eq.0) then
         dsepdiff=dsepspec(egev,hamwimp,sigvline,dsepsigvdnde,
     &     dsepeeuncut,dsepvvuncut,dsepwuncut)
      else if (how.eq.1.and.dsepdiffmethod.eq.1) then
         dsepdiff=dsepspec(egev,hamwimp,sigvline,dsepsigvdnde,
     &     dsepeecut,dsepvvcut,dsepwcut)
      else if (how.eq.2.and.dsepdiffmethod.eq.0) then
         dsepdiff=dsepspec(egev,hamwimp,sigvline,dsepipol,dsepeeuncut,
     &     dsepvvuncut,dsepwuncut)
      else
         dsepdiff=dsepspec(egev,hamwimp,sigvline,dsepipol,dsepeecut,
     &     dsepvvcut,dsepwcut)
      endif

c...Now rescale
      dsepdiff=(rhox/rho0)**2*dsepdiff

      return
      end

