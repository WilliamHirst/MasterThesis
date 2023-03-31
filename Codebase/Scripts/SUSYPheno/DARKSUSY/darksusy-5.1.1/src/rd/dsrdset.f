      subroutine dsrdset(key,value)
c...set parameters for relic density routines
c...  key - character string
c...  value - character string
c...parameters implemented:
c...  key='dof'
c...    value='help' - show a list of possible values
c...    value='default' - default
c...    value='1' - Gelmini-Gondolo 150MeV
c...    value='2' - Hindmarsch-Philipsen HP-A
c...    value='3' - Hindmarsch-Philipsen HP-B (default)
c...    value='4' - Hindmarsch-Philipsen HP-B2
c...    value='5' - Hindmarsch-Philipsen HP-B3
c...    value='6' - Hindmarsch-Philipsen HP-C
c...  key='help'
c...    any value - show a list of possible keys
c...author: paolo gondolo 2007-12-29
      implicit none
      include 'dsrdcom.h'
      character*(*) key,value

c...print list of keys
      if (key.eq.'help') then
         write (*,*) 'dsrdset: allowed keys are'
         write (*,*) '	dof	degrees of freedom in early universe'

c...default values
      else if (key.eq.'default') then
         if (dofcode.ne.1) rdinit=0
         dofcode=1

c...degrees of freedom
      else if (key.eq.'dof') then
         if (value.eq.'1') then 
            if (dofcode.ne.1) rdinit=0
            dofcode=1
         else if (value.eq.'2') then
            if (dofcode.ne.2) rdinit=0
            dofcode=2
         else if (value.eq.'3'.or.value.eq.'default') then
            if (dofcode.ne.3) rdinit=0
            dofcode=3
         else if (value.eq.'4') then
            if (dofcode.ne.4) rdinit=0
            dofcode=4
         else if (value.eq.'5') then
            if (dofcode.ne.5) rdinit=0
            dofcode=5
         else if (value.eq.'6') then
            if (dofcode.ne.6) rdinit=0
            dofcode=6
         else if (value.eq.'help') then
            write (*,*) 'dsrdset: key=dof can have values'
            write (*,*) '	''1''	Gelmini-Gondolo 150MeV (default)'
            write (*,*) '	''2''	Hindmarsch-Philipsen HP-A'
            write (*,*) '	''3''	Hindmarsch-Philipsen HP-B'
            write (*,*) '	''4''	Hindmarsch-Philipsen HP-B2'
            write (*,*) '	''5''	Hindmarsch-Philipsen HP-B3'
            write (*,*) '	''6''	Hindmarsch-Philipsen HP-C'
         else if (value.eq.'show') then
            write (*,*) 'dsrdset: key=dof has value ',dofcode
         else
            goto 1000
         endif

c...invalid choice
      else
         goto 1000
      endif

      return

 1000 continue
      write (*,*) 'dsrdset: unrecognized key and/or invalid value ',
     &     key,value
      stop
      end
