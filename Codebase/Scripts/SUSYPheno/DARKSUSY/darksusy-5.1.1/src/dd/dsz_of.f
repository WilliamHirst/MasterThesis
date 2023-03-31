      function dsz_of(elem)
c...return the charge Z of the chemical element with symbol elem
c...input:
c...  elem - character string specifying the chemical symbol
c...author: paolo gondolo 2008-02-17
      implicit none
      integer dsz_of
      character*(*) elem
      integer n,i
      parameter (n=118)
      character*3 c
      character*3 ctable(n)

      data ctable/
     &'  h',' he',
     &' li',' be','  b','  c','  n','  o','  f',' ne',
     &' na',' mg',' al',' si','  p','  s',' cl',' ar',
     &'  k',' ca',
     &      ' sc',' ti','  v',' cr',' mn',' fe',' co',' ni',' cu',' zn',
     &            ' ga',' ge',' as',' se',' br',' kr',
     &' rb',' sr',
     &      '  y',' zr',' nb',' mo',' tc',' ru',' rh',' pd',' ag',' cd',
     &            ' in',' sn',' sb',' te','  i',' xe',
     &' cs',' ba',
     &      ' la',' ce',' pr',' nd',' pm',' sm',' eu',' gd',
     &            ' tb',' dy',' hu',' er',' tm',' yb',' lu',
     &      ' hf',' ta','  w',' re',' os',' ir',' pt',' au',' hg',
     &            ' ti',' pb',' bi',' po',' at',' rn',
     &' fr',' ra',
     &      ' ac',' th',' pa','  u',' np',' pu',' am',' cm',
     &            ' bk',' cf',' es',' fm',' md',' no',' lr',
     &      ' rf',' db',' sg',' bh',' hs',' mt',' ds',' rg','uub',
     &            'uut','uuq','uup','uuh','uus','uuo'
     &     /

c...find index
      write (c,'(a3)') elem
c      write (*,*) 'dsz_of (dbg): >',c,'<'
      call dslowcase(c)
c      write (*,*) 'dsz_of (dbg): >',c,'<'
      do i=1,n
         if (ctable(i).eq.c) then
            dsz_of=i
            return
         endif
      enddo

c...invalid chemical symbol
      write (*,*) 'dsz_of: non-existent chemical element ''',elem,''''
      dsz_of=0
      return
      end
