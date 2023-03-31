      subroutine dsgalpropset(c)
c...set parameters for positron routines
c...  c - character string specifying choice to be made
c...author: joakim edsjo, 2006-02-21
      implicit none
      include 'dsgalpropcom.h'
      character*(*) c

c...
      if (c.eq.'default') then
         gpnumdecade=10
         gpnumin=5*gpnumdecade+1
         gpnumout=9*gpnumdecade+1
         gpgfread=.false.
         gpmodtag='RACC'
      else if (c.eq.'plain') then
         gpnumdecade=10
         gpnumin=5*gpnumdecade+1
         gpnumout=9*gpnumdecade+1
         gpgfread=.false.
         gpmodtag='DIFF'
      else if (c.eq.'test') then
         gpnumdecade=10
         gpnumin=5*gpnumdecade+1
         gpnumout=9*gpnumdecade+1
         gpgfread=.false.
         gpmodtag='TEST'
c...  invalid choice
      else
         write (*,*) 'dsgalpropset: unrecognized option ',c
         stop
      endif

      return
      end
