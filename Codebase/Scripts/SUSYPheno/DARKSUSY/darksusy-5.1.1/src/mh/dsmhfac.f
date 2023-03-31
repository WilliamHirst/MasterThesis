***********************************************************************
*** returns the factorial of an integer, based on NUMERICAL RECIPES
***
*** torsten bringmann (troms@physto.se), 2008-12-03
*** Modified: Joakim Edsjo, edsjo@fysik.su.se, 2011-05-25
***********************************************************************

      REAL*8 FUNCTION dsmhfac(n)
      implicit none
      INTEGER n
      INTEGER j,ntop
      REAL a(33) ! Table to be filled in only as required.
      real*8 dsmhgamma
      SAVE ntop,a
      DATA ntop,a(1)/0,1./ !Table initialized with 0! only.

      if (n.lt.0) then
         write (*,*) 'negative factorial in dsmhfac'
         stop
      else if (n.le.ntop) then ! Already in table.
        dsmhfac=a(n+1)
      else if (n.le.32) then ! Fill in table up to desired value.
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
  11    continue
        ntop=n
        dsmhfac=a(n+1)
      else ! Larger value than size of table is required. 
        dsmhfac=dsmhgamma(n+1d0)
      endif
      return
      END
