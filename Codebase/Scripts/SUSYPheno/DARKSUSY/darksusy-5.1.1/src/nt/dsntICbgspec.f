***********************************************************************
*** dsntICbgspec provides the interpolated probability distribution
*** function for the number of hit DOMs due to background events. 
***
*** Input:	nchan		number of hit DOMs
*** Output:                     pdf (chan^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Apr 24, 2011
***********************************************************************

      real*8 function dsntICbgspec(nchan)

      implicit none
      include 'dsntIC.h'

      integer nchan, nchan_index

      if (nchan .lt. nchan_min .or. nchan .gt. nchan_max) then
        write(*,*) 'Error in dsntICbgspec: nchan outside tabulated'
        write(*,*) 'range, nchan=',nchan,'.  Quitting...'
        stop
      endif
  
      nchan_index = nchan - nchan_min + 1 - nchan_hist2BGoffset
      if (BGnchandist_nchan(nchan_index) .ne. nchan) then
        write(*,*) BGnchandist_nchan(nchan_index), nchan, nchan_index
        stop 'Something is wrong with nchan_index in dsntICbgspec.'
      endif

      dsntICbgspec = BGnchandist_prob(nchan_index)

      end function dsntICbgspec
