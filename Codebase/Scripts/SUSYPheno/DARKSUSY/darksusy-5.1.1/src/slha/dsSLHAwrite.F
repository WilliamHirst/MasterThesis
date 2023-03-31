**********************************************************************
*** subroutine dsSLHAwrite writes SUSY Les Houches Accord files
*** with the help of SLHALIB by Thomas Hahn.
*** Input: file = filename of SLHA2 file to write to
***        opt = 1 - an SLHA2 file with full 6x6 squark and slepton
***                  mixings and 3x3 trilinear couplings is written.
***                  For minimal flavour violation
***                  this is a bit of an overkill, but more consistent
***                  with how things are stored internally in DarkSUSY
***              2 - a minimal flavour violation SLHA2 file is written
***                  with only third generation mixings and trilinear
***                  couplings. This option will in principle be an
***                  approximation to the internal setup in DarkSUSY,
***                  but most likely good enough for most scenarios.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-06-24
**********************************************************************
      subroutine dsSLHAwrite(file,opt)
      implicit none
c...Note: These include files need to be included with the precompiler
c...directive for include and not the fortran one. Hence, this file
c...has to end on .F instead of .f as well.

#include "SLHA.h"

#include "PDG.h"

      include 'dsidtag.h'
#include "dsslha.h"

      integer error,i,opt
      character*(*) file
      logical first
      data first/.true./
      save first

      if (first) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Info: You are now using the new DarkSUSY SLHA interface.'
         write(*,*) 
     &   'Please be advised that this interface will write SLHA2'
         write(*,*) 
     &   'files and that all possible options for how an SLHA2'
         write(*,*) 
     &   'file can look like are not in yet.'
         write(*,*) 
     &   'In particular, we have chosen some midway between dumping'
         write(*,*) 
     &   'all the details to the SLHA file instead of just the inputs.'
         write(*,*)
     &   'This choice has been made to make our files compatible with'
         write(*,*)
     &   'as many other codes as possible that read SLHA files.'
         write(*,*)
     &   'If these files do not work for you, please let us know'
         write(*,*)
     &   'and we will take a look at them.'
         write(*,*)
     &   'Contact Joakim Edsjo (edsjo@physto.se)',
     &   ' if you have questions.'
         write(*,*) ('-',i=1,70)
         first=.false.
      endif
      

c...Transfer DarkSUSY stuff to the SLHA data structure
      call dstoslha(opt)

c...Write it out
      call SLHAWrite(error,slhadata,file)

      if (error .ne. 0) then
         write(*,*) 'DS ERROR: Writing of SLHA file failed.'
         write(*,*) 'Error code from Hahn''s SLHALIB: ',error
         stop
      endif

      return
      end

      
