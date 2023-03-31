**********************************************************************
*** subroutine dsSLHAread reads SUSY Les Houches Accord files
*** with the help of SLHALIB by Thomas Hahn.
*** Input: file - filename
***        printlevel - 0 = print only errors
***                     1 = print errors and warnings
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-07-01
**********************************************************************
      subroutine dsSLHAread(file,printlevel)
      implicit none
c...Note: These include files need to be included with the precompiler
c...directive for include and not the fortran one. Hence, this file
c...has to end on .F instead of .f as well.

#include "SLHA.h"

#include "PDG.h"

      include 'dsidtag.h'
#include "dsslha.h"

      integer error,abort,i,printlevel
      character*(*) file
      logical first
      data first/.true./
      save first

      prl=printlevel ! transfer to common block

      if (first) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Info: You are now using the new DarkSUSY SLHA interface.'
         write(*,*) 
     &   'Please be advised that this interface is intended for SLHA2'
         write(*,*) 
     &   'files and that all possible options for how an SLHA2'
         write(*,*) 
     &   'file can look like are not in yet. If you are missing'
         write(*,*) 
     &   'a specific feature, please let us know and we will try'
         write(*,*) 
     &   'to implement it. With the current status of the code'
         write(*,*)
     &   'we will issue warnings when things are not exactly as'
         write(*,*)
     &   'expected, but probably OK (and proceed with the SLHA file).'
         write(*,*)
     &   'When options are found that are not yet implemented,',
     &   ' the program will stop.'
         write(*,*)
     &   'Contact Joakim Edsjo (edsjo@physto.se)',
     &   ' if you have questions.'
         write(*,*) ('-',i=1,70)
         first=.false.
      endif



c...Clear slhadata
      call SLHAClear(slhadata)

c...Read it in
      abort=0 ! write warning when reading, but continue
      call SLHARead(error,slhadata,file,abort)
      if (error .ne. 0) then
         write(*,*) 'DS ERROR: Reading of SLHA file failed.'
         write(*,*) 'Error code from Hahn''s SLHALIB: ',error
         stop
      endif

      idtag=file ! We wan't to get this from the file eventually

c...Transfer stuff from the SLHA data structure to DarkSUSY
      call dsfromslha

c...After this, we are at the same state as for an MSSM model after
c...a call to dsgive_model (or similar) and dsmodel_setup.
c...We should thus NOT call dsmodelsetup here, because then we
c...will override what we have read from the SLHA file.


      return
      end

      
