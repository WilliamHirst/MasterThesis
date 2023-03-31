************************************************************
*** dssetfromslha sets one complex parameter from the SLHA commong
*** block slhadata, if the value is set. If not, it will take
*** the DarkSUSY defualt
*** Input: slhaitemr - real part of array item in SLHADATA
***        slhaitemi - imaginary part of array item in SLHADATA
***        option - if 0, disregard missing item in the SLHA file
***                 if 1, print a warning that 'text' is missing
***                 if 2, print an error that 'text' is missing
***        text - explanatory text for what variable we are working with
*** Output: dsparam - changed to slhaitem if valid
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-07-01
************************************************************

      subroutine dssetfromslhac(slhaitem,dsparam,
     &  option,text)
      implicit none
      complex*16 slhaitem,dsparam
      integer option,i
      character*(*) text

#include "SLHA.h"
#include "dsslha.h"
      
c      if (slhaitemi.ne.invalid) then
c         aux=slhaitemi
c      else
c         aux=0.0d0
c      endif

      if (slhaitem.ne.invalid) then
         dsparam=slhaitem
      else
         if (option.eq.1.and.prl.ge.1) then
            write(*,*) ('-',i=1,70)
            write(*,*)
     &    'DS Warning: The parameter ',text
            write(*,*)
     &    'was not present in the SLHA file, but DarkSUSY expected it.'
            write(*,*) 
     &    'Will disregard that it was missing and proceed anyway.'
            write(*,*) ('-',i=1,70)
         elseif (option.eq.2) then
            write(*,*) ('-',i=1,70)
            write(*,*)
     &    'DS Error: The parameter ',text
            write(*,*)
     &    'was not present in the SLHA file, but DarkSUSY needs it.'
            write(*,*) 
     &    'The program will stop.',
     &    ' Either you have to supply the variable,'
            write(*,*) 
     &     'or rewrite the slha interface in DarkSUSY to',
     &     ' handle this case.'
            write(*,*) ('-',i=1,70)
            stop
         endif
      endif

      return
      end

