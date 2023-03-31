      program dstest
c
c     This program tests many DarkSUSY routines
c
c     In addition to the DarkSUSY programs and data files, this test
c     program needs the input file dstest.mod, provided with the
c     DarkSUSY distribution. The test program outputs one file for
c     checking purposes (in an actual application, the user should
c     decide what information to output and in which format):
c     dstest1.tmp is output for the first model. You should compare it
c     with dstest.tmp provided with the DarkSUSY distribution.
c     You will also find dstest.output with the distribution. This file is
c     the output (to the terminal) when running the test program. You
c     should compare this with your output as well. All numbers should
c     to within the numerical accuracy. Some numbers, like the Z gamma
c     cross sections may wary a few % when the cross sections are low
c     due to numerical differences.

c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      integer i

c     Here we include dsidtag.h which gives the possibility to tag
c     every model with a 12 character id tag (variable idtag). This
c     id tag is printed in case of errors for that model.
c
      include 'dsidtag.h'
      include 'dsgalpropcom.h'

c
c     This call initializes the DarkSUSY package. This call should
c     be the first call in any program using DarkSUSY. This call initializes
c     some global variables and calls various other modules to initialize
c     them with their default values. Check out src/ini/dsinit.f if you
c     want to see what it does.
c

      read (5,*) i
      call dsinit

      call dsgalpropset('plain')
      
      call dsgalprop_maketable_one(i)
      end
