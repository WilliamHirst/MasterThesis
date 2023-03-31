      program dsslha2slha
c
c     This program reads in an SLHA2 file (mSUGRA or low-energy MSSM)
c     and writes out a new SLHA file (in 6x6 FV setup, or 2x2 MFV setup).
c     If the input SLHA file is an mSUGRA input file, isasugra will be
c     run for this model.
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      integer opt
      character slhain*128,slhaout*128

c
c     Here we include the file dsmssm.h which contains common block
c     definitions for particle masses, susy parameters, and results.
c     We also include dsidtag.h which gives the possibility to tag
c     every model with a 12 character id tag (variable idtag). This
c     id tag is printed in case of errors for that model.
c
      include 'dsidtag.h'

c
c     This call initializes the DarkSUSY package. This call should
c     be the first call in any program using DarkSUSY. This call initializes
c     some global variables and calls various other modules to initialize
c     them with their default values. Check out src/ini/dsinit.f if you
c     want to see what it does.
c
      call dsinit

c
c     The amount of output onto standard output can be controlled by the
c     variable prtlevel. Setting prtlevel=0 suppresses all messages,
c     except fatal errors that terminate the program. Setting prtlevel=1
c     displays informational messages. Setting prtlevel to higher values
c     displays more and more information. In this example, we set
c     prtlevel=1 to have a message from dsreadpar. We reset it to 0
c     later. Had we set prtlevel=2 here, we would also have messages for
c     each variable which is (re)set in dsreadpar. Caution: if prtlevel
c     is reset in dsreadpar, the new value applies to subsequent
c     messages.
c
      prtlevel=0
c


c...Now read in the first SLHA file
      write(*,*) 'Enter the filename for your input SLHA2 file: '
      read(5,'(A)') slhain
      call dsSLHAread(slhain,0) ! 0=no warnings, 1=print warnings

      write(*,*) ' '
      write(*,*) '***** MODEL: ',idtag,' *****'

      write(*,*) 
     &  'What kind of SLHA file do you want to write?'
      write(*,*) '  1 = with full 6x6 sfermion mixing'
      write(*,*) '  2 = with minimal flavour violation',
     &  ' (2x2 sfermion mixing)'
      read(*,*) opt

      write(*,*) 'Give SLHA2 file name:'
      read(5,'(A)') slhaout
      call dsSLHAwrite(slhaout,opt) 

      write (*,*) 'Done.'
      stop
      end
