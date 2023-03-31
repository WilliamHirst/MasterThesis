      program dstest
c
c     This program tests some DarkSUSY routines
c
c     In addition to the DarkSUSY programs and data files, this test
c     program chooses a few mSUGRA models randomly and uses Isasugra to
c     calculate the low energy models. These are then fed into DarkSUSY
c     and the relic density is calculated.
c     If you want to calculate rates and stuff, use the calls as
c     outlined in dstest.f

c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      real*8 oh2,xf,dsrdomega                            ! relic density
      character*80 message,scr
      logical first
      data first/.true./
      integer ii, jj, idum,nfc,iwar,ierr,unphys,hwarning
      real*8 m0min,m0max,m5min,m5max,a0min,a0max,tbmin,tbmax,
     &  m0,mhf,a0,tgb,sgnmu
      real*8 dsrndlin,dsrndsgn,dsabsq
c
c     Here we include the file dsmssm.h which contains common block
c     definitions for particle masses, susy parameters, and results.
c     We also include dsidtag.h which gives the possibility to tag
c     every model with a 12 character id tag (variable idtag). This
c     id tag is printed in case of errors for that model.
c
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsidtag.h'

      idum = -363469   ! seed for random number generator

c
c     This call initializes the DarkSUSY package. This call should
c     be the first call in any program using DarkSUSY.
c
      call dsinit
      call aldata   ! to make sure isajet common block datas are set
                    ! on most compilers, this is not needed

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
      prtlevel=1
c
c     This is an example of redefining DarkSUSY parameters reading them
c     from an input file, in this case the file dsdefault.dat which
c     contains the default values (this is done for test purposes only,
c     because the default values are already included in the code). The
c     user can copy the file dsdefault.dat to a new file, change the
c     values of the variables, and read in the new file.  Not all
c     variables need to be assigned in the new file, but only those
c     whose value differs from the default value. Note: dsdefault.dat
c     resets prtlevel to 0, ie no messages.
c     Note: the dsdefault.dat is not upto date, but a call like the
c     following can be used to reset parameters to your own choices.
c
c      open (unit=10,file='dsdefault.dat',status='old')
c      call dsreadpar(10)
c      close (10)


c     Now we are ready to scan the supersymmetric parameter space. 
c     We here read some models from dstest-isasugra.mod, which
c     are a few of the Benchmark points in hep-ph/0306219
c
      

      open(unit=10,file='dstest-isasugra.mod',status='old',
     &  form='formatted')
      read(10,'(a)') scr ! read header lines
      read(10,'(a)') scr ! read header lines
      read(10,'(a)') scr ! read header lines


c     these are minimal SUGRA
 20   read(10,*,end=100) idtag,m0,mhf,a0,sgnmu,tgb
c 1001 format(1x,A12,3(1x,E14.8),1x,e7.1,5(1x,e14.8))

      write(*,*) ' '
      write(*,*) 'Model: ',idtag
      write(*,*) 'Parameters: ',m0,mhf,a0,sgnmu,tgb

c      Now define the model, see dsgive_model_isasugra how we transfer
c      these varibles to DarkSUSY common blocks
      call dsgive_model_isasugra(m0,mhf,a0,sgnmu,tgb)
         

c      Now set up the model by a call to dssusy_isasugra. dssusy_isasugra
c      calls isasugra and transfers the results to DarkSUSY common blocks.
c      DarkSUSY is then set-up for this model. Check the flags unphys
c      and hwarning which are non-zero if the model is not good.
      call dssusy_isasugra(unphys,hwarning)
       
      write(*,*) '   hwarning = ',hwarning
      write(*,*) '     unphys = ',unphys
      if (hwarning.eq.0.and.unphys.eq.0) then  ! Model OK

c
c     We now have all the masses and couplings calculated. Let's print
c     some out to see what they are. dsabsq is just the square of the
c     absolute value of it's complex*16 argument.

        write(*,*) '  Neutralino mass = ',mass(kn(1))
        write(*,*) '  Gaugino fraction = ',
     &    dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
        write(*,*) '  H1 mass =  ',mass(kh1),width(kh1)
        write(*,*) '  H2 mass =  ',mass(kh2),width(kh2)
        write(*,*) '  H3 mass =  ',mass(kh3),width(kh3)
        write(*,*) '  H+- mass = ',mass(khc),width(khc)

c...Now calculate the relic density with all coannihilations
        write(*,*) '   Calculating omega h^2, please be patient...'
        oh2=dsrdomega(1,1,xf,ierr,iwar,nfc)
        write(*,*) '   Oh2 = ',oh2

c...Go on calculating rates etc, as shown in dstest.f if we want to
      endif

      goto 20

 100  continue

      close(10)
      stop
      end
