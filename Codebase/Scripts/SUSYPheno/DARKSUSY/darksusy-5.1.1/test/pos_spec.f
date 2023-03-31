      program pos_spec

c     writes out positron spectra with and without FSR
c     Torsten Bringmann, 2007-10-20


      implicit none
      real*8 flux1,flux2,flux3 ! fluxes in e^+ cm^-2 s^-1 sr^-1 GeV^-1
      real*8 epos                              ! positron energy
      real*8 n1,n2,n3,n4,n5 ! number of positrons
      integer unphys,hwarning,iend,ierr,istat
      real*8 dsabsq,sv,dssigmav, dshaloyield,mygf,m0
      real*8 dsepspecm
      integer COUNT, COUNT2, COUNT_RATE, COUNT_MAX
      integer i,test
      character filenametrunk*15,modelfilename*50
      integer modelflag, skip
      real*8 dshayield,alpha2

      real*8 Cllww1, Crrww1, Clrww1, Clrww51,
     &       Cllww2, Crrww2, Clrww2, Clrww52


      include 'dsidtag.h'
      include 'dshrcom.h'
      include 'dsprep.h'
      include 'dsibcom.h'

      call dsinit
      call dshmset('isosm')
      write(*,*)
      prtlevel=0  ! suppress DS messages

      modelflag=2      ! 1-MSSM, 2-mSUGRA

      if (modelflag.eq.2) then
        filenametrunk='ted/ASNm03'    ! BM1
c        filenametrunk='ted/ASNm24'    ! BM2
c        filenametrunk='ted/FMCp03'    ! BM3
c        filenametrunk='ted/ASNp01'    ! BM4
        modelfilename='./'//
     -      filenametrunk(1:index(filenametrunk,' ')-1) //
     -      '.lemod'     ! low energy model file
        skip=748      ! BM1; skip the first <skip> entries in the .lemodfile
c        skip=543      ! BM2; skip the first <skip> entries in the .lemodfile
c        skip=765      ! BM3; skip the first <skip> entries in the .lemodfile
c        skip=2560     ! BM4; skip the first <skip> entries in the .lemodfile
      endif

      if (modelflag.eq.1) then
        filenametrunk='dstest_tb'
        modelfilename='./' //
     -      filenametrunk(1:index(filenametrunk,' ')-1) //
     -      '.mod'     ! test model file
        open (unit=10,file=modelfilename)
        read (10,*) ! this skips the header
        read (10,*)
c       read (10,*)
      endif

c  output file with spectra
      open (unit=30,file='pos_test.spec')    


      iend=0

      call SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)

c     begin loop to read in models
 1000 if (modelflag.eq.1) call read_model(10,0,iend,ierr)
      if (modelflag.eq.2) then
 1050   call read_mod(modelfilename,ierr)
        skip=skip-1
        if (skip.gt.0) goto 1050
        if (ierr.ne.0) goto 2000
      endif

      if (ierr.ne.0) then
         write (*,*) 'Error while reading model file'
         write (*,*) 'The positron test cannot continue'
         stop
      else if (iend.eq.1) then
         goto 2000 ! exit program
      else
c         call dswrite(0,0,'Model parameters read from file model file')
      endif

      if (modelflag.eq.1) call dssusy(unphys,hwarning)
      if (modelflag.eq.2) then
        call set_mod
        call call_dssusy(unphys,hwarning,1)
      endif   

c     make sure we computed the branching ratios
      sv=dssigmav(0)
      mygf=dsabsq(neunmx(1,1))+
     -     dsabsq(neunmx(1,2))
      m0=mass(kn(1))

c...set up couplings

      alpha2=alphem/s2thw
      Cllww1  = abs(gl(kw,kn(1),kcha(1)))**2/alpha2
      Crrww1  = abs(gr(kw,kn(1),kcha(1)))**2/alpha2
      Clrww1  = imag(gl(kw,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1))))
     - /alpha2
      Clrww51 = dble(gl(kw,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1))))
     - /alpha2
      Cllww2  = abs(gl(kw,kn(1),kcha(2)))**2/alpha2
      Crrww2  = abs(gr(kw,kn(1),kcha(2)))**2/alpha2
      Clrww2  = imag(gl(kw,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2))))
     - /alpha2
      Clrww52 = dble(gl(kw,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2))))
     - /alpha2


      write(*,*) 'Model: ',idtag
      write(*,*) '----------------------------'
      write(*,*) 'Neutralino mass:',m0
c      write(*,*) 'Stop masses:',mass(kst(1)),mass(kst(2))
c      write(*,*) 'Stau masses:',mass(kstau(1)),mass(kstau(2))
c      write(*,*) 'Selectron masses:',mass(kse(1)),mass(kse(2))
c      write(*,*) 'Smuon masses:',mass(ksmu(1)),mass(ksmu(2))
      write(*,*) 'Z_g/(1-Z_g): ', mygf/(1d0-mygf)
      write(*,*) 'total sigmav: ',sv
      write(*,*) 'W+W-: ',dssigmav(13)/sv  
c      write(*,*) 'light leptons: ',
c     -            (dssigmav(15)+dssigmav(17))/sv
c      write(*,*) 'tau: ',dssigmav(19)/sv
c      write(*,*) 'top: ',dssigmav(24)/sv 
c      write(*,*) 'bottom: ',dssigmav(25)/sv
c      write(*,*) 
c      write(*,*)  Cllww1, Crrww1, Clrww1, Clrww51,
c     &       Cllww2, Crrww2, Clrww2, Clrww52
      write(*,*)

      epos=m0*0.5d0

      call dsIBset('medium')
c      IBflag(7)=0
c      IBflag(8)=0
c      IBflag(9)=0
c      IBflag(10)=0
c      IBflag(11)=0
c      IBflag(12)=0


      IBflag(1)=0
      IBflag(4)=0
      IBflag(5)=0
      IBflag(6)=0
      n1=dshaloyield(epos,51,istat)            ! 151:differential, 51:integrated
      write(*,*) 'Number of positrons above ',
     -           '0.8 m_chi (without IB): ',n1

      IBflag(4)=1
      n2=dshaloyield(epos,51,istat)
      write(*,*) 'Number of positrons above ',
     -           '0.8 m_chi (only from eeg): ',n2-n1

      IBflag(5)=1
      n3=dshaloyield(epos,51,istat)
      write(*,*) 'Number of positrons above ',
     -           '0.8 m_chi (only from mumug): ',n3-n2

      IBflag(6)=1
      n4=dshaloyield(epos,51,istat)
      write(*,*) 'Number of positrons above ',
     -           '0.8 m_chi (only from tautaug): ',n4-n3

      IBflag(1)=1
      n5=dshaloyield(epos,51,istat)
      write(*,*) 'Number of positrons above ',
     -           '0.8 m_chi (only from WWg): ',n5-n4
      write(*,*) 'Number of positrons above ',
     -           '0.8 m_chi (total): ',n5

      write(*,*)

      IBflag(5)=0
      IBflag(6)=0
      n5=dshaloyield(epos,52,istat)
      IBflag(4)=0
c      IBflag(5)=0
c      IBflag(6)=0
      n4=dshaloyield(epos,52,istat)
      write(*,*) 'Photons: ',n4,n5
      write(*,*)   

      epos=m0/10000d0

c     loop to write out spectrum for a given model
c     format: E [GeV], flux without IB, total flux
      i=0
 1100 continue
 
      epos=epos+m0/50d0
 
      do 1509 i = 1, 12
        IBflag(i)=0
 1509 continue
      flux1=dsepspecm(epos,4)
 
      IBflag(1)=1
      IBflag(4)=1
      IBflag(5)=1
      IBflag(6)=1
      flux2=dsepspecm(epos,4)

      write(*,*) epos,flux1,flux2
      write(30,*) epos,flux1,flux2

c     end of loop to write spectrum for a given model      
 
      if (epos.lt.m0) goto 1100

c     end of loop to read in new models
      if (modelflag.eq.2) then
c        skip=skip+1
         goto 2000
      endif 
 1900 goto 1000

 2000 continue
  
      call SYSTEM_CLOCK(COUNT2, COUNT_RATE, COUNT_MAX)

      close (10)
        close (30)

      write (*,*) 'All spectra written out succesfully'
      write (*,*) 'Total time needed (in seconds): ',
     -     (COUNT2-COUNT)/(1d0*COUNT_RATE)
      stop
      end  






      subroutine read_model(lunit,nmodel,iend,ierr)
c
c     To read in a model from a file
c
c     NOTE: We here show how you can read in a model file and define all
c     DarkSUSY model parameters. For MSSM models, there exists a routine
c     src/su/dsgive_model.f that makes the definitions for you. 
c     We here read the model parameters and call that routine.

      implicit none
      include 'dsidtag.h'
      integer nmodel,lunit,iend,ierr
      real*8 at,ab,mqtild
      integer i
      character*40 message
 2000 format (1x,a12,7(1x,e14.8))
      ierr=0
c... When nmodel<0, skip -nmodel lines
      if (nmodel.lt.0) then
         do i=1,-nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
         return
      endif
c... If nmodel>0, read n-th model (assumes header line)
      if (nmodel.gt.0) then
         do i=1,nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
      endif
c... If nmodel=0, read next model
      read (lunit,2000,end=1000,err=3000) 
     &     idtag,mu,m2,ma,tanbe,mqtild,at,ab
c... modify/set additional parameters
      higloop=6  ! 5 = Full FeynHiggs;  6 = FeynHiggsFast
      ! W mass for unitarity of tree-level annihilation amplitudes
      call dsgive_model(mu,m2,ma,tanbe,mqtild,at,ab)
      return

 1000 continue
      iend=1
      ierr=0
      write (message,*) 'End of model file (unit=',lunit,')'
      call dswrite(1,0,message)
      return
 3000 continue
      ierr=1
      return
      end
