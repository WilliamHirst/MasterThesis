      subroutine dsntsunread

***********************************************************************
*** Reads in data about the solar model used and stores it in a
*** common block (as described in dssun.h).
*** Author: Joakim Edsjo
*** Date: 2003-11-25
*** Modified: 2004-01-28 (calculates potential instead of reading file)
***********************************************************************

      implicit none
      include 'dssun.h'
      include 'dsmpconst.h'
      include 'dsdirver.h'

      real*8 dsntsunpotint,dsntsuncdensint
      integer i,fl,l,m

      logical sdread
      data sdread/.false./
      save sdread

      character*200 file,filene,scr

      real*8 tmp1,tmp2,tmp3,totfr,totfrheavy

c...If already initialized, return
      if (sdread) then
        return
      endif

      sdread=.true.

c...Determine which file we should read in
c...This is from the Standard solar model, BP2000
c      file=dsinstall//'share/DarkSUSY/'//'bp2000stdmodel.dat'
c...This is from the Standard solar model, BS05(OP)
      file=dsinstall//'share/DarkSUSY/'//'bs05op.dat'
c...This is form the alternative Standard solar model, BS05(AGS,OP)
c...with new heavy element measurements (fits worse with helioseismology
c...so we don't use it as a defualt).
c      file=dsinstall//'share/DarkSUSY/'//'bs05_agsop.dat'
c...The electron density is from the Standard solar model, BS05(OP)
      filene=dsinstall//'share/DarkSUSY/'//'nele_bs05op.dat'

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 40     if (file(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            file(m:m)=file(m+1:m+1)
          enddo
          if (fl.eq.l) goto 50
          goto 40
        endif
      enddo
 50   continue

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 60     if (filene(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            filene(m:m)=filene(m+1:m+1)
          enddo
          if (fl.eq.l) goto 70
          goto 60
        endif
      enddo
 70   continue


c...First define the abundances of heavy elements not listed in the
c...Bahcall et al standard solar model file
c...The number fractions are given by n_i = n_H log10(A-12)
c...The abundances listed here are from N. Grevesse and A.J. Sauval,
c...Space Science Reviews 85 (1998) 161.
      sdabund(7)  = 8.08d0  ! Ne
      sdabund(8)  = 6.33d0  ! Na
      sdabund(9)  = 7.58d0  ! Mg
      sdabund(10) = 6.47d0  ! Al
      sdabund(11) = 7.55d0  ! Si
      sdabund(12) = 7.33d0  ! S
      sdabund(13) = 6.40d0  ! Ar
      sdabund(14) = 6.36d0  ! Ca
      sdabund(15) = 7.50d0  ! Fe
      sdabund(16) = 6.25d0  ! Ni
      

c...Mass numbers
      sdaa(1)  = 1.0d0  ! H
      sdaa(2)  = 4.0d0  ! He4
      sdaa(3)  = 3.0d0  ! He3
      sdaa(4)  = 12.0d0 ! C12
      sdaa(5)  = 14.0d0 ! N14
      sdaa(6)  = 16.0d0 ! O8
      sdaa(7)  = 20.0d0 ! Ne
      sdaa(8)  = 23.0d0 ! Na
      sdaa(9)  = 24.0d0 ! Mg
      sdaa(10) = 27.0d0 ! Al
      sdaa(11) = 28.0d0 ! Si
      sdaa(12) = 32.0d0 ! S
      sdaa(13) = 40.0d0 ! Ar
      sdaa(14) = 40.0d0 ! Ca
      sdaa(15) = 56.0d0 ! Fe
      sdaa(16) = 59.0d0 ! Ni

c...Masses
      do i=1,16
        sdma(i)=sdaa(i)*(m_p+m_n)/2.0d0
      enddo

c...Now read in data from solar model data file
      write(*,*) 'dsntsunread: Opening file ',file
      open(unit=13,file=file,
     &  form='formatted',status='old')
      sdn=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,23
        read(13,'(A)',end=110) scr
      enddo

c...Calculate the total mass fraction (relative to H) of elements
c...heavier than O. This is just for normalization of the fractions
c...as a function of radius below
      totfrheavy=0.0d0
      do i=7,16
        totfrheavy=totfrheavy+sdma(i)*10**(sdabund(i)-12.0d0)
      enddo

c...Read in table and calculate remaining mass fractions from
c...average abundances sdabund
 100  read(13,*,end=110) sdm(sdn),sdr(sdn),tmp1,sdrho(sdn),tmp2,tmp3,
     &  sdmfr(1,sdn),sdmfr(2,sdn),sdmfr(3,sdn),sdmfr(4,sdn),
     &  sdmfr(5,sdn),sdmfr(6,sdn)
      totfr=sdmfr(1,sdn)+sdmfr(2,sdn)+sdmfr(3,sdn)+sdmfr(4,sdn)
     &  +sdmfr(5,sdn)+sdmfr(6,sdn)
c...Add the heavy elements (>O16)
      do i=7,16
        sdmfr(i,sdn)=(1.0d0-totfr)*sdma(i)*10**(sdabund(i)-12.0d0)
     &    /totfrheavy      
      enddo
      sdn=sdn+1
      if (sdn.gt.sdmax-1) then  ! need one more entry for last line
        goto 110
      endif
      goto 100
 110  continue
      sdn=sdn-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdm(1)=0.0d0
      sdr(1)=0.0d0
      sdrho(1)=sdrho(2)
      do i=1,16
        sdmfr(i,1)=sdmfr(i,2)
      enddo

c...Now add the last line with r/r_sun=1
      sdn=sdn+1
      sdm(sdn)=1.0d0
      sdr(sdn)=1.0d0
      sdrho(sdn)=0.0d0
      do i=1,16
        sdmfr(i,sdn)=sdmfr(i,sdn-1)
      enddo


c...Electron density
c...Now read in data from solar model electron density file
      write(*,*) 'dsntsunread: Opening file ',filene
      open(unit=13,file=filene,
     &  form='formatted',status='old')
      sdnne=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,6
        read(13,'(A)',end=210) scr
      enddo

c...Read in table and calculate remaining mass fractions from
c...average abundances sdabund
 200  read(13,*,end=210) sdrne(sdnne),sdne(sdnne)
      sdnne=sdnne+1
      if (sdnne.gt.sdmax) then  ! need one more entry for last line
        write(*,*) 'ERROR in dsntsunread: array too small.'
        write(*,*) 'Increase sdmax.'
        goto 210
      endif
      goto 200
 210  continue
      sdnne=sdnne-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdne(1)=sdne(2)
      sdrne(1)=0.0d0

c...Sun parameters
      m_sun = 1.9889d30   ! solar mass in kg
      r_sun = 6.9598d8    ! solar radius in m


c=======================================================================
c=======================================================================
c=======================================================================

c...Now calculate and tabulate the potential inside the Sun

      write(*,*) 'dsntsunread: tabulating potential inside the Sun...'
      do i=1,sdn
        sdphi(i)=dsntsunpotint(sdr(i)*r_sun)
      enddo
      write(*,*) '  ...done'

      write(*,*) 
     &  'dsntsunread: tabulating column density inside the Sun...'
      do i=1,sdn
        sdcdens(i,0)=dsntsuncdensint(sdr(i)*r_sun,'N')
        sdcdens(i,1)=dsntsuncdensint(sdr(i)*r_sun,'p')
        sdcdens(i,2)=dsntsuncdensint(sdr(i)*r_sun,'n')
      enddo
      write(*,*) '  ...done'
      cd_sun(0)=sdcdens(sdn,0)  ! total (p+n) column density g/cm^2
      cd_sun(1)=sdcdens(sdn,1)  ! total proton column density g/cm^2
      cd_sun(2)=sdcdens(sdn,2)  ! total neutron column density g/cm^2

c      write(*,*) 'Ip = ',cd_sun(1)/r_sun/100.0d0
c      write(*,*) 'In = ',cd_sun(2)/r_sun/100.0d0

      return

      end


      

