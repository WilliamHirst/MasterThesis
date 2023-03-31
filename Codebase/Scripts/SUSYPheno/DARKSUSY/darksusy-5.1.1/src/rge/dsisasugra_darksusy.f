      subroutine dsisasugra_darksusy(valid)
c=======================================================================
c  interface between ISASUGRA and DarkSUSY common blocks
c  author: E.A.Baltz, 2001 eabaltz@alum.mit.edu
c  modified by J. Edsjo, 2002-03-19 to set alph3
c  updated to Isajet 7.74 by J. Edsjo and E.A. Baltz, 2006-02-20
c  modified by P. Ullio 02-07-10, 02-11-21
c=======================================================================
      implicit none

      real*8 zgm,dsabsq,aux,masstau,massb,masst,mwtmp
      integer i,j,g,valid,itmp
      complex*16 chaumxtmp(2,2),chavmxtmp(2,2)

      real*8 thetax,thetay,cotgl,cotgr,singl,singr,cosgl,cosgr
      real*8 cosbeta

      include 'dsmssm.h'
      include 'dsisasugra.h'

c...First check that Higgs masses are OK in case ISASUGRA has not
c...spotted this problem
      if (.not.mss(29).gt.0.0d0) valid=9  ! H2
      if (.not.mss(30).gt.0.0d0) valid=9  ! H1
      if (.not.mss(31).gt.0.0d0) valid=9  ! H3
      if (.not.mss(32).gt.0.0d0) valid=9  ! Hc
      if (valid.gt.0) return

c     now set all physical masses and mixings
c     ISASUGRA provides the neutralino and chargino masses and mixings
c     in the convention of signed masses and real matrices
c     We also get the stop, sbottom and stau mixing angles
c     Higgs sector from ISASUGRA also, for consistency

c...First set tan(beta)
      tanbe=tgbetavar

c...Then set mu. Added by edsjo 2002-06-26
      mu=musugra

c...And set m1, m2 and m3 to allow us to calcualte the chargino
c...and neutralino matrices from scratch. Added by J. Edsjo 2002-11-26
      m1=gss(7)
      m2=gss(8)
      m3=gss(9)

c...Now set up constants; we need this call already here
c...since we call dschasct below
      call dssuconst  

c======= Sfermion Masses "1" = "L", "2" = "R" ===================
      mass(ksu1)=mss(2)
      mass(ksu2)=mss(3)
      mass(ksd1)=mss(4)
      mass(ksd2)=mss(5)
      mass(kss1)=mss(6)
      mass(kss2)=mss(7)
      mass(ksc1)=mss(8)
      mass(ksc2)=mss(9)
      mass(ksb1)=mss(10)
      mass(ksb2)=mss(11)
      mass(kst1)=mss(12)
      mass(kst2)=mss(13)

      mass(kse1)=mss(17)
      mass(kse2)=mss(18)
      mass(ksmu1)=mss(19)
      mass(ksmu2)=mss(20)
      mass(kstau1)=mss(21)
      mass(kstau2)=mss(22)

      mass(ksnue)=mss(14)
      mass(ksnumu)=mss(15)
      mass(ksnutau)=mss(16)

c========== Sfermion Mixings =================      
      do g=1,3
         do i=1,3
            slulmx(i,g) = cmplx(0.d0,0.d0)
         enddo
         do i=1,6
            sldlmx(i,g) = cmplx(0.d0,0.d0)
            sldrmx(i,g) = cmplx(0.d0,0.d0)
            squlmx(i,g) = cmplx(0.d0,0.d0)
            squrmx(i,g) = cmplx(0.d0,0.d0)
            sqdlmx(i,g) = cmplx(0.d0,0.d0)
            sqdrmx(i,g) = cmplx(0.d0,0.d0)
         enddo
         slulmx(g,g) = cmplx(1.d0,0.d0)
         sldlmx(g,g) = cmplx(1.d0,0.d0)
         sldrmx(g+3,g) = cmplx(1.d0,0.d0)
         squlmx(g,g) = cmplx(1.d0,0.d0)
         squrmx(g+3,g) = cmplx(1.d0,0.d0)
         sqdlmx(g,g) = cmplx(1.d0,0.d0)
         sqdrmx(g+3,g) = cmplx(1.d0,0.d0)
      enddo

c======= the mixing angles ======
c...The sign of these changed by J. Edsjo 2002-11-27
      mix_stau=thetal
      mix_sbot=thetab
      mix_stop=thetat

c...The matrices we now have would be the same as the ones from
c...DarkSUSY, except that q~_1 and q~_2 are sometimes flipped, i.e.
c...the mixing matrices
c...are multiplied by +- ( 0  1)
c...                     (-1  0)  compared to if DarkSUSY would calculte these

c======= stau mixing ===============
      sldlmx(3,3) = cmplx(cos(mix_stau),0.d0)
      sldlmx(6,3) = cmplx(-sin(mix_stau),0.d0)
      sldrmx(3,3) = cmplx(sin(mix_stau),0.d0)
      sldrmx(6,3) = cmplx(cos(mix_stau),0.d0)

c========= stop mixing ===========
      squlmx(3,3) = cmplx(cos(mix_stop),0.d0)
      squlmx(6,3) = cmplx(-sin(mix_stop),0.d0)
      squrmx(3,3) = cmplx(sin(mix_stop),0.d0)
      squrmx(6,3) = cmplx(cos(mix_stop),0.d0)

c======= sbottom mixing=========
      sqdlmx(3,3) = cmplx(cos(mix_sbot),0.d0)
      sqdlmx(6,3) = cmplx(-sin(mix_sbot),0.d0)
      sqdrmx(3,3) = cmplx(sin(mix_sbot),0.d0)
      sqdrmx(6,3) = cmplx(cos(mix_sbot),0.d0)

c      write(*,*) 'ISASUGRA matrices:'
c      write(*,*) 'squlmx(3,3)=',squlmx(3,3)
c      write(*,*) 'squlmx(6,3)=',squlmx(6,3)
c      write(*,*) 'squrmx(3,3)=',squrmx(3,3)
c      write(*,*) 'squrmx(6,3)=',squrmx(6,3)

c======== Soft Terms ==========
      do i=1,2
         asofte(i)=0.0d0
         asoftu(i)=0.0d0
         asoftd(i)=0.0d0
      enddo
      asofte(3)=aal
      asoftu(3)=aat
      asoftd(3)=aab
      
c====== Gluino Mass ==========
      mass(kgluin)=mss(1)

c======= Higgs Masses=======
      mass(kh2)=mss(29)
      mass(kh1)=mss(30)
      mass(kh3)=mss(31)
      mass(khc)=mss(32)
c====== sign reversed =========
      alpha=-alfah

c====== Neutralino and Chargino Masses
c====== ISASUGRA orders the neutralino and chargino mass eigenvalues
      mass(kn1)=abs(amz1ss)
      mass(kn2)=abs(amz2ss)
      mass(kn3)=abs(amz3ss)
      mass(kn4)=abs(amz4ss)

c========= Neutralino Mixing ===========

      do 100 i=1,4
         do 200 j=1,4
            if (amziss(i).le.0.0d0) then
               neunmx(i,j)=cmplx(zmixss(5-j,i),0.0d0)
            else
               neunmx(i,j)=cmplx(0.0d0,zmixss(5-j,i))
            endif
 200     continue
         do 300 j=1,2
            neunmx(i,j)=-neunmx(i,j)
 300     continue
 100  continue

      zg = dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
      zgm = dsabsq(neunmx(1,3))+dsabsq(neunmx(1,4))
      lgzg = log10(zg)-log10(zgm)
      kln = 1
      lsp = kn(kln)

c======= Chargino Mixing ==========

c...This approach extracts the chargino matrices from ISASUGRA, but
c...with wrong signs of the eigenvalues sometimes
c...We will use the ISASUGRA values, but with the signs determined
c...from DarkSUSY

      cotgl=1.0d0/tan(gammal)
      cotgr=1.0d0/tan(gammar)
 
      singl=sin(gammal)
      singr=sin(gammar)
      cosgl=cos(gammal)
      cosgr=cos(gammar)

      chaumxtmp(1,2)=singl
      chavmxtmp(1,2)=singr
      chaumxtmp(1,1)=-cosgl
      chavmxtmp(1,1)=-cosgr
      chaumxtmp(2,2)=cosgl
      chavmxtmp(2,2)=cosgr
      chaumxtmp(2,1)=singl
      chavmxtmp(2,1)=singr

c...Use the m_W value of ISASUGRA and recalculate the charigno
c...matrices

      mwtmp=mass(kw) 
      mass(kw)=amw   ! Use mw=79.9522095 from isasugra to be consistent
      call dschasct
      mass(kw)=mwtmp

c...Now set use the isasugra matrices, but with the DarkSUSY signs
      do i=1,2
        do j=1,2
          chaumx(i,j)=abs(chaumxtmp(i,j))
     &      *sign(1.0d0,dreal(chaumx(i,j)))
          chavmx(i,j)=abs(chavmxtmp(i,j))
     &      *sign(1.0d0,dreal(chavmx(i,j)))
c          chaumx(i,j)=chaumxtmp(i,j)
c          chavmx(i,j)=chavmxtmp(i,j)
        enddo
      enddo

      mass(kcha2)=abs(amw1ss)
      mass(kcha1)=abs(amw2ss)
      
c----------------------------------------------------------------------
c...The things set below are not needed if only ISASUGRA is used for
c...the particle spectrum. However, if we want to calculate e.g. the Higgs
c...widths with HDECAY we need to have the low-energy mass parameters as
c...input and not only the mass matrices. This is in a sense doing
c...things double since HDECAY will recalculate things already calculated
c...elsewhere, but the interface to HDECAY is kept simple.
c...This addition is made by Joakim Edsjo, edsjo@physto.se, 02-09-11
      mass2l(1)=gss(16)
      mass2e(1)=gss(15)
      mass2q(1)=gss(19)
      mass2u(1)=gss(18)
      mass2d(1)=gss(17)

      mass2l(2)=mass2l(1)
      mass2e(2)=mass2e(1)
      mass2q(2)=mass2q(1)
      mass2u(2)=mass2u(1)
      mass2d(2)=mass2d(1)

      mass2l(3)=gss(21)
      mass2e(3)=gss(20)
      mass2q(3)=gss(24)
      mass2u(3)=gss(23)
      mass2d(3)=gss(22)

c      asofte(3)=gss(10)  ! already done above
c      asoftu(3)=gss(12)  ! already done above
c      asoftd(3)=gss(11)  ! already done above

c...Now check that this implementaion is OK. We have checked this
c...extensively, but to be on the safe side, let's check it for every
c...model. Note: With the loop corrections in isajet, this check
c...cannot longer be used due to all false alarms

c      call dsisasugra_check(valid)

c      write (*,*) 'AAAAA gorge = ',gorge

      return
      end
