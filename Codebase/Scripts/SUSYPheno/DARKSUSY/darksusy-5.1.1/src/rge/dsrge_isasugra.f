      subroutine dsrge_isasugra(unphys,valid)
c=======================================================================
c  interface to ISASUGRA (ISAJET 7.78) routines for SUSY spectra
c  author: E.A.Baltz, 2001 eabaltz@alum.mit.edu
c
c  if valid is non-zero, the model is no good
c  the valid flag is equal to the isasugra nogood flag:
c  valid  reason for model being bad
c  -----  --------------------------
c      1  TACHYONIC PARTICLES
c      2  NO EW SYMMETRY BREAKING
c      3  M(H_P)^2<0
c      4  YUKAWA>10
c      5  Z1SS NOT LSP
c      7  XT EWSB IS BAD
c      8  MHL^2<0
c      9  if in our check any Higgs mass is NaN
c  The following are not set, but can be set by uncommenting
c  the appropriate lines in dsisasugra_check.f
c  10-14  dsisasugra_check has reported a possible error in the interface
c         while checking that chargino, neutralino and sfermion mass
c         matrices are diagonlized by isasugra
c Updated to ISAJET 7.74 by J. Edsjo and E.A. Baltz, 2006-02-20
c Updated to ISAJET 7.78 by J. Edsjo, 2008-06-05.
c Updated to ISAJET 7.79 by J. Edsjo, 2010-03-08
c Updated to ISAJET 7.81 by P. Gondolo, 2011-09-04
c=======================================================================
      implicit none

      include 'dsmssm.h'
      include 'dsisasugra.h'
      integer unphys,valid,i,iallow,IMHL,IMHC
      real*8 mx
      real dssug(6)
      logical dsisnan

      unphys=0
      valid=0

c-----------------------------------------------------------------------
c     for simplicty, start with mSUGRA: m0,mhalf,tan(beta),A,sgn(mu)      
      dssug(1) = m0var
      dssug(2) = mhfvar
      dssug(3) = a0var
      dssug(4) = tgbetavar
      dssug(5) = sgnmuvar
      dssug(6) = mass(kt)
      
c-----------------------------------------------------------------------
c     here is the call to the ISASUGRA routine
c      write(*,*) dssug

c...Re-initialize isasugra
      GORGE = .false.  ! start-up value from aldata.f

c... Solve RGEs
      call sugra(dssug(1),dssug(2),dssug(3),
     +     dssug(4),dssug(5),dssug(6),1)      

c... Calculate masses from last iteration
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ dssug(6),iallow,1,IMHL,IMHC)

      valid=nogood
c      if (nogood.ne.0) then
c        write(*,*) 'dsrge_isasugra: nogood=',nogood
c      endif
      if (valid.gt.0) then
         write (*,*) 'Error in ISAJET. NOGOOD=',nogood
         return
      endif

c------------ transfer ISASUGRA to DarkSUSY common blocks 
      call dsisasugra_darksusy(valid)

      if (valid.gt.0) then
         write (*,*) 'Error in DarkSUSY-ISAJET interface. valid=',valid
         return
      endif
            
c---Check that the mass spectrum is OK
      mx=mass(kn(1))
      unphys=0
c...Check that neutralino is LSP
      do i=1,2
        if (mass(kcha(i)).lt.mx) unphys=ibset(unphys,0)
      enddo
      do i=1,3
        if (mass(ksnu(i)).lt.mx) unphys=ibset(unphys,1)
      enddo
      do i=1,6
        if (mass(ksl(i)).lt.mx) unphys=ibset(unphys,2)
        if (mass(ksqu(i)).lt.mx) unphys=ibset(unphys,3)
        if (mass(ksqd(i)).lt.mx) unphys=ibset(unphys,4)
      enddo
      if (mass(kgluin).lt.mx) unphys=ibset(unphys,5)

c...Check that masses are not NaN
      if (dsisnan(mass(kh1))) unphys=ibset(unphys,6)
      if (dsisnan(mass(kh2))) unphys=ibset(unphys,6)
      if (dsisnan(mass(kh3))) unphys=ibset(unphys,6)
      if (dsisnan(mass(khc))) unphys=ibset(unphys,6)
      do i=1,4
        if (dsisnan(mass(kn(i)))) unphys=ibset(unphys,6)
      enddo
      do i=1,2
        if (dsisnan(mass(kcha(i)))) unphys=ibset(unphys,6)
      enddo
      do i=1,3
        if (dsisnan(mass(ksnu(i)))) unphys=ibset(unphys,6)
      enddo
      do i=1,6
        if (dsisnan(mass(ksl(i)))) unphys=ibset(unphys,6)
        if (dsisnan(mass(ksqu(i)))) unphys=ibset(unphys,6)
        if (dsisnan(mass(ksqd(i)))) unphys=ibset(unphys,6)
      enddo
      if (dsisnan(mass(kgluin))) unphys=ibset(unphys,6)

c...Check if SSMSSM has warned about something else...
      if (unphys.eq.0.and.iallow.ne.0) then
        write(*,*) 'WARNING: iallow =',iallow,
     &    ' from SSMSSM (Isajet)'
        unphys=ibset(unphys,7)
      endif

      end
