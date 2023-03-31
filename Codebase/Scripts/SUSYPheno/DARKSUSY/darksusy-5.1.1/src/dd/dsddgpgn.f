      subroutine dsddgpgn(gps,gns,gpa,gna)
c_______________________________________________________________________
c  neutralino nucleon four-fermion couplings
c  common:
c    'dssusy.h' - file with susy common blocks
c  output:
c    gps, gns : proton and neutron scalar four-fermion couplings
c    gpa, gna : proton and neutron axial four-fermion couplings
c    units: GeV^-4
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c    2008-02-16 paolo gondolo double spin-dependent squark contribution
c    2008-02-16 paolo gondolo undo running quark masses in couplings
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsprep.h'
      include 'dsddcom.h'
      include 'dshacom.h'
      include 'dsmpconst.h'

      real*8 gps,gns,gpa,gna
      integer q,g,i,z,kx
      real*8 mq,aux,rmq,
     &     du,dd,ds,dsabsq,tmp1,tmp2,
     &     ftu(2),ftd(2),fts(2),ftg(2),asmt,asmb,
     &     dsdddn1,dsdddn2,dsdddn3,dsdddn4,dsddo
      real*8 f1q,df1q,fp,fn,aqh,aqd,aqs,fgh,fgd,fgs,agh,agd,ags
      ! fp,fn,f1q,df1q,du,dd,ds are twice those in drees-nojiri
      real*8 rmass(12) ! temporary array for lepton and quark running masses
                       ! currently extracted from yukawa 
                       ! one day we will have rmass(k,q)

      kx = kn(kln)
      mx = mass(kx)

c...... temporary fix for running masses (PG 2008-02-16)
      aux=dsqrt(2.d0)*mass(kw)/g2weak
      do i=1,3
         rmass(knu(i)) = aux*yukawa(knu(i))*sinbe
         rmass(kl(i)) = aux*yukawa(kl(i))*cosbe
         rmass(kqu(i)) = aux*yukawa(kqu(i))*sinbe
         rmass(kqd(i)) = aux*yukawa(kqd(i))*cosbe
      enddo
c...... end of temporary fix for running masses


c
c     SPIN INDEPENDENT
c

      fp = 0.d0
      fn = 0.d0
      if (dddn.eq.1) then
c         write (*,*) 'The Drees-Nojiri treatment of ',
c     &        'neutralino-nucleon scattering is not available'
c         write (*,*)  'DarkSUSY will stop'
c         stop
c
c     the drees-nojiri treatment
c

         ftu(1) = ftp(7)
         ftd(1) = ftp(8)
         fts(1) = ftp(10)
         ftg(1) = 2.d0/27.d0*(1.d0-ftu(1)-ftd(1)-fts(1))
         ftu(2) = ftn(7)
         ftd(2) = ftn(8)
         fts(2) = ftn(10)
         ftg(2) = 2.d0/27.d0*(1.d0-ftu(2)-ftd(2)-fts(2))
         asmt = alph3/(1.d0+(11.d0-10.d0/3.d0)/(4.d0*pi)*alph3
     &              * dlog(mass(kt)**2/mass(kz)**2))
         asmb = alph3/(1.d0+(11.d0-10.d0/3.d0)/(4.d0*pi)*alph3
     &              * dlog(mass(kb)**2/mass(kz)**2))
c     g
         fgh = 0.d0
         fgd = 0.d0
         fgs = 0.d0
         fgh = fgh              ! higgs: c,b,t triangle
     &        +dreal(gr(kh1,kx,kx)*gl(kh1,kc,kc))/
     &        mass(kh1)**2/rmass(kc)
     &        +dreal(gr(kh2,kx,kx)*gl(kh2,kc,kc))/
     &        mass(kh2)**2/rmass(kc)
     &        +dreal(gr(kh1,kx,kx)*gl(kh1,kb,kb))/
     &        mass(kh1)**2/rmass(kb)
     &        +dreal(gr(kh2,kx,kx)*gl(kh2,kb,kb))/
     &        mass(kh2)**2/rmass(kb)
     &        +dreal(gr(kh1,kx,kx)*gl(kh1,kt,kt))/
     &        mass(kh1)**2/rmass(kt)
     &        +dreal(gr(kh2,kx,kx)*gl(kh2,kt,kt))/
     &        mass(kh2)**2/rmass(kt)
         do i=1,6
            fgh = fgh+1.d0/8.d0*( ! higgs: squarks triangle
     &           dreal(gr(kh1,kx,kx)*
     &           gl(kh1,ksqu(i),ksqu(i)))/
     &           mass(kh1)**2/mass(ksqu(i))**2
     &           +dreal(gr(kh2,kx,kx)*
     &           gl(kh2,ksqu(i),ksqu(i)))/
     &           mass(kh2)**2/mass(ksqu(i))**2+
     &           dreal(gr(kh1,kx,kx)*
     &           gl(kh1,ksqd(i),ksqd(i)))/
     &           mass(kh1)**2/mass(ksqd(i))**2
     &           +dreal(gr(kh2,kx,kx)*
     &           gl(kh2,ksqd(i),ksqd(i)))/
     &           mass(kh2)**2/mass(ksqd(i))**2)
            fgd = fgd           ! squarks: c,b,t trace in a^2-b^2
     &           - 0.5d0/rmass(kc)*dreal(gl(ksqu(i),kx,kc)*
     &           conjg(gr(ksqu(i),kx,kc)))*
     &           dsdddn1(mass(ksqu(i)),mass(kc),mx)
     &           - 0.5d0/rmass(kb)*dreal(gl(ksqd(i),kx,kb)*
     &           conjg(gr(ksqd(i),kx,kb)))*
     &           dsdddn1(mass(ksqd(i)),mass(kb),mx)
     &           - 0.5d0/rmass(kt)*dreal(gl(ksqu(i),kx,kt)*
     &           conjg(gr(ksqu(i),kx,kt)))*
     &           dsdddn1(mass(ksqu(i)),mass(kt),mx)
            fgs = fgs           ! squarks: u,d,s,c,b,t trace in a^2+b^2
     &           + 0.125*mx*(dsabsq(gl(ksqu(i),kx,ku))+
     &           dsabsq(gr(ksqu(i),kx,ku)))*
     &           dsdddn2(mass(ksqu(i)),mass(ku),mx)
     &           + 0.125*mx*(dsabsq(gl(ksqd(i),kx,kd))+
     &           dsabsq(gr(ksqd(i),kx,kd)))*
     &           dsdddn2(mass(ksqd(i)),mass(kd),mx)
     &           + 0.125*mx*(dsabsq(gl(ksqd(i),kx,ks))+
     &           dsabsq(gr(ksqd(i),kx,ks)))*
     &           dsdddn2(mass(ksqd(i)),mass(ks),mx)
     &           + 0.125*mx*(dsabsq(gl(ksqu(i),kx,kc))+
     &           dsabsq(gr(ksqu(i),kx,kc)))*
     &           dsdddn2(mass(ksqu(i)),mass(kc),mx)
     &           + 0.125*mx*(dsabsq(gl(ksqd(i),kx,kb))+
     &           dsabsq(gr(ksqd(i),kx,kb)))*
     &           dsdddn2(mass(ksqd(i)),mass(kb),mx)
     &           + 0.125*mx*(dsabsq(gl(ksqu(i),kx,kt))+
     &           dsabsq(gr(ksqu(i),kx,kt)))*
     &           dsdddn2(mass(ksqu(i)),mass(kt),mx)
         enddo
         do z=1,2
            agh = fgh*ftg(z)
            agd = fgd*ftg(z)
            ags = fgs*ftg(z)
c     u,d,s
            aqh =               ! higgs: u,d,s tree level
     &           +dreal(gr(kh1,kx,kx)*gl(kh1,ku,ku))/
     &           mass(kh1)**2/rmass(ku)*ftu(z)
     &           +dreal(gr(kh2,kx,kx)*gl(kh2,ku,ku))/
     &           mass(kh2)**2/rmass(ku)*ftu(z)
     &           +dreal(gr(kh1,kx,kx)*gl(kh1,kd,kd))/
     &           mass(kh1)**2/rmass(kd)*ftd(z)
     &           +dreal(gr(kh2,kx,kx)*gl(kh2,kd,kd))/
     &           mass(kh2)**2/rmass(kd)*ftd(z)
     &           +dreal(gr(kh1,kx,kx)*gl(kh1,ks,ks))/
     &           mass(kh1)**2/rmass(ks)*fts(z)
     &           +dreal(gr(kh2,kx,kx)*gl(kh2,ks,ks))/
     &           mass(kh2)**2/rmass(ks)*fts(z)
            aqd = 0.d0
            aqs = 0.d0
            do i=1,6
               aqd = aqd        ! squarks: u,d,s tree a^2-b^2
     &              - 0.5d0*dreal(gl(ksqu(i),kx,ku)*
     &              conjg(gr(ksqu(i),kx,ku)))/rmass(ku)/
     &              (mass(ksqu(i))**2-(mx+mass(ku))**2)*ftu(z)
     &              - 0.5d0*dreal(gl(ksqd(i),kx,kd)*
     &              conjg(gr(ksqd(i),kx,kd)))/rmass(kd)/
     &              (mass(ksqd(i))**2-(mx+mass(kd))**2)*ftd(z)
     &              - 0.5d0*dreal(gl(ksqd(i),kx,ks)*
     &              conjg(gr(ksqd(i),kx,ks)))/rmass(ks)/
     &              (mass(ksqd(i))**2-(mx+mass(ks))**2)*fts(z)
               aqs = aqs        ! squarks: u,d,s tree & twist-2 a^2+b^2
     &              + 0.125*mx*(dsabsq(gl(ksqu(i),kx,ku))+
     &              dsabsq(gr(ksqu(i),kx,ku)))/
     &              (mass(ksqu(i))**2-(mx+mass(ku))**2)**2*(ftu(z)+
     &              3.d0*dsddo(1,z,mass(ksqu(i))**2-(mx+mass(ku))**2))
     &              + 0.125*mx*(dsabsq(gl(ksqd(i),kx,kd))+
     &              dsabsq(gr(ksqd(i),kx,kd)))/
     &              (mass(ksqd(i))**2-(mx+mass(kd))**2)**2*(ftd(z)+ 
     &              3.d0*dsddo(2,z,mass(ksqd(i))**2-(mx+mass(kd))**2))
     &              + 0.125*mx*(dsabsq(gl(ksqd(i),kx,ks))+
     &              dsabsq(gr(ksqd(i),kx,ks)))/
     &              (mass(ksqd(i))**2-(mx+mass(ks))**2)**2*(fts(z)+
     &              3.d0*dsddo(4,z,mass(ksqu(i))**2-(mx+mass(ks))**2))
               aqs = aqs        ! squarks: c twist-2 a^2+b^2
     &              + 0.125*mx*(dsabsq(gl(ksqu(i),kx,kc))+
     &              dsabsq(gr(ksqu(i),kx,kc)))/
     &              (mass(ksqu(i))**2-(mx+mass(kc))**2)**2*
     &              3.d0*dsddo(3,z,mass(ksqu(i))**2-(mx+mass(kc))**2)
               tmp1 =           ! squarks: b twist-2 a^2+b^2
     &              + 0.125*mx*(dsabsq(gl(ksqd(i),kx,kb))+
     &              dsabsq(gr(ksqd(i),kx,kb)))/
     &              (mass(ksqd(i))**2-(mx+mass(kb))**2)**2*
     &              3.d0*dsddo(6,z,mass(ksqd(i))**2-(mx+mass(kb))**2)
               tmp2 = asmb*dsddo(0,z,mass(kb)**2)/4.d0/pi*(
     &              +0.125*mx*(dsabsq(gl(ksqd(i),kx,kb))+
     &              dsabsq(gr(ksqd(i),kx,kb)))*
     &              dsdddn4(mass(ksqd(i)),mass(kb),mx)
     &              - 0.5d0*dreal(gl(ksqd(i),kx,kb)*
     &              conjg(gr(ksqd(i),kx,kb)))/rmass(kb)*
     &              dsdddn3(mass(ksqd(i)),mass(kb),mx))
               if (dabs(tmp1).lt.dabs(tmp2)) then
                  aqs = aqs + tmp1
               else
                  aqs = aqs + tmp2
               endif
               aqs = aqs        ! squarks: t twist-2 a^2+b^2
     &              +asmt*dsddo(0,z,mass(kt)**2)/4.d0/pi*(
     &              +0.125*mx*(dsabsq(gl(ksqu(i),kx,kt))+
     &              dsabsq(gr(ksqu(i),kx,kt)))*
     &              dsdddn4(mass(ksqu(i)),mass(kt),mx)
     &              - 0.5d0*dreal(gl(ksqu(i),kx,kt)*
     &              conjg(gr(ksqu(i),kx,kt)))/rmass(kt)*
     &              dsdddn3(mass(ksqu(i)),mass(kt),mx))
            enddo
            if (z.eq.1) then
               fp = fp + m_p * (aqh+aqd+aqs+agh+agd+ags)
               if (prtlevel.ge.1) then
                  write (*,*) '(p) aqh=',aqh,' aqd=',aqd, ' aqs=',aqs
                  write (*,*) '(p) agh=',agh,' agd=',agd, ' ags=',ags
               endif
            endif
            if (z.eq.2) then
               fn = fn + m_n * (aqh+aqd+aqs+agh+agd+ags)
               if (prtlevel.ge.1) then
                  write (*,*) '(n) aqh=',aqh,' aqd=',aqd, ' aqs=',aqs
                  write (*,*) '(n) agh=',agh,' agd=',agd, ' ags=',ags
               endif
            endif
         enddo
      else
c     
c     without poles, without twist-2 terms
c
         do g=1,3
c     up-quarks
            mq = mass(kqu(g))
            rmq = rmass(kqu(g))
            f1q =
     &           gr(kh1,kx,kx)*gl(kh1,kqu(g),kqu(g))/
     &           mass(kh1)**2
     &           +gr(kh2,kx,kx)*gl(kh2,kqu(g),kqu(g))/
     &           mass(kh2)**2
            do i=1,6
               if (ddpole.eq.1) then
                  df1q =
     &                 - 0.5*gl(ksqu(i),kx,kqu(g))*
     &                 conjg(gr(ksqu(i),kx,kqu(g)))/
     &                 (mass(ksqu(i))**2-(mx+mq)**2)
               else
                  df1q =
     &                 - 0.5*gl(ksqu(i),kx,kqu(g))*
     &                 conjg(gr(ksqu(i),kx,kqu(g)))/
     &                 (mass(ksqu(i))**2)
               endif
c     &        + 0.125*(glsnf(i,kln,q)**2+conjg(grsnf(i,kln,q))**2)/
c     &        (mass(ksqu(i))**2-(mx+mq)**2)**2
               f1q = f1q + df1q
            enddo
            q = 5+2*g
            fp = fp + f1q*ftp(q)*(m_p/rmq)
            fn = fn + f1q*ftn(q)*(m_n/rmq)
c     down-quarks
            mq = mass(kqd(g))
            rmq = rmass(kqd(g))
            f1q =
     &           gr(kh1,kx,kx)*gl(kh1,kqd(g),kqd(g))/
     &           mass(kh1)**2
     &           +gr(kh2,kx,kx)*gl(kh2,kqd(g),kqd(g))/
     &           mass(kh2)**2
            do i=1,6
               if (ddpole.eq.1) then
                  df1q =
     &                 - 0.5*gl(ksqd(i),kx,kqd(g))*
     &                 conjg(gr(ksqd(i),kx,kqd(g)))/
     &                 (mass(ksqd(i))**2-(mx+mq)**2)
               else
                  df1q =
     &                 - 0.5*gl(ksqd(i),kx,kqd(g))*
     &                 conjg(gr(ksqd(i),kx,kqd(g)))/
     &                 (mass(ksqd(i))**2)
               endif
c     &        + 0.125*(glsnf(i,kln,q)**2+conjg(grsnf(i,kln,q))**2)/
c     &        (mass(ksqd(i))**2-(mx+mq)**2)**2
               f1q = f1q + df1q
            enddo
            q = 6+2*g
            fp = fp + f1q*ftp(q)*(m_p/rmq)
            fn = fn + f1q*ftn(q)*(m_n/rmq)
         enddo
      endif
      gps = fp   ! in gev^-4
      gns = fn   ! in gev^-4

c
c     SPIN DEPENDENT
c

      du = 0.5*gr(kz,kx,kx)*(gl(kz,ku,ku)-gr(kz,ku,ku))/
     &     mass(kz)**2
      dd = 0.5*gr(kz,kx,kx)*(gl(kz,kd,kd)-gr(kz,kd,kd))/
     &     mass(kz)**2
      ds = 0.5*gr(kz,kx,kx)*(gl(kz,ks,ks)-gr(kz,ks,ks))/
     &     mass(kz)**2
      do i=1,6
        if (ddpole.eq.1.or.dddn.eq.1) then
          du = du
     &       + 0.25*(dsabsq(gl(ksqu(i),kx,ku))+
     &                dsabsq(gr(ksqu(i),kx,ku)))/
     &         (mass(ksqu(i))**2-(mx+mass(ku))**2)
          dd = dd
     &       + 0.25*(dsabsq(gl(ksqd(i),kx,kd))+
     &                dsabsq(gr(ksqd(i),kx,kd)))/
     &         (mass(ksqd(i))**2-(mx+mass(kd))**2)
          ds = ds
     &       + 0.25*(dsabsq(gl(ksqd(i),kx,ks))+
     &                dsabsq(gr(ksqd(i),kx,ks)))/
     &         (mass(ksqu(i))**2-(mx+mass(ks))**2)
        else
          du = du
     &       + 0.25*(dsabsq(gl(ksqu(i),kx,ku))+
     &                dsabsq(gr(ksqu(i),kx,ku)))/
     &         (mass(ksqu(i))**2)
          dd = dd
     &       + 0.25*(dsabsq(gl(ksqd(i),kx,kd))+
     &                dsabsq(gr(ksqd(i),kx,kd)))/
     &         (mass(ksqd(i))**2)
          ds = ds
     &       + 0.25*(dsabsq(gl(ksqd(i),kx,ks))+
     &                dsabsq(gr(ksqd(i),kx,ks)))/
     &         (mass(ksqu(i))**2)
         endif
      enddo
      gpa = du*delu+dd*deld+ds*dels ! in gev^-4
      gna = du*deld+dd*delu+ds*dels ! in gev^-4

      return
      end
