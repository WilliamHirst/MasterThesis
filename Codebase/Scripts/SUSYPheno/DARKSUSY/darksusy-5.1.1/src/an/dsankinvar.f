      subroutine dsankinvar(p,costheta,kp1,kp2,kpk,kp3,kp4)

c=======================================================================
c  calculate kinematical variables for s-, t-, and u-diagram
c  author: joakim edsjo, edsjo@physto.se
c  date: 97-01-09
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsandiacom.h'
      integer kp1,kp2,kpk,kp3,kp4
      integer kp1old,kp2old,kp3old,kp4old
      real*8 p,costheta
      real*8 pold,cthold
      real*8 sintheta,pp2,kk2,mp12,mp22,mp32,mp42
      logical cthch,pch,inich,finch
      data kp1old,kp2old,kp3old,kp4old/4*-1/
      data pold,cthold / -1.0d0,-5.0d0/
      save pold,cthold,kp1old,kp2old,kp3old,kp4old

c----------------------------------------- check changes since last call

      if (kp1.eq.0) then ! initialize for new point
        kp1old=-1
        kp2old=-1
        kp3old=-1
        kp4old=-1
        pold=-1.0d0
        cthold=-5.0d0
        return
      endif
      if (cthold.ne.costheta) then
        cthold=costheta
        cthch=.true.
      else
        cthch=.false.
      endif
      if (pold.ne.p) then
        pold=p
        pch=.true.
      else
        pch=.false.
      endif
      if (kp1.ne.kp1old.or.kp2.ne.kp2old) then
        kp1old=kp1
        kp2old=kp2
        inich=.true.
      else
        inich=.false.
      endif
      if (kp3.ne.kp3old.or.kp4.ne.kp4old) then
        kp3old=kp3
        kp4old=kp4
        finch=.true.
      else
        finch=.false.
      endif

c---------------------------------------------------- wigner d-functions
      if (cthch) then
        sintheta = sqrt(1-costheta)*sqrt(1+costheta)
        wd(1,1,1) = 0.5*(1+costheta)
        wd(1,1,0) = -sintheta/sqrt(2.0)
        wd(1,1,-1) = 0.5*(1-costheta)
        wd(1,0,1) = -wd(1,1,0)
        wd(1,0,0) = costheta
        wd(1,0,-1) = wd(1,1,0)
        wd(1,-1,1) = wd(1,1,-1)
        wd(1,-1,0) = wd(1,0,1)
        wd(1,-1,-1) = wd(1,1,1)
        wd(2,2,2) = (0.5*(1+costheta))**2
        wd(2,2,1) = -0.5*(1+costheta)*sintheta
        wd(2,2,0) = 0.25*sqrt(6.)*sintheta**2
        wd(2,2,-1) = -0.5*(1-costheta)*sintheta
        wd(2,2,-2) = (0.5*(1-costheta))**2
        wd(2,1,2) = -wd(2,2,1)
        wd(2,1,1) = 0.5*(1+costheta)*(2*costheta-1)
        wd(2,1,0) = -sqrt(3./2.)*sintheta*costheta
        wd(2,1,-1) = 0.5*(1-costheta)*(2*costheta+1)
        wd(2,1,-2) = wd(2,2,-1)
        wd(2,0,2) = wd(2,2,0)
        wd(2,0,1) = -wd(2,1,0)
        wd(2,0,0) = 1.5*costheta**2-0.5
        wd(2,0,-1) = wd(2,1,0)
        wd(2,0,-2) = wd(2,2,0)
        wd(2,-1,2) = -wd(2,2,-1)
        wd(2,-1,1) = wd(2,1,-1)
        wd(2,-1,0) = -wd(2,1,0)
        wd(2,-1,-1) = wd(2,1,1)
        wd(2,-1,-2) = wd(2,2,1)
        wd(2,-2,2) = wd(2,2,-2)
        wd(2,-2,1) = wd(2,-1,2)
        wd(2,-2,0) = wd(2,0,2)
        wd(2,-2,-1) = wd(2,1,2)
        wd(2,-2,-2) = wd(2,2,2)
      endif
c------------------------------------------------- kinematical variables

      mw=mass(kw)
      mz=mass(kz)
      if (inich.or.pch) then
        pp=p
        mp1=mass(kp1)
        mp2=mass(kp2)
        mp12=mp1**2
        mp22=mp2**2
        pp2=pp**2
        if (mp1.eq.mp2) then
          e1 = sqrt(mp12 + pp2)
          e2 = e1
          s = 4*e1**2
          ppl = pp
          pmi = 0.0d0
          epl = e1
          emi = mp1
        else
          e1=sqrt(mp12+pp2)
          e2=sqrt(mp22+pp2)
          s=mp12+mp22+2*e1*e2+2*pp2
          ppl=sqrt(max((e1*e2-mp1*mp2+pp2)/2.0d0,0.0d0))
          pmi=sign(1.0d0,mp1-mp2)*
     &      sqrt(max((e1*e2-mp1*mp2-pp2)/2.0d0,0.0d0))
          epl=sqrt(max((e1*e2+mp1*mp2+pp2)/2.0d0,0.0d0))
          emi=sqrt(max((e1*e2+mp1*mp2-pp2)/2.0d0,0.0d0))
        endif
      endif

      if (inich.or.finch.or.pch) then
        mp3 = mass(kp3)
        mp4 = mass(kp4)
        mp32 = mp3**2
        mp42 = mp4**2
        kk2 = (s-(mp3+mp4)**2)*(s-(mp3-mp4)**2)/(4*s)
        if (kk2.gt.0.0d0) then
          kk = sqrt(kk2)
          if (mp3.eq.mp4) then
            e3 = sqrt(mp3**2 + kk2)
            e4 = e3
            kpl = kk
            kmi = 0.0d0
            efpl = e3
            efmi = mp3
          else
            e3 = sqrt(mp3**2 + kk2)
            e4 = sqrt(mp4**2 + kk2)
            kpl = sqrt((e3*e4-mp3*mp4+kk2)/2.0d0)
            kmi = sign(1.0d0,mp3-mp4)*
     &        sqrt(max((e3*e4-mp3*mp4-kk2)/2.0d0,0.0d0))
            efpl = sqrt((e3*e4+mp3*mp4+kk2)/2.0d0)
            efmi = sqrt(max((e3*e4+mp3*mp4-kk2)/2.0d0,0.0d0))
          endif
        else
          kk=0.0d0
          e3=mp3
          e4=mp4
          kpl=0.0d0
          kmi=0.0d0
          efpl=sqrt(mp3*mp4)
          efmi=sqrt(mp3*mp4)
        endif
      endif

      t = mp1**2+mp3**2-2*e1*e3+2*pp*kk*costheta
      u = mp1**2+mp4**2-2*e1*e4-2*pp*kk*costheta
      mk = mass(kpk)

      end









