      integer function dsantures(kp1,kp2,kp3,kp4,p)

c=======================================================================
c  determine if t- or u-resonances can occur for a given 2 ->2
c  scattering. if t_max>0 or u_max>0, then tures=1, otherwise tures=0
c  author: joakim edsjo, edsjo@physto.se
c  date: 97-09-17
c=======================================================================
      implicit none
      include 'dsmssm.h'
      integer kp1,kp2,kp3,kp4
      real*8 mp1,mp2,mp3,mp4,mp12,mp22,mp32,mp42,
     &  pp,e1,e2,e3,e4,s,tmax,umax,pp2,kk2
      real*8 p,kk

c------------------------------------------------- kinematical variables


c... initial state
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
      else
        e1=sqrt(mp12+pp2)
        e2=sqrt(mp22+pp2)
        s=mp12+mp22+2*e1*e2+2*pp2
      endif

c... final state
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
        else
          e3 = sqrt(mp3**2 + kk2)
          e4 = sqrt(mp4**2 + kk2)
        endif
      else
        kk=0.0d0
        e3=mp3
        e4=mp4
      endif

c      t = mp12+mp32-2*e1*e3+2*pp*kk*costheta
c      u = mp12+mp42-2*e1*e4-2*pp*kk*costheta

      tmax = mp12+mp32-2*e1*e3+2*pp*kk
      umax = mp12+mp42-2*e1*e4+2*pp*kk

      if(tmax.gt.0.0d0.or.umax.gt.0.0d0) then
        dsantures=1
      else
        dsantures=0
      endif

      end











