      real*8 function dsrmq1loop(mscale,kpart)
      implicit none
      include 'dsmssm.h'
      real*8 mscale,dsralph31loop
      integer kpart
      real*8 mcmb,mcmt,mbmt
c
      if(kpart.eq.kc) then
        if(mscale.le.mcmc) then
          dsrmq1loop=mcmc 
        elseif(mscale.le.mbmb.and.mscale.gt.mcmc) then
          dsrmq1loop=mcmc 
     &      *(dsralph31loop(mscale)/dsralph31loop(mcmc))**(12.d0/25.d0)
        elseif(mscale.le.mtmt.and.mscale.gt.mbmb) then
          mcmb=mcmc 
     &      *(dsralph31loop(mbmb)/dsralph31loop(mcmc))**(12.d0/25.d0)
          dsrmq1loop=mcmb
     &      *(dsralph31loop(mscale)/dsralph31loop(mbmb))**(12.d0/23.d0)
        elseif(mscale.gt.mtmt) then
          mcmt=mcmc 
     &      *(dsralph31loop(mbmb)/dsralph31loop(mcmc))**(12.d0/25.d0)
     &      *(dsralph31loop(mtmt)/dsralph31loop(mbmb))**(12.d0/23.d0)
          dsrmq1loop=mcmt
     &      *(dsralph31loop(mscale)/dsralph31loop(mtmt))**(12.d0/21.d0)
        endif
      elseif(kpart.eq.kb) then
        if(mscale.le.mbmb) then
          dsrmq1loop=mbmb
        elseif(mscale.le.mtmt.and.mscale.gt.mbmb) then
          dsrmq1loop=mbmb
     &      *(dsralph31loop(mscale)/dsralph31loop(mbmb))**(12.d0/23.d0)
        elseif(mscale.gt.mtmt) then
          mbmt=mbmb
     &      *(dsralph31loop(mtmt)/dsralph31loop(mbmb))**(12.d0/23.d0)
          dsrmq1loop=mbmt
     &      *(dsralph31loop(mscale)/dsralph31loop(mtmt))**(12.d0/21.d0)
        endif
      elseif(kpart.eq.kt) then
        if(mscale.le.mtmt) then
          dsrmq1loop=mtmt
        else
          dsrmq1loop=mtmt
     &      *(dsralph31loop(mscale)/dsralph31loop(mtmt))**(12.d0/23.d0)
        endif
      else  
        write(*,*) 'dsrmq1loop called for wrong particle'
        write(*,*) 'particle =  ',pname(kpart)
        write(*,*) 'rather than c, b or t quark'
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end
