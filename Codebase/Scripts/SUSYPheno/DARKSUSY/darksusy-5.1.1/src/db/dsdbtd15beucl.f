**********************************************************************
*** function that gives the dbar diffusion time per unit volume
*** (units of 10^15 sec kpc^-3) for a dbar point source located
*** at rcl, zcl, thetacl (in the cylidrical framework with the sun 
*** located at r=r_0, z=0 theta=0) and some small "angular width"
*** deltathetacl which makes the routine converge much faster
*** rcl, zcl, thetacl and deltathetacl are in the dspb_clcom.h common 
*** blocks and must be before calling this routine. rcl and zcl are in
*** kpc, thetacl and deltathetacl in rad.
*** numerical convergence gets slower for rcl->0 or zcl->0
***
*** it assumes the diffusion model in:
***   bergstrom, edsjo & ullio, ajp 526 (1999) 215
*** inputs:
***     td - dbar kinetic energy (gev)
***
*** the conversion from this source function to the local dbar flux 
*** is the same as for dsdbtd15beu(td), except that dsdbtd15beucl(td)
*** must be multiplied by: 
***     int dV (rho_cl(\vec{x}_cl)/rho0)**2
***   where the integral is over the volume of the clump,
***   rho_cl(\vec{x}_cl) is the density profile in the clump
***   and the local halo density rho0 is the normalization scale used 
***   everywhere 
*** 
*** author: piero ullio (ullio@sissa.it)
*** date: 04-01-22
**********************************************************************

      real*8 function dsdbtd15beucl(td)
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h'
      include 'dspbprivate.h'
      integer vindexs(0:nbesselk)
      real*8 vepsis(0:nbesselk),vtermk(0:nbesselk)
     &  ,vterms(0:nbesselk,1:nzerojk),verrors(0:nbesselk)
      real*8 td,pp,ee,rig,nucleon
      real*8 kpc
      parameter (kpc=3.08567802d0)
      real*8 dspbkdiff,dsdbsigmavdbar,piover2,dsbessjw,dspbaddterm
      integer jkcheck,jscheck,ii,jk,js
      real*8 sumk,sums,sumscheck,sumkcheck,beta,dspbkdiffm,relprec,ratio
      real*8 jzero(nzerojk)
      integer incr
      parameter(incr=99)
      real*8 nusk,Jklocal,Jkplus1squared
      integer k,s,kk,kprevck
      logical rerun
c
      nucleon=2.d0
      ee=nucleon*td+m_d
      pp=dsqrt(dabs(ee**2-m_d**2))
      rig=pp/nucleon
      if(pbpropmodel.eq.2) then
        dh=dspbkdiff(rig,1)
        dg=dspbkdiff(rig,2)
      elseif(pbpropmodel.eq.3) then
        beta=pp/ee
        dh=dspbkdiffm(beta,rig,1)
        dg=dspbkdiffm(beta,rig,2)
      else
        write(*,*) 'dsdbtd15beu called with wrong pbpropmodel'
        write(*,*) 'pbpropmodel = ',pbpropmodel
        write(*,*) 'program stopped'
        stop
      endif
      axsec=dsdbsigmavdbar(ee)
c diffusion constant in units of 10^27 cm^2 s^-1
c axec in mb*10^10 cm s^-1
      gamh=pbcvel/2.d0/dh/10.d0
c
      jkcheck=25
      jscheck=25
      relprec=1.d-6 !about 4 digits precision
      kprevck=0
      do k=0,nbesselk
        vindexs(k)=0
        vepsis(k)=relprec
      enddo
 40   continue
      sumk=0.d0
      do k=0,nbesselk
      if(vindexs(k).eq.0) then   
        sums=0.d0
      else
        sums=vtermk(k) 
      endif
      do s=1+vindexs(k),nzerojk
 10   nusk=storage(k,s,1)
      if(nusk.lt.1.d-16) then
        CALL DBZEJY(dble(k),s+incr,1,1.d-16,jzero)
        do ii=s,s+incr,1
          storage(k,ii,1)=jzero(ii)
          storage(k,ii,2)=dsbessjw(k,jzero(ii)*pbr0/pbrh)
          storage(k,ii,3)=(dsbessjw(k+1,jzero(ii)))**2
        enddo
        goto 10
      endif
      Jklocal=storage(k,s,2)
      Jkplus1squared=storage(k,s,3)
      vterms(k,s)=dspbaddterm(k,nusk,Jklocal,Jkplus1squared)
      sums=sums+vterms(k,s)
      if(s.gt.jscheck) then
        sumscheck=0.d0 
        do js=0,jscheck-1 
          sumscheck=sumscheck+vterms(k,s-js)
        enddo 
        if(dabs(sumscheck/sums).lt.vepsis(k)) then
          vindexs(k)=s
          verrors(k)=sumscheck
          vtermk(k)=sums
          goto 20
        endif  
      endif
      enddo
      write(*,*) 'storage matrix has index s too small, k, s = ',k,s
      write(*,*) 'program stopped'
      stop 
 20   continue
      sumk=sumk+vtermk(k)
      if(k.ge.jkcheck) then
        sumkcheck=0.d0 
        do jk=0,jkcheck 
          sumkcheck=sumkcheck+vtermk(k-jk)
        enddo
c        write(*,*) k,sumk,vtermk(k)
        if(dabs(sumkcheck/sumk).lt.relprec.and.k.gt.kprevck) then
          kprevck=k 
          goto 30
        endif  
      endif
      enddo
      write(*,*) 'storage matrix has index k too small, k = ',k
      write(*,*) 'program stopped'
      stop 
 30   continue
      rerun=.false.
      do kk=0,k-1
        ratio=dabs(sumkcheck)/dble(k)/dabs(verrors(kk))
        if(ratio.lt.1.d0) then
c          write(*,*) kk,ratio,vepsis(kk)
          vepsis(kk)=ratio*0.9d0*vepsis(kk)
          rerun=.true.
        endif
      enddo
      if(rerun) goto 40
      piover2=2.d0*datan(1.d0)
c dspbtd15beucl in 1/10^6 s/cm kpc^-2
      dsdbtd15beucl=sumk/pbrh**2/piover2
c dsdbtd15beucl in 10^15 s kpc^-3
      dsdbtd15beucl=dsdbtd15beucl*kpc
      return
      end

