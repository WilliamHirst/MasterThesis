**********************************************************************
*** function which makes a tabulation of dspbtd15beucl as function
*** the distance between source and observer L, and neglecting the
*** weak dependence of dspbtd15beucl over the vertical coordinate for 
*** the source zcl
*** 
*** for every tp dspbtd15beuclsp is tabulated on first call in L, with        
*** L between: 
***    Lmin=0.9d0*(r_0-pbrcy) and 
***    Lmax=1.1d0*dsqrt((r_0+pbrcy)**2+pbzcy**2)
*** and stored in spline tables. 
***
*** pbrcy and pbzcy in kpc are passed through a common block in 
*** dspbcom.h and should be set before the calling this routine.
***
*** there is no internal check to verify whether between to consecutive 
*** calls, with the same tp, pbrcy and pbzcy, or halo parameters, or  
*** propagation parameters are changed. If this is done make sure, 
*** before calling this function, to reinitialize to zero the integer 
*** parameter clspset in the common block:
***
***      real*8 tpsetup
***      integer clspset
***      common/clspsetcom/tpsetup,clspset
*** 
*** input: L in kpc, tp in GeV
*** output in 10^15 s kpc^-3
***
*** author: piero ullio (ullio@sissa.it)
*** date: 04-01-22
**********************************************************************


      real*8 function dspbtd15beuclsp(L,tp)
      implicit none
      include 'dspbcom.h'
      include 'dspbprivate.h'
      include 'dshmcom.h'
      real*8 L,tp,dspbtd15beucl
      integer solarmod,npoints,k
      real*8 Lmin,Lmax,thetaL,Lint,Lhor,result
      real*8 xclsp(100),yclsp(100),yclsp2(100),reclspnpoints
      common/clspcom/xclsp,yclsp,yclsp2,reclspnpoints
      real*8 tpsetup
      integer clspset
      common/clspsetcom/tpsetup,clspset
      reclspnpoints=60.d0
      npoints=int(reclspnpoints)
      if(clspset.ne.123456.or.dabs((tp-tpsetup)/(tp+tpsetup))
     &   .gt.1.d-8) then
        Lmin=0.9d0*(r_0-pbrcy)
        Lmax=1.1d0*dsqrt((r_0+pbrcy)**2+pbzcy**2)
        zcl=0.05d0
        thetaL=2.d0*datan(1.d0)
        deltathetacl=datan(1.d0)/100.d0
c        write(*,*) 'dspbtd15beuclsp tabulation started for tp = ',tp
c     &    ,tpsetup  
        do k=1,npoints-2
          Lint=Lmin+(Lmax-Lmin)/dble(npoints-3)*(k-1)
          Lhor=dsqrt(Lint**2-zcl**2)
          rcl=dsqrt(Lhor**2+pbr0**2-2.d0*pbr0*Lhor*dcos(thetaL))
          thetacl=dasin(Lhor/rcl*dsin(thetaL))
          if(Lhor*dcos(thetaL).gt.pbr0.and.dcos(thetaL).gt.0.d0) then
            thetacl=4.d0*datan(1.d0)-thetacl
          endif
          xclsp(k+1)=Lint
          yclsp(k+1)=dlog(dspbtd15beucl(tp))
c          write(*,*) tp,xclsp(k+1),yclsp(k+1) 
        enddo
c        write(*,*) 'dspbtd15beuclsp tabulation is over'
        xclsp(1)=xclsp(2)*0.9999d0
        yclsp(1)=yclsp(2)*1.0001d0
        xclsp(npoints)=xclsp(npoints-1)*1.0001d0
        yclsp(npoints)=yclsp(npoints-1)*0.9999d0
        call dsspline(xclsp,yclsp,npoints,1.d31,1.d31,yclsp2)
        clspset=123456
        tpsetup=tp
      endif
      if(L.ge.xclsp(1).and.L.le.xclsp(npoints)) then
        call dssplint(xclsp,yclsp,yclsp2,npoints,L,result)
        dspbtd15beuclsp=dexp(result)
      else
        write(*,*) 'dspbtd15beuclsp called out of the allowed range :' 
        write(*,*) 'L = ',L 
        write(*,*) 'L range betweeen ',xclsp(1),' and ',xclsp(npoints) 
        write(*,*) 'program stopped' 
        stop
      endif
 1000 format(60(1x,e14.8))
      return
      end
