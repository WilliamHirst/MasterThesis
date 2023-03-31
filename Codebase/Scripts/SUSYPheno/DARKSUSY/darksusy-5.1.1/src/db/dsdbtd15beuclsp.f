**********************************************************************
*** function which makes a tabulation of dsdbtd15beucl as function
*** the distance between source and observer L, and neglecting the
*** weak dependence of dsdbtd15beucl over the vertical coordinate for 
*** the source zcl
*** 
*** for every td dsdbtd15beuclsp is tabulated on first call in L, with        
*** L between: 
***    Lmin=0.9d0*(r_0-pbrcy) and 
***    Lmax=1.1d0*dsqrt((r_0+pbrcy)**2+pbzcy**2)
*** and stored in spline tables. 
***
*** pbrcy and pbzcy in kpc are passed through a common block in 
*** dspb_clcom.h and should be set before the calling this routine.
***
*** there is no internal check to verify whether between to consecutive 
*** calls, with the same td, pbrcy and pbzcy, or halo parameters, or  
*** propagation parameters are changed. If this is done make sure, 
*** before calling this function, to reinitialize to zero the integer 
*** parameter clspset in the common block:
***
***      real*8 tdsetup
***      integer clspset
***      common/clspsetcom/tdsetup,clspset
*** 
*** input: L in kpc, td in GeV
*** output in 10^15 s kpc^-3
***
*** author: piero ullio (ullio@sissa.it)
*** date: 04-01-22
**********************************************************************


      real*8 function dsdbtd15beuclsp(L,td)
      implicit none
      include 'dspbcom.h'
      include 'dspbprivate.h'
      include 'dshmcom.h'
      real*8 L,td,dsdbtd15beucl
      integer solarmod,npoints,k
      real*8 Lmin,Lmax,thetaL,Lint,Lhor,result
      real*8 xclsp(100),yclsp(100),yclsp2(100),reclspnpoints
      common/dbclspcom/xclsp,yclsp,yclsp2,reclspnpoints
      real*8 tdsetup
      integer clspset
      common/dbclspsetcom/tdsetup,clspset
      reclspnpoints=60.d0
      npoints=int(reclspnpoints)
      if(clspset.ne.123456.or.dabs((td-tdsetup)/(td+tdsetup))
     &   .gt.1.d-8) then
        Lmin=0.9d0*(r_0-pbrcy)
        Lmax=1.1d0*dsqrt((r_0+pbrcy)**2+pbzcy**2)
        zcl=0.05d0
        thetaL=2.d0*datan(1.d0)
        deltathetacl=datan(1.d0)/100.d0
c        write(*,*) 'dsdbtd15beuclsp tabulation started for td = ',td
c     &    ,tdsetup  
        do k=1,npoints-2
          Lint=Lmin+(Lmax-Lmin)/dble(npoints-3)*(k-1)
          Lhor=dsqrt(Lint**2-zcl**2)
          rcl=dsqrt(Lhor**2+pbr0**2-2.d0*pbr0*Lhor*dcos(thetaL))
          thetacl=dasin(Lhor/rcl*dsin(thetaL))
          if(Lhor*dcos(thetaL).gt.pbr0.and.dcos(thetaL).gt.0.d0) then
            thetacl=4.d0*datan(1.d0)-thetacl
          endif
          xclsp(k+1)=Lint
          yclsp(k+1)=dlog(dsdbtd15beucl(td))
c          write(*,*) td,xclsp(k+1),yclsp(k+1) 
        enddo
c        write(*,*) 'dsdbtd15beuclsp tabulation is over'
        xclsp(1)=xclsp(2)*0.9999d0
        yclsp(1)=yclsp(2)*1.0001d0
        xclsp(npoints)=xclsp(npoints-1)*1.0001d0
        yclsp(npoints)=yclsp(npoints-1)*0.9999d0
        call dsspline(xclsp,yclsp,npoints,1.d31,1.d31,yclsp2)
        clspset=123456
        tdsetup=td
      endif
      if(L.ge.xclsp(1).and.L.le.xclsp(npoints)) then
        call dssplint(xclsp,yclsp,yclsp2,npoints,L,result)
        dsdbtd15beuclsp=dexp(result)
      else
        write(*,*) 'dsdbtd15beuclsp called out of the allowed range :' 
        write(*,*) 'L = ',L 
        write(*,*) 'L range betweeen ',xclsp(1),' and ',xclsp(npoints) 
        write(*,*) 'program stopped' 
        stop
      endif
 1000 format(60(1x,e14.8))
      return
      end
