**********************************************************************
*** function which computes the dbar diffusion time term corresponding 
*** to the axisymmetric diffuse source within a cylinder of radius
*** pbrcy and height 2* pbzcy. 
*** This routine assumes also that the Green function of
*** the diffusion equation dsdbtd15beuclsp(L,tp) does depend just
*** on kinetic energy tp and distance from the observer L, neglecting
*** a weak dependence on the cylindrical coordinate z.
*** For every tp, dsdbtd15beuclsp is tabulated on first call in L and
*** stored in spline tables. 
*** In this function and in dsdbtd15beuclsp, pbrcy and pbzcy in kpc are 
*** passed through a common block in dspbcom.h. There is no check 
*** in dsdbtd15beuclsp on whether, pbrcy and pbzcy which define the
*** interval of tabulation are changed. Check header dsdbtd15beuclsp
*** for more details on this and other warnings, and how to get the 
*** right implementation is such parameters are changed while running 
*** our own code
*** After the tabulation, the following integral is performed:
***
***  2 int_0^{pbzcy} int_0^{pbrcy} dr r  int_0^{2\pi} dphi
***     (dshmaxirho(r,zint)/rho0)^2 * dsdbtd15beuclsp(L(z,r,theta),tp)
***
*** The triple integral is splitted into a double integral on r and 
*** theta, this result is tabulated in z and then this integral is 
*** performed. The tabulation in z has at least 100 points on a 
*** regular grid between 0 and pbzcy (this is set by the parameter
*** incompnpoints in the dspbcompint1 function), however points are 
*** added as long as the values of the function in two nearest 
*** neighbour points differs more than 10% (this is set by the 
*** parameter reratio in the dspbcompint1 function)
***
*** input: scale in kpc, tp in GeV
*** output in 10^15 s
*** 
*** author: piero ullio (ullio@sissa.it)
*** date: 04-01-22
**********************************************************************

      real*8 function dsdbtd15comp(tp)
      implicit none
      include 'dspbcom.h'
      integer solarmod
      real*8 tp,dsdbcompint1,dummy
      real*8 up,low
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsdbcompint1
      real*8 tpint,zint,rint,rmaxint
      common/dsdbcompcom/tpint,zint,rint,rmaxint
      real*8 epsabs
      integer ndigits
      integer comp0set
      common/comp0setcom/epsabs,comp0set,ndigits
      limit=5000
      if(comp0set.ne.123456) then
        ndigits=4 
        epsabs=1.d-14     !numerical accuracy
        comp0set=123456
      endif  
      epsrel=10.d0**(-ndigits)
      rmaxint=pbrcy
      tpint=tp
      dummy=dsdbcompint1(0.d0)
      low=0.d0
      up=pbzcy
 10   call dqagse(dsdbcompint1,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 10
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      dsdbtd15comp=2.d0*result
 1000 format(a7,10(1x,e14.8))
      return
      end


      real*8 function dsdbcompint1(z)
      implicit none
      include 'dspbcom.h'
      real*8 z
      integer npoints,k,kk
      real*8 up,low,reratio
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsdbcompint2
      real*8 xcomp(1005),ycomp(1005),ycomp2(1005),recompnpoints
      common/compspcom/xcomp,ycomp,ycomp2,recompnpoints
      integer incompnpoints
      common/compspcom2/incompnpoints
      real*8 epsabs,tpsetup,rmaxsetup
      integer ndigits
      integer compset
      common/dbcompsetcom/epsabs,tpsetup,rmaxsetup,compset,ndigits
      real*8 tpint,zint,rint,rmaxint
      common/dsdbcompcom/tpint,zint,rint,rmaxint
      limit=5000
      npoints=int(recompnpoints)
      reratio=0.1d0
      if(compset.ne.123456.or.dabs((tpint-tpsetup)/(tpint+tpsetup))
     &   .gt.1.d-8.or.dabs((rmaxint-rmaxsetup)/(rmaxint+rmaxsetup))
     &   .gt.1.d-8) then
      ndigits=5 
      epsabs=1.d-14     !numerical accuracy
      epsrel=10.d0**(-ndigits)
c      write(*,*) 'dsdbcompint1 tabulation started for tp = ',tpint
c      write(*,*) '(tpsetup) = ',tpsetup
      incompnpoints=100
      npoints=incompnpoints
      do k=1,npoints-2
        zint=pbzcy/dble(npoints-3)*(k-1)
        low=0.d0
        up=rmaxint
 10     call dqagse(dsdbcompint2,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        if(dabs(result)/10.d0**ndigits.lt.epsabs) then
          epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
          goto 10
        elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
          epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
        endif
        xcomp(k+1)=zint
        ycomp(k+1)=result
      enddo
      xcomp(1)=xcomp(2)-0.001d0
      ycomp(1)=ycomp(2)
      xcomp(npoints)=xcomp(npoints-1)*1.0001d0
      ycomp(npoints)=ycomp(npoints-1)
 30   continue
      do k=1,npoints-1
        if(dabs(ycomp(k)-ycomp(k+1))
     &    .gt.reratio*max(dabs(ycomp(k)),dabs(ycomp(k+1)))) then
          zint=(xcomp(k)+xcomp(k+1))/2.d0
          low=0.d0
          up=rmaxint
 20       call dqagse(dsdbcompint2,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(result)/10.d0**ndigits.lt.epsabs) then
            epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
            goto 20
          elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
            epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
          endif
          do kk=npoints,k+1,-1
            xcomp(kk+1)=xcomp(kk)
            ycomp(kk+1)=ycomp(kk)
          enddo
          xcomp(k+1)=zint
          ycomp(k+1)=result
          npoints=npoints+1
          if(npoints.le.1005) then
            goto 30 
          else
            write(*,*) 'in dsdbcompint1 exceeded the maximum dim'
            write(*,*) 'allowed for vectors in the compspcom block'
            write(*,*) 'which is set equal to 1005' 
          endif
        endif  
      enddo
c      write(*,*) 'dsdbcompint1 tabulation is over with npoints =',
c     &  npoints
      recompnpoints=dble(npoints)
      call dsspline(xcomp,ycomp,npoints,1.d31,1.d31,ycomp2)
      compset=123456  
      tpsetup=tpint
      rmaxsetup=rmaxint
      endif
      if(z.ge.xcomp(1).and.z.le.xcomp(npoints)) then
        call dssplint(xcomp,ycomp,ycomp2,npoints,z,result) 
        dsdbcompint1=result
      else
        write(*,*) 'dsdbcompint1 called out of the allowed range :' 
        write(*,*) 'z = ',z
        write(*,*) 'z range betweeen ',xcomp(1),' and ',xcomp(npoints) 
        write(*,*) 'program stopped' 
        stop
      endif
      return
      end


      real*8 function dsdbcompint2(r)
      implicit none
      include 'dshmcom.h'
      real*8 tpint,zint,rint,rmaxint
      common/dsdbcompcom/tpint,zint,rint,rmaxint
      real*8 r,dshmaxirho
      real*8 up,low
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsdbcompint3
      integer ndigits
      integer comp2set
      real*8 epsabs
      common/dbcomp2setcom/epsabs,comp2set,ndigits
      limit=5000
      if(comp2set.ne.123456) then
        ndigits=5 
        epsabs=1.d-14     !numerical accuracy
        comp2set=123456
      endif
      epsrel=10.d0**(-ndigits)
      rint=r
      low=0.d0
      up=8.d0*datan(1.d0)
 10   call dqagseb(dsdbcompint3,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 10
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      dsdbcompint2=result*(dshmaxirho(rint,zint)/rho0)**2*rint
 1000 format(a7,10(1x,e14.8))
      return
      end


      real*8 function dsdbcompint3(theta)
      implicit none
      include 'dshmcom.h'
      real*8 theta,L,dsdbtd15beuclsp
      real*8 tpint,zint,rint,rmaxint
      common/dsdbcompcom/tpint,zint,rint,rmaxint
      L=dsqrt((r_0**2+rint**2-2.d0*r_0*rint*dcos(theta))+zint**2)
      dsdbcompint3=dsdbtd15beuclsp(L,tpint)
      return
      end

