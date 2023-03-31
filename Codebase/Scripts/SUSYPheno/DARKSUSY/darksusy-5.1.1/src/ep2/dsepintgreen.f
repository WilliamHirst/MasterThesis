      real*8 function dsepintgreen(DeltaV)
****************************************************************
***                                                          ***
*** function which gives the integral over volume of the     ***
*** positron green function times dsephaloterm (which is the ***
*** square of the density normalized to the local halo       ***
*** density for a smooth halo profile, i.e. for hclumpy = 1, ***
*** and the density probability of clumps normalized to the  ***
*** local halo density for a clumpy halo, i.e. hclumpy = 2)  ***
*** as a function of:                                        ***
***    DeltaV = 4 * K0 * tau_E * deltav in units of kpc**2   ***
*** i.e. of a given deltav = v(Eps)- v(Epsprime) this should ***
*** be called with:                                          ***
***    DeltaV=4.0d0*k27*tau16*deltav*10.d0/kpc**2            ***
*** where kpc=3.08567802d0 and the 10/kpc**2 converts from   ***
*** units of 10**43 cm**2 to units of kpc**2                 *** 
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-02-03                                         ***
****************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 DeltaV
      real*8 dv,rtilde2,Ltilde,rtildesun
      common/dsepzintcom/dv,rtilde2,Ltilde,rtildesun
ccc
      real*8 rint(0:20),low,up,par,recheck
      integer i,imax
      real*8 rint2(0:10),rstep,ref,dsephaloterm,ck
      integer ii,j,kk
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsepzint
      real*8 epsabs
      integer parrintset,ndigits
      common/parrintsetcom/epsabs,parrintset,ndigits
ccc
      if(parrintset.ne.123456) then
         ndigits=5
         epsabs=1.d-1          !numerical accuracy
         parrintset=123456
      endif
      epsrel=10.d0**(-ndigits)
      limit=5000
ccc
      dv=DeltaV
      Ltilde=l_h/dsqrt(dv)
      rtildesun=r_e/dsqrt(dv)
ccc
      rint(0)=0.d0
      rtildesun=r_e/dsqrt(dv)
      do i=1,5
        rint(i)=rtildesun/5.d0*i
      enddo
      imax=5
      rint2(0)=1.d-5/dsqrt(dv)
      rstep=rint2(0)
      ref=dsephaloterm(rstep*dsqrt(dv),0.d0)
      ii=0
      do i=1,5
        rstep=1.d-5/dsqrt(dv)*10.d0**i
        ck=dsephaloterm(rstep*dsqrt(dv),0.d0)
        if(ck.lt.1.d-1*ref) then
          ii=ii+1
          rint2(ii)=rstep
          ref=ck
        endif
      enddo
      do kk=0,ii
        do i=1,imax
          if(rint2(kk).gt.rint(i-1).and. rint2(kk).lt.rint(i)) then
            do j=imax,i,-1
              rint(j+1)= rint(j)
            enddo
            rint(i)=rint2(kk)
            imax=imax+1
            goto 12
          endif
        enddo
 12     continue
      enddo
      par=0.d0
      do i=imax,1,-1
      low=rint(i-1)
      up=rint(i)
 11   call dqagse(dsepzint,low,up,epsabs,epsrel,limit,
     &  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      recheck=par+result
      if(dabs(recheck)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(recheck)/10.d0**(ndigits+2) !too low accuracy
        goto 11
      elseif(dabs(recheck)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(recheck)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      par=par+result
      enddo
      do i=5,100
      low=rtildesun/5.d0*i
      up=rtildesun/5.d0*(i+1)
 13   call dqagse(dsepzint,low,up,epsabs,epsrel,limit,
     &  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      recheck=par+result
      if(dabs(recheck)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(recheck)/10.d0**(ndigits+2) !too low accuracy
        goto 13
      elseif(dabs(recheck)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(recheck)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      par=par+result
      if(dabs(result/par).lt.10.d0**(-ndigits)) goto 14 
      enddo
 14   continue
      dsepintgreen=2.d0*par/dsqrt(4.d0*datan(1.d0))
      return
      end



      real*8 function dsepzint(rtilde)
      implicit none
      include 'dsepcom.h'
      real*8 rtilde
      real*8 dv,rtilde2,Ltilde,rtildesun
      common/dsepzintcom/dv,rtilde2,Ltilde,rtildesun
ccc
      real*8 dsbessei0,zint(0:20),low,up,par
      integer i,imax
      real*8 zint2(0:10),zstep,ref,dsephaloterm,ck
      integer ii,j,kk
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsepzint2
      real*8 epsabs
      integer parzintset,ndigits
      common/parzintsetcom/epsabs,parzintset,ndigits
      if(parzintset.ne.123456) then
         ndigits=5
         epsabs=1.d-1          !numerical accuracy
         parzintset=123456
      endif
      epsrel=10.d0**(-ndigits)
      limit=5000
ccc
      rtilde2=rtilde
      zint(0)=0.d0
      do i=1,10
        zint(i)=1.d0/Ltilde*i
        if(zint(i).gt.1.d0) then
          imax=i
          zint(imax)=1.d0
          goto 10
        endif
      enddo
      imax=11
      zint(imax)=1.d0
 10   continue
      zint2(0)=1.d-5/l_h
      zstep=zint2(0)
      ref=dsephaloterm(rtilde2*dsqrt(dv),l_h*zstep)
      ii=0
      do i=1,5
        zstep=1.d-5/l_h*10.d0**i
        ck=dsephaloterm(rtilde2*dsqrt(dv),l_h*zstep)
        if(ck.lt.1.d-1*ref) then
          ii=ii+1
          zint2(ii)=zstep
          ref=ck
        endif
      enddo
      do kk=0,ii
        do i=1,imax
          if(zint2(kk).gt.zint(i-1).and. zint2(kk).lt.zint(i)) then
            do j=imax,i,-1
              zint(j+1)= zint(j)
            enddo
            zint(i)=zint2(kk)
            imax=imax+1
            goto 12
          endif
        enddo
 12     continue
      enddo
      par=0.d0
      do i=1,imax
      low=zint(i-1)
      up=zint(i)
 11   call dqagseb(dsepzint2,low,up,epsabs,epsrel,limit,
     &  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 11
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      par=par+result
c      write(*,*) 'int ',low,up,epsabs,result,par
      enddo
      dsepzint=2.d0*Ltilde*par*rtilde
     &  *dsbessei0(2.d0*rtilde*rtildesun)*dexp(-(rtilde-rtildesun)**2)
      return
      end
      


      real*8 function dsepzint2(zhat)
      implicit none
      include 'dsepcom.h'
      include 'dshmcom.h'
      real*8 dv,rtilde2,Ltilde,rtildesun
      common/dsepzintcom/dv,rtilde2,Ltilde,rtildesun
      real*8 zhat,dsephaloterm
      integer n
      real*8 sign,plus,minus,term,sum
      n=0
      sum=0.d0
      sign=1.d0
      plus=Ltilde**2*(-sign*zhat-2.d0*n)**2
      term=sign*exp(-plus)
      sum=sum+term
 10   n=n+1
      sign=sign*(-1.d0)
      plus=Ltilde**2*(sign*zhat+2.d0*n)**2
      minus=Ltilde**2*(sign*zhat-2.d0*n)**2
      term=sign*(exp(-plus)+exp(-minus))
      sum=sum+term
      if (abs(term).gt.1.d-5*abs(sum).or.n.lt.5) goto 10
ccc
      dsepzint2=sum*dsephaloterm(dsqrt(dv)*rtilde2,l_h*zhat)
      return
      end



      real*8 function dsephaloterm(radcoord,vertcoord)
      implicit none
      include 'dshmcom.h'
      real*8 radcoord,vertcoord,val,dshmaxirho,dshmaxiprob
ccc

      if(hclumpy.eq.1) then                     ! smooth halo
        val=dshmaxirho(radcoord,vertcoord)
        val=val/rho0
        val=val**2
      elseif(hclumpy.eq.2) then                 ! clumpy halo
        val=dshmaxiprob(radcoord,vertcoord)
        val=val/rho0
      else
        write(*,*) 'dsephaloterm called with wrong hclumpy=',hclumpy
        stop
      endif
      dsephaloterm=val
      return
      end
