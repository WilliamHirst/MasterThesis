      function dsanwxint(p,a,b)
c_______________________________________________________________________
c  neutralino self-annihilation invariant rate integrated between
c  cos(theta)=a and cos(theta)=b.
c  input:
c    p - initial cm momentum (real) for lsp annihilations
c    integration limits a and b
c  called by wx.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c  modified slightly by joakim edsjo (edsjo@physto.se) 96-04-10
c=======================================================================
      implicit none
      include 'dsio.h'
      include 'dsidtag.h'
      real*8 dsanwxint,p,dsandwdcos,a,b,tot,eps,st,os,ost,del,sum,x
      integer jmax,it,l,j,nfcn,jdid
      external dsandwdcos
      parameter (eps=1.0d-2,jmax=20)  ! je change in eps
      dsanwxint=0.d0
      del=b-a
      ost=0.5*del*(dsandwdcos(p,a)+dsandwdcos(p,b))
      x=0.5*(b+a)
      st=0.5*(ost+del*dsandwdcos(p,x))
      os=(4.0*st-ost)/3.0
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=a+0.5*del
        sum=0.0
        do l=1,it
          sum=sum+dsandwdcos(p,x)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (abs(tot-os).le.eps*abs(os)) then
           dsanwxint=tot
           return
        endif
     	os=tot
        ost=st
      enddo
      call dswrite(0,1,'too many steps in dsanwxint')
      if (prtlevel.ge.2) then
         write(luerr,*) ' p = ',p
         write(luerr,*) ' a = ',a
         write(luerr,*) ' b = ',b
         write(luerr,*) ' model: ',idtag
c         if (prtlevel.ge.4) then
c            call dsrdwdwdcos(p,2**(jmax-1))
c            stop
c         endif
      endif
      end






