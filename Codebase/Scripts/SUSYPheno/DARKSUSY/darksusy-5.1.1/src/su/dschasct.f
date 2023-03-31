      subroutine dschasct
c_______________________________________________________________________
c  chargino masses and mixings.
c  called by susyin or mrkin.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c     940407 correction to chargino mixing
c     990724 drop v1,v2 (pg)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      real*8 auxa,auxb,auxh,auxf,e1,e2,tp,tm,cp,cm,sp,sm,
     & ap,am,bp,bm,dp,dm,temp
      complex*16 chawmx(2,2)
      auxh=sqrt(2.d0)*mass(kw)*cosbe
      auxf=sqrt(2.d0)*mass(kw)*sinbe
      if (prtlevel.ge.2) then
         write (*, '('' chargino mass matrix ='',2(f10.2))') m2,auxf
         write (*, '(''                       '',2(f10.2))') auxh,mu
      endif
      ap=m2**2+auxh**2
      bp=m2*auxf+mu*auxh
      dp=mu**2+auxf**2
      am=m2**2+auxf**2
      bm=m2*auxh+mu*auxf
      dm=mu**2+auxh**2
c---------------------------------------------------------------- masses
      auxa=dp+ap
      auxb=sqrt(4*bp**2+(dp-ap)**2)
      mass(kcha(1)) = sqrt(0.5d0*(auxa+auxb))
      temp = 0.5d0*(auxa-auxb)
      if (temp.gt.0.d0) then
         mass(kcha(2)) = sqrt(temp)
      else
         mass(kcha(2)) = 0.d0
      endif
c---------------------------------------------------------- eigenvectors
      if(bp.eq.0.d0) then
         cp=0.d0
         sp=1.d0
      else
         auxa=dp-ap
         if(auxa.ge.0.d0) then
            tp=(auxa+auxb)/(2*bp)
         else
            tp=2*bp/(abs(auxa)+auxb)
         endif
         cp=1.d0/sqrt(1.d0+tp**2)
         sp=tp*cp
      endif
      if(bm.eq.0.d0) then
         cm=0.d0
         sm=1.d0
      else
         auxa=dm-am
         if(auxa.ge.0.d0) then
            tm=(auxa+auxb)/(2*bm)
         else
            tm=2*bm/(abs(auxa)+auxb)
         endif
         cm=1.d0/sqrt(1.d0+tm**2)
         sm=tm*cm
      endif
      chaumx(1,1)=dcmplx(cm,0.d0)
      chaumx(1,2)=dcmplx(sm,0.d0)
      chaumx(2,1)=dcmplx(-sm,0.d0)
      chaumx(2,2)=dcmplx(cm,0.d0)
      chawmx(1,1)=dcmplx(cp,0.d0)
      chawmx(1,2)=dcmplx(sp,0.d0)
      chawmx(2,1)=dcmplx(-sp,0.d0)
      chawmx(2,2)=dcmplx(cp,0.d0)
c----------------------------------------------------- mixing convention
      if (cm*(m2*cp+auxf*sp)+sm*(auxh*cp+mu*sp).lt.0.d0) then
         e1=-1
      else
         e1=+1
      endif
      if (-sm*(-m2*sp+auxf*cp)+cm*(-auxh*sp+mu*cp).lt.0.d0) then
         e2=-1
      else
         e2=+1
      endif

c               ! this is negative-mass convention
c      mass(kcha(1)) = e1*mass(kcha(1))
c      mass(kcha(2)) = e2*mass(kcha(2))
c      e1=1
c      e2=1

      chavmx(1,1)=dcmplx(e1*cp,0.d0)
      chavmx(1,2)=dcmplx(e1*sp,0.d0)
      chavmx(2,1)=dcmplx(-e2*sp,0.d0)
      chavmx(2,2)=dcmplx(e2*cp,0.d0)

      end
