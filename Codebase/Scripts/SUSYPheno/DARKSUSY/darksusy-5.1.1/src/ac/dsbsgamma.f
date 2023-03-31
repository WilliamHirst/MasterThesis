      subroutine dsbsgamma(ratio,flag) 
c_______________________________________________________________________
c  b -> s + gamma branching ratio
c  input:
c    flag : 0 no qcd correction -- 1 qcd corrections
c  output:
c    ratio : branching ratio
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c     28-nov-94 formulas in bertolini et al, nucl phys b353 (1991) 591
c     modified: joakim edsjo, 2000-09-03, vertices from dsvertx.f
c     correctly implemented
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      real*8 ratio
      integer flag
      real*8 dsbsgf1,dsbsgf2,dsbsgf3,dsbsgf4,eu,ed,xtw,xth,xcu,xgd,xnd
      real*8 eta,lambda,rho,f,r,gfermi2b,alphemb,alphwb,alph3b,
     &     gexp,gbsg,gbcenu,dsabsq,pi
      parameter (pi=3.141592653589793238d0)
      complex*16 aw,ah,ac(2),ag,an(4),aa,a0
      integer i,k,c,j,g
      complex*16 coeff,temp1,temp2
      complex*16 glscf,grscf,gjll,gjlr,gjrl,gjrr
      complex*16 clsuchipb(6,2),clsuchips(6,2),crsuchipb(6,2)
      complex*16 clsdgluinb(6),clsdgluins(6),crsdgluinb(6)
      complex*16 clsdchi0b(6,4),clsdchi0s(6,4),crsdchi0b(6,4)
      real*8 alph3w,g3strow,alphemw,yuka(50),masstw,massbw,
     & g2weakw,gyweakw,sinbeta,cosbeta,aux

      logical first
      data first/.true./
      save first

      if (first) then
        write(*,*) '*************************************************'
        write(*,*) 'WARNING in dsbsgamma: This routine is obsolete.'
        write(*,*) 'It is only included here for testing purposes.'
        write(*,*) 'You should use the routine dsbsgammafull instead.'
        write(*,*) 'This warning is only printed once.'
        write(*,*) '*************************************************'
        first=.false.
      endif


      ratio = 0.d0
      if (dsabsq(ckm(2,3)).eq.0.d0) return

c at the mz scale
      alph3w=0.118d0
      g3strow=sqrt(4.d0*pi*alph3w)
      alphemw = 0.0078186083d0 ! =1/127.9
      g2weakw=sqrt(4.0d0*pi*alphem/s2thw)
      gyweakw=g2weakw*sqrt(s2thw)/sqrt(1.0d0-s2thw)
      
      massbw=2.83d0 ! Kuhn & Steinhauser
      masstw=167.d0 ! mt(mt) with pole mass = 175 and 1-loop relation
      cosbeta = 1.0d0 / sqrt(1.0d0+tanbe*tanbe)
      sinbeta = tanbe * cosbeta
      aux = g2weakw/dsqrt(2.d0)/mass(kw)
      do i=1,2
         yuka(kqu(i)) = aux*mass(kqu(i))/sinbeta
         yuka(kqd(i)) = aux*mass(kqd(i))/cosbeta
      enddo
      yuka(kqu(3)) = aux*masstw/sinbeta
      yuka(kqd(3)) = aux*massbw/cosbeta


c at the b mass scale: 
      alphemb = 1.d0/137.d0
      alphwb = alphemb/s2thw
      eta = 1.d0/0.548d0        ! eta = alph3(mb)/alph3(mw)
      alph3b = alph3*eta
      gfermi2b = 2.d0*(4.d0*pi*alphwb)**2/mass(kw)**4/64.d0

      eu = echarg(ku)
      ed = echarg(kd)


c------ up-squark chargino bottom or strange ------
      do c=1,2
         do k=1,6
            g=3
            glscf = dcmplx(0.d0,0.d0)
            grscf = dcmplx(0.d0,0.d0)
            do j=1,3
              glscf = glscf -
     &          (conjg(chavmx(c,2))*squrmx(k,j)*yuka(kqu(j))-
     &           g2weakw*conjg(chavmx(c,1))*squlmx(k,j))*ckm(j,g)
              grscf = grscf -
     &          chaumx(c,2)*squlmx(k,j)*ckm(j,g)*yuka(kqd(g))
            enddo
            clsuchipb(k,c)=glscf
            crsuchipb(k,c)=grscf
            g=2
            glscf = dcmplx(0.d0,0.d0)
            grscf = dcmplx(0.d0,0.d0)
            do j=1,3
              glscf = glscf -
     &          (conjg(chavmx(c,2))*squrmx(k,j)*yuka(kqu(j))-
     &           g2weakw*conjg(chavmx(c,1))*squlmx(k,j))*ckm(j,g)
              grscf = grscf -
     &          chaumx(c,2)*squlmx(k,j)*ckm(j,g)*yuka(kqd(g))
            enddo
            clsuchips(k,c) = glscf
         enddo
      enddo

c------ down-squark gluino bottom or strange ------

      do k=1,6
         g=3
         clsdgluinb(k)= -sqrt(2.d0)*g3strow*sqdlmx(k,g)
         crsdgluinb(k)= +sqrt(2.d0)*g3strow*sqdrmx(k,g)
         g=2
         clsdgluins(k)= -sqrt(2.d0)*g3strow*sqdlmx(k,g)
      enddo


c------ down-squark neutralino bottom or strange ------

      do j=1,4
         do k=1,6
            g=3
            gjll=-(-g2weakw*neunmx(j,2)+
     &       gyweakw*neunmx(j,1)/3.0d0)/sqrt(2.0d0)
            gjlr=-yuka(kqd(g))*neunmx(j,3)
            gjrl=-yuka(kqd(g))*neunmx(j,3)
            gjrr=sqrt(2.0d0)*echarg(kd)*gyweakw*neunmx(j,1)
            clsdchi0b(k,j)=
     &           conjg(gjll)*sqdlmx(k,g)+conjg(gjrl)*sqdrmx(k,g)
            crsdchi0b(k,j)=
     &           gjlr*sqdlmx(k,g)+gjrr*sqdrmx(k,g)
            g=2
            gjll=-(-g2weakw*neunmx(j,2)+
     &       gyweakw*neunmx(j,1)/3.0d0)/sqrt(2.0d0)
            gjlr=-yuka(kqd(g))*neunmx(j,3)
            gjrl=-yuka(kqd(g))*neunmx(j,3)
            gjrr=sqrt(2.0d0)*echarg(kd)*gyweakw*neunmx(j,1)
            clsdchi0s(k,j)=
     &           conjg(gjll)*sqdlmx(k,g)+conjg(gjrl)*sqdrmx(k,g)
        enddo
      enddo   

c-------------------------------------------------------- b -> c e nubar
      gexp = 0.107d0
      r = mass(kc)/mass(kb)
      rho = 1.d0-8.d0*r*r+8.d0*r**6-r**8-24.d0*r**4*dlog(r)
      f = 2.41d0
      if (flag.ne.0) then
         lambda = 1.d0-(2.d0/3.d0)*f*alph3*eta/pi
      else
         lambda = 1.d0
      endif
      lambda = 1.d0-(2.d0/3.d0)*f*alph3b/pi
      gbcenu = gfermi2b/(192.d0*pi**3)*rho*lambda* ! without mb**5
     &     (ckm(2,3)*conjg(ckm(2,3)))

      coeff = (alphwb/2.d0)*dsqrt(alphemb/pi) * 1/mass(kw)**2 *
     &     conjg(ckm(3,2))*ckm(3,3)

c-------------------------------------------------------------------- w+
      xtw = (mass(kt)/mass(kw))**2
      aw = (alphwb/4.d0)*dsqrt(alphemb/pi) * 1.d0/mass(kw)**2 *
     &     conjg(ckm(3,2))*ckm(3,3) *
     &     3*xtw*(eu*dsbsgf1(xtw)+dsbsgf2(xtw))

c--------------------------------------------------------- charged higgs
      xth = (mass(kt)/mass(khc))**2
      ah = (alphwb/4.d0)*dsqrt(alphemb/pi) * 1.d0/mass(kw)**2 *
     &     conjg(ckm(3,2))*ckm(3,3) *
     &     xth*((eu*dsbsgf1(xth)+dsbsgf2(xth))/tanbe**2
     &         +(eu*dsbsgf3(xth)+dsbsgf4(xth)))

c-------------------------------------------------------------- chargino
   
      do c=1,2
         ac(c) = cmplx(0.d0,0.d0)
         do k=1,6
            xcu = (mass(kcha(c))/mass(ksqu(k)))**2
            temp1 =  - (alphwb/2.d0)*dsqrt(alphemb/pi) *
     &           1.d0/mass(ksqu(k))**2 * (1.d0/g2weakw**2) *
     &           (clsuchipb(k,c)*conjg(clsuchips(k,c))*
     &                    (dsbsgf1(xcu)+eu*dsbsgf2(xcu)))
            temp2 =  - (alphwb/2.d0)*dsqrt(alphemb/pi) *
     &           1.d0/mass(ksqu(k))**2 * (1.d0/g2weakw**2) *
     &           (crsuchipb(k,c)*conjg(clsuchips(k,c))*
     &                    (mass(kcha(c))/mass(kb))*
     &                    (dsbsgf3(xcu)+eu*dsbsgf4(xcu)))
            ac(c) = ac(c) + temp1 + temp2
         enddo
      enddo

c---------------------------------------------------------------- gluino

      ag = cmplx(0.d0,0.d0)
      do k=1,6
         xgd = (mass(kgluin)/mass(ksqd(k)))**2
         ag = ag - alph3b*dsqrt(alphemb/pi) * (4.d0/3.d0) * ed *
     &        1.d0/(mass(ksqd(k)))**2 * (1.d0/(2.d0*g3strow**2)) *
     &        (clsdgluinb(k)*conjg(clsdgluins(k))*
     &                             dsbsgf2(xgd)
     &        +crsdgluinb(k)*conjg(clsdgluins(k))*
     &                (mass(kgluin)/mass(kb))*
     &                             dsbsgf4(xgd))
      enddo

c----------------------------------------------------------- neutralino

      do j=1,4
         an(j) = cmplx(0.d0,0.d0)
         do k=1,6
            xnd = (mass(kn(j))/mass(ksqd(k)))**2
            an(j) = an(j) - (alphwb/2.d0)*dsqrt(alphemb/pi) * ed *
     &           1.d0/(mass(ksqd(k)))**2 * (1.d0/g2weakw**2) *
     &           (clsdchi0b(k,j)*conjg(clsdchi0s(k,j))*
     &                             dsbsgf2(xnd)
     &           +clsdchi0b(k,j)*conjg(clsdchi0s(k,j))*
     &                 (mass(kn(j))/mass(kb))*
     &                             dsbsgf4(xnd))
         enddo
      enddo

c---------------------------------------------------------- b -> s gamma
      aa = aw+ah+ac(1)+ac(2)+ag+an(1)+an(2)+an(3)+an(4)
      if (prtlevel.ge.2) then
         write (*,*) 'bsgamma: aw=',aw/coeff
         write (*,*) '       : ah=',ah/coeff
         write (*,*) '       : ac(1)=',ac(1)/coeff
         write (*,*) '       : ac(2)=',ac(2)/coeff
         write (*,*) '       : ag=',ag/coeff
         write (*,*) '       : an(1)=',an(1)/coeff
         write (*,*) '       : an(2)=',an(2)/coeff
         write (*,*) '       : an(3)=',an(3)/coeff
         write (*,*) '       : an(4)=',an(4)/coeff
      endif
      if (flag.ne.0) then
         a0 = (alphwb/4.d0)*dsqrt(alphemb/pi) * 1.d0/mass(kw)**2 *
     &     conjg(ckm(3,2))*ckm(3,3)
         aa = eta**(-16.d0/23.d0) * ( aa + a0 *
     &        (116.d0/135.d0*(eta**(10.d0/23.d0)-1.d0)
     &        +58.d0/189.d0*(eta**(28.d0/35.d0)-1.d0)) )
      endif
      gbsg = 1.d0/(16.d0*pi) * (aa*conjg(aa)) ! without mb**5

c---------------------- normalization of b -> s gamma to b -> c e nubar
      ratio = gexp * (gbsg/gbcenu)

      end
