c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.34, October 8, 2001, edsjo@physto.se)
c....Template file for dsaschicasea begins here

**************************************************************
*** SUBROUTINE dsaschicasea                                ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** sfermion(i) + neutralino(j)/chargino^+(j)              ***
*** -> gauge-boson + fermion                               ***
***                                                        ***
*** The sfermion must be the first mentioned               ***
*** particle (kp1) and the neutralino/chargino             ***
*** the other (kp2) -- not the opposite.                   ***
*** For the final state the gauge boson must be mentioned  ***
*** first (i.e. kp3) and next the fermion (kp4) --         ***
*** not the opposite.                                      ***                 
***                                                        ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-05                                         ***
*** QCD included: 02-03-21                                 ***
*** comment added by Piero Ullio, 02-07-01                 ***
*** added flavour changing charged exchange f the 2nd case ***
*** 3rd, 4th and the 5th case                              ***
*** added by Mia Schelke June 2006                         ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
*** bug with sfermion particle codes fixed by JE 090403    ***
**************************************************************

***** Note that it is assumed that coupling constants that do 
***** not exist have already been set to zero!!!!!
***** Thus, many of the coefficients defined in this code 
***** simplify when the diagrams contain sneutrinos 
***** or neutrinos.

      subroutine dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      complex*16 dsaschicaseaas,tmp1,tmp2
      real*8  msfer2,mchi2,mgb2,mfer2
      real*8  massorga
      complex*16  zftt1l,zftt1r,zftt2l,zftt2r
      complex*16  zfuu1l,zfuu1r,zfuu2l,zfuu2r
      complex*16  zfuu3l,zfuu3r,zfuu4l,zfuu4r,zfuu5l,zfuu5r
      complex*16  zfuu6l,zfuu6r,zfuu7l,zfuu7r,zfuu8l,zfuu8r
      complex*16  zfss,zfss1l,zfss1r,zfss2l,zfss2r 
      complex*16  zfss3l,zfss3r,zfss4l,zfss4r,zfss5l,zfss5r
      complex*16  zfss6l,zfss6r,zfss7l,zfss7r,zfss8l,zfss8r
      complex*16  zftu1l,zftu1r,zftu2l,zftu2r
      complex*16  zftu3l,zftu3r,zftu4l,zftu4r
      complex*16  zfts,zfts1l,zfts1r,zfts2l,zfts2r
      complex*16  zfts3l,zfts3r,zfts4l,zfts4r
      complex*16  zfsu,zfsu1l,zfsu1r,zfsu2l,zfsu2r 
      complex*16  zfsu3l,zfsu3r,zfsu4l,zfsu4r,zfsu5l,zfsu5r
      complex*16  zfsu6l,zfsu6r,zfsu7l,zfsu7r,zfsu8l,zfsu8r 
      real*8      ctt,cuu,css,ctu,cts,csu
      integer i,k,l
      integer kgb,kfer
      integer ksfert(6),kchiu(4)
      integer kfers(3)
      integer nsfert,nchiu,nfers
      real*8 pi
      parameter (pi=3.141592653589793238d0)


c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        par=0.d0
        return
      endif      
***** set particles in final state:
      kgb=kp3
      kfer=kp4

***** masses in final state:  
      mass3=mass(kgb)  
      mass4=mass(kfer)  
***** define the kinematic variables  
      call dsaskinset2
   
c....symmetry factor, for non-identical final state particles
      s34=1.d0

c....mass symbols used in form-code  
      msfer2=mass(kp1)**2
      mchi2=mass(kp2)**2
      mgb2=mass(kgb)**2
      mfer2=mass(kfer)**2
      s=Svar  
      t=Tvar  
      u=Uvar  


***** now set the color factors to 1.d0 for a slepton in
***** the initial state and to 3.d0 for a squark in 
***** the initial state
***** if the gauge boson in the final state is a gluon, then the 
***** color factor is changed -- see case nine and ten below 
***** the color factors enter the form generated expression for 
***** the squared amplitude

      ctt=1.d0
      cuu=1.d0
      css=1.d0
      ctu=1.d0
      cts=1.d0
      csu=1.d0

      if(iifam(1).ge.ivfam(ksu1).and.iifam(1).le.ivfam(ksb2)) then
        ctt=3.d0
        cuu=3.d0
        css=3.d0
        ctu=3.d0
        cts=3.d0
        csu=3.d0
      endif

***** initially setting the array of exchange sfermions 
***** to be empty
      do i=1,6
        ksfert(i)=0
      enddo	
      do i=1,3
        kfers(i)=0
      enddo

*****       
***** 
***** the first case    
***** sfermion(i) + neutralino(j) -> Z + fermion
***** the fermion must be a supersymmetric partner to the sfermion 
      if(icase.eq.1) then
        massorga=1.d0/mass(kgb)**2
***** in the completeness relation for the polarization vectors 
***** massorga is a factor that we have defined for the 
***** term that is inversely proportional to the mass squared
***** massorga is mass(kgb)**-2 for the massive gauge bosons
***** and 0 for the photon
***** now set particles in the intermediate states 
***** in the t-channel (1 or 2 sfermions):  
        nsfert=ncsfert
        do k=1,nsfert
          ksfert(k)=kcsfertn(k)
        enddo
***** in the u-channel (4 neutralinos):
        do k=1,4
          kchiu(k)=kn(k)
        enddo
        nchiu=4
***** in the s-channel (1 fermion):
        kfers(1)=kfer
        nfers=1
        goto 100
      endif
*****
***** the second case
***** down-type-sfermion(i) + neutralino(j) 
***** -> W^- + up-type-fermion
      if(icase.eq.2) then
        massorga=1.d0/mass(kgb)**2 
***** set particles in the intermediate states 
***** in the t-channel (1 or 2 sfermions):
        nsfert=ncsfertc
        do i=1,nsfert ! JE Correction 09-04-03
          ksfert(i)=kcsfertc(i)
        enddo
***** in the u-channel (2 charginos):
        do k=1,2
          kchiu(k)=kcha(k)
          kchiu(k+2)=0
        enddo
        nchiu=2
***** in the s-channel (1 fermion):
        kfers(1)=kcfers
        nfers=1
        goto 100
      endif
*****
***** the third case
***** down-type-sfermion(i) + chargino^+(j) 
***** -> Z + up-type-fermion
      if(icase.eq.3) then
        massorga=1.d0/mass(kgb)**2 
***** in the t-channel (2 sfermions):  
        nsfert=ncsfert
        do i=1,nsfert
          ksfert(i)=kcsfertn(i)
        enddo
***** in the u-channel (2 charginos):
        do k=1,2
          kchiu(k)=kcha(k)
          kchiu(k+2)=0
        enddo
        nchiu=2
***** in the s-channel (1 fermion):
        kfers(1)=kfer
        nfers=1
        goto 100
      endif
*****
***** the fourth case
***** up-type-sfermion(i) + chargino^+(j) 
***** -> W^+ + up-type-fermion
      if(icase.eq.4) then
        massorga=1.d0/mass(kgb)**2
***** in the t-channel, 6 squarks or 2 sleptons ********************** 
        nsfert=ncsfertc
        do i=1,nsfert
          ksfert(i)=kcsfertc(i)
        enddo
***** in the u-channel (4 neutralinos):
***** only when initial/final state (s)fermions are of same family
        if(iifam(1).eq.ivfam(kfer)) then
          do k=1,4
            kchiu(k)=kn(k)
          enddo
          nchiu=4
	else
          nchiu=0
        endif
***** no s-channel:
        nfers=0
        goto 100
      endif
*****
***** the fifth case
***** down-type-sfermion(i) + chargino^+(j) 
***** -> W^+ + down-type-fermion 
      if(icase.eq.5) then
        massorga=1.d0/mass(kgb)**2
***** no t-channel:  
        nsfert=0 
***** in the u-channel (4 neutralinos):
***** only when initial/final state (s)fermions are of same family
        if(iifam(1).eq.ivfam(kfer)) then
          do k=1,4
            kchiu(k)=kn(k)
          enddo
          nchiu=4
        else
          nchiu=0
        endif
***** in the s-channel (1 lepton or 3 quarks):
        nfers=ncferd
        do i=1,nfers
          kfers(i)=kcferd(i)
        enddo
        goto 100
      endif

      write(*,*) 'DS: dsaschicasea called with wrong icase : ',icase   
      write(*,*) 'DS: initial or final states : '  
      write(*,*)   pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  
      
 100  continue


c....specify the coefficients used in the form-code
c....first for the t-channel(sfermion exchange)
      zftt1l=dcmplx(0.d0,0.d0)
      zftt1r=dcmplx(0.d0,0.d0)
      zftt2l=dcmplx(0.d0,0.d0)
      zftt2r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0) then
      do k=1,nsfert
        do l=1,nsfert 
         zftt1l=zftt1l
     &        +dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,ksfert(k),kp1)*conjg(gl(kgb,ksfert(l),kp1))
     &        *gr(ksfert(k),kfer,kp2)*conjg(gr(ksfert(l),kfer,kp2))
         zftt1r=zftt1r
     &        +dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,ksfert(k),kp1)*conjg(gl(kgb,ksfert(l),kp1))
     &        *gl(ksfert(k),kfer,kp2)*conjg(gl(ksfert(l),kfer,kp2))
         zftt2l=zftt2l+mass(kfer)*mass(kp2)
     &        *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,ksfert(k),kp1)*conjg(gl(kgb,ksfert(l),kp1))
     &        *gl(ksfert(k),kfer,kp2)*conjg(gr(ksfert(l),kfer,kp2))
         zftt2r=zftt2r+mass(kfer)*mass(kp2)
     &        *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,ksfert(k),kp1)*conjg(gl(kgb,ksfert(l),kp1))
     &        *gr(ksfert(k),kfer,kp2)*conjg(gl(ksfert(l),kfer,kp2))
        enddo
      enddo
      endif



c....then the u-channel(neutralino/chargino exchange)
      zfuu1l=dcmplx(0.d0,0.d0)
      zfuu1r=dcmplx(0.d0,0.d0)
      zfuu2l=dcmplx(0.d0,0.d0)
      zfuu2r=dcmplx(0.d0,0.d0)
      zfuu3l=dcmplx(0.d0,0.d0)
      zfuu3r=dcmplx(0.d0,0.d0)
      zfuu4l=dcmplx(0.d0,0.d0)
      zfuu4r=dcmplx(0.d0,0.d0)
      zfuu5l=dcmplx(0.d0,0.d0)
      zfuu5r=dcmplx(0.d0,0.d0)
      zfuu6l=dcmplx(0.d0,0.d0)
      zfuu6r=dcmplx(0.d0,0.d0)
      zfuu7l=dcmplx(0.d0,0.d0)
      zfuu7r=dcmplx(0.d0,0.d0)
      zfuu8l=dcmplx(0.d0,0.d0)
      zfuu8r=dcmplx(0.d0,0.d0)
      if(nchiu.gt.0) then
      do k=1,nchiu
        do l=1,nchiu
          zfuu1l=zfuu1l+dsasdepro(u,kchiu(k))
     &         *conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu1r=zfuu1r+dsasdepro(u,kchiu(k))
     &         *conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu2l=zfuu2l+mass(kp2)*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu2r=zfuu2r+mass(kp2)*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu3l=zfuu3l+mass(kchiu(k))*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu3r=zfuu3r+mass(kchiu(k))*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu4l=zfuu4l+mass(kchiu(k))*mass(kp2)
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu4r=zfuu4r+mass(kchiu(k))*mass(kp2)
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu5l=zfuu5l+mass(kfer)*mass(kchiu(k))
     &         *mass(kp2)*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu5r=zfuu5r+mass(kfer)*mass(kchiu(k))
     &         *mass(kp2)*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu6l=zfuu6l+mass(kfer)*mass(kchiu(k))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l))) 
          zfuu6r=zfuu6r+mass(kfer)*mass(kchiu(k))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu7l=zfuu7l+mass(kfer)*mass(kp2)
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu7r=zfuu7r+mass(kfer)*mass(kp2)
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
          zfuu8l=zfuu8l+mass(kfer)*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gl(kp1,kfer,kchiu(k))*gl(kgb,kchiu(k),kp2)
     &         *conjg(gl(kgb,kchiu(l),kp2))
     &         *conjg(gr(kp1,kfer,kchiu(l)))
          zfuu8r=zfuu8r+mass(kfer)*mass(kchiu(l))
     &         *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &         *gr(kp1,kfer,kchiu(k))*gr(kgb,kchiu(k),kp2)
     &         *conjg(gr(kgb,kchiu(l),kp2))
     &         *conjg(gl(kp1,kfer,kchiu(l)))
        enddo
      enddo
      endif       



c....then the s-channel(fermion-exchange)
      zfss=dcmplx(0.d0,0.d0)
      zfss1l=dcmplx(0.d0,0.d0)
      zfss1r=dcmplx(0.d0,0.d0)
      zfss2l=dcmplx(0.d0,0.d0)
      zfss2r=dcmplx(0.d0,0.d0) 
      zfss3l=dcmplx(0.d0,0.d0)
      zfss3r=dcmplx(0.d0,0.d0)
      zfss4l=dcmplx(0.d0,0.d0)
      zfss4r=dcmplx(0.d0,0.d0) 
      zfss5l=dcmplx(0.d0,0.d0)
      zfss5r=dcmplx(0.d0,0.d0)
      zfss6l=dcmplx(0.d0,0.d0)
      zfss6r=dcmplx(0.d0,0.d0) 
      zfss7l=dcmplx(0.d0,0.d0)
      zfss7r=dcmplx(0.d0,0.d0)
      zfss8l=dcmplx(0.d0,0.d0)
      zfss8r=dcmplx(0.d0,0.d0)
      if(nfers.gt.0) then
***** the expressions below has ben changed June 2006 
***** before we had only kfers now kfers(l), kfers(k)	  
      do k=1,nfers
        do l=1,nfers	  
         zfss=dsasdepro(s,kfers(k))*conjg(dsasdepro(s,kfers(l)))
         zfss1l=zfss1l+zfss*gl(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l))) 
         zfss1r=zfss1r+zfss*gr(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l))) 
         zfss2l=zfss2l+zfss*mass(kp2)*mass(kfers(l))
     &          *gl(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l))) 
         zfss2r=zfss2r+zfss*mass(kp2)*mass(kfers(l))
     &          *gr(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l)))  
         zfss3l=zfss3l+zfss*mass(kfers(k))*mass(kfers(l))
     &          *gl(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l)))  
         zfss3r=zfss3r+zfss*mass(kfers(k))*mass(kfers(l))
     &          *gr(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l))) 
         zfss4l=zfss4l+zfss*mass(kfers(k))*mass(kp2)
     &          *gl(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l))) 
         zfss4r=zfss4r+zfss*mass(kfers(k))*mass(kp2)
     &          *gr(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l)))  
         zfss5l=zfss5l+zfss*mass(kfer)*mass(kfers(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gr(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l)))  
         zfss5r=zfss5r+zfss*mass(kfer)*mass(kfers(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gl(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l)))
         zfss6l=zfss6l+zfss*mass(kfer)*mass(kfers(k))
     &          *gr(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l)))
         zfss6r=zfss6r+zfss*mass(kfer)*mass(kfers(k))
     &          *gl(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l))) 
         zfss7l=zfss7l+zfss*mass(kfer)*mass(kp2)
     &          *gr(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l))) 
         zfss7r=zfss7r+zfss*mass(kfer)*mass(kp2)
     &          *gl(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l))) 
         zfss8l=zfss8l+zfss*mass(kfer)*mass(kfers(l))
     &          *gr(kgb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(kgb,kfer,kfers(l))) 
         zfss8r=zfss8r+zfss*mass(kfer)*mass(kfers(l))
     &          *gl(kgb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(kgb,kfer,kfers(l))) 
        enddo	 
	  enddo
	  zfss=dcmplx(1.d0,0.d0)
      endif


c....then t-channel amplitude multiplied by
c....hermitian conjugated u-channel amplitude
      zftu1l=dcmplx(0.d0,0.d0)
      zftu1r=dcmplx(0.d0,0.d0)
      zftu2l=dcmplx(0.d0,0.d0)
      zftu2r=dcmplx(0.d0,0.d0)
      zftu3l=dcmplx(0.d0,0.d0)
      zftu3r=dcmplx(0.d0,0.d0)
      zftu4l=dcmplx(0.d0,0.d0)
      zftu4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nchiu.gt.0) then
      do k=1,nsfert
        do l=1,nchiu
         zftu1l=zftu1l+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *gl(kgb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &        *conjg(gr(kgb,kchiu(l),kp2))
     &        *conjg(gr(kp1,kfer,kchiu(l)))
         zftu1r=zftu1r+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *gl(kgb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &        *conjg(gl(kgb,kchiu(l),kp2))
     &        *conjg(gl(kp1,kfer,kchiu(l)))
         zftu2l=zftu2l+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *mass(kp2)*mass(kchiu(l))
     &        *gl(kgb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &        *conjg(gl(kgb,kchiu(l),kp2))
     &        *conjg(gr(kp1,kfer,kchiu(l)))
         zftu2r=zftu2r+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *mass(kp2)*mass(kchiu(l))
     &        *gl(kgb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &        *conjg(gr(kgb,kchiu(l),kp2))
     &        *conjg(gl(kp1,kfer,kchiu(l)))
         zftu3l=zftu3l+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *mass(kfer)*mass(kchiu(l))
     &        *gl(kgb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &        *conjg(gl(kgb,kchiu(l),kp2))
     &        *conjg(gr(kp1,kfer,kchiu(l))) 
         zftu3r=zftu3r+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *mass(kfer)*mass(kchiu(l))
     &        *gl(kgb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &        *conjg(gr(kgb,kchiu(l),kp2))
     &        *conjg(gl(kp1,kfer,kchiu(l)))
         zftu4l=zftu4l+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *mass(kfer)*mass(kp2)
     &        *gl(kgb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &        *conjg(gr(kgb,kchiu(l),kp2))
     &        *conjg(gr(kp1,kfer,kchiu(l))) 
         zftu4r=zftu4r+dsasdepro(t,ksfert(k))
     &        *conjg(dsasdepro(u,kchiu(l)))
     &        *mass(kfer)*mass(kp2)
     &        *gl(kgb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &        *conjg(gl(kgb,kchiu(l),kp2))
     &        *conjg(gl(kp1,kfer,kchiu(l)))
        enddo
      enddo
      endif

c....s-channel amplitude multiplied by
c....hermitian conjugated t-channel amplitude 
      zfts=dcmplx(0.d0,0.d0)
      zfts1l=dcmplx(0.d0,0.d0)
      zfts1r=dcmplx(0.d0,0.d0)
      zfts2l=dcmplx(0.d0,0.d0)
      zfts2r=dcmplx(0.d0,0.d0)
      zfts3l=dcmplx(0.d0,0.d0)
      zfts3r=dcmplx(0.d0,0.d0)
      zfts4l=dcmplx(0.d0,0.d0)
      zfts4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nfers.gt.0) then
*       zfts=dsasdepro(s,kfers)
       do k=1,nsfert
* below has been changed June 2006 to include kfers(l) 	   
	    do l=1,nfers
         zfts=dsasdepro(s,kfers(l))
         zfts1l=zfts1l+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *conjg(gl(kgb,ksfert(k),kp1))*gl(kgb,kfer,kfers(l))
     &          *gr(kp1,kfers(l),kp2)*conjg(gr(ksfert(k),kfer,kp2))
         zfts1r=zfts1r+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *conjg(gl(kgb,ksfert(k),kp1))*gr(kgb,kfer,kfers(l))
     &          *gl(kp1,kfers(l),kp2)*conjg(gl(ksfert(k),kfer,kp2))
         zfts2l=zfts2l+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *mass(kfers(l))*mass(kp2)	 
     &          *conjg(gl(kgb,ksfert(k),kp1))*gl(kgb,kfer,kfers(l))
     &          *gl(kp1,kfers(l),kp2)*conjg(gr(ksfert(k),kfer,kp2))
         zfts2r=zfts2r+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *mass(kfers(l))*mass(kp2)
     &          *conjg(gl(kgb,ksfert(k),kp1))*gr(kgb,kfer,kfers(l))
     &          *gr(kp1,kfers(l),kp2)*conjg(gl(ksfert(k),kfer,kp2))
          zfts3l=zfts3l+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *mass(kfer)*mass(kp2)
     &          *conjg(gl(kgb,ksfert(k),kp1))*gr(kgb,kfer,kfers(l))
     &          *gl(kp1,kfers(l),kp2)*conjg(gr(ksfert(k),kfer,kp2))
          zfts3r=zfts3r+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *mass(kfer)*mass(kp2)
     &          *conjg(gl(kgb,ksfert(k),kp1))*gl(kgb,kfer,kfers(l))
     &          *gr(kp1,kfers(l),kp2)*conjg(gl(ksfert(k),kfer,kp2))
          zfts4l=zfts4l+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *mass(kfer)*mass(kfers(l))
     &          *conjg(gl(kgb,ksfert(k),kp1))*gr(kgb,kfer,kfers(l))
     &          *gr(kp1,kfers(l),kp2)*conjg(gr(ksfert(k),kfer,kp2))
          zfts4r=zfts4r+zfts*conjg(dsasdepro(t,ksfert(k)))
     &          *mass(kfer)*mass(kfers(l))
     &          *conjg(gl(kgb,ksfert(k),kp1))*gl(kgb,kfer,kfers(l))
     &          *gl(kp1,kfers(l),kp2)*conjg(gl(ksfert(k),kfer,kp2))
        enddo
	   enddo
	   zfts=dcmplx(1.d0,0.d0)
      endif

c....s-channel amplitude multiplied by
c....hermitian conjugated u-channel amplitude
      zfsu=dcmplx(0.d0,0.d0)
      zfsu1l=dcmplx(0.d0,0.d0)
      zfsu1r=dcmplx(0.d0,0.d0)
      zfsu2l=dcmplx(0.d0,0.d0)
      zfsu2r=dcmplx(0.d0,0.d0) 
      zfsu3l=dcmplx(0.d0,0.d0)
      zfsu3r=dcmplx(0.d0,0.d0)
      zfsu4l=dcmplx(0.d0,0.d0)
      zfsu4r=dcmplx(0.d0,0.d0)
      zfsu5l=dcmplx(0.d0,0.d0)
      zfsu5r=dcmplx(0.d0,0.d0)
      zfsu6l=dcmplx(0.d0,0.d0)
      zfsu6r=dcmplx(0.d0,0.d0)
      zfsu7l=dcmplx(0.d0,0.d0)
      zfsu7r=dcmplx(0.d0,0.d0)
      zfsu8l=dcmplx(0.d0,0.d0)
      zfsu8r=dcmplx(0.d0,0.d0) 
      if(nfers.gt.0.and.nchiu.gt.0) then
*       zfsu=dsasdepro(s,kfers) 
       do k=1,nchiu 
* below changed June 2006 to allow for kfers(l)	   
	    do l=1,nfers
         zfsu=dsasdepro(s,kfers(l)) 
         zfsu1l=zfsu1l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *gl(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu1r=zfsu1r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *gr(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))
         zfsu2l=zfsu2l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kp2)*mass(kchiu(k))
     &           *gl(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu2r=zfsu2r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kp2)*mass(kchiu(k))
     &           *gr(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k))) 
         zfsu3l=zfsu3l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kchiu(k))
     &           *gl(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu3r=zfsu3r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kchiu(k))
     &           *gr(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))
         zfsu4l=zfsu4l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kp2)
     &           *gl(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu4r=zfsu4r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kp2)
     &           *gr(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))
         zfsu5l=zfsu5l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))*mass(kp2)*mass(kchiu(k))
     &           *gr(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu5r=zfsu5r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))*mass(kp2)*mass(kchiu(k))
     &           *gl(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))
         zfsu6l=zfsu6l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))
     &           *gr(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu6r=zfsu6r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))
     &           *gl(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))
         zfsu7l=zfsu7l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kp2)
     &           *gr(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu7r=zfsu7r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kp2)
     &           *gl(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))
         zfsu8l=zfsu8l+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kchiu(k))
     &           *gr(kgb,kfer,kfers(l))*gl(kp1,kfers(l),kp2)
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kfer,kchiu(k)))
         zfsu8r=zfsu8r+zfsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kchiu(k))
     &           *gl(kgb,kfer,kfers(l))*gr(kp1,kfers(l),kp2)
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kfer,kchiu(k)))  
	    enddo
       enddo
	   zfsu=dcmplx(1.d0,0.d0)
      endif

***** After the template file follows the form expression of the 
***** ('complex') amplitude squared: dsaschicaseaas
***** which is (M_tM_t^\dagger + M_sM_s^\dagger + M_uM_u^\dagger
***** + 2M_tM_u^\dagger + 2M_sM_t^\dagger + 2M_sM_u^\dagger)
***** As M_tM_u^\dagger+M_uM_t^\dagger=2Re(M_tM_u^\dagger) 
***** we have to take the real part of the interference terms,
***** in order to obtain the true amplitude squared.
***** This is most easily done by taking the real part of the
***** whole expression dsaschicaseaas.
***** This is done at the very end of the code.

c....Template file for dsaschicasea ends here










      tmp1=
     &  + msfer2*mchi2*mgb2*massorga*zfuu1l*cuu
     &  + msfer2*mchi2*mgb2*massorga*zfuu1r*cuu
     &  + msfer2*mchi2*mfer2*massorga*zfuu1l*cuu
     &  + msfer2*mchi2*mfer2*massorga*zfuu1r*cuu
     &  - 2.d0*msfer2*mchi2*mfer2*massorga*zfts*zfts1l*cts
     &  - 2.d0*msfer2*mchi2*mfer2*massorga*zfts*zfts1r*cts
     &  + 2.d0*msfer2*mchi2*mfer2*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*msfer2*mchi2*mfer2*massorga*zfsu*zfsu1r*csu
     &  - msfer2*mchi2*s*massorga*zfuu1l*cuu
     &  - msfer2*mchi2*s*massorga*zfuu1r*cuu
     &  + 2.d0*msfer2*mchi2*s*massorga*zfts*zfts1l*cts
     &  + 2.d0*msfer2*mchi2*s*massorga*zfts*zfts1r*cts
     &  - 2.d0*msfer2*mchi2*s*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*msfer2*mchi2*s*massorga*zfsu*zfsu1r*csu
     &  - 2.d0*msfer2*mchi2*t*massorga*zftt1l*ctt
     &  - 2.d0*msfer2*mchi2*t*massorga*zftt1r*ctt
     &  - 2.d0*msfer2*mchi2*t*massorga*zfuu1l*cuu
     &  - 2.d0*msfer2*mchi2*t*massorga*zfuu1r*cuu
     &  + 4.d0*msfer2*mchi2*t*massorga*zftu1l*ctu
     &  + 4.d0*msfer2*mchi2*t*massorga*zftu1r*ctu
     &  - 2.d0*msfer2*mchi2*zftt1l*ctt
     &  - 2.d0*msfer2*mchi2*zftt1r*ctt
     &  + 2.d0*msfer2*mchi2*zftu1l*ctu
     &  + 2.d0*msfer2*mchi2*zftu1r*ctu
     &  - 4.d0*msfer2*mchi2*zfts*zfts1l*cts
     &  - 4.d0*msfer2*mchi2*zfts*zfts1r*cts
     &  + 4.d0*msfer2*mchi2*zfsu*zfsu1l*csu
     &  + 4.d0*msfer2*mchi2*zfsu*zfsu1r*csu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfuu1l*cuu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfuu1r*cuu
     &  + msfer2*mgb2*mfer2*massorga*zfss*zfss1l*css
     &  + msfer2*mgb2*mfer2*massorga*zfss*zfss1r*css
     &  - 2.d0*msfer2*mgb2*mfer2*massorga*zftu1l*ctu
     &  - 2.d0*msfer2*mgb2*mfer2*massorga*zftu1r*ctu
     &  - 2.d0*msfer2*mgb2*mfer2*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*msfer2*mgb2*mfer2*massorga*zfsu*zfsu1r*csu
     &  - 2.d0*msfer2*mgb2*t*massorga*zfuu1l*cuu
     &  - 2.d0*msfer2*mgb2*t*massorga*zfuu1r*cuu
     &  + 2.d0*msfer2*mgb2*t*massorga*zftu1l*ctu
     &  + 2.d0*msfer2*mgb2*t*massorga*zftu1r*ctu
     &  + msfer2*mgb2*massorga*zfuu3l*cuu
     &  + msfer2*mgb2*massorga*zfuu3r*cuu
     &  + msfer2*mgb2*massorga*zfuu6l*cuu
     &  + msfer2*mgb2*massorga*zfuu6r*cuu
     &  + 2.d0*msfer2*mgb2*massorga*zfuu7l*cuu
     &  + 2.d0*msfer2*mgb2*massorga*zfuu7r*cuu
     &  + msfer2*mgb2*massorga*zfuu8l*cuu
     &  + msfer2*mgb2*massorga*zfuu8r*cuu
     &  + msfer2*mgb2*massorga*zfss*zfss3l*css
     &  + msfer2*mgb2*massorga*zfss*zfss3r*css
     &  - msfer2*mgb2*massorga*zfss*zfss6l*css
     &  - msfer2*mgb2*massorga*zfss*zfss6r*css
     &  - msfer2*mgb2*massorga*zfss*zfss8l*css
     &  - msfer2*mgb2*massorga*zfss*zfss8r*css
     &  - 2.d0*msfer2*mgb2*massorga*zftu2l*ctu
     &  - 2.d0*msfer2*mgb2*massorga*zftu2r*ctu
     &  - 4.d0*msfer2*mgb2*massorga*zftu4l*ctu
     &  - 4.d0*msfer2*mgb2*massorga*zftu4r*ctu
     &  - 2.d0*msfer2*mgb2*massorga*zfts*zfts2l*cts
     &  - 2.d0*msfer2*mgb2*massorga*zfts*zfts2r*cts
     &  + 2.d0*msfer2*mgb2*massorga*zfts*zfts3l*cts
     &  + 2.d0*msfer2*mgb2*massorga*zfts*zfts3r*cts
     &  + 2.d0*msfer2*mgb2*massorga*zfsu*zfsu3l*csu
     &  + 2.d0*msfer2*mgb2*massorga*zfsu*zfsu3r*csu
     &  + 2.d0*msfer2*mgb2*massorga*zfsu*zfsu6l*csu
     &  + 2.d0*msfer2*mgb2*massorga*zfsu*zfsu6r*csu
     &  - 2.d0*msfer2*mgb2*massorga*zfsu*zfsu8l*csu
     &  - 2.d0*msfer2*mgb2*massorga*zfsu*zfsu8r*csu
     &  + 2.d0*msfer2*mgb2*zfuu1l*cuu
     &  + 2.d0*msfer2*mgb2*zfuu1r*cuu
     &  + 2.d0*msfer2*mgb2*zfss*zfss1l*css
     &  + 2.d0*msfer2*mgb2*zfss*zfss1r*css
     &  + 4.d0*msfer2*mgb2*zfsu*zfsu1l*csu
     &  + 4.d0*msfer2*mgb2*zfsu*zfsu1r*csu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfuu1l*cuu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfuu1r*cuu
     &  + msfer2*mfer2*s*massorga*zfss*zfss1l*css
     &  + msfer2*mfer2*s*massorga*zfss*zfss1r*css
     &  + 2.d0*msfer2*mfer2*s*massorga*zftu1l*ctu
     &  + 2.d0*msfer2*mfer2*s*massorga*zftu1r*ctu
     &  + 2.d0*msfer2*mfer2*s*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*msfer2*mfer2*s*massorga*zfsu*zfsu1r*csu
     &  - 2.d0*msfer2*mfer2*t*massorga*zftt1l*ctt
     &  - 2.d0*msfer2*mfer2*t*massorga*zftt1r*ctt
     &  - 4.d0*msfer2*mfer2*t*massorga*zfuu1l*cuu
     &  - 4.d0*msfer2*mfer2*t*massorga*zfuu1r*cuu
     &  + 6.d0*msfer2*mfer2*t*massorga*zftu1l*ctu
     &  + 6.d0*msfer2*mfer2*t*massorga*zftu1r*ctu
     &  - 2.d0*msfer2*mfer2*t*massorga*zfts*zfts1l*cts
     &  - 2.d0*msfer2*mfer2*t*massorga*zfts*zfts1r*cts
     &  + 2.d0*msfer2*mfer2*t*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*msfer2*mfer2*t*massorga*zfsu*zfsu1r*csu
     &  + msfer2*mfer2*massorga*zfuu3l*cuu
     &  + msfer2*mfer2*massorga*zfuu3r*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfuu6l*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfuu6r*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfuu8l*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfuu8r*cuu
     &  + msfer2*mfer2*massorga*zfss*zfss3l*css
     &  + msfer2*mfer2*massorga*zfss*zfss3r*css
     &  - 2.d0*msfer2*mfer2*massorga*zftu2l*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zftu2r*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zftu3l*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zftu3r*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zftu4l*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zftu4r*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zfts*zfts2l*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfts*zfts2r*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfts*zfts3l*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfts*zfts3r*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfts*zfts4l*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfts*zfts4r*cts
     &  + 2.d0*msfer2*mfer2*massorga*zfsu*zfsu3l*csu
     &  + 2.d0*msfer2*mfer2*massorga*zfsu*zfsu3r*csu
     &  + 4.d0*msfer2*mfer2*massorga*zfsu*zfsu6l*csu
     &  + 4.d0*msfer2*mfer2*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*msfer2*mfer2*massorga*zfsu*zfsu7l*csu
     &  + 2.d0*msfer2*mfer2*massorga*zfsu*zfsu7r*csu
     &  - 2.d0*msfer2*mfer2*zftt1l*ctt
     &  - 2.d0*msfer2*mfer2*zftt1r*ctt
     &  + 2.d0*msfer2*mfer2*zfuu1l*cuu
     &  + 2.d0*msfer2*mfer2*zfuu1r*cuu
     &  - 2.d0*msfer2*mfer2*zfss*zfss1l*css
     &  - 2.d0*msfer2*mfer2*zfss*zfss1r*css
     &  + 6.d0*msfer2*mfer2*zftu1l*ctu
     &  + 6.d0*msfer2*mfer2*zftu1r*ctu
     &  - 2.d0*msfer2*mfer2*zfts*zfts1l*cts
     &  - 2.d0*msfer2*mfer2*zfts*zfts1r*cts
     &  + 8.d0*msfer2*mfer2*zfsu*zfsu1l*csu
     &  + 8.d0*msfer2*mfer2*zfsu*zfsu1r*csu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfuu1l*cuu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfuu1r*cuu
     &  - msfer2*mfer2**2*massorga*zfss*zfss1l*css
     &  - msfer2*mfer2**2*massorga*zfss*zfss1r*css
     &  - 2.d0*msfer2*mfer2**2*massorga*zftu1l*ctu
     &  - 2.d0*msfer2*mfer2**2*massorga*zftu1r*ctu
     &  - 2.d0*msfer2*mfer2**2*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*msfer2*mfer2**2*massorga*zfsu*zfsu1r*csu
     &  + 2.d0*msfer2*s*t*massorga*zfuu1l*cuu
     &  + 2.d0*msfer2*s*t*massorga*zfuu1r*cuu
     &  - 2.d0*msfer2*s*t*massorga*zftu1l*ctu
     &  - 2.d0*msfer2*s*t*massorga*zftu1r*ctu
     &  - 2.d0*msfer2*s*t*massorga*zfts*zfts1l*cts
     &  - 2.d0*msfer2*s*t*massorga*zfts*zfts1r*cts
     &  + 2.d0*msfer2*s*t*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*msfer2*s*t*massorga*zfsu*zfsu1r*csu
     &  - msfer2*s*massorga*zfuu3l*cuu
     &  - msfer2*s*massorga*zfuu3r*cuu
     &  - 2.d0*msfer2*s*massorga*zfuu6l*cuu
     &  - 2.d0*msfer2*s*massorga*zfuu6r*cuu
     &  - 2.d0*msfer2*s*massorga*zfuu8l*cuu
     &  - 2.d0*msfer2*s*massorga*zfuu8r*cuu
     &  - msfer2*s*massorga*zfss*zfss3l*css
     &  - msfer2*s*massorga*zfss*zfss3r*css
     &  + 2.d0*msfer2*s*massorga*zftu2l*ctu
     &  + 2.d0*msfer2*s*massorga*zftu2r*ctu
     &  + 2.d0*msfer2*s*massorga*zftu3l*ctu
     &  + 2.d0*msfer2*s*massorga*zftu3r*ctu
     &  + 2.d0*msfer2*s*massorga*zftu4l*ctu
     &  + 2.d0*msfer2*s*massorga*zftu4r*ctu
     &  + 2.d0*msfer2*s*massorga*zfts*zfts2l*cts
     &  + 2.d0*msfer2*s*massorga*zfts*zfts2r*cts
     &  + 2.d0*msfer2*s*massorga*zfts*zfts3l*cts
     &  + 2.d0*msfer2*s*massorga*zfts*zfts3r*cts
     &  + 2.d0*msfer2*s*massorga*zfts*zfts4l*cts
     &  + 2.d0*msfer2*s*massorga*zfts*zfts4r*cts
     &  - 2.d0*msfer2*s*massorga*zfsu*zfsu3l*csu
     &  - 2.d0*msfer2*s*massorga*zfsu*zfsu3r*csu
     &  - 4.d0*msfer2*s*massorga*zfsu*zfsu6l*csu
     &  - 4.d0*msfer2*s*massorga*zfsu*zfsu6r*csu
     &  - 2.d0*msfer2*s*massorga*zfsu*zfsu7l*csu
     &  - 2.d0*msfer2*s*massorga*zfsu*zfsu7r*csu
     &  - 2.d0*msfer2*s*zfuu1l*cuu
     &  - 2.d0*msfer2*s*zfuu1r*cuu
     &  - 2.d0*msfer2*s*zfss*zfss1l*css
     &  - 2.d0*msfer2*s*zfss*zfss1r*css
     &  - 4.d0*msfer2*s*zfsu*zfsu1l*csu
     &  - 4.d0*msfer2*s*zfsu*zfsu1r*csu
     &  - 4.d0*msfer2*t*massorga*zftt2l*ctt
     &  - 4.d0*msfer2*t*massorga*zftt2r*ctt
     &  - 2.d0*msfer2*t*massorga*zfuu6l*cuu
     &  - 2.d0*msfer2*t*massorga*zfuu6r*cuu
     &  - 2.d0*msfer2*t*massorga*zfuu8l*cuu
     &  - 2.d0*msfer2*t*massorga*zfuu8r*cuu
     &  + 4.d0*msfer2*t*massorga*zftu3l*ctu
     &  + 4.d0*msfer2*t*massorga*zftu3r*ctu
     &  + 4.d0*msfer2*t*massorga*zftu4l*ctu
     &  + 4.d0*msfer2*t*massorga*zftu4r*ctu
     &  + 4.d0*msfer2*t*massorga*zfts*zfts4l*cts
     &  + 4.d0*msfer2*t*massorga*zfts*zfts4r*cts
     &  - 4.d0*msfer2*t*massorga*zfsu*zfsu6l*csu
     &  - 4.d0*msfer2*t*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*msfer2*t*zftt1l*ctt
     &  + 2.d0*msfer2*t*zftt1r*ctt
     &  - 2.d0*msfer2*t*zftu1l*ctu
     &  - 2.d0*msfer2*t*zftu1r*ctu
     &  + 4.d0*msfer2*t*zfts*zfts1l*cts
     &  + 4.d0*msfer2*t*zfts*zfts1r*cts
     &  - 4.d0*msfer2*t*zfsu*zfsu1l*csu
     &  - 4.d0*msfer2*t*zfsu*zfsu1r*csu
     &  + 2.d0*msfer2*t**2*massorga*zftt1l*ctt
     &  + 2.d0*msfer2*t**2*massorga*zftt1r*ctt
     &  + 2.d0*msfer2*t**2*massorga*zfuu1l*cuu
     &  + 2.d0*msfer2*t**2*massorga*zfuu1r*cuu
     &  - 4.d0*msfer2*t**2*massorga*zftu1l*ctu
     &  - 4.d0*msfer2*t**2*massorga*zftu1r*ctu
     &  - 4.d0*msfer2*zftt2l*ctt
     &  - 4.d0*msfer2*zftt2r*ctt
     &  + 2.d0*msfer2*zfuu6l*cuu
     &  + 2.d0*msfer2*zfuu6r*cuu
     &  - 8.d0*msfer2*zfuu7l*cuu
     &  - 8.d0*msfer2*zfuu7r*cuu
     &  + 2.d0*msfer2*zfuu8l*cuu
     &  + 2.d0*msfer2*zfuu8r*cuu
     &  + 4.d0*msfer2*zfss*zfss6l*css
     &  + 4.d0*msfer2*zfss*zfss6r*css
     &  + 4.d0*msfer2*zfss*zfss8l*css
     &  + 4.d0*msfer2*zfss*zfss8r*css
     &  + 2.d0*msfer2*zftu3l*ctu
     &  + 2.d0*msfer2*zftu3r*ctu
     &  + 6.d0*msfer2*zftu4l*ctu
     &  + 6.d0*msfer2*zftu4r*ctu
     &  - 4.d0*msfer2*zfts*zfts3l*cts
     &  - 4.d0*msfer2*zfts*zfts3r*cts
     &  + 2.d0*msfer2*zfts*zfts4l*cts
     &  + 2.d0*msfer2*zfts*zfts4r*cts
     &  + 4.d0*msfer2*zfsu*zfsu6l*csu
     &  + 4.d0*msfer2*zfsu*zfsu6r*csu
     &  - 4.d0*msfer2*zfsu*zfsu7l*csu
     &  - 4.d0*msfer2*zfsu*zfsu7r*csu
     &  + 8.d0*msfer2*zfsu*zfsu8l*csu
     &  + 8.d0*msfer2*zfsu*zfsu8r*csu
     &  + msfer2**2*mchi2*massorga*zftt1l*ctt
     &  + msfer2**2*mchi2*massorga*zftt1r*ctt
     &  + msfer2**2*mchi2*massorga*zfuu1l*cuu
     &  + msfer2**2*mchi2*massorga*zfuu1r*cuu
     &  - 2.d0*msfer2**2*mchi2*massorga*zftu1l*ctu
     &  - 2.d0*msfer2**2*mchi2*massorga*zftu1r*ctu
     &  + msfer2**2*mfer2*massorga*zftt1l*ctt
     &  + msfer2**2*mfer2*massorga*zftt1r*ctt
     &  + msfer2**2*mfer2*massorga*zfuu1l*cuu
     &  + msfer2**2*mfer2*massorga*zfuu1r*cuu
     &  - 2.d0*msfer2**2*mfer2*massorga*zftu1l*ctu
     &  - 2.d0*msfer2**2*mfer2*massorga*zftu1r*ctu
     &  + 2.d0*msfer2**2*mfer2*massorga*zfts*zfts1l*cts
     &  + 2.d0*msfer2**2*mfer2*massorga*zfts*zfts1r*cts
     &  - 2.d0*msfer2**2*mfer2*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*msfer2**2*mfer2*massorga*zfsu*zfsu1r*csu
     &  - msfer2**2*t*massorga*zftt1l*ctt
     &  - msfer2**2*t*massorga*zftt1r*ctt
     &  - msfer2**2*t*massorga*zfuu1l*cuu
     &  - msfer2**2*t*massorga*zfuu1r*cuu
     &  + 2.d0*msfer2**2*t*massorga*zftu1l*ctu
     &  + 2.d0*msfer2**2*t*massorga*zftu1r*ctu
     &  + 2.d0*msfer2**2*massorga*zftt2l*ctt
     &  + 2.d0*msfer2**2*massorga*zftt2r*ctt
     &  + msfer2**2*massorga*zfuu6l*cuu
     &  + msfer2**2*massorga*zfuu6r*cuu
     &  + msfer2**2*massorga*zfuu8l*cuu
     &  + msfer2**2*massorga*zfuu8r*cuu
     &  - 2.d0*msfer2**2*massorga*zftu3l*ctu
     &  - 2.d0*msfer2**2*massorga*zftu3r*ctu
     &  - 2.d0*msfer2**2*massorga*zftu4l*ctu
     &  - 2.d0*msfer2**2*massorga*zftu4r*ctu
     &  - 2.d0*msfer2**2*massorga*zfts*zfts4l*cts
     &  - 2.d0*msfer2**2*massorga*zfts*zfts4r*cts
     &  + 2.d0*msfer2**2*massorga*zfsu*zfsu6l*csu
     &  + 2.d0*msfer2**2*massorga*zfsu*zfsu6r*csu
     &  - 2.d0*mchi2*mgb2*mfer2*massorga*zfuu1l*cuu
     &  - 2.d0*mchi2*mgb2*mfer2*massorga*zfuu1r*cuu
     &  - mchi2*mgb2*mfer2*massorga*zfss*zfss1l*css
     &  - mchi2*mgb2*mfer2*massorga*zfss*zfss1r*css
     &  + 4.d0*mchi2*mgb2*mfer2*massorga*zfsu*zfsu1l*csu
     &  + 4.d0*mchi2*mgb2*mfer2*massorga*zfsu*zfsu1r*csu
     &  + mchi2*mgb2*s*massorga*zfuu1l*cuu
     &  + mchi2*mgb2*s*massorga*zfuu1r*cuu
     &  + mchi2*mgb2*massorga*zfuu2l*cuu
     &  + mchi2*mgb2*massorga*zfuu2r*cuu
     &  - mchi2*mgb2*massorga*zfuu3l*cuu
     &  - mchi2*mgb2*massorga*zfuu3r*cuu
     &  + mchi2*mgb2*massorga*zfuu4l*cuu
     &  + mchi2*mgb2*massorga*zfuu4r*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfuu6l*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfuu6r*cuu
     &  + 2.d0*mchi2*mgb2*massorga*zfuu7l*cuu
     &  + 2.d0*mchi2*mgb2*massorga*zfuu7r*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfuu8l*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfuu8r*cuu
     &  - mchi2*mgb2*massorga*zfss*zfss3l*css
     &  - mchi2*mgb2*massorga*zfss*zfss3r*css
     &  + mchi2*mgb2*massorga*zfss*zfss6l*css
     &  + mchi2*mgb2*massorga*zfss*zfss6r*css
     &  + mchi2*mgb2*massorga*zfss*zfss8l*css
     &  + mchi2*mgb2*massorga*zfss*zfss8r*css
     &  - 2.d0*mchi2*mgb2*massorga*zfsu*zfsu3l*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfsu*zfsu3r*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfsu*zfsu4l*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfsu*zfsu4r*csu
     &  - 4.d0*mchi2*mgb2*massorga*zfsu*zfsu6l*csu
     &  - 4.d0*mchi2*mgb2*massorga*zfsu*zfsu6r*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfsu*zfsu7l*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfsu*zfsu7r*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfsu*zfsu8l*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfsu*zfsu8r*csu
     &  + mchi2*mgb2*zftt1l*ctt
     &  + mchi2*mgb2*zftt1r*ctt
     &  + 2.d0*mchi2*mgb2*zfuu1l*cuu
     &  + 2.d0*mchi2*mgb2*zfuu1r*cuu
     &  - 2.d0*mchi2*mgb2*zfss*zfss1l*css
     &  - 2.d0*mchi2*mgb2*zfss*zfss1r*css
     &  + 4.d0*mchi2*mgb2*zfts*zfts1l*cts
     &  + 4.d0*mchi2*mgb2*zfts*zfts1r*cts
     &  - 4.d0*mchi2*mgb2*zfsu*zfsu1l*csu
     &  - 4.d0*mchi2*mgb2*zfsu*zfsu1r*csu
     &  - mchi2*mgb2**2*massorga*zfuu1l*cuu
     &  - mchi2*mgb2**2*massorga*zfuu1r*cuu
     &  - 2.d0*mchi2*mfer2*s*massorga*zfss*zfss1l*css
     &  - 2.d0*mchi2*mfer2*s*massorga*zfss*zfss1r*css
     &  - mchi2*mfer2*t*massorga*zfuu1l*cuu
     &  - mchi2*mfer2*t*massorga*zfuu1r*cuu
     &  + 2.d0*mchi2*mfer2*t*massorga*zfts*zfts1l*cts
     &  + 2.d0*mchi2*mfer2*t*massorga*zfts*zfts1r*cts
     &  - 2.d0*mchi2*mfer2*t*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*mchi2*mfer2*t*massorga*zfsu*zfsu1r*csu
     &  + 6.d0*mchi2*mfer2*zfuu1l*cuu
     &  + 6.d0*mchi2*mfer2*zfuu1r*cuu
     &  + 2.d0*mchi2*mfer2*zfss*zfss1l*css
     &  + 2.d0*mchi2*mfer2*zfss*zfss1r*css
     &  + 4.d0*mchi2*mfer2*zftu1l*ctu
     &  + 4.d0*mchi2*mfer2*zftu1r*ctu
     &  - 2.d0*mchi2*mfer2*zfts*zfts1l*cts
     &  - 2.d0*mchi2*mfer2*zfts*zfts1r*cts
     &  - 8.d0*mchi2*mfer2*zfsu*zfsu1l*csu
     &  - 8.d0*mchi2*mfer2*zfsu*zfsu1r*csu
     &  + mchi2*mfer2**2*massorga*zfss*zfss1l*css
     &  + mchi2*mfer2**2*massorga*zfss*zfss1r*css
     &  + mchi2*s*t*massorga*zfuu1l*cuu
     &  + mchi2*s*t*massorga*zfuu1r*cuu
     &  - 2.d0*mchi2*s*t*massorga*zfts*zfts1l*cts
     &  - 2.d0*mchi2*s*t*massorga*zfts*zfts1r*cts
     &  + 2.d0*mchi2*s*t*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*mchi2*s*t*massorga*zfsu*zfsu1r*csu
     &  - 4.d0*mchi2*s*zfuu1l*cuu
     &  - 4.d0*mchi2*s*zfuu1r*cuu
     &  - 2.d0*mchi2*s*zfts*zfts1l*cts
     &  - 2.d0*mchi2*s*zfts*zfts1r*cts
     &  + mchi2*s**2*massorga*zfss*zfss1l*css
     &  + mchi2*s**2*massorga*zfss*zfss1r*css
     &  - 2.d0*mchi2*t*zftt1l*ctt
     &  - 2.d0*mchi2*t*zftt1r*ctt
     &  - 2.d0*mchi2*t*zfuu1l*cuu
     &  - 2.d0*mchi2*t*zfuu1r*cuu
     &  + 2.d0*mchi2*t*zftu1l*ctu
     &  + 2.d0*mchi2*t*zftu1r*ctu
     &  - 4.d0*mchi2*t*zfts*zfts1l*cts
     &  - 4.d0*mchi2*t*zfts*zfts1r*cts
     &  + 4.d0*mchi2*t*zfsu*zfsu1l*csu
     &  + 4.d0*mchi2*t*zfsu*zfsu1r*csu
     &  + mchi2*t**2*massorga*zftt1l*ctt
     &  + mchi2*t**2*massorga*zftt1r*ctt
     &  + mchi2*t**2*massorga*zfuu1l*cuu
     &  + mchi2*t**2*massorga*zfuu1r*cuu
     &  - 2.d0*mchi2*t**2*massorga*zftu1l*ctu
     &  - 2.d0*mchi2*t**2*massorga*zftu1r*ctu
     &  - 4.d0*mchi2*zfuu2l*cuu
     &  - 4.d0*mchi2*zfuu2r*cuu
     &  + 2.d0*mchi2*zfuu3l*cuu
     &  + 2.d0*mchi2*zfuu3r*cuu
     &  - 4.d0*mchi2*zfuu4l*cuu
     &  - 4.d0*mchi2*zfuu4r*cuu
     &  + 4.d0*mchi2*zfuu6l*cuu
     &  + 4.d0*mchi2*zfuu6r*cuu
     &  - 8.d0*mchi2*zfuu7l*cuu
     &  - 8.d0*mchi2*zfuu7r*cuu
     &  + 4.d0*mchi2*zfuu8l*cuu
     &  + 4.d0*mchi2*zfuu8r*cuu
     &  + 2.d0*mchi2*zfss*zfss3l*css
     &  + 2.d0*mchi2*zfss*zfss3r*css
     &  - 4.d0*mchi2*zfss*zfss6l*css
     &  - 4.d0*mchi2*zfss*zfss6r*css
     &  - 4.d0*mchi2*zfss*zfss8l*css
     &  - 4.d0*mchi2*zfss*zfss8r*css
     &  + 4.d0*mchi2*zftu2l*ctu
     &  + 4.d0*mchi2*zftu2r*ctu
     &  + 4.d0*mchi2*zftu3l*ctu
     &  + 4.d0*mchi2*zftu3r*ctu
     &  + 4.d0*mchi2*zftu4l*ctu
     &  + 4.d0*mchi2*zftu4r*ctu
     &  + 4.d0*mchi2*zfts*zfts2l*cts
     &  + 4.d0*mchi2*zfts*zfts2r*cts
     &  + 4.d0*mchi2*zfts*zfts3l*cts
     &  + 4.d0*mchi2*zfts*zfts3r*cts
     &  + 4.d0*mchi2*zfts*zfts4l*cts
     &  + 4.d0*mchi2*zfts*zfts4r*cts
     &  + 4.d0*mchi2*zfsu*zfsu3l*csu
     &  + 4.d0*mchi2*zfsu*zfsu3r*csu
     &  - 8.d0*mchi2*zfsu*zfsu4l*csu
     &  - 8.d0*mchi2*zfsu*zfsu4r*csu
     &  + 8.d0*mchi2*zfsu*zfsu6l*csu
     &  + 8.d0*mchi2*zfsu*zfsu6r*csu
     &  + 4.d0*mchi2*zfsu*zfsu7l*csu
     &  + 4.d0*mchi2*zfsu*zfsu7r*csu
     &  - 8.d0*mchi2*zfsu*zfsu8l*csu
     &  - 8.d0*mchi2*zfsu*zfsu8r*csu
     &  - mchi2**2*mgb2*massorga*zfuu1l*cuu
     &  - mchi2**2*mgb2*massorga*zfuu1r*cuu
     &  + 2.d0*mchi2**2*zfuu1l*cuu
     &  + 2.d0*mchi2**2*zfuu1r*cuu
     &  + 4.d0*mchi2**2*zfts*zfts1l*cts
     &  + 4.d0*mchi2**2*zfts*zfts1r*cts
     &  - 4.d0*mchi2**2*zfsu*zfsu1l*csu
     &  - 4.d0*mchi2**2*zfsu*zfsu1r*csu
     &  - 2.d0*mgb2*mfer2*s*massorga*zfuu1l*cuu
     &  - 2.d0*mgb2*mfer2*s*massorga*zfuu1r*cuu
     &  - mgb2*mfer2*s*massorga*zfss*zfss1l*css
     &  - mgb2*mfer2*s*massorga*zfss*zfss1r*css
     &  - 4.d0*mgb2*mfer2*t*massorga*zfuu1l*cuu
     &  - 4.d0*mgb2*mfer2*t*massorga*zfuu1r*cuu
     &  + 2.d0*mgb2*mfer2*t*massorga*zftu1l*ctu
     &  + 2.d0*mgb2*mfer2*t*massorga*zftu1r*ctu
     &  + 2.d0*mgb2*mfer2*massorga*zfuu2l*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfuu2r*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfuu4l*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfuu4r*cuu
     &  + mgb2*mfer2*massorga*zfuu6l*cuu
     &  + mgb2*mfer2*massorga*zfuu6r*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfuu7l*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfuu7r*cuu
     &  + mgb2*mfer2*massorga*zfuu8l*cuu
     &  + mgb2*mfer2*massorga*zfuu8r*cuu
     &  - mgb2*mfer2*massorga*zfss*zfss2l*css
     &  - mgb2*mfer2*massorga*zfss*zfss2r*css
     &  - mgb2*mfer2*massorga*zfss*zfss4l*css
     &  - mgb2*mfer2*massorga*zfss*zfss4r*css
     &  - 2.d0*mgb2*mfer2*massorga*zfsu*zfsu2l*csu
     &  - 2.d0*mgb2*mfer2*massorga*zfsu*zfsu2r*csu
     &  + 4.d0*mgb2*mfer2*massorga*zfsu*zfsu4l*csu
     &  + 4.d0*mgb2*mfer2*massorga*zfsu*zfsu4r*csu
     &  + 2.d0*mgb2*mfer2*massorga*zfsu*zfsu6l*csu
     &  + 2.d0*mgb2*mfer2*massorga*zfsu*zfsu6r*csu
     &  + mgb2*mfer2*zftt1l*ctt
     &  + mgb2*mfer2*zftt1r*ctt
     &  - 2.d0*mgb2*mfer2*zftu1l*ctu
     &  - 2.d0*mgb2*mfer2*zftu1r*ctu
     &  + 2.d0*mgb2*mfer2**2*massorga*zfuu1l*cuu
     &  + 2.d0*mgb2*mfer2**2*massorga*zfuu1r*cuu
     &  + 2.d0*mgb2*s*t*massorga*zfuu1l*cuu
     &  + 2.d0*mgb2*s*t*massorga*zfuu1r*cuu
     &  + 2.d0*mgb2*s*t*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*mgb2*s*t*massorga*zfsu*zfsu1r*csu
     &  - mgb2*s*massorga*zfuu2l*cuu
     &  - mgb2*s*massorga*zfuu2r*cuu
     &  - mgb2*s*massorga*zfuu3l*cuu
     &  - mgb2*s*massorga*zfuu3r*cuu
     &  - mgb2*s*massorga*zfuu4l*cuu
     &  - mgb2*s*massorga*zfuu4r*cuu
     &  - mgb2*s*massorga*zfuu6l*cuu
     &  - mgb2*s*massorga*zfuu6r*cuu
     &  - 2.d0*mgb2*s*massorga*zfuu7l*cuu
     &  - 2.d0*mgb2*s*massorga*zfuu7r*cuu
     &  - mgb2*s*massorga*zfuu8l*cuu
     &  - mgb2*s*massorga*zfuu8r*cuu
     &  - mgb2*s*massorga*zfss*zfss2l*css
     &  - mgb2*s*massorga*zfss*zfss2r*css
     &  - mgb2*s*massorga*zfss*zfss3l*css
     &  - mgb2*s*massorga*zfss*zfss3r*css
     &  - mgb2*s*massorga*zfss*zfss4l*css
     &  - mgb2*s*massorga*zfss*zfss4r*css
     &  + mgb2*s*massorga*zfss*zfss6l*css
     &  + mgb2*s*massorga*zfss*zfss6r*css
     &  + 2.d0*mgb2*s*massorga*zfss*zfss7l*css
     &  + 2.d0*mgb2*s*massorga*zfss*zfss7r*css
     &  + mgb2*s*massorga*zfss*zfss8l*css
     &  + mgb2*s*massorga*zfss*zfss8r*css
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu2l*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu2r*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu3l*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu3r*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu4l*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu4r*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu6l*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu6r*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu7l*csu
     &  - 2.d0*mgb2*s*massorga*zfsu*zfsu7r*csu
     &  + 2.d0*mgb2*s*massorga*zfsu*zfsu8l*csu
     &  + 2.d0*mgb2*s*massorga*zfsu*zfsu8r*csu
     &  - 2.d0*mgb2*s*zfuu1l*cuu
     &  - 2.d0*mgb2*s*zfuu1r*cuu
     &  - 2.d0*mgb2*s*zfss*zfss1l*css
     &  - 2.d0*mgb2*s*zfss*zfss1r*css
     &  - 4.d0*mgb2*s*zfsu*zfsu1l*csu
     &  - 4.d0*mgb2*s*zfsu*zfsu1r*csu
     &  - mgb2*t*massorga*zfuu2l*cuu
     &  - mgb2*t*massorga*zfuu2r*cuu
     &  - mgb2*t*massorga*zfuu4l*cuu
     &  - mgb2*t*massorga*zfuu4r*cuu
     &  - mgb2*t*massorga*zfuu6l*cuu
     &  - mgb2*t*massorga*zfuu6r*cuu
     &  - 2.d0*mgb2*t*massorga*zfuu7l*cuu
     &  - 2.d0*mgb2*t*massorga*zfuu7r*cuu
     &  - mgb2*t*massorga*zfuu8l*cuu
     &  - mgb2*t*massorga*zfuu8r*cuu
     &  + 2.d0*mgb2*t*massorga*zftu2l*ctu
     &  + 2.d0*mgb2*t*massorga*zftu2r*ctu
     &  + 4.d0*mgb2*t*massorga*zftu4l*ctu
     &  + 4.d0*mgb2*t*massorga*zftu4r*ctu
     &  + 2.d0*mgb2*t*massorga*zfts*zfts2l*cts
     &  + 2.d0*mgb2*t*massorga*zfts*zfts2r*cts
      tmp2=
     &  - 2.d0*mgb2*t*massorga*zfts*zfts3l*cts
     &  - 2.d0*mgb2*t*massorga*zfts*zfts3r*cts
     &  - 2.d0*mgb2*t*massorga*zfsu*zfsu4l*csu
     &  - 2.d0*mgb2*t*massorga*zfsu*zfsu4r*csu
     &  - 2.d0*mgb2*t*massorga*zfsu*zfsu6l*csu
     &  - 2.d0*mgb2*t*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*mgb2*t*massorga*zfsu*zfsu7l*csu
     &  + 2.d0*mgb2*t*massorga*zfsu*zfsu7r*csu
     &  - mgb2*t*zftt1l*ctt
     &  - mgb2*t*zftt1r*ctt
     &  + 2.d0*mgb2*t*zftu1l*ctu
     &  + 2.d0*mgb2*t*zftu1r*ctu
     &  + 2.d0*mgb2*t**2*massorga*zfuu1l*cuu
     &  + 2.d0*mgb2*t**2*massorga*zfuu1r*cuu
     &  - 2.d0*mgb2*t**2*massorga*zftu1l*ctu
     &  - 2.d0*mgb2*t**2*massorga*zftu1r*ctu
     &  + 2.d0*mgb2*massorga*zfuu5l*cuu
     &  + 2.d0*mgb2*massorga*zfuu5r*cuu
     &  + 2.d0*mgb2*massorga*zfss*zfss5l*css
     &  + 2.d0*mgb2*massorga*zfss*zfss5r*css
     &  + 4.d0*mgb2*massorga*zfsu*zfsu5l*csu
     &  + 4.d0*mgb2*massorga*zfsu*zfsu5r*csu
     &  + 2.d0*mgb2*zftt2l*ctt
     &  + 2.d0*mgb2*zftt2r*ctt
     &  - 4.d0*mgb2*zfuu2l*cuu
     &  - 4.d0*mgb2*zfuu2r*cuu
     &  - 4.d0*mgb2*zfuu4l*cuu
     &  - 4.d0*mgb2*zfuu4r*cuu
     &  - 8.d0*mgb2*zfuu7l*cuu
     &  - 8.d0*mgb2*zfuu7r*cuu
     &  - 2.d0*mgb2*zfss*zfss2l*css
     &  - 2.d0*mgb2*zfss*zfss2r*css
     &  - 2.d0*mgb2*zfss*zfss4l*css
     &  - 2.d0*mgb2*zfss*zfss4r*css
     &  + 2.d0*mgb2*zftu2l*ctu
     &  + 2.d0*mgb2*zftu2r*ctu
     &  + 2.d0*mgb2*zfts*zfts2l*cts
     &  + 2.d0*mgb2*zfts*zfts2r*cts
     &  + 2.d0*mgb2*zfts*zfts3l*cts
     &  + 2.d0*mgb2*zfts*zfts3r*cts
     &  - 4.d0*mgb2*zfsu*zfsu2l*csu
     &  - 4.d0*mgb2*zfsu*zfsu2r*csu
     &  - 8.d0*mgb2*zfsu*zfsu4l*csu
     &  - 8.d0*mgb2*zfsu*zfsu4r*csu
     &  - 4.d0*mgb2*zfsu*zfsu7l*csu
     &  - 4.d0*mgb2*zfsu*zfsu7r*csu
     &  + mgb2**2*mfer2*massorga*zfuu1l*cuu
     &  + mgb2**2*mfer2*massorga*zfuu1r*cuu
     &  - mgb2**2*t*massorga*zfuu1l*cuu
     &  - mgb2**2*t*massorga*zfuu1r*cuu
     &  + mgb2**2*massorga*zfuu2l*cuu
     &  + mgb2**2*massorga*zfuu2r*cuu
     &  + mgb2**2*massorga*zfuu4l*cuu
     &  + mgb2**2*massorga*zfuu4r*cuu
     &  + 2.d0*mgb2**2*massorga*zfuu7l*cuu
     &  + 2.d0*mgb2**2*massorga*zfuu7r*cuu
     &  + 2.d0*mgb2**2*massorga*zfsu*zfsu4l*csu
     &  + 2.d0*mgb2**2*massorga*zfsu*zfsu4r*csu
     &  - 2.d0*mgb2**2*massorga*zfsu*zfsu7l*csu
     &  - 2.d0*mgb2**2*massorga*zfsu*zfsu7r*csu
     &  + 4.d0*mfer2*s*t*massorga*zfuu1l*cuu
     &  + 4.d0*mfer2*s*t*massorga*zfuu1r*cuu
     &  + mfer2*s*t*massorga*zfss*zfss1l*css
     &  + mfer2*s*t*massorga*zfss*zfss1r*css
     &  - 2.d0*mfer2*s*t*massorga*zftu1l*ctu
     &  - 2.d0*mfer2*s*t*massorga*zftu1r*ctu
     &  + 2.d0*mfer2*s*t*massorga*zfsu*zfsu1l*csu
     &  + 2.d0*mfer2*s*t*massorga*zfsu*zfsu1r*csu
     &  - 2.d0*mfer2*s*massorga*zfuu3l*cuu
     &  - 2.d0*mfer2*s*massorga*zfuu3r*cuu
     &  - 2.d0*mfer2*s*massorga*zfuu6l*cuu
     &  - 2.d0*mfer2*s*massorga*zfuu6r*cuu
     &  - 2.d0*mfer2*s*massorga*zfuu8l*cuu
     &  - 2.d0*mfer2*s*massorga*zfuu8r*cuu
     &  - 2.d0*mfer2*s*massorga*zfss*zfss2l*css
     &  - 2.d0*mfer2*s*massorga*zfss*zfss2r*css
     &  - 2.d0*mfer2*s*massorga*zfss*zfss3l*css
     &  - 2.d0*mfer2*s*massorga*zfss*zfss3r*css
     &  - 2.d0*mfer2*s*massorga*zfss*zfss4l*css
     &  - 2.d0*mfer2*s*massorga*zfss*zfss4r*css
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu2l*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu2r*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu3l*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu3r*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu6l*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu6r*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu7l*csu
     &  - 4.d0*mfer2*s*massorga*zfsu*zfsu7r*csu
     &  - 4.d0*mfer2*s*zfuu1l*cuu
     &  - 4.d0*mfer2*s*zfuu1r*cuu
     &  - 2.d0*mfer2*s*zftu1l*ctu
     &  - 2.d0*mfer2*s*zftu1r*ctu
     &  - 8.d0*mfer2*s*zfsu*zfsu1l*csu
     &  - 8.d0*mfer2*s*zfsu*zfsu1r*csu
     &  + mfer2*s**2*massorga*zfuu1l*cuu
     &  + mfer2*s**2*massorga*zfuu1r*cuu
     &  - mfer2*t*massorga*zfuu3l*cuu
     &  - mfer2*t*massorga*zfuu3r*cuu
     &  - 2.d0*mfer2*t*massorga*zfuu6l*cuu
     &  - 2.d0*mfer2*t*massorga*zfuu6r*cuu
     &  - 2.d0*mfer2*t*massorga*zfuu8l*cuu
     &  - 2.d0*mfer2*t*massorga*zfuu8r*cuu
     &  - mfer2*t*massorga*zfss*zfss3l*css
     &  - mfer2*t*massorga*zfss*zfss3r*css
     &  + 2.d0*mfer2*t*massorga*zftu2l*ctu
     &  + 2.d0*mfer2*t*massorga*zftu2r*ctu
     &  + 2.d0*mfer2*t*massorga*zftu3l*ctu
     &  + 2.d0*mfer2*t*massorga*zftu3r*ctu
     &  + 2.d0*mfer2*t*massorga*zftu4l*ctu
     &  + 2.d0*mfer2*t*massorga*zftu4r*ctu
     &  + 2.d0*mfer2*t*massorga*zfts*zfts2l*cts
     &  + 2.d0*mfer2*t*massorga*zfts*zfts2r*cts
     &  + 2.d0*mfer2*t*massorga*zfts*zfts3l*cts
     &  + 2.d0*mfer2*t*massorga*zfts*zfts3r*cts
     &  + 2.d0*mfer2*t*massorga*zfts*zfts4l*cts
     &  + 2.d0*mfer2*t*massorga*zfts*zfts4r*cts
     &  - 2.d0*mfer2*t*massorga*zfsu*zfsu3l*csu
     &  - 2.d0*mfer2*t*massorga*zfsu*zfsu3r*csu
     &  - 4.d0*mfer2*t*massorga*zfsu*zfsu6l*csu
     &  - 4.d0*mfer2*t*massorga*zfsu*zfsu6r*csu
     &  - 2.d0*mfer2*t*massorga*zfsu*zfsu7l*csu
     &  - 2.d0*mfer2*t*massorga*zfsu*zfsu7r*csu
     &  - 2.d0*mfer2*t*zftt1l*ctt
     &  - 2.d0*mfer2*t*zftt1r*ctt
     &  - 2.d0*mfer2*t*zfuu1l*cuu
     &  - 2.d0*mfer2*t*zfuu1r*cuu
     &  + 4.d0*mfer2*t*zftu1l*ctu
     &  + 4.d0*mfer2*t*zftu1r*ctu
     &  + mfer2*t**2*massorga*zftt1l*ctt
     &  + mfer2*t**2*massorga*zftt1r*ctt
     &  + 3.d0*mfer2*t**2*massorga*zfuu1l*cuu
     &  + 3.d0*mfer2*t**2*massorga*zfuu1r*cuu
     &  - 4.d0*mfer2*t**2*massorga*zftu1l*ctu
     &  - 4.d0*mfer2*t**2*massorga*zftu1r*ctu
     &  - 8.d0*mfer2*zfuu2l*cuu
     &  - 8.d0*mfer2*zfuu2r*cuu
     &  + 2.d0*mfer2*zfuu3l*cuu
     &  + 2.d0*mfer2*zfuu3r*cuu
     &  - 8.d0*mfer2*zfuu4l*cuu
     &  - 8.d0*mfer2*zfuu4r*cuu
     &  + 2.d0*mfer2*zfuu6l*cuu
     &  + 2.d0*mfer2*zfuu6r*cuu
     &  - 8.d0*mfer2*zfuu7l*cuu
     &  - 8.d0*mfer2*zfuu7r*cuu
     &  + 2.d0*mfer2*zfuu8l*cuu
     &  + 2.d0*mfer2*zfuu8r*cuu
     &  + 2.d0*mfer2*zfss*zfss2l*css
     &  + 2.d0*mfer2*zfss*zfss2r*css
     &  + 2.d0*mfer2*zfss*zfss3l*css
     &  + 2.d0*mfer2*zfss*zfss3r*css
     &  + 2.d0*mfer2*zfss*zfss4l*css
     &  + 2.d0*mfer2*zfss*zfss4r*css
     &  - 2.d0*mfer2*zftu2l*ctu
     &  - 2.d0*mfer2*zftu2r*ctu
     &  - 2.d0*mfer2*zftu3l*ctu
     &  - 2.d0*mfer2*zftu3r*ctu
     &  - 2.d0*mfer2*zftu4l*ctu
     &  - 2.d0*mfer2*zftu4r*ctu
     &  - 2.d0*mfer2*zfts*zfts2l*cts
     &  - 2.d0*mfer2*zfts*zfts2r*cts
     &  - 2.d0*mfer2*zfts*zfts3l*cts
     &  - 2.d0*mfer2*zfts*zfts3r*cts
     &  - 2.d0*mfer2*zfts*zfts4l*cts
     &  - 2.d0*mfer2*zfts*zfts4r*cts
     &  + 4.d0*mfer2*zfsu*zfsu2l*csu
     &  + 4.d0*mfer2*zfsu*zfsu2r*csu
     &  + 4.d0*mfer2*zfsu*zfsu3l*csu
     &  + 4.d0*mfer2*zfsu*zfsu3r*csu
     &  - 16.d0*mfer2*zfsu*zfsu4l*csu
     &  - 16.d0*mfer2*zfsu*zfsu4r*csu
     &  + 4.d0*mfer2*zfsu*zfsu6l*csu
     &  + 4.d0*mfer2*zfsu*zfsu6r*csu
     &  + 4.d0*mfer2*zfsu*zfsu7l*csu
     &  + 4.d0*mfer2*zfsu*zfsu7r*csu
     &  - 2.d0*mfer2**2*s*massorga*zfuu1l*cuu
     &  - 2.d0*mfer2**2*s*massorga*zfuu1r*cuu
     &  - 3.d0*mfer2**2*t*massorga*zfuu1l*cuu
     &  - 3.d0*mfer2**2*t*massorga*zfuu1r*cuu
     &  + 2.d0*mfer2**2*t*massorga*zftu1l*ctu
     &  + 2.d0*mfer2**2*t*massorga*zftu1r*ctu
     &  + mfer2**2*massorga*zfuu3l*cuu
     &  + mfer2**2*massorga*zfuu3r*cuu
     &  + mfer2**2*massorga*zfuu6l*cuu
     &  + mfer2**2*massorga*zfuu6r*cuu
     &  + mfer2**2*massorga*zfuu8l*cuu
     &  + mfer2**2*massorga*zfuu8r*cuu
     &  + mfer2**2*massorga*zfss*zfss2l*css
     &  + mfer2**2*massorga*zfss*zfss2r*css
     &  + mfer2**2*massorga*zfss*zfss3l*css
     &  + mfer2**2*massorga*zfss*zfss3r*css
     &  + mfer2**2*massorga*zfss*zfss4l*css
     &  + mfer2**2*massorga*zfss*zfss4r*css
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu2l*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu2r*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu3l*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu3r*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu6l*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu7l*csu
     &  + 2.d0*mfer2**2*massorga*zfsu*zfsu7r*csu
     &  + 2.d0*mfer2**2*zfuu1l*cuu
     &  + 2.d0*mfer2**2*zfuu1r*cuu
     &  - 2.d0*mfer2**2*zftu1l*ctu
     &  - 2.d0*mfer2**2*zftu1r*ctu
     &  + mfer2**3*massorga*zfuu1l*cuu
     &  + mfer2**3*massorga*zfuu1r*cuu
     &  + s*t*massorga*zfuu3l*cuu
     &  + s*t*massorga*zfuu3r*cuu
     &  + 2.d0*s*t*massorga*zfuu6l*cuu
     &  + 2.d0*s*t*massorga*zfuu6r*cuu
     &  + 2.d0*s*t*massorga*zfuu8l*cuu
     &  + 2.d0*s*t*massorga*zfuu8r*cuu
     &  + s*t*massorga*zfss*zfss3l*css
     &  + s*t*massorga*zfss*zfss3r*css
     &  - 2.d0*s*t*massorga*zftu2l*ctu
     &  - 2.d0*s*t*massorga*zftu2r*ctu
     &  - 2.d0*s*t*massorga*zftu3l*ctu
     &  - 2.d0*s*t*massorga*zftu3r*ctu
     &  - 2.d0*s*t*massorga*zftu4l*ctu
     &  - 2.d0*s*t*massorga*zftu4r*ctu
     &  - 2.d0*s*t*massorga*zfts*zfts2l*cts
     &  - 2.d0*s*t*massorga*zfts*zfts2r*cts
     &  - 2.d0*s*t*massorga*zfts*zfts3l*cts
     &  - 2.d0*s*t*massorga*zfts*zfts3r*cts
     &  - 2.d0*s*t*massorga*zfts*zfts4l*cts
     &  - 2.d0*s*t*massorga*zfts*zfts4r*cts
     &  + 2.d0*s*t*massorga*zfsu*zfsu3l*csu
     &  + 2.d0*s*t*massorga*zfsu*zfsu3r*csu
     &  + 4.d0*s*t*massorga*zfsu*zfsu6l*csu
     &  + 4.d0*s*t*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*s*t*massorga*zfsu*zfsu7l*csu
     &  + 2.d0*s*t*massorga*zfsu*zfsu7r*csu
     &  + 2.d0*s*t*zfuu1l*cuu
     &  + 2.d0*s*t*zfuu1r*cuu
     &  + 2.d0*s*t*zfss*zfss1l*css
     &  + 2.d0*s*t*zfss*zfss1r*css
     &  - 2.d0*s*t*zftu1l*ctu
     &  - 2.d0*s*t*zftu1r*ctu
     &  - 2.d0*s*t*zfts*zfts1l*cts
     &  - 2.d0*s*t*zfts*zfts1r*cts
     &  + 4.d0*s*t*zfsu*zfsu1l*csu
     &  + 4.d0*s*t*zfsu*zfsu1r*csu
     &  - 2.d0*s*t**2*massorga*zfuu1l*cuu
     &  - 2.d0*s*t**2*massorga*zfuu1r*cuu
     &  + 2.d0*s*t**2*massorga*zftu1l*ctu
     &  + 2.d0*s*t**2*massorga*zftu1r*ctu
     &  + 2.d0*s*t**2*massorga*zfts*zfts1l*cts
     &  + 2.d0*s*t**2*massorga*zfts*zfts1r*cts
     &  - 2.d0*s*t**2*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*s*t**2*massorga*zfsu*zfsu1r*csu
     &  + 4.d0*s*zfuu2l*cuu
     &  + 4.d0*s*zfuu2r*cuu
     &  + 4.d0*s*zfuu4l*cuu
     &  + 4.d0*s*zfuu4r*cuu
     &  - 2.d0*s*zfuu6l*cuu
     &  - 2.d0*s*zfuu6r*cuu
     &  + 8.d0*s*zfuu7l*cuu
     &  + 8.d0*s*zfuu7r*cuu
     &  - 2.d0*s*zfuu8l*cuu
     &  - 2.d0*s*zfuu8r*cuu
     &  + 2.d0*s*zfss*zfss2l*css
     &  + 2.d0*s*zfss*zfss2r*css
     &  + 2.d0*s*zfss*zfss4l*css
     &  + 2.d0*s*zfss*zfss4r*css
     &  - 4.d0*s*zfss*zfss6l*css
     &  - 4.d0*s*zfss*zfss6r*css
     &  - 8.d0*s*zfss*zfss7l*css
     &  - 8.d0*s*zfss*zfss7r*css
     &  - 4.d0*s*zfss*zfss8l*css
     &  - 4.d0*s*zfss*zfss8r*css
     &  - 2.d0*s*zftu2l*ctu
     &  - 2.d0*s*zftu2r*ctu
     &  - 2.d0*s*zftu3l*ctu
     &  - 2.d0*s*zftu3r*ctu
     &  - 2.d0*s*zftu4l*ctu
     &  - 2.d0*s*zftu4r*ctu
     &  - 2.d0*s*zfts*zfts2l*cts
     &  - 2.d0*s*zfts*zfts2r*cts
     &  - 2.d0*s*zfts*zfts3l*cts
     &  - 2.d0*s*zfts*zfts3r*cts
     &  - 2.d0*s*zfts*zfts4l*cts
     &  - 2.d0*s*zfts*zfts4r*cts
     &  + 4.d0*s*zfsu*zfsu2l*csu
     &  + 4.d0*s*zfsu*zfsu2r*csu
     &  + 8.d0*s*zfsu*zfsu4l*csu
     &  + 8.d0*s*zfsu*zfsu4r*csu
     &  - 4.d0*s*zfsu*zfsu6l*csu
     &  - 4.d0*s*zfsu*zfsu6r*csu
     &  - 8.d0*s*zfsu*zfsu8l*csu
     &  - 8.d0*s*zfsu*zfsu8r*csu
     &  - s**2*t*massorga*zfuu1l*cuu
     &  - s**2*t*massorga*zfuu1r*cuu
     &  - s**2*t*massorga*zfss*zfss1l*css
     &  - s**2*t*massorga*zfss*zfss1r*css
     &  - 2.d0*s**2*t*massorga*zfsu*zfsu1l*csu
     &  - 2.d0*s**2*t*massorga*zfsu*zfsu1r*csu
     &  + s**2*massorga*zfuu3l*cuu
     &  + s**2*massorga*zfuu3r*cuu
     &  + s**2*massorga*zfuu6l*cuu
     &  + s**2*massorga*zfuu6r*cuu
     &  + s**2*massorga*zfuu8l*cuu
     &  + s**2*massorga*zfuu8r*cuu
     &  + s**2*massorga*zfss*zfss2l*css
     &  + s**2*massorga*zfss*zfss2r*css
     &  + s**2*massorga*zfss*zfss3l*css
     &  + s**2*massorga*zfss*zfss3r*css
     &  + s**2*massorga*zfss*zfss4l*css
     &  + s**2*massorga*zfss*zfss4r*css
     &  + 2.d0*s**2*massorga*zfsu*zfsu2l*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu2r*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu3l*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu3r*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu6l*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu7l*csu
     &  + 2.d0*s**2*massorga*zfsu*zfsu7r*csu
     &  + 2.d0*s**2*zfuu1l*cuu
     &  + 2.d0*s**2*zfuu1r*cuu
     &  + 2.d0*s**2*zfss*zfss1l*css
     &  + 2.d0*s**2*zfss*zfss1r*css
     &  + 4.d0*s**2*zfsu*zfsu1l*csu
     &  + 4.d0*s**2*zfsu*zfsu1r*csu
     &  - 4.d0*t*zftt2l*ctt
     &  - 4.d0*t*zftt2r*ctt
     &  + 4.d0*t*zfuu2l*cuu
     &  + 4.d0*t*zfuu2r*cuu
     &  - 2.d0*t*zfuu3l*cuu
     &  - 2.d0*t*zfuu3r*cuu
     &  + 4.d0*t*zfuu4l*cuu
     &  + 4.d0*t*zfuu4r*cuu
     &  - 2.d0*t*zfuu6l*cuu
     &  - 2.d0*t*zfuu6r*cuu
     &  + 8.d0*t*zfuu7l*cuu
     &  + 8.d0*t*zfuu7r*cuu
     &  - 2.d0*t*zfuu8l*cuu
     &  - 2.d0*t*zfuu8r*cuu
     &  - 2.d0*t*zfss*zfss3l*css
     &  - 2.d0*t*zfss*zfss3r*css
     &  - 4.d0*t*zftu2l*ctu
     &  - 4.d0*t*zftu2r*ctu
     &  + 2.d0*t*zftu3l*ctu
     &  + 2.d0*t*zftu3r*ctu
     &  - 2.d0*t*zftu4l*ctu
     &  - 2.d0*t*zftu4r*ctu
     &  - 4.d0*t*zfts*zfts2l*cts
     &  - 4.d0*t*zfts*zfts2r*cts
     &  + 2.d0*t*zfts*zfts4l*cts
     &  + 2.d0*t*zfts*zfts4r*cts
     &  - 4.d0*t*zfsu*zfsu3l*csu
     &  - 4.d0*t*zfsu*zfsu3r*csu
     &  + 8.d0*t*zfsu*zfsu4l*csu
     &  + 8.d0*t*zfsu*zfsu4r*csu
     &  - 4.d0*t*zfsu*zfsu6l*csu
     &  - 4.d0*t*zfsu*zfsu6r*csu
     &  + 2.d0*t**2*massorga*zftt2l*ctt
     &  + 2.d0*t**2*massorga*zftt2r*ctt
     &  + t**2*massorga*zfuu6l*cuu
     &  + t**2*massorga*zfuu6r*cuu
     &  + t**2*massorga*zfuu8l*cuu
     &  + t**2*massorga*zfuu8r*cuu
     &  - 2.d0*t**2*massorga*zftu3l*ctu
     &  - 2.d0*t**2*massorga*zftu3r*ctu
     &  - 2.d0*t**2*massorga*zftu4l*ctu
     &  - 2.d0*t**2*massorga*zftu4r*ctu
     &  - 2.d0*t**2*massorga*zfts*zfts4l*cts
     &  - 2.d0*t**2*massorga*zfts*zfts4r*cts
     &  + 2.d0*t**2*massorga*zfsu*zfsu6l*csu
     &  + 2.d0*t**2*massorga*zfsu*zfsu6r*csu
     &  + 2.d0*t**2*zftt1l*ctt
     &  + 2.d0*t**2*zftt1r*ctt
     &  - 2.d0*t**2*zftu1l*ctu
     &  - 2.d0*t**2*zftu1r*ctu
     &  - t**3*massorga*zftt1l*ctt
     &  - t**3*massorga*zftt1r*ctt
     &  - t**3*massorga*zfuu1l*cuu
     &  - t**3*massorga*zfuu1r*cuu
     &  + 2.d0*t**3*massorga*zftu1l*ctu
     &  + 2.d0*t**3*massorga*zftu1r*ctu
     &  - 8.d0*zfuu5l*cuu
     &  - 8.d0*zfuu5r*cuu
     &  - 8.d0*zfss*zfss5l*css
     &  - 8.d0*zfss*zfss5r*css
     &  - 16.d0*zfsu*zfsu5l*csu
     &  - 16.d0*zfsu*zfsu5r*csu

      dsaschicaseaas=tmp1+tmp2
      if (dble(dsaschicaseaas).lt.0.0d0) then     
        if(aszeroprint) then
        write(*,*) ' '
        write(*,*) 
     &    'DS: ERROR IN dsaschicasea with negative cross section:'
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
        write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
        write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
        write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
        write(*,*) 'DS: dsaschicaseaas = ',dsaschicaseaas
        write(*,*) 'DS: mass1=',mass1
        write(*,*) 'DS: mass2=',mass2
        write(*,*) 'DS: mass3=',mass3
        write(*,*) 'DS: mass4=',mass4
        write(*,*) 'DS: s=',s
        write(*,*) 'DS: t=',t
        write(*,*) 'DS: u=',u
        write(*,*) 'DS: p12=',p12
        write(*,*) 'DS: costheta=',costheta
        endif
        dsaschicaseaas=dcmplx(0.0d0,0.0d0)
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsaschicaseaas)*k34/(8.0d0*pi*gg1c*gg2c*s34*dsqrt(s))
      return
      end
