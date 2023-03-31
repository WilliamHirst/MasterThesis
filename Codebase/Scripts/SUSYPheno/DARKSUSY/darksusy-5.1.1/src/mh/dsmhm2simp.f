***********************************************************************
*** dsmhm2simp returns the scattering amplitude at zero 
*** momentum transfer squared, averaged over the initial spin states of 
*** the DM particle and summed over all other spin states, in the limit 
*** of relativistic scattering partners with small energies omega. In 
*** this limit, it can be expanded as
***
***    |M|**2 = cn*(omega/m0)**n + O( (omega/m0)**(n+1) )
*** 
***  input: m0       - DM mass (in GeV)
***         SMtype   - SM scattering partners:
***                    7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      subroutine dsmhm2simp(m0,SMtype,cn,n)
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'
      include 'dsidtag.h'

      integer n, SMtype,kf,nsf,ksf(2)
      real*8  m0,cn
      integer i,j

      real*8 tmpres,sv,dssigmav

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/


      cn=0d0
      n=0
      tmpres=0d0



c... neutralino DM     
      if (mhtype.eq.1.or.mhtype.eq.2) then
        n=2
        if (mhtype.eq.1) sv=dssigmav(0) ! make sure all vertices 
                                        ! and masses are calculated...

c... set up sparticles/KK-particles
        if (SMtype.le.3) then
          nsf=1
          kf=knu(SMtype)
        endif
        if (SMtype.ge.4) then
          nsf=2
          kf=kl(SMtype-3)
        endif
        if (SMtype.ge.7) return
        do 10 i=1,nsf
          if (SMtype.le.3) ksf(i)=ksnu(SMtype)
          if (SMtype.ge.4) ksf(i)=ksl(SMtype+3*i-6)
  10    continue

c... neutralino DM
        if (mhtype.eq.1) then
          tmpres=8*(m0/mass(kz))**4
     -            *abs(gl(kz,kn(1),kn(1)))**2
     -            *(abs(gl(kz,kf,kf))**2
     -              +abs(gr(kz,kf,kf))**2)
          do 50 i=1,nsf
           do 40 j=1,nsf
            tmpres=tmpres+2*
     -     (abs(gl(ksf(i),kf,kn(1)))**2
     -      *abs(gl(ksf(j),kf,kn(1)))**2
     -      +abs(gr(ksf(i),kf,kn(1)))**2
     -      *abs(gr(ksf(j),kf,kn(1)))**2
     -      +4*dble(conjg(gl(ksf(i),kf,kn(1)))
     -               *gr(ksf(i),kf,kn(1))
     -               *gl(ksf(j),kf,kn(1))
     -               *conjg(gr(ksf(j),kf,kn(1)))))
     -     /((mass(ksf(i))/m0)**2-1d0)
     -     /((mass(ksf(j))/m0)**2-1d0)
   40      continue
   50     continue
          cn=2*nsf*tmpres 
        endif

c... KK DM
        if (mhtype.eq.2) then
         do 60 i=1,nsf
            tmpres=tmpres+16/3.*gyweak**4*(i/2.)**4
     -     /((mass(ksf(i))/m0)**2-1.)**2
   60     continue
          cn=2*nsf*tmpres
        endif

c... user defined
      elseif (mhtype.eq.3) then

        write (*,*) 'WARNING: Choice ''user'' for dsmhset ', 
     -              'not yet implemented in dsmhm2simp.f!'
      endif


      if (cn.lt.0d0) then
        if (memory.ne.idtag) then
          write(*,*) 'WARNING: negative |M|^2 in dsmhm2simp.f for model'
     &               ,idtag, ' !'
c          write(*,*) 'SMtype, cn = ',SMtype,cn
          memory=idtag
        endif
        cn=0d0
      endif


      return

      end





        

