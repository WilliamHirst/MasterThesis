**************************************************************
*** SUBROUTINE dsasfercol                                  ***
*** computes dW_{ij}/dcostheta                             ***
*** sfermion(i) + antisfermion(j)                          ***
*** -> fermion(k1) + antifermion(k2)                       ***
*** version to be used if i or j are both squarks          ***
***                                                        *** 
*** input askin variables: p12,costheta                    ***
*** kpk1 and kpk2 are the fermion code                     ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-08-10                                         ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      SUBROUTINE dsasfercol(kpi,kpj,kpk1,kpk2,icase,result)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kpi,kpj,kpk1,kpk2,icase
      integer iii,jjj
      real*8 result
      real*8 ampl2
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c....check if the final state is kinematically allowed
      if((Svar-(mass(kpk1)+mass(kpk2))**2)/Svar.lt.thstep) then
        result=0.d0
        return
      endif      
***** type for kpk1 and kpk2:
      itype(3)=ivtype(kpk1)
      itype(4)=ivtype(kpk2)
***** non equal particle final state: 
      s34=1.d0
***** fermion type in final state:
      q3=2
      q4=1
***** masses in final state:
      mass3=mass(kpk1)
      mass4=mass(kpk2)
***** define the kinematic variables
      call dsaskinset2
      call dsaskinset3
***** initialize ASx and ASy to zero
      do jjj=1,6
      do iii=1,4
      ASxplc(jjj,iii)=(0.d0,0.d0)
      ASxprc(jjj,iii)=(0.d0,0.d0)
      enddo
      ASylc(jjj)=(0.d0,0.d0)
      ASyrc(jjj)=(0.d0,0.d0)
      enddo
***** then fill them in:
*****
*************** case 1 ****************************************
*****
***** up-type + up-type*, or down-type + down-type*  **********
      if(icase.eq.1) then
***** if particle - antiparticle in f.s. add s-channels *******
      if(kpk1.eq.kpk2) then 
***** Z in s-channel **+***************************************
        call dsasscscsSVffbcol(kpi,kpj,kpk1,kpk2,kz)
***** if kpk1 is not a neutrino: ******************************
        if(dabs(echarg(kpk1)).gt.0.d0) then
***** H01 in s-channel ****************************************
          call dsasscscsSHffbcol(kpi,kpj,kpk1,kpk2,kh1)
***** H02 in s-channel ****************************************
          call dsasscscsSHffbcol(kpi,kpj,kpk1,kpk2,kh2)
***** H03 in s-channel ****************************************
          call dsasscscsSHffbcol(kpi,kpj,kpk1,kpk2,kh3)
***** gamma in s-channel **************************************
          call dsasscscsSVffbcol(kpi,kpj,kpk1,kpk2,kgamma)
        endif  
***** gluon in s-channel **************************************
        if(ncolor(kpk1).gt.2.d0) then
          call dsasscscsSVffbcol(kpi,kpj,kpk1,kpk2,kgluon)
        endif 
      endif 
***** t-channels: *********************************************
      if(itype(1).eq.itype(3).and.itype(2).eq.itype(4)) then
***** neutralinos in t-channel ********************************
        do iii=1,4
          call dsasscscsTCffbcol(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo
***** gluino in t-channel *************************************
        call dsasscscsTCffbcol(kpi,kpj,kpk1,kpk2,kgluin)
      endif  
      if(abs(itype(1)-itype(3)).eq.1
     &   .and.abs(itype(2)-itype(4)).eq.1) then
***** charginos in t-channel **********************************
        do iii=1,2
          call dsasscscsTCffbcol(kpi,kpj,kpk1,kpk2,kcha(iii))
        enddo  
      endif
*****
*************** case 2 ****************************************
*****
***** up-type + down-type*, or down-type + up-type*  **********
      elseif(icase.eq.2) then
***** up-type + down-type-bar, or down-type + up-type-bar *****
***** in f.s. matching the i.s.: add s-channels ***************
      if(itype(3)-itype(4).eq.itype(1)-itype(2)) then
***** W in s-channel **+***************************************
        call dsasscscsSVffbcol(kpi,kpj,kpk1,kpk2,kw)  
***** H^+ in s-channel ****************************************
        call dsasscscsSHffbcol(kpi,kpj,kpk1,kpk2,khc)
      endif
***** t-channels: *********************************************
      if(itype(1).eq.itype(3).and.itype(2).eq.itype(4)) then
***** neutralinos in t-channel ********************************
        do iii=1,4
          call dsasscscsTCffbcol(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo
***** gluino in t-channel *************************************
        call dsasscscsTCffbcol(kpi,kpj,kpk1,kpk2,kgluin)
      endif  
*****
*************** error *****************************************
*****
      else  
        write(*,*) 'DS: dsasfercol called with wrong particles'
        write(*,*) 'DS: in the initial state :'
        write(*,*) pname(kpi),pname(kpj)
        stop        
      endif  
*****
*****
***** compute the amplitude squared 
      call dsasffcol(ampl2)
***** computes dW_{ij}/dcostheta
      result=ampl2*k34/(8.d0*pi*gg1*gg2*s34*dsqrt(Svar))
      return
      end



