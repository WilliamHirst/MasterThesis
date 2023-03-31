**************************************************************
*** SUBROUTINE dsasfere                                    ***
*** computes dW_{ij}/dcostheta                             ***
*** sfermion(i) + sfermion(j)                              ***
*** -> fermion(k1) + fermion(k2)                           ***
*** version to be used if i or j has non-zero lepton #     ***
***                                                        *** 
*** input askin variables: p12,costheta                    ***
*** kpk1 and kpk2 are the fermion code                     ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-08-10                                         ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      SUBROUTINE dsasfere(kpi,kpj,kpk1,kpk2,icase,result)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kpi,kpj,kpk1,kpk2,icase
      integer iii
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
***** equal / non equal particle final state: 
      if(kpk1.eq.kpk2) then
        s34=2.d0
      else  
        s34=1.d0
      endif  
***** fermion type in final state:
      q3=2
      q4=1
***** masses in final state:
      mass3=mass(kpk1)
      mass4=mass(kpk2)
***** define the kinematic variables  
      call dsaskinset2
      call dsaskinset3
***** check if the final state is allowed
      if((Svar-(mass3+mass4)**2)/Svar.lt.1.d-16) then
        result=0.d0
        return
      endif 
***** initialize ASx and ASy to zero
      do iii=1,4
      ASxpl(iii)=(0.d0,0.d0)
      ASxpr(iii)=(0.d0,0.d0)
      enddo
      ASyl=(0.d0,0.d0)
      ASyr=(0.d0,0.d0)
***** then fill them in:
*****
*************** case 1 ****************************************
*****
***** slepton + slepton or sneutrino + sneutrino **************
***** both in the same family *********************************
***** i. s.: particle - particle ******************************
      if(icase.eq.1) then
***** f.s. particle consistent? *******************************
        if(itype(3).ne.itype(1).or.itype(4).ne.itype(1)) then
          write(*,*) 'DS: dsasfere called with wrong particle'
          write(*,*) 'DS: in the final or initial state :'
          write(*,*) pname(kpi),pname(kpj)
     &               ,'and ',pname(kpk1),pname(kpk2)
          stop        
        endif  
***** neutralinos in t-channel **+*****************************
        do iii=1,4
          call dsasscscTCff(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo
***** neutralinos in u-channel **+*****************************
        do iii=1,4
          call dsasscscUCff(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo
*****
*************** case 2 ****************************************
*****
***** sneutrino + slepton* in the same family *****************
      elseif(icase.eq.2) then
***** up-type + down-type-bar, or down-type + up-type-bar *****
***** t-channels: *********************************************
      if(itype(1).eq.itype(3).and.itype(2).eq.itype(4)) then
***** neutralinos in t-channel **+*****************************
        do iii=1,4
          call dsasscscTCff(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo  
***** charginos in u-channel **********************************
        do iii=1,2
            call dsasscscUCff(kpi,kpj,kpk1,kpk2,kcha(iii))
          enddo
      elseif(itype(1).eq.itype(4).and.itype(2).eq.itype(3)) then 
***** neutralinos in u-channel **+*****************************
        do iii=1,4
          call dsasscscUCff(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo  
***** charginos in t-channel **********************************
        do iii=1,2
          call dsasscscTCff(kpi,kpj,kpk1,kpk2,kcha(iii))
        enddo
      else  
        write(*,*) 'DS: dsasfere called with wrong particle'
        write(*,*) 'DS: in the final or initial state :'
        write(*,*) pname(kpi),pname(kpj)
     &               ,'and ',pname(kpk1),pname(kpk2)
        stop        
      endif  
*****
*************** case 3 ****************************************
*****
***** slepton1 + slepton2 (different families) ****************
***** or squark + slepton *************************************
      elseif(icase.eq.3) then  
        if(itype(1).eq.itype(3).and.itype(2).eq.itype(4)) then  
***** neutralinos in t-channel **+*****************************
          do iii=1,4
            call dsasscscTCff(kpi,kpj,kpk1,kpk2,kn(iii))
          enddo
        elseif(abs(itype(1)-itype(3)).eq.1.and.
     &         abs(itype(2)-itype(4)).eq.1) then
***** charginos in t-channel **********************************
          do iii=1,2
            call dsasscscTCff(kpi,kpj,kpk1,kpk2,kcha(iii))
          enddo
        else
***** f. s.: not the same families as i. s. *******************
          write(*,*) 'DS: dsassferfer called with wrong particle'
          write(*,*) 'DS: in the final or initial state :'
          write(*,*) pname(kpi),pname(kpj)
     &               ,'and ',pname(kpk1),pname(kpk2)
          stop        
        endif 
      endif  
*****
*****
***** compute the amplitude squared 
      call dsasff(ampl2)
***** computes dW_{ij}/dcostheta
      result=ampl2*k34/(8.d0*pi*gg1*gg2*s34*dsqrt(Svar))
      return
      end

