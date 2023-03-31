C***********************************************************************
C*=======================                                              *
       SUBROUTINE BSDATE                                                
C*=======================                                              *
C*((Purpose))                                                          *
C*    Changethe format of the time stamp.                              *
C*    This program should be modified according to the machine.        *
C*((Author))                                                           *
C*    S.Kawabata  Nov. '91 at KEK                                      *
C*    For HP      Jul. '92 at KEK                                      *
C***********************************************************************
       COMMON /BDATE/ IDATE(3),ITIME(2)                                 
       COMMON /SLATE/ IS(40)                                            
*            IDATE(1) : year        ITIME(1) : hour                     
*            IDATE(2) : month       ITIME(2) : minute                   
*            IDATE(3) : day                                             
                                                                        
*      CALL UXDATE(IY,IM,ID,IHH,IMM)                                    
       IDATE(1) = MOD(IDATE(1),1900)                                    
*      IDATE(1) = IY                                                    
*      IDATE(2) = IM                                                    
*      IDATE(3) = ID                                                    
*      ITIME(1) = IHH                                                   
*      ITIME(2) = IMM                                                   
       RETURN                                                           
       END                                                              
