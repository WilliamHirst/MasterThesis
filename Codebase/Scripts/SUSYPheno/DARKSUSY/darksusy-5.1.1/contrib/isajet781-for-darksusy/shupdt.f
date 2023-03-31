************************************************************************
*    ====================                                              *
      SUBROUTINE SHUPDT                                                 
*    ====================                                              *
* ((Function))                                                         *
*     To update histograms and scatter plots with unit weight.         *
*   The bin number to be updated is marked by XHFILL and DHFILL.       *
* ((Author))                                                           *
*     S.Kawabata  June '90 at KEK                                      *
*                                                                      *
************************************************************************
                                                                        
      PARAMETER ( NHS = 50, NSC = 50 )                                  
      INTEGER*4 XHASH,DHASH,NHIST,MAPL,IFBASE,NSCAT,MAPD                
      COMMON/PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),        
     .              NHIST, MAPL(4,NHS),                                 
     .              NSCAT, MAPD(4,NSC),                                 
     .              NW                                                  
      COMMON /PLOTB/ IBUF( 281*NHS + 2527*NSC )                         
      REAL*4         BUFF( 281*NHS + 2527*NSC )                         
      EQUIVALENCE (IBUF(1),BUFF(1))                                     
                                                                        
      IF( NHIST .GT. 0 ) THEN                                           
         DO 150   IHIST   = 1, NHIST                                    
            IP3       = MAPL(4,IHIST)                                   
            IX        = IBUF(IP3)                                       
            IF( IX .GE. 0 ) THEN                                        
                IP       = IX + MAPL(3,IHIST) + 52                      
                BUFF(IP) = BUFF(IP) + 1.                                
                                                                        
                IBUF(IP3)  = -1                                         
            ENDIF                                                       
  150    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      IF( NSCAT .GT. 0 ) THEN                                           
         DO 250   ISCAT   = 1, NSCAT                                    
            IP3         = MAPD(4,ISCAT)                                 
            IX          = IBUF(IP3)                                     
            IF( IX .GT. 0 ) THEN                                        
                IP1   = MAPD(2,ISCAT)                                   
                MXBIN = IBUF(IP1+2)                                     
                MYBIN = IBUF(IP1+6)                                     
                IP2       = MAPD(3,ISCAT)                               
                IBUF(IP2) = IBUF(IP2) + 1                               
                IY        = IBUF(IP3+1)                                 
                IF( IX .GT. 0 .AND. IX .LE. MXBIN .AND.                 
     .              IY .GT. 0 .AND. IY .LE. MYBIN ) THEN                
                    IP       = IX + MXBIN*(IY-1) + IP2                  
                    BUFF(IP) = BUFF(IP) + 1.0                           
                ENDIF                                                   
                IBUF(IP3)   =  0                                        
                IBUF(IP3+1) =  0                                        
           ENDIF                                                        
C                                                                       
  250    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
