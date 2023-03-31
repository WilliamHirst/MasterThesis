************************************************************************
*     ==========================                                       *
       SUBROUTINE BHINIT( LUN )                                         
*     ==========================                                       *
*                                                                      *
* ((Purpose))                                                          *
*    Initialization program for  histograms and scatter plots.         *
*    This program is called by USERIN.                                 *
* ((Arguments))                                                        *
*    LUN    : logical unit number for message print                    *
* !! Caution!!                                                         *
*    When LUN is set equal to 0, the message print is suppressed.      *
* (( Common /PLOTH/ ))                                                 *
*                                                                      *
*    NW                     : Total number of words of used buffer     *
*                                                                      *
*    NHIST                  : Number of Histograms                     *
*    NSCAT                  : Number of Scat_Plots                     *
*                                                                      *
*   -----------------                                                  *
*     Hashing Table                                                    *
*   -----------------                                                  *
*                                                                      *
*     XHASH(   1,i)      : NH Number of histograms for the i-th class  *
*     XHASH(   2,i) = K  : Serial number of histograms                 *
*              :                     :                                 *
*     XHASH(NH+1,i) = K  : Serial number of histograms                 *
*                     |                                                *
*              MAPL(1,K) = ID  : Histogram ID                          *
*              MAPL(2,K) = IP1 : the 1st pointer to the K-th buffer    *
*              MAPL(3,K) = IP2 : the 2nd pointer to the K-th buffer    *
*              MAPL(4,K) = IP3 : the 3rd pointer to the K-th buffer    *
*                                                                      *
* (( Common /PLOTB/ ))                                                 *
*                                                                      *
*   --------------------                                               *
*     Histogram buffer                                                 *
*   --------------------                                               *
*                                                                      *
*    IP1  = NW + 1                                                     *
*           NW = NW + 281    : Updated NW                              *
*       BUFF( IP1 )          = Xmin                                    *
*       BUFF( IP1 + 1 )      = Xmax                                    *
*       IBUF( IP1 + 2 )      = No. of bins                             *
*       BUFF( IP1 + 3 )      = Bin width                               *
*    IP2  = IP1 + 4                                                    *
*       IBUF(   IP2       )                                            *
*          => IBUF( +  51 )  = No. of sampling points                  *
*       BUFF(   IP2 +  52)                                             *
*          => BUFF( + 103 )  = Sum of Fi for the current IT            *
*       BUFF(   IP2 + 104)                                             *
*          => BUFF( + 155 )  = Sum of Fi**2 for the current IT         *
*       BUFF(   IP2 + 156)                                             *
*          => BUFF( + 207 )  = Sum of Fi for total                     *
*       BUFF(   IP2 + 208)                                             *
*          => BUFF( + 259 )  = Sum of Fi**2 for total                  *
*    IP3  = IP1 + 264                                                  *
*       IBUF( IP3 )          = Tag for spring                          *
*       IBUF( IP3   +  1 )                                             *
*          => IBUF( + 16 )   = Title of this histogram                 *
*                                                                      *
*   --------------------                                               *
*     Scat_Plot buffer                                                 *
*   --------------------                                               *
*                                                                      *
* IP1   = NW + 1                                                       *
*         NW  = NW + 2527                                              *
*       BUFF( IP1 )          = Xmin                                    *
*       BUFF( IP1 + 1 )      = Xmax                                    *
*       IBUF( IP1 + 2 )      = No. of bins for X                       *
*       BUFF( IP1 + 3 )      = Bin width for X                         *
*       BUFF( IP1 + 4 )      = Ymin                                    *
*       BUFF( IP1 + 5 )      = Ymax                                    *
*       IBUF( IP1 + 6 )      = No. of bins for Y                       *
*       BUFF( IP1 + 7 )      = Bin width for Y                         *
* IP2   = IP1 + 8                                                      *
*       BUFF(   IP2       )  = No. of sampling points                  *
*       BUFF(   IP2 +   1 )                                            *
*          => BUFF( +2500 )  = Sum of Fi                               *
* IP3   = IP1 + 2509                                                   *
*       IBUF( IP3 )          = X-Tag for spring                        *
*       IBUF( IP3   +  1 )   = Y-Tag for spring                        *
*       IBUF( IP3   +  2 )                                             *
*          => IBUF( + 17 )   = Title of this histogram                 *
*                                                                      *
*  ((Author))                                                          *
*    S.Kawabata    June '90 at KEK                                     *
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
                                                                        
      COMMON/PLOTLU/ LU                                                 
*                                                                      *
*--------------------------- Entry point ------------------------------*
*                                                                      *
         LU   = LUN                                                     
                                                                        
         NW     = 0                                                     
                                                                        
         DO 50 I = 1, 13                                                
           XHASH(1,I) = 0                                               
           DHASH(1,I) = 0                                               
   50    CONTINUE                                                       
         NHIST    = 0                                                   
         NSCAT    = 0                                                   
         DO 100 I = 1, NHS                                              
           MAPL(1,I)= 0                                                 
  100    CONTINUE                                                       
         DO 200 I = 1, NSC                                              
           MAPD(1,I)= 0                                                 
  200    CONTINUE                                                       
C                                                                       
      RETURN                                                            
      END                                                               
