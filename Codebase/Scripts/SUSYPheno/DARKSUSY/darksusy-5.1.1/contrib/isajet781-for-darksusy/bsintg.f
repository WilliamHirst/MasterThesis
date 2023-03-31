*********************************************************************** 
*                                                                     * 
*    ==========================                                       * 
      SUBROUTINE BSINTG( FXN )                                          
*    ==========================                                       * 
*((Function))                                                         * 
*    Subroutine performs N-dimensional Monte Carlo integration        * 
*    for four vector generation of simulated events                   * 
*                                                                     * 
*       JFLAG = 0 ; First Trial of Defining Grid                      * 
*       JFLAG = 1 ; First Trial of Data Accumulation                  * 
*       JFLAG = 2 ; Second Trial of Defining Grid                     * 
*       JFLAG = 3 ; Second Trial of Data Accumulation                 * 
*                                                                     * 
*    Coded   by S.Kawabata    July 1980 at DESY, Hamburg              * 
*    Last update              March 1994                              * 
*                                                                     * 
*********************************************************************** 
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
                                                                        
      EXTERNAL FXN                                                      
      PARAMETER (MXDIM = 50, NDMX = 50, LENG = 32768)                   
      COMMON /BASE0/ JFLAG,IBASES                                       
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                    
     .               IG(MXDIM),NCALL                                    
      COMMON /BASE2/ ACC1,ACC2,ITMX1,ITMX2                              
      COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT                          
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),      
     .               ND,NG,NPG,MA(MXDIM)                                
      PARAMETER (ITM = 50)                                              
      REAL*4 TIME, EFF, WRONG, TRSLT, TSTD, PCNT                        
      COMMON /BASE5/ ITRAT(ITM,0:1),TIME(ITM,0:2),EFF(ITM,0:1),         
     .               WRONG(ITM,0:1),RESLT(ITM,0:1),ACSTD(ITM,0:1),      
     .               TRSLT(ITM,0:1),TSTD(ITM,0:1),PCNT(ITM,0:1)         
      COMMON /BASE6/ D(NDMX,MXDIM),                                     
     .               ALPH,XSAVE(NDMX,MXDIM),XTI,XTSI,XACC,ITSX          
      REAL*4 STIME                                                      
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF                        
      CHARACTER*80 ERROR                                                
      COMMON /BWARN1/ NERROR                                            
      COMMON /BWARN2/ ERROR(3,3)                                        
*                                                                       
*        INTV = ( 0 / 1 / any ) = ( Batch / Batch(Unix) / Interactive ) 
*        IPNT = ( 0 / any ) = ( IBM Type / Ascii printer )              
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
                                                                        
      REAL*8  X(MXDIM)                                                  
      INTEGER KG(MXDIM),IA(MXDIM)                                       
                                                                        
      COMMON/NINFO/ NODEID, NUMNOD                                      
      REAL*4 TIMEBS,TIMINT,TIMESP,TIME0,RTIME,TIMEB1,TIMEB2,TIMES1      
      COMMON /BTIME1/ TIME0,RTIME,TIMEB1,TIMEB2,TIMES1                  
      COMMON /BTIME2/ TIMEBS(0:2),TIMINT,TIMESP(0:2)                    
*                                                                       
*     Parameters for checking convergency                               
*                                                                       
      DATA ACLMT,FC / 25.0D0, 5.0D0 /                                   
                                                                        
C     REAL*8  TX(2)                                                     
      INTEGER NCNODE(2,512),NPNODE(2,512)                               
C     INTEGER NEFF(2)                                                   
                                                                        
      DATA  ONE/ 1.0D0/, ZERO/0.0D0/, LU / 6/                           
      DATA  N0/0/, N1/1/, HUNDRT/100.0D0/                               
                                                                        
************************************************************************
*                       Initialization Part                             
************************************************************************
*=======================================================================
*          Determine the number of hypercubes NSP                       
*=======================================================================
                                                                        
      XND     = ND                                                      
      NSP     = NG**NWILD                                               
      XJAC    = 1.0D0                                                   
      DO  5 I = 1, NDIM                                                 
         XJAC = XJAC*DX(I)                                              
    5 CONTINUE                                                          
      CALLS   = NSP*NPG                                                 
      DXG     = 1.0D0/NG                                                
      DV2G    = DXG**(2*NWILD)/NPG/NPG/(NPG-1)                          
      DXG     = DXG*XND                                                 
                                                                        
      IF( NSP .EQ. 1 ) THEN                                             
*=======================================================================
*           Determination of the number of sampling points              
*               per node in the single hypercube case                   
*=======================================================================
          MEX     = MOD(NPG,NUMNOD)                                     
          NPERCP  = NPG/NUMNOD                                          
          NPGT    = 0                                                   
          DO  12 NODEX = 1,NUMNOD                                       
             NPGS  = NPGT + 1                                           
             NPGT  = NPGT + NPERCP                                      
             IF( NODEX .LE. MEX ) NPGT = NPGT + 1                       
             NCNODE(1,NODEX) = 1                                        
             NCNODE(2,NODEX) = 1                                        
             NPNODE(1,NODEX) = NPGS                                     
             NPNODE(2,NODEX) = NPGT                                     
   12     CONTINUE                                                      
      ELSE                                                              
*=======================================================================
*          Determination of the number of hypercubes                    
*              per node in many hypercubes case                         
*=======================================================================
          MEX     = MOD(NSP,NUMNOD)                                     
          NPERCP  = NSP/NUMNOD                                          
          NSPT    = 0                                                   
          DO  15 NODEX = 1,NUMNOD                                       
             NSPS  = NSPT + 1                                           
             NSPT  = NSPT + NPERCP                                      
             IF( NODEX .LE. MEX ) NSPT = NSPT + 1                       
             NCNODE(1,NODEX) = NSPS                                     
             NCNODE(2,NODEX) = NSPT                                     
             NPNODE(1,NODEX) = 1                                        
             NPNODE(2,NODEX) = NPG                                      
   15     CONTINUE                                                      
      ENDIF                                                             
*=======================================================================
      NEND    = N0                                                      
      ATACC   = ZERO                                                    
      NERROR  = N0                                                      
      NER1    = N0                                                      
      NER2    = N0                                                      
      NER3    = N0                                                      
      SUMTI   = ZERO                                                    
      SUMTSI  = ZERO                                                    
                                                                        
      IF(JFLAG .EQ. N0 .OR. JFLAG .EQ. N1 ) THEN                        
*-----------------------------------------------------------------------
*        JFLAG = 0  : The first trial of the grid optim. step           
*        JFLAG = 1  : The first trial of the integration step           
*-----------------------------------------------------------------------
         DO 10 J  = N1,NSP                                              
           DXD(J) = ZERO                                                
           DXP(J) = ZERO                                                
   10    CONTINUE                                                       
*       -----------------                                               
         ISTEP   = JFLAG                                                
*       -----------------                                               
         IT1   = N1                                                     
         SI    = ZERO                                                   
         SI2   = ZERO                                                   
         SWGT  = ZERO                                                   
         SCHI  = ZERO                                                   
*       =============                                                   
         CALL BHRSET                                                    
*       =============                                                   
         NSU     = N0                                                   
         SCALLS= ZERO                                                   
      ELSE                                                              
*-----------------------------------------------------------------------
*        JFLAG = 2  : The continuation of the grid optim. step          
*        JFLAG = 3  : The continuation of the integration step          
*-----------------------------------------------------------------------
C        IF( JFLAG .EQ. 2 ) THEN                                        
*           -------------                                               
C            ISTEP  = N0                                                
*           -------------                                               
C        ELSE                                                           
C    .   IF( JFLAG .EQ. 3 ) THEN                                        
*           -------------                                               
C            ISTEP  = N1                                                
*           -------------                                               
C        ELSE                                                           
C                *****************                                      
C                      STOP                                             
C                *****************                                      
C         ENDIF                                                         
C                                                                       
C         IT1   = IT + 1                                                
      ENDIF                                                             
                                                                        
*------- Set the expected accuracy and the max. iteration number -------
                                                                        
      ITMX   = ITMX1                                                    
      ACC    = ACC1*0.01D0                                              
      IF( ISTEP .EQ. N1 ) THEN                                          
         ITMX = ITMX2                                                   
         ACC  = ACC2*0.01D0                                             
      ENDIF                                                             
                                                                        
*-------- Print the title of the convergency behavior table ----------- 
*                  in the interactive mode                              
      IF( INTV .GT. 1 ) THEN                                            
*         -----------------------------------                           
           CALL BSPRNT( LU, 5, ISTEP, IDUM2 )                           
*         -----------------------------------                           
      ENDIF                                                             
      NEGFLG     = 0                                                    
                                                                        
*    =====================                                              
      CALL BSUTIM( 0, 2 )                                               
*    =====================                                              
                                                                        
*********************************************************************   
*               Main Integration Loop                                   
*********************************************************************   
*    ========                                                           
      DO 500  IT = IT1,ITMX                                             
*    ========                                                           
*=======================================================================
*                 Initialization for the iteration                      
*=======================================================================
                                                                        
         SCALLS  = SCALLS + CALLS                                       
         NGOOD   = N0                                                   
         NEGTIV  = N0                                                   
         TI      = ZERO                                                 
         TSI     = TI                                                   
                                                                        
         IF( ISTEP .EQ. N0 ) THEN                                       
             DO 200 J= N1,NDIM                                          
             DO 200 I=1,ND                                              
                D(I,J)= TI                                              
  200        CONTINUE                                                   
         ENDIF                                                          
                                                                        
         NODEX  = NODEID                                                
         IF( NODEID .EQ. 0 )  NODEX = NUMNOD                            
                                                                        
*---------------------------------------------------------------------  
*        Distributing hyper cubes to NumNode nodes                      
*           NCNODE(1,NODEX)   : 1st cube number for the node NODEX      
*           NCNODE(2,NODEX)   : Last cube number for the node NODEX     
*                    NODEX    : node number 1 => NumNode(=0)            
*                    NODEX    : node number 1 => NumNode(=0)            
*---------------------------------------------------------------------  
                                                                        
         NSP1  = NCNODE(1,NODEX)                                        
         NSP2  = NCNODE(2,NODEX)                                        
*                                 Dummy loopfor a parallel processor    
C                                 IF( NSP1 .GT. 1 ) THEN                
C                                     CALL DRLOOP( NDIM*NPG*(NSP1-1) )  
C                                 ENDIF                                 
                                                                        
*=====================================================================  
*      Loop for hypercube from NSP1 to NSP2 in the NodeX-th node        
*=====================================================================  
*       ========                                                        
         DO 400 NCB = NSP1, NSP2                                        
*       ========                                                        
            FB      = 0.0                                               
            F2B     = 0.0                                               
            NP      = NCB - 1                                           
            IF( NWILD .GT. 1 ) THEN                                     
                DO 210 J = 1,NWILD-1                                    
                   NUM   = MOD(NP,MA(J+1))                              
                   KG(J) = NUM/MA(J) + 1                                
  210           CONTINUE                                                
            ENDIF                                                       
            KG(NWILD)     = NP/MA(NWILD) + 1                            
                                                                        
*---------------------------------------------------------------------  
*       If number of hypercubes is only one,                            
*        Distributing sampling points to NumNode nodes                  
*           NPNODE(1,NODEX)   : 1st sample point for the node NODEX     
*           NPNODE(2,NODEX)   : Last sample point for the node NODEX    
*                    NODEX    : node number 1 => NumNode(=0)            
*---------------------------------------------------------------------  
                                                                        
            NPG1  = NPNODE(1,NODEX)                                     
            NPG2  = NPNODE(2,NODEX)                                     
*                                 Dummy loop for a parallel processor   
C                                 IF( NPG1 .GT. 1 ) THEN                
C                                     CALL DRLOOP( NDIM*(NPG1-1) )      
C                                 ENDIF                                 
                                                                        
*=====================================================================  
*          Loop for sampling points from NPG1 to NPG2                   
*                in the single hypercube case                           
*=====================================================================  
*          ========                                                     
            DO 300 NTY = NPG1,NPG2                                      
*          ========                                                     
*---------------------------------------------------------------------  
*        Determine the integration variables by random numbers          
*---------------------------------------------------------------------  
                                                                        
               WGT   = XJAC                                             
               DO 250 J= 1,NDIM                                         
                  IF( J .LE. NWILD ) THEN                               
                      XN  = (KG(J)-DRN(IDUMY))*DXG+1.D0                 
                  ELSE                                                  
                      XN  = ND*DRN(IDUMY)+1.D0                           
                  ENDIF                                                 
                  IA(J)   = XN                                          
                  IAJ     = IA(J)                                       
                  IF( IAJ .EQ. 1) THEN                                  
                      XO  = XI(IAJ,J)                                   
                      RC  = (XN-IA(J))*XO                               
                  ELSE                                                  
                      XO  = XI(IAJ,J)-XI(IAJ-1,J)                       
                      RC  = XI(IAJ-1,J)+(XN-IAJ)*XO                     
                  ENDIF                                                 
                  X(J)    = XL(J)+RC*DX(J)                              
                  WGT     = WGT*XO*XND                                  
  250          CONTINUE                                                 
*-----------------------------------------------------------------------
*                     =======                                           
               FXG  =  FXN(X)*WGT                                       
*                     =======                                           
*-----------------------------------------------------------------------
*             Check the value of the integrand                          
*-----------------------------------------------------------------------
                                                                        
               IF( FXG .NE. 0.0 ) THEN                                  
                   NGOOD = NGOOD + 1                                    
                   IF( ISTEP .EQ. 1 ) THEN                              
                       DXD(NCB) = DXD(NCB) + FXG                        
                       IF( FXG .GT. DXP(NCB) ) DXP(NCB) = FXG           
                   ENDIF                                                
                   IF( FXG .LT. 0.0 ) THEN                              
                       NEGTIV= NEGTIV+ 1                                
                       IF( NEGFLG .EQ. 0 ) THEN                         
                          WRITE(6,9200) IT,NODEID                       
 9200                     FORMAT(1X,                                    
     .                       '******* WARNING FROM BASES ********',     
     .                       '***********',                             
     .                       /1X,'*  Negative FUNCTION at IT =',I3,1X,  
     .                       ', node = ',I3,1X,'*',                     
     .                       /1X,'***********************************', 
     .                       '***********')                             
                          NEGFLG  = 1                                   
                       ENDIF                                            
                   ENDIF                                                
               ENDIF                                                    
                                                                        
*-----------------------------------------------------------------------
*              Accumulation of FXG and FXG*FXG                          
*-----------------------------------------------------------------------
                                                                        
               F2    = FXG*FXG                                          
               FB    = FB + FXG                                         
               F2B   = F2B + F2                                         
                                                                        
               IF( ISTEP .EQ. 0 ) THEN                                  
                   DO 260  J = 1,NDIM                                   
                      D(IA(J),J)= D(IA(J),J)+F2                         
  260              CONTINUE                                             
               ENDIF                                                    
*======                                                                 
  300       CONTINUE                                                    
*======                                                                 
*------------------------------------------- for a parallel processor   
*                                 Dummy loop for a parallel processor   
C                                 IF( NPG2 .LT. NPG ) THEN              
C                                     CALL DRLOOP(NDIM*(NPG-NPG1))      
C                                 ENDIF                                 
*                                 Global sum of FB and F2B              
C                                 IF( NSP .EQ. 1 ) THEN                 
C                                     CALL BSDSUM(  FB, 1 )             
C                                     CALL BSDSUM( F2B, 1 )             
C                                 ENDIF                                 
*-----------------------------------------------------------------------
                                                                        
*-----------------------------------------------------------------------
*         Calculate the estimate and variance in the hypercube          
*-----------------------------------------------------------------------
                                                                        
            F2B   = DSQRT(F2B*NPG)                                      
            F2S   = (F2B-FB)*(F2B+FB)                                   
            TI    = TI+FB                                               
            TSI   = TSI + F2S                                           
                                                                        
*======                                                                 
  400    CONTINUE                                                       
*======                                                                 
*------------------------------------------- for a parallel processor   
*                                 Dummy loop                            
C                                 IF( NSP2 .LT. NSP ) THEN              
C                                     CALL DRLOOP(NDIM*NPG*(NSP-NSP2))  
C                                 ENDIF                                 
                                                                        
*                                 Global sum of efficiency and frequency
*                                     of negative valued function       
C                                 NEFF(1) = NGOOD                       
C                                 NEFF(2) = NEGTIV                      
C                                 CALL BSISUM( NEFF, 2 )                
                                                                        
C                                 TX(1) = TI                            
C                                 TX(2) = TSI                           
C                                 IF( NSP .EQ. 1 ) THEN                 
C                                     CALL BSDSUM(   TX, 2 )            
C                                 ENDIF                                 
                                                                        
*                                 Global sum of grid information        
C                                 IF( ISTEP .EQ. 0 ) THEN               
C                                     NOWORK = NDMX*NDIM                
C                                     CALL BSDSUM(    D, NOWORK )       
C                                 ENDIF                                 
                                                                        
*=====================================================================  
*           Compute Result of this Iteration                            
*=====================================================================  
*--------------------------------------------------------------------   
*           Accumulate the histogram entries                            
*--------------------------------------------------------------------   
*       -------------                                                   
         CALL BHSAVE                                                    
*       -------------                                                   
*--------------------------------------------------------------------   
                                                                        
C        TI     = TX(1)                                                 
C        TSI    = TX(2)                                                 
C        NGOOD  = NEFF(1)                                               
C        NEGTIV = NEFF(2)                                               
                                                                        
         TI    = TI/CALLS                                               
         TSI   = TSI*DV2G                                               
**                                                                      
         IF( TSI .LE. 1.0D-37 ) TSI = 1.0D-37                           
**                                                                      
         TI2   = TI*TI                                                  
                                                                        
         IF( NGOOD .LE. 10 ) THEN                                       
*           --------------------------------                            
             CALL BSPRNT( LU, 9, IDUM1, IDUM2 )                         
*           --------------------------------                            
*            *****************                                          
                   STOP                                                 
*            *****************                                          
                                                                        
         ENDIF                                                          
                                                                        
*--------------------------------------------------------------------   
*               Calculate the cumulative result                         
*--------------------------------------------------------------------   
                                                                        
         WGT   = ONE/TSI                                                
         SI    = SI+TI*WGT                                              
         SWGT  = SWGT+WGT                                               
         SCHI  = SCHI+TI2*WGT                                           
         AVGI  = SI/SWGT                                                
         CHI2A = ZERO                                                   
         IF(IT .GT. N1 ) CHI2A = (SCHI - SI*AVGI)/(IT-.999D0)           
         SD    = DSQRT(ONE/SWGT)                                        
                                                                        
*---------------------------------------------------------------------  
*             Save the results in the buffer                            
*---------------------------------------------------------------------  
                                                                        
         TSI   = DSQRT(TSI)                                             
         ITX         = MOD( IT, ITM)                                    
         IF( ITX .EQ. 0 ) ITX = ITM                                     
         ITRAT(ITX,ISTEP)  = IT                                         
         EFF  (ITX,ISTEP)  = NGOOD/CALLS*HUNDRT                         
         WRONG(ITX,ISTEP)  = NEGTIV/CALLS*HUNDRT                        
         RESLT(ITX,ISTEP)  = AVGI                                       
         ACSTD(ITX,ISTEP)  = SD                                         
         TRSLT(ITX,ISTEP)  = TI                                         
         TACC              = ABS(TSI/TI*HUNDRT)                         
         TSTD (ITX,ISTEP)  = TACC                                       
         PCNT (ITX,ISTEP)  = ABS(SD/AVGI*HUNDRT)                        
                                                                        
*---------------------------------------------------------------------- 
*                  Check cumulative accuracy                            
*---------------------------------------------------------------------- 
                                                                        
         IF( NODEID .EQ. 0 ) THEN                                       
                                                                        
*-------------------  Check cumulative accuracy ----------------------- 
                                                                        
             SDAV  = SD/AVGI                                            
             IF((ABS(SDAV) .LE. ACC)) NEND = N1                         
                                                                        
             IF( ISTEP .EQ. N1 ) THEN                                   
                 IF( TACC .GT. ACLMT ) THEN                             
                     IF( NER1 .EQ. 0 ) THEN                             
                         NERROR = NERROR + 1                            
                         WRITE(ERROR(1,NERROR),9900) NERROR,IT,ACLMT    
 9900                    FORMAT('* (',I1,') Temp. accuracy of it-#',    
     .                         I3,' is too large comparing to',         
     .                         F6.2,' percent.',6X,'*')                 
                         WRITE(ERROR(2,NERROR),9901) TACC,ACLMT         
 9901                    FORMAT('*',8X,'Temp. accuracy (',              
     .                         F7.4,' % )  >>   (',                     
     .                         F7.4,' % )',23X,'*')                     
                         WRITE(ERROR(3,NERROR),9902)                    
 9902                    FORMAT('*',77X,'*')                            
                         NER1  = 1                                      
                     ENDIF                                              
                 ENDIF                                                  
                 IF( IT .GT. 1 ) THEN                                   
                     IF(( TI .GT. AVTI+FDEVI ) .OR.                     
     .                  ( TI .LT. AVTI-FDEVI )      ) THEN              
                          IF( NER2 .EQ. 0 ) THEN                        
                              NERROR = NERROR + 1                       
                              WRITE(ERROR(1,NERROR),9910) NERROR,IT,FC  
 9910                         FORMAT('* (',I1,') Temp. estimate of ',   
     .                        'it-#',I3,' fluctuates more than ',       
     .                               F4.1,'*average-sigma.',6X,'*')     
                              RE = TI                                   
                              CALL BSORDR( RE, FX2, ORDER, IORDR )      
                              RE = TI/ORDER                             
                              RE1 = AVTI                                
                              AC  = FDEVI                               
                              IF( RE1 .GE. AC ) THEN                    
                                  CALL BSORDR( RE1, FX2, ORDR1, IORDR1) 
                              ELSE                                      
                                  CALL BSORDR( AC, FX2, ORDR1, IORDR1)  
                              ENDIF                                     
                              RE1 = AVTI/ORDR1                          
                              AC  = AC/ORDR1                            
                              WRITE(ERROR(2,NERROR),9911) RE,IORDR,     
     .                                          RE1,AC,IORDR1           
 9911                         FORMAT('*        Temp. Estimate (',       
     .                         F10.6,' E',I3,')  >  (',F10.6,'+',F8.6,  
     .                         ' ) E',I3,', or',1X,'*')                 
                              WRITE(ERROR(3,NERROR),9912) RE,IORDR,     
     .                                          RE1,AC,IORDR1           
 9912                         FORMAT('*        Temp. Estimate (',       
     .                         F10.6,' E',I3,')  <  (',F10.6,'-',F8.6,  
     .                         ' ) E',I3,5X,'*')                        
                              NER2 = 1                                  
                          ENDIF                                         
                     ENDIF                                              
                     IF( TSI .GT. FDEVI ) THEN                          
                         IF( NER3 .EQ. 0 ) THEN                         
                             NERROR = NERROR + 1                        
                             WRITE(ERROR(1,NERROR),9920) NERROR,IT,FC   
 9920                        FORMAT('* (',I1,') Error of it-#',         
     .                              I3,' fluctuates more than',F4.1,    
     .                              '*average-sigma.',16X,'*')          
                             RE1 = TSI                                  
                             CALL BSORDR( RE1, FX2, ORDER, IORDR)       
                             RE1 = TSI/ORDER                            
                             AC  = FDEVI                                
                             CALL BSORDR( AC, FX2, ORDR1, IORDR1)       
                             AC  = AC/ORDR1                             
                             WRITE(ERROR(2,NERROR),9921) RE1,IORDR,     
     .                                         AC,IORDR1                
 9921                        FORMAT('*        Temp. Error (',           
     .                         F10.6,' E',I3,')  >  (',F10.6,           
     .                         ' E',I3,')',18X,'*')                     
                             WRITE(ERROR(3,NERROR),9902)                
                             NER3  = 1                                  
                         ENDIF                                          
                     ENDIF                                              
                 ENDIF                                                  
                 SUMTSI = SUMTSI + TSI                                  
                 SUMTI  = SUMTI  + TI                                   
                 AVTSI  = SUMTSI/FLOAT(IT)                              
                 AVTI   = SUMTI/FLOAT(IT)                               
                 FDEVI  = FC*AVTSI                                      
             ENDIF                                                      
         ENDIF                                                          
                                                                        
*------------------------------------------- for a parallel processor   
                                                                        
*                                  Broadcast                            
C                                  CALL BSCAST( NEND, 1 )               
                                                                        
*---------------------------------------------------------------------- 
*        Smoothing the Distribution D(I,J) and refine the grids         
*---------------------------------------------------------------------- 
                                                                        
         IF( ISTEP .LE. N0 ) THEN                                       
             IF( IT .EQ. ITMX ) NEND = N1                               
*           ---------------------                                       
             CALL BSETGV( NEND )                                        
*           ---------------------                                       
         ENDIF                                                          
*       ==========================                                      
         CALL BSUTIM( 0, ISTEP )                                        
*       ==========================                                      
                                                                        
         TIME (ITX,ISTEP)  = TIMINT                                     
         STIME             = TIMINT                                     
                                                                        
*---- Print the convergency behavior table in the interactive mode ---- 
         IF( INTV .GT. 1 ) THEN                                         
*            ---------------------------------                          
              CALL BSPRNT ( LU, 6, ISTEP, IDUM2 )                       
*            ---------------------------------                          
         ENDIF                                                          
                                                                        
         IF( NEND .EQ. N1 ) GO TO 600                                   
                                                                        
*       ======================                                          
         CALL BSUTIM( 0, 2 )                                            
*       ======================                                          
*======                                                                 
  500 CONTINUE                                                          
*======                                                                 
      IT    = IT - N1                                                   
      NEND  = N1                                                        
                                                                        
*********************************************************************** 
*                   Termination of BASES                                
*********************************************************************** 
*======                                                                 
  600 CONTINUE                                                          
*======                                                                 
*---------------------------------------------- For a parallel computer 
                                                                        
*                                 Global sum of histograms              
C                                 CALL BHSUM                            
*                                 Global sum of probabilities           
C                                 CALL BSDSUM(  DXD, NSP )              
*                                 Global sum of the max.value in each HC
C                                 CALL BSDSUM(  DXP, NSP )              
                                                                        
                                                                        
*======================= End of the step ? ============================ 
                                                                        
      IF( NEND .EQ. N1 ) THEN                                           
          IF( INTV .GT. 1 ) THEN                                        
*            ---------------------------------                          
              CALL BSPRNT ( LU, 7, IDUM1, IDUM2 )                       
*            ---------------------------------                          
          ENDIF                                                         
          IF( ISTEP .EQ. N0) THEN                                       
              JFLAG   = N1                                              
              ITG     = IT                                              
          ELSE                                                          
              JFLAG   = N0                                              
              ITF     = IT                                              
          ENDIF                                                         
      ENDIF                                                             
*    ======================                                             
       CALL BSUTIM( 0, 2 )                                              
*    ======================                                             
                                                                        
      RETURN                                                            
      END                                                               
