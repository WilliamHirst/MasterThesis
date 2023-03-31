CDECK  ID>, STRADD.
C
      FUNCTION STRADD(PREFIX,EXT)
C
C     This function adds character string EXT to string PREFIX
C     to make new character string STRADD
C
      CHARACTER*128 PREFIX,STRADD
      CHARACTER*6 EXT
      INTEGER NPRE,I,J,JJ
      NPRE=0
      DO 100 I=1,LEN(STRADD)
        STRADD(I:I)=' '
        IF(PREFIX(I:I).NE.' ') THEN
          STRADD(I:I)=PREFIX(I:I)
          NPRE=NPRE+1
        ENDIF
100   CONTINUE
C      write(6,*) 'EXT=',EXT
      DO 110 J=1,LEN(EXT)
        JJ=NPRE+J
        STRADD(JJ:JJ)=EXT(J:J)
110   CONTINUE
      RETURN
      END
