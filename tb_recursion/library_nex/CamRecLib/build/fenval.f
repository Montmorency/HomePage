      FUNCTION FENVAL(AN,A,B2,ND,LL,NC,ERR,EPS,EB,WK,NW,IW,IFT)

      DIMENSION A(ND,NC),B2(ND,NC),EB(2),WK(NW,4),IW(NW),FEB(2)

      ALP=0.5

      ITMX=IFT

      IFT=0

      DO 1 I=1,2

      FEB(I)=-AN

      DO 1 J=1,NC

      ANT=DENQD(EB(I),EB(I),A(1,J),B2(1,J),LL,ALP,EPS,WK,NW,NQ,NE,IW)

      IF(NE.LE.0)RETURN

      FEB(I)=FEB(I)-WK(NE,2)*(1.0-ALP)

      DO 1 K=1,NE

1     FEB(I)=FEB(I)+WK(K,2)

      IFT=-1

      IF(FEB(1)*FEB(2).GT.0)RETURN

      IFT=0

      EF=(EB(1)+EB(2))*0.5

      DO 6 IT=1,ITMX

      IF(ABS(EB(2)-EB(1)).LT.ERR)GOTO 7

      ANTL=-AN

      DANTL=0.0

      DO 3 J=1,NC

      ANT=DENQD(EF,EF,A(1,J),B2(1,J),LL,ALP,EPS,WK,NW,NQ,NE,IW)

      IF(NE.LT.0)RETURN

      DANTL=DANTL+ANT

      ANTL=ANTL-WK(NE,2)*(1.0-ALP)

      DO 3 K=1,NE

3     ANTL=ANTL+WK(K,2)

      IF(ANTL*FEB(1).GT.0.0)GOTO 4

      EB(2)=EF

      FEB(2)=ANTL

      GOTO 5

4     EB(1)=EF

      FEB(1)=ANTL

5     IF(ABS(DANTL).LT.EPS*ABS(ANTL))DANTL=EPS*ANTL

      DE=ANTL/DANTL

      EF=EF-DE

      IF(ABS(DE).LT.ERR)GOTO 7

      IF(EF.LT.EB(1))EF=(EB(1)+EB(2))*0.5

      IF(EF.GT.EB(2))EF=(EB(1)+EB(2))*0.5

6     CONTINUE

      IT=ITMX+1

7     IFT=IT

      FENVAL=EF

      RETURN

      END
