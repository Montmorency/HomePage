      SUBROUTINE PEEL(CRD,NDIM,NAT,NN,ND,NM,IST,NS,IZP,IZERO,NSH,IW)

      INTEGER*2 NN,IST,IZP,IZERO,IW

      DIMENSION NN(ND,NM),IST(NS),IZP(NAT),IZERO(NAT),IW(NAT)

     1 ,CRD(NDIM,NAT)

      NA=0

      DO 1 I=1,NAT

      IW(I)=0

      IF(IZERO(I).GT.NSH)GOTO 1

      NA=NA+1

      IW(I)=NA

1     CONTINUE

      DO 4 I=1,NAT

      II=IW(I)

      IF(II.EQ.0)GOTO 4

      DO 2 J=1,3

2     CRD(J,II)=CRD(J,I)

      IZERO(II)=IZERO(I)

      IZP(II)=IZP(I)

      JJ=NN(I,1)

      KV=1

      IF(JJ.GE.2) GOTO 38

      JJ=2

      GOTO 37

38    DO 3 J=2,JJ

      K=NN(I,J)

      KK=IW(K)

      IF(KK.EQ.0)GOTO 3

      KV=KV+1

      NN(II,KV)=KK

3     CONTINUE

      JJ=KV+1

      IF(JJ.GT.NM)GOTO 36

37    DO 35 J=JJ,NM

35    NN(II,J)=0

36    NN(II,1)=KV

4     CONTINUE

      NAT=NA

      DO 5 I=1,NS

      J=IST(I)

5     IST(I)=IW(J)

      RETURN

      END
