            FUNCTION DENQD(E,EMX,A,B2,LL,ALP,EPS,TB,NT,NQ,NE,IWK)

      DIMENSION A(LL),B2(LL),TB(NT,4),IWK(LL)

      DOUBLE PRECISION P(2),PP(2),PPP(2),FI(2),FIP(2),XW,SC,DABS

      DOUBLE PRECISION P1,P0,DUM,AN

      N=-1

      DUM=RECWT(E,A,B2,LL,EPS,N,TB,1)

      NN=LL+N

      NR=0

      CALL RECRTS(A,B2,NN,EPS,EMX+EPS,NR,TB,IWK,TB(1,2),TB(1,3))

      IC=1

      NE=0

      DO 5 NQ=1,NR

      IF(IWK(NQ).NE.1)RETURN

      IF(TB(NQ,1).GT.EMX+EPS)GOTO 6

      XV=TB(NQ,1)

      P(1)=1.0D0

      P(2)=XV-A(1)

      PP(1)=0.0D0

      PP(2)=1.0D0

      PPP(1)=0.0D0

      PPP(2)=0.0D0

      FI(1)=0.0D0

      FI(2)=B2(1)

      FIP(1)=0.0D0

      FIP(2)=0.0D0

      I=2

      DO 4 J=2,NN

      XW=XV-A(J)

      I1=3-I

      P(I1)=XW*P(I)-B2(J)*P(I1)

      PP(I1)=XW*PP(I)-B2(J)*PP(I1)+P(I)

      PPP(I1)=XW*PPP(I)-B2(J)*PPP(I1)+2.0*PP(I)

      FI(I1)=XW*FI(I)-B2(J)*FI(I1)

      FIP(I1)=XW*FIP(I)-B2(J)*FIP(I1)+FI(I)

      SC=DABS(P(I1))+DABS(PP(I1))+DABS(PPP(I1))+DABS(FI(I1))

     1 +DABS(FIP(I1))

      DO 3 K=1,2

      P(K)=P(K)/SC

      PP(K)=PP(K)/SC

      PPP(K)=PPP(K)/SC

      FI(K)=FI(K)/SC

3     FIP(K)=FIP(K)/SC

4     I=I1

      I2=3-I

      TB(NQ,2)=FI(I)/PP(I)

      P2=PP(I)*PP(I)

      TB(NQ,4)=P(I2)/PP(I)

5     TB(NQ,3)=(FI(I)*PP(I2)-PP(I)*FI(I2)+TB(NQ,4)

     1 *(PP(I)*FIP(I)-FI(I)*PPP(I)))/P2

      NQ=NR+1

6     NQ=NQ-1

      DEV=ABS(TB(1,1)-E)

      DEN=0.0

      I=2

      IF(NQ.EQ.1)GOTO 8

      DO 7 I=2,NQ

      WD=ABS(TB(I,1)-E)

      IF(DEV.LT.WD)GOTO 8

      DEV=WD

7     DEN=DEN+TB(I-1,3)

      I=NQ+1

8     NE=I-1

      DENQD=(DEN+TB(NE,3)*ALP)/TB(NE,4)

      IF(DEV.GT.10.0*EPS)NE=-NE

      RETURN

      END
