      FUNCTION DENCRS(E,A,B2,LL,AA,RNG,WB,NB,AM,BM2)

      DIMENSION A(LL),B2(LL),AA(NB),RNG(NB),WB(NB),AM(LL),BM2(LL)

      DIMENSION P(2),Q(2)

      PI=3.14159265359

      AIF=0.0

      RF=0.0

      DENCRS=0.0

      DO 2 I=1,NB

      DISC=(E-AA(I))*(E-AA(I)-RNG(I))

      WT=8.0*WB(I)/(RNG(I)*RNG(I))

      IF (DISC.LT.0.0)GOTO 1

      DISC=SQRT(DISC)

      IF(E.LT.AA(I))DISC=-DISC

      RF=RF+WT*(E-AA(I)-RNG(I)*0.5-DISC)

      GOTO 2

1     DISC=SQRT(-DISC)

      RF=RF+WT*(E-AA(I)-RNG(I)*0.5)

      AIF=AIF-WT*DISC

2     CONTINUE

      IF(AIF.GE.0.0)RETURN

      CALL PLYVAL(E,AM,BM2,LL-1,P,Q)

      RD=Q(1)-RF*P(1)

      AID=-AIF*P(1)

      DENOM=BM2(LL)*(RD*RD+AID*AID)

      RT=((Q(2)-RF*P(2))*RD-AIF*P(2)*AID)/DENOM

      AIT=(-RD*AIF*P(2)-(Q(2)-RF*P(2))*AID)/DENOM

      CALL PLYVAL(E,A,B2,LL-1,P,Q)

      T1=P(2)-RT*P(1)*B2(LL)

      T2=AIT*P(1)*B2(LL)

      DENCRS=-AIT*B2(LL)*(Q(2)*P(1)-Q(1)*P(2))/((T1*T1+T2*T2)*PI)

      RETURN

      END
