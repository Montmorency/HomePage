      SUBROUTINE BNDWT(AA,RNG,WB,NB,A,B2,LL,EPS,WK,NW,IWK)

      DIMENSION AA(NB),RNG(NB),WB(NB),A(LL),B2(LL),WK(NW,4),IWK(NW)

      IF(NB.LT.2)GOTO 23

      DO 22 I=2,NB

      E=AA(I-1)+RNG(I-1)

      ID=0

      WB(I-1)=0.0

      DO 22 K=1,2

      DUM=DENQD(E,E,A,B2,LL,0.5,EPS,WK,NW,NQ,NE,IWK)

      IF(NE.LE.0)GOTO 35

      JL=NE-ID

      IF(JL.LT.1)GOTO 21

      DO 9 J=1,JL

9     WB(I-1)=WB(I-1)+WK(J,2)

21    ID=1

      E=AA(I)

22    CONTINUE

23    WB(NB)=B2(1)

      IF(NB.LT.2)RETURN

      DO 24 II=2,NB

      I=NB-II+2

      WB(I-1)=WB(I-1)*0.5

24    WB(I)=WB(I)-WB(I-1)

      RETURN

35    NB=0

      RETURN

      END
