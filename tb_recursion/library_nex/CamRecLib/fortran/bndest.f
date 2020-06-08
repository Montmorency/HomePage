      SUBROUTINE BNDEST(A,B2,LL,EPS,AA,RNG,NB,EV,FEV,IC,NET,DE,WK,NW)

      DIMENSION A(LL),B2(LL),AA(NB),RNG(NB),EV(NET),FEV(NET),IC(NET)

     1 ,WK(NW,4),P(2,3)

      NPTS=NW

      ALP=1.0

      LM1=LL-1

      EMN=AMIN1(A(1)-B2(2),A(LM1)-1.0)

      EMX=AMAX1(A(1)+B2(2),A(LM1)+1.0)

      NID=LM1-1

      IF(NID.LT.2)GOTO 2

      DO 1 I=2,NID

      EMN=AMIN1(EMN,A(I)-1.0-B2(I+1))

1     EMX=AMAX1(EMX,A(I)+1.0+B2(I+1))

2     DE=EMX-EMN

      EMN=EMN-DE/FLOAT(LL)

      EMX=EMX+DE/FLOAT(LL)

      DE=(EMX-EMN)/FLOAT(NPTS-1)

      E0=EMN

      DO 3 NE=1,NPTS

      WK(NE,1)=E0+FLOAT(NE-1)*DE

3     WK(NE,2)=RECWT(WK(NE,1),A,B2,LL,EPS,1,P,1)

      THU=ALP*B2(1)/FLOAT(LL)

      THL=B2(1)/(6.0*FLOAT(LL)**3)

      DO 4 I=1,5

      NET=NW

      CALL TABAN(WK,WK(1,2),NPTS,THU,THL,WK(1,3),IC,WK(1,4),NET)

      IF(NET.GT.0)GOTO 5

4     THU=(THU+THL)*0.5

      RETURN

5     CALL BNDCRD(WK(1,3),IC,WK(1,4),NET,NB,AA,RNG)

      IF(NB.LE.0)RETURN

      IB=0

      J=0

      DO 7 I=1,NET

      IF(IABS(IC(I)).GT.1)GOTO 7

      J=J+1

      IF(IC(I).EQ.0)GOTO 6

      IC(J)=-IC(I)

      EV(J)=WK(I,3)

      FEV(J)=RECWT(EV(J),A,B2,LL,EPS,1,P,1)

      GOTO 7

6     IC(J)=0

      E=WK(I,3)

      IT=30

      CALL WTMIN(E-DE,E+DE,A,B2,LL,EPS,EPS*10.0,IT,EV(J),FEV(J))

7     CONTINUE

      NET=J

      RETURN

      END
