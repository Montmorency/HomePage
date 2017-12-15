      SUBROUTINE CFGPGN(AA,RNG,WB,NBP1,IC,EPS,A,B2,LM1,WK,NW)

      DIMENSION AA(NBP1),RNG(NBP1),WB(NBP1),A(LM1),B2(LM1)

      DIMENSION WK(NW,2,NBP1)

      DOUBLE PRECISION SA,SB,SAJ,SBJ,DSQRT,DUM

      PI=3.14159265359

      NB=NBP1-1

      IF(IC-2)10,20,1

10    DO 11 I=1,LM1

      TH=FLOAT(LM1-I+1)*PI/FLOAT(LM1+1)

      WK(I,1,1)=COS(TH)

      WK(I,2,1)=2.0*(1.0-WK(I,1,1)*WK(I,1,1))/FLOAT(LM1+1)

11    WK(I,1,1)=(WK(I,1,1)+1.0)*0.5

      GOTO 1

20    DO 21 L=1,LM1

      A(L)=0.5

      AL2=FLOAT((L-1)*(L-1))

21    B2(L)=AL2/(16.0*AL2-4.0)

      B2(1)=1.0

      CALL RECQD(A,B2,LM1,WK,WK(1,2,1),LO,EPS,WK(1,1,2),WK(1,2,2))

      IF(LO.EQ.LM1)GOTO 1

      IC=0

      RETURN

1     SA=0.0D0

      SB=0.0D0

      DO 3 J=1,NB

      SAJ=0.0D0

      JP=J+1

      DO 2 I=1,LM1

      SAJ=SAJ+(WK(I,1,1)*RNG(J)+AA(J))*WK(I,2,1)

      WK(I,1,JP)=0.0

2     WK(I,2,JP)=1.0

      SB=SB+WB(J)

3     SA=SA+SAJ*WB(J)

      B2(1)=SB

      A(1)=SA/SB

      IF(LM1.LE.1)RETURN

      DO 6 L=2,LM1

      ANRM=DSQRT(SB)

      SA=0.0D0

      SB=0.0D0

      DO 5 J=1,NB

      JP=J+1

      SAJ=0.0D0

      SBJ=0.0D0

      DO 4 I=1,LM1

      XX=WK(I,1,1)*RNG(J)+AA(J)

      DUM=WK(I,2,JP)

      WK(I,2,JP)=((XX-A(L-1))*WK(I,2,JP)-B2(L-1)*WK(I,1,JP))/ANRM

      WK(I,1,JP)=DUM/ANRM

      DUM=WK(I,2,JP)*WK(I,2,JP)*WK(I,2,1)

      SBJ=SBJ+DUM

4     SAJ=SAJ+XX*DUM

      SA=SA+SAJ*WB(J)

5     SB=SB+SBJ*WB(J)

      B2(L)=SB

      IF(SB.LT.EPS)GOTO 7

6     A(L)=SA/SB

      L=LM1+1

7     LM1=L-1

      RETURN

      END
