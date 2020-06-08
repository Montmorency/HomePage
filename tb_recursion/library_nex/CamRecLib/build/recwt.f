      FUNCTION RECWT(E,A,B2,LL,EPS,N,P,NS)

      DIMENSION A(LL),B2(LL),P(2,3)

      RECWT=B2(1)

      IF(LL.EQ.1)RETURN

      K=NS

      IF(K.GE.2)GOTO 1

      P(1,1)=1.0

      P(2,1)=E-A(1)

      P(1,2)=0.0

      P(2,2)=1.0

      P(1,3)=0.0

      P(2,3)=B2(1)

      K=2

1     LLIM=LL-IABS(N)

      IF(LLIM.LT.K)GOTO 4

      DO 3 L=K,LLIM

      SC=ABS(P(1,1))+ABS(P(1,3))

      DO 2 J=1,3

      DUM=P(2,J)

      P(2,J)=((E-A(L))*DUM-B2(L)*P(1,J))/SC

2     P(1,J)=DUM/SC

3     P(2,2)=P(2,2)+P(1,1)

4     IF(N)6,5,6

5     RECWT=P(2,3)/P(2,2)

      GOTO 7

6     RECWT=(P(1,1)*P(2,3)-P(2,1)*P(1,3))/(P(1,1)*P(2,2)-P(1,2)

     1 *P(2,1)+P(2,1)*P(2,1)/B2(LL))

7     IF(RECWT.LT.0.0)RECWT=0.0

      IF(N.GE.0)RETURN

      IF(ABS(P(2,1)).LE.EPS*ABS(P(1,1)))RETURN

      N=0

      A(LL)=E-B2(LL)*P(1,1)/P(2,1)

      RETURN

      END
