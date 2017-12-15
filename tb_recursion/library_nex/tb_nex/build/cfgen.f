      SUBROUTINE CFGEN(X,W,N,EPS,A,B2,LM1,WK)

      DIMENSION X(N),W(N),A(LM1),B2(LM1),WK(N,2)

      DOUBLE PRECISION SA,SB,DSQRT,DUM

      SA=0.0D0

      SB=0.0D0

      DO 1 I=1,N

      SB=SB+W(I)

      SA=SA+X(I)*W(I)

      WK(I,1)=0.0

1     WK(I,2)=1.0

      B2(1)=SB

      A(1)=SA/SB

      IF(LM1.LT.2)RETURN

      DO 3 L=2,LM1

      ANRM=DSQRT(SB)

      SA=0.0D0

      SB=0.0D0

      DO 2 I=1,N

      DUM=WK(I,2)

      WK(I,2)=((X(I)-A(L-1))*WK(I,2)-B2(L-1)*WK(I,1))/ANRM

      WK(I,1)=DUM/ANRM

      DUM=WK(I,2)*WK(I,2)*W(I)

      SB=SB+DUM

2     SA=SA+X(I)*DUM

      B2(L)=SB

      IF(SB.LT.EPS)GOTO 4

3     A(L)=SA/SB

      L=LM1+1

4     LM1=L-1

      RETURN

      END
