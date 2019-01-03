      SUBROUTINE RECPER(HOP,VOP,W1,W0,A,B,NW,LLIM,NA,NL,AMAT)

      DIMENSION W1(NW,LLIM),W0(NW,LLIM),A(NA,LLIM),B(NA,LLIM)

     1,AMAT(LLIM,LLIM)

      DOUBLE PRECISION SUM,DSQRT

      NLM1=NL-1

      FAC=1.0/SQRT(B(1,1))

      DO 16 N=1,NLM1

      DO 2 LL=1,LLIM

      L=LLIM-LL+1

      DO 2 I=1,NW

      SUM=0.0

      DO 1 M=1,L

1     SUM=SUM+B(N,M)*W0(I,L-M+1)

2     W0(I,L)=SUM

      CALL HOP(W1,W0,AMAT,NW,LLIM,LLIM)

      A(N,1)=AMAT(1,1)

      IF(LLIM.LE.1)GOTO 6

      DO 3 L=2,LLIM

      A(N,L)=AMAT(1,L)

      DO 3 M=2,L

3     A(N,L)=A(N,L)+AMAT(M,L-M+1)

      CALL VOP(W1,W0(1,2),AMAT,NW,LLIM,LLIM-1)

      DO 5 L=2,LLIM

      LL=L-1

      DO 4 M=1,LL

4     A(N,L)=A(N,L)+AMAT(M,L-M)

5     A(N,L)=A(N,L)/B(N,1)

6     A(N,1)=A(N,1)/B(N,1)

      DO 7 L=1,LLIM

      DO 7 M=1,L

      DUM=A(N,M)

      LL=L-M+1

      DO 7 I=1,NW

7     W0(I,L)=W0(I,L)-DUM*W1(I,LL)

      DO 17 L=1,LLIM

      DO 17 I=1,NW

      DUM=W0(I,L)*FAC

      W0(I,L)=-W1(I,L)*FAC

17    W1(I,L)=DUM

      DO 9 L=1,LLIM

      SUM=0.0D0

      DO 8 M=1,L

      LL=L-M+1

      DO 8 I=1,NW

8     SUM=SUM+W1(I,M)*W1(I,LL)

      IF(L.EQ.1)FAC=1.0/DSQRT(SUM)

9     B(N+1,L)=SUM

      IF(LLIM.LE.1)GOTO 16

      DO 11 L=2,LLIM

      LL=L-1

      IF(LL.LT.2)GOTO 11

      DO 10 M=2,LL

10    B(N+1,L)=B(N+1,L)-B(N+1,M)*B(N+1,L-M+1)/B(N+1,1)

11    B(N+1,L)=B(N+1,L)*0.5

      DO 14 L=2,LLIM

      LL=L-1

      DO 14 M=1,LL

      DUM=B(N+1,L-M+1)/B(N+1,1)

      DO 14 I=1,NW

14    W1(I,L)=W1(I,L)-DUM*W1(I,M)

16    CONTINUE

      RETURN

      END
