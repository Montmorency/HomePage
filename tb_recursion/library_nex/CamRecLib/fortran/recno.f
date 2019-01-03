      SUBROUTINE RECNO(HOP,SOP,U,M,NIT,LS,LL,A,B2,PSI,PMN,EMACH)

      DIMENSION U(M),PSI(M),PMN(M),A(LL),B2(LL)

      DOUBLE PRECISION SUM

      COMMON /BLKNNN/SUM

      NNIT=IABS(NIT)

      LM1=LL-1

      IF(LS.GT.1)GOTO 3

      DO 1 I=1,M

1     PMN(I)=0.0

      IF (NIT) 11,17,15

11    DO 12 I=1,M

12    U(I)=PSI(I)

      DO 13 IT=1,NNIT

13    CALL SOP(U,PSI,U)

      GOTO 17

15    CALL SOP(U,PMN,PSI)

      DO 16 I=1,M

16    PSI(I)=U(I)-PSI(I)

17    SUM=0.0D0

      DO 2 I=1,M

2     SUM=SUM+U(I)*PSI(I)

      B2(1)=SUM

      LS=1

3     DO 8 L=LS,LM1

      DO 4 I=1,M

4     PMN(I)=-B2(L)*PMN(I)

      CALL HOP(U,PMN,A(L))

      A(L)=A(L)/SUM

      ANORM=1.0/DSQRT(SUM)

      DO 5 I=1,M

      DUM=PMN(I)*ANORM

      PMN(I)=PSI(I)*ANORM

5     PSI(I)=DUM-A(L)*PMN(I)

      IF(NNIT.EQ.0)GOTO 10

      DO 6 J=1,NNIT

6     CALL SOP(U,PSI,U)

10    SUM=0.0D0

      DO 7 I=1,M

7     SUM=SUM+U(I)*PSI(I)

      IF(SUM.LT.EMACH) GOTO 9

8     B2(L+1)=SUM

      RETURN

9     LL=-L

      RETURN

      END
