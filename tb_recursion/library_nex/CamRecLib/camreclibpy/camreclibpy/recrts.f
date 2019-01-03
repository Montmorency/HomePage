      SUBROUTINE RECRTS(A,B2,LM1,EPS,XLIM,N,X,MULT,BI,NI)

      DIMENSION A(LM1),B2(LM1),X(LM1),MULT(LM1),BI(LM1),NI(LM1)

      IF(LM1.GT.1)GOTO 1

      IF(N.EQ.0.AND.A(1).GT.XLIM)RETURN

      X(1)=A(1)

      MULT(1)=1

      N=1

      RETURN

1     MULT(1)=0

C

C ESTIMATE LIMITS OF ROOTS

C

      CN=AMIN1(A(1)-B2(2),A(LM1)-1.0)

      NID=LM1-1

      IF(NID.LT.2)GOTO 20

      DO 2 I=2,NID

2     CN=AMIN1(CN,A(I)-1.0-B2(I+1))

20    AA=CN

      IF(N.EQ.0)GOTO 4

      CN=AMAX1(A(1)+B2(2),A(LM1)+1.0)

      IF(NID.LT.2)GOTO 21

      DO 3 I=2,NID

3     CN=AMAX1(CN,A(I)+1.0+B2(I+1))

21    B=CN

4     IF(N.EQ.0) B=XLIM

      NMAX=AMAX1(ALOG(EPS/(B-AA))/ALOG(0.5)*10,70.0)

      NA=NUMC(A,B2,AA,LM1)

      NB=NUMC(A,B2,B,LM1)

5     NL=NA

      NM=NB

      N=0

      IF (NL.LE.NM) GOTO 19

      NE=NL-NM

      DO 22 I=1,NE

      X(I)=0.0

22    MULT(I)=0

      I=1

      NI(1)=NB

      BI(1)=B

      NCOUNT=0

C

C START BISECTION SEARCH

C

9     C=(AA+B)*0.5

      NC=NUMC(A,B2,C,LM1)

      NCOUNT=NCOUNT+1

      IF(NCOUNT.GE.NMAX)GOTO 16

      IF(NC.LT.NL) GOTO 11

      AA=C

      NA=NC

      GOTO 15

11    IF(NC.LE.NM)GOTO 14

      IF(NI(I).EQ.NC) GOTO 14

      I=I+1

      BI(I)=C

      NI(I)=NC

14    B=C

      NB=NC

15    IF( (NA-NB).GT.1 ) GOTO 36

C DO A FEW NEWTON-RAPHSON STEPS

      C=(AA+B)*0.5

      DO 34 ITR=1,10

      CALL NUMD(A,B2,C,LM1,DE)

      C=C-DE

      IF(C.GT.B .OR. C.LT.AA) GOTO 9

      IF( ABS(DE).LE.EPS) GOTO 35

34    CONTINUE

      GOTO 9

35    XX=AMAX1(C-EPS,AA)

      NXX=NUMC(A,B2,XX,LM1)

      IF( NXX.NE.NA) GOTO 9

      XX=AMIN1(C+EPS,B)

      NXX=NUMC(A,B2,XX,LM1)

      IF( NXX.NE.NB) GOTO 9

      N=N+1

      X(N)=C

      MULT(N)=1

      GOTO 37

36    IF(ABS(AA-B).GT.EPS)GOTO 9

16    N=N+1

      X(N)=(AA+B)*0.5

      MULT(N)=NA-NB

      IF(NCOUNT.GE.NMAX)MULT(N)=-MULT(N)

37    IF(NB.LE.NM)RETURN

      IF(I.LE.1)GOTO 19

18    AA=BI(I)

      NA=NI(I)

      I=I-1

      B=BI(I)

      NB=NI(I)

      NL=NL-IABS(MULT(N))

      NCOUNT=0

      GOTO 9

19    MULT(1)=-MULT(1)

      RETURN

      END
