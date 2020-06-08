      SUBROUTINE RECSUM(AC,BC,NA,NN,NC,A,B,EPS,WK,NW)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION AC(NA,NC),BC(NA,NC),WK(NW),A(NA),B(NA)

      N=IABS(NN)

      IF(NN.LT.0)GOTO 12

      DO 11 I=1,NC

11    AC(N,I)=AC(N-1,I)

12    L1=NC*N

      L2=2*L1

      L3=L2+N

      NL=0

      IC=0

      DO 2 I=1,NC

      CALL RECQD(AC(1,I),BC(1,I),N,WK(IC+1),WK(L1+IC+1),M,

     1EPS,WK(L2+1),WK(L3+1))

      IF(M.NE.N)GOTO 5

      NN=0

      DO 1 J=1,N

      IF(WK(L1+IC+J).GT.0.0)GOTO 1

      WK(L1+IC+J)=0.0

      NN=NN+1

1     CONTINUE

      IC=IC+N

2     NL=MAX0(NL,NN)

      NN=N-NL

      CALL CFGEN(WK,WK(L1+1),IC,EPS,A,B,NN,WK(L2+1))

      RETURN

5     NN=-M

      RETURN

      END
