      SUBROUTINE SCAN(NN,ND,NNMX,N0,NAT,NON,SUB)

      INTEGER*2 NN

      DIMENSION NN(ND,NNMX),IA(3),NA(3)

      DO 2 I=N0,NAT

      NA(1)=1

      IA(1)=I

      CALL SUB(IA,NA,1)

      NJ=NN(I,1)

      DO 2 JJ=2,NJ

      NA(2)=JJ

      J=NN(I,JJ)

      IA(2)=J

      CALL SUB(IA,NA,2)

      IF(NON.EQ.1)GOTO 2

      NK=NN(J,1)

      DO 1 KK=2,NK

      NA(3)=KK

      K=NN(J,KK)

      IA(3)=K

      IF(K.EQ.I)GOTO 1

      CALL SUB(IA,NA,3)

1     CONTINUE

2     CONTINUE

      RETURN

      END
