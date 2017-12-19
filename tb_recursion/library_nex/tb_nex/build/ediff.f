      FUNCTION EDIFF(EF,A,B2,LL,EPS,WK,MU)

      DIMENSION A(LL),B2(LL),WK(LL,4),MU(LL,2)

      DIMENSION P(2),PP(2),PPP(2),FI(2),FIP(2)

      DOUBLE PRECISION P,PP,PPP,FI,FIP,XW,SC,DABS

      DOUBLE PRECISION P1,P0,DUM,AN

      NN=LL-1

      P1=1.0D0

      P0=0.0D0

      DO 1 I=2,NN

      DUM=(EF-A(I))*P1-B2(I)*P0

      AN=DABS(DUM)+DABS(P1)

      P0=P1/AN

1     P1=DUM/AN

      IF(DABS(P1).LT.EPS*DABS(P0))GOTO 2

      A(LL)=EF-B2(LL)*P0/P1

      NN=LL

2     CONTINUE

      NR=0

      XLIM=EF+EPS

      CALL RECRTS(A(2),B2(2),NN-1,EPS,XLIM,NR,WK,MU,WK(1,2),WK(1,3))

      DO 5 I=1,NR

      IF(MU(I,1).NE.1)WRITE(6,98)WK(I,1),MU(I,1),EF

98    FORMAT(' ***  WARNING ROOT AT ',E16.8,' IS OF MULTIPLICITY',

     1 I6,' IN THE CALL OF EDIFF AT EF=',E16.8/,

     2 ' INCREASING THE ACCURACY IN THE CALL TO EDIFF MAY CURE THIS')

5     CONTINUE

      IC=0

      CALL RECRTS(A,B2,NN,EPS,EF-EPS,IC,WK(1,2),MU(1,2),WK(1,3),WK(1,4))

      DO 6 I=1,IC

      IF(MU(I,2).NE.1)WRITE(6,98)WK(I,2),MU(I,2),EF

6     CONTINUE

      IF(IC.EQ.NR)GOTO 7

      WRITE(6,44)NR,IC

44    FORMAT(' TROUBLE IN EDIFF ',I4,' ZEROS AND',I4,' POLES - TRY',

     1 ' INCREASING THE ACCURACY . EIGENVALUES AS FOLLOWS')

      WRITE(6,45)(WK(I,1),I=1,NR)

      WRITE(6,46)(WK(I,2),I=1,IC)

45    FORMAT(4E26.5)

46    FORMAT(E13.5,3E26.5)

      STOP

7     EDIFF=0.0

      DO 8 J=1,NR

8     EDIFF=EDIFF-WK(J,1)+WK(J,2)

      RETURN

      END
