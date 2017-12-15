            FUNCTION NUMC(A,B2,ALAM,LM1)

      DIMENSION A(LM1),B2(LM1)

7     NU=0

      P0=0.0

      P1=1.0

      DO 6 I=1,LM1

      P2=(A(I)-ALAM)*P1-B2(I)*P0

      SC=ABS(P1)+ABS(P2)

      IF(P2*P1)5,14,4

14    IF(P0*P2)4,5,5

4     NU=NU+1

5     P0=P1/SC

6     P1=P2/SC

      NUMC=NU

      RETURN

      END
