      SUBROUTINE NUMD(A,B2,ALAM,LM1,DE)

      DIMENSION A(LM1),B2(LM1)

      P0=0.0

      PP0=0.0

      P1=1.0

      PP1=0.0

      DO 6 I=1,LM1

      P2=(A(I)-ALAM)*P1-B2(I)*P0

      PP2=(A(I)-ALAM)*PP1-B2(I)*PP0-P1

      SC=ABS(P1)+ABS(P2)

      P0=P1/SC

      PP0=PP1/SC

      PP1=PP2/SC

6     P1=P2/SC

      DE=P1/PP1

      RETURN

      END
