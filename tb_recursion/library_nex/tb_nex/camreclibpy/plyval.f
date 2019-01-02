      SUBROUTINE PLYVAL(E,A,B2,LM1,P,Q)

      DIMENSION A(LM1),B2(LM1),P(2),Q(2)

      P(1)=1.0

      P(2)=E-A(1)

      Q(1)=0.0

      Q(2)=B2(1)

      IF(LM1.LT.2)RETURN

      DO 1 L=2,LM1

      SC=1.0/(ABS(P(2))+ABS(Q(2)))

      DUM=P(2)

      P(2)=((E-A(L))*P(2)-B2(L)*P(1))*SC

      P(1)=DUM*SC

      DUM=Q(2)

      Q(2)=((E-A(L))*Q(2)-B2(L)*Q(1))*SC

1     Q(1)=DUM*SC

      RETURN

      END
