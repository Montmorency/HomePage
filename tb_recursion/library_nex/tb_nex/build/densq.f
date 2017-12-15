            FUNCTION DENSQ(E,A,B2,LL,EI)

      DIMENSION A(LL),B2(LL),EI(2),P(2),Q(2)

      PI=3.14159265359

      DENSQ=0.0

      DISC=(E-EI(1))*(E-EI(2))

      IF(DISC.GT.0.0)RETURN

      BI2=(EI(2)-EI(1))*0.25

      BI2=BI2*BI2*2.0

      AIT=-SQRT(-DISC)/BI2

      RT=(E-(EI(1)+EI(2))*0.5)/BI2

      CALL PLYVAL(E,A,B2,LL-1,P,Q)

      T1=P(2)-RT*P(1)*B2(LL)

      T2=AIT*P(1)*B2(LL)

      DENSQ=-AIT*B2(LL)*(Q(2)*P(1)-Q(1)*P(2))/((T1*T1+T2*T2)*PI)

      RETURN

      END
