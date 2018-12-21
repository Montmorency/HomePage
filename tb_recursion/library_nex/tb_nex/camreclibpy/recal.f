      SUBROUTINE RECAL(HOP,PSI,PMN,M,A,B2,N)

      DIMENSION PSI(M),PMN(M),A(N),B2(N)

      DOUBLE PRECISION SUM,DSQRT

      EXTERNAL HOP

      SUM=B2(1)

      NM1=N-1

      DO 2 J=1,NM1

      B2(J)=SUM

      S=1.0/DSQRT(SUM)

      CALL HOP(PSI,PMN,A(J))

      A(J)=A(J)/B2(J)

      SUM=0.0D0

      DO 1 I=1,M

      DUM=PSI(I)

      PSI(I)=(PMN(I)-A(J)*PSI(I))*S

      PMN(I)=DUM

1     SUM=SUM+PSI(I)*PSI(I)

      S=S*SUM

      DO 2 I=1,M

2     PMN(I)=-PMN(I)*S

      B2(N)=SUM

      RETURN

      END
