      SUBROUTINE WTMIN(AA,BB,A,B2,LL,EPS,ACC,ITMX,EM,FEM)

      DIMENSION A(LL),B2(LL),P(2,3),E(4),FE(4)

      E(1)=AA

      E(4)=BB

      FE(1)=RECWT(AA,A,B2,LL,EPS,1,P,1)

      FE(4)=RECWT(BB,A,B2,LL,EPS,1,P,1)

      E(2)=(AA+BB)*0.5

      FE(2)=RECWT(E(2),A,B2,LL,EPS,1,P,1)

      IF(FE(2).GT.FE(1))GOTO 6

      IF(FE(2).GT.FE(4))GOTO 6

      IC=3

      DO 4 I=1,ITMX

      E(IC)=(E(IC+1)+E(IC-1))*0.5

      FE(IC)=RECWT(E(IC),A,B2,LL,EPS,1,P,1)

      IF(FE(2).GT.FE(3))GOTO 2

      E(4)=E(3)

      FE(4)=FE(3)

      IF(IC.EQ.2)GOTO 3

      E(3)=E(2)

      FE(3)=FE(2)

      GOTO 3

2     E(1)=E(2)

      FE(1)=FE(2)

      IF(IC.EQ.3)GOTO 3

      E(2)=E(3)

      FE(2)=FE(3)

3     ERR=ABS(E(4)-E(1))

      IF(ERR.LT.ACC)GOTO 5

4     IC=5-IC

      IT=0

5     EM=(E(4)+E(1))*0.5

      FEM=RECWT(EM,A,B2,LL,EPS,1,P,1)

      ITMX=IT

      RETURN

6     IT=1

      IF(FE(1).GT.FE(4))IT=4

      EM=E(IT)

      FEM=FE(IT)

      ITMX=-2+IT/4

      RETURN

      END
