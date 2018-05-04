      SUBROUTINE SLKODE(DUM,I,J,EM,NGBR)

      DIMENSION DUM(3),EM(5,5)

      DIMENSION DD(3),X(6),X2(6),DM(20),E(20,3)

COMMENT D ELECTRON HAMILTONIAN CALCULATION NGBR comes in as bccbfe

      R2=0.0

      DO 1 L=1,3

      X(L)=DUM(L)

      X2(L)=X(L)*X(L)

1     R2=R2+X2(L)

      IF(R2.LT.1.0E-4)GOTO 3

      R2I=1.0/R2

      RI=SQRT(R2I)

      DO 2 L=1,3

      X(L)=X(L)*RI

      X(L+3)=X(L)

      X2(L)=X2(L)*R2I

2     X2(L+3)=X2(L)

      CALL SKDD(X,X2,I,J,R2,NGBR,EM,5)

      RETURN

C GOTO 3 returns 0s along diagonal of self energy mat.

3     CALL SELFD(I,J,R2,NGBR,EM,5)

      RETURN

      END
