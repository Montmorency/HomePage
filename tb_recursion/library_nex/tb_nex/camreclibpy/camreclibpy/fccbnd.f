      INTEGER FUNCTION FCCBND(I,J,R2,DD)

      DIMENSION DD(13)

      IF(R2-3.0)1,3,3

3     FCCBND=0

      RETURN

1     FCCBND=1

      IF(R2.LT.1.0E-4)GOTO 2

      DD(1)=-0.027784

      DD(2)=0.012535

      DD(3)=-0.001554

      RETURN

2     DD(11)=0.0

      RETURN

      END
