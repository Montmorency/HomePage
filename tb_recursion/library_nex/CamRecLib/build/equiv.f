      LOGICAL FUNCTION EQUIV(V,W)

      REAL DIMENSION V(3),W(3)

      EQUIV=.FALSE.

      DO 1 I=1,3

      IF(ABS(V(I)-W(I)).GT.1.0E-4)RETURN

1     CONTINUE

      EQUIV=.TRUE.

      RETURN

      END
