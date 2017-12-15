      INTEGER FUNCTION BCCBFE(I,J,R2,DD)

C

C  MICHAEL YOU'S IBONDS FOR A BCC LATTICE

C  FIRST AND SECOND NEAREST NEIGHBOUR DISTANCES

C  AND THE PETTIFOR PARAMETERS DD.

C

      DIMENSION DD(13)

      IF(R2-4.5)1,3,3

3     BCCBFE=0

      RETURN

1     BCCBFE=1

      IF(R2-3.5)2,4,4

4     DD(1)=-0.03195

      DD(2)=0.02130

      DD(3)=-0.00537

      RETURN

2     IF(R2.LT.1.0E-4)GOTO 5

      DD(1)=-0.06560

      DD(2)=0.04373

      DD(3)=-0.01093

      RETURN

5     DD(11)=0.0

      RETURN

      END
