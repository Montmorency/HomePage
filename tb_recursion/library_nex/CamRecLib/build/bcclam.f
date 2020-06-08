      INTEGER FUNCTION BCCLAM(I,J,R2,DD)

C    HL: Rescaling Masuda's first and second nearest neighbour distances and parameters

C    J. Phys. F 14 47-53 (1984).

C    The scaling of the hopping integrals transforms like a power law |R_t-R_s|^-Q

C    With Q being chosen for each metal:

C    Fe:5, Nb:3.3151, Mo:3.57, W:3.255 

C    So that they fit the LSDA vasp calculations I can scale those parameters. But probably simplest to leave them

      REAL LSDA

      DIMENSION DD(13)
    
      LSDA = 1.00

      IF(R2-4.5)1,3,3

3     BCCLAM=0

      RETURN

1     BCCLAM=1

      IF(R2-3.5)2,4,4

4     DD(1)=-0.02253*LSDA

      DD(2)=0.01502*LSDA

      DD(3)=-0.00376*LSDA

      RETURN

2     IF(R2.LT.1.0E-4)GOTO 5

      DD(1)=-0.04624*LSDA

      DD(2)=0.03083*LSDA

      DD(3)=-0.00771*LSDA

      RETURN

5     DD(11)=0.0

      RETURN

      END
