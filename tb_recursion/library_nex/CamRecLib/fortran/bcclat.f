      SUBROUTINE BCCLAT(CRD,NDIM,IZP,NAT,NX,NY,NZ,NTYPE)

      DIMENSION CRD(NDIM,NAT),IZP(NAT)

      INTEGER*2 IZP

      NA=NAT/2

      N=0

      DO 1 I=2,NX,2

      DO 1 J=2,NY,2

      DO 1 K=2,NZ,2

      N=N+1

      IF(N.GT.NA)GOTO 2

      IZP(N)=NTYPE

      CRD(1,N)=FLOAT(I)-1.0

      CRD(2,N)=FLOAT(J)-1.0

      CRD(3,N)=FLOAT(K)-1.0

1     CONTINUE

      NAT=N+N

      GOTO 4

2     WRITE(6,3)I,J,K

      NAT=0

3     FORMAT(' INCREASE NAT IN THE CALL TO BCCLAT - LATTICE GENERATED'

     1 ,' AS FAR AS ',3I4)

4     DO 5 I=1,N

      K=I+N

      IZP(K)=NTYPE

      DO 5 J=1,3

5     CRD(J,K)=CRD(J,I)+1.0

      RETURN

      END
