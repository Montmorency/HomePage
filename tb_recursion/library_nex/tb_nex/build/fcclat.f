      SUBROUTINE FCCLAT(CRD,NDIM,IZP,NAT,NX,NY,NZ,NTYPE)

      DIMENSION CRD(NDIM,NAT),IZP(NAT)

      INTEGER*2 IZP

      N=0

      KK=2

      DO 5 K=1,NZ

      II=KK

      RK=FLOAT(K)

      DO 2 J=1,NY

      RJ=FLOAT(J)

      DO 1 I=II,NX,2

      RI=FLOAT(I)

      N=N+1

      IF(N.GT.NAT)GOTO 3

      IZP(N)=NTYPE

      CRD(1,N)=RI

      CRD(2,N)=RJ

1     CRD(3,N)=RK

2     II=3-II

5     KK=3-KK

      NAT=N

      RETURN

3     NAT=0

      WRITE(6,4)I,J,K

4     FORMAT(' NAT TOO SMALL IN FCCLAT- LATTICE GENERATED UP TO',3I4)

      RETURN

      END
