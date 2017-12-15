      SUBROUTINE SETUP(CRD,ND,NAT,EV,NTYPE,IZP,MM,NN,NND,NM,HCAL,NGBR

     1,IOVPAR,EE,NP,NED,NE,VEC,IW)

      INTEGER*2 MM(NND,NM),NN(NND,NM),IZP(NAT),IW(2,NED)

      LOGICAL EV

      DIMENSION CRD(ND,NAT),VEC(3,NED),EE(NP,NP,NED)

      EXTERNAL NGBR,IOVPAR,EV

C

C COMPUTE NEIGHBOUR MAP

C

      CALL NNCAL(CRD,ND,NAT,IZP,NN,NND,NM,NGBR)

C

C CONSTRUCT HAMILTONIAN MAP MM AND EE

C

      NE=NED

      CALL MMCAL(CRD,ND,NAT,NN,NND,NM,EV,IZP,NE,MM,VEC,IW)

      IF(NE.EQ.0)GOTO 4

      DO 1 I=1,NAT

      II=IZP(I)+NE

1     MM(I,1)=II

C

C  LOAD HAMILTONIAN MATRICES GENERATED BY EACH DISTINCT VECTOR

C

      DO 2 K=1,NE

      II=IW(1,K)

      IJ=IW(2,K)

2     CALL HCAL(VEC(1,K),II,IJ,EE(1,1,K),IOVPAR)

      IF(NE+NTYPE.GT.NED)GOTO 4

      DO 3 KK=1,NTYPE

      K=NE+KK

C

C COPY 'SELF ENERGY' MATRICES

C

      DO 8 I=1,3

8     VEC(I,K)=0.0

      IW(1,K)=KK

      IW(2,K)=KK

      CALL HCAL(VEC(1,K),KK,KK,EE(1,1,K),IOVPAR)

3     CONTINUE

      NE=NE+NTYPE

      RETURN

4     WRITE(6,5)NE,KK

5     FORMAT(' INCREASE DIMENSONS OF EE ETC. AS',I4,'+',I4,

     1' INTERATION MATRICES ARE REQUIRED')

      STOP

      END
