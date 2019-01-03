      SUBROUTINE ADDAT(CRD,ND,NAT,EV,IZP,MM,NN,NND,NM,NGBR,NE,EE,NP

     1,VEC,IW,NED,OVPAR,HCAL)

C      INTEGER*2 IZP,NN,IW,MM

      INTEGER*2 MM(NND,NM),NN(NND,NM),IZP(NAT),IW(2,NED)

      LOGICAL EV

      DIMENSION CRD(ND,NAT),VEC(3,NED),V(3),EE(NP,NP,NED)
C     1,IZP(NAT),NN(NND,NM),IW(2,NED),MM(NND,NM)
CH
C     COMMON /BLKNNM/NNMAT

      EXTERNAL EV,NGBR,OVPAR

      CALL NNCAL(CRD,ND,NAT,IZP,NN,NND,NM,NGBR)

      NE1=1+NE
      WRITE(6,*) 'THIS IS H NM:', NM
      CALL MMCAL(CRD,ND,NAT,NN,NND,NM,EV,IZP,NE,MM,VEC,IW)

      IF(NE1.GT.NE)RETURN

      WRITE(6,*)'NE1',' NE', NE1, NE
      DO 1 K = NE1, NE

      II=IW(1,K)

      IJ=IW(2,K)

1     CALL HCAL(VEC(1,K),II,IJ,EE(1,1,K),OVPAR)

      DO 2 I=1,3

2     V(I)=0.0

      DO 5 I=1,NAT

      IF(IZP(I).GT.0)GOTO 5

      IT=-IZP(I)

      DO 3 K=1,NE

      IF(.NOT.EV(V,VEC(1,K)))GOTO 3

      IF(IT.EQ.IW(1,K))GOTO 4

3     CONTINUE

      NE=NE+1

      K=NE

      DO 6 J=1,3

6     VEC(J,K)=0.0

      WRITE(6,*) 'Doing diagonal part'

      IW(1,K)=IT

      IW(2,K)=IT

      CALL HCAL(VEC(1,K),IT,IT,EE(1,1,K),OVPAR)

4     MM(I,1)=K

5     CONTINUE

      RETURN

      END
