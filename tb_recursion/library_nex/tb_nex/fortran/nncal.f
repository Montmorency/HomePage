      SUBROUTINE NNCAL(CRD,NDIM,NAT,IZP,NN,ND,NM,NGBR)

      INTEGER*2 IZP,NN

      DIMENSION CRD(NDIM,NAT),DUM(13),IZP(NAT),NN(ND,NM)

      NNMAX=0

      IADD=1

      IF(IZP(1).LT.0)NN(1,1)=1

      DO 11 II=1,NAT

      IF(IZP(II).LT.0)GOTO 12

11    CONTINUE

      NN(1,1)=1

      DO 1 I=1,NAT

      DO 1 J=2,NM

1     NN(I,J)=0

      IADD=0

      II=2

12    IF(II.LE.1)II=2

50    FORMAT(' FROM NNCAL')

      DO 5 I=II,NAT

      IF(IZP(I)*IADD.GT.0)GOTO 5

      NN(I,1)=1

      IIP=IZP(I)

      IIP=IABS(IIP)

      ILJ=I-1

      DO 4 J=1,ILJ

      JJP=IZP(J)

      JJP=IABS(JJP)

      R2=0.0

      DO 2 L=1,3

      DUM(L)=CRD(L,I)-CRD(L,J)

2     R2=R2+DUM(L)*DUM(L)

      ID=NGBR(IIP,JJP,R2,DUM)

      IF(ID.EQ.0)GOTO 4

      ID=NN(I,1)+1

      NN(I,1)=ID

      NN(I,ID)=J

      NNMAX=MAX0(NNMAX,ID)

      ID=NN(J,1)+1

      NN(J,1)=ID

      NN(J,ID)=I

      NNMAX=MAX0(NNMAX,ID)

      IF(NNMAX.GT.NM)GOTO 33

3     FORMAT(20H TOO MANY NEIGHBOURS)

4     CONTINUE

5     CONTINUE

      NM=NNMAX

      RETURN

33    WRITE(6,3)

      WRITE(6,38)I

38    FORMAT(24H NEIGHBOUR MAP AS FAR AS,I6,7HTH SITE)

      DO 39 L=1,NAT

      MMM=NN(L,1)

39    WRITE(6,40)L,(NN(L,M),M=1,MMM)

40    FORMAT(21I4)

      STOP

      END
