      SUBROUTINE ORPEEL(NSTRT,NORB,NO,MM,NN,ND,ID,EE,NP,NE,NED,MEM)

      INTEBGER*2 MM,NN,MEM

      DIMENSION EE(NP,NP,NED),MM(ND,ID),NN(ND,ID),MEM(3,ID)

      I=NN(NSTRT,1)

      IF(NO.EQ.NP)GOTO 8

      IF(NO.NE.1)GOTO 5

      MEM(2,1)=NE

      DO 1 L=1,I

      II=MM(NSTRT,L)

      NE=NE+1

      IF(NE.GT.NED)GOTO 11

      MEM(1,L)=II

      MM(NSTRT,L)=NE

      DO 1 M=1,NP

      DO 1 N=1,NP

1     EE(M,N,NE)=EE(M,N,II)

      DO 4 L=2,I

      J=NN(NSTRT,L)

      KK=NN(J,1)

      DO 2 K=2,KK

      JJ=NN(J,K)

      IF(JJ.EQ.NSTRT)GOTO 3

2     CONTINUE

      GOTO 13

3     II=MM(J,K)

      NE=NE+1

      IF(NE.GT.NED)GOTO 11

      MEM(2,L)=II

      MEM(3,L)=K

      MM(J,K)=NE

      DO 4 M=1,NP

      DO 4 N=1,NP

4     EE(M,N,NE)=EE(M,N,II)

5     II=MEM(2,1)+1

      DO 55 M=1,NP

55    EE(M,NORB,II)=0.0

      DO 6 L=1,I

      II=MEM(2,1)+L

      DO 6 M=1,NP

6     EE(NORB,M,II)=0.0

      DO 7 L=2,I

      II=II+1

      DO 7 M=1,NP

7     EE(M,NORB,II)=0.0

      RETURN

8     DO 9 L=1,I

9     MM(NSTRT,L)=MEM(1,L)

      DO 10 L=2,I

      J=NN(NSTRT,L)

      K=MEM(3,L)

10    MM(J,K)=MEM(2,L)

      NE=MEM(2,1)

      RETURN

11    WRITE(6,12)

12    FORMAT(' TOO MANY INTERACTION MATRICES FOR STORE IN ORPEEL :

     1 INCREASE DIMENSION AND NED IN CALLING ROUTINE')

      STOP

13    WRITE(6,14)NSTRT,J,NSTRT

14    FORMAT(' INCONSISTENCY IN NEIGHBOUR MAP : ATOM',I6,

     1' IS NOT A NGBR OF ATOM',I6,' WHICH IS A NGBR OF ATOM',I6)

      STOP

      END
