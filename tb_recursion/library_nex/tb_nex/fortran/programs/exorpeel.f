      PROGRAM EXORPEEL

      INTEGER*2 MM,NN,IZERO,IZP,MEM,IW

      DIMENSION VEC(3,60),CRD(3,100),IZP(100),MEM(3,15),IW(2,60)

      DIMENSION PSI(5,100),PMN(5,100),A(34),B2(34)

      COMMON /BLKREC/NN(100,15),MM(100,15),NAT,NP,EE(5,5,60)

     1,IZERO(100)

      COMMON DUM(3)

      EXTERNAL SLKODE,HOP,FCCBND,EQUIV,PADDOV,PADBND,PADD

C

C THIS PROGRAM COMPUTES THE RECURSION COEFFICIENTS FOR THE LOCAL D.O.S.

C FOR THE 3 P-ORBITALS ON A SURFACE ATOM OF A SMALL F.C.C. CLUSTER OF

C ATOMS INTERACTING VIA THEIR D-ELECTRONS.

C

      NNDIM=100

      INTDIM=15

      NED=60

      LL=11

      NORP=3

      NP=5

      NTYPE=1

      NX=5

      NY=5

      NZ=6

      NPOSN=2

C

C  NNDIM   IS THE FIRST DIMENSION OF ARRAYS NN & MM AND LAST OF CRD

C  INTDIM  IS THE SECOND DIMENSION OF ARRAYS NN & MM

C  NED     IS THE LAST DIMENSION OF IW,VEC AND EE

C  LL     IS THE REQUIRED LENGTH OF THE CONTINUED FRACTION

C  NORP  THE NUMBER OF ORBITALS TO BE PEELED

C  NP IS THE NUMBER OF ORBITALS PER ATOM SITE

C  NTYPE IS THE NUMBER OF 'TYPES' OF ATOM SITE INITIALLY

C  FCCLAT GENERATES A F.C.C. LATTICE OF SIZE NX X NY X NZ

C  NPOSN IS THE NUMBER OF ADD ATOM POSTIONS REQUIRED

C

      NAT=NNDIM

      CALL FCCLAT(CRD,3,IZP,NAT,NX,NY,NZ,NTYPE)

      IF(NAT.EQ.0)STOP

        WRITE(6,4)

4       FORMAT(' EXAMPLE OF RECURSION WITH ORBITAL PEELING',

     1 /20H LATTICE COORDINATES)

        WRITE(6,5)(I,(CRD(J,I),J=1,3),I=1,NAT)

5       FORMAT(4(I5,3F4.1))

        CALL CLCK4X

C

C SETUP ASSEMBLES THE HAMILTONIAN MATRIX USING THE SUPPLIED ROUTINES

C EQUIV,SLKODE, AND FCCBFE, AND THE LIBRARY ROUTINES NNCAL,MMCAL .

C THE INFORMATION IS STORED IN THE COMMON BLOCK /BLKREC/ FOR USE BY HOP

C

C  NN   IS THE NEIGHBOUR MAP GENERATED BY NNCAL

C  MM   IS THE INTERACTION MAP GENERATED BY MMCAL

C  NAT  IS THE NO. OF ATOMS IN THE CLUSTER

C  NP   IS THE NO. OF 'ORBITALS' ON EACH SITE

C  EE   IS THE LIST OF SUBMATRICES

C IZERO IS SUCH THAT IZERO(I)=0 IMPLIES PSI(.,I)=0.0

C

C FOR A CLUSTER OF MORE THAN 100 ATOMS ONLY THE OBVIOUS DIMENSIONS

C ABOVE NEED TO BE CHANGED , LIKEWISE FOR MORE THAN 15

C CONNECTED ATOMS ; THE ONLY OTHER CHANGES  ARE THE CORRESPONDING

C DIMENSION INFORMATION NNDIM AND INTDIM

C

C SETUP IS CALLED TO GENERATE THE D-ELECTRON INTERACTIONS OF THE

C FCC CLUSTER

C

      NM=INTDIM

      CALL SETUP(CRD,3,NAT,EQUIV,NTYPE,IZP,MM,NN,NNDIM,NM,SLKODE

     1,FCCBND,FCCBND,EE,NP,NED,NE,VEC,IW)

      CALL CLCK3F(IT)

      WRITE(6,6)IT

6       FORMAT(' TIME IN SETUP',I6)

      CALL OUT(6,1,1,IZP,IW,VEC,NED,NE,NAT,MM,NN,NNDIM,INTDIM,EE,NP)

C

C ADD ONE ATOM TO THE CLUSTER AND COMPUTE THE P-ELECTRON INTERACTIONS

C USING THE ROUTINE ADDAT

C

      NAT=NAT+1

      CRD(3,NAT)=0.0

      CRD(2,NAT)=2.0

      CRD(1,NAT)=2.0

      IZP(NAT)=-2

      CALL CLCK4X

      NE1=NE

      NM=INTDIM

      CALL ADDAT(CRD,3,NAT,EQUIV,IZP,MM,NN,NNDIM,NM,PADBND,

     1NE,EE,NP,VEC,IW,NED,PADDOV,PADD)

      CALL CLCK3F(IT)

      WRITE(6,7)IT

7     FORMAT(' TIME IN ADDAT',I6)

      CALL OUT(6,NAT,NE1,IZP,IW,VEC,NED,NE,NAT,MM,NN,NNDIM,INTDIM,EE,NP)

      IZP(NAT)=-IZP(NAT)

      WRITE(1,12)LL

      WRITE(1,10)NAT,NPOSN,NORP

10    FORMAT(' ORBITAL PEELING ON A ',I6,' ATOM FCC CLUSTER',/I4

     1,' POSTIONS',I3,' ORBITALS PER POSITION')

C

C ADD A SECOND ATOM TO THE SURFACE IN VARIOUS PLACES

C

      DO 17 NADD=1,NPOSN

      NAT=NAT+1

      CRD(3,NAT)=0.0

      CRD(2,NAT)=2.0

      CRD(1,NAT)=FLOAT(NADD)*0.5+2.5

      IZP(NAT)=-2

      NM=INTDIM

      CALL  ADDAT(CRD,3,NAT,EQUIV,IZP,MM,NN,NNDIM,NM,PADBND,

     1NE,EE,NP,VEC,IW,NED,PADDOV,PADD)

      CALL OUT(6,NAT,0,IZP,IW,VEC,NED,NE,NAT,MM,NN,NNDIM,INTDIM,EE,NP)

C

C  CONSTRUCT INITIAL WAVE FUNCTION NONZERO ON LAST ATOM AND

C  COMPUTE CONTINUED FRACTION

C

      NSTRT=NAT

      WRITE(6,44) NSTRT

44    FORMAT(I6,' ATOM TAKEN AS STARTING POINT'//)

C

C COMPUTES LOCAL D.O.S. RECURSION CFTS. FOR  THE THREE ORBITALS ON THE

C STARTING ATOM. IZERO IS AN ARRAY SET UP TO AVIOD UNNECCESSARY MATRIX

C MULTIPLICATIONS: IZERO(I)=0 IMPLIES PSI(.,I)=0.0  AND

C THE ROUTINE HOP UPDATES THIS WITH EACH APPLICATION OF THE MATRIX.

C

      WRITE(1,8)(CRD(K,NSTRT),K=1,3)

8     FORMAT(' ADD ATOM POSITION',3E13.5)

9     FORMAT(' ORBITAL',I3)

      DO 15 III=1,NORP

        CALL CLCK4X

      DO 11 I=1,NAT

      IZERO(I)=0

      DO 11 J=1,NP

      PSI(J,I)=0.0

11    PMN(J,I)=0.0

      B2(1)=1.0

      A(LL)=0.0

      PSI(III,NSTRT)=1.0

      IZERO(NSTRT)=1

      CALL RECAL(HOP,PSI,PMN,NAT*NP,A,B2,LL)

      WRITE(1,9)III

12    FORMAT(48H CONTINUED FRACTION COEFFICIENTS A(I),B2(I),I=1,,I3)

      WRITE(1,13)(A(I),B2(I),I=1,LL)

13    FORMAT(2E14.6)

        CALL CLCK3F(IT)

        WRITE(6,14)IT

14      FORMAT(35H TIME TO COMPUTE TRIDIAGONALISATION,I6)

      ICODE=III

      IF(III.EQ.NORP)ICODE=NP

C

C 'PEEL' THE ORBITAL FOR WHICH WE HAVE COMPUTED THE CONTINUED FRACTION

C

      CALL ORPEEL(NSTRT,III,ICODE,MM,NN,NNDIM,INTDIM,EE,NP,NE,NED,MEM)

15    CONTINUE

C

C DETACH ATOM READY FOR A RUN WITH A NEW POSITION

C

      MEM(1,1)=1

      DO 16 I=1,NAT

16    IZERO(I)=IZP(I)

      IZERO(NAT)=3

      CALL PEEL(CRD,3,NAT,NN,NNDIM,INTDIM,MEM,1,IZP,IZERO,2,PSI)

17    CONTINUE

      STOP

      END

      SUBROUTINE HOP(PSI,PMN,AN)

      INTEGER*2 MM,NN,IZERO,IDUM

      DOUBLE PRECISION SUM,DUM

      DIMENSION PSI(5,100),PMN(5,100),DUM(5)

      COMMON /BLKREC/NN(100,15),MM(100,15),NAT,NP,EE(5,5,60),IZERO(100)

      COMMON IDUM(100)

C

C HOP SATISFIES THE SPEC. IN RECAL. UNLABELLED COMMON IS USED FOR

C AN INTEGER*2 WORK ARRAY.ALL INFORMATION FOR THE MATRIX DEFINITION

C IS HELD IN THE COMMON BLOCK /BLKREC/ WHICH WAS ASSIGNED BY A CALL TO

C SETUP

C

C   IZERO(I)=0 INDICATES PSI(.,I)=0.0 AND HENCE THAT CERTAIN LINES

C MAY BE SKIPPED TO SAVE TIME.

C

      SUM=0.0

      DO 4 I=1,NAT

C

C THE MATRIX MULTIPLICATION PROCEEDS A LINE AT A TIME BY COMPUTING

C AT EACH SITE I THE EFFECT OF THE HAMILTONIAN OPERATING AT EACH

C OF ITS NEIGHBOURING SITES, WHICH ARE DEFINED BY NN(I,.)

C THE PRODUCT VECTOR IS ACCUMULATED IN DUM AND THE INNER PRODUCT

C <PSI H PSI> IN SUM.

C

      IDUM(I)=IZERO(I)

      DO 11 L=1,NP

11    DUM(L)=0.0

      JSZ=NN(I,1)

      IE=MM(I,1)

      IF(IZERO(I).EQ.0)GOTO 10

C

C MULTIPLY PSI(I) BY THE 'SELF ENERGY' MATRIX (OR DIAGONAL SUBMATRIX)

C

      DO 1 M=1,NP

      DO 1 L=1,NP

1     DUM(L)=DUM(L)+EE(L,M,IE)*PSI(M,I)

10    IF(JSZ.LT.2)GOTO 3

      DO 12 J=2,JSZ

      JJ=NN(I,J)

      IF(IZERO(JJ).EQ.0)GOTO 12

      IE=MM(I,J)

C

C MULTIPLY PSI(.,JJ) BY THE APPROPRIATE SUBMATRIX TO CALCULATE THE

C EFFECT AT SITE I OF THE HAMILTONIAN OPERATING AT SITE JJ.NOTE

C THE PSI(.,I) IS NOW NON-ZERO AND THIS IS RECORDED IN IDUM FOR

C LATER COPYING TO IZERO

C

      DO 2 M=1,NP

      DO 2 L=1,NP

2     DUM(L)=DUM(L)+EE(L,M,IE)*PSI(M,JJ)

      IDUM(I)=1

12    CONTINUE

C

C ACCUMULATE THE REQUIRED INNER PRODUCT AND OVERWRITE PMN

C

3     DO 4 L=1,NP

      SUM=SUM+DUM(L)*PSI(L,I)

4     PMN(L,I)=DUM(L)+PMN(L,I)

      AN=SUM

      DO 14 I=1,NAT

14    IZERO(I)=IDUM(I)

      RETURN

      END

      INTEGER FUNCTION PADBND(I,J,R2,DD)

C

C  FUNCTION DEFINING 'NEIGHBOUR' FOR THE ADD ATOM :

C  TAKES THE VALUE 1 IF I & J ARE NEIGHBOURS , 0 OTHERWISE

C

      PADBND=0

      IF(R2.GT.4.0)RETURN

      PADBND=1

      RETURN

      END

      SUBROUTINE PADD(DUM,I,J,EM,OVP)

      DIMENSION DUM(3),X(6),X2(6),EM(5,5),E(5,5)

      EXTERNAL OVP

C

C  COMPUTES THE OVERLAP MATRICES FOR THE ADD ATOM : P-D AND PP

C  ELECTRONS . SATISFIES THE SPEC OF HCAL FOR SETUP  AND ADDAT.

C





      R2=0.0

      DO 1 L=1,3

      X(L)=DUM(L)

      X2(L)=X(L)*X(L)

1     R2=R2+X2(L)

      DO 6 L=1,5

      DO 6 M=1,5

6     EM(L,M)=0.0

      IF(R2.LT.1.0E-4)GOTO 5

      R2I=1.0/R2

      RI=SQRT(R2I)

      DO 2 L=1,3

      X(L)=X(L)*RI

      X(L+3)=X(L)

      X2(L)=X2(L)*R2I

2     X2(L+3)=X2(L)

      IF(I.EQ.J)GOTO 4

      IF(J.EQ.1)GOTO 3

      CALL SKPD(X,X2,I,J,R2,OVP,E,EM,5)

      RETURN

3     CALL SKPD(X,X2,I,J,R2,OVP,EM,E,5)

      RETURN

4     CALL SKPP(X,X2,I,J,R2,OVP,EM,5)

      RETURN

5     CALL SELFP(I,J,R2,OVP,EM,5)

      RETURN

      END

      INTEGER FUNCTION PADDOV(I,J,R2,PAR)

      DIMENSION PAR(13)

C

C RETURNS THE VALUES OF THE OVERLAP PARAMETERS FOR THE ADD ATOMS.

C

C  PAR(4)= P D SIGMA

C  PAR(5)= P D PI

C  PAR(6)= P P SIGMA

C  PAR(7)= P P PI

C  PAR(12) = SELF ENERGY  OF P ELECTRON

C

      PADDOV=1

      PAR(4)=0.04

      PAR(5)=-0.003

      PAR(6)=0.2

      PAR(7)=-0.002

      PAR(12)=0.5

      RETURN

      END


