      PROGRAM EXRECAL

      INTEGER*2 MM,NN,IZERO,IZP,IW

      DIMENSION CRD(3,1500),VEC(3,20),DUM(3),PSI(5,1500),PMN(5,1500)
     1,A(34), B2(34),IZP(1500),IW(2,20)

      COMMON /BLKREC/NN(1500,20),MM(1500,20),NAT,NP,EE(5,5,20)
     1,IZERO(1500)

      EQUIVALENCE (IZP(1),PMN(1,1))

      EQUIVALENCE (PSI(1,1),CRD(1,1))

      COMMON VEC,IW

      EXTERNAL SLKODE,HOP,BCCBFE,EQUIV, SADDOV, SADBND, SADD

C

C THIS PROGRAM COMPUTES THE RECURSION COEFFICIENTS FOR THE LOCAL D.O.S.

C FOR THE 5 D-ORBITALS ON A CENTRAL ATOM OF A SMALL B.C.C. CLUSTER OF

C ATOMS. THIS SHOULD PROVIDE A 'BLUEPRINT' FOR OTHER APPLICATIONS OF THE

C RECURSION LIBRARY; IN PARTICULAR THE ROUTINES HOP AND SET  NEED VERY

C LITTLE MODIFICATION TO BE APPROPRIATE FOR MANY ELECTRONIC D.O.S

C CALCULATIONS, AND CHANGES SHOULD BE ABLE TO BE RESTRICTED TO THE

C DIMENSION AND OUTPUT STATEMENTS.AS THE ARRAYS PSI AND PMN ARE REFERRED

C TO WITH A SINGLE SUBSCRIPT IN RECAL, THEIR FIRST DIMENSION

C SHOULD EQUAL THE NUMBER OF ORBITALS PER SITE, OTHERWISE THE BASIC

C PROGRAM IS SUITABLE FOR ANY NUMBER OF ORBITALS PER SITE.

C IN THIS EXAMPLE THE ARRAYS CRD AND IZP ARE NOT NEEDED ONCE THE

C HAMILTONIAN HAS BEEN SET UP IN /BLKREC/ AND THEY ARE THEREFORE

C EQUIVALENCED TO PSI AND PMN TO SAVE SPACE. THE 'TYPE' ARRAY IZP IS

C PROVIDED TO ENABLE INTERACTIONS BETWEEN GEOMETRICALLY EQUIVALENT

C SITES TO BE DIFFERENT, AS MAY BE REQUIRED IN LATTICES WITH TWO OR MORE

C DIFFERENT ATOM TYPES,OR IN ALLOWING SURFACE INTERACTIONS TO BE

C DIFFERENT FROM THOSE IN THE BULK.

C

      LL=11

      NED=20

      NNDIM=1500

C      NNMX=15
      NNMX=20

      NP=5

      NTYPE=2

      NX=6

      NY=6

      NZ=6

C

C  LL IS THE REQUIRED LENGTH OF THE CONTINUED FRACTION

C  NED IS THE DIMENSION OF ARRAYS IW,VEC AND EE

C  NNDIM IS THE SECOND DIMENSION OF ARRAYS NN & MM AND PSI,PMN & CRD

C  NNMX IS THE MAXIMUM NUMBER OF NEIGHBOUR INTERACTIONS

C  NP IS THE NUMBER OR ORBITALS PER SITE AND FIRST 2 DIMENSIONS OF EE

C  NTYPE IS THE NUMBER OF 'TYPES' OF ATOM SITE

C  BCCLAT GENERATES A B.C.C. LATTICE OF SIZE NX X NY X NZ

C

      NAT=NNDIM

C     CALL BCCLAT(CRD,3,IZP,NAT,NX,NY,NZ,NTYPE)
      CALL BCCLAT(CRD,3,IZP,NAT,NX,NY,NZ,1)

      IF(NAT.EQ.0)STOP

        WRITE(6,4)

4       FORMAT(' OUTPUT FROM EXAMPLE RUN OF RECURSION LIBRARY - RECAL',

     1 //20H LATTICE COORDINATES)

        WRITE(6,5)(I,(CRD(J,I),J=1,3),I=1,NAT)

5       FORMAT(4(I7,3F7.1))

CH        CALL CLCK4X

C

C SETUP GENERATES THE COMMON BLOCK /BLKREC/ DEFINING THE

C HAMILTONIAN MATRIX

C

C  THE VARIABLES GENERATED AND PUT IN THE COMMON /BLKREC/

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

C CONNECTED ATOMS

      NM=NNMX

      CALL SETUP(CRD,3,NAT,EQUIV,NTYPE,IZP,MM,NN,NNDIM,NM,SLKODE,BCCBFE

     1,BCCBFE,EE,NP,NED,NE,VEC,IW)

CH      CALL CLCK3F(IT)

        WRITE(6,6)IT

6       FORMAT(14H TIME IN SETUP,I6)

      CALL OUT(6,1,1,IZP,IW,VEC,NED,NE,NAT,MM,NN,NNDIM,NM,EE,NP)

C HERE GOES ADDING ADD ATOMS (COPYING EXORPEEL)

      NAT=NAT+1

      CRD(3,NAT)=2.0

      CRD(2,NAT)=2.25

      CRD(1,NAT)=2.5

      IZP(NAT)=-2

      NE1 = NE

      NM = NNMX

      NORP = 1

      NHYD = 1

      CALL ADDAT(CRD,3,NAT,EQUIV,IZP,MM,NN,NNDIM,NM,SADBND,
     1NE,EE,NP,VEC,IW,NED,SADDOV,SADD)


      CALL OUT(6,NAT,NE1,IZP,IW,VEC,NED,NE,NAT,MM,NN,NNDIM,NM,EE,NP)

      IZP(NAT)=-IZP(NAT)

C     WRITE(1,10)NAT,NPOSN,NORP

C  10    FORMAT(' 1S ORBITAL IN A ',I6,' ATOM BCC CLUSTER',/I4

C    1,' POSTIONS',I3,' ORBITALS PER POSITION')

C


C

C  CONSTRUCT INITIAL WAVE FUNCTION NONZERO ON ATOM 4,3,3 AND

C  COMPUTE CONTINUED FRACTION

C

      DUM(1)=2.0

      DUM(2)=2.0

      DUM(3)=2.0

      DO 41 NSTRT=1,NAT

      IF(CRD(1,NSTRT).EQ.DUM(1).AND.CRD(2,NSTRT).EQ.DUM(2).AND.

     1  CRD(3,NSTRT).EQ.DUM(3))GOTO 43

41    CONTINUE

      WRITE(6,42)

42    FORMAT(' NO SUCH ATOM TO START')

      STOP

43    CONTINUE

      WRITE(6,44) NSTRT

44    FORMAT(I6,' ATOM TAKEN AS STARTING POINT'//)

C

C COMPUTES LOCAL D.O.S. RECURSION CFTS. FOR  THE FIVE ORBITALS ON THE

C STARTING ATOM. IZERO IS AN ARRAY SET UP TO AVIOD UNNECCESSARY MATRIX

C MULTIPLICATIONS: IZERO(I)=0 IMPLIES PSI(.,I)=0.0  AND

C THE ROUTINE HOP UPDATES THIS WITH EACH APPLICATION OF THE MATRIX.

C IN THIS CASE (BECAUSE OF SYMMETRY) IT WOULD HAVE BEEN SUFFICIENT

C TO COMPUTE ONLY THE 2 DISTINCT DENSITY FUNCTIONS. THE NUMERICAL

C DISCREPENCIES IN THE ANSWERS ARE CAUSED BY LACK OF FULL SYMMETRY

C IN THE CLUSTER AND STARTING STATES.FOR THE LOCAL D.O.S. SUMMED OVER

C ALL ORBITALS ON THE STARTING ATOM , THE STARTING STATE (1,1,1,1,1)

C ON THAT ATOM COULD HAVE BEEN USED, AND B2(1) SET EQUAL TO 5.0

C

      WRITE(1,10)

10    FORMAT(' DOS FOR A SMALL BCC CLUSTER')

      DO 15 III=1,NP

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

      WRITE(1,12)LL

12    FORMAT(5X,48H CONTINUED FRACTION COEFFICIENTS A(I),B2(I),I=1,,I3)

      WRITE(1,13)(A(I),B2(I),I=1,LL)

13    FORMAT(5X,2E14.6)

CH        CALL CLCK3F(IT)

        WRITE(6,14)IT

14      FORMAT(35H TIME TO COMPUTE TRIDIAGONALISATION,I6)

15    CONTINUE

      STOP

      END

      SUBROUTINE HOP(PSI,PMN,AN)

      INTEGER*2 MM,NN,IZERO,IDUM

      DOUBLE PRECISION SUM,DUM

      DIMENSION PSI(5,1500),PMN(5,1500),DUM(5)

      COMMON /BLKREC/NN(1500,20),MM(1500,20),NAT,NP,EE(5,5,20)
     1,IZERO(1500)

      COMMON IDUM(1500)

C

C HOP SATISFIES THE SPEC. IN RECAL. UNLABELLED COMMON IS USED FOR

C AN INTEGER*2 WORK ARRAY.ALL INFORMATION FOR THE MATRIX DEFINITION

C IS HELD IN THE COMMON BLOCK /BLKREC/ WHICH IS ASSIGNED IN SET

C WHERE A DESCRIPTION MAY ALSO BE FOUND.

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

      INTEGER FUNCTION SADBND(I,J,R2,DD)

C

C  FUNCTION DEFINING 'NEIGHBOUR' FOR THE ADD ATOM :

C  TAKES THE VALUE 1 IF I & J ARE NEIGHBOURS , 0 OTHERWISE

C

      SADBND=0

      IF(R2.GT.0.5)RETURN

      SADBND=1

      RETURN

      END


      SUBROUTINE SADD(DUM,I,J,EM,OVP)

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

      CALL SKSD(X,X2,I,J,R2,OVP,E,EM,5)

      RETURN

3     CALL SKSD(X,X2,I,J,R2,OVP,EM,E,5)

      RETURN

4     CALL SKSS(X,X2,I,J,R2,OVP,EM,5)

      RETURN

5     CALL SELFS(I,J,R2,OVP,EM,5)

      RETURN

      END

      INTEGER FUNCTION SADDOV(I,J,R2,PAR)

      DIMENSION PAR(13)

C

C RETURNS THE VALUES OF THE OVERLAP PARAMETERS FOR THE ADD ATOMS

C HL ADAPATED FROM EXORPEEL

C

C  PAR(8) is this sd\sigma i think so?

C  PAR(10) looks like ss selfenergy.

C  PAR(4)= P D SIGMA

C  PAR(5)= P D PI

C  PAR(6)= P P SIGMA

C  PAR(7)= P P PI

C  PAR(12) = SELF ENERGY  OF P ELECTRON

C

      SADDOV=1

      PAR(8) = -0.2

      PAR(10) = -0.5

      PAR(13) = -0.5

C     PAR(4)=0.04

C     PAR(5)=-0.003

C     PAR(6)=0.2

C     PAR(7)=-0.002

C     PAR(12)=-0.3

      RETURN

      END


