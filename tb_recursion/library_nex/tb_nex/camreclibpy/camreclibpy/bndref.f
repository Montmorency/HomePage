      SUBROUTINE BNDREF(DEL,AM,BM2,LL,EPS,AA,RNG,NB,BWK,NBD,IC,NET)

      DIMENSION AM(LL),BM2(LL),AA(NB),RNG(NB),BWK(NBD,3),IC(NBD)

     1 ,P(2,3)

      J=0

      DO 33 I=1,NET

      J=J+1

      IF(IC(I).EQ.0)GOTO 32

      VAL=RECWT(BWK(I,1),AM,BM2,LL,EPS,1,P,1)

      BWK(J,3)=(BWK(I,2)-VAL)*FLOAT(IC(I))

      GOTO 33

32    IT=30

      ERR=ABS(BWK(J,3))

      AE=BWK(I,1)-ERR

      BE=BWK(I,1)+ERR

      ERR=ERR*0.01

      CALL WTMIN(AE,BE,AM,BM2,LL,EPS,ERR,IT,EM,FEM)

      BWK(J,3)=BWK(I,1)-EM

      J=J+1

      BWK(J,3)=FEM-BWK(I,2)

33    CONTINUE

      IB=0

      J=0

      DO 44 I=1,NET

      J=J+1

      IF(IC(I))41,43,42

41    IB=IB+1

      DA=SIGN(DEL,BWK(J,3))

      AA(IB)=AA(IB)+DA

      GOTO 44

42    RNG(IB)=RNG(IB)+SIGN(DEL,BWK(J,3))-DA

      GOTO 44

43    D1=SIGN(DEL,BWK(J,3))

      J=J+1

      D2=SIGN(DEL,BWK(J,3))*0.5

      RNG(IB)=RNG(IB)+D1-D2-DA

      IB=IB+1

      DA=D1+D2

      AA(IB)=AA(IB)+DA

44    CONTINUE

      RETURN

      END
