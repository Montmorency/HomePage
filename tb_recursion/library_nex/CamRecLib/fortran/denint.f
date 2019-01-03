            FUNCTION DENINT(E,A,B2,NA,NP,LL,ALP,EPS,WK,IWK,ICODE)

      DIMENSION A(NA,NP),B2(NA,NP),WK(LL,4),IWK(LL)

      SUM=0.0

      DO 2 II=1,NP

      DUM=DENQD(E,E,A(1,II),B2(1,II),LL,ALP,EPS,WK,LL,NQ,NE,IWK)

      IF(NE.LT.1) GOTO 3

      SUM=SUM+(ALP-1.0)*WK(NE,2)

      DO 1 I=1,NE

1     SUM=SUM+WK(I,2)

2     CONTINUE

      ICODE=0

      DENINT=SUM

      RETURN

3     ICODE=-1

      RETURN

      END
