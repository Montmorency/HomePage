      SUBROUTINE ONION(NN,ND,NM,IZERO,NAT,IST,NNS,IW)

      INTEGER*2 NN,IZERO,IST,IW

      DIMENSION NN(ND,NM),IZERO(NAT),IST(NNS),IW(NAT)

      NS=NNS

      DO 1 I=1,NAT

1     IZERO(I)=0

      DO 2 I=1,NS

      II=IST(I)

      IW(I)=II

2     IZERO(II)=1

      NF=1

      ISH=2

      NCD=NS

3     DO 4 I=NF,NS

      II=IW(I)

      JL=NN(II,1)

      DO 4 J=2,JL

      JJ=NN(II,J)

      IF(IZERO(JJ).NE.0)GOTO 4

      IZERO(JJ)=ISH

      NCD=NCD+1

      IW(NCD)=JJ

      IF(NCD.EQ.NAT)GOTO 5

4     CONTINUE

      ISH=ISH+1

      NF=NS+1

      NS=NCD

      GOTO 3

5     CONTINUE

      RETURN

      END
