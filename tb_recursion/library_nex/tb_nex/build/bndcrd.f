      SUBROUTINE BNDCRD(ET,IC,WTT,NET,NB,AA,RNG)

      DIMENSION ET(NET),IC(NET),WTT(NET),AA(NB),RNG(NB)

      DIMENSION NC(9)

      DO 31 I=1,9

31    NC(I)=0

      WMX=0.0

      DO 32 I=1,NET

      IF(IC(I).EQ.4)WMX=AMAX1(WMX,WTT(I))

      ICC=IC(I)+5

32    NC(ICC)=NC(ICC)+1

      NL=NC(6)+NC(8)

      NR=NC(2)+NC(4)

      IF(NL.NE.NR)GOTO 51

      IF(NL.GT.NB)GOTO 52

      IF(NL.EQ.NB)GOTO 37

      NM=NB-NL

      IF(NM.LT.NC(9))GOTO 34

      DO 33 I=1,NET

      IF(IC(I).EQ.4)IC(I)=0

33    CONTINUE

      GOTO 37

34    DO 36 I=1,NM

      WM=WMX

      DO 35 J=1,NET

      IF(IC(J).NE.4)GOTO 35

      IF(WTT(J).GT.WM)GOTO 35

      WM=WTT(J)

      JJ=J

35    CONTINUE

36    IC(JJ)=0

37    NB=0

      DO 41 II=1,NET

      ICC=IC(II)+5

      GOTO(41,39,41,39,40,38,41,38,41),ICC

38    NB=NB+1

      AA(NB)=ET(II)

      GOTO 41

39    RNG(NB)=ET(II)-AA(NB)

      GOTO 41

40    RNG(NB)=(ET(II)+ET(II-1))*0.5-AA(NB)

      NB=NB+1

      AA(NB)=(ET(II+1)+ET(II))*0.5

41    CONTINUE

      RETURN

51    NB=0

      RETURN

52    NB=-NL

      RETURN

      END
