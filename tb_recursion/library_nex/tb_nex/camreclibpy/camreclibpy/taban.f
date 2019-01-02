      SUBROUTINE TABAN(E,WT,NPTS,THU,THL,ET,IC,WTT,NET)

      DIMENSION E(NPTS),WT(NPTS),ET(NET),IC(NET),WTT(NET)

C

C CONTROL IS MODERATED THROUGHOUT BY THE VARIABLE IREG WHICH

C INDICATES IN WHICH REGION THE PREVIOUS POINT LIES.

C IREG=1   FUNCTION VALUE LESS THAN LOWER THRESHOLD, THL

C IREG=2   FUNCTION VALUE LIES BETWEEN THL AND THU (IN WINDOW)

C IREG=3   FUNCTION VALUE LIES ABOVE UPPER THRESHOLD, THU

C

      I=0

      IREG=2

      IF(WT(1).GT.THU)IREG=3

      IF(WT(1).LT.THL)IREG=1

      INCR=0

      W0=2.0*WT(1)-WT(2)

      IF(WT(1).GT.W0)INCR=4

      IF(WT(1).LT.W0)INCR=-4

      IF(INCR.EQ.0)EFL=E(1)

      DO 10 II=1,NPTS

      GOTO(1,3,7),IREG

1     IF(WT(II).LT.THL)GOTO 9

      IF(WT(II).LT.THU)GOTO 2

      IREG=3

      I=I+1

      ET(I)=(E(II)+E0)*0.5

      IC(I)=3

      GOTO 9

2     IREG=2

      I=I+1

      ET(I)=(E(II)+E0)*0.5

      IC(I)=1

      INCR=4

      GOTO 85

3     NCR=0

      IF(WT(II).LT.W0)NCR=-4

      IF(WT(II).GT.W0)NCR=4

      IF(NCR.EQ.0)GOTO 9

      IF(NCR.EQ.INCR)GOTO 4

      I=I+1

      ET(I)=(E0+EFL)*0.5

      WTT(I)=W0

      IC(I)=NCR

      IF(I.GE.NET)GOTO 11

4     INCR=NCR

      IF(WT(II).LT.THU)GOTO 5

      IREG=3

      I=I+1

      ET(I)=(E(II)+E0)*0.5

      IC(I)=2

      GOTO 9

5     IF(WT(II).GT.THL)GOTO 85

      IREG=1

      I=I+1

      ET(I)=(E(II)+E0)*0.5

      IC(I)=-1

      GOTO 9

7     IF(WT(II).GT.THU)GOTO 9

      IF(WT(II).GT.THL)GOTO 8

      IREG=1

      I=I+1

      ET(I)=(E(II)+E0)*0.5

      IC(I)=-3

      GOTO 9

8     IREG=2

      I=I+1

      ET(I)=(E(II)+E0)*0.5

      IC(I)=-2

      INCR=-4

85    W0=WT(II)

      EFL=E(II)

9     IF(I.GE.NET)GOTO 11

      E0=E(II)

10    CONTINUE

      II=-I

11    NET=-II

      RETURN

      END
