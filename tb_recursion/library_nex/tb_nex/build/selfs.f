      SUBROUTINE SELFS(I,J,R2,IOVPAR,EM,NE)

      DIMENSION PARM(13),EM(NE,1)

      ID=IOVPAR(I,J,R2,PARM)

      WRITE(6,*) "HL Applying self energy for s."
      WRITE(6,*) PARM(13)

      EM(1,1)=PARM(13)

      RETURN

      END
