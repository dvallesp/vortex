***************************************************************
       SUBROUTINE NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3)
***************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*9 FILNOM1,FILNOM2
       CHARACTER*10 FILNOM3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILNOM1='clus'//NOM
       FILNOM2='cldm'//NOM
       FILNOM3='grids'//NOM

       RETURN
       END

*****************************************************************
       SUBROUTINE NOMFILE2(ITER,FILE5)
*****************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*14 FILE5
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILE5='velocity'//NOM

       RETURN
       END
