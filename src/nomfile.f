***************************************************************
       SUBROUTINE NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3)
***************************************************************
*     Filenames for the input files
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
*     Filenames for the output files
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

*****************************************************************
      SUBROUTINE NOMFILEMACH5(ITER,FILE5)
*****************************************************************
*     Filenames for the Mach number files.
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

      FILE5='MachNum_'//NOM

      RETURN
      END

*****************************************************************
       SUBROUTINE NOMFILE_FILTLEN(ITER,FILE5)
*****************************************************************
*     Filenames for the turbulent outer scale files
*****************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*13 FILE5
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILE5='filtlen_'//NOM

       RETURN
       END


*****************************************************************
       SUBROUTINE NOMFILE_PARTICLES_ERR(ITER,FILE5)
*****************************************************************
*     Filename for the erro statistics
*****************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*21 FILE5
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILE5='error-particles_'//NOM

       RETURN
       END

***************************************************************
       SUBROUTINE NOMFILE_GRIDVARS(ITER,FILNOM1,FILNOM2)
***************************************************************
*     Filenames for the grid outputs
***************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*13 FILNOM1
       CHARACTER*10 FILNOM2
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILNOM1='gridvars'//NOM
       FILNOM2='grids'//NOM

       RETURN
       END

*****************************************************************
       SUBROUTINE NOMFILE3(ITER,FILE5)
*****************************************************************
*     Filenames for the output files (particles)
*****************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*24 FILE5
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILE5='velocity-particles'//NOM

       RETURN
       END
