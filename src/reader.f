***********************************************************************
       SUBROUTINE LEER(ITER,NX,NY,NZ,T,ZETA,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************

       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ,LOW1, LOW2
       real T,AAA,BBB,CCC,MAP,ZETA
       INTEGER I,J,K,IX,NL,IR,IRR,N1,N2,N3

       INTEGER FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &         FLAG_W_VELOCITIES
       COMMON /FLAGS/ FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &                FLAG_W_VELOCITIES

*      R4 VARIABLES
       real*4, ALLOCATABLE::SCR4(:,:,:)
C      INTEGER, ALLOCATABLE::SCR4_INT(:,:,:)

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       CHARACTER*9 FILNOM1,FILNOM2, FILNOM4
       CHARACTER*10 FILNOM3

       CHARACTER*24 FIL1,FIL2,FIL4
       CHARACTER*25 FIL3

       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC/ U2,U3,U4,U12,U13,U14

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

*      READING DATA
       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3)
       WRITE(*,*) 'Reading iteration: ',ITER,' ',FILNOM1, ' ', FILNOM2, ' ',
     &            FILNOM3

       FIL1='simu_masclet/'//FILNOM1
       FIL2='simu_masclet/'//FILNOM2
       FIL3='simu_masclet/'//FILNOM3

       WRITE(*,*) FIL1
       OPEN (33,FILE=FIL3,STATUS='UNKNOWN',ACTION='READ')
       OPEN (31,FILE=FIL1,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')


*      pointers and general data (grids file)
       READ(33,*) IRR,T,NL,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ
       DO IR=1,NL
       READ(33,*) IRR,NPATCH(IR)
       READ(33,*)
       WRITE(*,*) 'IR,NPATCH(IR)=', IR,NPATCH(IR)

       IF (IR.NE.IRR) WRITE(*,*)'Warning: fail in restart'
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        READ(33,*) PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        READ(33,*) PATCHX(I),PATCHY(I),PATCHZ(I)
        READ(33,*) AAA,BBB,CCC
        PATCHRX(I)=AAA
        PATCHRY(I)=BBB
        PATCHRZ(I)=CCC
        READ(33,*) PARE(I)
       END DO
       IF (FLAG_VERBOSE.EQ.1) THEN
         WRITE(*,*) 'IR, PATCHX', IR, MINVAL(PATCHX(LOW1:LOW2)),
     &                                MAXVAL(PATCHX(LOW1:LOW2))
         WRITE(*,*) 'IR, PATCHY', IR, MINVAL(PATCHY(LOW1:LOW2)),
     &                                MAXVAL(PATCHY(LOW1:LOW2))
         WRITE(*,*) 'IR, PATCHZ', IR, MINVAL(PATCHZ(LOW1:LOW2)),
     &                                MAXVAL(PATCHZ(LOW1:LOW2))

       END IF

       END DO
       CLOSE(33)

*      BARYONIC (clus file)
       READ(31)
       IR=0
       ALLOCATE(SCR4(0:NMAX+1,0:NMAY+1,0:NMAZ+1))
       SCR4=0.0
C       ALLOCATE(SCR4_INT(0:NMAX+1,0:NMAY+1,0:NMAZ+1))
C       SCR4_INT=0
       N1=NX
       N2=NY
       N3=NZ
       READ(31) !(((U1(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       READ(31) (((SCR4(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       U2(1:NX,1:NY,1:NZ)=SCR4(1:NX,1:NY,1:NZ)
       READ(31) (((SCR4(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       U3(1:NX,1:NY,1:NZ)=SCR4(1:NX,1:NY,1:NZ)
       READ(31) (((SCR4(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       U4(1:NX,1:NY,1:NZ)=SCR4(1:NX,1:NY,1:NZ)
       READ(31) !(((PRES(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       READ(31) !(((POT(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       READ(31) !OPOT
       READ(31) !(((TTT(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       READ(31) !!new: metalicity!! depends on MASCLET version!: TRACER
       READ(31) !!(((SCR4_INT(I,J,K),I=1,N1),J=1,N2),K=1,N3)
c        CR0AMR(1:NX,1:NY,1:NZ)=SCR4_INT(1:NX,1:NY,1:NZ)

       DEALLOCATE(SCR4)
C        DEALLOCATE(SCR4_INT)


       DO IR=1,NL
         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))
         DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))
          SCR4=0.0
C        ALLOCATE(SCR4_INT(NAMRX,NAMRY,NAMRZ))
C        SCR4_INT=0
          READ(31) !(((U11(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
          READ(31) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
             U12(1:N1,1:N2,1:N3,I)=SCR4(1:N1,1:N2,1:N3)
          READ(31) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
             U13(1:N1,1:N2,1:N3,I)=SCR4(1:N1,1:N2,1:N3)
          READ(31) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
             U14(1:N1,1:N2,1:N3,I)=SCR4(1:N1,1:N2,1:N3)
          READ(31) !(((PRES21(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
          READ(31) !(((POT1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
          READ(31) !OPOT
          READ(31) !(((TTT1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
          READ(31) !!new: metalicity!! depends on MASCLET version! TRACER
          READ(31) !!(((SCR4_INT(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
c           CR0AMR1(:,:,:,I)=SCR4_INT(:,:,:)
          READ(31) !!(((SCR4_INT(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
c           SOLAP(:,:,:,I)=SCR4_INT(:,:,:)
          DEALLOCATE(SCR4)
C        DEALLOCATE(SCR4_INT)
         END DO
       END DO

       CLOSE(31)

       RETURN
       END
