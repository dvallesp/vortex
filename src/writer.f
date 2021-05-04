***********************************************************************
       SUBROUTINE WRITE_DIVROT(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the divergence and each component of the rotational
*     of the velocity field to a file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      CHARACTER*30 FILERR5
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /DIVERGENCE/ DIVER0, DIVER

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)
      real*4 scrvar1, scrvar2

*     OPEN THE OUTPUT FILE
      OPEN(25,FILE=FILERR5,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &     POSITION='APPEND')

*     WRITE GENERAL DATA
      scrvar1 = T
      scrvar2 = ZETA
      WRITE(25) ITER, scrvar1, scrvar2, NL


*     WRITE THE DIVERGENCE FIELD
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ) = DIVER0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=DIVER(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO

      DEALLOCATE(SCR4)

*     WRITE THE ROTATIONAL FIELDS (x, y, z)
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ) = ROTAX_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ) = ROTAY_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ) = ROTAZ_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=ROTAX_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=ROTAY_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=ROTAZ_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO

      DEALLOCATE(SCR4)

      CLOSE(25)

      END


***********************************************************************
       SUBROUTINE WRITE_POTENTIALS(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,
     &            NPATCH,PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the scalar and the vector potentials to a file
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      CHARACTER*30 FILERR5
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /DIVERGENCE/ DIVER0, DIVER

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)

*     OPEN THE OUTPUT FILE
      OPEN(25,FILE=FILERR5,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &     POSITION='APPEND')

      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

*     WRITE THE SCALAR POTENTIAL
      SCR4(1:NX,1:NY,1:NZ) = DIVER0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=DIVER(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO

      DEALLOCATE(SCR4)

*     WRITE THE VECTOR POTENTIAL (x, y, z)
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ) = ROTAX_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ) = ROTAY_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ) = ROTAZ_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRX,NAMRX))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=ROTAX_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=ROTAY_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=ROTAZ_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO


      DEALLOCATE(SCR4)

      CLOSE(25)

      END


***********************************************************************
       SUBROUTINE WRITE_TOTALVELOCITY(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,
     &                             NPATCH,PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the initial velocity field to a file (DEPRECATED)
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      CHARACTER*30 FILERR5
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES

      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)

*     OPEN THE OUTPUT FILE
      OPEN(25,FILE=FILERR5,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &     POSITION='APPEND')

*     WRITE THE TOTAL VELOCITY FIELD (x, y, z)
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ)=U2(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=U3(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=U4(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=U12(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=U13(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=U14(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO

      DEALLOCATE(SCR4)

      CLOSE(25)

      END



***********************************************************************
       SUBROUTINE WRITE_VELOCITIES(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,
     &                             NPATCH,PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the total, the compressive and the rotational velocities
*     to a file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      CHARACTER*30 FILERR5
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES
*     ROTA has been recycled as the rotational velocity
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P

*      original, total velocity
       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC_ORIGINAL/ U2,U3,U4,U12,U13,U14

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)

*     OPEN THE OUTPUT FILE
      OPEN(25,FILE=FILERR5,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &     POSITION='APPEND')

*     WRITE THE TOTAL VELOCITY
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ)=U2(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=U3(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=U4(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         SCR4(1:N1,1:N2,1:N3)=U12(1:N1,1:N2,1:N3,I)
         WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
         SCR4(1:N1,1:N2,1:N3)=U13(1:N1,1:N2,1:N3,I)
         WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
         SCR4(1:N1,1:N2,1:N3)=U14(1:N1,1:N2,1:N3,I)
         WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

       END DO
      END DO

      DEALLOCATE(SCR4)

*     WRITE THE COMPRESSIONAL VELOCITY
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ)=U2P(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=U3P(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=U4P(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=U12P(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=U13P(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=U14P(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO

      DEALLOCATE(SCR4)

*     WRITE THE ROTATIONAL VELOCITY
      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))

      SCR4(1:NX,1:NY,1:NZ)=ROTAX_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=ROTAY_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      SCR4(1:NX,1:NY,1:NZ)=ROTAZ_0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))

      DO IR=1,NL

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          SCR4(1:N1,1:N2,1:N3)=ROTAX_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=ROTAY_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          SCR4(1:N1,1:N2,1:N3)=ROTAZ_1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)

        END DO
      END DO

      DEALLOCATE(SCR4)

      CLOSE(25)

      END


**********************************************************************
       SUBROUTINE WRITE_FILTLEN(FILERR5,NX,NY,NZ,ITER,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ,L0,L1)
***********************************************************************
*     Writes the filter length to a separate file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      CHARACTER*30 FILERR5
      INTEGER NX, NY, NZ, NL, ITER
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      real L0(1:NMAX,1:NMAY,1:NMAZ)
      real L1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)

*     OPEN THE OUTPUT FILE
      OPEN(25,FILE=FILERR5,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &     POSITION='APPEND')

*     WRITE THE 'COHERENCE' LENGTH

      ALLOCATE(SCR4(NMAX,NMAY,NMAZ))
      SCR4(1:NX,1:NY,1:NZ) = L0(1:NX,1:NY,1:NZ)
      WRITE(25) (((SCR4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      DEALLOCATE(SCR4)

      ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))
      DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          SCR4(1:N1,1:N2,1:N3)=L1(1:N1,1:N2,1:N3,I)
          WRITE(25) (((SCR4(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
        END DO
      END DO
      DEALLOCATE(SCR4)

      CLOSE(25)

      END

**********************************************************************
      SUBROUTINE WRITE_GRID_PARTICLES(NL,NX,NY,NZ,NPATCH,PATCHNX,
     &                          PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &                          PATCHRX,PATCHRY,PATCHRZ,PARE,CR0AMR,
     &                          CR0AMR1,SOLAP)
***********************************************************************
*     Writes the GRIDS and CR0AMR/SOLAP variables for the created AMR
*     structure
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'vortex_parameters.dat'

      INTEGER NL,NX,NY,NZ
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER cr0amr(NMAX,NMAY,NMAZ)
      INTEGER cr0amr1(NAMRX,NAMRY,NAMRZ,NPALEV)
      INTEGER solap(NAMRX,NAMRY,NAMRZ,NPALEV)

      real u1(1:NMAX,1:NMAY,1:NMAZ)
      real u11(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /dens/ u1,u11

      INTEGER IX,JY,KZ,I,LOW1,LOW2,IR,IPATCH,J,K,N1,N2,N3

      CHARACTER*13 FILNOMGRIDVARS
      CHARACTER*10 FILNOMGRIDS

      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      CALL NOMFILE_GRIDVARS(ITER,FILNOMGRIDVARS,FILNOMGRIDS)

      write(*,*) 'Writing grids data', filnomgrids
      OPEN (23,FILE='output_files/'//FILNOMGRIDS,STATUS='UNKNOWN')
      WRITE(23,*) ITER,' ',0.0,' ',NL,' ',0.0,' ',0.0
      WRITE(23,*) 0.0
      IR=0
      WRITE(23,*) IR,0,0,NX,NY,NZ
      DO IR=1,NL
       WRITE(23,*) IR,' ',NPATCH(IR),' ',0,' ',0,' ',0
       WRITE(23,*) '----------------- within level=',IR,' -----------'
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        WRITE(23,*) PATCHNX(I),' ',PATCHNY(I),' ',PATCHNZ(I)
        WRITE(23,*) PATCHX(I),' ',PATCHY(I),' ',PATCHZ(I)
        WRITE(23,*) PATCHRX(I),' ',PATCHRY(I),' ',PATCHRZ(I)
        WRITE(23,*) PARE(I)
        END DO
       END DO
       CLOSE(23)

       OPEN(99,FILE='output_files/'//FILNOMGRIDVARS,STATUS='UNKNOWN',
     &     FORM='UNFORMATTED')

       write(99) (((u1(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       write(99) (((cr0amr(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       do ipatch=1,sum(npatch(0:nl))
        n1=patchnx(ipatch)
        n2=patchny(ipatch)
        n3=patchnz(ipatch)
        write(99) (((u11(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
        write(99) (((cr0amr1(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
        write(99) (((solap(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
       end do
       CLOSE(99)

       RETURN
       END


*********************************************************************
      SUBROUTINE WRITE_PARTICLES(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                           PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                           PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &                           MASAP,U2DM,U3DM,U4DM,NPART,LADO0)
***********************************************************************
*     Writes the GRIDS and CR0AMR/SOLAP variables for the created AMR
*     structure
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'vortex_parameters.dat'

      INTEGER NL,NX,NY,NZ
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV),LADO0
      INTEGER PARE(NPALEV)

      REAL*4 RXPA(NDM),RYPA(NDM),RZPA(NDM),
     &        U2DM(NDM),U3DM(NDM),U4DM(NDM),MASAP(NDM)

      INTEGER LIHAL(NDM),LIHAL_IX(NDM),LIHAL_JY(NDM),LIHAL_KZ(NDM)
      REAL SCRPART(NDM)

*     INPUTS FROM COMMON MODULES
*     ROTA has been recycled as the rotational velocity
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P

*      original, total velocity
      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_ORIGINAL/ U2,U3,U4,U12,U13,U14

      INTEGER IX,JY,KZ,I,LOW1,LOW2,IR,IP,NPARTTOT

      CHARACTER*24 FILNOM

      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      NPARTTOT=SUM(NPART(0:NLEVELS))
      CALL NOMFILE3(ITER,FILNOM)
      OPEN(99,FILE='output_files/'//FILNOM,STATUS='UNKNOWN',
     &     FORM='UNFORMATTED')

      WRITE(99) NPARTTOT
      WRITE(99) (RXPA(I),I=1,NPARTTOT)
      WRITE(99) (RYPA(I),I=1,NPARTTOT)
      WRITE(99) (RZPA(I),I=1,NPARTTOT)
      WRITE(99) (U2DM(I),I=1,NPARTTOT)
      WRITE(99) (U3DM(I),I=1,NPARTTOT)
      WRITE(99) (U4DM(I),I=1,NPARTTOT)
      WRITE(99) (MASAP(I),I=1,NPARTTOT)

      CALL PLACE_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ)


      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            U2(1:NX,1:NY,1:NZ),U12(1:NAMRX,1:NAMRY,1:NAMRZ,:),
     &            SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            U3(1:NX,1:NY,1:NZ),U13(1:NAMRX,1:NAMRY,1:NAMRZ,:),
     &            SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            U4(1:NX,1:NY,1:NZ),U14(1:NAMRX,1:NAMRY,1:NAMRZ,:),
     &            SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            U2P(1:NX,1:NY,1:NZ),U12P(1:NAMRX,1:NAMRY,1:NAMRZ,:),
     &            SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            U3P(1:NX,1:NY,1:NZ),U13P(1:NAMRX,1:NAMRY,1:NAMRZ,:),
     &            SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            U4P(1:NX,1:NY,1:NZ),U14P(1:NAMRX,1:NAMRY,1:NAMRZ,:),
     &            SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            ROTAX_0(1:NX,1:NY,1:NZ),
     &            ROTAX_1(1:NAMRX,1:NAMRY,1:NAMRZ,:),SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            ROTAY_0(1:NX,1:NY,1:NZ),
     &            ROTAY_1(1:NAMRX,1:NAMRY,1:NAMRZ,:),SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,
     &            ROTAZ_0(1:NX,1:NY,1:NZ),
     &            ROTAZ_1(1:NAMRX,1:NAMRY,1:NAMRZ,:),SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CLOSE(99)

      RETURN
      END
