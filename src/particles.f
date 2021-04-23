************************************************************************
      SUBROUTINE CREATE_MESH(NX,NY,NZ,NL_MESH,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,
     &            U4DM,MASAP,NPART,LADO0)
************************************************************************
*     Creats a mesh hierarchy for the given particle distribution
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL_MESH
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL*4 RXPA(NDM),RYPA(NDM),RZPA(NDM),
     &       U2DM(NDM),U3DM(NDM),U4DM(NDM),MASAP(NDM)
      REAL LADO0

*     COMMON VARIABLES
      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

*     LOCAL VARIABLES
      INTEGER PLEV(NDM)
      !REAL,ALLOCATABLE::U1(:,:,:)
      !REAL,ALLOCATABLE::U11(:,:,:,:)
      INTEGER,ALLOCATABLE::CR0(:,:,:)
      INTEGER,ALLOCATABLE::CR01(:,:,:,:)
      INTEGER,ALLOCATABLE::CONTA1(:,:,:)
      INTEGER,ALLOCATABLE::CONTA11(:,:,:,:)
      REAL MAP,XL,YL,ZL,DXPA,DYPA,DZPA
      INTEGER I,IX,JY,KZ,REFINE_THR,REFINE_COUNT,BOR,MIN_PATCHSIZE
      INTEGER INI_EXTENSION,NBIS,IRPA,BORAMR,LOW1,LOW2
      INTEGER INMAX(3),I1,I2,J1,J2,K1,K2,N1,N2,N3,IR,MARCA,IPATCH

!     hard-coded parameters (for now, at least)
      REFINE_THR=2
      BOR=6
      BORAMR=0
      INI_EXTENSION=2 !initial extension of a patch around a cell (on each direction)
      MIN_PATCHSIZE=18 !minimum size (child cells) to be accepted

      MAP=MAXVAL(MASAP(1:SUM(NPART)))

      PLEV=0
!$OMP PARALLEL DO SHARED(NPART,PLEV,MAP,MASAP),PRIVATE(I),DEFAULT(NONE)
      DO I=1,SUM(NPART)
       PLEV(I)=LOG(MAP/MASAP(I)+.5)/LOG(8.)
      END DO
      WRITE(*,*) 'Particle levels: min and max values:', MINVAL(PLEV),
     &           MAXVAL(PLEV)

      XL=-FLOAT(NMAX)*DX/2.
      YL=-FLOAT(NMAY)*DY/2.
      ZL=-FLOAT(NMAZ)*DZ/2.

*     FIRST LEVEL OF REFINEMENT ========================================
      IR=1
      DXPA=DX/(2.0**IR)
      DYPA=DY/(2.0**IR)
      DZPA=DZ/(2.0**IR)
      !ALLOCATE(U1(NMAX,NMAY,NMAZ))
      ALLOCATE(CONTA1(NMAX,NMAY,NMAZ))
      ALLOCATE(CR0(NMAX,NMAY,NMAZ))

!$OMP PARALLEL DO SHARED(CONTA1,NX,NY,NZ),PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       CONTA1(IX,JY,KZ)=0
      END DO
      END DO
      END DO

!$OMP PARALLEL DO SHARED(NPART,RXPA,RYPA,RZPA,XL,YL,ZL,DX,DY,DZ,
!$OMP+                   NX,NY,NZ,PLEV,REFINE_THR),
!$OMP+            PRIVATE(I,IX,JY,KZ), DEFAULT(NONE)
!$OMP+            REDUCTION(+: CONTA1)
      DO I=1,SUM(NPART)
       IX=INT((RXPA(I)-XL)/DX)+1
       JY=INT((RYPA(I)-YL)/DY)+1
       KZ=INT((RZPA(I)-ZL)/DZ)+1
       IF (IX.LT.1) IX=1
       IF (IX.GT.NX) IX=NX
       IF (JY.LT.1) JY=1
       IF (JY.GT.NY) JY=NY
       IF (KZ.LT.1) KZ=1
       IF (KZ.GT.NZ) KZ=NZ

       !U1(IX,JY,KZ)=U1(IX,JY,KZ)+MASAP(I)

       IF (PLEV(I).EQ.0) THEN
        CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+1
       ELSE
        CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+REFINE_THR
       END IF
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,BOR,CONTA1,CR0),
!$OMP+            PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       IF(IX.LE.BOR.OR.IX.GE.NX-BOR+1.OR.
     &    JY.LE.BOR.OR.JY.GE.NY-BOR+1.OR.
     &    KZ.LE.BOR.OR.KZ.GE.NZ-BOR+1) THEN
         CONTA1(IX,JY,KZ)=0
       END IF
       CR0(IX,JY,KZ)=CONTA1(IX,JY,KZ)
      END DO
      END DO
      END DO

      !WRITE(*,*) 'TOTAL DM MASS: ',SUM(U1*9.18E18)
      !WRITE(*,*) 'PARTICLE COUNT CHECK: ',SUM(CONTA1),SUM(NPART)
      REFINE_COUNT=COUNT(CR0.GE.REFINE_THR)
      WRITE(*,*) 'REFINABLE CELLS:', REFINE_COUNT

      IPATCH=0

      DO WHILE (REFINE_COUNT.GT.0.AND.IPATCH.LT.NPALEV) !--------------
       INMAX=MAXLOC(CR0)
       IX=INMAX(1)
       JY=INMAX(2)
       KZ=INMAX(3)
       !IF (CONTA1(IX,JY,KZ).LT.REFINE_THR) EXIT

       I1=IX-INI_EXTENSION
       I2=IX+INI_EXTENSION
       J1=JY-INI_EXTENSION
       J2=JY+INI_EXTENSION
       K1=KZ-INI_EXTENSION
       K2=KZ+INI_EXTENSION

       N1=2*(I2-I1+1)
       N2=2*(J2-J1+1)
       N3=2*(K2-K1+1)
       !NBAS=MAXVAL(N1,N2,N3)
       !NBIS=MINVAL(N1,N2,N3)

       MARCA = 1
       DO WHILE (MARCA.EQ.1) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MARCA=0
        IF (N1.LE.NAMRX-2.AND.I1.GT.BOR+1) THEN
         IF (COUNT(CONTA1(I1-1,J1:J2,K1:K2).GE.REFINE_THR).GT.0) THEN
          I1=I1-1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N1.LE.NAMRX-2.AND.I2.LT.NX-BOR) THEN
         !IF (IPATCH.EQ.153) WRITE(*,*) IX,JY,KZ,I1,I2,J1,J2,K1,K2
         IF (COUNT(CONTA1(I2+1,J1:J2,K1:K2).GE.REFINE_THR).GT.0) THEN
          I2=I2+1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J1.GT.BOR+1) THEN
         IF (COUNT(CONTA1(I1:I2,J1-1,K1:K2).GE.REFINE_THR).GT.0) THEN
          J1=J1-1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J2.LT.NY-BOR) THEN
         IF (COUNT(CONTA1(I1:I2,J2+1,K1:K2).GE.REFINE_THR).GT.0) THEN
          J2=J2+1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K1.GT.BOR+1) THEN
         IF (COUNT(CONTA1(I1:I2,J1:J2,K1-1).GE.REFINE_THR).GT.0) THEN
          K1=K1-1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K2.LT.NZ-BOR) THEN
         IF (COUNT(CONTA1(I1:I2,J1:J2,K2+1).GE.REFINE_THR).GT.0) THEN
          K2=K2+1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

       END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       NBIS=MIN(N1,N2,N3)
       IF (NBIS.LE.MIN_PATCHSIZE) THEN
        CR0(IX,JY,KZ)=0
       ELSE
        IPATCH=IPATCH+1
*       WRITE(*,*) IPATCH,N1,N2,N3,
*     &             COUNT(CONTA1(I1:I2,J1:J2,K1:K2).GE.REFINE_THR)
        CONTA1(I1:I2,J1:J2,K1:K2)=0
        CR0(I1:I2,J1:J2,K1:K2)=-1

        PATCHNX(IPATCH)=N1
        PATCHNY(IPATCH)=N2
        PATCHNZ(IPATCH)=N3

        PATCHX(IPATCH)=I1
        PATCHY(IPATCH)=J1
        PATCHZ(IPATCH)=K1

        PATCHRX(IPATCH)=RADX(I1)
        PATCHRY(IPATCH)=RADY(J1)
        PATCHRZ(IPATCH)=RADZ(K1)

        PARE(IPATCH)=0
       END IF

       REFINE_COUNT=COUNT(CR0.GE.REFINE_THR)
       !WRITE(*,*) REFINE_COUNT


      END DO  !-----------------------------------------------------

      NPATCH(IR)=IPATCH

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CR0,CR0AMR),
!$OMP+            PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       IF (CR0(IX,JY,KZ).EQ.-1) THEN
        CR0AMR(IX,JY,KZ)=0
       ELSE
        CR0AMR(IX,JY,KZ)=1
       END IF
      END DO
      END DO
      END DO

      !DEALLOCATE(U1)
      DEALLOCATE(CONTA1)
      DEALLOCATE(CR0)

      WRITE(*,*) 'l=1 patches, l=0 cells refined', NPATCH(IR),
     &           COUNT(CR0AMR.EQ.0)
*     END FIRST LEVEL OF REFINEMENT ====================================

*     START SUBSEQUENT LEVELS OF REFINEMENT ============================
      pare_levels: DO IRPA=1,NL_MESH-1 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       IF (NPATCH(IRPA).EQ.0) THEN
         WRITE(*,*) 'Mesh building stops at level: ', IRPA
         WRITE(*,*) 'There are no more candidate patches'
         EXIT pare_levels
       END IF

       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       !WRITE(*,*) IRPA, LOW1,LOW2

       ALLOCATE(CR01(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))
       ALLOCATE(CONTA11(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))

! PARALLELIZE!!!!!!!!!!!!!!!!!!
       DO IPATCH=LOW1,LOW2 != = = = = = = = = = = = = = = = = = = = = =
        !WRITE(*,*) IPATCH, LOW2
        XL=PATCHRX(IPATCH)-DXPA
        YL=PATCHRY(IPATCH)-DYPA
        ZL=PATCHRZ(IPATCH)-DZPA

        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        DO I=1,SUM(NPART)
         IX=INT((RXPA(I)-XL)/DXPA)+1
         JY=INT((RYPA(I)-YL)/DYPA)+1
         KZ=INT((RZPA(I)-ZL)/DZPA)+1
         IF (IX.GE.1.AND.IX.LE.N1.AND.
     &       JY.GE.1.AND.JY.LE.N2.AND.
     &       KZ.GE.1.AND.KZ.LE.N3) THEN !*****************************
          !U1(IX,JY,KZ)=U1(IX,JY,KZ)+MASAP(I)
          IF (PLEV(I).LE.IRPA) THEN
           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+1
          ELSE
           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+REFINE_THR
          END IF
         END IF !*****************************************************
        END DO

        DO IX=1,N1
        DO JY=1,N2
        DO KZ=1,N3
         IF(IX.LE.BORAMR.OR.IX.GE.N1-BORAMR+1.OR.
     &      JY.LE.BORAMR.OR.JY.GE.N2-BORAMR+1.OR.
     &      KZ.LE.BORAMR.OR.KZ.GE.N3-BORAMR+1) THEN
           CONTA11(IX,JY,KZ,IPATCH)=0
         END IF
         CR01(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)
        END DO
        END DO
        END DO

       END DO != = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       WRITE(*,*) 'Max particles at a cell at l=',IRPA,
     &            MAXVAL(CR01(:,:,:,LOW1:LOW2))

       !CALL VEINSGRID_REDUCED

       DEALLOCATE(CR01)
       DEALLOCATE(CONTA11)

      END DO pare_levels !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
*     END SUBSEQUENT LEVELS OF REFINEMENT
      STOP

      RETURN
      END
