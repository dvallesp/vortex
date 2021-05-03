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
      INTEGER INI_EXTENSION,NBIS,IRPA,BORAMR,LOW1,LOW2,IPATCH,IPARE
      INTEGER INMAX(3),INMAX2(2),I1,I2,J1,J2,K1,K2,N1,N2,N3,IR,MARCA
      INTEGER NP1,NP2,NP3,BASINT,NPALEV3

      INTEGER,ALLOCATABLE::LNPATCH(:)
      INTEGER,ALLOCATABLE::LPATCHNX(:,:),LPATCHNY(:,:),LPATCHNZ(:,:)
      INTEGER,ALLOCATABLE::LPATCHX(:,:),LPATCHY(:,:),LPATCHZ(:,:)
      REAL,ALLOCATABLE::LPATCHRX(:,:),LPATCHRY(:,:),LPATCHRZ(:,:)
      INTEGER,ALLOCATABLE::LVAL(:,:)

!     hard-coded parameters (for now, at least)
      REFINE_THR=3
      BOR=8
      BORAMR=3
      INI_EXTENSION=2 !initial extension of a patch around a cell (on each direction)
      MIN_PATCHSIZE=16 !minimum size (child cells) to be accepted
      NPALEV3=(INT(NAMRX/5)**3)+1
      write(*,*) 'NPALEV3=',NPALEV3

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

!$OMP PARALLEL DO SHARED(CONTA1,CR0,NX,NY,NZ),PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       CONTA1(IX,JY,KZ)=0
       CR0(IX,JY,KZ)=0
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

       I1=MAX(IX-INI_EXTENSION,BOR+1)
       I2=MIN(IX+INI_EXTENSION,NX-BOR)
       J1=MAX(JY-INI_EXTENSION,BOR+1)
       J2=MIN(JY+INI_EXTENSION,NY-BOR)
       K1=MAX(KZ-INI_EXTENSION,BOR+1)
       K2=MIN(KZ+INI_EXTENSION,NZ-BOR)

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
        CR0(I1:I2,J1:J2,K1:K2)=0
        CONTA1(I1:I2,J1:J2,K1:K2)=0
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

      WRITE(*,*) 'At l=',1,', patches:', NPATCH(IR)
      WRITE(*,*) '  --> l=',0,' cells refined:', COUNT(CR0AMR.EQ.0)
*     END FIRST LEVEL OF REFINEMENT ====================================

*     START SUBSEQUENT LEVELS OF REFINEMENT ============================
      pare_levels: DO IRPA=1,NL_MESH-1 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       IF (NPATCH(IRPA).EQ.0) THEN
         WRITE(*,*) 'Mesh building stops at level: ', IRPA
         WRITE(*,*) 'There are no more candidate patches'
         EXIT pare_levels
       END IF

       DXPA=DX/(2.0**IRPA)
       DYPA=DY/(2.0**IRPA)
       DZPA=DZ/(2.0**IRPA)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       !WRITE(*,*) IRPA, LOW1,LOW2

       ALLOCATE(CR01(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))
       ALLOCATE(CONTA11(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))

!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR01,CONTA11),
!$OMP+            PRIVATE(IPATCH,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        DO IX=1,NAMRX
        DO JY=1,NAMRY
        DO KZ=1,NAMRZ
         CR01(IX,JY,KZ,IPATCH)=0
         CONTA11(IX,JY,KZ,IPATCH)=0
        END DO
        END DO
        END DO
       END DO

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHRX,PATCHRY,PATCHRZ,DXPA,DYPA,
!$OMP+                   DZPA,PATCHNX,PATCHNY,PATCHNZ,NPART,RXPA,RYPA,
!$OMP+                   RZPA,CONTA11,REFINE_THR,IRPA,PLEV,BORAMR,CR01),
!$OMP+            PRIVATE(IPATCH,XL,YL,ZL,N1,N2,N3,I,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
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

       WRITE(*,*) '  --> Max particles at a cell:',
     &            MAXVAL(CR01(:,:,:,LOW1:LOW2))

c       REFINE_COUNT=COUNT(CR01(:,:,:,LOW1:LOW2).GE.REFINE_THR)
c       WRITE(*,*) 'Refinable cells BEFORE cleaning at l=',IRPA,
c     &            REFINE_COUNT

       CALL VEINSGRID_REDUCED(IRPA,NPATCH,PARE,PATCHNX,PATCHNY,
     &      PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,CR01,
     &      CONTA11,LOW1,LOW2)

       REFINE_COUNT=COUNT(CR01(:,:,:,LOW1:LOW2).GE.REFINE_THR)
       WRITE(*,*) '  --> Refinable cells AFTER cleaning:',REFINE_COUNT


       ! mesh creation at the next level
       IR=IRPA+1
       !DXPA=DX/(2.0**IR)
       !DYPA=DY/(2.0**IR)
       !DZPA=DZ/(2.0**IR)

       ALLOCATE(LNPATCH(LOW1:LOW2))
       ALLOCATE(LPATCHNX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHNY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHNZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LVAL(NPALEV3,LOW1:LOW2))

       LNPATCH(:)=0
       LVAL(:,:)=0

c       WRITE(*,*) 'REFINABLE CELLS:', REFINE_COUNT

!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR01,REFINE_THR,NPALEV3,PATCHNX,
!$OMP+                   PATCHNY,PATCHNZ,INI_EXTENSION,BORAMR,CONTA11,
!$OMP+                   MIN_PATCHSIZE,DXPA,DYPA,DZPA,LPATCHNX,
!$OMP+                   LPATCHNY,LPATCHNZ,LPATCHX,LPATCHY,LPATCHZ,
!$OMP+                   LPATCHRX,LPATCHRY,LPATCHRZ,LVAL,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ),
!$OMP+            PRIVATE(IPARE,REFINE_COUNT,IPATCH,INMAX,IX,JY,KZ,
!$OMP+                    BASINT,NP1,NP2,NP3,I1,I2,J1,J2,K1,K2,N1,N2,N3,
!$OMP+                    MARCA,NBIS),
!$OMP+            DEFAULT(NONE)
       DO IPARE=LOW1,LOW2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REFINE_COUNT=COUNT(CR01(:,:,:,IPARE).GE.REFINE_THR)
        IPATCH=0
        DO WHILE (REFINE_COUNT.GT.0.AND.IPATCH.LT.NPALEV3) !------------
         INMAX=MAXLOC(CR01(:,:,:,IPARE))
         IX=INMAX(1)
         JY=INMAX(2)
         KZ=INMAX(3)
         BASINT=CR01(IX,JY,KZ,IPARE)
         !IF (CONTA1(IX,JY,KZ).LT.REFINE_THR) EXIT

         NP1=PATCHNX(IPARE)
         NP2=PATCHNY(IPARE)
         NP3=PATCHNZ(IPARE)

         I1=MAX(IX-INI_EXTENSION,BORAMR+1)
         I2=MIN(IX+INI_EXTENSION,NP1-BORAMR)
         J1=MAX(JY-INI_EXTENSION,BORAMR+1)
         J2=MIN(JY+INI_EXTENSION,NP2-BORAMR)
         K1=MAX(KZ-INI_EXTENSION,BORAMR+1)
         K2=MIN(KZ+INI_EXTENSION,NP3-BORAMR)

         N1=2*(I2-I1+1)
         N2=2*(J2-J1+1)
         N3=2*(K2-K1+1)
         !NBAS=MAXVAL(N1,N2,N3)
         !NBIS=MINVAL(N1,N2,N3)

         MARCA = 1
         DO WHILE (MARCA.EQ.1) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          MARCA=0
          IF (N1.LE.NAMRX-2.AND.I1.GT.BORAMR+1) THEN
           IF (COUNT(CONTA11(I1-1,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            I1=I1-1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N1.LE.NAMRX-2.AND.I2.LT.NP1-BORAMR) THEN
           IF (COUNT(CONTA11(I2+1,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            I2=I2+1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J1.GT.BORAMR+1) THEN
           IF (COUNT(CONTA11(I1:I2,J1-1,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            J1=J1-1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J2.LT.NP2-BORAMR) THEN
           IF (COUNT(CONTA11(I1:I2,J2+1,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            J2=J2+1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K1.GT.BORAMR+1) THEN
           IF (COUNT(CONTA11(I1:I2,J1:J2,K1-1,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            K1=K1-1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K2.LT.NP3-BORAMR) THEN
           IF (COUNT(CONTA11(I1:I2,J1:J2,K2+1,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            K2=K2+1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

         END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         NBIS=MIN(N1,N2,N3)
         IF (NBIS.LE.MIN_PATCHSIZE) THEN
          CR01(I1:I2,J1:J2,K1:K2,IPARE)=0
          CONTA11(I1:I2,J1:J2,K1:K2,IPARE)=0
         ELSE
          IPATCH=IPATCH+1
c          WRITE(*,*) 'new,pare:',IPATCH,IPARE
c          WRITE(*,*) 'N1,N2,N3,refinable:',N1,N2,N3,
c     &             COUNT(CONTA11(I1:I2,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
c          write(*,*) 'x,y,z',i1,j1,k1

          CONTA11(I1:I2,J1:J2,K1:K2,IPARE)=0
          CR01(I1:I2,J1:J2,K1:K2,IPARE)=-1

          LPATCHNX(IPATCH,IPARE)=N1
          LPATCHNY(IPATCH,IPARE)=N2
          LPATCHNZ(IPATCH,IPARE)=N3

          LPATCHX(IPATCH,IPARE)=I1
          LPATCHY(IPATCH,IPARE)=J1
          LPATCHZ(IPATCH,IPARE)=K1

          ! remember that dxpa is the cellsize of the parent!!!
          LPATCHRX(IPATCH,IPARE)=PATCHRX(IPARE)+(I1-1.5)*DXPA
          LPATCHRY(IPATCH,IPARE)=PATCHRY(IPARE)+(J1-1.5)*DYPA
          LPATCHRZ(IPATCH,IPARE)=PATCHRZ(IPARE)+(K1-1.5)*DZPA

          LVAL(IPATCH,IPARE)=BASINT
         END IF

         REFINE_COUNT=COUNT(CR01(:,:,:,IPARE).GE.REFINE_THR)
         !WRITE(*,*) REFINE_COUNT
        END DO  !-------------------------------------------------------
       END DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       IPATCH=LOW2
       DO WHILE(COUNT(LVAL(:,:).GT.0).GT.0.AND.IPATCH.LT.NPALEV)
        INMAX2=MAXLOC(LVAL)
        I=INMAX2(1)
        IPARE=LOW1-1+INMAX2(2)
C        WRITE(*,*) LVAL(I,IPARE)
        LVAL(I,IPARE)=0

        IPATCH=IPATCH+1
        IF (IPATCH.GT.NPALEV) EXIT

        PATCHNX(IPATCH)=LPATCHNX(I,IPARE)
        PATCHNY(IPATCH)=LPATCHNY(I,IPARE)
        PATCHNZ(IPATCH)=LPATCHNZ(I,IPARE)

        PATCHX(IPATCH)=LPATCHX(I,IPARE)
        PATCHY(IPATCH)=LPATCHY(I,IPARE)
        PATCHZ(IPATCH)=LPATCHZ(I,IPARE)

        PATCHRX(IPATCH)=LPATCHRX(I,IPARE)
        PATCHRY(IPATCH)=LPATCHRY(I,IPARE)
        PATCHRZ(IPATCH)=LPATCHRZ(I,IPARE)

        PARE(IPATCH)=IPARE
       END DO

       NPATCH(IR)=IPATCH-SUM(NPATCH(0:IR-1))
       IF (SUM(NPATCH).GE.NPALEV) STOP 'NPALEV too small'

       DEALLOCATE(LNPATCH,LPATCHNX,LPATCHNY,LPATCHNZ,LPATCHRX,LPATCHRY,
     &            LPATCHRZ,LPATCHX,LPATCHY,LPATCHZ,LVAL)

       ! still needing to compute cr0amr1(:,:,:,low1:low2)
       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR0AMR1),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        DO IX=1,NAMRX
        DO JY=1,NAMRY
        DO KZ=1,NAMRZ
         CR0AMR1(IX,JY,KZ,IPATCH)=1
        END DO
        END DO
        END DO
       END DO

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PATCHX,
!$OMP+                   PATCHY,PATCHZ,PARE,CR0AMR1),
!$OMP+            PRIVATE(IPATCH,IPARE,I1,I2,J1,J2,K1,K2,N1,N2,N3,
!$OMP+                    IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        IPARE=PARE(IPATCH)
        I1=PATCHX(IPATCH)
        J1=PATCHY(IPATCH)
        K1=PATCHZ(IPATCH)
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        I2=I1+N1/2
        J2=J1+N2/2
        K2=K1+N3/2
        DO IX=I1,I2
        DO JY=J1,J2
        DO KZ=K1,K2
         CR0AMR1(IX,JY,KZ,IPARE)=0
        END DO
        END DO
        END DO
       END DO

       CALL VEINSGRID_CR0AMR(IRPA,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                       PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                       PATCHRZ)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       WRITE(*,*) 'At l=',IR,', patches:', NPATCH(IR)
       WRITE(*,*) '  --> l=',IRPA,' cells refined:',
     &            COUNT(CR0AMR1(:,:,:,LOW1:LOW2).EQ.0)

       DEALLOCATE(CR01, CONTA11)

      END DO pare_levels !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
*     END SUBSEQUENT LEVELS OF REFINEMENT

      IR=IR-1
      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,N1,N2,N3,IX,JY,KZ,CR0AMR1),
!$OMP+            PRIVATE(IPATCH,PATCHNX,PATCHNY,PATCHNZ)
      DO IPATCH=LOW1,LOW2
       N1=PATCHNX(IPATCH)
       N2=PATCHNY(IPATCH)
       N3=PATCHNZ(IPATCH)
       DO IX=1,N1
       DO JY=1,N2
       DO KZ=1,N3
        CR0AMR1(IX,JY,KZ,IPATCH)=1
       END DO
       END DO
       END DO
      END DO

      RETURN
      END



************************************************************************
      SUBROUTINE INTERPOLATE_VELOCITIES(NX,NY,NZ,NL,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,
     &            U4DM,MASAP,NPART,LADO0,BUF,BUFAMR)
************************************************************************
*     Compute the velocity field on the grid
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL*4 RXPA(NDM),RYPA(NDM),RZPA(NDM),
     &       U2DM(NDM),U3DM(NDM),U4DM(NDM),MASAP(NDM)
      REAL LADO0
      INTEGER BUF,BUFAMR

*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      REAL RX(-2:NAMRX+3,NPALEV)
      REAL RY(-2:NAMRX+3,NPALEV)
      REAL RZ(-2:NAMRX+3,NPALEV)
      REAL RMX(-2:NAMRX+3,NPALEV)
      REAL RMY(-2:NAMRX+3,NPALEV)
      REAL RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      INTEGER cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      INTEGER cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      COMMON /cr0/ cr0amr, cr0amr1
      INTEGER solap(NAMRX,NAMRY,NAMRZ,NPALEV)

      REAL U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

      real u1(1:NMAX,1:NMAY,1:NMAZ)
      real u11(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /dens/ u1,u11

      REAL OMEGA0,ACHE,FDM,RHOB0
      COMMON /COSMO/ OMEGA0,ACHE,FDM

      REAL L0(NMAX,NMAY,NMAZ)
      REAL L1(NAMRX,NAMRY,NAMRZ,NPALEV)

      INTEGER NPART0(-BUF+1:NMAX+BUF,-BUF+1:NMAY+BUF,-BUF+1:NMAZ+BUF)
      INTEGER NPART1(-BUFAMR+1:NAMRX+BUFAMR,-BUFAMR+1:NAMRY+BUFAMR,
     &               -BUFAMR+1:NAMRZ+BUFAMR,NPALEV)

      INTEGER IX,JY,KZ,IR,I,IPATCH,LOW1,LOW2,MARCA,CONTA,KPARTICLES
      INTEGER N1,N2,N3,RINT,MINI,MINJ,MINK,MAXI,MAXJ,MAXK,II,JJ,KK
      INTEGER INSI,INSJ,INSK,BASINT,CONTA2,CONTA_PA,I1,I2,J1,J2,K1,K2
      INTEGER CONTA3,JPATCH,J,K,SCR_PAR(NMAX),SCR_IDX(NMAX),IXPAR
      REAL DXPA,DYPA,DZPA,BASX,BASY,BASZ,BAS,RBAS,BASXX,BASYY,BASZZ
      REAL XL,YL,ZL,XR,YR,ZR,MEDIOLADO0,WBAS,CONSTA_DENS,PI
      REAL,ALLOCATABLE::DIST(:),MINS(:),U2INS(:),U3INS(:),U4INS(:)
      REAL,ALLOCATABLE::RXPA_PA(:),RYPA_PA(:),RZPA_PA(:),U2DM_PA(:),
     &                  U3DM_PA(:),U4DM_PA(:),MASAP_PA(:)
      INTEGER,ALLOCATABLE::INDICES_PA(:),IXPA_PA(:),JYPA_PA(:),
     &                     KZPA_PA(:)
      REAL FUIN,U(2,2,2),UW(2,2,2)

      !INTEGER,ALLOCATABLE::SCRINT(:,:,:)
      real t1,t2
      INTEGER IXPA(NDM),JYPA(NDM),KZPA(NDM)

      KPARTICLES=32 !NUMBER OF PARTICLES INSIDE THE KERNEL

      MEDIOLADO0=0.5*LADO0
      PI=ACOS(-1.0)
      RHOB0=MAXVAL(MASAP)/FDM/DX**3
      ! to get overdensity from mass in a sphere (multiply by mass in
      ! code units, divide by radius squared)
      CONSTA_DENS=1/(4*PI/3)/RHOB0
      WRITE(*,*) RHOB0,CONSTA_DENS

      CALL VEINSGRID_ALL_L(NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                     PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &                     SOLAP)


*    1. Find number of particles per l=0 cell, and lengths including at
*       least KPARTICLES particles

!$OMP PARALLEL DO SHARED(NX,NY,NZ,NPART0,BUF),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IX=-BUF+1,NX+BUF
      DO JY=-BUF+1,NY+BUF
      DO KZ=-BUF+1,NZ+BUF
       NPART0(IX,JY,KZ)=0
      END DO
      END DO
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,L0),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       L0(IX,JY,KZ)=0.0
      END DO
      END DO
      END DO

      XL=-LADO0/2
      YL=-LADO0/2
      ZL=-LADO0/2
!$OMP PARALLEL DO SHARED(NPART,NL,RXPA,RYPA,RZPA,XL,YL,ZL,DX,DY,DZ,
!$OMP+                   NX,NY,NZ,IXPA,JYPA,KZPA),
!$OMP+            PRIVATE(I,IX,JY,KZ),
!$OMP+            DEFAULT(NONE),
!$OMP+            REDUCTION(+:NPART0)
      DO I=1,SUM(NPART(0:NL))
       IX=INT((RXPA(I)-XL)/DX)+1
       JY=INT((RYPA(I)-YL)/DY)+1
       KZ=INT((RZPA(I)-ZL)/DZ)+1
       IF (IX.LT.1) IX=1
       IF (IX.GT.NX) IX=NX
       IF (JY.LT.1) JY=1
       IF (JY.GT.NY) JY=NY
       IF (KZ.LT.1) KZ=1
       IF (KZ.GT.NZ) KZ=NZ
       NPART0(IX,JY,KZ)=NPART0(IX,JY,KZ)+1
       IXPA(I)=IX
       JYPA(I)=JY
       KZPA(I)=KZ
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,NPART0,BUF),
!$OMP+            PRIVATE(IX,JY,KZ,II,JJ,KK)
      DO IX=-BUF+1,NX+BUF
      DO JY=-BUF+1,NY+BUF
      DO KZ=-BUF+1,NZ+BUF
       IF (IX.LT.1.OR.IX.GT.NX.OR.
     &     JY.LT.1.OR.JY.GT.NY.OR.
     &     KZ.LT.1.OR.KZ.GT.NZ) THEN
        II=IX
        JJ=JY
        KK=KZ
        IF (II.LT.1) II=II+NX
        IF (II.GT.NX) II=II-NX
        IF (JJ.LT.1) JJ=JJ+NY
        IF (JJ.GT.NY) JJ=JJ-NY
        IF (KK.LT.1) KK=KK+NZ
        IF (KK.GT.NZ) KK=KK-NZ
        NPART0(IX,JY,KZ)=NPART0(II,JJ,KK)
       END IF
      END DO
      END DO
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CR0AMR,KPARTICLES,NPART0,L0),
!$OMP+            PRIVATE(IX,JY,KZ,MARCA,RINT,MINI,MAXI,MINJ,MAXJ,MINK,
!$OMP+                    MAXK),
!$OMP+            DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       IF (CR0AMR(IX,JY,KZ).EQ.1) THEN
        MARCA=0
        RINT=0
        IF (NPART0(IX,JY,KZ).GE.KPARTICLES) MARCA=1
        DO WHILE (MARCA.EQ.0)
         RINT=RINT+1
         MINI=IX-RINT
         MAXI=IX+RINT
         MINJ=JY-RINT
         MAXJ=JY+RINT
         MINK=KZ-RINT
         MAXK=KZ+RINT
         IF (SUM(NPART0(MINI:MAXI,MINJ:MAXJ,MINK:MAXK)).GE.
     &       KPARTICLES) MARCA=1
        END DO
        L0(IX,JY,KZ)=FLOAT(RINT)
       END IF
      END DO
      END DO
      END DO
      WRITE(*,*) 0, MAXVAL(NPART0), MAXVAL(L0)

!$OMP PARALLEL DO SHARED(NPATCH,NL,NPART1,BUFAMR),
!$OMP+            PRIVATE(IPATCH,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IPATCH=1,SUM(NPATCH(0:NL))
       DO IX=-BUFAMR+1,NAMRX+BUFAMR
       DO JY=-BUFAMR+1,NAMRY+BUFAMR
       DO KZ=-BUFAMR+1,NAMRZ+BUFAMR
        NPART1(IX,JY,KZ,IPATCH)=0
       END DO
       END DO
       END DO
      END DO

!$OMP PARALLEL DO SHARED(NPATCH,NL,L1),
!$OMP+            PRIVATE(IPATCH,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IPATCH=1,SUM(NPATCH(0:NL))
       DO IX=1,NAMRX
       DO JY=1,NAMRY
       DO KZ=1,NAMRZ
        L1(IX,JY,KZ,IPATCH)=0.0
       END DO
       END DO
       END DO
      END DO

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ,DXPA,DYPA,DZPA,NPART,NL,RXPA,
!$OMP+                   RYPA,RZPA,NPART1,BUFAMR,CR0AMR1,SOLAP,
!$OMP+                   KPARTICLES,L1,IR),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,XL,YL,ZL,IX,JY,KZ,XR,YR,ZR,
!$OMP+                    MARCA,RINT,MINI,MINJ,MINK,MAXI,MAXJ,MAXK),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        XL=PATCHRX(IPATCH)-(BUFAMR+1)*DXPA
        YL=PATCHRY(IPATCH)-(BUFAMR+1)*DYPA
        ZL=PATCHRZ(IPATCH)-(BUFAMR+1)*DZPA
        XR=XL+(N1+2*BUFAMR)*DXPA
        YR=XL+(N2+2*BUFAMR)*DYPA
        ZR=XL+(N3+2*BUFAMR)*DZPA
        DO I=1,SUM(NPART(0:NL))
         IX=INT((RXPA(I)-XL)/DXPA)+1-BUFAMR
         JY=INT((RYPA(I)-YL)/DYPA)+1-BUFAMR
         KZ=INT((RZPA(I)-ZL)/DZPA)+1-BUFAMR
         IF (IX.GT.-BUFAMR.AND.IX.LE.N1+BUFAMR.AND.
     &       JY.GT.-BUFAMR.AND.JY.LE.N2+BUFAMR.AND.
     &       KZ.GT.-BUFAMR.AND.KZ.LE.N3+BUFAMR) THEN
          NPART1(IX,JY,KZ,IPATCH)=NPART1(IX,JY,KZ,IPATCH)+1
         END IF
        END DO

        DO IX=1,N1
        DO JY=1,N2
        DO KZ=1,N3
         IF (CR0AMR1(IX,JY,KZ,IPATCH).EQ.1.AND.
     &       SOLAP(IX,JY,KZ,IPATCH).EQ.1) THEN
          MARCA=0
          RINT=0
          IF (NPART1(IX,JY,KZ,IPATCH).GE.KPARTICLES) MARCA=1
          DO WHILE (MARCA.EQ.0)
           RINT=RINT+1
           MINI=IX-RINT
           MAXI=IX+RINT
           MINJ=JY-RINT
           MAXJ=JY+RINT
           MINK=KZ-RINT
           MAXK=KZ+RINT
           IF (SUM(NPART1(MINI:MAXI,MINJ:MAXJ,MINK:MAXK,IPATCH)).GE.
     &         KPARTICLES) MARCA=1

           IF (MIN(MINI,MINJ,MINK).EQ.-BUFAMR+1.OR.
     &         MAX(MAXI-N1,MAXJ-N2,MAXK-N3).EQ.BUFAMR) THEN
            IF (SUM(NPART1(MINI:MAXI,MINJ:MAXJ,MINK:MAXK,IPATCH)).GE.
     &          KPARTICLES/4) THEN
             MARCA=1
            ELSE
             WRITE(*,*) IX,JY,KZ,IPATCH,'.',RINT,'.',
     &               SUM(NPART1(MINI:MAXI,MINJ:MAXJ,MINK:MAXK,IPATCH))
            END IF
           END IF
          END DO
          L1(IX,JY,KZ,IPATCH)=FLOAT(RINT)
         END IF
        END DO
        END DO
        END DO

       END DO
       WRITE(*,*) IR, MAXVAL(NPART1(:,:,:,LOW1:LOW2)),
     &            MAXVAL(L1(:,:,:,LOW1:LOW2))
      END DO

      DO IX=1,NX
       SCR_PAR(IX)=COUNT(CR0AMR(IX,:,:).EQ.1)
      END DO
      CALL INDEXX(NX,SCR_PAR,SCR_IDX)

*     3. Now, finally, interpolate the velocity field
!$OMP PARALLEL DO SHARED(NX,NY,NZ,CR0AMR,RADX,RADY,RADZ,L0,NPART0,
!$OMP+                   NPART,NL,IXPA,JYPA,KZPA,MEDIOLADO0,LADO0,
!$OMP+                   KPARTICLES,U1,U2,U3,U4,RXPA,RYPA,RZPA,MASAP,
!$OMP+                   U2DM,U3DM,U4DM,DX,CONSTA_DENS,SCR_IDX,scr_par),
!$OMP+            PRIVATE(IXPAR,IX,JY,KZ,BASX,BASY,BASZ,RINT,MINI,MAXI,
!$OMP+                    MINJ,MAXJ,MINK,MAXK,BASINT,CONTA,DIST,MINS,
!$OMP+                    U2INS,U3INS,U4INS,CONTA2,I,II,JJ,KK,
!$OMP+                    BASXX,BASYY,BASZZ,BAS,RBAS,WBAS,CONTA3),
!$OMP+            DEFAULT(NONE),
!$OMP+            SCHEDULE(DYNAMIC)
      DO IXPAR=1,NX !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IX=SCR_IDX(NX+1-IXPAR)
       write(*,*) ixpar,'-->',IX,scr_par(ix)
      DO JY=1,NY
      DO KZ=1,NZ
       IF (CR0AMR(IX,JY,KZ).EQ.1) THEN
        BASX=RADX(IX)
        BASY=RADY(JY)
        BASZ=RADZ(KZ)

        RINT=INT(L0(IX,JY,KZ))
        IF (RINT.EQ.0) RINT=1

        MINI=IX-RINT
        MAXI=IX+RINT
        MINJ=JY-RINT
        MAXJ=JY+RINT
        MINK=KZ-RINT
        MAXK=KZ+RINT
        BASINT=2*RINT
        CONTA=SUM(NPART0(MINI:MAXI,MINJ:MAXJ,MINK:MAXK))

        ALLOCATE(DIST(CONTA),MINS(CONTA),U2INS(CONTA),
     &           U3INS(CONTA),U4INS(CONTA))
        CONTA2=0
        DO I=1,SUM(NPART(0:NL))
         II=MOD(IXPA(I)-MINI,NMAX)
         IF (II.LT.0) II=II+NMAX
         IF (II.LE.BASINT) THEN
          JJ=MOD(JYPA(I)-MINJ,NMAY)
          IF (JJ.LT.0) JJ=JJ+NMAY
          IF (JJ.LE.BASINT) THEN
           KK=MOD(KZPA(I)-MINK,NMAZ)
           IF (KK.LT.0) KK=KK+NMAZ
           IF (KK.LE.BASINT) THEN
            CONTA2=CONTA2+1
            BASXX=ABS(BASX-RXPA(I))
            BASYY=ABS(BASY-RYPA(I))
            BASZZ=ABS(BASZ-RZPA(I))
            IF (BASXX.GT.MEDIOLADO0) BASXX=BASXX-LADO0
            IF (BASYY.GT.MEDIOLADO0) BASYY=BASYY-LADO0
            IF (BASZZ.GT.MEDIOLADO0) BASZZ=BASZZ-LADO0
            DIST(CONTA2)=SQRT(BASXX**2+BASYY**2+BASZZ**2)
            MINS(CONTA2)=MASAP(I)
            U2INS(CONTA2)=U2DM(I)
            U3INS(CONTA2)=U3DM(I)
            U4INS(CONTA2)=U4DM(I)
           END IF
          END IF
         END IF
        END DO

        IF (CONTA2.LE.KPARTICLES) THEN
         BAS=0.5*MAXVAL(DIST)+0.00001
        ELSE
         CALL SELECT(KPARTICLES,DIST,CONTA2,BAS)
         BAS=0.5*MAX(1.5*DX,BAS)+0.00001
        END IF

        L0(IX,JY,KZ)=BAS
        CALL KERNEL_CUBICSPLINE(CONTA,CONTA2,BAS,DIST)

        RBAS=0.0
        BASXX=0.0
        BASX=0.0
        BASY=0.0
        BASZ=0.0
        CONTA3=0
        DO I=1,CONTA2
         !! INTERPOLATE MASS WITHOUT KERNEL
         !RBAS=RBAS+MINS(I) ! THIS IS THE MASS SUM
         !! INTERPOLATE MASS WITH KERNEL
         IF (DIST(I).GT.0.0) THEN
          CONTA3=CONTA3+1
          RBAS=RBAS+DIST(I)
          WBAS=DIST(I)*MINS(I)
          BASXX=BASXX+WBAS ! THIS IS THE MASS-WEIGHTED KERNEL SUM
          BASX=BASX+WBAS*U2INS(I)
          BASY=BASY+WBAS*U3INS(I)
          BASZ=BASZ+WBAS*U4INS(I)
         END IF
        END DO

        BASX=BASX/BASXX
        BASY=BASY/BASXX
        BASZ=BASZ/BASXX
        !! INTERPOLATE MASS WITHOUT KERNEL
        !BASXX=CONSTA_DENS*RBAS/BAS**3
        !! INTERPOLATE MASS WITH KERNEL
        BASXX=BASXX/RBAS*CONTA3
        BASXX=CONSTA_DENS*BASXX/(2.0*BAS)**3
        !WRITE(*,*) IX,JY,KZ,'->',BAS,BASX,BASY,BASZ,BASXX

        U1(IX,JY,KZ)=BASXX
        U2(IX,JY,KZ)=BASX
        U3(IX,JY,KZ)=BASY
        U4(IX,JY,KZ)=BASZ

        DEALLOCATE(DIST,MINS,U2INS,U3INS,U4INS)
       END IF
      END DO
      END DO
      END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      write(*,*) minval(u1(1:nx,1:ny,1:nz)),maxval(u1(1:nx,1:ny,1:nz))
      write(*,*) minval(u2(1:nx,1:ny,1:nz)),maxval(u2(1:nx,1:ny,1:nz))
      write(*,*) minval(u3(1:nx,1:ny,1:nz)),maxval(u3(1:nx,1:ny,1:nz))
      write(*,*) minval(u4(1:nx,1:ny,1:nz)),maxval(u4(1:nx,1:ny,1:nz))

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ,BUFAMR,DXPA,DYPA,DZPA,
!$OMP+                   MEDIOLADO0,DX,DY,DZ,NPART0,NPART,NL,RXPA,
!$OMP+                   RYPA,RZPA,CR0AMR1,SOLAP,RX,RY,RZ,L1,
!$OMP+                   KPARTICLES,U11,U12,U13,U14,U2DM,U3DM,U4DM,
!$OMP+                   MASAP,NPART1,CONSTA_DENS),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,XL,YL,ZL,XR,YR,ZR,I1,I2,J1,J2,
!$OMP+                    K1,K2,CONTA_PA,INDICES_PA,CONTA,I,BASX,BASY,
!$OMP+                    BASZ,RXPA_PA,RYPA_PA,RZPA_PA,U2DM_PA,U3DM_PA,
!$OMP+                    U4DM_PA,MASAP_PA,IXPA_PA,JYPA_PA,KZPA_PA,IX,
!$OMP+                    JY,KZ,RINT,MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
!$OMP+                    BASINT,DIST,MINS,U2INS,U3INS,U4INS,CONTA2,
!$OMP+                    II,JJ,KK,BASXX,BASYY,BASZZ,BAS,RBAS,WBAS,
!$OMP+                    CONTA3),
!$OMP+            DEFAULT(NONE), SCHEDULE(DYNAMIC)
       DO IPATCH=LOW1,LOW2 !--------------------------------------------
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        write(*,*) ipatch,count(solap(1:n1,1:n2,1:n3,ipatch).eq.1.and.
     &                          cr0amr1(1:n1,1:n2,1:n3,ipatch).eq.1)

        !I1=PATCHX(IPATCH)
        !J1=PATCHY(IPATCH)
        !K1=PATCHZ(IPATCH)
        !I2=I1+N1/2
        !J2=J1+N2/2
        !K2=K1+N3/2

        XL=PATCHRX(IPATCH)-(BUFAMR+1)*DXPA
        YL=PATCHRY(IPATCH)-(BUFAMR+1)*DYPA
        ZL=PATCHRZ(IPATCH)-(BUFAMR+1)*DZPA
        XR=XL+(N1+2*BUFAMR)*DXPA
        YR=YL+(N2+2*BUFAMR)*DYPA
        ZR=ZL+(N3+2*BUFAMR)*DZPA

        I1=INT((XL+MEDIOLADO0)/DX)+1
        I2=INT((XR+MEDIOLADO0)/DX)+2
        J1=INT((YL+MEDIOLADO0)/DY)+1
        J2=INT((YR+MEDIOLADO0)/DY)+2
        K1=INT((ZL+MEDIOLADO0)/DZ)+1
        K2=INT((ZR+MEDIOLADO0)/DZ)+2

        CONTA_PA=SUM(NPART0(I1:I2,J1:J2,K1:K2))

        !write(*,*) ipatch,':',xl,yl,zl,i1,j1,k1
        !write(*,*) ipatch,':',xr,yr,zr,i2,j2,k2
        !write(*,*) ipatch,':',conta_pa
        !stop

        ALLOCATE(INDICES_PA(CONTA_PA))

        CONTA=0
        DO I=1,SUM(NPART(0:NL))
         BASX=RXPA(I)
         IF (XL.LT.BASX.AND.BASX.LT.XR) THEN
          BASY=RYPA(I)
          IF (YL.LT.BASY.AND.BASY.LT.YR) THEN
           BASZ=RZPA(I)
           IF (ZL.LT.BASZ.AND.BASZ.LT.ZR) THEN
            CONTA=CONTA+1
            INDICES_PA(CONTA)=I
           END IF
          END IF
         END IF
        END DO

        !write(*,*) 'conta:',conta

        CONTA_PA=CONTA

        ALLOCATE(RXPA_PA(CONTA_PA),RYPA_PA(CONTA_PA),RZPA_PA(CONTA_PA),
     &           U2DM_PA(CONTA_PA),U3DM_PA(CONTA_PA),U4DM_PA(CONTA_PA),
     &           MASAP_PA(CONTA_PA),IXPA_PA(CONTA_PA),JYPA_PA(CONTA_PA),
     &           KZPA_PA(CONTA_PA))

        DO I=1,CONTA_PA
         BASX=RXPA(INDICES_PA(I))
         BASY=RYPA(INDICES_PA(I))
         BASZ=RZPA(INDICES_PA(I))

         RXPA_PA(I)=BASX
         RYPA_PA(I)=BASY
         RZPA_PA(I)=BASZ
         U2DM_PA(I)=U2DM(INDICES_PA(I))
         U3DM_PA(I)=U3DM(INDICES_PA(I))
         U4DM_PA(I)=U4DM(INDICES_PA(I))
         MASAP_PA(I)=MASAP(INDICES_PA(I))

         IXPA_PA(I)=INT((BASX-XL)/DXPA)+1-BUFAMR
         JYPA_PA(I)=INT((BASY-YL)/DYPA)+1-BUFAMR
         KZPA_PA(I)=INT((BASZ-ZL)/DZPA)+1-BUFAMR
        END DO

        DEALLOCATE(INDICES_PA)

C        write(*,*) 'pre-a.',minval(ixpa_pa),maxval(ixpa_pa),
C     & minval(jypa_pa),maxval(jypa_pa),minval(kzpa_pa),maxval(kzpa_pa)
C        write(*,*) 'rx,nx,xl,xr',PATCHRX(IPATCH),n1,xl,xr
C        write(*,*) 'ry,ny,yl,yr',PATCHRY(IPATCH),n2,yl,yr
C        write(*,*) 'rz,nz,zl,zr',PATCHRZ(IPATCH),n3,zl,zr
c        stop

        DO IX=1,N1 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO JY=1,N2
        DO KZ=1,N3
         IF (CR0AMR1(IX,JY,KZ,IPATCH).EQ.1.AND.
     &       SOLAP(IX,JY,KZ,IPATCH).EQ.1) THEN
          BASX=RX(IX,IPATCH)
          BASY=RY(JY,IPATCH)
          BASZ=RZ(KZ,IPATCH)

          RINT=INT(L1(IX,JY,KZ,IPATCH))
          IF (RINT.EQ.0) RINT=1

          MINI=IX-RINT
          MAXI=IX+RINT
          MINJ=JY-RINT
          MAXJ=JY+RINT
          MINK=KZ-RINT
          MAXK=KZ+RINT
          BASINT=2*RINT
          CONTA=SUM(NPART1(MINI:MAXI,MINJ:MAXJ,MINK:MAXK,IPATCH))
          ALLOCATE(DIST(CONTA),MINS(CONTA),U2INS(CONTA),
     &             U3INS(CONTA),U4INS(CONTA))

c          write(*,*) '.a.',basx,basy,basz,L1(ix,jy,kz,ipatch),rint
c          write(*,*) '.b.',mini,maxi,minj,maxj,mink,maxk
c          write(*,*) '.c.',conta

          CONTA2=0
          conta_pa_loop: DO I=1,CONTA_PA
           II=IXPA_PA(I)-MINI
           IF (II.GE.0.AND.II.LE.BASINT) THEN
            JJ=JYPA_PA(I)-MINJ
            IF (JJ.GE.0.AND.JJ.LE.BASINT) THEN
             KK=KZPA_PA(I)-MINK
             IF (KK.GE.0.AND.KK.LE.BASINT) THEN
              CONTA2=CONTA2+1
              IF (CONTA2.GT.CONTA) THEN
               WRITE(*,*) 'WARNING TOO MUCH PARTICLES', IX, JY, KZ,
     &                    CONTA, CONTA2
               EXIT conta_pa_loop
              END IF
              BASXX=ABS(BASX-RXPA_PA(I))
              BASYY=ABS(BASY-RYPA_PA(I))
              BASZZ=ABS(BASZ-RZPA_PA(I))
              DIST(CONTA2)=SQRT(BASXX**2+BASYY**2+BASZZ**2)
              MINS(CONTA2)=MASAP_PA(I)
              U2INS(CONTA2)=U2DM_PA(I)
              U3INS(CONTA2)=U3DM_PA(I)
              U4INS(CONTA2)=U4DM_PA(I)
             END IF
            END IF
           END IF
          END DO conta_pa_loop

c          write(*,*) '.d.',conta2

          IF (CONTA2.LE.KPARTICLES) THEN
           BAS=0.5*MAXVAL(DIST)+0.00001
          ELSE
           CALL SELECT(KPARTICLES,DIST,CONTA2,BAS)
           BAS=0.5*MAX(1.5*DXPA,BAS)+0.00001
          END IF

c          write(*,*) '.e.',bas

          L1(IX,JY,KZ,IPATCH)=BAS
          CALL KERNEL_CUBICSPLINE(CONTA,CONTA2,BAS,DIST)

c          write(*,*) '.f.',minval(dist),maxval(dist)

          RBAS=0.0
          BASXX=0.0
          BASX=0.0
          BASY=0.0
          BASZ=0.0
          CONTA3=0
          DO I=1,CONTA2
           !! INTERPOLATE MASS WITHOUT KERNEL
           !RBAS=RBAS+MINS(I) ! THIS IS THE MASS SUM
           !! INTERPOLATE MASS WITH KERNEL
           IF (DIST(I).GT.0.0) THEN
            CONTA3=CONTA3+1
            RBAS=RBAS+DIST(I)
            WBAS=DIST(I)*MINS(I)
            BASXX=BASXX+WBAS ! THIS IS THE MASS-WEIGHTED KERNEL SUM
            BASX=BASX+WBAS*U2INS(I)
            BASY=BASY+WBAS*U3INS(I)
            BASZ=BASZ+WBAS*U4INS(I)
            !write(*,*) u2ins(i),u3ins(i),u4ins(i)
           END IF
          END DO

C          write(*,*) '.g.',conta3,rbas,basxx,basx,basy,basz

          BASX=BASX/BASXX
          BASY=BASY/BASXX
          BASZ=BASZ/BASXX
          !! INTERPOLATE MASS WITHOUT KERNEL
          !BASXX=CONSTA_DENS*RBAS/BAS**3
          !! INTERPOLATE MASS WITH KERNEL
          BASXX=BASXX/RBAS*CONTA3
          BASXX=CONSTA_DENS*BASXX/(2.0*BAS)**3
          !WRITE(*,*) IX,JY,KZ,IPATCH,'->',BAS,BASX,BASY,BASZ,BASXX

          U11(IX,JY,KZ,IPATCH)=BASXX
          U12(IX,JY,KZ,IPATCH)=BASX
          U13(IX,JY,KZ,IPATCH)=BASY
          U14(IX,JY,KZ,IPATCH)=BASZ

          DEALLOCATE(DIST,MINS,U2INS,U3INS,U4INS)

         END IF
        END DO
        END DO
        END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        DEALLOCATE(RXPA_PA,RYPA_PA,RZPA_PA,U2DM_PA,U3DM_PA,U4DM_PA,
     &             MASAP_PA,IXPA_PA,JYPA_PA,KZPA_PA)

       END DO
c       write(*,*) 'End level', IR
c       write(*,*) 'density min, max',minval(u11(:,:,:,low1:low2)),
c     &                               maxval(u11(:,:,:,low1:low2))
c       write(*,*) 'vx min, max',minval(u12(:,:,:,low1:low2)),
c     &                               maxval(u12(:,:,:,low1:low2))
c       write(*,*) 'vy min, max',minval(u13(:,:,:,low1:low2)),
c     &                               maxval(u13(:,:,:,low1:low2))
c       write(*,*) 'vz min, max',minval(u14(:,:,:,low1:low2)),
c     &                               maxval(u14(:,:,:,low1:low2))
      END DO

*     refill refined and overlapping cells
      DO IR=NL,1,-1
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U11,NL)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    U12(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    U13(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    U14(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO ipatch=LOW1,LOW2
          !WRITE(*,*) 'FINISHING PATCH', IPATCH
          N1 = PATCHNX(IPATCH)
          N2 = PATCHNY(IPATCH)
          N3 = PATCHNZ(IPATCH)
          JPATCH = PARE(IPATCH)
          DO I=1,N1,2
          DO J=1,N2,2
          DO K=1,N3,2
            II = PATCHX(ipatch) + int((I-1)/2)
            JJ = PATCHY(ipatch) + int((J-1)/2)
            KK = PATCHZ(ipatch) + int((K-1)/2)
            if (jpatch.ne.0) then
             uw(1:2,1:2,1:2) = 1.0
             u(1:2,1:2,1:2) = u11(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u11(II,JJ,KK,JPATCH) = FUIN

             uw(1:2,1:2,1:2) = u(1:2,1:2,1:2)
             u(1:2,1:2,1:2) = u12(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u12(II,JJ,KK,JPATCH) = FUIN

             u(1:2,1:2,1:2) = u13(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u13(II,JJ,KK,JPATCH) = FUIN

             u(1:2,1:2,1:2) = u14(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u14(II,JJ,KK,JPATCH) = FUIN
            else
             uw(1:2,1:2,1:2) = 1.0
             u(1:2,1:2,1:2) = u11(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u1(II,JJ,KK) = FUIN

             uw(1:2,1:2,1:2) = u(1:2,1:2,1:2)
             u(1:2,1:2,1:2) = u12(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u2(II,JJ,KK) = FUIN

             u(1:2,1:2,1:2) = u13(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u3(II,JJ,KK) = FUIN

             u(1:2,1:2,1:2) = u14(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u4(II,JJ,KK) = FUIN
            end if
          END DO
          END DO
          END DO
        END DO
      END DO !IR=NL,1,-1

      write(*,*) 'At level', 0
      write(*,*) 'density min, max',minval(u1),maxval(u1)
      write(*,*) 'vx min,max',minval(u2),maxval(u2)
      write(*,*) 'vy min,max',minval(u3),maxval(u3)
      write(*,*) 'vz min,max',minval(u4),maxval(u4)

      DO IR=1,NL
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
       write(*,*) 'density min, max',minval(u11(:,:,:,low1:low2)),
     &                               maxval(u11(:,:,:,low1:low2))
       write(*,*) 'vx min, max',minval(u12(:,:,:,low1:low2)),
     &                          maxval(u12(:,:,:,low1:low2))
       write(*,*) 'vy min, max',minval(u13(:,:,:,low1:low2)),
     &                          maxval(u13(:,:,:,low1:low2))
       write(*,*) 'vz min, max',minval(u14(:,:,:,low1:low2)),
     &                          maxval(u14(:,:,:,low1:low2))
      END DO

      OPEN(99,FILE='output_files/particle-grid',STATUS='UNKNOWN',
     &     FORM='UNFORMATTED')

      write(99) nl
      write(99) (npatch(i),i=0,nl)
      write(99) (patchnx(i),i=1,sum(npatch(0:nl)))
      write(99) (patchny(i),i=1,sum(npatch(0:nl)))
      write(99) (patchnz(i),i=1,sum(npatch(0:nl)))
      write(99) (patchx(i),i=1,sum(npatch(0:nl)))
      write(99) (patchy(i),i=1,sum(npatch(0:nl)))
      write(99) (patchz(i),i=1,sum(npatch(0:nl)))
      write(99) (patchrx(i),i=1,sum(npatch(0:nl)))
      write(99) (patchry(i),i=1,sum(npatch(0:nl)))
      write(99) (patchrz(i),i=1,sum(npatch(0:nl)))
      write(99) (pare(i),i=1,sum(npatch(0:nl)))
      write(99) (((u1(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      write(99) (((u2(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      write(99) (((u3(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      write(99) (((u4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      write(99) (((cr0amr(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
      do ipatch=1,sum(npatch(0:nl))
       n1=patchnx(ipatch)
       n2=patchny(ipatch)
       n3=patchnz(ipatch)
       write(99) (((u11(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
       write(99) (((u12(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
       write(99) (((u13(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
       write(99) (((u14(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
       write(99) (((cr0amr1(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
       write(99) (((solap(I,J,K,ipatch),I=1,n1),J=1,n2),K=1,n3)
      end do

      CLOSE(99)
      stop

      RETURN
      END




************************************************************************
      SUBROUTINE KERNEL_CUBICSPLINE(N,N2,W,DIST)
************************************************************************
*     DIST contains initially the distance (particle to cell), and it is
*     updated with the (unnormalised) value of the kernel
      INTEGER N,N2 ! N is the dimension of the array dist;
                   ! N2, the actual number of particles filled in
      REAL W,DIST(N)

      REAL DISTS
      INTEGER I

      DO I=1,N2
       DISTS=DIST(I)/W
       IF (DISTS.LE.1.0) THEN
        DIST(I)=1.0-1.5*DISTS**2*(1-0.5*DISTS)
       ELSE IF (DISTS.LE.2.0) THEN
        DIST(I)=0.25*(2-DISTS)**3
       ELSE
        DIST(I)=0.0
       END IF
      END DO

      RETURN
      END



************************************************************************
      SUBROUTINE VEINSGRID_REDUCED(IR,NPATCH,PARE,
     &      PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &      PATCHRY,PATCHRZ,CR01,CONTA11,LOW1,LOW2)
************************************************************************
*     small fraction of patches with rare geometry were not detected
*     when overlapping

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NPALEV2

*     U11(PATCHNX,PATCHNY,PATCHNZ,NLEVEL,NPALEV)
*     PATCHNX,PATCHNY,PATCHNZ patches dimensions
*     IPATCH number of patches per level
*     NLEVELS total number of levels

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      REAL*4  PATCHRX(NPALEV)
      REAL*4  PATCHRY(NPALEV)
      REAL*4  PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK,I2,J2
      INTEGER N1,N2,N3,L1,L2,L3,NN1,NN2,NN3,LL1,LL2,LL3
      INTEGER KK2,JJ2,II2,KZ2,JY2,IX2
      INTEGER NV,A2,B2,C2,K

      INTEGER LOW1,LOW2
      INTEGER SOLAP_PATCH(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)
      INTEGER CR01(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)
      INTEGER CONTA11(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)

      REAL*4 A1,B1,C1,RIV1,RIV2,RIV3
      INTEGER CONTROL
      INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
      INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
      INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
      REAL*4 RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
      REAL*4 RZ1,RZZ1,RZ2,RZZ2,ORXX1,ORYY1,ORZZ1

      REAL*4 DXPA,DYPA,DZPA
      REAL*4 DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL*4  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      INTEGER,ALLOCATABLE::VECINO(:,:)
      INTEGER,ALLOCATABLE::NVECI(:)

      INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
      REAL*4 RXFIX,RYFIX,RZFIX

      REAL*4 OVERLAP
*
       NPALEV2=MAX(100,INT(NPALEV/5))

       ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
       ALLOCATE(NVECI(NPATCH(IR)))

       DXPA=DX/(2.**IR)
       DYPA=DY/(2.**IR)
       DZPA=DZ/(2.**IR)

*      built auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

!$OMP   PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI,
!$OMP+        DXPA,DYPA,DZPA,PATCHRX,PATCHRY,PATCHRZ,
!$OMP+        SOLAP_PATCH,LOW1,LOW2,RXFIX,RYFIX,RZFIX),
!$OMP+  PRIVATE(I,N1,N2,N3,NV,J,NN1,NN2,NN3,
!$OMP+          RX1,RY1,RZ1,RXX1,RYY1,RZZ1,I2,
!$OMP+          IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4),
!$OMP+  DEFAULT(NONE)
       DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         SOLAP_PATCH(:,:,:,I)=0

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         NV=0

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA

         IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
         JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
         KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

         IG2=IG1 + N1 - 1
         JG2=JG1 + N2 - 1
         KG2=KG1 + N3 - 1

         DO J=LOW1,LOW2
          IF (J.NE.I) THEN

          NN1=PATCHNX(J)
          NN2=PATCHNY(J)
          NN3=PATCHNZ(J)

          RXX1=PATCHRX(J)-0.5*DXPA
          RYY1=PATCHRY(J)-0.5*DYPA
          RZZ1=PATCHRZ(J)-0.5*DZPA

          IG3=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1

          IG4=IG3 + NN1 - 1
          JG4=JG3 + NN2 - 1
          KG4=KG3 + NN3 - 1

          IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &        JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &        KG1.LE.KG4.AND.KG3.LE.KG2) THEN
           NV=NV+1
           VECINO(NV,I2)=J
          END IF

          END IF

         END DO
         NVECI(I2)=NV
       END DO


       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) THEN
         WRITE(*,*) 'ERROR: gvecino ST second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))
         STOP
       END IF


       DO I=LOW1,LOW2

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA

         I2=I-LOW1+1

         DO K=1,NVECI(I2)
         J=VECINO(K,I2)
         J2=J-LOW1+1

         LL1=PATCHX(J)
         LL2=PATCHY(J)
         LL3=PATCHZ(J)

         NN1=PATCHNX(J)
         NN2=PATCHNY(J)
         NN3=PATCHNZ(J)

         RXX1=PATCHRX(J)-0.5*DXPA
         RYY1=PATCHRY(J)-0.5*DYPA
         RZZ1=PATCHRZ(J)-0.5*DZPA
         RXX2=PATCHRX(J)-0.5*DXPA+(NN1-1)*DXPA
         RYY2=PATCHRY(J)-0.5*DYPA+(NN2-1)*DYPA
         RZZ2=PATCHRZ(J)-0.5*DZPA+(NN3-1)*DZPA

*        X
         IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=1
            CORNXX2=NN1
         END IF
         IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=N1
            CORNXX1=1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) +1
         END IF
         IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
            CORNX1=1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=NN1
         END IF
         IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=1
            CORNX2=N1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) + 1
         END IF

*        Y
         IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=1
            CORNYY2=NN2
         END IF
         IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=N2
            CORNYY1=1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) +1
         END IF
         IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
            CORNY1=1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=NN2
         END IF
         IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=1
            CORNY2=N2
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) + 1
         END IF

*        Z
         IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=N3
            CORNZZ1=1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) +1
         END IF
         IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
            CORNZ1=1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=1
            CORNZ2=N3
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) + 1
         END IF

*        overlap is the fraction of volumen overlapping
*        the idea is to fix a thershold of overlapping that below
*        this therhold for instance 10%, the cells are not marked as overlap
         OVERLAP=(CORNZZ2-CORNZZ1+1)*(CORNYY2-CORNYY1+1)*
     &           (CORNXX2-CORNXX1+1)
         OVERLAP=ABS(OVERLAP)/PATCHNX(J)/PATCHNY(J)/PATCHNZ(J)

**       celdas madre del nivel inferior
         IF (OVERLAP.GT.0.1) THEN
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
            IX=II-CORNXX1+CORNX1
            JY=JJ-CORNYY1+CORNY1
            KZ=KK-CORNZZ1+CORNZ1
            IF (SOLAP_PATCH(IX,JY,KZ,I).EQ.0) THEN
*            the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i
             SOLAP_PATCH(II,JJ,KK,J)=1
             CR01(II,JJ,KK,J)=0
             CONTA11(II,JJ,KK,J)=0
            END IF
           END DO
           END DO
           END DO
          END IF
       END DO
       END DO

      DEALLOCATE(VECINO)
      DEALLOCATE(NVECI)

      RETURN
      END


************************************************************************
      SUBROUTINE VEINSGRID_CR0AMR(IR,NPATCH,PARE,
     &      PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &      PATCHRY,PATCHRZ)
************************************************************************
*     small fraction of patches with rare geometry were not detected
*     when overlapping

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NPALEV2

*     U11(PATCHNX,PATCHNY,PATCHNZ,NLEVEL,NPALEV)
*     PATCHNX,PATCHNY,PATCHNZ patches dimensions
*     IPATCH number of patches per level
*     NLEVELS total number of levels

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      REAL*4  PATCHRX(NPALEV)
      REAL*4  PATCHRY(NPALEV)
      REAL*4  PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK,I2,J2
      INTEGER N1,N2,N3,L1,L2,L3,NN1,NN2,NN3,LL1,LL2,LL3
      INTEGER KK2,JJ2,II2,KZ2,JY2,IX2
      INTEGER NV,A2,B2,C2,K

      INTEGER LOW1,LOW2
      INTEGER SOLAP_PATCH(NAMRX,NAMRY,NAMRZ,NPATCH(IR))

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

      REAL*4 A1,B1,C1,RIV1,RIV2,RIV3
      INTEGER CONTROL
      INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
      INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
      INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
      REAL*4 RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
      REAL*4 RZ1,RZZ1,RZ2,RZZ2,ORXX1,ORYY1,ORZZ1

      REAL*4 DXPA,DYPA,DZPA
      REAL*4 DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL*4  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      INTEGER,ALLOCATABLE::VECINO(:,:)
      INTEGER,ALLOCATABLE::NVECI(:)

      INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
      REAL*4 RXFIX,RYFIX,RZFIX

      REAL*4 OVERLAP
*
       NPALEV2=MAX(100,INT(NPALEV/5))

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

       ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
       ALLOCATE(NVECI(NPATCH(IR)))

       DXPA=DX/(2.**IR)
       DYPA=DY/(2.**IR)
       DZPA=DZ/(2.**IR)

*      built auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

!$OMP   PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI,
!$OMP+        DXPA,DYPA,DZPA,PATCHRX,PATCHRY,PATCHRZ,
!$OMP+        SOLAP_PATCH,LOW1,LOW2,RXFIX,RYFIX,RZFIX),
!$OMP+  PRIVATE(I,N1,N2,N3,NV,J,NN1,NN2,NN3,
!$OMP+          RX1,RY1,RZ1,RXX1,RYY1,RZZ1,I2,
!$OMP+          IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4),
!$OMP+  DEFAULT(NONE)
       DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         SOLAP_PATCH(:,:,:,I2)=0

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         NV=0

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA

         IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
         JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
         KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

         IG2=IG1 + N1 - 1
         JG2=JG1 + N2 - 1
         KG2=KG1 + N3 - 1

         DO J=LOW1,LOW2
          IF (J.NE.I) THEN

          NN1=PATCHNX(J)
          NN2=PATCHNY(J)
          NN3=PATCHNZ(J)

          RXX1=PATCHRX(J)-0.5*DXPA
          RYY1=PATCHRY(J)-0.5*DYPA
          RZZ1=PATCHRZ(J)-0.5*DZPA

          IG3=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1

          IG4=IG3 + NN1 - 1
          JG4=JG3 + NN2 - 1
          KG4=KG3 + NN3 - 1

          IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &        JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &        KG1.LE.KG4.AND.KG3.LE.KG2) THEN
           NV=NV+1
           VECINO(NV,I2)=J
          END IF

          END IF

         END DO
         NVECI(I2)=NV
       END DO


       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) THEN
         WRITE(*,*) 'ERROR: gvecino ST second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))
         STOP
       END IF


       DO I=LOW1,LOW2

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA

         I2=I-LOW1+1

         DO K=1,NVECI(I2)
         J=VECINO(K,I2)
         J2=J-LOW1+1

         LL1=PATCHX(J)
         LL2=PATCHY(J)
         LL3=PATCHZ(J)

         NN1=PATCHNX(J)
         NN2=PATCHNY(J)
         NN3=PATCHNZ(J)

         RXX1=PATCHRX(J)-0.5*DXPA
         RYY1=PATCHRY(J)-0.5*DYPA
         RZZ1=PATCHRZ(J)-0.5*DZPA
         RXX2=PATCHRX(J)-0.5*DXPA+(NN1-1)*DXPA
         RYY2=PATCHRY(J)-0.5*DYPA+(NN2-1)*DYPA
         RZZ2=PATCHRZ(J)-0.5*DZPA+(NN3-1)*DZPA

*        X
         IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=1
            CORNXX2=NN1
         END IF
         IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=N1
            CORNXX1=1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) +1
         END IF
         IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
            CORNX1=1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=NN1
         END IF
         IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=1
            CORNX2=N1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) + 1
         END IF

*        Y
         IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=1
            CORNYY2=NN2
         END IF
         IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=N2
            CORNYY1=1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) +1
         END IF
         IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
            CORNY1=1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=NN2
         END IF
         IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=1
            CORNY2=N2
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) + 1
         END IF

*        Z
         IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=N3
            CORNZZ1=1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) +1
         END IF
         IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
            CORNZ1=1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=1
            CORNZ2=N3
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) + 1
         END IF

**       celdas madre del nivel inferior
          DO KK=CORNZZ1,CORNZZ2
          DO JJ=CORNYY1,CORNYY2
          DO II=CORNXX1,CORNXX2
           IX=II-CORNXX1+CORNX1
           JY=JJ-CORNYY1+CORNY1
           KZ=KK-CORNZZ1+CORNZ1
           IF (SOLAP_PATCH(IX,JY,KZ,I2).EQ.0) THEN
*            the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i
            IF (CR0AMR1(IX,JY,KZ,I).EQ.0) THEN
             CR0AMR1(II,JJ,KK,J)=0
            ELSE IF (CR0AMR1(II,JJ,KK,J).EQ.0) THEN
             CR0AMR1(IX,JY,KZ,I)=0
            END IF
           END IF
          END DO
          END DO
          END DO
       END DO
       END DO

      DEALLOCATE(VECINO)
      DEALLOCATE(NVECI)

      RETURN
      END

***************************************************************************
      SUBROUTINE indexx(n,arr,indx)
***************************************************************************
* From Press,Teukoslky,Vetterling & Flannery, Numerical Recipes in Fortran90,
* Cambridge University Press
***************************************************************************

      INTEGER n,indx(n),M,NSTACK
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
         write(*,*) 'NSTACK too small in indexx'
         stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
