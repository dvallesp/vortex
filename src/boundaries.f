***********************************************************************
      SUBROUTINE CORRECT_SOURCE_BOUNDARIES(NL,NX,NY,NZ,NPATCH,
     &           PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     Function parameters
      INTEGER NL, NX, NY, NZ
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     GLOBAL VARIABLES
*     Grid
      REAL*4 RX(-2:NAMRX+3,NPALEV)
      REAL*4 RY(-2:NAMRX+3,NPALEV)
      REAL*4 RZ(-2:NAMRX+3,NPALEV)
      REAL*4 RMX(-2:NAMRX+3,NPALEV)
      REAL*4 RMY(-2:NAMRX+3,NPALEV)
      REAL*4 RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      REAL*4  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/ RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

*     Differential operators (rotational and divergence)
      REAL*4 ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 ROTAX_1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL*4 ROTAY_1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL*4 ROTAZ_1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      REAL*4 DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 DIVER(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /DIVERGENCE/ DIVER0, DIVER

*     private variables
      INTEGER IR,I,N1,N2,N3,LOW1,LOW2,IX,JY,KZ
      REAL*4 UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),FUIN
      REAL*4 AAA, BBB, CCC
      INTEGER KARE,KR1,KR2,KR3

*     Level 1 patches: interpolate from the base grid
      IR=1

!$OMP PARALLEL DO SHARED(IR,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+         CR3AMR1X,CR3AMR1Y,CR3AMR1Z,
!$OMP+         ROTAX_0, ROTAY_0, ROTAZ_0, DIVER0, ROTAX_1, ROTAY_1,
!$OMP+         ROTAZ_1, DIVER),
!$OMP+      PRIVATE(I,N1,N2,N3,IX,JY,KZ,UBAS,FUIN,KR1,KR2,KR3)
      DO I=1,NPATCH(IR)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1

*       #######################################################
        IF (IX.EQ.1.OR.IX.EQ.N1.OR.JY.EQ.1.OR.JY.EQ.N2.OR.
     &      KZ.EQ.1.OR.KZ.EQ.N3) THEN
*       #######################################################

          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          UBAS(1:3,1:3,1:3)=ROTAX_0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          ROTAX_1(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=ROTAY_0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          ROTAY_1(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=ROTAZ_0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          ROTAZ_1(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=DIVER0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          DIVER(IX,JY,KZ,I)=FUIN

*       #######################################################
        END IF
*       #######################################################

       END DO
       END DO
       END DO

      END DO

*     Other levels (l>1)
      DO IR=2,NL

      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(IR,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+         LOW1,LOW2,RX,RY,RZ,RADX,RADY,RADZ,
!$OMP+         CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,
!$OMP+         ROTAX_0, ROTAY_0, ROTAZ_0, DIVER0, ROTAX_1, ROTAY_1,
!$OMP+         ROTAZ_1, DIVER),
!$OMP+     PRIVATE(I,N1,N2,N3,IX,JY,KZ,KARE,KR1,KR2,KR3,UBAS,RXBAS,
!$OMP+             RYBAS,RZBAS,FUIN,AAA,BBB,CCC)
      DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1

*       #######################################################
        IF (IX.EQ.1.OR.IX.EQ.N1.OR.JY.EQ.1.OR.JY.EQ.N2.OR.
     &      KZ.EQ.1.OR.KZ.EQ.N3) THEN
*       #######################################################

          KARE=CR3AMR1(IX,JY,KZ,I)
          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          AAA=RX(IX,I)
          BBB=RY(JY,I)
          CCC=RZ(KZ,I)

          IF (KARE.GT.0) THEN
            RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
            RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
            RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)
          ELSE
            RXBAS(1:3)=RADX(KR1-1:KR1+1)
            RYBAS(1:3)=RADY(KR2-1:KR2+1)
            RZBAS(1:3)=RADZ(KR3-1:KR3+1)
          ENDIF

*         ROTAX
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        ROTAX_1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        ROTAX_0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          ROTAX_1(IX,JY,KZ,I)=FUIN

*         ROTAY
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        ROTAY_1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        ROTAY_0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          ROTAY_1(IX,JY,KZ,I)=FUIN

*         ROTAZ
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        ROTAZ_1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        ROTAZ_0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          ROTAZ_1(IX,JY,KZ,I)=FUIN

*         DIVER
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        DIVER(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        DIVER0(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          DIVER(IX,JY,KZ,I)=FUIN

*       #######################################################
        END IF
*       #######################################################

       END DO
       END DO
       END DO

      END DO ! I = LOW1, LOW2
      END DO ! IR = 2, NL

      RETURN
      END


***********************************************************************
      SUBROUTINE CORRECT_VELOCITY_BOUNDARIES(NL,NX,NY,NZ,NPATCH,
     &           PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     Function parameters
      INTEGER NL, NX, NY, NZ
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     GLOBAL VARIABLES
*     Grid
      REAL*4 RX(-2:NAMRX+3,NPALEV)
      REAL*4 RY(-2:NAMRX+3,NPALEV)
      REAL*4 RZ(-2:NAMRX+3,NPALEV)
      REAL*4 RMX(-2:NAMRX+3,NPALEV)
      REAL*4 RMY(-2:NAMRX+3,NPALEV)
      REAL*4 RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      REAL*4  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/ RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

*     VELOCITY COMPONENTS
*     compressive velocity
      REAL*4 U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL*4 U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL*4 U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P
*     rotational velocity (ROTS variables were reused)
      REAL*4 U2R(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 U3R(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 U4R(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL*4 U12R(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL*4 U13R(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL*4 U14R(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /ROTS/ U2R,U3R,U4R,U12R,U13R,U14R

*     private variables
      INTEGER IR,I,N1,N2,N3,LOW1,LOW2,IX,JY,KZ
      REAL*4 UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),FUIN
      REAL*4 AAA, BBB, CCC
      INTEGER KARE,KR1,KR2,KR3

*     Level 1 patches: interpolate from the base grid
      IR=1

!$OMP PARALLEL DO SHARED(IR,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+         CR3AMR1X,CR3AMR1Y,CR3AMR1Z,
!$OMP+         U2P, U3P, U4P, U2R, U3R, U4R, U12P, U13P, U14P, U12R,
!$OMP+         U13R, U14R),
!$OMP+      PRIVATE(I,N1,N2,N3,IX,JY,KZ,UBAS,FUIN,KR1,KR2,KR3)
      DO I=1,NPATCH(IR)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1

*       #######################################################
        IF (IX.EQ.1.OR.IX.EQ.N1.OR.JY.EQ.1.OR.JY.EQ.N2.OR.
     &      KZ.EQ.1.OR.KZ.EQ.N3) THEN
*       #######################################################

          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          UBAS(1:3,1:3,1:3)=U2P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U12P(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=U3P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U13P(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=U4P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U14P(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=U2R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U12R(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=U3R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U13R(IX,JY,KZ,I)=FUIN

          UBAS(1:3,1:3,1:3)=U4R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U14R(IX,JY,KZ,I)=FUIN

*       #######################################################
        END IF
*       #######################################################

       END DO
       END DO
       END DO

      END DO

*     Other levels (l>1)
      DO IR=2,NL

      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(IR,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+         LOW1,LOW2,RX,RY,RZ,RADX,RADY,RADZ,
!$OMP+         CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,
!$OMP+         U2P, U3P, U4P, U2R, U3R, U4R, U12P, U13P, U14P, U12R,
!$OMP+         U13R, U14R),
!$OMP+     PRIVATE(I,N1,N2,N3,IX,JY,KZ,KARE,KR1,KR2,KR3,UBAS,RXBAS,
!$OMP+             RYBAS,RZBAS,FUIN,AAA,BBB,CCC)
      DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1

*       #######################################################
        IF (IX.EQ.1.OR.IX.EQ.N1.OR.JY.EQ.1.OR.JY.EQ.N2.OR.
     &      KZ.EQ.1.OR.KZ.EQ.N3) THEN
*       #######################################################

          KARE=CR3AMR1(IX,JY,KZ,I)
          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          AAA=RX(IX,I)
          BBB=RY(JY,I)
          CCC=RZ(KZ,I)

          IF (KARE.GT.0) THEN
            RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
            RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
            RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)
          ELSE
            RXBAS(1:3)=RADX(KR1-1:KR1+1)
            RYBAS(1:3)=RADY(KR2-1:KR2+1)
            RZBAS(1:3)=RADZ(KR3-1:KR3+1)
          ENDIF

*         U2P
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U12P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U2P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          U12P(IX,JY,KZ,I)=FUIN

*         U3P
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U13P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U3P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          U13P(IX,JY,KZ,I)=FUIN

*         U4P
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U14P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U4P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          U14P(IX,JY,KZ,I)=FUIN

*         U2R
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U12R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U2R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          U12R(IX,JY,KZ,I)=FUIN

*         U3R
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U13R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U3R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          U13R(IX,JY,KZ,I)=FUIN

*         U4R
          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U14R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U4R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          ENDIF

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
          U14R(IX,JY,KZ,I)=FUIN

*       #######################################################
        END IF
*       #######################################################

       END DO
       END DO
       END DO

      END DO ! I = LOW1, LOW2
      END DO ! IR = 2, NL

      RETURN
      END
