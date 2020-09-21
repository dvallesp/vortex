***********************************************************************
      SUBROUTINE MALLA(NX,NY,NZ,LADO)
***********************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,I,NY,J,NZ,K,PA4
      real A,B,C,LADO

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ


*     GENERAL INITIAL CONDITIONS
*     GRID LIMITS
       A=-LADO/2.0
       B=LADO/2.0


*     GRID

*     X-AXIS
      C=(B-A)/(NX-1)
      RADX(1)=A
      DO I=2,NX
        RADX(I)=RADX(I-1)+C
      END DO

*     FICTICIUS CELLS
      RADX(0)=RADX(1)-C
      RADX(NX+1)=RADX(NX)+C

*     Y-AXIS
      C=(B-A)/(NY-1)
      RADY(1)=A
      DO J=2,NY
        RADY(J)=RADY(J-1)+C
      END DO

*     FICTICIUS CELLS
      RADY(0)=RADY(1)-C
      RADY(NY+1)=RADY(NY)+C
*     Z-AXIS
      C=(B-A)/(NZ-1)
      RADZ(1)=A
      DO K=2,NZ
        RADZ(K)=RADZ(K-1)+C
      END DO

*     FICTICIUS CELLS
      RADZ(0)=RADZ(1)-C
      RADZ(NZ+1)=RADZ(NZ)+C


*     COORDINATE FOR INTERFACES ***************************************
      DO I=0,NX
        RADMX(I) = (RADX(I)+RADX(I+1))/2.D0
      END DO
      DO J=0,NY
        RADMY(J) = (RADY(J)+RADY(J+1))/2.D0
      END DO
      DO K=0,NZ
        RADMZ(K) = (RADZ(K)+RADZ(K+1))/2.D0
      END DO

      DX=RADX(2)-RADX(1)
      DY=RADY(2)-RADY(1)
      DZ=RADZ(2)-RADZ(1)


      RETURN

      END

************************************************************************
      SUBROUTINE GRIDAMR(NX,NY,NZ,NL,NPATCH,
     &                   PATCHNX,PATCHNY,PATCHNZ,
     &                   PATCHX,PATCHY,PATCHZ,
     &                   PATCHRX,PATCHRY,PATCHRZ,PARE)
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,NY,NZ,NL1,NL

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PARE(NPALEV)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      real  PATCHRX(NPALEV)
      real  PATCHRY(NPALEV)
      real  PATCHRZ(NPALEV)

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

      INTEGER I,J,K,IX,JY,KZ,I1,J1,K1,IR,IPALE,BOR
      INTEGER NEF,NCELL,PAX1,PAX2,PAY1,PAY2,PAZ1,PAZ2
      INTEGER IPATCH,II,JJ,KK,N1,N2,N3
      INTEGER INMAX(3)
      INTEGER NBAS,IPA2,NPX,NPY,NPZ
      real DXPA,DYPA,DZPA

      real RX(-2:NAMRX+3,NPALEV)
      real RY(-2:NAMRX+3,NPALEV)
      real RZ(-2:NAMRX+3,NPALEV)
      real RMX(-2:NAMRX+3,NPALEV)
      real RMY(-2:NAMRX+3,NPALEV)
      real RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      real XX,YY,ZZ,XX1,YY1,ZZ1,XX2,YY2,ZZ2,BAS1,BAS2,BAS3

      INTEGER CONTROL,IP,P1,P2,P3,NP1,NP2,NP3
      INTEGER L1,L2,L3,LL1,LL2,LL3,CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER LN1,LN2,LN3,LNN1,LNN2,LNN3
      INTEGER LOW1,LOW2,LOW3,LOW4
      INTEGER KR1,KR2,KR3,ABUELO,MARCA

*      ---PARALLEL---
      INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
      COMMON /PROCESADORES/ NUM


*     modificado 10/5/2012
*     marca celdas al borde en niveles recursivos
*     solo hasta el nivel 2 porque el nivel 1 nunca se acerca
*     a las caras de la malla base


      DO IR=NL,2,-1
      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(IR,NL,LOW1,LOW2,PATCHX,PATCHY,PATCHZ,
!$OMP+            PATCHNX,PATCHNY,PATCHNZ,CR3AMR1,NX,NY,NZ,
!$OMP+            CR3AMR1X,CR3AMR1Y,CR3AMR1Z,PARE),
!$OMP+ PRIVATE(I,L1,L2,L3,IX,JY,KZ,KR1,KR2,KR3,MARCA,
!$OMP+         CR1,CR2,CR3,ABUELO)
      DO I=LOW1,LOW2
       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)
       DO KZ=-2,PATCHNZ(I)+3
       DO JY=-2,PATCHNY(I)+3
       DO IX=-2,PATCHNX(I)+3

       MARCA=0

       CR1=L1-1+INT((IX+1)/2)
       CR2=L2-1+INT((JY+1)/2)
       CR3=L3-1+INT((KZ+1)/2)

       IF (IX.LT.-1) CR1=INT((IX+1)/2)+L1-2
       IF (JY.LT.-1) CR2=INT((JY+1)/2)+L2-2
       IF (KZ.LT.-1) CR3=INT((KZ+1)/2)+L3-2

       KR1=CR1
       KR2=CR2
       KR3=CR3

       ABUELO=PARE(I)

       DO WHILE (MARCA.EQ.0)  !...................
        IF (CR1.LT.2.OR.CR1.GT.PATCHNX(ABUELO)-1.OR.   !----------------
     &      CR2.LT.2.OR.CR2.GT.PATCHNY(ABUELO)-1.OR.
     &      CR3.LT.2.OR.CR3.GT.PATCHNZ(ABUELO)-1) THEN
*        la celda madre de la celda ix,jy,kz
*        esta en el borde de su parche, hay que buscar
*        recursivamente que celda no esta en la cara

         CR1=PATCHX(ABUELO)-1+INT((KR1+1)/2)
         CR2=PATCHY(ABUELO)-1+INT((KR2+1)/2)
         CR3=PATCHZ(ABUELO)-1+INT((KR3+1)/2)

         IF (KR1.LT.-1) CR1=PATCHX(ABUELO)-2+INT((KR1+1)/2)
         IF (KR2.LT.-1) CR2=PATCHY(ABUELO)-2+INT((KR2+1)/2)
         IF (KR3.LT.-1) CR3=PATCHZ(ABUELO)-2+INT((KR3+1)/2)

         KR1=CR1
         KR2=CR2
         KR3=CR3

         ABUELO=PARE(ABUELO)

         IF (ABUELO.EQ.0) THEN
          MARCA=1
          CR3AMR1(IX,JY,KZ,I)=0
          CR3AMR1X(IX,JY,KZ,I)=KR1
          CR3AMR1Y(IX,JY,KZ,I)=KR2
          CR3AMR1Z(IX,JY,KZ,I)=KR3

*          IF (KR1.LT.1) WRITE(*,*) 'X OUT <', ABUELO,IR,I,KR1
*          IF (KR2.LT.1) WRITE(*,*) 'Y OUT <', ABUELO,IR,I,KR2
*          IF (KR3.LT.1) WRITE(*,*) 'Z OUT <', ABUELO,IR,I,KR3
*
*          IF (KR1.GT.NX) WRITE(*,*) 'X OUT >', ABUELO,IR,I,KR1
*          IF (KR2.GT.NY) WRITE(*,*) 'Y OUT >', ABUELO,IR,I,KR2
*          IF (KR3.GT.NZ) WRITE(*,*) 'Z OUT >', ABUELO,IR,I,KR3

         END IF
        ELSE                                           !----------------
         CR3AMR1(IX,JY,KZ,I)=ABUELO
         CR3AMR1X(IX,JY,KZ,I)=KR1
         CR3AMR1Y(IX,JY,KZ,I)=KR2
         CR3AMR1Z(IX,JY,KZ,I)=KR3

*         IF (KR1.LT.1) WRITE(*,*) 'X OUT L <', ABUELO,IR,I,KR1
*         IF (KR2.LT.1) WRITE(*,*) 'Y OUT L <', ABUELO,IR,I,KR2
*         IF (KR3.LT.1) WRITE(*,*) 'Z OUT L <', ABUELO,IR,I,KR3
*
*         IF (KR1.GT.PATCHNX(ABUELO))
*     &     WRITE(*,*) 'X OUT L >', ABUELO,IR,I,KR1,PATCHNX(ABUELO)
*         IF (KR2.GT.PATCHNY(ABUELO))
*     &     WRITE(*,*) 'Y OUT L >', ABUELO,IR,I,KR2,PATCHNY(ABUELO)
*         IF (KR3.GT.PATCHNZ(ABUELO))
*     &     WRITE(*,*) 'Z OUT L >', ABUELO,IR,I,KR3,PATCHNZ(ABUELO)

*        esta celda si esta bien dentro de su parche padre
         MARCA=1                                       !----------------
        END IF
       END DO                 !...................

       END DO
       END DO
       END DO
      END DO
      END DO


      IR=1
!$OMP PARALLEL DO SHARED(IR,NL,PATCHX,PATCHY,PATCHZ,
!$OMP+            PATCHNX,PATCHNY,PATCHNZ,CR3AMR1,NX,NY,NZ,
!$OMP+            CR3AMR1X,CR3AMR1Y,CR3AMR1Z,PARE),
!$OMP+ PRIVATE(I,L1,L2,L3,IX,JY,KZ,KR1,KR2,KR3,CR1,CR2,CR3)
      DO I=1,NPATCH(IR)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       DO KZ=-2,PATCHNZ(I)+3
       DO JY=-2,PATCHNY(I)+3
       DO IX=-2,PATCHNX(I)+3


       CR1=L1-1+INT((IX+1)/2)
       CR2=L2-1+INT((JY+1)/2)
       CR3=L3-1+INT((KZ+1)/2)

       IF (IX.LT.-1) CR1=INT((IX+1)/2)+L1-2
       IF (JY.LT.-1) CR2=INT((JY+1)/2)+L2-2
       IF (KZ.LT.-1) CR3=INT((KZ+1)/2)+L3-2

       KR1=CR1
       KR2=CR2
       KR3=CR3


       CR3AMR1(IX,JY,KZ,I)=0
       CR3AMR1X(IX,JY,KZ,I)=KR1
       CR3AMR1Y(IX,JY,KZ,I)=KR2
       CR3AMR1Z(IX,JY,KZ,I)=KR3

*       IF (KR1.LT.1) WRITE(*,*) '0 X OUT <', ABUELO,IR,I,KR1
*       IF (KR2.LT.1) WRITE(*,*) '0 Y OUT <', ABUELO,IR,I,KR2
*       IF (KR3.LT.1) WRITE(*,*) '0 Z OUT <', ABUELO,IR,I,KR3
*
*       IF (KR1.GT.NX) WRITE(*,*) '0 X OUT >', ABUELO,IR,I,KR1
*       IF (KR2.GT.NY) WRITE(*,*) '0 Y OUT >', ABUELO,IR,I,KR2
*       IF (KR3.GT.NZ) WRITE(*,*) '0 Z OUT >', ABUELO,IR,I,KR3

       END DO
       END DO
       END DO
      END DO




*     MINI GRIDS
      RX=0.0
      RY=0.0
      RZ=0.0
      RMX=0.0
      RMY=0.0
      RMZ=0.0
      DO IR=1,NL
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

        CALL MINIMALLA(N1,N2,N3,DXPA,DYPA,DZPA,PATCHRX(I),
     &        PATCHRY(I),PATCHRZ(I),
     &        RX(-2:N1+3,I),RY(-2:N2+3,I),RZ(-2:N3+3,I),
     &        RMX(-2:N1+3,I),RMY(-2:N2+3,I),RMZ(-2:N3+3,I))

       END DO
      END DO


      RETURN
      END

**********************************************************************
      SUBROUTINE MINIMALLA(N1,N2,N3,DX,DY,DZ,RPAX,RPAY,RPAZ,
     &                     RX,RY,RZ,RMX,RMY,RMZ)
**********************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER N1,N2,N3,I
      real RX(-2:N1+3),RY(-2:N2+3),RZ(-2:N3+3)
      real RMX(-2:N1+3),RMY(-2:N2+3),RMZ(-2:N3+3)
      real RPAX,RPAY,RPAZ,DX,DY,DZ

*     NEW DX
      DO I=1,N1
        RX(I)=RPAX-(DX/2.0)+(I-1)*DX
      END DO
      RX(0)=RX(1)-DX
      RX(N1+1)=RX(N1)+DX
      RX(-1)=RX(0)-DX
      RX(N1+2)=RX(N1+1)+DX
      RX(-2)=RX(-1)-DX
      RX(N1+3)=RX(N1+2)+DX

      DO I=-2,N1+2
        RMX(I)=(RX(I)+RX(I+1))/2.0
      END DO

      DO I=1,N2
        RY(I)=RPAY-(DY/2.0)+(I-1)*DY
      END DO
      RY(0)=RY(1)-DY
      RY(N2+1)=RY(N2)+DY
      RY(-1)=RY(0)-DY
      RY(N2+2)=RY(N2+1)+DY
      RY(-2)=RY(-1)-DY
      RY(N2+3)=RY(N2+2)+DY

      DO I=-2,N2+2
        RMY(I)=(RY(I)+RY(I+1))/2.0
      END DO

      DO I=1,N3
        RZ(I)=RPAZ-(DZ/2.0)+(I-1)*DZ
      END DO
      RZ(0)=RZ(1)-DZ
      RZ(N3+1)=RZ(N3)+DZ
      RZ(-1)=RZ(0)-DZ
      RZ(N3+2)=RZ(N3+1)+DZ
      RZ(-2)=RZ(-1)-DZ
      RZ(N3+3)=RZ(N3+2)+DZ

      DO I=-2,N3+2
        RMZ(I)=(RZ(I)+RZ(I+1))/2.0
      END DO


      RETURN
      END

***********************************************************************
       SUBROUTINE EXTEND_VAR(NX,NY,NZ,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************

       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ,I,J,K,IR,IX,JY,KZ
       INTEGER NL,N1,N2,N3,L1,L2,L3
       INTEGER II,JJ,KK,LOW1, LOW2

       real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &         RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &         RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC/ U2,U3,U4,U12,U13,U14

       real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

       real DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       real DXPA,DYPA,DZPA,XXX1,YYY1,ZZZ1
       INTEGER CR1,CR2,CR3
       INTEGER MARK, ABUELO, KR1, KR2, KR3, KARE,IR_ABUE

       real UBAS(1:3,1:3,1:3),FUIN,RXBAS(3),RYBAS(3),RZBAS(3)
       real AAA,BBB,CCC

       INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

       real RX(-2:NAMRX+3,NPALEV)
       real RY(-2:NAMRX+3,NPALEV)
       real RZ(-2:NAMRX+3,NPALEV)
       real RMX(-2:NAMRX+3,NPALEV)
       real RMY(-2:NAMRX+3,NPALEV)
       real RMZ(-2:NAMRX+3,NPALEV)
       COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ


*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM


*      EXPANSION DE LAS MATRICES POR INTERPOLACION

      IR=1
      DXPA=DX/(2.0**IR)
      DYPA=DY/(2.0**IR)
      DZPA=DZ/(2.0**IR)

      DO I=1,NPATCH(IR)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=0,N3+1
       DO JY=0,N2+1
       DO IX=0,N1+1

*      #######################################################
       IF (IX.LT.1.OR.IX.GT.N1.OR.JY.LT.1.OR.JY.GT.N2.OR.
     &     KZ.LT.1.OR.KZ.GT.N3) THEN
*      #######################################################

        KR1=CR3AMR1X(IX,JY,KZ,I)
        KR2=CR3AMR1Y(IX,JY,KZ,I)
        KR3=CR3AMR1Z(IX,JY,KZ,I)

        !VX
        UBAS(1:3,1:3,1:3)=U2(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
        CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
        U12(IX,JY,KZ,I)=FUIN

        !VY
        UBAS(1:3,1:3,1:3)=U3(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
        CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
        U13(IX,JY,KZ,I)=FUIN

        !VZ
        UBAS(1:3,1:3,1:3)=U4(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
        CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
        U14(IX,JY,KZ,I)=FUIN

*      #######################################################
        END IF
*      #######################################################

        END DO
        END DO
        END DO
       END DO


       DO IR=2,NL
        DXPA=DX/(2.0**IR)
        DYPA=DY/(2.0**IR)
        DZPA=DZ/(2.0**IR)

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

         DO KZ=0,N3+1
         DO JY=0,N2+1
         DO IX=0,N1+1

         IF (IX.LT.1.OR.IX.GT.N1.OR.
     &       JY.LT.1.OR.JY.GT.N2.OR.
     &       KZ.LT.1.OR.KZ.GT.N3) THEN

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

             !Vx:
              UBAS(1:3,1:3,1:3)=
     &             U12(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
              CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                             RXBAS,RYBAS,RZBAS,UBAS,FUIN)
              U12(IX,JY,KZ,I)=FUIN

              !Vy:
              UBAS(1:3,1:3,1:3)=
     &             U13(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
              CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
              U13(IX,JY,KZ,I)=FUIN

              !Vz:
              UBAS(1:3,1:3,1:3)=
     &             U14(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
              CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
              U14(IX,JY,KZ,I)=FUIN


             ELSE
              RXBAS(1:3)=RADX(KR1-1:KR1+1)
              RYBAS(1:3)=RADY(KR2-1:KR2+1)
              RZBAS(1:3)=RADZ(KR3-1:KR3+1)

              !Vx
              UBAS(1:3,1:3,1:3)=
     &             U2(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
              CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
              U12(IX,JY,KZ,I)=FUIN

              !Vy
              UBAS(1:3,1:3,1:3)=
     &             U3(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
              CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
              U13(IX,JY,KZ,I)=FUIN

              !Vz
              UBAS(1:3,1:3,1:3)=
     &            U4(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
              CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
              U14(IX,JY,KZ,I)=FUIN

             ENDIF
c          ELSE
c              U12(IX,JY,KZ,I)=U12_Z(IX,JY,KZ,I)
c              U13(IX,JY,KZ,I)=U13_Z(IX,JY,KZ,I)
c              U14(IX,JY,KZ,I)=U14_Z(IX,JY,KZ,I)
          ENDIF
         END DO
         END DO
         END DO

        END DO
       END DO


* IR=0
       DO K=0,NZ+1
       DO J=0,NY+1
       DO I=0,NX+1

        IX=I
        JY=J
        KZ=K
        IF (I.LT.1) IX=I+NX
        IF (J.LT.1) JY=J+NY
        IF (K.LT.1) KZ=K+NZ
        IF (I.GT.NX) IX=I-NX
        IF (J.GT.NY) JY=J-NY
        IF (K.GT.NZ) KZ=K-NZ

        U2(I,J,K)=U2(IX,JY,KZ)
        U3(I,J,K)=U3(IX,JY,KZ)
        U4(I,J,K)=U4(IX,JY,KZ)
       END DO
       END DO
       END DO

***    ALL VARIALBES ARE EXTENDED ONE CELL

       RETURN
       END



***********************************************************************
*******   END PROGRAM  ************************************************
***********************************************************************
