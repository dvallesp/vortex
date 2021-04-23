***********************************************************************
      SUBROUTINE CORRECT_OUTLIERS(NL,NX,NY,NZ,NPATCH,
     &           PATCHNX,PATCHNY,PATCHNZ, ERR_THR)
***********************************************************************
*     Cells with large relative error get their values obtained by
*     interpolation from a coarser grid.
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     Function parameters
      INTEGER NL, NX, NY, NZ
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      real err_thr

*     GLOBAL VARIABLES
*     Grid
      real RX(-2:NAMRX+3,NPALEV)
      real RY(-2:NAMRX+3,NPALEV)
      real RZ(-2:NAMRX+3,NPALEV)
      real RMX(-2:NAMRX+3,NPALEV)
      real RMY(-2:NAMRX+3,NPALEV)
      real RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
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
      real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P
*     rotational velocity (ROTS variables were reused)
      real U2R(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3R(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4R(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12R(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real U13R(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real U14R(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ U2R,U3R,U4R,U12R,U13R,U14R
*      original, total velocity
      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_ORIGINAL/ U2,U3,U4,U12,U13,U14

*     relative error (we assume it has been computed before)
      real ERR0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ERR1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /ERROR/ ERR0, ERR1

*     private variables
      INTEGER IR,I,N1,N2,N3,LOW1,LOW2,IX,JY,KZ
      real UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),FUIN
      real AAA, BBB, CCC, ERRBAS, BAS2, BAS3, BAS4, BAS
      real U2PBAS, U3PBAS, U4PBAS, U2RBAS, U3RBAS, U4RBAS, EU2, EU3, EU4
      INTEGER KARE,KR1,KR2,KR3

      integer c1, c2, c3

      CALL CELLWISE_ERROR(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ)

*     Level 1 patches: interpolate from the base grid
      IR=1
      c1 = 0
      c2 = 0
      c3 = 0

!$OMP PARALLEL DO SHARED(IR,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+         CR3AMR1X,CR3AMR1Y,CR3AMR1Z,U2P,U3P,U4P,U2R,U3R,U4R,
!$OMP+         U12P,U13P,U14P,U12R, U13R,U14R,U12,U13,U14,ERR1,ERR_THR,
!$OMP+         c1,c2,c3)
!$OMP+      PRIVATE(I,N1,N2,N3,IX,JY,KZ,UBAS,FUIN,KR1,KR2,KR3, ERRBAS,
!$OMP+              U2PBAS,U3PBAS,U4PBAS,U2RBAS,U3RBAS,U4RBAS,
!$OMP+              EU2,EU3,EU4,BAS2,BAS3,BAS4,BAS),
!$OMP+      DEFAULT(NONE)
      DO I=1,NPATCH(IR)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1
         c1 = c1+1

*       #######################################################
        IF (ERR1(IX,JY,KZ,I).GE.ERR_THR) THEN
*       #######################################################
          c2 = c2+1
          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          UBAS(1:3,1:3,1:3)=U2P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U2PBAS=FUIN

          UBAS(1:3,1:3,1:3)=U3P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U3PBAS=FUIN

          UBAS(1:3,1:3,1:3)=U4P(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U4PBAS=FUIN

          UBAS(1:3,1:3,1:3)=U2R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U2RBAS=FUIN

          UBAS(1:3,1:3,1:3)=U3R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U3RBAS=FUIN

          UBAS(1:3,1:3,1:3)=U4R(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U4RBAS=FUIN

          eu2 =  (u2pbas+u2rbas) / u12(ix,jy,kz,i) - 1.0
          eu3 =  (u3pbas+u3rbas) / u13(ix,jy,kz,i) - 1.0
          eu4 =  (u4pbas+u4rbas) / u14(ix,jy,kz,i) - 1.0

          bas2 = u12(ix,jy,kz,i) ** 2
          bas3 = u13(ix,jy,kz,i) ** 2
          bas4 = u14(ix,jy,kz,i) ** 2
          bas = bas2 + bas3 + bas4

          bas2 = (bas2 / bas * eu2) ** 2
          bas3 = (bas3 / bas * eu3) ** 2
          bas4 = (bas4 / bas * eu4) ** 2

          errbas = sqrt(bas2 + bas3 + bas4)

          if (errbas.lt.ERR1(IX,JY,KZ,I)) then
            u12p(ix,jy,kz,i) = u2pbas
            u13p(ix,jy,kz,i) = u3pbas
            u14p(ix,jy,kz,i) = u4pbas
            u12r(ix,jy,kz,i) = u2rbas
            u13r(ix,jy,kz,i) = u3rbas
            u14r(ix,jy,kz,i) = u4rbas
            err1(ix,jy,kz,i) = errbas
            c3 = c3+1
          end if

*       #######################################################
        END IF
*       #######################################################

       END DO
       END DO
       END DO

      END DO

      write(*,*) ir,c1, c2,c3

*     Other levels (l>1)
      DO IR=2,NL
        c1 = 0
        c2 = 0
        c3 = 0

      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(IR,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+         LOW1,LOW2,RX,RY,RZ,RADX,RADY,RADZ,CR3AMR1,CR3AMR1X,
!$OMP+         CR3AMR1Y,CR3AMR1Z,U2P,U3P,U4P,U2R,U3R,U4R,U12P,U13P,
!$OMP+         U14P,U12R,U13R,U14R,U12,U13,U14,ERR1,ERR_THR,c1,c2,c3),
!$OMP+     PRIVATE(I,N1,N2,N3,IX,JY,KZ,KARE,KR1,KR2,KR3,UBAS,RXBAS,
!$OMP+             RYBAS,RZBAS,FUIN,AAA,BBB,CCC,ERRBAS,U2PBAS,U3PBAS,
!$OMP+             U4PBAS,U2RBAS,U3RBAS,U4RBAS,EU2,EU3,EU4,BAS2,BAS3,
!$OMP+             BAS4,BAS),
!$OMP+     DEFAULT(NONE)
      DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1
         c1 = c1 + 1

*       #######################################################
        IF (ERR1(IX,JY,KZ,I).GE.ERR_THR) THEN
*       #######################################################
          c2 = c2 + 1
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
          U2PBAS=FUIN

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
          U3PBAS=FUIN

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
          U4PBAS=FUIN

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
          U2RBAS=FUIN

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
          U3RBAS=FUIN

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
          U4RBAS=FUIN

          eu2 =  (u2pbas+u2rbas) / u12(ix,jy,kz,i) - 1.0
          eu3 =  (u3pbas+u3rbas) / u13(ix,jy,kz,i) - 1.0
          eu4 =  (u4pbas+u4rbas) / u14(ix,jy,kz,i) - 1.0

          bas2 = u12(ix,jy,kz,i) ** 2
          bas3 = u13(ix,jy,kz,i) ** 2
          bas4 = u14(ix,jy,kz,i) ** 2
          bas = bas2 + bas3 + bas4

          bas2 = (bas2 / bas * eu2) ** 2
          bas3 = (bas3 / bas * eu3) ** 2
          bas4 = (bas4 / bas * eu4) ** 2

          errbas = sqrt(bas2 + bas3 + bas4)

          if (errbas.lt.ERR1(IX,JY,KZ,I)) then
            u12p(ix,jy,kz,i) = u2pbas
            u13p(ix,jy,kz,i) = u3pbas
            u14p(ix,jy,kz,i) = u4pbas
            u12r(ix,jy,kz,i) = u2rbas
            u13r(ix,jy,kz,i) = u3rbas
            u14r(ix,jy,kz,i) = u4rbas
            err1(ix,jy,kz,i) = errbas
             c3 = c3 + 1
          end if

*       #######################################################
        END IF
*       #######################################################

       END DO
       END DO
       END DO

      END DO ! I = LOW1, LOW2
       write(*,*) ir, c1, c2, c3
      END DO ! IR = 2, NL

      RETURN
      END
