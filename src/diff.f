***********************************************************************
       SUBROUTINE ROTARY(NX,NY,NZ,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************
*     Computes the rotational of the velocity field (U2,U3,U4,U12,U13,
*     U14) using a first-order centered scheme.
************************************************************************

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
       real  PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       real DXPA,DYPA,DZPA,XXX1,YYY1,ZZZ1
       real BAS21,BAS32,BAS33,BAS43,AAA,BBB,CCC
       INTEGER CR1,CR2,CR3
       INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       DO IR=1,NL

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,U12,U13,U14,
!$OMP+                   ROTAX_1,ROTAY_1,ROTAZ_1,DXPA,DYPA,DZPA),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ,BAS21,BAS32,BAS43),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

        BAS21=0.0
        BAS32=0.0
        BAS43=0.0
* X
        BAS21=U14(IX,JY+1,KZ,I)-U14(IX,JY-1,KZ,I)
        BAS32=U13(IX,JY,KZ+1,I)-U13(IX,JY,KZ-1,I)

        ROTAX_1(IX,JY,KZ,I)=(BAS21-BAS32)/(2.0*DXPA)

* Y
        BAS21=U12(IX,JY,KZ+1,I)-U12(IX,JY,KZ-1,I)
        BAS32=U14(IX+1,JY,KZ,I)-U14(IX-1,JY,KZ,I)

        ROTAY_1(IX,JY,KZ,I)=(BAS21-BAS32)/(2.0*DYPA)

* Z
        BAS21=U13(IX+1,JY,KZ,I)-U13(IX-1,JY,KZ,I)
        BAS32=U12(IX,JY+1,KZ,I)-U12(IX,JY-1,KZ,I)

        ROTAZ_1(IX,JY,KZ,I)=(BAS21-BAS32)/(2.0*DZPA)

       END DO
       END DO
       END DO

       END DO
       END DO

*-------------------------------*
*      COARSE LEVEL
*-------------------------------*

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2,U3,U4,DX,DY,DZ,ROTAX_0,ROTAY_0,
!$OMP+                   ROTAZ_0),
!$OMP+            PRIVATE(IX,JY,KZ,BAS21,BAS32),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX

        BAS21=0.0
        BAS32=0.0

* X
        BAS21=U4(IX,JY+1,KZ)-U4(IX,JY-1,KZ)
        BAS32=U3(IX,JY,KZ+1)-U3(IX,JY,KZ-1)

        ROTAX_0(IX,JY,KZ)=(BAS21-BAS32)/(2.0*DX)

* Y
        BAS21=U2(IX,JY,KZ+1)-U2(IX,JY,KZ-1)
        BAS32=U4(IX+1,JY,KZ)-U4(IX-1,JY,KZ)

        ROTAY_0(IX,JY,KZ)=(BAS21-BAS32)/(2.0*DY)

* Z
        BAS21=U3(IX+1,JY,KZ)-U3(IX-1,JY,KZ)
        BAS32=U2(IX,JY+1,KZ)-U2(IX,JY-1,KZ)

        ROTAZ_0(IX,JY,KZ)=(BAS21-BAS32)/(2.0*DZ)

       END DO
       END DO
       END DO

        RETURN
       END

***********************************************************************
       SUBROUTINE DIVER_FINA(NX,NY,NZ,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************
*     Computes the divergence of the velocity field (U2,U3,U4,U12,U13,
*     U14) using a first-order centered scheme.
************************************************************************

       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ
       INTEGER I,J,K,IR,IX,JY,KZ
       INTEGER NL,N1,N2,N3,L1,L2,L3
       INTEGER II,JJ,KK, LOW1, LOW2

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

       real DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       COMMON /DIVERGENCE/ DIVER0, DIVER

       real DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       real DXPA,DYPA,DZPA,XXX1,YYY1,ZZZ1
       real BAS21,BAS32,BAS33,BAS43,AAA
       INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE
       INTEGER CR1,CR2,CR3

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       DO IR=1,NL

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,U12,U13,U14,
!$OMP+                   DIVER,DXPA),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ,BAS21,BAS32,BAS43),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

        BAS21=0.0
        BAS32=0.0
        BAS43=0.0

        BAS21=U12(IX+1,JY,KZ,I)-U12(IX-1,JY,KZ,I)
        BAS32=U13(IX,JY+1,KZ,I)-U13(IX,JY-1,KZ,I)
        BAS43=U14(IX,JY,KZ+1,I)-U14(IX,JY,KZ-1,I)

        DIVER(IX,JY,KZ,I)=(BAS21+BAS32+BAS43)/(2.0*DXPA)

       END DO
       END DO
       END DO

       END DO
       END DO

*-------------------------------*
*      COARSE LEVEL
*-------------------------------*

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2,U3,U4,DIVER0,DX),
!$OMP+            PRIVATE(IX,JY,KZ,BAS21,BAS32,BAS43),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX

        BAS21=0.0
        BAS32=0.0
        BAS43=0.0

        BAS21=U2(IX+1,JY,KZ)-U2(IX-1,JY,KZ)
        BAS32=U3(IX,JY+1,KZ)-U3(IX,JY-1,KZ)
        BAS43=U4(IX,JY,KZ+1)-U4(IX,JY,KZ-1)

        DIVER0(IX,JY,KZ)=(BAS21+BAS32+BAS43)/(2.0*DX)

       END DO
       END DO
       END DO

        RETURN
       END


***********************************************************************
       SUBROUTINE GRADIENTE(NX,NY,NZ,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************
*     Computes the gradient of the scalar potential (stored in DIVER0,
*     DIVER) using a first-order scheme.
***********************************************************************

       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ
       INTEGER I,J,K,IR,IX,JY,KZ
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

       real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P

       real DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       COMMON /DIVERGENCE/ DIVER0, DIVER

       real DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real  PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       real DXPA,DYPA,DZPA,XXX1,YYY1,ZZZ1
       real BAS21,BAS32,BAS33,BAS43,AAA
       INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE,CR1,CR2,CR3
       real GRAD_P_X,GRAD_P_Y,GRAD_P_Z


*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

      GRAD_P_X=0.0
      GRAD_P_Y=0.0
      GRAD_P_Z=0.0


       DO IR=1,NL

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0

       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,DXPA,DYPA,
!$OMP+                   DZPA,U12P,U13P,U14P,DIVER),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ,BAS21,BAS32,BAS43,
!$OMP+                    GRAD_P_X,GRAD_P_Y,GRAD_P_Z),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

        BAS21=0.0
        BAS32=0.0
        BAS43=0.0

* X
        GRAD_P_X=DIVER(IX+1,JY,KZ,I)-DIVER(IX-1,JY,KZ,I)
        GRAD_P_X=GRAD_P_X/(2.0*DXPA)
        U12P(IX,JY,KZ,I)=GRAD_P_X

* Y
        GRAD_P_Y=DIVER(IX,JY+1,KZ,I)-DIVER(IX,JY-1,KZ,I)
        GRAD_P_Y=GRAD_P_Y/(2.0*DYPA)
        U13P(IX,JY,KZ,I)=GRAD_P_Y

* Z
        GRAD_P_Z=DIVER(IX,JY,KZ+1,I)-DIVER(IX,JY,KZ-1,I)
        GRAD_P_Z=GRAD_P_Z/(2.0*DZPA)
        U14P(IX,JY,KZ,I)=GRAD_P_Z

       END DO
       END DO
       END DO

       END DO
       END DO

*-------------------------------*
*      COARSE LEVEL
*-------------------------------*

!$OMP PARALLEL DO SHARED(NX,NY,NZ,DX,DY,DZ,DIVER0,U2P,U3P,U4P),
!$OMP+            PRIVATE(IX,JY,KZ,BAS21,BAS32,GRAD_P_X,GRAD_P_Y,
!$OMP+                    GRAD_P_Z),
!$OMP+            DEFAULT(NONE)
       DO KZ=1, NZ
       DO JY=1, NY
       DO IX=1, NX

        BAS21=0.0
        BAS32=0.0

        GRAD_P_X=0.0
        GRAD_P_Y=0.0
        GRAD_P_Z=0.0

* X
        GRAD_P_X=DIVER0(IX+1,JY,KZ)-DIVER0(IX-1,JY,KZ)
        GRAD_P_X=GRAD_P_X/(2.0*DX)
        U2P(IX,JY,KZ)=GRAD_P_X

* Y
        GRAD_P_Y=DIVER0(IX,JY+1,KZ)-DIVER0(IX,JY-1,KZ)
        GRAD_P_Y=GRAD_P_Y/(2.0*DY)
        U3P(IX,JY,KZ)=GRAD_P_Y

* Z
        GRAD_P_Z=DIVER0(IX,JY,KZ+1)-DIVER0(IX,JY,KZ-1)
        GRAD_P_Z=GRAD_P_Z/(2.0*DZ)
        U4P(IX,JY,KZ)=GRAD_P_Z

       END DO
       END DO
       END DO


        RETURN
       END
