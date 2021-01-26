************************************************************************
      SUBROUTINE CELLWISE_ERROR(NX,NY,NZ,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ)
************************************************************************
*     Computes the cell-wise relative error in reconstructing the
*     velocity field, as ABS((UP + UR) / U - 1), being U, UP and UR
*     the total, compressive and solenoidal velocity components,
*     respectively. This is done for each component and the total
*     error is found by a weighted mean (corresponds to variance
*     propagation to sqrt(u2**2 + u3**2 + u4**2)).
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX, NY, NZ, NL
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     global variables
*     original velocity
      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_ORIGINAL/ U2,U3,U4,U12,U13,U14
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

*     variable to store the relative error
      real ERR0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ERR1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /ERROR/ ERR0, ERR1

*     Auxiliary variables
      real EU2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real EU3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real EU4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real EU12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real EU13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real EU14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)

*     private VARIABLES
      INTEGER IX, JY, KZ, I, J, K, LOW1, LOW2, N1, N2, N3, IR
      real BAS, BAS1, BAS2, BAS3


*     UP will contain the relative error per component!!!!

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2P,U3P,U4P,U2R,U3R,U4R,
!$OMP+                   U2,U3,U4,EU2,EU3,EU4),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
       DO K=1, NZ
       DO J=1, NY
       DO I=1, NX
         IF(U2(I,J,K).NE.0.0) EU2(I,J,K) = (U2P(I,J,K) + U2R(I,J,K)) /
     &                                     U2(I,J,K) - 1.0
         IF(U3(I,J,K).NE.0.0) EU3(I,J,K)= (U3P(I,J,K) + U3R(I,J,K)) /
     &                                     U3(I,J,K) - 1.0
         IF(U4(I,J,K).NE.0.0) EU4(I,J,K)= (U4P(I,J,K) + U4R(I,J,K)) /
     &                                     U4(I,J,K) - 1.0
       END DO
       END DO
       END DO

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U12P,U13P,U14P,U12R,U13R,U14R,U12,U13,U14,
!$OMP+                   EU12,EU13,EU14),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
         DO KZ=1, N3
         DO JY=1, N2
         DO IX=1, N1
           IF(U12(IX,JY,KZ,I).NE.0.0) THEN
             EU12(IX,JY,KZ,I) = (U12P(IX,JY,KZ,I) + U12R(IX,JY,KZ,I)) /
     &                          U12(IX,JY,KZ,I) - 1.0
           END IF
           IF(U13(IX,JY,KZ,I).NE.0.0) THEN
             EU13(IX,JY,KZ,I) = (U13P(IX,JY,KZ,I) + U13R(IX,JY,KZ,I)) /
     &                          U13(IX,JY,KZ,I) - 1.0
           END IF
           IF(U14(IX,JY,KZ,I).NE.0.0) THEN
             EU14(IX,JY,KZ,I) = (U14P(IX,JY,KZ,I) + U14R(IX,JY,KZ,I)) /
     &                          U14(IX,JY,KZ,I) - 1.0
           END IF
         END DO
         END DO
         END DO

        END DO
        END DO

*       Compute the relative error

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2P,U3P,U4P,U2,U3,U4,ERR0,
!$OMP+                   EU2,EU3,EU4),
!$OMP+            PRIVATE(I,J,K,BAS,BAS1,BAS2,BAS3),
!$OMP+            DEFAULT(NONE)
       DO K=1, NZ
       DO J=1, NY
       DO I=1, NX
         BAS1 = U2(I,J,K)**2
         BAS2 = U3(I,J,K)**2
         BAS3 = U4(I,J,K)**2
         BAS = BAS1 + BAS2 + BAS3

         IF(BAS.NE.0.0) THEN
           BAS1 = (BAS1 / BAS * EU2(I,J,K)) ** 2
           BAS2 = (BAS2 / BAS * EU3(I,J,K)) ** 2
           BAS3 = (BAS3 / BAS * EU4(I,J,K)) ** 2
         ELSE
           BAS1 = 0
           BAS2 = 0
           BAS3 = 0
         END IF

         ERR0(I,J,K) = SQRT(BAS1 + BAS2 + BAS3)
       END DO
       END DO
       END DO

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U12P,U13P,U14P,U12,U13,U14,ERR1,
!$OMP+                   EU12,EU13,EU14),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I,BAS,BAS1,BAS2,BAS3),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
         DO KZ=1, N3
         DO JY=1, N2
         DO IX=1, N1
           BAS1 = U12(IX,JY,KZ,I)**2
           BAS2 = U13(IX,JY,KZ,I)**2
           BAS3 = U14(IX,JY,KZ,I)**2
           BAS = BAS1 + BAS2 + BAS3

           IF(BAS.NE.0) THEN
             BAS1 = (BAS1 / BAS * EU12(IX,JY,KZ,I)) ** 2
             BAS2 = (BAS2 / BAS * EU13(IX,JY,KZ,I)) ** 2
             BAS3 = (BAS3 / BAS * EU14(IX,JY,KZ,I)) ** 2
           ELSE
             BAS1 = 0
             BAS2 = 0
             BAS3 = 0
           END IF

           ERR1(IX,JY,KZ,I) = SQRT(BAS1 + BAS2 + BAS3)
         END DO
         END DO
         END DO

        END DO
        END DO

      RETURN
      END


************************************************************************
       SUBROUTINE SYNC_AMR_VELOCITIES(NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)
************************************************************************
*     Ensures overlapping cells at a given refinement levels have
*     the same values of the reconstructed velocities. This
*     subroutine is based on the last version of the "VEINSGRID"
*     subroutine in MASCLET to detect the overlaps.
************************************************************************

*     modified VQ 23-4-2020, DV 23-4-2020

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV)
       INTEGER PATCHNY(NPALEV)
       INTEGER PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV)
       INTEGER PATCHY(NPALEV)
       INTEGER PATCHZ(NPALEV)
       real  PATCHRX(NPALEV)
       real  PATCHRY(NPALEV)
       real  PATCHRZ(NPALEV)
       INTEGER PARE(NPALEV)

       INTEGER CR1,CR2,CR3,CR4,CR5,CR6
       INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK
       INTEGER N1,N2,N3,L1,L2,L3, NL
       INTEGER NN1,NN2,NN3,LL1,LL2,LL3
       INTEGER KZ2,JY2,IX2,I2
       INTEGER, ALLOCATABLE::VECINO(:,:)
       INTEGER, ALLOCATABLE::NVECI(:)
       INTEGER NV,A2,B2,C2,K,LOW1, LOW2

       INTEGER MARCA(NAMRX,NAMRY,NAMRZ,NPALEV)

       real A1,B1,C1,RIV1,RIV2,RIV3
       INTEGER CONTROL(3)
       INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
       INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
       INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
       real RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
       real RZ1,RZZ1,RZ2,RZZ2

       INTEGER NX,NY,NZ,ITER
       COMMON /ITERI/ NX,NY,NZ,ITER

       real DXPA,DYPA,DZPA
       real DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &         RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &         RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADMX,RADY,RADMY,RADZ,RADMZ

       INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
       real RXFIX,RYFIX,RZFIX
       INTEGER NPALEV2

*      variable to store the relative error
       real ERR0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ERR1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /ERROR/ ERR0, ERR1

*      compressive velocity
       real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P
*      rotational velocity (ROTS variables were reused)
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

       write(*,*) 'Computing relative errors...'
*      compute the err0, err1 variables!
       CALL CELLWISE_ERROR(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ)

       write(*,*) minval(err0), maxval(err0)
       write(*,*) minval(err1), maxval(err1)

       DO IR=1,NL
        NPALEV2=MAX(100,INT(NPALEV/10))
        ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
        ALLOCATE(NVECI(NPATCH(IR)))

        DXPA=DX/(2.0**IR)
        DYPA=DY/(2.0**IR)
        DZPA=DZ/(2.0**IR)

*       build auxiliar grid for comparison
        RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
        RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
        RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

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

*         Valles: modified 23-04-2020, find if there's overlap without
*         the need of using loops

          IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &        JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &        KG1.LE.KG4.AND.KG3.LE.KG2) THEN
            NV=NV+1
            VECINO(NV,I2)=J
          END IF

         END IF       ! if J.NE.I
         END DO
         NVECI(I2)=NV

       END DO


       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) WRITE(*,*)
     &    'ERROR: gvecino second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))


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


         CORNX1=0
         CORNX2=0
         CORNXX1=0
         CORNXX2=0
         CORNY1=0
         CORNY2=0
         CORNYY1=0
         CORNYY2=0
         CORNZ1=0
         CORNZ2=0
         CORNZZ1=0
         CORNZZ2=0


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
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
              IX=II-CORNXX1+CORNX1
              JY=JJ-CORNYY1+CORNY1
              KZ=KK-CORNZZ1+CORNZ1
              IF (ERR1(IX,JY,KZ,I).LE.ERR1(II,JJ,KK,J)) THEN
                U12P(II,JJ,KK,J) = U12P(IX,JY,KZ,I)
                U13P(II,JJ,KK,J) = U13P(IX,JY,KZ,I)
                U14P(II,JJ,KK,J) = U14P(IX,JY,KZ,I)
                U12R(II,JJ,KK,J) = U12R(IX,JY,KZ,I)
                U13R(II,JJ,KK,J) = U13R(IX,JY,KZ,I)
                U14R(II,JJ,KK,J) = U14R(IX,JY,KZ,I)
                U12(II,JJ,KK,J) = U12(IX,JY,KZ,I)
                U13(II,JJ,KK,J) = U13(IX,JY,KZ,I)
                U14(II,JJ,KK,J) = U14(IX,JY,KZ,I)
              ELSE
                U12P(IX,JY,KZ,I) = U12P(II,JJ,KK,J)
                U13P(IX,JY,KZ,I) = U13P(II,JJ,KK,J)
                U14P(IX,JY,KZ,I) = U14P(II,JJ,KK,J)
                U12R(IX,JY,KZ,I) = U12R(II,JJ,KK,J)
                U13R(IX,JY,KZ,I) = U13R(II,JJ,KK,J)
                U14R(IX,JY,KZ,I) = U14R(II,JJ,KK,J)
                U12(IX,JY,KZ,I) = U12(II,JJ,KK,J)
                U13(IX,JY,KZ,I) = U13(II,JJ,KK,J)
                U14(IX,JY,KZ,I) = U14(II,JJ,KK,J)
              END IF
           END DO
           END DO
           END DO

       END DO
       END DO

       DEALLOCATE(VECINO)
       DEALLOCATE(NVECI)

      END DO !IR=1,NL

      RETURN
      END
