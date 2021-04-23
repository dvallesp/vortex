***********************************************************************
       SUBROUTINE LEER(ITER,NX,NY,NZ,T,ZETA,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************
*     Reads the outputs of the simulation. This subroutine may be
*     changed for different simulation codes. In particular, in any
*     case we'll need to feed the main code with:
****  GRIDS info:
*     NPATCH: number of patches per refinement level
*     PATCHNX, PATCHNY, PATCHNZ: cell extensions of each refinement patch
*     PATCHX, PATCHY, PATCHZ: grid coordinates of the leftmost cell of each patch
*     PATCHRX, PATCHRY, PATCHRZ: origin position of each patch (position of the leftmost "mother" cell)
*     PARE: coarser patch a given patch is embedded in
****  CLUS info:
*     U2, U3, U4: initial velocity field (base level, i.e coarse grid)
*     U12, U13, U14: initial velocity field (refinement patches)
*     CR0AMR: whether a cell is refined (=0) or it isn't (=1)
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
       INTEGER, ALLOCATABLE::SCR4_INT(:,:,:)

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       CHARACTER*9 FILNOM1,FILNOM2, FILNOM4
       CHARACTER*10 FILNOM3

       CHARACTER*24 FIL1,FIL2,FIL4
       CHARACTER*25 FIL3

       real U1(1:NMAX,1:NMAY,1:NMAZ)
       real U11(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       common /dens/ u1,u11

       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC/ U2,U3,U4,U12,U13,U14

       integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
       integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       common /cr0/ cr0amr, cr0amr1
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

       IF(SUM(NPATCH(0:NLEVELS)).GT.NPALEV) THEN
         write(*,*) "NPALEV too small!! Should at least be",
     &               SUM(NPATCH(0:NLEVELS))
       END IF

*      BARYONIC (clus file)
       READ(31)
       IR=0
       ALLOCATE(SCR4(0:NMAX+1,0:NMAY+1,0:NMAZ+1))
       SCR4=0.0
       ALLOCATE(SCR4_INT(0:NMAX+1,0:NMAY+1,0:NMAZ+1))
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
       READ(31) (((SCR4_INT(I,J,K),I=1,N1),J=1,N2),K=1,N3)
       CR0AMR(1:NX,1:NY,1:NZ)=SCR4_INT(1:NX,1:NY,1:NZ)

       DEALLOCATE(SCR4)
       DEALLOCATE(SCR4_INT)


       DO IR=1,NL
         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))
         DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          ALLOCATE(SCR4(NAMRX,NAMRY,NAMRZ))
          SCR4=0.0
         ALLOCATE(SCR4_INT(NAMRX,NAMRY,NAMRZ))
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
          READ(31) (((SCR4_INT(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
          CR0AMR1(:,:,:,I)=SCR4_INT(:,:,:)
          READ(31) !!(((SCR4_INT(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
c           SOLAP(:,:,:,I)=SCR4_INT(:,:,:)
          DEALLOCATE(SCR4)
          DEALLOCATE(SCR4_INT)
         END DO
       END DO

       CLOSE(31)

       RETURN
       END


************************************************************************
      SUBROUTINE LEE_MACH(ITER,NPATCH,PARE,PATCHNX,PATCHNY,
     &            PATCHNZ,SHOCK0,SHOCK1)
************************************************************************
*     Only used for the filter. Reads the outputs of a shock finder,
*     and identifies the shocked cells (M >= 1.3)
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     Parameters: input
      INTEGER ITER, NPATCH(0:NLEVELS), PARE(NPALEV), PATCHNX(NPALEV),
     &        PATCHNY(NPALEV), PATCHNZ(NPALEV)

*     Parameters: output
      INTEGER SHOCK0(NMAX,NMAY,NMAZ), SHOCK1(NAMRX,NAMRY,NAMRZ,NPALEV)

*     Private variables
      INTEGER I,J,K,IPATCH,N1,N2,N3
      REAL BAS, THR
      CHARACTER*20 FILNOM,FIL1
      real*4, allocatable::scr4(:,:,:)

      thr = 1.3 ! mach no. threshold (>thr --> shocked, <thr --> unshocked)

      CALL NOMFILEMACH5(ITER,FILNOM)
      FIL1='shocks/'//FILNOM
      OPEN (31,FILE=FIL1,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

      shock0 = 0
      allocate(scr4(nmax,nmay,nmaz))
      read(31) (((scr4(i,j,k),i=1,n2),j=1,n2),k=1,n3)
      n1 = nmax
      n2 = nmay
      n3 = nmaz
      do i=1,n1
        do j=1,n2
          do k=1,n3
            if (scr4(i,j,k).ge.thr) shock0(i,j,k) = 1
          end do
        end do
      end do
      deallocate(scr4)

      allocate(scr4(namrx,namry,namrz))
      do ipatch=1,sum(npatch)
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        shock1(:,:,:,ipatch) = 0
        do i=1,n1
          do j=1,n2
            do k=1,n3
              if (scr4(i,j,k).ge.thr) shock1(i,j,k,ipatch) = 1
            end do
          end do
        end do
      end do
      deallocate(scr4)

      RETURN
      END
