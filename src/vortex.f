*-------------------------------------------------------------------*
*********************************************************************
       PROGRAM VORTEX
*********************************************************************
*-------------------------------------------------------------------*
*      AUTHORS:  David Vallés-Pérez, Susana Planelles and Vicent Quilis
*      'vortex' has been developed at the Departament d'Astronomia
*      i Astrofísica of the Universitat de València, in the
*      Computational Cosmology group. This project has been supported
*      by the Spanish Ministerio de Ciencia e Innovación (MICINN,
*      grants AYA2016-77237-C3-3-P and PID2019-107427GB-C33) and by
*      the Generalitat Valenciana (grant PROMETEO/2019/071).
*-------------------------------------------------------------------*
*      'vortex' is a code which implements the Helmholtz-Hodge
*      decomposition for an AMR velocity field.
*      It has been designed to be coupled to the outputs of the
*      cosmological code MASCLET (Quilis 2004), although it can be
*      straightforwardly applied to any block-based AMR code or even
*      to particle-based outputs by means of a smoothing scheme.
*---------------------GENERAL CONSIDERATIONS------------------------*
*      Besides the source code, the following files are needed:
*      1) vortex_parameters.dat. This file dimensions the arrays,
*         and therefore the code needs to be compiled when these
*         parameters are changed.
*      2) vortex.dat. This file contains runtime parameters. They
*         can be changed once the code has been compiled.
*      3) simulation data. By default, we read the simulation data
*         in a folder simu_masclet, which contains the "gas" files.
*         In order to use vortex on other code's outputs, the
*         functions in "reader.f" (actual reader of the outputs)
*         and "nomfile.f" (names of the simulation data files) need
*         to be adapted.
*
*      The outputs of the code are written, by default, inside a
*      folder "output_files". This folder needs to be created before
*      running vortex. The output file will be saved as
*      "velocitiesXXXXX" (XXXXX is the iteration number). This
*      behaviour can be changed in "nomfile.f". As an additional
*      safety measure, the code stops if a file with the same name is
*      already in the folder.
*
*      The code is parallelised according to the OpenMP standard
*      directives. For the code to run in parallel, it has to be
*      compiled with the flag -fopenmp (gfortran), and the environment
*      variable OMP_NUM_THREADS needs to be set to the number of cores
*      set to run the code.
*********************************************************************
*-------------------------------------------------------------------*

       IMPLICIT NONE

*      COMPILATION-TIME PARAMETERS
       INCLUDE 'vortex_parameters.dat'

*      GLOBAL VARIABLES (COMMON MODULES)
       INTEGER NX,NY,NZ,ITER
       COMMON /ITERI/ NX,NY,NZ,ITER

       real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &         RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &         RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
       COMMON /GRID/  RADX,RADMX,RADY,RADMY,RADZ,RADMZ

       real DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       real U1(0:NMAX+1,0:NMAY+1,0:NMAZ+1)     ! source in poisson equation
       real POT(0:NMAX+1,0:NMAY+1,0:NMAZ+1)    ! field to solve
       COMMON /BASE/ U1,POT

       real  U11(-1:NAMRX+2,-1:NAMRY+2,-1:NAMRZ+2,NPALEV)    !source in Poisson eq.
       real  POT1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)   !field to solve
       COMMON /UAMR/ U11,POT1

       real dens0(1:NMAX,1:NMAY,1:NMAZ)
       real dens1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       common /dens/ dens0,dens1

       integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
       integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       common /cr0/ cr0amr, cr0amr1

       ! original velocities (reused afterwards)
       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC/ U2,U3,U4,U12,U13,U14

       ! compressive velocities will be saved here
       real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P

       ! to save the original velocities and reuse U2, U3, ...
       real UORI2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real UORI3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real UORI4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real UORI12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real UORI13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real UORI14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC_ORIGINAL/ UORI2,UORI3,UORI4,UORI12,UORI13,UORI14

       ! differential operators (reused as potentials afterwards)
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

       ! runtime flags
       INTEGER FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &         FLAG_W_VELOCITIES, FLAG_FILTER, FLAG_W_FILTLEN,
     &         FLAG_PARTICLES
       COMMON /FLAGS/ FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &                FLAG_W_VELOCITIES

       ! AMR grid parent cells
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

       ! SOR precision parameter and max num of iterations
       real PRECIS
       INTEGER MAXIT
       COMMON /SOR/ PRECIS,MAXIT

*      LOCAL VARIABLES
       INTEGER I,J,K,LOW1,LOW2,II,JJ,IX,JY,KZ,NL,IR,N1,N2,N3,FILT_MAXIT
       INTEGER NFILE,FIRST,EVERY,IFI,LAST,BOR,NL_PARTICLE_GRID
       real ZI,LADO,LADO0,ZETA,LIM,ERR_THR,T,FILT_TOL,FILT_STEP
       REAL OMEGA0,ACHE,FDM
       COMMON /COSMO/ OMEGA0,ACHE,FDM
       LOGICAL FILE_EXISTS
       CHARACTER*14 FILE5
       CHARACTER*30 FILERR5

       ! For the particles version
       REAL*4 RXPA(NDM),RYPA(NDM),RZPA(NDM),
     &        U2DM(NDM),U3DM(NDM),U4DM(NDM),MASAP(NDM)

       ! grids
       INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real  PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       ! base-grid sources backup (they get overwritten)
       real ROTAX_0BKP(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAY_0BKP(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAZ_0BKP(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real DIVER0BKP(0:NMAX+1,0:NMAY+1,0:NMAZ+1)

       real KKK(NMAX,NMAY,NMAZ)    !KKK coeficients of Fourier series

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM


****************************************************
*      READING INITIAL DATA                        *
****************************************************
       OPEN(1,FILE='vortex.dat',STATUS='UNKNOWN',ACTION='READ')

       READ(1,*)
       READ(1,*)
       READ(1,*) FIRST,LAST,EVERY
       READ(1,*)
       READ(1,*) NX,NY,NZ
       READ(1,*)
       READ(1,*) ZI,LADO0
       READ(1,*)
       READ(1,*) OMEGA0,ACHE,FDM
       READ(1,*)
       READ(1,*) NL
       READ(1,*)
       READ(1,*) PRECIS, MAXIT, BOR
       READ(1,*)
       READ(1,*) ERR_THR
       READ(1,*)
       READ(1,*) FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &           FLAG_W_VELOCITIES
       READ(1,*)
       READ(1,*) FLAG_FILTER, FLAG_W_FILTLEN, FILT_TOL, FILT_STEP,
     &           FILT_MAXIT
       READ(1,*)
       READ(1,*) FLAG_PARTICLES, NL_PARTICLE_GRID

       CLOSE(1)


**************************************************************
*     ...PARALLEL RUNNING...
*     NUMBER OF PROCESSORS
      NUM=4
!$OMP PARALLEL SHARED(NUM)
!$OMP SINGLE
!$      NUM=OMP_GET_NUM_THREADS()
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
**************************************************************


* ===========  this is global for a given output ============================
*     GRID BUILDER
      LADO=LADO0-(LADO0/NX) ! from leftmost center to rightmost center

*     coarse grid:
      CALL MALLA(NX,NY,NZ,LADO)

*     KKK coeficients of Fourier series (coarse grid)
      CALL MOMENTO(DX,NX,NY,NZ,KKK) ! once in the code
*============================================================================


       NFILE=INT((LAST-FIRST)/EVERY) + 1
       WRITE(*,*) 'NFILE=',NFILE

* === first, we check that no output files will be overwritten
      DO IFI=1,NFILE
        ITER=FIRST+EVERY*(IFI-1)
        CALL NOMFILE2(ITER,FILE5)
        FILERR5='./output_files/'//FILE5
        INQUIRE(FILE=FILERR5,EXIST=FILE_EXISTS)
        IF (FILE_EXISTS) THEN
          WRITE(*,*) 'The file ', FILERR5, 'already exists'
          STOP 'Program will terminate'
        END IF
      END DO
* === end of check ===========================================


*////////////////////////////////////
       DO IFI=1,NFILE
*////////////////////////////////////

       ITER=FIRST+EVERY*(IFI-1)

*      define the output file name
       CALL NOMFILE2(ITER,FILE5)
       FILERR5='./output_files/'//FILE5

* ===========  READ DATA FROM THE SIMULATION ============================

       NX=NMAX
       NY=NMAY
       NZ=NMAZ

*      READING DATA
       IF (FLAG_PARTICLES.EQ.0) THEN
        CALL LEER(ITER,NX,NY,NZ,T,ZETA,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
       ELSE
        ! HERE WE WILL CALL THE LEER_PARTICLES(), CREATE_MESH(),
        ! AND INTERPOLATE_PARTICLES() ROUTINES.
        CALL LEER_PARTICLES(ITER,NX,NY,NZ,T,ZETA,NL,NL_PARTICLE_GRID,
     &            NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,
     &            PATCHZ,PATCHRX,PATCHRY,PATCHRZ,LADO0,NPART,RXPA,
     &            RYPA,RZPA,MASAP,U2DM,U3DM,U4DM)
       END IF


*      INITIALIZE VARIABLES TO ZERO
*      BASE LEVEL

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2P,U3P,U4P,U1,POT,DIVER0,
!$OMP+                   ROTAX_0,ROTAY_0,ROTAZ_0),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
        U2P(I,J,K)=0.0
        U3P(I,J,K)=0.0
        U4P(I,J,K)=0.0

        U1(I,J,K)=0.0
        POT(I,J,K)=0.0
        DIVER0(I,J,K)=0.0
        ROTAX_0(I,J,K)=0.0
        ROTAY_0(I,J,K)=0.0
        ROTAZ_0(I,J,K)=0.0
      END DO
      END DO
      END DO

*     REFINEMENT LEVELS
      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U11),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2

          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          DO KZ=1,N3
          DO JY=1,N2
          DO IX=1,N1
             U11(IX,JY,KZ,I)=0.0
          END DO
          END DO
          END DO
      END DO
      END DO


      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U12,U13,U14,U12P,U13P,U14P,POT1,DIVER,
!$OMP+                   ROTAX_1,ROTAY_1,ROTAZ_1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
        DO KZ=0,N3+1
        DO JY=0,N2+1
        DO IX=0,N1+1
          U12P(IX,JY,KZ,I)=0.0
          U13P(IX,JY,KZ,I)=0.0
          U14P(IX,JY,KZ,I)=0.0

          POT1(IX,JY,KZ,I)=0.0
          DIVER(IX,JY,KZ,I)=0.0
          ROTAX_1(IX,JY,KZ,I)=0.0
          ROTAY_1(IX,JY,KZ,I)=0.0
          ROTAZ_1(IX,JY,KZ,I)=0.0
       END DO
       END DO
       END DO
      END DO
      END DO

       IF (ZETA.LT.0.0) ZETA=0.0


*     Build the AMR grid:
      CALL GRIDAMR(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &      PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,PARE)

*     Runtime check: velocities have been read properly
      IF (FLAG_VERBOSE.EQ.1) THEN
        write(*,*) 'velocity: min and max values'
        call p_minmax(u2,u12,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                npatch)
        call p_minmax(u3,u13,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                npatch)
        call p_minmax(u4,u14,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                npatch)
      END IF

*     Filter velocities (if specified to do so in vortex.dat)
      IF (flag_filter.eq.1) then
        IF (flag_verbose.eq.1) write(*,*) 'Applying multiscale filter'
        call MULTISCALE_FILTER(NX,NY,NZ,NL,NPATCH,pare,
     &            PATCHNX,PATCHNY,PATCHNZ,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,DX,ITER,FLAG_W_FILTLEN,
     &            FILT_TOL,FILT_STEP, FILT_MAXIT)
        IF (FLAG_VERBOSE.EQ.1) THEN
         write(*,*) 'Computation ended!'
         write(*,*) 'filtered velocity: min and max values'
         call p_minmax(u2,u12,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                npatch)
         call p_minmax(u3,u13,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                npatch)
         call p_minmax(u4,u14,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                npatch)
        END IF
      END IF

*     All patches are extended with one extra cell per direction
      CALL EXTEND_VAR(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)

*------------------------------------------------------------------*
*      We compute the velocity divergence and curl
*------------------------------------------------------------------*

        WRITE(*,*) 'Computing the velocity rotational...'

        CALL ROTARY(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &              PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)

        WRITE(*,*) 'Computation ends!'

        WRITE(*,*) 'Computing velocity divergence...'

        CALL DIVER_FINA(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,
     &                  PATCHNY,PATCHNZ,PATCHX,PATCHY,
     &                  PATCHZ,PATCHRX,PATCHRY,PATCHRZ)

        WRITE(*,*) 'Computation ends!'

        IF (FLAG_VERBOSE.EQ.1) THEN
          write(*,*) 'rotational: min and max values'
          call p_minmax(rotax_0,rotax_1,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(rotay_0,rotay_1,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(rotaz_0,rotaz_1,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)

          write(*,*) 'divergence: min and max values'
          call p_minmax(diver0,diver,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)

        END IF

        ! we back-up the sources for later use (they will get overwritten
        ! once the base grid is solved)
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DIVER0,ROTAX_0,ROTAY_0,ROTAZ_0,
!$OMP+                   DIVER0BKP,ROTAX_0BKP,ROTAY_0BKP,ROTAZ_0BKP),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
        DO K=1,NZ
        DO J=1,NY
        DO I=1,NX
         DIVER0BKP(I,J,K)=DIVER0(I,J,K)
         ROTAX_0BKP(I,J,K)=ROTAX_0(I,J,K)
         ROTAY_0BKP(I,J,K)=ROTAY_0(I,J,K)
         ROTAZ_0BKP(I,J,K)=ROTAZ_0(I,J,K)
        END DO
        END DO
        END DO

        IF (FLAG_W_DIVROT.EQ.1) THEN
          CALL WRITE_DIVROT(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,NPATCH,
     &                      PATCHNX,PATCHNY,PATCHNZ)
        END IF


* >>>>>>>>>>>>>    THIS IS FOR EACH FIELD TO BE SOLVED <<<<<<<<<<<<<<<<<<<<<<<<<<<<
*      COARSE GRID:
*
*      NABLA(FIELD)=SOURCE
*      SOURCE: SOURCE OF POISSON EQUATION, ARRAY(0:NX+1,0:NY+1,0:NZ+1)
*      FIELD: SOLUTION OF POISSON EQUATION, ARRAY(0:NX+1,0:NY+1,0:NZ+1)
*      ONE FICTICIOUS CELL IN EACH DIRECTION ASSUMING PERIODIC BOUNDARY CONDITIONS
*      real POT(0:NMAX+1,0:NMAY+1,0:NMAZ+1)    ! field to solve
*      real U1(0:NMAX+1,0:NMAY+1,0:NMAZ+1)     ! source in poisson equation
*      COMMON /BASE/ U1,POT
*
*      WE 'MUST' SOLVE 4 DIFFERENT POISSON EQUATIONS --> 4 CALLS:
*       --> IR=0, i.e. base level: we call to POFFT3D
*       --> IR>0, i.e. refinement patches: we call to POTAMR
*
*      To save memory, we then copy the potentials to the source arrays.
*      So after this calculations, DIVER0/DIVER will contain the scalar
*      potential, and ROTA(X,Y,Z)_0/1 will contain the vector potential.

      WRITE(*,*) 'Solving Poisson eqns. Base level'
      WRITE(*,*) 'Scalar potential'

**     SOURCE=-DIVER0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,DIVER0),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
          U1(I,J,K)=-1.0*DIVER0(I,J,K)
      END DO
      END DO
      END DO

      CALL POFFT3D(NX,NY,NZ,KKK)    ! returns field POT --> PHI

!$OMP PARALLEL DO SHARED(NX,NY,NZ,DIVER0,POT),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
         DIVER0(I,J,K)=POT(I,J,K)
      END DO
      END DO
      END DO

      WRITE(*,*) 'Vector potential: x component'
**     SOURCE=-ROTAX_0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,ROTAX_0),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
        U1(I,J,K)=-1.0*ROTAX_0(I,J,K)
      END DO
      END DO
      END DO

      CALL POFFT3D(NX,NY,NZ,KKK)    ! returns field POT --> W_x

!$OMP PARALLEL DO SHARED(NX,NY,NZ,ROTAX_0,POT),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
         ROTAX_0(I,J,K)=POT(I,J,K)
      END DO
      END DO
      END DO

      WRITE(*,*) 'Vector potential: y component'
**     SOURCE=-ROTAY_0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,ROTAY_0),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
        U1(I,J,K)=-1.0*ROTAY_0(I,J,K)
      END DO
      END DO
      END DO

      CALL POFFT3D(NX,NY,NZ,KKK)    ! returns field POT --> W_y

!$OMP PARALLEL DO SHARED(NX,NY,NZ,ROTAY_0,POT),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
        ROTAY_0(I,J,K)=POT(I,J,K)
      END DO
      END DO
      END DO

      WRITE(*,*) 'Vector potential: z component'
**     SOURCE=-ROTAZ_0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,ROTAZ_0),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
         U1(I,J,K)=-1.0*ROTAZ_0(I,J,K)
      END DO
      END DO
      END DO

      CALL POFFT3D(NX,NY,NZ,KKK)    ! returns field POT --> W_z

!$OMP PARALLEL DO SHARED(NX,NY,NZ,ROTAZ_0,POT),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
        ROTAZ_0(I,J,K)=POT(I,J,K)
      END DO
      END DO
      END DO

      WRITE(*,*) 'Solving Poisson eqns. AMR levels'
*      AMR levels:
*
*      Source and field are too large, must be sent by COMMON
*      real  U11(NAMRX,NAMRY,NAMRZ,NPALEV)
*      real  POT1(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
*      COMMON /UAMR/U11,POT1
*
*     WARNING: POTAMR WORKS ON THE WHOLE GRID!!!
*

** 1ST CALL:
      WRITE(*,*) 'Scalar potential'
!$OMP PARALLEL DO SHARED(NX,NY,NZ,POT,DIVER0,U1,DIVER0BKP),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
          !Initial guess...
          POT(I,J,K)=DIVER0(I,J,K)
          U1(I,J,K)=-1.0*DIVER0BKP(I,J,K)
      END DO
      END DO
      END DO

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U11,DIVER),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1
*      SOURCE=-DIVER
         U11(IX,JY,KZ,I)=-1.0*DIVER(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

       CALL POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,
     &             PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,BOR)

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   POT1,DIVER),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=-2, N3+3
       DO JY=-2, N2+3
       DO IX=-2, N1+3
          DIVER(IX,JY,KZ,I)=POT1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

** 2ND CALL:
       WRITE(*,*) 'Vector potential: x component'
!$OMP PARALLEL DO SHARED(NX,NY,NZ,POT,ROTAX_0,U1,ROTAX_0BKP),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
       DO K=0, NZ+1
       DO J=0, NY+1
       DO I=0, NX+1
            !Initial guess...
            POT(I,J,K)=ROTAX_0(I,J,K)
            U1(I,J,K)=-1.0*ROTAX_0BKP(I,J,K)
       END DO
       END DO
       END DO

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U11,ROTAX_1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1
*      SOURCE=-ROTAX_1
          U11(IX,JY,KZ,I)=-1.0*ROTAX_1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

       CALL POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,
     &             PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,BOR)

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   ROTAX_1,POT1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=-2, N3+3
       DO JY=-2, N2+3
       DO IX=-2, N1+3
         ROTAX_1(IX,JY,KZ,I)=POT1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

* 3RD CALL:
       WRITE(*,*) 'Vector potential: y component'
!$OMP PARALLEL DO SHARED(NX,NY,NZ,POT,ROTAY_0,U1,ROTAY_0BKP),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
       DO K=0, NZ+1
       DO J=0, NY+1
       DO I=0, NX+1
!Initial guess...
        POT(I,J,K)=ROTAY_0(I,J,K)
        U1(I,J,K)=-1.0*ROTAY_0BKP(I,J,K)
       END DO
       END DO
       END DO

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U11,ROTAY_1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1
*      SOURCE=-ROTAY_1
          U11(IX,JY,KZ,I)=-1.0*ROTAY_1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

       CALL POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,
     &             PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,BOR)


       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   ROTAY_1,POT1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=-2, N3+3
       DO JY=-2, N2+3
       DO IX=-2, N1+3
         ROTAY_1(IX,JY,KZ,I)=POT1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

* 4TH CALL:
       WRITE(*,*) 'Vector potential: z component'
!$OMP PARALLEL DO SHARED(NX,NY,NZ,POT,ROTAZ_0,U1,ROTAZ_0BKP),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
       DO K=0, NZ+1
       DO J=0, NY+1
       DO I=0, NX+1
!Initial guess...
         POT(I,J,K)=ROTAZ_0(I,J,K)
         U1(I,J,K)=-1.0*ROTAZ_0BKP(I,J,K)
       END DO
       END DO
       END DO

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U11,ROTAZ_1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1
*      SOURCE=-ROTAZ_1
          U11(IX,JY,KZ,I)=-1.0*ROTAZ_1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

       CALL POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,
     &             PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,BOR)

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   ROTAZ_1,POT1),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       DO KZ=-2, N3+3
       DO JY=-2, N2+3
       DO IX=-2, N1+3
          ROTAZ_1(IX,JY,KZ,I)=POT1(IX,JY,KZ,I)
       END DO
       END DO
       END DO

       END DO
       END DO

       WRITE(*,*) 'Computation ended!'

**** WARNING: NOW DIVER AND (ROTAX,ROTAY,ROTAZ) ARE SOLUTIONS OF POISSON EQ.
       IF (FLAG_W_POTENTIALS.EQ.1) THEN
         CALL WRITE_POTENTIALS(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,NPATCH,
     &                         PATCHNX,PATCHNY,PATCHNZ)
       END IF

**** ---> WE NEED TO COMPUTE: -GRAD(DIVER) AND ROT(ROTAX,ROTAY,ROTAZ)

        IF (FLAG_VERBOSE.EQ.1) THEN
          WRITE(*,*) '...Total velocity...'
          call p_minmax(u2,u12,1,1,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(u3,u13,1,1,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(u4,u14,1,1,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
        END IF

*       We backup the original velocities in UORI

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2,U3,U4,UORI2,UORI3,UORI4),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
       DO K=1, NZ
       DO J=1, NY
       DO I=1, NX
         UORI2(I,J,K)=U2(I,J,K)
         UORI3(I,J,K)=U3(I,J,K)
         UORI4(I,J,K)=U4(I,J,K)
       END DO
       END DO
       END DO

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U12,U13,U14,UORI12,UORI13,UORI14),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1, N3
        DO JY=1, N2
        DO IX=1, N1
           UORI12(IX,JY,KZ,I)=U12(IX,JY,KZ,I)
           UORI13(IX,JY,KZ,I)=U13(IX,JY,KZ,I)
           UORI14(IX,JY,KZ,I)=U14(IX,JY,KZ,I)
        END DO
        END DO
        END DO

        END DO
        END DO
*       END backuping original velocities

        WRITE(*,*) '...Differencing the potentials...'
*     We compute the -grad(PHI)  ---> we get (U2P,U3P,U4P)
*     PHI NOW IS IN DIVER0 AND DIVER!!!
       CALL GRADIENTE(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)

**       !!!! Por convenio hay que multiplicar los UP por -1

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2P,U3P,U4P),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
       DO K=1, NZ
       DO J=1, NY
       DO I=1, NX
         U2P(I,J,K)=-1.0*U2P(I,J,K)
         U3P(I,J,K)=-1.0*U3P(I,J,K)
         U4P(I,J,K)=-1.0*U4P(I,J,K)
       END DO
       END DO
       END DO

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,
!$OMP+                   U12P,U13P,U14P),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1, N3
        DO JY=1, N2
        DO IX=1, N1
           U12P(IX,JY,KZ,I)=-1.0*U12P(IX,JY,KZ,I)
           U13P(IX,JY,KZ,I)=-1.0*U13P(IX,JY,KZ,I)
           U14P(IX,JY,KZ,I)=-1.0*U14P(IX,JY,KZ,I)
        END DO
        END DO
        END DO

        END DO
        END DO

**     We compute the rotational of (rotax,rotay,rotaz) ---> we get (U2R,U3R,U4R)
*      Note that we lose the original velocities (U2, U3, U4), as we overwrite them

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2,U3,U4,ROTAX_0,ROTAY_0,ROTAZ_0),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=0, NZ+1
      DO J=0, NY+1
      DO I=0, NX+1
       U2(I,J,K)=ROTAX_0(I,J,K)
       U3(I,J,K)=ROTAY_0(I,J,K)
       U4(I,J,K)=ROTAZ_0(I,J,K)
      END DO
      END DO
      END DO

       DO IR=1,NL
         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))
         DO I=LOW1,LOW2
           N1=PATCHNX(I)
           N2=PATCHNY(I)
           N3=PATCHNZ(I)
           U12(0:N1+1,0:N2+1,0:N3+1,I)=ROTAX_1(0:N1+1,0:N2+1,0:N3+1,I)
           U13(0:N1+1,0:N2+1,0:N3+1,I)=ROTAY_1(0:N1+1,0:N2+1,0:N3+1,I)
           U14(0:N1+1,0:N2+1,0:N3+1,I)=ROTAZ_1(0:N1+1,0:N2+1,0:N3+1,I)
        END DO
       END DO

       CALL ROTARY_2(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &             PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)



*      Ensure that overlapping cells have the same value of the velocity
*      after the calculations
       CALL SYNC_AMR_VELOCITIES(NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                       PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                       PATCHRZ)

*      Outliers (cells where differentiation noise has produced large
*      relative errors) get corrected by interpolation from coarser grids
       CALL CORRECT_OUTLIERS(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                       ERR_THR)

        IF (FLAG_VERBOSE.EQ.1) THEN
          WRITE(*,*) '...Compressional velocity...'
          call p_minmax(u2p,u12p,1,1,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(u3p,u13p,1,1,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(u4p,u14p,1,1,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          WRITE(*,*) '...Rotational velocity...'
          call p_minmax(rotax_0,rotax_1,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(rotay_0,rotay_1,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
          call p_minmax(rotaz_0,rotaz_1,1,3,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch)
        END IF

        WRITE(*,*) 'Computation ended!'

        IF (FLAG_W_VELOCITIES.EQ.1) THEN
          CALL WRITE_VELOCITIES(FILERR5,NX,NY,NZ,ITER,T,ZETA,NL,
     &                          NPATCH, PATCHNX,PATCHNY,PATCHNZ)
          IF (FLAG_PARTICLES.EQ.1) THEN
           CALL WRITE_PARTICLES(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                          PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &                          MASAP,U2DM,U3DM,U4DM,NPART,LADO0)
          END IF
        END IF

*//////////////////////////////////// ! DO IFI=1,NFILE
       END DO
*////////////////////////////////////

*********************************************************************
*********************************************************************
       END
*********************************************************************
*********************************************************************

      !!!! FUNCTIONS in EXTERNAL FILES
*     Differential operators
      INCLUDE 'diff_ho.f' ! diff.f for first order
*     Filenames
      INCLUDE 'nomfile.f'
*     Build the base and AMR grids
      INCLUDE 'grids.f'
*     Linear interpolation routines
      INCLUDE 'interp.f'
*     Solve elliptic equations at base and refined levels
      INCLUDE 'poisson.f'
*     Read the input data
      INCLUDE 'reader.f'
*     Write the outputs
      INCLUDE 'writer.f'
*     Handle the overlaps between different patches
      INCLUDE 'overlaps.f'
*     Handle the boundaries of AMR patches
      INCLUDE 'boundaries.f'
*     Detect cells with large errors and interpolate from coarser levels
      INCLUDE 'outliers.f'
*     Multiscale filter as in (Vazza, 2012) to extract turbulent field
      INCLUDE 'filter.f'
*     Routines for working with particles
      INCLUDE 'particles.f'
*     Routines from 'Numerical Recipes in Fortran90', Press, Teukoslky et al.
      INCLUDE 'nr.f'
