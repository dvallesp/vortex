***********************************************************************
       SUBROUTINE ROTARY(NX,NY,NZ,NL,NPATCH,
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
       real  PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       real DXPA,DYPA,DZPA,XXX1,YYY1,ZZZ1
       real bas23,bas34,bas24,bas32,bas43,bas42,AAA,BBB,CCC
       INTEGER CR1,CR2,CR3
       INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE

       real coef4(-4:4), coef3(-3:3), coef2(-2:2), coef1(-1:1)
       integer idx

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

*-------------------------------------
*      Divergencia fina  (DIVER)
*            (IR: 1 ---> NL)
*            (celdas: 2 ----> NX-1) (excluimos los bordes)
*-------------------------------------

       coef4 = (/1.0/280,-4.0/105,1.0/5,-4.0/5,0.0,4.0/5,-1.0/5,4.0/105,
     &           -1.0/280 /)
       coef3 = (/-1.0/60,3.0/20,-3.0/4,0.0,3.0/4,-3.0/20,1.0/60/)
       coef2 = (/1.0/12,-2.0/3,0.0,2.0/3,-1.0/12/)
       coef1 = (/-1.0/2,0.0,1.0/2/)


       DO IR=1,NL

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

* X
        bas32 = 0.0
        bas42 = 0.0
        if (ix.ge.4.and.ix.le.n1-3) then
         do idx=-4,4
           bas32 = bas32 + coef4(idx) * u13(ix+idx,jy,kz,i)
           bas42 = bas42 + coef4(idx) * u14(ix+idx,jy,kz,i)
         end do
        else if (ix.ge.3.and.ix.le.n1-2) then
         do idx=-3,3
           bas32 = bas32 + coef3(idx) * u13(ix+idx,jy,kz,i)
           bas42 = bas42 + coef3(idx) * u14(ix+idx,jy,kz,i)
         end do
        else if (ix.ge.2.and.ix.le.n1-1) then
         do idx=-2,2
           bas32 = bas32 + coef2(idx) * u13(ix+idx,jy,kz,i)
           bas42 = bas42 + coef2(idx) * u14(ix+idx,jy,kz,i)
         end do
        else
         do idx=-1,1
           bas32 = bas32 + coef1(idx) * u13(ix+idx,jy,kz,i)
           bas42 = bas42 + coef1(idx) * u14(ix+idx,jy,kz,i)
         end do
        end if
        bas32 = bas32/dxpa
        bas42 = bas42/dxpa

* Y
        bas23 = 0.0
        bas43 = 0.0
        if (jy.ge.4.and.jy.le.n2-3) then
         do idx=-4,4
           bas23 = bas23 + coef4(idx) * u12(ix,jy+idx,kz,i)
           bas43 = bas43 + coef4(idx) * u14(ix,jy+idx,kz,i)
         end do
        else if (jy.ge.3.and.jy.le.n2-2) then
         do idx=-3,3
           bas23 = bas23 + coef3(idx) * u12(ix,jy+idx,kz,i)
           bas43 = bas43 + coef3(idx) * u14(ix,jy+idx,kz,i)
         end do
        else if (jy.ge.2.and.jy.le.n2-1) then
         do idx=-2,2
           bas23 = bas23 + coef2(idx) * u12(ix,jy+idx,kz,i)
           bas43 = bas43 + coef2(idx) * u14(ix,jy+idx,kz,i)
         end do
        else
         do idx=-1,1
           bas23 = bas23 + coef1(idx) * u12(ix,jy+idx,kz,i)
           bas43 = bas43 + coef1(idx) * u14(ix,jy+idx,kz,i)
         end do
        end if
        bas23 = bas23/dypa
        bas43 = bas43/dypa

* Z
        bas24 = 0.0
        bas34 = 0.0
        if (kz.ge.4.and.kz.le.n3-3) then
         do idx=-4,4
           bas24 = bas24 + coef4(idx) * u12(ix,jy,kz+idx,i)
           bas34 = bas34 + coef4(idx) * u13(ix,jy,kz+idx,i)
         end do
        else if (kz.ge.3.and.kz.le.n3-2) then
         do idx=-3,3
           bas24 = bas24 + coef3(idx) * u12(ix,jy,kz+idx,i)
           bas34 = bas34 + coef3(idx) * u13(ix,jy,kz+idx,i)
         end do
        else if (kz.ge.2.and.kz.le.n3-1) then
         do idx=-2,2
           bas24 = bas24 + coef2(idx) * u12(ix,jy,kz+idx,i)
           bas34 = bas34 + coef2(idx) * u13(ix,jy,kz+idx,i)
         end do
        else
         do idx=-1,1
           bas24 = bas24 + coef1(idx) * u12(ix,jy,kz+idx,i)
           bas34 = bas34 + coef1(idx) * u13(ix,jy,kz+idx,i)
         end do
        end if
        bas24 = bas24/dzpa
        bas34 = bas34/dzpa

        ROTAX_1(IX,JY,KZ,I)=bas43-bas34
        ROTAY_1(IX,JY,KZ,I)=bas24-bas42
        ROTAZ_1(IX,JY,KZ,I)=bas32-bas23

       END DO
       END DO
       END DO


       END DO
       END DO

*-------------------------------*
*      COARSE LEVEL
*-------------------------------*


       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX

* X
       bas32 = 0.0
       bas42 = 0.0
       if (ix.ge.4.and.ix.le.nx-4) then
        do idx=-4,4
          bas32 = bas32 + coef4(idx) * u3(ix+idx,jy,kz)
          bas42 = bas42 + coef4(idx) * u4(ix+idx,jy,kz)
        end do
       else if (ix.ge.3.and.ix.le.nx-3) then
        do idx=-3,3
          bas32 = bas32 + coef3(idx) * u3(ix+idx,jy,kz)
          bas42 = bas42 + coef3(idx) * u4(ix+idx,jy,kz)
        end do
       else if (ix.ge.2.and.ix.le.nx-2) then
        do idx=-2,2
          bas32 = bas32 + coef2(idx) * u3(ix+idx,jy,kz)
          bas42 = bas42 + coef2(idx) * u4(ix+idx,jy,kz)
        end do
       else
        do idx=-1,1
          bas32 = bas32 + coef1(idx) * u3(ix+idx,jy,kz)
          bas42 = bas42 + coef1(idx) * u4(ix+idx,jy,kz)
        end do
       end if
       bas32 = bas32/dx
       bas42 = bas42/dx

* Y
       bas23 = 0.0
       bas43 = 0.0
       if (jy.ge.4.and.jy.le.ny-4) then
        do idx=-4,4
          bas23 = bas23 + coef4(idx) * u2(ix,jy+idx,kz)
          bas43 = bas43 + coef4(idx) * u4(ix,jy+idx,kz)
        end do
       else if (jy.ge.3.and.jy.le.ny-3) then
        do idx=-3,3
          bas23 = bas23 + coef3(idx) * u2(ix,jy+idx,kz)
          bas43 = bas43 + coef3(idx) * u4(ix,jy+idx,kz)
        end do
       else if (jy.ge.2.and.jy.le.ny-2) then
        do idx=-2,2
          bas23 = bas23 + coef2(idx) * u2(ix,jy+idx,kz)
          bas43 = bas43 + coef2(idx) * u4(ix,jy+idx,kz)
        end do
       else
        do idx=-1,1
          bas23 = bas23 + coef1(idx) * u2(ix,jy+idx,kz)
          bas43 = bas43 + coef1(idx) * u4(ix,jy+idx,kz)
        end do
       end if
       bas23 = bas23/dy
       bas43 = bas43/dy

* Z
       bas24 = 0.0
       bas34 = 0.0
       if (kz.ge.4.and.kz.le.nz-4) then
        do idx=-4,4
          bas24 = bas24 + coef4(idx) * u2(ix,jy,kz+idx)
          bas34 = bas34 + coef4(idx) * u3(ix,jy,kz+idx)
        end do
      else if (kz.ge.3.and.kz.le.nz-3) then
        do idx=-3,3
          bas24 = bas24 + coef3(idx) * u2(ix,jy,kz+idx)
          bas34 = bas34 + coef3(idx) * u3(ix,jy,kz+idx)
        end do
      else if (kz.ge.2.and.kz.le.nz-2) then
        do idx=-2,2
          bas24 = bas24 + coef2(idx) * u2(ix,jy,kz+idx)
          bas34 = bas34 + coef2(idx) * u3(ix,jy,kz+idx)
        end do
       else
        do idx=-1,1
          bas24 = bas24 + coef1(idx) * u2(ix,jy,kz+idx)
          bas34 = bas34 + coef1(idx) * u3(ix,jy,kz+idx)
        end do
       end if
       bas24 = bas24/dz
       bas34 = bas34/dz

       ROTAX_0(IX,JY,KZ)=bas43-bas34
       ROTAY_0(IX,JY,KZ)=bas24-bas42
       ROTAZ_0(IX,JY,KZ)=bas32-bas23

       END DO
       END DO
       END DO


        RETURN
       END

***********************************************************************
      SUBROUTINE ROTARY_2(NX,NY,NZ,NL,NPATCH,
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
      real U12(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real U13(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real U14(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
*      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

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
      real bas23,bas34,bas24,bas32,bas43,bas42,AAA,BBB,CCC
      INTEGER CR1,CR2,CR3
      INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE

      real coef4(-4:4), coef3(-3:3), coef2(-2:2), coef1(-1:1)
      integer idx

*      ---PARALLEL---
      INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
      COMMON /PROCESADORES/ NUM

*-------------------------------------
*      Divergencia fina  (DIVER)
*            (IR: 1 ---> NL)
*            (celdas: 2 ----> NX-1) (excluimos los bordes)
*-------------------------------------

      coef4 = (/1.0/280,-4.0/105,1.0/5,-4.0/5,0.0,4.0/5,-1.0/5,4.0/105,
     &           -1.0/280 /)
      coef3 = (/-1.0/60,3.0/20,-3.0/4,0.0,3.0/4,-3.0/20,1.0/60/)
      coef2 = (/1.0/12,-2.0/3,0.0,2.0/3,-1.0/12/)
      coef1 = (/-1.0/2,0.0,1.0/2/)


      DO IR=1,NL

      DXPA=0.0
      DYPA=0.0
      DZPA=0.0
      DXPA=DX/(2.0**IR)
      DYPA=DY/(2.0**IR)
      DZPA=DZ/(2.0**IR)

      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
      DO I=LOW1,LOW2

      N1=PATCHNX(I)
      N2=PATCHNY(I)
      N3=PATCHNZ(I)

      U12(-2:N1+3,-2:N2+3,-2:N3+3,I) =
     &        ROTAX_1(-2:N1+3,-2:N2+3,-2:N3+3,I)
      U13(-2:N1+3,-2:N2+3,-2:N3+3,I) =
     &        ROTAY_1(-2:N1+3,-2:N2+3,-2:N3+3,I)
      U14(-2:N1+3,-2:N2+3,-2:N3+3,I) =
     &        ROTAZ_1(-2:N1+3,-2:N2+3,-2:N3+3,I)

      DO KZ=1, N3
      DO JY=1, N2
      DO IX=1, N1

* X
       bas32 = 0.0
       bas42 = 0.0
       if (ix.ge.2.and.ix.le.n1-1) then
        do idx=-4,4
          bas32 = bas32 + coef4(idx) * u13(ix+idx,jy,kz,i)
          bas42 = bas42 + coef4(idx) * u14(ix+idx,jy,kz,i)
        end do
       else if (ix.ge.1.and.ix.le.n1) then
        do idx=-3,3
          bas32 = bas32 + coef3(idx) * u13(ix+idx,jy,kz,i)
          bas42 = bas42 + coef3(idx) * u14(ix+idx,jy,kz,i)
        end do
       end if
       bas32 = bas32/dxpa
       bas42 = bas42/dxpa

* Y
       bas23 = 0.0
       bas43 = 0.0
       if (jy.ge.2.and.jy.le.n2-1) then
        do idx=-4,4
          bas23 = bas23 + coef4(idx) * u12(ix,jy+idx,kz,i)
          bas43 = bas43 + coef4(idx) * u14(ix,jy+idx,kz,i)
        end do
       else if (jy.ge.1.and.jy.le.n2) then
        do idx=-3,3
          bas23 = bas23 + coef3(idx) * u12(ix,jy+idx,kz,i)
          bas43 = bas43 + coef3(idx) * u14(ix,jy+idx,kz,i)
        end do
       end if
       bas23 = bas23/dypa
       bas43 = bas43/dypa

* Z
       bas24 = 0.0
       bas34 = 0.0
       if (kz.ge.2.and.kz.le.n3-1) then
        do idx=-4,4
          bas24 = bas24 + coef4(idx) * u12(ix,jy,kz+idx,i)
          bas34 = bas34 + coef4(idx) * u13(ix,jy,kz+idx,i)
        end do
       else if (kz.ge.1.and.kz.le.n3) then
        do idx=-3,3
          bas24 = bas24 + coef3(idx) * u12(ix,jy,kz+idx,i)
          bas34 = bas34 + coef3(idx) * u13(ix,jy,kz+idx,i)
        end do
       end if
       bas24 = bas24/dzpa
       bas34 = bas34/dzpa

       ROTAX_1(IX,JY,KZ,I)=bas43-bas34
       ROTAY_1(IX,JY,KZ,I)=bas24-bas42
       ROTAZ_1(IX,JY,KZ,I)=bas32-bas23

      END DO
      END DO
      END DO


      END DO
      END DO

*-------------------------------*
*      COARSE LEVEL
*-------------------------------*

      U2(0:NX+1,0:NY+1,0:NZ+1) = ROTAX_0(0:NX+1,0:NY+1,0:NZ+1)
      U3(0:NX+1,0:NY+1,0:NZ+1) = ROTAY_0(0:NX+1,0:NY+1,0:NZ+1)
      U4(0:NX+1,0:NY+1,0:NZ+1) = ROTAZ_0(0:NX+1,0:NY+1,0:NZ+1)

      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX

* X
      bas32 = 0.0
      bas42 = 0.0
      if (ix.ge.4.and.ix.le.nx-4) then
       do idx=-4,4
         bas32 = bas32 + coef4(idx) * u3(ix+idx,jy,kz)
         bas42 = bas42 + coef4(idx) * u4(ix+idx,jy,kz)
       end do
      else if (ix.ge.3.and.ix.le.nx-3) then
       do idx=-3,3
         bas32 = bas32 + coef3(idx) * u3(ix+idx,jy,kz)
         bas42 = bas42 + coef3(idx) * u4(ix+idx,jy,kz)
       end do
      else if (ix.ge.2.and.ix.le.nx-2) then
       do idx=-2,2
         bas32 = bas32 + coef2(idx) * u3(ix+idx,jy,kz)
         bas42 = bas42 + coef2(idx) * u4(ix+idx,jy,kz)
       end do
      else
       do idx=-1,1
         bas32 = bas32 + coef1(idx) * u3(ix+idx,jy,kz)
         bas42 = bas42 + coef1(idx) * u4(ix+idx,jy,kz)
       end do
      end if
      bas32 = bas32/dx
      bas42 = bas42/dx

* Y
      bas23 = 0.0
      bas43 = 0.0
      if (jy.ge.4.and.jy.le.ny-4) then
       do idx=-4,4
         bas23 = bas23 + coef4(idx) * u2(ix,jy+idx,kz)
         bas43 = bas43 + coef4(idx) * u4(ix,jy+idx,kz)
       end do
      else if (jy.ge.3.and.jy.le.ny-3) then
       do idx=-3,3
         bas23 = bas23 + coef3(idx) * u2(ix,jy+idx,kz)
         bas43 = bas43 + coef3(idx) * u4(ix,jy+idx,kz)
       end do
      else if (jy.ge.2.and.jy.le.ny-2) then
       do idx=-2,2
         bas23 = bas23 + coef2(idx) * u2(ix,jy+idx,kz)
         bas43 = bas43 + coef2(idx) * u4(ix,jy+idx,kz)
       end do
      else
       do idx=-1,1
         bas23 = bas23 + coef1(idx) * u2(ix,jy+idx,kz)
         bas43 = bas43 + coef1(idx) * u4(ix,jy+idx,kz)
       end do
      end if
      bas23 = bas23/dy
      bas43 = bas43/dy

* Z
      bas24 = 0.0
      bas34 = 0.0
      if (kz.ge.4.and.kz.le.nz-4) then
       do idx=-4,4
         bas24 = bas24 + coef4(idx) * u2(ix,jy,kz+idx)
         bas34 = bas34 + coef4(idx) * u3(ix,jy,kz+idx)
       end do
      else if (kz.ge.3.and.kz.le.nz-3) then
       do idx=-3,3
         bas24 = bas24 + coef3(idx) * u2(ix,jy,kz+idx)
         bas34 = bas34 + coef3(idx) * u3(ix,jy,kz+idx)
       end do
      else if (kz.ge.2.and.kz.le.nz-2) then
       do idx=-2,2
         bas24 = bas24 + coef2(idx) * u2(ix,jy,kz+idx)
         bas34 = bas34 + coef2(idx) * u3(ix,jy,kz+idx)
       end do
      else
       do idx=-1,1
         bas24 = bas24 + coef1(idx) * u2(ix,jy,kz+idx)
         bas34 = bas34 + coef1(idx) * u3(ix,jy,kz+idx)
       end do
      end if
      bas24 = bas24/dz
      bas34 = bas34/dz

      ROTAX_0(IX,JY,KZ)=bas43-bas34
      ROTAY_0(IX,JY,KZ)=bas24-bas42
      ROTAZ_0(IX,JY,KZ)=bas32-bas23

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
       real bas22,bas33,bas44,AAA,BBB,CCC
       INTEGER CR1,CR2,CR3
       INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE

       real coef4(-4:4), coef3(-3:3), coef2(-2:2), coef1(-1:1)
       integer idx

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

*-------------------------------------
*      Divergencia fina  (DIVER)
*            (IR: 1 ---> NL)
*            (celdas: 2 ----> NX-1)
*-------------------------------------

       coef4 = (/1.0/280,-4.0/105,1.0/5,-4.0/5,0.0,4.0/5,-1.0/5,4.0/105,
     &           -1.0/280 /)
       coef3 = (/-1.0/60,3.0/20,-3.0/4,0.0,3.0/4,-3.0/20,1.0/60/)
       coef2 = (/1.0/12,-2.0/3,0.0,2.0/3,-1.0/12/)
       coef1 = (/-1.0/2,0.0,1.0/2/)

       DO IR=1,NL

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

* X
       bas22 = 0.0
       if (ix.ge.4.and.ix.le.n1-3) then
        do idx=-4,4
          bas22 = bas22 + coef4(idx) * u12(ix+idx,jy,kz,i)
        end do
       else if (ix.ge.3.and.ix.le.n1-2) then
        do idx=-3,3
          bas22 = bas22 + coef3(idx) * u12(ix+idx,jy,kz,i)
        end do
       else if (ix.ge.2.and.ix.le.n1-1) then
        do idx=-2,2
          bas22 = bas22 + coef2(idx) * u12(ix+idx,jy,kz,i)
        end do
       else
        do idx=-1,1
          bas22 = bas22 + coef1(idx) * u12(ix+idx,jy,kz,i)
        end do
       end if
       bas22 = bas22 / dxpa

* Y
       bas33 = 0.0
       if (jy.ge.4.and.jy.le.n2-3) then
        do idx=-4,4
          bas33 = bas33 + coef4(idx) * u13(ix,jy+idx,kz,i)
        end do
      else if (jy.ge.3.and.jy.le.n2-2) then
        do idx=-3,3
          bas33 = bas33 + coef3(idx) * u13(ix,jy+idx,kz,i)
        end do
      else if (jy.ge.2.and.jy.le.n2-1) then
        do idx=-2,2
          bas33 = bas33 + coef2(idx) * u13(ix,jy+idx,kz,i)
        end do
       else
        do idx=-1,1
          bas33 = bas33 + coef1(idx) * u13(ix,jy+idx,kz,i)
        end do
       end if
       bas33 = bas33 / dypa

* Z
       bas44 = 0.0
       if (kz.ge.4.and.kz.le.n3-3) then
        do idx=-4,4
          bas44 = bas44 + coef4(idx) * u14(ix,jy,kz+idx,i)
        end do
       else if (kz.ge.3.and.kz.le.n3-2) then
        do idx=-3,3
          bas44 = bas44 + coef3(idx) * u14(ix,jy,kz+idx,i)
        end do
       else if (kz.ge.2.and.kz.le.n3-1) then
        do idx=-2,2
          bas44 = bas44 + coef2(idx) * u14(ix,jy,kz+idx,i)
        end do
       else
        do idx=-1,1
          bas44 = bas44 + coef1(idx) * u14(ix,jy,kz+idx,i)
        end do
       end if
       bas44 = bas44/dzpa

       DIVER(IX,JY,KZ,I) = bas22 + bas33 + bas44

       END DO
       END DO
       END DO

       END DO
       END DO


*-------------------------------*
*      COARSE LEVEL
*-------------------------------*
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX

* X
        bas22 = 0.0
        if (ix.ge.4.and.ix.le.nx-4) then
         do idx=-4,4
           bas22 = bas22 + coef4(idx) * u2(ix+idx,jy,kz)
         end do
        else if (ix.ge.3.and.ix.le.nx-3) then
         do idx=-3,3
           bas22 = bas22 + coef3(idx) * u2(ix+idx,jy,kz)
         end do
        else if (ix.ge.2.and.ix.le.nx-2) then
         do idx=-2,2
           bas22 = bas22 + coef2(idx) * u2(ix+idx,jy,kz)
         end do
        else
         do idx=-1,1
           bas22 = bas22 + coef1(idx) * u2(ix+idx,jy,kz)
         end do
        end if
        bas22 = bas22 / dx

* Y
        bas33 = 0.0
        if (jy.ge.4.and.jy.le.ny-4) then
         do idx=-4,4
           bas33 = bas33 + coef4(idx) * u3(ix,jy+idx,kz)
         end do
        else if (jy.ge.3.and.jy.le.ny-3) then
         do idx=-3,3
           bas33 = bas33 + coef3(idx) * u3(ix,jy+idx,kz)
         end do
        else if (jy.ge.2.and.jy.le.ny-2) then
         do idx=-2,2
           bas33 = bas33 + coef2(idx) * u3(ix,jy+idx,kz)
         end do
        else
         do idx=-1,1
           bas33 = bas33 + coef1(idx) * u3(ix,jy+idx,kz)
         end do
        end if
        bas33 = bas33 / dy

* Z
        bas44 = 0.0
        if (kz.ge.4.and.kz.le.nz-4) then
         do idx=-4,4
           bas44 = bas44 + coef4(idx) * u4(ix,jy,kz+idx)
         end do
        else if (kz.ge.3.and.kz.le.nz-3) then
         do idx=-3,3
           bas44 = bas44 + coef3(idx) * u4(ix,jy,kz+idx)
         end do
        else if (kz.ge.2.and.kz.le.nz-2) then
         do idx=-2,2
           bas44 = bas44 + coef2(idx) * u4(ix,jy,kz+idx)
         end do
        else
         do idx=-1,1
           bas44 = bas44 + coef1(idx) * u4(ix,jy,kz+idx)
         end do
        end if
        bas44 = bas44/dz

        DIVER0(IX,JY,KZ) = bas22 + bas33 + bas44

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
       real bas22,bas33,bas44,AAA,BBB,CCC
       INTEGER MARK,ABUELO,KR1,KR2,KR3,IR_ABUE,CR1,CR2,CR3
       real GRAD_P_X,GRAD_P_Y,GRAD_P_Z

       real coef4(-4:4), coef3(-3:3), coef2(-2:2), coef1(-1:1)
       integer idx


*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

*-------------------------------------
*      Divergencia fina  (DIVER)
*            (IR: 1 ---> NL)
*            (celdas: 2 ----> NX-1)
*-------------------------------------

       coef4 = (/1.0/280,-4.0/105,1.0/5,-4.0/5,0.0,4.0/5,-1.0/5,4.0/105,
     &           -1.0/280 /)
       coef3 = (/-1.0/60,3.0/20,-3.0/4,0.0,3.0/4,-3.0/20,1.0/60/)
       coef2 = (/1.0/12,-2.0/3,0.0,2.0/3,-1.0/12/)
       coef1 = (/-1.0/2,0.0,1.0/2/)

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
       DO I=LOW1,LOW2

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

* X
        bas22 = 0.0
        if (ix.ge.2.and.ix.le.n1-1) then
         do idx=-4,4
           bas22 = bas22 + coef4(idx) * diver(ix+idx,jy,kz,i)
         end do
        else if (ix.ge.1.and.ix.le.n1) then
         do idx=-3,3
           bas22 = bas22 + coef3(idx) * diver(ix+idx,jy,kz,i)
         end do
        end if
        bas22 = bas22 / dxpa

* Y
        bas33 = 0.0
        if (jy.ge.2.and.jy.le.n2-1) then
         do idx=-4,4
           bas33 = bas33 + coef4(idx) * diver(ix,jy+idx,kz,i)
         end do
        else if (jy.ge.1.and.jy.le.n2) then
         do idx=-3,3
           bas33 = bas33 + coef3(idx) * diver(ix,jy+idx,kz,i)
         end do
        end if
        bas33 = bas33 / dypa

* Z
        bas44 = 0.0
        if (kz.ge.2.and.kz.le.n3-1) then
         do idx=-4,4
           bas44 = bas44 + coef4(idx) * diver(ix,jy,kz+idx,i)
         end do
        else if (kz.ge.1.and.kz.le.n3) then
         do idx=-3,3
           bas44 = bas44 + coef3(idx) * diver(ix,jy,kz+idx,i)
         end do
        end if
        bas44 = bas44 / dzpa

        U12P(IX,JY,KZ,I)=bas22
        U13P(IX,JY,KZ,I)=bas33
        U14P(IX,JY,KZ,I)=bas44

       END DO
       END DO
       END DO

       END DO
       END DO

*-------------------------------*
*      COARSE LEVEL
*-------------------------------*

       DO KZ=1, NZ
       DO JY=1, NY
       DO IX=1, NX

* X
        bas22 = 0.0
        if (ix.ge.4.and.ix.le.nx-4) then
         do idx=-4,4
           bas22 = bas22 + coef4(idx) * diver0(ix+idx,jy,kz)
         end do
        else if (ix.ge.3.and.ix.le.nx-3) then
         do idx=-3,3
           bas22 = bas22 + coef3(idx) * diver0(ix+idx,jy,kz)
         end do
        else if (ix.ge.2.and.ix.le.nx-2) then
         do idx=-2,2
           bas22 = bas22 + coef2(idx) * diver0(ix+idx,jy,kz)
         end do
        else
         do idx=-1,1
           bas22 = bas22 + coef1(idx) * diver0(ix+idx,jy,kz)
         end do
        end if
        bas22 = bas22 / dx

*  Y
        bas33 = 0.0
        if (jy.ge.4.and.jy.le.ny-4) then
         do idx=-4,4
          bas33 = bas33 + coef4(idx) * diver0(ix,jy+idx,kz)
         end do
        else if (jy.ge.3.and.jy.le.ny-3) then
         do idx=-3,3
           bas33 = bas33 + coef3(idx) * diver0(ix,jy+idx,kz)
         end do
        else if (jy.ge.2.and.jy.le.ny-2) then
         do idx=-2,2
           bas33 = bas33 + coef2(idx) * diver0(ix,jy+idx,kz)
         end do
        else
         do idx=-1,1
           bas33 = bas33 + coef1(idx) * diver0(ix,jy+idx,kz)
         end do
        end if
        bas33 = bas33 / dy

* Z
        bas44 = 0.0
        if (kz.ge.4.and.kz.le.nz-4) then
         do idx=-4,4
           bas44 = bas44 + coef4(idx) * diver0(ix,jy,kz+idx)
         end do
        else if (kz.ge.3.and.kz.le.nz-3) then
         do idx=-3,3
           bas44 = bas44 + coef3(idx) * diver0(ix,jy,kz+idx)
         end do
        else if (kz.ge.2.and.kz.le.nz-2) then
         do idx=-2,2
           bas44 = bas44 + coef2(idx) * diver0(ix,jy,kz+idx)
         end do
        else
         do idx=-1,1
           bas44 = bas44 + coef1(idx) * diver0(ix,jy,kz+idx)
         end do
        end if
        bas44 = bas44 / dz

        U2P(IX,JY,KZ)=bas22
        U3P(IX,JY,KZ)=bas33
        U4P(IX,JY,KZ)=bas44

       END DO
       END DO
       END DO


        RETURN
       END
