************************************************************************
      SUBROUTINE MULTISCALE_FILTER(NX,NY,NZ,NL,NPATCH,pare,
     &            PATCHNX,PATCHNY,PATCHNZ,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,DX,output_iter,flag_w_filtlen,
     &            tol,step,maxit)
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX, NY, NZ, NL
      INTEGER NPATCH(0:NLEVELS), pare(npalev)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      real PATCHrX(NPALEV),PATCHrY(NPALEV),PATCHrZ(NPALEV)
      REAL DX
      INTEGER output_iter, flag_w_filtlen, maxit
      real tol, step

*     global variables
*     original velocity: at the end, the filtered velocity will be
*     stored here
      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

*     filtering scales
      real L0(1:NMAX,1:NMAY,1:NMAZ)
      real L1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

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
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      real dens0(1:NMAX,1:NMAY,1:NMAZ)
      real dens1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /dens/ dens0,dens1

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

      integer solap(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

      integer shock0(1:NMAX,1:NMAY,1:NMAZ)
      integer shock1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

*     Auxiliary variables
      real u2bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u3bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u4bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u12bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real u13bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real u14bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

*     private VARIABLES
      INTEGER IX, JY, KZ, I, J, K, LOW1, LOW2, N1, N2, N3, IR, irr
      integer ii, jj, kk, iixx, jjyy, kkzz, ipatch, jpatch, llow1, llow2
      integer nn1, nn2, nn3
      integer marca, iter
      real BAS1, BAS2, BAS3, bas4, DXPA, L, l2, err, dxpa_i
      real thisx, thisy, thisz, dv2, dv3, dv4, dv2prev, dv3prev, dv4prev
      real lado0
      integer mini, maxi, minj, maxj, mink, maxk
      real rx1, rx2, ry1, ry2, rz1, rz2, rxx1, rxx2, ryy1, ryy2,
     &     rzz1, rzz2
      real u(2,2,2), fuin, uw(2,2,2)
      character*13 filenom
      character*30 filerr5

      integer exectime, time


      lado0 = nx * dx

      dens0 = dens0 + 1.0
      dens1 = dens1 + 1.0

*     DENS0, DENS1 PROPORTIONAL TO CELLS MASSES!!!
      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA = DX / (2.0**IR)
       DO I=LOW1,LOW2
         DENS1(:,:,:,I) = DENS1(:,:,:,I) / 8.0**IR
       END DO
      END DO

      call veinsgrid_all_l(NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &            PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,solap)

      call lee_mach(output_iter,NPATCH,PARE,PATCHNX,PATCHNY,
     &              PATCHNZ,SHOCK0,SHOCK1)

!!!!! for each cell, we ought to find the optimum coherence length
      ! we first initialize the lengths
      L0 = 2.0 * DX
      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA = DX / (2.0**IR)
       DO I=LOW1,LOW2
         L1(:,:,:,I) = 2.0 * DXPA
       END DO
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,L0,RADX,RADY,RADZ,LADO0,CR0AMR,
!$OMP+                   DENS0,U2,U3,U4,NL,NPATCH,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,CR0AMR1,SOLAP,RX,RY,RZ,DENS1,
!$OMP+                   U12,U13,U14,TOL,U2BULK,U3BULK,U4BULK,STEP,DX,
!$OMP+                   SHOCK0,SHOCK1),
!$OMP+            PRIVATE(I,J,K,MARCA,ITER,L,THISX,THISY,THISZ,BAS1,
!$OMP+                    BAS2,BAS3,BAS4,L2,II,JJ,KK,IRR,LLOW1,LLOW2,
!$OMP+                    JPATCH,NN1,NN2,NN3,IIXX,JJYY,KKZZ,ERR,
!$OMP+                    DV2,DV3,DV4,DV2PREV,DV3PREV,DV4PREV,
!$OMP+                    MINI,MAXI,MINJ,MAXJ,MINK,MAXK,rx1,rx2,ry1,
!$OMP+                    ry2,rz1,rz2,rxx1,rxx2,ryy1,ryy2,rzz1,rzz2,
!$OMP+                    dxpa),
!$OMP+            schedule(dynamic)
      do i=1,nx
        do j=1,ny
          do k=1,nz
            marca = 0
            if (cr0amr(i,j,k).eq.1) marca = 1
            iter = 0
            L = L0(i,j,k)
            thisx = radx(i)
            thisy = rady(j)
            thisz = radz(k)

            if (l.gt.min(thisx+0.5*lado0,0.5*lado0-thisx,
     &                     thisy+0.5*lado0,0.5*lado0-thisy,
     &                     thisz+0.5*lado0,0.5*lado0-thisz)) marca=0

            iter_while_c: do while (marca.eq.1)
              bas1 = 0.0
              bas2 = 0.0
              bas3 = 0.0
              bas4 = 0.0
              l2 = l**2

              mini = int(((thisx - l) / lado0 + 0.5) * nx) + 1
              maxi = int(((thisx + l) / lado0 + 0.5) * nx) + 1
              minj = int(((thisy - l) / lado0 + 0.5) * ny) + 1
              maxj = int(((thisy + l) / lado0 + 0.5) * ny) + 1
              mink = int(((thisz - l) / lado0 + 0.5) * nz) + 1
              maxk = int(((thisz + l) / lado0 + 0.5) * nz) + 1
              if (mini.lt.0) mini=1
              if (maxi.gt.nx) maxi=nx
              if (minj.lt.0) minj=1
              if (maxj.gt.ny) maxj=ny
              if (mink.lt.0) mink=1
              if (maxk.gt.nz) maxk=nz

              outer0_c: do ii=mini,maxi
                do jj=minj,maxj
                  do kk=mink,maxk
                    if (cr0amr(ii,jj,kk).eq.1) then
                    if ((radx(ii)-thisx)**2+(rady(jj)-thisy)**2+
     &                  (radz(kk)-thisz)**2.le.l2) then
                      bas1 = bas1 + dens0(ii,jj,kk)
                      bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
                      bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
                      bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
                      if (shock0(ii,jj,kk).eq.1) then
                        if (iter.ge.1) then
                          marca = 0
                          exit outer0_c
                        end if
                      end if
                    end if
                    end if
                  end do
                end do
              end do outer0_c

              if (marca.eq.0) exit iter_while_c

              outer1_c: DO irr=1,NL
               LLOW1=SUM(NPATCH(0:IRR-1))+1
               LLOW2=SUM(NPATCH(0:IRR))
               dxpa = dx / 2.0 ** irr
               DO jpatch=LLOW1,LLOW2
                 nn1 = patchnx(jpatch)
                 nn2 = patchny(jpatch)
                 nn3 = patchnz(jpatch)

                 RX1=PATCHRX(jpatch)-0.5*dxpa
                 RY1=PATCHRY(jpatch)-0.5*dxpa
                 RZ1=PATCHRZ(jpatch)-0.5*dxpa
                 RX2=PATCHRX(jpatch)-0.5*dxpa+(nn1-1)*dxpa
                 RY2=PATCHRY(jpatch)-0.5*dxpa+(nn2-1)*dxpa
                 RZ2=PATCHRZ(jpatch)-0.5*dxpa+(nn3-1)*dxpa

                 RXX1 = thisx - l
                 RXX2 = thisx + l
                 RYY1 = thisy - l
                 RYY2 = thisy + l
                 RZZ1 = thisz - l
                 RZZ2 = thisz + l

                 IF (rxx1.le.rx2.AND.rx1.le.rxx2.AND.
     &               ryy1.le.ry2.AND.ry1.le.ryy2.AND.
     &               rzz1.le.rz2.AND.rz1.le.rzz2) then

                  !X
                  IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
                     mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                     maxi=INT(((RXX2-RX1)/DXPA)) + 1
                  END IF
                  IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
                     mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                     maxi=nn1
                  END IF
                  IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
                     mini=1
                     maxi=INT(((RXX2-RX1)/DXPA)) + 1
                  END IF
                  IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
                     mini=1
                     maxi=nn1
                  END IF

                  !Y
                  IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
                     minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                     maxj=INT(((RYY2-RY1)/dxpa)) + 1
                  END IF
                  IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
                     minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                     maxj=nn2
                  END IF
                  IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
                     minj=1
                     maxj=INT(((RYY2-RY1)/dxpa)) + 1
                  END IF
                  IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
                     minj=1
                     maxj=nn2
                  END IF

                  !Z
                  IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
                     mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                     maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                  END IF
                  IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
                     mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                     maxk=nn3
                  END IF
                  IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
                     mink=1
                     maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                  END IF
                  IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
                     mink=1
                     maxk=nn3
                  END IF

                   do iixx = mini,maxi
                   do jjyy = minj,maxj
                   do kkzz = mink,maxk
                     if (cr0amr1(iixx,jjyy,kkzz,jpatch).eq.1.and.
     &                 solap(iixx,jjyy,kkzz,jpatch).eq.1) then
                     if ((rx(iixx,jpatch)-thisx)**2+
     &                 (ry(jjyy,jpatch)-thisy)**2+
     &                 (rz(kkzz,jpatch)-thisz)**2.le.l2) then
                       bas1 = bas1 + dens1(iixx,jjyy,kkzz,jpatch)
                       bas2 = bas2 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u12(iixx,jjyy,kkzz,jpatch)
                       bas3 = bas3 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u13(iixx,jjyy,kkzz,jpatch)
                       bas4 = bas4 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u14(iixx,jjyy,kkzz,jpatch)
                       if (shock1(iixx,jjyy,kkzz,jpatch).eq.1) then
                         if (iter.ge.1) then
                           marca = 0
                           exit outer1_c
                         end if
                       end if
                     end if
                     end if
                   end do
                   end do
                   end do
                END IF
               END DO
             end do outer1_c

              if (marca.eq.0) exit iter_while_c

              bas2 = bas2 / bas1
              bas3 = bas3 / bas1
              bas4 = bas4 / bas1

              !!! 2. ERROR CALCULATION
              if (iter.eq.1) then
                err = 2.0 * tol
              else
                dv2prev = u2(i,j,k) - u2bulk(i,j,k)
                dv3prev = u3(i,j,k) - u3bulk(i,j,k)
                dv4prev = u4(i,j,k) - u4bulk(i,j,k)
                dv2 = u2(i,j,k) - bas2
                dv3 = u3(i,j,k) - bas3
                dv4 = u4(i,j,k) - bas4
                err = max(abs((dv2/dv2prev-1.0)), abs(dv3/dv3prev-1.0),
     &                    abs(dv4/dv4prev-1.0))
              end if

              u2bulk(i,j,k) = bas2
              u3bulk(i,j,k) = bas3
              u4bulk(i,j,k) = bas4

*             stop condition: min error or max num of its
              if (err.le.tol.or.iter.gt.maxit) marca = 0
*             stop when the growing sphere touches the domain boundary
*             (so we do not worry about boundary condition)
              if (l.gt.min(thisx+0.5*lado0,0.5*lado0-thisx,
     &                     thisy+0.5*lado0,0.5*lado0-thisy,
     &                     thisz+0.5*lado0,0.5*lado0-thisz)) marca=0

*             the sphere grows
              if (marca.eq.1) l = max(l*step, l+dx)
              iter = iter+1
            end do iter_while_c
            l0(i,j,k) = l
            !if (cr0amr(i,j,k).eq.1) write(*,*) i,j,k,iter,l,err !DEBUGGING
          end do
        end do
        !write(*,*) i, 'done'
      end do ! do i=1,nx


 !!! REFINED CELLS
      LOW1=1
      LOW2=SUM(NPATCH)
!$OMP PARALLEL DO SHARED(L1,radx,rady,radz,LADO0,CR0AMR,
!$OMP+                   DENS0,U2,U3,U4,NL,NPATCH,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,CR0AMR1,SOLAP,RX,RY,RZ,DENS1,
!$OMP+                   U12,U13,U14,TOL,U12BULK,U13BULK,U14BULK,STEP,
!$OMP+                   DX,SHOCK0,SHOCK1,LOW1,LOW2),
!$OMP+            PRIVATE(I,J,K,MARCA,ITER,L,THISX,THISY,THISZ,BAS1,
!$OMP+                    BAS2,BAS3,BAS4,L2,II,JJ,KK,IRR,LLOW1,LLOW2,
!$OMP+                    JPATCH,NN1,NN2,NN3,IIXX,JJYY,KKZZ,ERR,
!$OMP+                    DV2,DV3,DV4,DV2PREV,DV3PREV,DV4PREV,
!$OMP+                    MINI,MAXI,MINJ,MAXJ,MINK,MAXK,rx1,rx2,ry1,
!$OMP+                    ry2,rz1,rz2,rxx1,rxx2,ryy1,ryy2,rzz1,rzz2,
!$OMP+                    ipatch,n1,n2,n3,exectime,dxpa,dxpa_i,ir),
!$OMP+            schedule(dynamic)
      DO ipatch=LOW1,LOW2
        exectime = time()
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        ! get the level
        DO IR=1,NL
          IF (IPATCH.LE.SUM(NPATCH(0:IR))) EXIT
        END DO
        dxpa_i = dx/(2.0**ir)

c        write(*,*) ipatch, sum(cr0amr1(1:n1,1:n2,1:n3,ipatch) *
c     &                         solap(1:n1,1:n2,1:n3,ipatch))
        do i=1,n1
        do j=1,n2
        do k=1,n3
          marca = 0
          if (cr0amr1(i,j,k,ipatch).eq.1.and.
     &          solap(i,j,k,ipatch).eq.1) marca = 1
          iter = 0
          L = L1(i,j,k,ipatch)
          thisx = rx(i,ipatch)
          thisy = ry(j,ipatch)
          thisz = rz(k,ipatch)

          if (l.gt.min(thisx+0.5*lado0,0.5*lado0-thisx,
     &                     thisy+0.5*lado0,0.5*lado0-thisy,
     &                     thisz+0.5*lado0,0.5*lado0-thisz)) marca=0

          iter_while: do while (marca.eq.1)
            bas1 = 0.0
            bas2 = 0.0
            bas3 = 0.0
            bas4 = 0.0
            l2 = l**2

            mini = int(((thisx - l) / lado0 + 0.5) * nx) + 1
            maxi = int(((thisx + l) / lado0 + 0.5) * nx) + 1
            minj = int(((thisy - l) / lado0 + 0.5) * ny) + 1
            maxj = int(((thisy + l) / lado0 + 0.5) * ny) + 1
            mink = int(((thisz - l) / lado0 + 0.5) * nz) + 1
            maxk = int(((thisz + l) / lado0 + 0.5) * nz) + 1
            if (mini.lt.0) mini=1
            if (maxi.gt.nx) maxi=nx
            if (minj.lt.0) minj=1
            if (maxj.gt.ny) maxj=ny
            if (mink.lt.0) mink=1
            if (maxk.gt.nz) maxk=nz

            outer0: do ii=mini,maxi
              do jj=minj,maxj
                do kk=mink,maxk
                  if (cr0amr(ii,jj,kk).eq.1) then
                  if ((radx(ii)-thisx)**2+(rady(jj)-thisy)**2+
     &                  (radz(kk)-thisz)**2.le.l2) then
                    bas1 = bas1 + dens0(ii,jj,kk)
                    bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
                    bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
                    bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
                    if (shock0(ii,jj,kk).eq.1) then
                      if (iter.ge.1) then
                        marca=0
                        exit outer0
                      end if
                    end if
                  end if
                  end if
                end do
              end do
            end do outer0

            if (marca.eq.0) exit iter_while

            outer1: DO irr=1,NL
             LLOW1=SUM(NPATCH(0:IRR-1))+1
             LLOW2=SUM(NPATCH(0:IRR))
             dxpa = dx / (2.0 ** irr)
             DO jpatch=LLOW1,LLOW2
               nn1 = patchnx(jpatch)
               nn2 = patchny(jpatch)
               nn3 = patchnz(jpatch)

               RX1=PATCHRX(jpatch)-0.5*dxpa
               RY1=PATCHRY(jpatch)-0.5*dxpa
               RZ1=PATCHRZ(jpatch)-0.5*dxpa
               RX2=PATCHRX(jpatch)-0.5*dxpa+(nn1-1)*dxpa
               RY2=PATCHRY(jpatch)-0.5*dxpa+(nn2-1)*dxpa
               RZ2=PATCHRZ(jpatch)-0.5*dxpa+(nn3-1)*dxpa

               RXX1 = thisx - l
               RXX2 = thisx + l
               RYY1 = thisy - l
               RYY2 = thisy + l
               RZZ1 = thisz - l
               RZZ2 = thisz + l

               IF (rxx1.le.rx2.AND.rx1.le.rxx2.AND.
     &               ryy1.le.ry2.AND.ry1.le.ryy2.AND.
     &               rzz1.le.rz2.AND.rz1.le.rzz2) then

                !X
                IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
                   mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                   maxi=INT(((RXX2-RX1)/DXPA)) + 1
                END IF
                IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
                   mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                   maxi=nn1
                END IF
                IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
                   mini=1
                   maxi=INT(((RXX2-RX1)/DXPA)) + 1
                END IF
                IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
                   mini=1
                   maxi=nn1
                END IF

                !Y
                IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
                   minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                   maxj=INT(((RYY2-RY1)/dxpa)) + 1
                END IF
                IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
                   minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                   maxj=nn2
                END IF
                IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
                   minj=1
                   maxj=INT(((RYY2-RY1)/dxpa)) + 1
                END IF
                IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
                   minj=1
                   maxj=nn2
                END IF

                !Z
                IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
                   mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                   maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                END IF
                IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
                   mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                   maxk=nn3
                END IF
                IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
                   mink=1
                   maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                END IF
                IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
                   mink=1
                   maxk=nn3
                END IF

                 do iixx = mini,maxi
                 do jjyy = minj,maxj
                 do kkzz = mink,maxk
                   if (cr0amr1(iixx,jjyy,kkzz,jpatch).eq.1.and.
     &                 solap(iixx,jjyy,kkzz,jpatch).eq.1) then
                   if ((rx(iixx,jpatch)-thisx)**2+
     &                 (ry(jjyy,jpatch)-thisy)**2+
     &                 (rz(kkzz,jpatch)-thisz)**2.le.l2) then
                     bas1 = bas1 + dens1(iixx,jjyy,kkzz,jpatch)
                     bas2 = bas2 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u12(iixx,jjyy,kkzz,jpatch)
                     bas3 = bas3 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u13(iixx,jjyy,kkzz,jpatch)
                     bas4 = bas4 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u14(iixx,jjyy,kkzz,jpatch)
                     if (shock1(iixx,jjyy,kkzz,jpatch).eq.1) then
                       if (iter.ge.1) then
                         marca=0
                         exit outer1
                       end if
                     end if
                   end if
                   end if
                 end do
                 end do
                 end do
              END IF
             END DO
           end do outer1

            if (marca.eq.0) exit iter_while

            bas2 = bas2 / bas1
            bas3 = bas3 / bas1
            bas4 = bas4 / bas1

            !!! 2. ERROR CALCULATION
            if (iter.eq.1) then
              err = 2.0 * tol
            else
              dv2prev = u12(i,j,k,ipatch) - u12bulk(i,j,k,ipatch)
              dv3prev = u13(i,j,k,ipatch) - u13bulk(i,j,k,ipatch)
              dv4prev = u14(i,j,k,ipatch) - u14bulk(i,j,k,ipatch)
              dv2 = u12(i,j,k,ipatch) - bas2
              dv3 = u13(i,j,k,ipatch) - bas3
              dv4 = u14(i,j,k,ipatch) - bas4
              err = max(abs((dv2/dv2prev-1.0)), abs(dv3/dv3prev-1.0),
     &                    abs(dv4/dv4prev-1.0))
            end if

            u12bulk(i,j,k,ipatch) = bas2
            u13bulk(i,j,k,ipatch) = bas3
            u14bulk(i,j,k,ipatch) = bas4

*             stop condition: min error or max num of its
            if (err.le.tol.or.iter.gt.maxit) marca = 0
*             stop when the growing sphere touches the domain boundary
*             (so we do not worry about boundary condition)
            if (l.gt.min(thisx+0.5*lado0,0.5*lado0-thisx,
     &                     thisy+0.5*lado0,0.5*lado0-thisy,
     &                     thisz+0.5*lado0,0.5*lado0-thisz)) marca=0

*             the sphere grows
            if (marca.eq.1) l = max(l*step, l+dxpa_i)
            iter = iter+1
          end do iter_while
          l1(i,j,k,ipatch) = l

          !DEBUGGING
C            if (cr0amr1(i,j,k,ipatch).eq.1.and.
C     &          solap(i,j,k,ipatch).eq.1) write(*,*) 'amr:',ipatch,i,j,
C     &                                                k,iter,l,err
          !END DEBUGGING
        end do
        end do
        end do ! do i=1,n1
        exectime = time() - exectime
        !write(*,*) ipatch, sum(cr0amr1(1:n1,1:n2,1:n3,ipatch) *
     &             solap(1:n1,1:n2,1:n3,ipatch)), exectime, dx, dxpa_i
      END DO


*     refill refined and overlapping cells
      DO IR=NL,1,-1
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,L1)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U12BULK)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U13BULK)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U14BULK)

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
            II = PATCHX(JPATCH) + int((I-1)/2)
            JJ = PATCHY(JPATCH) + int((J-1)/2)
            KK = PATCHZ(JPATCH) + int((K-1)/2)
            if (jpatch.ne.0) then
              uw(1:2,1:2,1:2) = dens1(I:I+1,J:J+1,K:K+1,IPATCH)

              u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              L1(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u12bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u12bulk(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u13bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u13bulk(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u14bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u14bulk(II,JJ,KK,JPATCH) = FUIN
            else
              uw(1:2,1:2,1:2) = dens1(I:I+1,J:J+1,K:K+1,IPATCH)

              u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              L0(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u12bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              U2bulk(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u13bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u3bulk(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u14bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u4bulk(II,JJ,KK) = FUIN
            end if
          END DO
          END DO
          END DO
        END DO
      END DO !IR=NL,1,-1

      ! U2,U3,U4,U12,U13,U14 gets updated with the values of the
      ! velocity fluctuation
      DO ipatch=1,sum(npatch)
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        do i=1,n1
        do j=1,n2
        do k=1,n3
          u12(i,j,k,ipatch)=u12(i,j,k,ipatch)-u12bulk(i,j,k,ipatch)
          u13(i,j,k,ipatch)=u13(i,j,k,ipatch)-u13bulk(i,j,k,ipatch)
          u14(i,j,k,ipatch)=u14(i,j,k,ipatch)-u14bulk(i,j,k,ipatch)
        end do
        end do
        end do
      end do

      do i=1,nx
        do j=1,ny
          do k=1,nz
            u2(i,j,k) = u2(i,j,k) - u2bulk(i,j,k)
            u3(i,j,k) = u3(i,j,k) - u3bulk(i,j,k)
            u4(i,j,k) = u4(i,j,k) - u4bulk(i,j,k)
          end do
        end do
      end do

      if (flag_w_filtlen.eq.1) then
        CALL NOMFILE_FILTLEN(output_iter,FILENOM)
        FILERR5 = './output_files/'//FILENOM
        CALL WRITE_FILTLEN(FILERR5,NX,NY,NZ,ITER,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ,L0,L1)
      end if

      RETURN
      END


************************************************************************
       SUBROUTINE finer_to_coarser(u,uw,fuin)
************************************************************************
        REAL U(2,2,2)
        REAL UW(2,2,2) ! weights!!
        REAL FUIN, WSUM
        INTEGER I,J,K

        FUIN = 0.0
        WSUM = 0.0
        DO I=1,2
          DO J=1,2
            DO K=1,2
              FUIN = FUIN + U(I,J,K)*UW(I,J,K)
              WSUM = WSUM + UW(I,J,K)
            END DO
          END DO
        END DO

        FUIN = FUIN/WSUM

        RETURN
       END


************************************************************************
       SUBROUTINE veinsgrid_all_l(NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &            solap)
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

       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
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

*!$OMP PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
*!$OMP+    PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI),
*!$OMP+  PRIVATE(I,L1,L2,L3,N1,N2,N3,CR1,CR2,CR3,NV,J,LL1,
*!$OMP+         LL2,LL3,NN1,NN2,NN3,CR4,CR5,CR6,A1,A2,B1,B2,
*!$OMP+         C1,C2)
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         SOLAP(:,:,:,I)=1

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
              IF (SOLAP(IX,JY,KZ,I).EQ.1) SOLAP(II,JJ,KK,J)=0
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


************************************************************************
       SUBROUTINE SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &            VARIABLE)
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
       INTEGER N1,N2,N3,L1,L2,L3
       INTEGER NN1,NN2,NN3,LL1,LL2,LL3
       INTEGER KZ2,JY2,IX2,I2
       INTEGER, ALLOCATABLE::VECINO(:,:)
       INTEGER, ALLOCATABLE::NVECI(:)
       INTEGER NV,A2,B2,C2,K,LOW1, LOW2

       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER MARCA(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL VARIABLE(NAMRX,NAMRY,NAMRZ,NPALEV) !PARAMETER

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

         SOLAP(:,:,:,I)=1

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
              IF (SOLAP(IX,JY,KZ,I).EQ.1) THEN
                SOLAP(II,JJ,KK,J)=0
                VARIABLE(II,JJ,KK,J) = VARIABLE(IX,JY,KZ,I)
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
