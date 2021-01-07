      subroutine CROSS_4HE(xb,d2g,dq2,ebeam,eh,phideg,dlambda,crosstot)
      IMPLICIT REAL*8(A-H,O-Z)
    
cccccc xb is the experimental Bjorken variable defined as xb = q^2/(2 M_NUCLEON nu)
cccccc dq2 is the virtuality of the initial photon. It in expressed in GeV^2
ccccccd2g is the momentum trasferred to the nuclear system. It is negative and expressed in GeV^2.
      

      phidegg = (180.d0-phideg)
      CALL IMAGINARYCFF(xb,d2g,dq2,ris)
      CALL REALCFF(xb,dq2,d2g,ris2)
      
      
      CALL COHDVCSXSECTION(xb,dq2,d2g,ebeam,eh,dlambda,
     1     phidegg,ris, ris2,crosstot)



 1000 format(4(d14.5,1x))
      return
      END


CCCCCC
CCCCCC THE FOLLOWING SUBROUTINES REPRODUCES THE RESULTS PUBLISHED IN PRC 98, 015203 (2018)   [1].
CCCCCC


cccccc
cccccc Imaginary part of teh compton Form Factor
cccccc

      SUBROUTINE IMAGINARYCFF(XB,D2G,DQ2,DIMGPD)
CCCCC D2G IS NEGATIVE
      IMPLICIT REAl*8(A-H,O-Z)

      DIMENSION gpdexp(2),gpdb(2)
      dimension nz(3),az(4),np(3),ap(4),aph(4),nph(3)
      DIMENSION WW(100),DPG(100),z(100),wwz(100)
      DIMENSION wff(100),xx(100),ff(2),xp(100),wx(100)
      parameter(pi=dacos(-1.d0),dc=1.d0/197.31d0)
      common/nk/dnk(100),dnkg(100),dk(100)
      common/ind/nn
      common/ph/phi(100),wphi(100)
      common/cose/dma,xi,dks,dma1,dma1ex,dm,ea,xiacheck,dxiamax
      common/delta/de0,d3,dy,dx,dvec2,dp,d2


CCCCCCC
CCCCCCC     DMA=MASS OF THE 4HE NUCLEUS EXPRESSED IN fm^-1
CCCCCCC     DMP=PROTON MASS EXPRESSED IN fm^-1
CCCCCCC     DMN=NEUTRON MASS EXPRESSED IN fm^-1
CCCCCCC     IA=MASS NUMBER
CCCCCCC

      IA = 4
      DMP = 938.3d0*DC
      DMN = 939.6D0*DC
      EA = - 28.295D0*DC
      DMA = 2.D0*(DMP+DMN) + EA

      open(unit=9,file='data/nk.dat')
      open(unit=19,file='data/nk0.dat')
      do ik = 1, 41   !READ THE MOMENTUM DISTRIBUTION FROM EXTERNAL FILES
         read(9,*) dk(ik), dnk(ik)
         read(19,*) dk(ik), dnkg(ik)
      enddo
      close(unit=9)
      close(unit=19)

      ISR=1      !FIX THE LAB REFERENCE FRAME (ISR=2 IS THE BREIT FRAME. IMPORTANT TO CALCULATE THE COMPONENTS OF DELTA^2)

c      CALL GAULEG(0.d0,2.d0*pi,PHI,WPHI,40)

ccccccc
       mph = 1
       aph(1) = 0.d0
       aph(2) = 2.d0*pi
       nph(1) = 20
       CALL DMPP(Mph,NPOINTphi,Nph,Aph,phi,wphi)
cccccc      
         D2 = D2G * 1.d6 * dc*dc
         dks = 1.d0+0.68d0/(1.d0+0.52d0*dlog(dq2/4.d0))


         do nn = 1, 2  ! LOOP OVER THE NUCLEONS (NN)

ccccccc
C     EMED IS THE AVERAGE ENERGY FOR THE EXCITED A-1 BODY STATE
C     DMA1 IS THE MASS OF THE A-1 BODY SYSTEM
C     E3 IS THE BINDING ENERGY OF THE A-1 BODY SYSTEM (3H OR 3HE)
C     DMAIEX IS THE MASS OF EXCITED THE A-1 BODY SYSTEM
ccccccc

            emed = (105.69d0-26.d0)*DC

            IF(NN.Eq.1) then
               DM = DMP
               E3 = -8.4820D0*dc
               DMA1 = 2.d0*DMN + DMP +E3

            else
               DM = DMN
               E3 = -7.718D0*dc
               DMA1 = 2.d0*DMP + DMN + E3

            endif

            DMA1EX = DMA1 + EMED
ccccccc
C XB=Q^2/(2*M*\NU) RANGES IN [0,A]; XI RANGES IN [0;1]
C XB IS THE EXPERIMENTAL X BJORKEN FOR THE BOUND PROTON.
C dm/dma*XB IS THE  BJORKEN VARIABLE FOR THE WHOLE NUCLEUS
ccccccc
            xi = dmp/dma*xb / ( 2.d0 -  dmp/dma*xb )

            XIAcheck = XI* (DMA/DMP)

            DXIAMAX = (DMA)/(DMP)*dSQRT(-D2)/(DSQRT(4.D0*DMA**2 - D2))
            IF(XIAcheck.GT.DXIAMAX) then
               WRITE (6,*) 'prob immmm'
            endif
               if(ISR.eq.1)then

                  DE0 = -D2/(2.d0*DMA)
                  DP = -DSQRT(2.D0)*DMA*XI/(XI+1.D0)
                  D3 = DSQRT(2.D0)*DP-DE0
                  DME = (DE0-D3)/DSQRT(2.D0)
                  DPERP = DSQRT(DE0**2-D2-D3**2)
                  DY = 0.d0
                  DX= DSQRT(DPERP**2-DY**2)

               else

                  D3= DMA*XI*2.d0/(XI+1.D0)
                  DE0 = 0.d0
                  dp= D3/DSQRT(2.D0)
                  DME=-D3/ DSQRT(2.D0)
                  DPERP=DSQRT(-D2-d3**2)
                  DY=0.D0
                  DX=dsqrt(DPERP**2-DY*DY)

               endif

               DVEC2 = DE0**2-D2

      if(dvec2.lt.0.d0) return

       do isign = 1, 2                    !LOOP OVER THE SIGN OF XI TO EVALUATE THE IMAGINARY PART OF CFF (ISIGN)


               if(isign.eq.1) then

                  x = xi

               else

                  x = - xi

               endif

                TILDM= DM*(DMA+DP/DSQRT(2.D0))/(DMA)

                    ZMIN = dabs(x)

                    zmax = 0.625d0

                    dintz=0.d0

c          call gauleg(zmin,zmax,z,wwz,100)

cccccccc
cccccccc                    
                mz = 2
                az(1) = zmin
                dintz=0.d0

               if(zmin.le.0.24d0) then
                  az(2) = 0.25d0
               else
                  az(2) = (zmin+zmax)/2.d0
               endif
               az(3) = zmax
               nz(1) = 24
               nz(2) = 24
               
               CALL DMPP(MZ,NPOINTZ,NZ,AZ,Z,wWZ)
cccccccc
cccccccc                    

                  do iz = 1, npointz   !LOOP OVER THE LIGHTCONE MOMENTUM FRACTION (Z)

                     TILDZ = dma/dm*(z(iz) + xi)

                     xg = x/z(iz)
                     xig = xi/z(iz)

                     CALL GPDNUCLEON(XG,XIG,D2G,1,dq2,nn,GPDU)
                     CALL GPDNUCLEON(XG,XIG,D2G,2,dq2,nn,GPDD)
                     CALL GPDNUCLEON(XG,XIG,D2G,3,dq2,nn,GPDS)
                     gpdsu= dks*gpds

                    if(nn.eq.1)then
                       GPDT=2.d0*((4.D0/9.d0*(gpdu+gpdsu)+1.d0/9.d0*
     *                       (gpdd+gpdsu)+1.d0/9.d0*gpds))
                     else
                       GPDT=2.d0*(( 1.D0/9.d0*(gpdu+gpdsu)+4.d0/9.d0*
     *                       (gpdd+ gpdsu)+1.d0/9.d0*gpds))
                    endif

                     DLCMDG = 0.d0
                     DLCMDEX = 0.d0

                     do iphi = 1, 40      !LOOP OVER THE AZIMUTHAL ANGLE OF THE INNER NUCLEON (PHI)

                       dcosphi = dcos(PHI(iphi))
                       dsinphi = dsqrt(1.d0-dcosphi*dcosphi)

                        do ig = 1, 2    ! LOOP OVER THE GROUND & THE ECITED PART OF THE SPECTRAL FUNCTION (IG)

                           if(ig.eq.1) then

                     DNUM = DMA1**2 - DMA*DMA*(1.d0-TILDM*TILDZ/DMA)**2
                     cont1=0.d0

                           else

                    DNUM = DMA1EX**2 - DMA*DMA*(1.d0-TILDM*TILDZ/DMA)**2
                    cont2=0.d0

                           endif

                    PMIN = 0.5d0*DABS(DNUM/(DMA*(1.d0-TILDM*TILDZ/DMA)))


                           if(pmin.gt.8.d0) then

                              cont1 = 0.d0
                              cont2 = 0.d0
                              go to 500


                           endif

ccccccccc
                        mp = 1
                        ap(1) = pmin
                        ap(2) = 8.d0
                        np(1) = 24
                        CALL DMPP(Mp,NPOINTp,Np,Ap,dpg,wW)

                           
ccccccccc
c                           CALL GAULEG(PMIN,8.d0,DPG,WW,90)
ccccccc
C
C For details about our model for then non-diagonal spectral function see Eqs.(13)-(14) in [1].
C
cccccccc
C     DPG IS THE INTEGRATION VARIABLE P
C     DPDG IS (P+DELTA)
C     DNPG E DNPEX ARE THE MOMENTUM DISTRIBUTION (G=GROUND, EX=EXCITED)
C     PROTON AND NEUTRON CONTRIBUTE IN THE SAME WAY TO THE MOM. DISTRIBUTION
C     THETA IS THE ANGLE BRETWEEN P AND THE Z AXIS
C     TLDZ IS THE VARIABLE DEFINED AS TLDZ=Z+(MA/M)XI
C     CONT1 AND CONT2 ARE THE GROUND AND THE EXCITED LIGHT-CONE MOMENTUM DISTRIBUTION, RESPECTIVELY.
cccccccc
                           DO ip = 1 ,npointp     ! LOOP OVER THE MOMENTUM OF THE INNER NUCLEON (P)

                              if(ig.eq.1) then

                               dnpg=dnkgr(DPG(IP))
                               DP0  = DMA - DSQRT(DMA1**2 + DPG(IP)**2)
                               DCOST=(TILDM*tildz-dp0)/DPg(ip)
                               DSINT=dsqrt(1.d0-dcost*dcost)
                               if(dabs(dcost).gt.1.d0) then
                               write(6,*) dcost,'WARNING:cosine
     1                             out of range'
                                  stop
                               endif
                               DPDG=DSQRT(DPG(IP)**2 + DVEC2 + 2.d0*
     1                         (  DPG(IP)*DSINT*DCOSPHI*DX +
     1                              DPG(IP)*DSINT*DSINPHI*DY + DPG(IP)*
     *                              DCOST*D3  ))
                                if(dpdg.gt.8.d0)then
                                  dpdg = 7.99d0
                                  endif

                               DNPDG = DNKGR(DPDG)
                               DNTOTG=DSQRT(DNPG*DNPDG)

                         CONT1 = CONT1 + TILDM*DPG(IP)*DNTOTG * WW(IP)

                              else

                               DNPEX = dnktot(DPG(IP))-dnkgr(DPG(IP))
                               DP0 = DMA - DSQRT(DMA1EX**2 + DPG(IP)**2)
                               DCOST=(TILDM*tildz-dp0)/DPg(ip)
                               DSINT=dsqrt(1.d0-dcost*dcost)
                               if(dabs(dcost).gt.1.d0) then
                                  write(6,*) 'WARNING:cosine
     1                             out of range',dcost
                                  stop
                               endif
                               DPDG=DSQRT(DPG(IP)**2 + DVEC2 +2.d0*
     1                              (  DPG(IP)*DSINT*DCOSPHI*DX +
     1                              DPG(IP)*DSINT*DSINPHI*DY + DPG(IP)*
     *                              DCOST*D3  ))
                               if(dpdg.gt.8.d0)then
                                  dpdg = 7.99d0
                                  endif
                               DNPDEX = dnktot(DPDG)-dnkgr(DPDG)
                               DNTOTEX=DSQRT(DNPEX * DNPDEX)

                             CONT2 = CONT2 + TILDM*DPDG*DNTOTEX *WW(IP)

                          endif

                           ENDDO    ! END LOOP OVER P


                        ENDDO    !END LOOP OVER IG

ccccccc light cone momentum distribution G= ground EX= excited

 500                    DLCMDG = DLCMDG + CONT1 * WPHI(IPHI)
                        DLCMDEX = DLCMDEX + CONT2 * WPHI(IPHI)

                     enddo    !END LOOP OVER PHI

                     DLCTOT = (DLCMDG + DLCMDEX)*GPDT/z(iz)

cccccc
C    dintz is the GPD H for the 4He
C    SEE EQ.(10) IN [1]
cccccc

                     dintz= dintz+ dlctot*wwz(iz)

                  enddo      !END LOOP OVER Z

               gpdb(isign)=dintz

            enddo       !END LOOP OVER ISIGN

            gpdexp(nn) =(gpdb(1) - gpdb(2))*dma/dm

         enddo       !END LOOP OVER IN

          DIMGPD= (gpdexp(1) + gpdexp(2))


      RETURN
      END


c*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

cccccc
cccccc SUBROUTINE THAT CALCULATES THE REAL PART OF THE COMPTON FORM FACTOR OF THE 4HE NUCLEUS. DCAUCH CALLS A CERN SUBROUTINE TO CALCULATE THE PART  WITH THE SINGULARITY (1/(X -XI)) IN THE INTEGRAL OVER X.
cccccc
      subroutine REALCFF(xb,dq2,d2g,risretot)
      IMPLICIT REAl*8(A-H,O-Z)

      dimension risu(2),risuu(2)
      dimension xp(150),ww(150)
      external f,g
      dimension  ax(3), nx(2) 
      parameter(pi=dacos(-1.d0),dc=1.d0/197.31d0)
      common/nk/dnk(100),dnkg(100),dk(100)
      common/ind/nn
      common/cose/dma,xi,dks,dma1,dma1ex,dm,ea,xiacheck,dxiamax
      common/kin/xbb,dq22,d2gg

      xbb = xb
      dq22 = dq2
      d2gg = d2g

      IF(XIAcheck.GT.DXIAMAX) then
      write(6,*)'prob re'
      endif
      call system_clock ( j_s, j_r )
          do nn = 1, 2
ccccc
ccccc PARAMETERS REQUIRED BY THE CERN SUBROUTINE DCAUCH
ccccc
      call system_clock ( j_f, j_r )
c      write(*,*) 'Time a: ', nn, ' ', (j_f-j_s)/real(j_r), ' seconds'

            aa = 0.d0
            bb = 0.625d0
            dl = xi
            eps = 0.01d0
            ris = dcauch(f,aa,bb,dl,eps)
            risu(nn) = ris*dma/dm


      call system_clock ( j_f, j_r )
c      write(*,*) 'Time b: ', nn, ' ', (j_f-j_s)/real(j_r), ' seconds'

c     CALL GAULEG(0.d0,0.625d0,xp,ww,150)

            mx = 1
            ax(1) = aa
            ax(3) = bb
                if(aa.le.0.24d0) then
                  ax(2) = 0.25d0
               else
                  ax(2) = (aa+bb)/2.d0
               endif
            nx(1)=44
            nx(2)=44
            
            CALL DMPP(mx,NPOINTxpr,Nx,Ax,xp,ww)

         dintx=0.d0
ccccccc
ccccccc CALCULATE THE PART OF THE INTEGRAL OVER X WITH KERNEL 1/(X+XI)
ccccccc
         do ipr = 1,npointxpr

               xpr = xp(ipr)
               dker = ( 1.d0/(xpr+xi))
               dintx = dintx +g(xpr)*dker*ww(ipr)

         enddo
         risuu(nn)= dintx*dma/dm

      enddo                             !END LOOP OVER NUCLEONS

         rispr=( risu(1) + risu(2))
         risre =( risuu(1) + risuu(2))
         risretot = rispr + risre

      RETURN
      END

C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
cccccc
cccccc INTEGRAND REQUEST BY THE PREVIOUS SUBROUTINE
cccccc
      double precision function g(xpr)

      implicit real*8(a-h,o-z)

      DIMENSION WW(500),DPG(100),z(100),wwz(100)
      dimension nz(3),az(4),np(3),ap(4),aph(2),nph(3)
      parameter(pi=dacos(-1.d0),dc=1.d0/197.31d0)
      common/nk/dnk(100),dnkg(100),dk(100)
      common/ph/phi(100),wphi(100)
      common/ind/nn
      common/cose/dma,xi,dks,dma1,dma1ex,dm,ea,xiacheck,dxiamax
      common/kin/xbb,dq22,d2gg
      common/delta/de0,d3,dy,dx,dvec2,dp,d2

       

       mph = 1
       aph(1) = 0.d0
       aph(2) = 2.d0*pi
       nph(1) = 20
       CALL DMPP(Mph,NPOINTphi,Nph,Aph,phi,wphi)
      
      
       xb = xbb
       dq2 = dq22
       d2g =  d2gg


            TILDM= DM*(DMA+DP/DSQRT(2.D0))/(DMA)

cc            CALL GAULEG(0.d0,2.d0*pi,PHI,WPHI,20)

               dintz=0.d0

               ZMIN = dabs(xpr)
               zmax = 0.625d0

               mz = 2
               az(1) = zmin
               if(zmin.le.0.24d0) then
                  az(2) = 0.25d0
               else
                  az(2) = (zmin+zmax)/2.d0
               endif
               az(3) = zmax
               nz(1) = 24
               nz(2) = 24
               CALL DMPP(MZ,NPOINTZ,NZ,AZ,Z,wWZ)

c               call gauleg(zmin,zmax,z,wwz,50)
c               do iz=1,50

                    
               do iz = 1,npointz

c               call gauleg(zmin,zmax,z,wwz,50)
c
c               do iz = 1,50

                  TILDZ = dma/dm*(z(iz) + xi)

                  xg = xpr/z(iz)
                  xig = (xi)/z(iz)

                  CALL GPDNUCLEON(XG,XIG,D2G,1,dq2,nn,GPDU)
                  CALL GPDNUCLEON(XG,XIG,D2G,2,dq2,nn,GPDD)
                  CALL GPDNUCLEON(XG,XIG,D2G,3,dq2,nn,GPDS)
                  CALL GPDNUCLEON(-XG,XIG,D2G,1,dq2,nn,GPDU1)
                  CALL GPDNUCLEON(-XG,XIG,D2G,2,dq2,nn,GPDD1)
                  CALL GPDNUCLEON(-XG,XIG,D2G,3,dq2,nn,GPDS1)

                  gpdsu= dks*gpds
                  gpdsu1 = dks*gpds1

                  gpdut = gpdu-gpdu1
                  gpddt = gpdd-gpdd1
                  gpdst = gpds-gpds1
                  gpdsut = gpdsu-gpdsu1

                  if(nn.eq.1)then
                     GPDT=2.d0*((4.D0/9.d0*(gpdut+gpdsut)+1.d0/9.d0*
     *                    (gpddt+gpdsut)+1.d0/9.d0*gpdst))
                  else

                     GPDT=2.d0*(( 1.D0/9.d0*(gpdut+gpdsut)+4.d0/9.d0*
     *                       (gpddt+ gpdsut)+1.d0/9.d0*gpdst))
                  endif

                  DLCMDG = 0.d0
                  DLCMDEX = 0.d0

                  do iphi = 1, 20

                     dcosphi = dcos(PHI(iphi))
                     dsinphi = dsqrt(1.d0-dcosphi*dcosphi)

                     do ig = 1, 2

                        if(ig.eq.1) then

                     DNUM = DMA1**2 - DMA*DMA*(1.d0-TILDM*TILDZ/DMA)**2
                     cont1=0.d0

                        else

                    DNUM = DMA1EX**2 - DMA*DMA*(1.d0-TILDM*TILDZ/DMA)**2
                    cont2=0.d0

                        endif

                    PMIN = 0.5d0*DABS(DNUM/(DMA*(1.d0-TILDM*TILDZ/DMA)))

                        if(pmin.gt.8.d0) then

                           cont1 = 0.d0
                           cont2 = 0.d0
                           go to 500


                        endif

c                        CALL GAULEG(PMIN,8.d0,DPG,WW,50)


                        mp = 1
                        ap(1) = pmin
                        ap(2) = 8.d0
                        np(1) = 32
                        CALL DMPP(Mp,NPOINTp,Np,Ap,dpg,wW)

                        do ip = 1, npointp

cc                        DO IP = 1 , 50

                           if(ig.eq.1) then

                              dnpg=dnkgr(DPG(IP))
                              DP0  = DMA - DSQRT(DMA1**2 + DPG(IP)**2)
                              DCOST=(TILDM*tildz-dp0)/DPg(ip)
                              DSINT=dsqrt(1.d0-dcost*dcost)
                              if(dabs(dcost).gt.1.d0) then
                         write(6,*) dcost,'Warning:cosine out of range'
                                 stop
                              endif
                              DPDG=DSQRT(DPG(IP)**2 + DVEC2 + 2.d0*
     1                             (  DPG(IP)*DSINT*DCOSPHI*DX +
     1                   DPG(IP)*DSINT*DSINPHI*DY + DPG(IP)*DCOST*D3  ))
                              if(dpdg.ge.8.d0) then
                                 dpdg = 7.99d0
                              endif
                              DNPDG = DNKGR(DPDG)
                              DNTOTG=DSQRT(DNPG*DNPDG)

                         CONT1 = CONT1 + TILDM*DPG(IP)*DNTOTG * WW(IP)


                            else

                               DNPEX = dnktot(DPG(IP))-dnkgr(DPG(IP))
                               DP0 = DMA - DSQRT(DMA1EX**2 + DPG(IP)**2)
                               DCOST=(TILDM*tildz-dp0)/DPg(ip)
                               DSINT=dsqrt(1.d0-dcost*dcost)
                               if(dabs(dcost).gt.1.d0) then
                         write(6,*) 'Warning: cosine out of range',dcost
                                  stop
                               endif
                               DPDG=DSQRT(DPG(IP)**2 + DVEC2 +2.d0*
     1                              (  DPG(IP)*DSINT*DCOSPHI*DX +
     1               DPG(IP)*DSINT*DSINPHI*DY + DPG(IP)*DCOST*D3  ))
                               if(dpdg.ge.8.d0) then
                                  dpdg = 7.99d0
                               endif
                               DNPDEX = dnktot(DPDG)-dnkgr(DPDG)
                               DNTOTEX=DSQRT(DNPEX * DNPDEX)

                              CONT2 = CONT2 + TILDM*DPDG*DNTOTEX *WW(IP)
c                        write(6,*) dnpex, dnpdex
                            endif
                      ENDDO

                      ENDDO


 500                  DLCMDG = DLCMDG + CONT1 * WPHI(IPHI)
                      DLCMDEX = DLCMDEX + CONT2 * WPHI(IPHI)

                   enddo

                   DLCTOT = (DLCMDG + DLCMDEX)*gpdt/z(iz)
                   dintz= dintz+ dlctot*wwz(iz)

                enddo

               g = dintz

             return
             end      


c     *-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      
      double precision function f(xpr)

      implicit real*8(a-h,o-z)
      common/cose/dma,xi,dks,dma1,dma1ex,dm,ea,xiacheck,dxiamax

      f = g(xpr)/(xpr-xi)

      return
      end

C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
ccccccc
ccccccc GAUSS LEGENDRE INTEGRATION METHOD
ccccccc

C  X1, X2 = estremi di integrazione
C  N= numbers of integration point
C  X(N)   = points for the integration
C  W(N)   = weights for each point

      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(1),W(1)
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
      Z=DCOS(3.141592653589793D0*(DFLOAT(I)-.25D0)/(DFLOAT(N)+.5D0))
1     CONTINUE
      P1=1.D0
      P2=0.D0
      DO 11 J=1,N
       P3=P2
       P2=P1
      P1=((2.D0*DFLOAT(J)-1.D0)*Z*P2-(DFLOAT(J)-1.D0)*P3)/DFLOAT(J)
11    CONTINUE
      PP=DFLOAT(N)*(Z*P1-P2)/(Z*Z-1.D0)
      Z1=Z
      Z=Z1-P1/PP
      IF (ABS(Z-Z1).GT.EPS)GO TO 1
      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END


c*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
cccccc
cccccc LEADING ORDER CROSS SECTION FOR A DVCS PROCESS OFF A SPINLESS TARGET
cccccc

      subroutine COHDVCSXSECTION(xb,dq2,d2g,ebeam,ehad,dlambda,phideg,
     1    dimcff,recff,crosstot)
      implicit real*8(a-h,o-z)
      parameter(pi=dacos(-1.d0),dc=1.d0/197.31d0)
      common/cose/dma,xi,dks,dma1,dma1ex,dm,ea,xiacheck,dxiamax

     
c      CALL IMAGINARYCFF(XB,D2G,DQ2,DIMCFF)
c      CALL REALCFF(XB,DQ2,D2G,RECFF)
    
cccccc
cccccc dlambda is the helicity of the beam; It is =+1 if the spin is aligned with the direction of the lepton tree-momentum.
cccccc
      
cccccc      
cccccc4HE MASS IN GEV
cccccc

      DMA = 2.D0*(0.9383d0+0.9396d0)-28.295d-3
      ia = 4

cccccc
cccccc Parameters to define the cross section in  [nb/GeV^4]
cccccc      
      gevtobarn = (1.D0/DC)**2*10.d0
      alphaelm = 1.d0/(137.04d0)
      grtorad = 2.d0*pi/360.d0

ccccc
ccccc PARAMETRIZATION OF THE 4HE FORM FACTOR FROM C. R. Ottermann et al.
ccccc  Nucl. Phys. A436, 688 (1985).
ccccc

      a = 0.316d0
      b = 0.681d0
      rq2 = -dc*dc*d2g*1.d6
      base = (1.d0-(a*a*rq2)**6)
      esp = -b*b*rq2
      ff4he = 2.d0*base*dexp(esp)
    

cccccc
cccccc x ranges in [0,1] (it's x_A); x_b is the experiemental x_Bjorken ranging in [0,A=4]
cccccc
      
      q = dsqrt(dq2)
      x  = xb/dfloat(ia)
      eps = 2.d0*dma*x/q
      s = dma*dma +2.d0*ebeam*(ehad+dsqrt(ehad**2-dma**2))
      y = dq2/(s-dma**2)/x

      dkincros = x*y**2*alphaelm**3/(q**4*dsqrt(1.d0+eps**2)*
     -      16.d0*pi**2)


cccccc
c THE RELATION BEETWEEN THE EXPERIMENTAL AND THE THEORETICAL AZIMUTHAL ANGLE IS PHI_EXP = 180° - PHI_THE (SEE E.G. arXiv:1810.02110v1 [hep-ex] 2018 AND PRD 79, 014017 (2009)). THE FOLLOWING HARMONICS ARE FUNCTIONS OF PHI_THE.
cccccc

      phi = phideg*grtorad

cccccc
cccccc Formulas from Ref. A. V. Belitsky and D. Müller, Phys. Rev. D 79, 014017 (2009).
cccccc

      d2gmin = -q**2*( 2.D0*(1.D0-X)*(1.D0-DSQRT(1.D0+eps**2))+eps**2) /
     /          (4.D0*X*(1.D0-X)+eps**2)


      if(d2g.gt.d2gmin) write(6,*)'prob'
      
      delt= (D2G-D2GMIN)/q**2

      dk2=-DELT*(1.D0-X)*(1.D0-y-(y**2*eps**2)/4.D0 ) *
     *        ( DSQRT(1.D0+eps**2) + DELT*( 4.D0*X*(1.D0-X) + eps**2 )/
     /        ( 4.D0*(1-X)) )

      DK = DSQRT(DK2)

cc      write(6,*) dk2,dk

      dj = ( 1.D0 - y - (y*eps**2 )/2.D0 )*(1.D0+D2G/DQ2)-
     -     (1.D0-X)*(2.D0-y)*D2G/DQ2

      P1 = - ( DJ + 2.D0 * DK * DCOS(PHI) )/(y*(1.D0+eps**2) )


      P2 = 1.D0 + D2G/DQ2 - P1
      propa= p1*p2
      dkt =DK*Q/DSQRT( 1.D0-y+eps**2*y**2/4.D0 )

CCCCCC
C THE FOLLOWING AMPLITUDES ARE EXPRESSED IN GEV^(-2)
CCCCCC

cccccc BETHE HEITLER AMPLITUDE

      c0bh= ( ( (2.D0 - y)**2 + y*y*(1.D0+eps**2)**2 )*
     *       (eps**2*DQ2/D2G + 4.D0*(1.D0-X) +(4.D0*X+eps**2)*D2G/DQ2)+
     +       2.D0*eps**2*(4.D0*(1.D0-y)*(3.D0+2.D0*eps**2)+
     +       y*y*(2.D0-eps**4) ) - 4.D0*X*X*(2.D0-y)**2*
     *       (2.D0+eps**2)*D2G/DQ2 + 8.D0*DK2*eps**2*DQ2/D2G )*ff4he**2

      c1bh= -8.D0*(2.D0-y)*DK*(2.D0*X+eps**2-(eps**2*DQ2)/D2G)*ff4he**2
      c2bh= 8.D0*DK2*eps**2*DQ2*ff4he**2/D2G

      tbh = 1.d0/((1.d0+eps**2)**2*x**2*y**2*d2g*propa)*
     * (c0bh+c1bh*dcos(phi)+
     +     c2bh*dcos(2.d0*phi) )

       
cccccc DVCS AMPLITUDE

      c0dvcs=2.d0*(2.d0-2.d0*y +y**2+eps**2*y**2/2.d0)/(1.d0+eps**2)*
     *     (dimcff**2 + recff**2)
      c1dvcs=0.
      c2dvcs=0.
      s1dvcs=0.
      s2dvcs=0.

      tdvcs = 1.d0/(y**2*Q**2)*( c0dvcs+c1dvcs*dcos(phi)+
     +     c2dvcs*dcos(2.d0*phi) + dlambda*(s1dvcs*dsin(phi) +
     +     s2dvcs*dsin(2.d0*phi)) )

cccccc INTERFERENCE BH-DVCS AMPLITUDE

	c0pp = -(4.D0*(2.D0-y)*(1.D0+DSQRT(1.d0 + eps**2)))/(1.D0+eps**2)**2*
     *   ((DKT*DKT*(2.D0-y)**2)
     /   /(q**2*DSQRT(1.d0+eps**2))+D2G/q**2*(1.D0-y-eps**2*y*y/4.D0)*
     *   (2.D0-X)*(1.D0+( 2.D0*X*(2.D0-X + eps**2/(2.D0*X)+
     +   (DSQRT(1.d0+eps**2)-1.D0)/2.D0) + eps**2)/( (2.D0-X)*
     +   ( 1.D0+DSQRT(1.D0+eps**2) ) )) )

        c0bhdvcs= recff*ff4he*c0pp


	c1pp = ((-16.D0*DK*(1.D0-y+Eps**2*y**2/4.D0)/
     /     (1.D0+eps**2)**(5/2))*((1.D0+(1.D0-X)*
     *     (DSQRT(1.D0+eps**2)-1.D0)/(2.D0*X)+eps**2/(4.D0*X))*
     *      X*D2G/Q**2-(3.D0*eps**2)/4.D0)-
     -     4.D0*dk*(2.D0-2.D0*y+y**2+y**2*eps**2/2.D0)*((1.D0+
     +     DSQRT(1.d0+eps**2)-eps**2)/(1.D0+eps**2)**(5/2)) *
     *       (1.D0-(1.D0-3.D0*X)*D2G/Q**2+(1.D0-DSQRT(1.d0+eps**2)+
     +      3.D0*eps**2)/
     /     (1.D0+DSQRT(1.D0+eps**2)-eps**2)*X*D2G/q**2))


      c1bhdvcs= ff4he*recff*c1pp

      c2bhdvcs=0.d0

      spp=(8.D0*DK*(2.D0-y)*y)/(1.d0+eps**2)*
     *        ( 1.D0+ delt * (1.d0-X+
     +        ( DSQRT(1.d0+eps**2)-1.D0)/2.D0  ) /
     *        (1.D0+eps**2))

      s1bhdvcs=ff4he*dimcff*spp

      s2bhdvcs=0.d0

       tbhdvcs=1.d0/(x*y**3*d2g*propa)*(c0bhdvcs + c1bhdvcs*dcos(phi) +
     +     c2bhdvcs*dcos(2.d0*phi)+dlambda*(s1bhdvcs*dsin(phi) +
     +     s2bhdvcs*dsin(2.d0*phi)) )


CCC   DIFFERENTIAL CROSS SECTION E+4HE-> E+GAMMA+4HE

       cros = dkincros*(tbh+tdvcs+tbhdvcs)*gevtobarn
       crosbh = dkincros*(tbh)*gevtobarn
       crosdvcs = dkincros*(tdvcs)*gevtobarn
       crosbhdvcs= dkincros*(tbhdvcs)*gevtobarn

       crosstot = cros
       
      return
      end



C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
cccccc
cccccc INTERPOLATION FOR THE TOTAL MOMENTUM DISTRIBUTION
cccccc
        double precision function dnktot(dkk)

      implicit real*8 (a-h,o-z)

      dimension dkn(1),dnkn(1)
      common/nk/dnk(100),dnkg(100),dk(100)

      pi = dacos(-1.d0)
      dnor = 1./(2.d0*pi)**3


      dkn(1) = dkk

      call lagm(dnk,1,41,0.d0,0.2d0,dkn,1,1,dnkn,1)

      dnktot = dnkn(1) * dnor

      return
      end

      double precision function dnkgr(dkk)
CCCCCCC
CCCCCCC  INTERPATION FOR THE GROUND MOMENTUM DISTRIBUTION
CCCCCCC
      implicit real*8 (a-h,o-z)

      dimension dkn(1),dnkn(1)
      common/nk/dnk(100),dnkg(100),dk(100)

      pi = dacos(-1.d0)
      dnor = 1./(2.d0*pi)**3

      dkn(1) = dkk

      call lagm(dnkg,1,41,0.d0,0.2d0,dkn,1,1,dnkn,1)

      dnkgr = dnkn(1) * dnor

      return
      end

c*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-
cccccc
cccccc SUBROUTINE FOR INTERPOLATING POINTS GIVEN A NUMERICAL TABLE
cccccc

      SUBROUTINE LAGM(G,J0,JF,R0,H,RIN,IN0,NIN,GIN,I0)

CCCCCCC THE FUCNTION  G IS TABULATED  WITH STEP H FROM G(J0)=G(R0)
CCCCCCC TO G(JF) WITH (JF-J0)> 4 .GIN IS THE INTERPOLATED FUNCTION IN POINTS RINCCC  WHOSE FIRST VALUE IS
CCCCCCC RIN(IN0) CORRESPONDING TO GIN(I0). RIN(IN0),..
CCCCCCC ..,RIN(IN0+NIN-1)MUST RANGE BEETWEEN R0 AND R0+H*(JF-J0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION G(1),RIN(1),GIN(1)
      N=JF-J0+1
      IF(N-6) 1,10,10
    1 WRITE(6,1000)
      STOP
   10 R=R0+H*(N-1)
cc      WRITE(6,*) RIN(IN0),R0,RIN(IN0+NIN-1),R
      IF(RIN(IN0).LT.R0.OR.RIN(IN0+NIN-1).GT.R) GO TO 100
      DO 80 I=1,NIN
      AMR=(RIN(IN0+I-1)-R0)/H
      MR=AMR+1
      MR=AMIN0(MR,N-3)
      IF(MR-3) 20,20,30
   20 P=AMR-2.D0
      G0=G(J0)
      L=J0+1
      GO TO 50
   30 L=MR+J0-2
      P=AMR-DFLOAT(L-J0+1)
      G0=G(L-1)
   50 P3=P-3.D0
      P2=P-2.D0
      P1=P-1.D0
      P4=P+1.D0
      P5=P+2.D0
      P23=P3*P5
      P12=P4*P2
      PP1=P*P1
      GIN(I0+I-1)=(PP1*(.1D0*P12*(P5*G(L+4)-P3*G0)+.5D0*P23*(P2*G(L)-
     1            P4*G(L+3)))+P12*P23*(-P1*G(L+1)+P*G(L+2)))/12.D0
   80 CONTINUE
      RETURN
  100 WRITE(6,1100)
      STOP
C 1000 FORMAT(2X,'N TROPPO PICCOLO IN LAGM')
 1000 FORMAT(2X,'PROBLEMS IN INTERPOLATING SUBROUTINE')
 1100 FORMAT(2X,'POINTS OUT OF RANGE IN LAGM')
      END

C*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

CCCCCC
CCCCCC  SUBROUTINE THAT GIVES THE GOLOSKOKOV MODEL FOR THE GPD H. SEE REF. EPJ J C 50, 829-842 (2007).
CCCCCC

      subroutine gpdnucleon(X,XI,D2G,I,dq2,nn,GPD)

       IMPLICIT REAl*8(A-H,O-Z)
       DIMENSION CONT(3),coe(3,4),HIJ(10),CONTAt(3),HIJI(10),hj(4)
       COMMON/PAR/D,A,B
       COMMON/MATR/C(4)


       CALL COEFF(I,dq2)
       CALL PARAREGGE(I,nn,dq2)

CCCCCCC II=1 DGLAP
CCCCCCC II=2 ERBL

       IF(I.EQ.1.or.i.eq.2) THEN

          if(xi.lt.0.0001d0)then

             IF(X.LE.0.D0)THEN
                WRITE(6,*)'Warning! Negative values for x at \xi=0'
                STOP
             ENDIF

             conta=0.d0

            do J=1,4

           hj(J)= c(j)*x**(DFLOAT(J-1)/2.d0)

          conta = CONTA+ HJ(J)

      ENDDO

      RIS = CONTA*x**(-D)*(1.d0-x)**3

      GPD = RIS*DEXP((B+A*DLOG(1.d0/X))*D2G)

         else
            IF(X.GE.-1.D0.AND.X.LT.-XI)then
               GPD=0.D0
            else

            X1=(X+XI)/(1.D0+XI)
            X2=(X-XI)/(1.D0-XI)
            X3=(X-XI)/(1.D0+XI)
            X4=(-X+XI)/(1.D0+XI)

       DO II=1,2

            IF(II.EQ.1) THEN
               IF(X.GE.XI) THEN
                  CONT(1) = 0.D0
               ELSE
                  CONT(1) =0.D0
                  GO TO 12
               ENDIF
            ELSE
               IF(X.GE.(-XI).AND.X.LT.XI)THEN
                  CONT(2)=0.D0
               ELSE
                  CONT(2)=0.D0
                  GO TO 12
               ENDIF
            ENDIF


            DO J=1,4

           DMIJ= 2.D0 + DFLOAT(J-1)/2.D0 - D-(A*D2G)
         IF (II.EQ.1) THEN
       DGRA = (XI**2-X)*(X1**DMIJ - X2**DMIJ)+DMIJ*XI*(1.D0 -X)*
     *       (X1**DMIJ+x2**DMIJ)

      HIJ(J)=(3.D0*DGAMMA(DMIJ-1.D0)*DGRA*C(J))/(2.D0*XI**3.D0*
     *       DGAMMA(DMIJ+2.D0))
         ELSE
CCCCCCC HVALJ ERBL

        DGRA= XI**2-X + DMIJ*XI*(1.D0-X)
      HIJ(J)=(3.D0*DGAMMA(DMIJ-1.D0)*DGRA*C(J)*X1**DMIJ)/(2.D0*XI**3.D0*
     *         DGAMMA(DMIJ+2.D0))
          ENDIF


        cont(ii)=cont(ii)+hij(j)
      enddo

 12     ENDDO

        GPD=(CONT(1)+CONT(2))*DEXP(B*D2G)
      endif
      endif
      ELSE

          xt=dabs(x)
          if(xi.lt.0.00001d0)then

             conta=0.d0

            do J=1,4

           hj(J)= c(j)*xt**(DFLOAT(J-1)/2.d0)

          conta = CONTA+ HJ(J)

      ENDDO

      RIS = CONTA*xt**(-(D+1.d0))*(1.d0-xt)**5

      GPD = RIS*DEXP((B+A*DLOG(1.d0/Xt))*D2G)
      else
            X1=(XT+XI)/(1.D0+XI)
            X2=(XT-XI)/(1.D0-XI)
            X3=(XT-XI)/(1.D0+XI)
            X4=(-XT+XI)/(1.D0+XI)

          DO III=1,2

            IF(III.EQ.1) THEN
               IF(XT.GE.XI) THEN
                  CONTAT(1) = 0.D0
               ELSE
                  CONTAT(1) =0.D0
                  GO TO 66
               ENDIF
            ELSE
               IF(XT.LT.XI)THEN
                  CONTAT(2)=0.D0
               ELSE
                  CONTAT(2)=0.D0
                  GO TO 66
               ENDIF
            ENDIF
              DO J=1,4

      DMIJ= 3.d0 + dfloat(J-1)/2.D0-D-(A*D2G)-1.d0

          IF(III.EQ.1) THEN
        DGRA =((DMIJ**2+2.D0)*(XI**2-XT)**2-(DMIJ**2-1.D0)*(1.D0-XI**2)*
     -        (XT**2-XI**2))* (X1**DMIJ -X2**DMIJ)+3.D0*DMIJ*XI*
     *        (1.D0-XT)*(XI**2-XT)*(X1**DMIJ+X2**DMIJ)

        HIJI(J) = (15.D0*DGAMMA(DMIJ-2.D0)*DGRA*C(J))/(2.D0*XI**5*
     *         DGAMMA(DMIJ+3.D0))

        else

CCCCCCC  HSEAJ ERBL

        DGRA=(X1**DMIJ)*((DMIJ**2 +2.D0)*(XI**2-XT)**2.D0 + 3.D0*DMIJ*
     *      XI*(1.D0-XT)*(XI**2-XT)-(DMIJ**2-1.D0)*(1.D0-XI**2)*
     *      (XT**2-XI**2))-(X4**DMIJ)*((DMIJ**2+2.D0)*(XI**2+XT)**2.D0 +
     +    3.D0*DMIJ*XI*(1.D0+XT)*(XI**2+XT)-(DMIJ**2-1.D0)*(1.D0-XI**2)*
     *      (XT**2-XI**2))

         HIJI(J) = (15.D0*DGAMMA(DMIJ-2.D0)*DGRA*C(J))/(2.D0*XI**5.D0*
     *       DGAMMA(DMIJ+3.D0))


      ENDIF


        CONTAT(III)= CONTAT(III)+HIJI(J)

        ENDDO


 66    ENDDO

       GPD = (CONTAT(1)+CONTAT(2))*DEXP(B*D2G)
       endif
       if (x.lt.0) GPD= -(CONTAT(1)+CONTAT(2))*DEXP(B*D2G)



       endif
      RETURN
      END

      SUBROUTINE COEFF(I,dq2)
      IMPLICIT REAL*8(A-H,O-Z)
      common/Matr/c(4)

       DL = DLOG(DQ2/4.D0)
       IF(I.EQ.1)GO TO 30
       IF(I.EQ.2)GO TO 50
       IF(I.EQ.3)GO TO 20
30     C(1) = 1.52D0 + 0.248*DL
        C(2) = 2.88D0 - 0.940*DL
        C(3) = -0.095D0*DL
        C(4) = 0.D0
        RETURN
50       C(1) = 0.76D0 + 0.248*DL
         C(2) = 3.11D0 - 1.36*DL
         C(3) = -3.99D0 + 1.15*DL
         C(4) = 0.D0
      RETURN
 20     C(1) = 0.123D0 +0.0003D0*DL
        C(2) = -0.327D0 - 0.004D0*DL
        C(3) = 0.692D0 - 0.068D0*DL
        C(4) = -0.486D0 + 0.038D0*DL
      END

      SUBROUTINE PARAREGGE(I,nn,dq2)
      IMPLICIT REAL*8(A-H,O-Z)

       COMMON/PAR/D,A,B
       COMMON/MATR/C(4)

       if(nn.eq.1)then
          DM = 0.9383d0
       else
          DM = 0.9396D0
       endif

       DL = DLOG(DQ2/4.D0)

       IF(I.EQ.1.OR.I.EQ.2)GO TO 1
       IF(I.EQ.3)GO TO 33
 1       D = 0.48D0
         A = 0.9D0
         B = 0.D0
         RETURN
 33      D = 0.10D0 + 0.06*DL
         A = 0.15D0
         B = 2.58D0+0.25D0*DLOG(DM**2/(DQ2+DM**2))
         RETURN

      END
c     *--*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

        SUBROUTINE DMPP(MX,NPX,NX,AX,X,WX)                                 DGA00010
      IMPLICIT REAL*8 (A-H,O-Z)                                         DGA00020
CCCCCCC                                                                 DGA00030
CCCCCCC              IMPOSTA I PUNTI DI INTEGRAZIONE E I PESI           DGA00040
CCCCCCC                        NELLA VARIABILE X                        DGA00050
CCCCCCC                                                                 DGA00060
CCCCCCC              MX = NUMERO DEI SOTTOINTERVALLI                    DGA00070
CCCCCCC              NX(MX) = NUMERO DEI PUNTI NEI SOTTOINTERVALLI      DGA00080
CCCCCCC              AX(MX+1) = ESTREMI DEI SOTTOINTERVALLI             DGA00090
CCCCCCC              X = VARIABILE DI INTEGRAZIONE                      DGA00100
CCCCCCC              WX = PESI NELLA VARIABILE X                        DGA00110
CCCCCCC                                                                 DGA00120
      DIMENSION NX(1),AX(1),X(1),WX(1)                                  DGA00130
      DIMENSION POINT(60),WEIGHT(60)                                    DGA00140
      DATA Z/-1.D0/                                                     DGA00150
      IF(MX.GT.0) GO TO 5                                               DGA00160
      WRITE(6,1000)                                                     DGA00170
 1000 FORMAT(//,10X,'WARNING FROM DMP: MX MUST BE POSITIVE',//)         DGA00180
      CALL EXIT                                                         DGA00190
    5 L0=0                                                              DGA00200
      DO 30 I=1,MX                                                      DGA00210
      NPX=NX(I)                                                         DGA00220
      IF(NPX.LE.0) GO TO 30                                             DGA00230
      CALL DVKN(NPX)                                                    DGA00240
      CALL GPT(POINT,NPX)                                               DGA00250
      CALL GWT(WEIGHT,NPX)                                              DGA00260
      DF=0.5D0*(AX(I+1)+Z*AX(I))                                        DGA00270

      IF(DABS(DF).GE.1.D-08) GO TO 10                                   DGA00280
c     write(6,*) DF
c      WRITE(26,2000)                                                   DGA00290
c2000 FORMAT(//,10X,'WARNING FROM DMP: VALUE OF |DF| NOT ALLOWED',//)   DGA00300
cc     CALL EXIT                                                        DGA00310
c      write(26,*) df
      DS=0.5D0*(AX(I+1)+AX(I))                                          DGA00320
      DO 22 L=1,NPX                                                     DGA00330
      WX(L0+L)=0.d0                                                     DGA00340
   22 X(L0+L)=DS+DF*POINT(L)                                            DGA00350
      L0=L0+NPX                                                         DGA00360
      return

   10 DS=0.5D0*(AX(I+1)+AX(I))                                          DGA00320
      DO 20 L=1,NPX                                                     DGA00330
      WX(L0+L)=DF*WEIGHT(L)                                             DGA00340
   20 X(L0+L)=DS+DF*POINT(L)                                            DGA00350
      L0=L0+NPX                                                         DGA00360
   30 CONTINUE                                                          DGA00370
      NPX=L0                                                            DGA00380
      IF(L0.GT.0) RETURN                                                DGA00390
      WRITE(6,3000)                                                     DGA00400
 3000 FORMAT(//,10X,'WARNING FROM DMP: VERIFY INTEGRATION POINTS',//)   DGA00410
      CALL EXIT                                                         DGA00420
      RETURN                                                            DGA00430
      END                                                               DGA00440
      SUBROUTINE DVKN(N)                                                DGA00450
CCCCCCC                                                                 DGA00460
CCCCCCC                 CHECK SUL VALORE DI N                           DGA00470
CCCCCCC       (N=4,8,12,16,20,24,28,32,36,40,44,48,52,56,60)            DGA00480
CCCCCCC                                                                 DGA00490
      M=N/4                                                             DGA00500
      IF(M.LE.0) GO TO 1                                                DGA00510
      IF((N-4*M).NE.0) GO TO 1                                          DGA00520
      IF(M.GT.15) GO TO 1                                               DGA00530
      RETURN                                                            DGA00540
    1 WRITE(6,1000)                                                     DGA00550
 1000 FORMAT(//,10X,'WARNING FROM DVKN: VALUE OF N NOT ALLOWED',//)     DGA00560
      CALL EXIT                                                         DGA00570
      RETURN                                                            DGA00580
      END                                                               DGA00590
      SUBROUTINE GPT(P,N)                                               DGA01250
CCCCCCC                                                                 DGA01260
CCCCCCC        IMPOSTA I PUNTI DI GAUSS-LEGENDRE NEL VETTORE P          DGA01270
CCCCCCC        (N=4,8,12,16,20,24,28,32,36,40,44,48,52,56,60)           DGA01280
CCCCCCC                        ( -1<P<1 )                               DGA01290
CCCCCCC                                                                 DGA01300
      DOUBLE PRECISION P(60)                                            DGA01310
      DATA JP/-1/,Z/-1.D0/                                              DGA01320
      L=N/4                                                             DGA01330
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),L                     DGA01340
    1 CONTINUE                                                          DGA01350
      P(1)=0.8611363115940526D0                                         DGA01360
      P(2)=0.3399810435848563D0                                         DGA01370
      GO TO 40                                                          DGA01380
    2 CONTINUE                                                          DGA01390
      P(1)=0.9602898564975362D0                                         DGA01400
      P(2)=0.7966664774136267D0                                         DGA01410
      P(3)=0.5255324099163290D0                                         DGA01420
      P(4)=0.1834346424956498D0                                         DGA01430
      GO TO 40                                                          DGA01440
    3 CONTINUE                                                          DGA01450
      P(1)=0.9815606342467193D0                                         DGA01460
      P(2)=0.9041172563704749D0                                         DGA01470
      P(3)=0.7699026741943047D0                                         DGA01480
      P(4)=0.5873179542866175D0                                         DGA01490
      P(5)=0.3678314989981802D0                                         DGA01500
      P(6)=0.1252334085114689D0                                         DGA01510
      GO TO 40                                                          DGA01520
    4 CONTINUE                                                          DGA01530
      P(1)=0.9894009349916499D0                                         DGA01540
      P(2)=0.9445750230732326D0                                         DGA01550
      P(3)=0.8656312023878317D0                                         DGA01560
      P(4)=0.7554044083550030D0                                         DGA01570
      P(5)=0.6178762444026438D0                                         DGA01580
      P(6)=0.4580167776572274D0                                         DGA01590
      P(7)=0.2816035507792589D0                                         DGA01600
      P(8)=0.9501250983763744D-01                                       DGA01610
      GO TO 40                                                          DGA01620
    5 CONTINUE                                                          DGA01630
      P(1)=0.9931285991850949D0                                         DGA01640
      P(2)=0.9639719272779138D0                                         DGA01650
      P(3)=0.9122344282513259D0                                         DGA01660
      P(4)=0.8391169718222188D0                                         DGA01670
      P(5)=0.7463319064601508D0                                         DGA01680
      P(6)=0.6360536807265150D0                                         DGA01690
      P(7)=0.5108670019508271D0                                         DGA01700
      P(8)=0.3737060887154196D0                                         DGA01710
      P(9)=0.2277858511416451D0                                         DGA01720
      P(10)=0.7652652113349733D-01                                      DGA01730
      GO TO 40                                                          DGA01740
    6 CONTINUE                                                          DGA01750
      P(1)=0.9951872199970214D0                                         DGA01760
      P(2)=0.9747285559713095D0                                         DGA01770
      P(3)=0.9382745520027328D0                                         DGA01780
      P(4)=0.8864155270044010D0                                         DGA01790
      P(5)=0.8200019859739029D0                                         DGA01800
      P(6)=0.7401241915785544D0                                         DGA01810
      P(7)=0.6480936519369756D0                                         DGA01820
      P(8)=0.5454214713888395D0                                         DGA01830
      P(9)=0.4337935076260451D0                                         DGA01840
      P(10)=0.3150426796961634D0                                        DGA01850
      P(11)=0.1911188674736163D0                                        DGA01860
      P(12)=0.6405689286260563D-01                                      DGA01870
      GO TO 40                                                          DGA01880
    7 CONTINUE                                                          DGA01890
      P(1)=0.9964424975739544D0                                         DGA01900
      P(2)=0.9813031653708728D0                                         DGA01910
      P(3)=0.9542592806289382D0                                         DGA01920
      P(4)=0.9156330263921321D0                                         DGA01930
      P(5)=0.8658925225743950D0                                         DGA01940
      P(6)=0.8056413709171792D0                                         DGA01950
      P(7)=0.7356108780136318D0                                         DGA01960
      P(8)=0.6566510940388650D0                                         DGA01970
      P(9)=0.5697204718114017D0                                         DGA01980
      P(10)=0.4758742249551183D0                                        DGA01990
      P(11)=0.3762515160890787D0                                        DGA02000
      P(12)=0.2720616276351781D0                                        DGA02010
      P(13)=0.1645692821333808D0                                        DGA02020
      P(14)=0.5507928988403427D-01                                      DGA02030
      GO TO 40                                                          DGA02040
    8 CONTINUE                                                          DGA02050
      P(1)=0.9972638618494816D0                                         DGA02060
      P(2)=0.9856115115452683D0                                         DGA02070
      P(3)=0.9647622555875064D0                                         DGA02080
      P(4)=0.9349060759377397D0                                         DGA02090
      P(5)=0.8963211557660521D0                                         DGA02100
      P(6)=0.8493676137325700D0                                         DGA02110
      P(7)=0.7944837959679424D0                                         DGA02120
      P(8)=0.7321821187402897D0                                         DGA02130
      P(9)=0.6630442669302152D0                                         DGA02140
      P(10)=0.5877157572407623D0                                        DGA02150
      P(11)=0.5068999089322294D0                                        DGA02160
      P(12)=0.4213512761306353D0                                        DGA02170
      P(13)=0.3318686022821276D0                                        DGA02180
      P(14)=0.2392873622521371D0                                        DGA02190
      P(15)=0.1444719615827965D0                                        DGA02200
      P(16)=0.4830766568773832D-01                                      DGA02210
      GO TO 40                                                          DGA02220
    9 CONTINUE                                                          DGA02230
      P(1)=0.9978304624840858D0                                         DGA02240
      P(2)=0.9885864789022122D0                                         DGA02250
      P(3)=0.9720276910496980D0                                         DGA02260
      P(4)=0.9482729843995076D0                                         DGA02270
      P(5)=0.9174977745156591D0                                         DGA02280
      P(6)=0.8799298008903971D0                                         DGA02290
      P(7)=0.8358471669924753D0                                         DGA02300
      P(8)=0.7855762301322065D0                                         DGA02310
      P(9)=0.7294891715935566D0                                         DGA02320
      P(10)=0.6680012365855211D0                                        DGA02330
      P(11)=0.6015676581359805D0                                        DGA02340
      P(12)=0.5306802859262452D0                                        DGA02350
      P(13)=0.4558639444334203D0                                        DGA02360
      P(14)=0.3776725471196892D0                                        DGA02370
      P(15)=0.2966849953440283D0                                        DGA02380
      P(16)=0.2135008923168656D0                                        DGA02390
      P(17)=0.1287361038093848D0                                        DGA02400
      P(18)=0.4301819847370861D-01                                      DGA02410
      GO TO 40                                                          DGA02420
   10 CONTINUE                                                          DGA02430
      P(1)=0.9982377097105592D0                                         DGA02440
      P(2)=0.9907262386994570D0                                         DGA02450
      P(3)=0.9772599499837743D0                                         DGA02460
      P(4)=0.9579168192137917D0                                         DGA02470
      P(5)=0.9328128082786765D0                                         DGA02480
      P(6)=0.9020988069688743D0                                         DGA02490
      P(7)=0.8659595032122595D0                                         DGA02500
      P(8)=0.8246122308333117D0                                         DGA02510
      P(9)=0.7783056514265194D0                                         DGA02520
      P(10)=0.7273182551899271D0                                        DGA02530
      P(11)=0.6719566846141795D0                                        DGA02540
      P(12)=0.6125538896679802D0                                        DGA02550
      P(13)=0.5494671250951282D0                                        DGA02560
      P(14)=0.4830758016861787D0                                        DGA02570
      P(15)=0.4137792043716050D0                                        DGA02580
      P(16)=0.3419940908257585D0                                        DGA02590
      P(17)=0.2681521850072537D0                                        DGA02600
      P(18)=0.1926975807013711D0                                        DGA02610
      P(19)=0.1160840706752552D0                                        DGA02620
      P(20)=0.3877241750605082D-01                                      DGA02630
      GO TO 40                                                          DGA02640
   11 CONTINUE                                                          DGA02650
      P(1)=0.9985402006367742D0                                         DGA02660
      P(2)=0.9923163921385158D0                                         DGA02670
      P(3)=0.9811518330779140D0                                         DGA02680
      P(4)=0.9650996504224931D0                                         DGA02690
      P(5)=0.9442395091181941D0                                         DGA02700
      P(6)=0.9186752599841758D0                                         DGA02710
      P(7)=0.8885342382860432D0                                         DGA02720
      P(8)=0.8539665950047104D0                                         DGA02730
      P(9)=0.8151445396451350D0                                         DGA02740
      P(10)=0.7722614792487559D0                                        DGA02750
      P(11)=0.7255310536607170D0                                        DGA02760
      P(12)=0.6751860706661224D0                                        DGA02770
      P(13)=0.6214773459035758D0                                        DGA02780
      P(14)=0.5646724531854708D0                                        DGA02790
      P(15)=0.5050543913882023D0                                        DGA02800
      P(16)=0.4429201745254115D0                                        DGA02810
      P(17)=0.3785793520147071D0                                        DGA02820
      P(18)=0.3123524665027858D0                                        DGA02830
      P(19)=0.2445694569282013D0                                        DGA02840
      P(20)=0.1755680147755168D0                                        DGA02850
      P(21)=0.1056919017086532D0                                        DGA02860
      P(22)=0.3528923696413536D-01                                      DGA02870
      GO TO 40                                                          DGA02880
   12 CONTINUE                                                          DGA02890
      P(1)=0.9987710072524261D0                                         DGA02900
      P(2)=0.9935301722663508D0                                         DGA02910
      P(3)=0.9841245837228269D0                                         DGA02920
      P(4)=0.9705915925462473D0                                         DGA02930
      P(5)=0.9529877031604309D0                                         DGA02940
      P(6)=0.9313866907065543D0                                         DGA02950
      P(7)=0.9058791367155697D0                                         DGA02960
      P(8)=0.8765720202742479D0                                         DGA02970
      P(9)=0.8435882616243935D0                                         DGA02980
      P(10)=0.8070662040294426D0                                        DGA02990
      P(11)=0.7671590325157403D0                                        DGA03000
      P(12)=0.7240341309238147D0                                        DGA03010
      P(13)=0.6778723796326639D0                                        DGA03020
      P(14)=0.6288673967765136D0                                        DGA03030
      P(15)=0.5772247260839727D0                                        DGA03040
      P(16)=0.5231609747222330D0                                        DGA03050
      P(17)=0.4669029047509584D0                                        DGA03060
      P(18)=0.4086864819907167D0                                        DGA03070
      P(19)=0.3487558862921607D0                                        DGA03080
      P(20)=0.2873624873554556D0                                        DGA03090
      P(21)=0.2247637903946891D0                                        DGA03100
      P(22)=0.1612223560688917D0                                        DGA03110
      P(23)=0.9700469920946270D-01                                      DGA03120
      P(24)=0.3238017096286936D-01                                      DGA03130
      GO TO 40                                                          DGA03140
   13 CONTINUE                                                          DGA03150
      P(1)=0.9989511111039503D+00                                       DGA03160
      P(2)=0.9944775909292160D+00                                       DGA03170
      P(3)=0.9864461956515498D+00                                       DGA03180
      P(4)=0.9748838842217445D+00                                       DGA03190
      P(5)=0.9598318269330866D+00                                       DGA03200
      P(6)=0.9413438536413591D+00                                       DGA03210
      P(7)=0.9194861289164245D+00                                       DGA03220
      P(8)=0.8943368905344953D+00                                       DGA03230
      P(9)=0.8659861628460676D+00                                       DGA03240
      P(10)=0.8345354323267345D+00                                      DGA03250
      P(11)=0.8000972834304683D+00                                      DGA03260
      P(12)=0.7627949951937450D+00                                      DGA03270
      P(13)=0.7227620997499832D+00                                      DGA03280
      P(14)=0.6801419042271677D+00                                      DGA03290
      P(15)=0.6350869776952459D+00                                      DGA03300
      P(16)=0.5877586049795791D+00                                      DGA03310
      P(17)=0.5383262092858274D+00                                      DGA03320
      P(18)=0.4869667456980961D+00                                      DGA03330
      P(19)=0.4338640677187617D+00                                      DGA03340
      P(20)=0.3792082691160937D+00                                      DGA03350
      P(21)=0.3231950034348078D+00                                      DGA03360
      P(22)=0.2660247836050018D+00                                      DGA03370
      P(23)=0.2079022641563661D+00                                      DGA03380
      P(24)=0.1490355086069492D+00                                      DGA03390
      P(25)=0.8963524464890057D-01                                      DGA03400
      P(26)=0.2991410979733877D-01                                      DGA03410
      GO TO 40                                                          DGA03420
   14 CONTINUE                                                          DGA03430
      P(1)=0.9990943438014656D+00                                       DGA03440
      P(2)=0.9952312260810697D+00                                       DGA03450
      P(3)=0.9882937155401615D+00                                       DGA03460
      P(4)=0.9783017091402564D+00                                       DGA03470
      P(5)=0.9652859019054902D+00                                       DGA03480
      P(6)=0.9492864795619626D+00                                       DGA03490
      P(7)=0.9303528802474963D+00                                       DGA03500
      P(8)=0.9085436204206555D+00                                       DGA03510
      P(9)=0.8839261083278275D+00                                       DGA03520
      P(10)=0.8565764337627486D+00                                      DGA03530
      P(11)=0.8265791321428817D+00                                      DGA03540
      P(12)=0.7940269228938665D+00                                      DGA03550
      P(13)=0.7590204227051289D+00                                      DGA03560
      P(14)=0.7216678344501881D+00                                      DGA03570
      P(15)=0.6820846126944705D+00                                      DGA03580
      P(16)=0.6403931068070069D+00                                      DGA03590
      P(17)=0.5967221827706633D+00                                      DGA03600
      P(18)=0.5512068248555346D+00                                      DGA03610
      P(19)=0.5039877183843817D+00                                      DGA03620
      P(20)=0.4552108148784596D+00                                      DGA03630
      P(21)=0.4050268809270913D+00                                      DGA03640
      P(22)=0.3535910321749545D+00                                      DGA03650
      P(23)=0.3010622538672207D+00                                      DGA03660
      P(24)=0.2476029094343372D+00                                      DGA03670
      P(25)=0.1933782386352753D+00                                      DGA03680
      P(26)=0.1385558468103762D+00                                      DGA03690
      P(27)=0.8330518682243537D-01                                      DGA03700
      P(28)=0.2779703528727544D-01                                      DGA03710
      GO TO 40                                                          DGA03720
   15 CONTINUE                                                          DGA03730
      P(1)=0.9992101232274360D+00                                       DGA03740
      P(2)=0.9958405251188382D+00                                       DGA03750
      P(3)=0.9897878952222217D+00                                       DGA03760
      P(4)=0.9810672017525982D+00                                       DGA03770
      P(5)=0.9697017887650527D+00                                       DGA03780
      P(6)=0.9557222558399961D+00                                       DGA03790
      P(7)=0.9391662761164232D+00                                       DGA03800
      P(8)=0.9200784761776276D+00                                       DGA03810
      P(9)=0.8985103108100459D+00                                       DGA03820
      P(10)=0.8745199226468983D+00                                      DGA03830
      P(11)=0.8481719847859296D+00                                      DGA03840
      P(12)=0.8195375261621458D+00                                      DGA03850
      P(13)=0.7886937399322641D+00                                      DGA03860
      P(14)=0.7557237753065857D+00                                      DGA03870
      P(15)=0.7207165133557304D+00                                      DGA03880
      P(16)=0.6837663273813554D+00                                      DGA03890
      P(17)=0.6449728284894771D+00                                      DGA03900
      P(18)=0.6044405970485104D+00                                      DGA03910
      P(19)=0.5622789007539445D+00                                      DGA03920
      P(20)=0.5186014000585697D+00                                      DGA03930
      P(21)=0.4735258417617071D+00                                      DGA03940
      P(22)=0.4271737415830784D+00                                      DGA03950
      P(23)=0.3796700565767980D+00                                      DGA03960
      P(24)=0.3311428482684482D+00                                      DGA03970
      P(25)=0.2817229374232617D+00                                      DGA03980
      P(26)=0.2315435513760293D+00                                      DGA03990
      P(27)=0.1807399648734254D+00                                      DGA04000
      P(28)=0.1294491353969450D+00                                      DGA04010
      P(29)=0.7780933394953657D-01                                      DGA04020
      P(30)=0.2595977230124780D-01                                      DGA04030
   40 M=N/2+1                                                           DGA04040
      J=M                                                               DGA04050
      DO 50 I=M,N                                                       DGA04060
      J=J+JP                                                            DGA04070
      P(I)=P(J)                                                         DGA04080
   50 P(J)=Z*P(J)                                                       DGA04090
      RETURN                                                            DGA04100
      END                                                               DGA04110
      SUBROUTINE GWT(W,N)                                               DGA04340
CCCCCCC                                                                 DGA04350
CCCCCCC        IMPOSTA I PESI DI GAUSS-LEGENDRE NEL VETTORE W           DGA04360
CCCCCCC        (N=4,8,12,16,20,24,28,32,36,40,44,48,52,56,60)           DGA04370
CCCCCCC                                                                 DGA04380
      DOUBLE PRECISION W(60)                                            DGA04390
      DATA JP/-1/                                                       DGA04400
      L=N/4                                                             DGA04410
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),L                     DGA04420
    1 CONTINUE                                                          DGA04430
      W(1)=0.3478548451374539D0                                         DGA04440
      W(2)=0.6521451548625461D0                                         DGA04450
      GO TO 40                                                          DGA04460
    2 CONTINUE                                                          DGA04470
      W(1)=0.1012285362903763D0                                         DGA04480
      W(2)=0.2223810344533745D0                                         DGA04490
      W(3)=0.3137066458778873D0                                         DGA04500
      W(4)=0.3626837833783620D0                                         DGA04510
      GO TO 40                                                          DGA04520
    3 CONTINUE                                                          DGA04530
      W(1)=0.4717533638651183D-01                                       DGA04540
      W(2)=0.1069393259953184D0                                         DGA04550
      W(3)=0.1600783285433462D0                                         DGA04560
      W(4)=0.2031674267230659D0                                         DGA04570
      W(5)=0.2334925365383548D0                                         DGA04580
      W(6)=0.2491470458134028D0                                         DGA04590
      GO TO 40                                                          DGA04600
    4 CONTINUE                                                          DGA04610
      W(1)=0.2715245941175410D-01                                       DGA04620
      W(2)=0.6225352393864789D-01                                       DGA04630
      W(3)=0.9515851168249279D-01                                       DGA04640
      W(4)=0.1246289712555339D0                                         DGA04650
      W(5)=0.1495959888165767D0                                         DGA04660
      W(6)=0.1691565193950025D0                                         DGA04670
      W(7)=0.1826034150449236D0                                         DGA04680
      W(8)=0.1894506104550685D0                                         DGA04690
      GO TO 40                                                          DGA04700
    5 CONTINUE                                                          DGA04710
      W(1)=0.1761400713915212D-01                                       DGA04720
      W(2)=0.4060142980038694D-01                                       DGA04730
      W(3)=0.6267204833410906D-01                                       DGA04740
      W(4)=0.8327674157670475D-01                                       DGA04750
      W(5)=0.1019301198172404D0                                         DGA04760
      W(6)=0.1181945319615184D0                                         DGA04770
      W(7)=0.1316886384491766D0                                         DGA04780
      W(8)=0.1420961093183821D0                                         DGA04790
      W(9)=0.1491729864726037D0                                         DGA04800
      W(10)=0.1527533871307259D0                                        DGA04810
      GO TO 40                                                          DGA04820
    6 CONTINUE                                                          DGA04830
      W(1)=0.1234122979998720D-01                                       DGA04840
      W(2)=0.2853138862893367D-01                                       DGA04850
      W(3)=0.4427743881741981D-01                                       DGA04860
      W(4)=0.5929858491543678D-01                                       DGA04870
      W(5)=0.7334648141108031D-01                                       DGA04880
      W(6)=0.8619016153195328D-01                                       DGA04890
      W(7)=0.9761865210411389D-01                                       DGA04900
      W(8)=0.1074442701159656D0                                         DGA04910
      W(9)=0.1155056680537256D0                                         DGA04920
      W(10)=0.1216704729278034D0                                        DGA04930
      W(11)=0.1258374563468283D0                                        DGA04940
      W(12)=0.1279381953467522D0                                        DGA04950
      GO TO 40                                                          DGA04960
    7 CONTINUE                                                          DGA04970
      W(1)=0.9124282593094518D-02                                       DGA04980
      W(2)=0.2113211259277126D-01                                       DGA04990
      W(3)=0.3290142778230438D-01                                       DGA05000
      W(4)=0.4427293475900423D-01                                       DGA05010
      W(5)=0.5510734567571675D-01                                       DGA05020
      W(6)=0.6527292396699960D-01                                       DGA05030
      W(7)=0.7464621423456878D-01                                       DGA05040
      W(8)=0.8311341722890122D-01                                       DGA05050
      W(9)=0.9057174439303284D-01                                       DGA05060
      W(10)=0.9693065799792992D-01                                      DGA05070
      W(11)=0.1021129675780608D0                                        DGA05080
      W(12)=0.1060557659228464D0                                        DGA05090
      W(13)=0.1087111922582941D0                                        DGA05100
      W(14)=0.1100470130164752D0                                        DGA05110
      GO TO 40                                                          DGA05120
    8 CONTINUE                                                          DGA05130
      W(1)=0.7018610009470097D-02                                       DGA05140
      W(2)=0.1627439473090567D-01                                       DGA05150
      W(3)=0.2539206530926206D-01                                       DGA05160
      W(4)=0.3427386291302143D-01                                       DGA05170
      W(5)=0.4283589802222668D-01                                       DGA05180
      W(6)=0.5099805926237618D-01                                       DGA05190
      W(7)=0.5868409347853555D-01                                       DGA05200
      W(8)=0.6582222277636185D-01                                       DGA05210
      W(9)=0.7234579410884851D-01                                       DGA05220
      W(10)=0.7819389578707031D-01                                      DGA05230
      W(11)=0.8331192422694676D-01                                      DGA05240
      W(12)=0.8765209300440381D-01                                      DGA05250
      W(13)=0.9117387869576388D-01                                      DGA05260
      W(14)=0.9384439908080457D-01                                      DGA05270
      W(15)=0.9563872007927486D-01                                      DGA05280
      W(16)=0.9654008851472780D-01                                      DGA05290
      GO TO 40                                                          DGA05300
    9 CONTINUE                                                          DGA05310
      W(1)=0.5565719664245045D-02                                       DGA05320
      W(2)=0.1291594728406557D-01                                       DGA05330
      W(3)=0.2018151529773547D-01                                       DGA05340
      W(4)=0.2729862149856878D-01                                       DGA05350
      W(5)=0.3421381077030723D-01                                       DGA05360
      W(6)=0.4087575092364490D-01                                       DGA05370
      W(7)=0.4723508349026598D-01                                       DGA05380
      W(8)=0.5324471397775992D-01                                       DGA05390
      W(9)=0.5886014424532482D-01                                       DGA05400
      W(10)=0.6403979735501549D-01                                      DGA05410
      W(11)=0.6874532383573644D-01                                      DGA05420
      W(12)=0.7294188500565306D-01                                      DGA05430
      W(13)=0.7659841064587068D-01                                      DGA05440
      W(14)=0.7968782891207160D-01                                      DGA05450
      W(15)=0.8218726670433971D-01                                      DGA05460
      W(16)=0.8407821897966194D-01                                      DGA05470
      W(17)=0.8534668573933863D-01                                      DGA05480
      W(18)=0.8598327567039475D-01                                      DGA05490
      GO TO 40                                                          DGA05500
   10 CONTINUE                                                          DGA05510
      W(1)=0.4521277098533191D-02                                       DGA05520
      W(2)=0.1049828453115281D-01                                       DGA05530
      W(3)=0.1642105838190789D-01                                       DGA05540
      W(4)=0.2224584919416696D-01                                       DGA05550
      W(5)=0.2793700698002340D-01                                       DGA05560
      W(6)=0.3346019528254785D-01                                       DGA05570
      W(7)=0.3878216797447202D-01                                       DGA05580
      W(8)=0.4387090818567327D-01                                       DGA05590
      W(9)=0.4869580763507223D-01                                       DGA05600
      W(10)=0.5322784698393682D-01                                      DGA05610
      W(11)=0.5743976909939155D-01                                      DGA05620
      W(12)=0.6130624249292894D-01                                      DGA05630
      W(13)=0.6480401345660104D-01                                      DGA05640
      W(14)=0.6791204581523390D-01                                      DGA05650
      W(15)=0.7061164739128678D-01                                      DGA05660
      W(16)=0.7288658239580406D-01                                      DGA05670
      W(17)=0.7472316905796826D-01                                      DGA05680
      W(18)=0.7611036190062624D-01                                      DGA05690
      W(19)=0.7703981816424797D-01                                      DGA05700
      W(20)=0.7750594797842481D-01                                      DGA05710
      GO TO 40                                                          DGA05720
   11 CONTINUE                                                          DGA05730
      W(1)=0.3745404803112778D-02                                       DGA05740
      W(2)=0.8700481367524844D-02                                       DGA05750
      W(3)=0.1361958675557999D-01                                       DGA05760
      W(4)=0.1847148173681475D-01                                       DGA05770
      W(5)=0.2323148190201921D-01                                       DGA05780
      W(6)=0.2787578282128101D-01                                       DGA05790
      W(7)=0.3238122281206982D-01                                       DGA05800
      W(8)=0.3672534781380887D-01                                       DGA05810
      W(9)=0.4088651231034622D-01                                       DGA05820
      W(10)=0.4484398408197003D-01                                      DGA05830
      W(11)=0.4857804644835204D-01                                      DGA05840
      W(12)=0.5207009609170446D-01                                      DGA05850
      W(13)=0.5530273556372805D-01                                      DGA05860
      W(14)=0.5825985987759550D-01                                      DGA05870
      W(15)=0.6092673670156197D-01                                      DGA05880
      W(16)=0.6329007973320385D-01                                      DGA05890
      W(17)=0.6533811487918143D-01                                      DGA05900
      W(18)=0.6706063890629365D-01                                      DGA05910
      W(19)=0.6844907026936666D-01                                      DGA05920
      W(20)=0.6949649186157258D-01                                      DGA05930
      W(21)=0.7019768547355821D-01                                      DGA05940
      W(22)=0.7054915778935407D-01                                      DGA05950
      GO TO 40                                                          DGA05960
   12 CONTINUE                                                          DGA05970
      W(1)=0.3153346052305839D-02                                       DGA05980
      W(2)=0.7327553901276262D-02                                       DGA05990
      W(3)=0.1147723457923454D-01                                       DGA06000
      W(4)=0.1557931572294385D-01                                       DGA06010
      W(5)=0.1961616045735553D-01                                       DGA06020
      W(6)=0.2357076083932438D-01                                       DGA06030
      W(7)=0.2742650970835695D-01                                       DGA06040
      W(8)=0.3116722783279809D-01                                       DGA06050
      W(9)=0.3477722256477044D-01                                       DGA06060
      W(10)=0.3824135106583071D-01                                      DGA06070
      W(11)=0.4154508294346475D-01                                      DGA06080
      W(12)=0.4467456085669428D-01                                      DGA06090
      W(13)=0.4761665849249047D-01                                      DGA06100
      W(14)=0.5035903555385447D-01                                      DGA06110
      W(15)=0.5289018948519367D-01                                      DGA06120
      W(16)=0.5519950369998416D-01                                      DGA06130
      W(17)=0.5727729210040322D-01                                      DGA06140
      W(18)=0.5911483969839564D-01                                      DGA06150
      W(19)=0.6070443916589388D-01                                      DGA06160
      W(20)=0.6203942315989266D-01                                      DGA06170
      W(21)=0.6311419228625403D-01                                      DGA06180
      W(22)=0.6392423858464819D-01                                      DGA06190
      W(23)=0.6446616443595008D-01                                      DGA06200
      W(24)=0.6473769681268392D-01                                      DGA06210
      GO TO 40                                                          DGA06220
   13 CONTINUE                                                          DGA06230
      W(1)=0.2691316950047111D-02                                       DGA06240
      W(2)=0.6255523962973277D-02                                       DGA06250
      W(3)=0.9802634579462752D-02                                       DGA06260
      W(4)=0.1331511498234096D-01                                       DGA06270
      W(5)=0.1678002339630074D-01                                       DGA06280
      W(6)=0.2018489150798079D-01                                       DGA06290
      W(7)=0.2351751355398446D-01                                       DGA06300
      W(8)=0.2676595374650401D-01                                       DGA06310
      W(9)=0.2991858114714395D-01                                       DGA06320
      W(10)=0.3296410908971880D-01                                      DGA06330
      W(11)=0.3589163483509723D-01                                      DGA06340
      W(12)=0.3869067831042398D-01                                      DGA06350
      W(13)=0.4135121950056027D-01                                      DGA06360
      W(14)=0.4386373425900041D-01                                      DGA06370
      W(15)=0.4621922837278479D-01                                      DGA06380
      W(16)=0.4840926974407490D-01                                      DGA06390
      W(17)=0.5042601856634238D-01                                      DGA06400
      W(18)=0.5226225538390699D-01                                      DGA06410
      W(19)=0.5391140693275726D-01                                      DGA06420
      W(20)=0.5536756966930265D-01                                      DGA06430
      W(21)=0.5662553090236860D-01                                      DGA06440
      W(22)=0.5768078745252683D-01                                      DGA06450
      W(23)=0.5852956177181387D-01                                      DGA06460
      W(24)=0.5916881546604297D-01                                      DGA06470
      W(25)=0.5959626017124816D-01                                      DGA06480
      W(26)=0.5981036574529186D-01                                      DGA06490
      GO TO 40                                                          DGA06500
   14 CONTINUE                                                          DGA06510
      W(1)=0.2323855375773216D-02                                       DGA06520
      W(2)=0.5402522246015338D-02                                       DGA06530
      W(3)=0.8469063163307888D-02                                       DGA06540
      W(4)=0.1150982434038338D-01                                       DGA06550
      W(5)=0.1451508927802147D-01                                       DGA06560
      W(6)=0.1747551291140095D-01                                       DGA06570
      W(7)=0.2038192988240257D-01                                       DGA06580
      W(8)=0.2322535156256532D-01                                       DGA06590
      W(9)=0.2599698705839195D-01                                       DGA06600
      W(10)=0.2868826847382274D-01                                      DGA06610
      W(11)=0.3129087674731045D-01                                      DGA06620
      W(12)=0.3379676711561176D-01                                      DGA06630
      W(13)=0.3619819387231519D-01                                      DGA06640
      W(14)=0.3848773425924766D-01                                      DGA06650
      W(15)=0.4065831138474452D-01                                      DGA06660
      W(16)=0.4270321608466709D-01                                      DGA06670
      W(17)=0.4461612765269228D-01                                      DGA06680
      W(18)=0.4639113337300190D-01                                      DGA06690
      W(19)=0.4802274679360026D-01                                      DGA06700
      W(20)=0.4950592468304758D-01                                      DGA06710
      W(21)=0.5083608261779848D-01                                      DGA06720
      W(22)=0.5200910915174140D-01                                      DGA06730
      W(23)=0.5302137852401076D-01                                      DGA06740
      W(24)=0.5386976186571449D-01                                      DGA06750
      W(25)=0.5455163687088942D-01                                      DGA06760
      W(26)=0.5506489590176243D-01                                      DGA06770
      W(27)=0.5540795250324512D-01                                      DGA06780
      W(28)=0.5557974630651440D-01                                      DGA06790
      GO TO 40                                                          DGA06800
   15 CONTINUE                                                          DGA06810
      W(1)=0.2026811968873758D-02                                       DGA06820
      W(2)=0.4712729926953569D-02                                       DGA06830
      W(3)=0.7389931163345456D-02                                       DGA06840
      W(4)=0.1004755718228798D-01                                       DGA06850
      W(5)=0.1267816647681596D-01                                       DGA06860
      W(6)=0.1527461859678480D-01                                       DGA06870
      W(7)=0.1782990101420772D-01                                       DGA06880
      W(8)=0.2033712072945729D-01                                       DGA06890
      W(9)=0.2278951694399782D-01                                       DGA06900
      W(10)=0.2518047762152125D-01                                      DGA06910
      W(11)=0.2750355674992479D-01                                      DGA06920
      W(12)=0.2975249150078895D-01                                      DGA06930
      W(13)=0.3192121901929633D-01                                      DGA06940
      W(14)=0.3400389272494642D-01                                      DGA06950
      W(15)=0.3599489805108450D-01                                      DGA06960
      W(16)=0.3788886756924344D-01                                      DGA06970
      W(17)=0.3968069545238080D-01                                      DGA06980
      W(18)=0.4136555123558476D-01                                      DGA06990
      W(19)=0.4293889283593564D-01                                      DGA07000
      W(20)=0.4439647879578711D-01                                      DGA07010
      W(21)=0.4573437971611449D-01                                      DGA07020
      W(22)=0.4694898884891220D-01                                      DGA07030
      W(23)=0.4803703181997118D-01                                      DGA07040
      W(24)=0.4899557545575684D-01                                      DGA07050
      W(25)=0.4982203569055018D-01                                      DGA07060
      W(26)=0.5051418453250937D-01                                      DGA07070
      W(27)=0.5107015606985563D-01                                      DGA07080
      W(28)=0.5148845150098093D-01                                      DGA07090
      W(29)=0.5176794317491019D-01                                      DGA07100
      W(30)=0.5190787763122064D-01                                      DGA07110
   40 M=N/2+1                                                           DGA07120
      J=M                                                               DGA07130
      DO 50 I=M,N                                                       DGA07140
      J=J+JP                                                            DGA07150
   50 W(I)=W(J)                                                         DGA07160
      RETURN                                                            DGA07170
      END                                                               DGA07180
 
