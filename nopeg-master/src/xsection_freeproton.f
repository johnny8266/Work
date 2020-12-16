      subroutine xsection_freeproton(xb,d2g,dq2,ebeam,eh,phideg,dlambda,
     1        crosstotfree)
      implicit real*8(a-h,o-z)
      parameter(pi=dacos(-1.d0),alphaelm=1.d0/137.04d0,dm=0.9383d0)
     
      
      phidegg = 180.d0- phideg
      grtorad = 2.d0*pi/360.d0 
      phin = phidegg*grtorad   
      gev2nanobarn = (197.31d0**2*10.d0)
      s = dm*dm +2.d0*ebeam*(eh+dsqrt(eh**2-dm**2))
      y = dq2/(s-dm**2)/xb
      EPS = 2.d0*dm*xb/dsqrt(dq2)
      call REALCFFproton(xb,dq2,d2g,1,reh)
      call REALCFFproton(xb,dq2,d2g,2,ree)
    
      CALL INTERFERENCEMULLER(xb,d2g,dq2,ebeam,eh,phin,dlambda,reh,ree,
     1     DINTREST)
c      write(6,*) reh, ree
      CALL CBH(xb,d2g,dq2,ebeam,eh,phin,DCOREST)
      
      CALL DVCSMULLER(xb,d2g,dq2,ebeam,eh,PHIN,dlambda,REH,REE,
     1    DVCSREST)
    
      dkinrest = xb*y**2/(8.d0*pi*dq2**2*dsqrt(1.d0+eps**2))
      crosssectbhfermo = dkinrest*alphaelm**3*dcorest*gev2nanobarn
      crosssectintfermo = dkinrest*alphaelm**3*dintrest*gev2nanobarn
      crosssectdvcsfermo = dkinrest*alphaelm**3*gev2nanobarn*DVCSREST
      crosstotfree = crosssectbhfermo + crosssectintfermo +
     1     crosssectdvcsfermo
      RETURN
      END  

      SUBROUTINE CBH(xb,d2g,dq2,ebeam,ehad,phI,DCTOT)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (pi=dacos(-1.d0),dm=0.9383d0) 
      common/nn/in
      in = 1
      q = dsqrt(dq2)
      s = dm*dm +2.d0*ebeam*(ehad+dsqrt(ehad**2-dm**2))
      y = dq2/(s-dm**2)/xb
      Q2F = D2G*(1.D6/(197.31d0)**2)
      eps = 2.d0*dm*xb/q
    
      DK = DSQRT((4. *(1.D0  - xB)*xB + EPS**2)*
     -     (-1.D0  + y + (y**2*EPS**2)/4.D0)*(D2G +
     -     (2.D0  * Q**2*(1.D0  - xB) + Q**2*EPS**2 -
     -     2.D0  * Q**2*(1.D0  - xB) * DSQRT(1.D0  + EPS**2))/
     -     (4.D0  * (1.D0  - xB)*xB + EPS**2))*
     -     (D2G + (2.D0 *Q**2*(1.D0  - xB) +Q**2*EPS**2 +
     -     2.D0  *Q**2*(1.D0  - xB)*DSQRT(1.D0  +EPS**2))/
     -     (4.D0  * (1.D0  - xB)*xB + EPS**2)))/(2.D0 *Q**2)

      CALL FORMF(IN,5,-Q2F,F1,F2,GEP,gmP)
      
      A= (F1+ F2)**2
      B = F1**2 -F2**2*(D2G/(4.D0 * DM**2))

    
      DC0 =(  (-16.D0 *A*D2G**2*Q**2 - 
     -     2.D0 *A*D2G**2*(1.D0  - D2G/Q**2)*Q**2 + 
     -     2.D0 *B*D2G**2*(1.D0  - D2G/Q**2)*Q**2 + 
     -     2.D0 *A*D2G**2*(1.D0  + D2G/Q**2)*Q**2 + 
     -     2.D0 *B*D2G**2*(1.D0  + D2G/Q**2)*Q**2 - 
     -     6.D0 *A*D2G*(1.D0  - D2G/Q**2)*Q**4 + 
     -     2.D0 *B*D2G*(1.D0  - D2G/Q**2)*Q**4 - 
     -     2.D0 *A*D2G*(1.D0  + D2G/Q**2)*Q**4 + 
     -     2.D0 *B*D2G*(1.D0  + D2G/Q**2)*Q**4 + 
     -     6.D0 *A*D2G*(1.D0  - D2G/Q**2)*(1.D0  + D2G/Q**2)*Q**4 - 
     -     2.D0 *B*D2G*(1.D0  - D2G/Q**2)*(1.D0  + D2G/Q**2)*Q**4 + 
     -     8.D0 *A*D2G*(1.D0  + D2G/Q**2)**2*Q**4 - 
     -     4.D0 *B*D2G*(1.D0  + D2G/Q**2)**2*Q**4 + 
     -     2.D0 *A*(1.D0  + D2G/Q**2)**2*Q**6 - 
     -     2.D0 *B*(1.D0  + D2G/Q**2)**2*Q**6 - 
     -     2.D0 *A*(1.D0  + D2G/Q**2)**3*Q**6 + 
     -     2.D0 *B*(1.D0  + D2G/Q**2)**3*Q**6 +
     -     8.D0 *A*D2G*Q**2*DM**2 - 
     -     8.D0 *B*D2G*Q**2*DM**2 - 
     -     4.D0 *A*(1.D0  + D2G/Q**2)**2*Q**4*DM**2 + 
     -     4.D0 *B*(1.D0  + D2G/Q**2)**2*Q**4*DM**2 - 
     -     8.D0 *A*D2G*Q**2*(1.D0 - D2G/(2.D0*DM**2))*DM**2 - 
     -     8.D0 *B*D2G*Q**2*(1.D0 - D2G/(2.D0*DM**2))*DM**2 + 
     -     4.D0 *A*(1.D0  + D2G/Q**2)**2*Q**4*(1.D0 - D2G/(2.D0*DM**2))*
     -     DM**2 + 4.D0 *B*(1.D0  + D2G/Q**2)**2*Q**4*
     -     (1.D0  - D2G/(2.D0*DM**2))*DM**2 - (4.D0 *B*D2G*Q**4)/XB**2 - 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*Q**6)/XB**2 + 
     -     (4.D0 *B*(1.D0  + D2G/Q**2)*Q**6)/XB**2 + 
     -     (4.D0 *B*D2G*Q**4)/XB + 
     -     (2.D0 *A*D2G*(1.D0  - D2G/Q**2)*Q**4)/XB + 
     -     (2.D0 *B*D2G*(1.D0  - D2G/Q**2)*Q**4)/XB + 
     -     (2.D0 *A*D2G*(1.D0  + D2G/Q**2)*Q**4)/XB - 
     -     (2.D0 *B*D2G*(1.D0  + D2G/Q**2)*Q**4)/XB + 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*Q**6)/XB - 
     -     (4.D0 *B*(1.D0  + D2G/Q**2)*Q**6)/XB - 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*(1.D0  + D2G/Q**2)*Q**6)/XB - 
     -     (2.D0 *A*(1.D0  + D2G/Q**2)**2*Q**6)/XB + 
     -     (2.D0 *B*(1.D0  + D2G/Q**2)**2*Q**6)/XB - 
     -     (4.D0 *B*D2G*Q**4*(1.D0  + (D2G*XB)/Q**2))/XB**2 - 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*Q**6*(1.D0  + (D2G*XB)/Q**2))/
     -     XB**2 + (4.D0 *B*(1.D0  + D2G/Q**2)*Q**6*
     -     (1.D0  + (D2G*XB)/Q**2))/XB**2 - 
     -     (2.D0 *A*D2G*(1.D0  - D2G/Q**2)*Q**4*
     -     (1.D0  + (D2G*XB)/Q**2))/
     -     XB - (2.D0*B*D2G*(1.D0  - D2G/Q**2)*Q**4*
     -     (1.D0  + (D2G*XB)/Q**2))/XB - 
     -     (2.D0 *A*D2G*(1.D0  + D2G/Q**2)*Q**4*
     -     (1.D0  + (D2G*XB)/Q**2))/
     -     XB - (2.D0 *B*D2G*(1.D0  + D2G/Q**2)*Q**4*
     -     (1.D0  + (D2G*XB)/Q**2))/XB + 
     -     (2.D0 *A*(1.D0  + D2G/Q**2)**2*Q**6*(1.D0  + (D2G*XB)/Q**2))/
     -     XB + (2.D0 *B*(1.D0  + D2G/Q**2)**2*Q**6*
     -     (1.D0  + (D2G*XB)/Q**2))/XB + 
     -     (32.D0 *A*D2G*DK**2*Q**4)/((1.D0  + EPS**2)**2*y**2) - 
     -     (16.D0 *B*D2G*DK**2*Q**4)/((1.D0  + EPS**2)**2*y**2) + 
     -     (8.D0 *A*DK**2*Q**6)/((1.D0  + EPS**2)**2*y**2) - 
     -     (8.D0 *B*DK**2*Q**6)/((1.D0  + EPS**2)**2*y**2) + 
     -     (4.D0 *A*DK**2*(1.D0  - D2G/Q**2)*Q**6)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (4.D0 *B*DK**2*(1.D0  - D2G/Q**2)*Q**6)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (12.D0 *A*DK**2*(1.D0  + D2G/Q**2)*Q**6)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (12.D0 *B*DK**2*(1.D0  + D2G/Q**2)*Q**6)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (16.D0 *A*DK**2*Q**4*DM**2)/((1.D0  + EPS**2)**2*y**2) + 
     -     (16.D0 *B*DK**2*Q**4*DM**2)/((1.D0  + EPS**2)**2*y**2) + 
     -     (16.D0 *A*DK**2*Q**4*(1.D0 - D2G/(2.D0*DM**2))*DM**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (16.D0 *B*DK**2*Q**4*(1.D0 - D2G/(2.D0*DM**2))*DM**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (8.D0 *B*(1.D0  - D2G/Q**2)*Q**6)/(XB**2*y**2) + 
     -     (8.D0 *B*(1.D0  + D2G/Q**2)*Q**6)/(XB**2*y**2) - 
     -     (8.D0 *A*DK**2*Q**6)/((1.D0  + EPS**2)**2*XB*y**2) - 
     -     (8.D0 *B*DK**2*Q**6)/((1.D0  + EPS**2)**2*XB*y**2) + 
     -     (8.D0 *A*DK**2*Q**6*(1.D0  + (D2G*XB)/Q**2))/
     -     ((1.D0  + EPS**2)**2*XB*y**2) + 
     -     (8.D0 *B*DK**2*Q**6*(1.D0  + (D2G*XB)/Q**2))/
     -     ((1.D0  + EPS**2)**2*XB*y**2) + 
     -     (8.D0 *B*D2G*Q**4)/(XB**2*y) + 
     -     (12.D0 *B*(1.D0  - D2G/Q**2)*Q**6)/(XB**2*y) - 
     -     (8.D0 *B*(1.D0  + D2G/Q**2)*Q**6)/(XB**2*y) - 
     -     (8.D0 *B*D2G*Q**4)/(XB*y) + 
     -     (8.D0 *B*D2G*(1.D0  + D2G/Q**2)*Q**4)/(XB*y) - 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*Q**6)/(XB*y) + 
     -     (8.D0 *B*(1.D0  + D2G/Q**2)*Q**6)/(XB*y) + 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*(1.D0  + D2G/Q**2)*Q**6)/(XB*y) - 
     -     (8.D0 *B*(1.D0  + D2G/Q**2)**2*Q**6)/(XB*y) + 
     -     (8.D0 *B*D2G*Q**4*(1.D0  + (D2G*XB)/Q**2))/(XB**2*y) + 
     -     (4.D0 *B*(1.D0  - D2G/Q**2)*Q**6*(1.D0  + (D2G*XB)/Q**2))/
     -     (XB**2*y) - (8*B*(1.D0  + D2G/Q**2)*Q**6*
     -     (1.D0  + (D2G*XB)/Q**2))/(XB**2*y) + 
     -     (8.D0 *B*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB**2*y**2) + 
     -     (16.D0 *B*D2G*Q**4*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y**2) + 
     -     (8.D0 *B*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y**2) + 
     -     (8.D0 *B*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y**2) - 
     -     (16.D0 *B*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y**2) - 
     -     (8*B*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB**2*y**2) + 
     -     (8.D0 *A*D2G**2*Q**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (4.D0 *A*D2G*Q**4*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (4.D0 *B*D2G*Q**4*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (4.D0 *A*D2G*(1.D0  - D2G/Q**2)*Q**4*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (12.D0 *A*D2G*(1.D0  + D2G/Q**2)*Q**4*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (8.D0 *B*D2G*(1.D0  + D2G/Q**2)*Q**4*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (2.D0 *A*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (2.D0 *B*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (6.D0 *A*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (6.D0 *B*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (2.D0 *A*(1.D0  - D2G/Q**2)*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (2.D0 *B*(1.D0  - D2G/Q**2)*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (6.D0 *A*(1.D0  + D2G/Q**2)**2*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (6.D0 *B*(1.D0  + D2G/Q**2)**2*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (8.D0 *A*(1.D0  + D2G/Q**2)*Q**4*DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (8.D0 *B*(1.D0  + D2G/Q**2)*Q**4*DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (8.D0 *A*(1.D0  + D2G/Q**2)*Q**4*(1.D0 - D2G/(2.D0*DM**2))*
     -     DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) + 
     -     (8.D0 *B*(1.D0  + D2G/Q**2)*Q**4*(1.D0 - D2G/(2.D0*DM**2))* 
     -     DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*y) - 
     -     (4.D0 *B*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB**2*y) + 
     -     (4.D0 *A*D2G*Q**4*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (12.D0 *B*D2G*Q**4*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (4.D0 *B*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) + 
     -     (2.D0 *A*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (6.D0 *B*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (6.D0 *A*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) + 
     -     (6.D0 *B*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) + 
     -     (4.D0 *B*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB**2*y) - 
     -     (4.D0 *A*D2G*Q**4*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (4.D0 *B*D2G*Q**4*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (2.D0 *A*(1.D0  - D2G/Q**2)*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) - 
     -     (2.D0 *B*(1.D0  - D2G/Q**2)*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) + 
     -     (6.D0 *A*(1.D0  + D2G/Q**2)*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) + 
     -     (6.D0 *B*(1.D0  + D2G/Q**2)*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0)))/
     -     ((1.D0  + EPS**2)*XB*y) + 
     -     (16.D0 *A*D2G*Q**4*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (8.D0 *B*D2G*Q**4*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (4.D0 *A*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (4.D0 *B*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (2.D0 *A*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (2.D0 *B*(1.D0  - D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (6.D0 *A*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (6.D0 *B*(1.D0  + D2G/Q**2)*Q**6*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (8.D0 *A*Q**4*DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (8.D0 *B*Q**4*DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (8.D0 *A*Q**4*(1.D0 - D2G/(2.D0*DM**2))*DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) + 
     -     (8.D0 *B*Q**4*(1.D0 - D2G/(2.D0*DM**2)) *DM**2*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*y**2) - 
     -     (4.D0 *A*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*XB*y**2) - 
     -     (4.D0 *B*Q**6*(-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*XB*y**2) + 
     -     (4.D0 *A*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*XB*y**2) + 
     -     (4.D0 *B*Q**6*(1.D0  + (D2G*XB)/Q**2)*
     -     (-((D2G*(1.D0  - XB)*(2.D0  - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (EPS**2*y)/2.D0))**2)/
     -     ((1.D0  + EPS**2)**2*XB*y**2)) )/(d2g**2*Q**4)
      
      propa =   -(((-((D2G*(1.D0  - XB)*(2.D0 - y))/Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (y*EPS**2)/2.D0 )+
     -     2.D0 *DK*DCOS(PHI))*
     -     (1.D0  + D2G/Q**2 + (-((D2G*(1.D0  - XB)*(2.D0  - y))/
     -     Q**2) + 
     -     (1.D0  + D2G/Q**2)*(1.D0  - y - (y*EPS**2)/2.D0 ) +
     -     2.D0  *DK*DCOS(PHI))/
     -     (y*(1.D0  +EPS**2))))/
     -     (y*(1.D0  +EPS**2)))
      
      dc2 = ( (16.D0 *DK**2*(A*D2G +
     -     2.D0 *B*DM**2)*Q**4)/(y**2*(1.D0  +EPS**2)**2) )/(d2g**2*
     1     Q**4)
     
      
      dc1 =  ( (-16.D0 *DK*Q**2*(-(B*D2G*Q**2) -
     -     B*D2G*EPS**2*Q**2 -
     -     A*D2G**2*XB - 
     -     2.D0 *B*D2G*DM**2*XB + A*D2G*Q**2*XB + 2.D0 *B*DM**2*
     -     Q **2*XB + 2.D0 *A*D2G**2*XB**2 + 4.D0 *B*D2G*DM**2*XB**2)*
     -     (-2.D0  + y))/((1.D0  + EPS**2)**2*XB*y**2) )/(d2g**2*q**4)


      DCTOT = ( DC0 + DC1*DCOS(PHI) + DC2*DCOS(2.D0 *PHI) ) / PROPA
      RETURN
      END

      SUBROUTINE IMCFFF1rest(xb,d2g,dq2,ikind,DIMF)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI = dacos(-1.d0),dm =0.938d0)
      DIMENSION  GPDIM(2)
      COMMON/PAR/D,A,B,DKS
      common/skwe/xi
      common/nn/in
      
      xi = xb/(2.d0-xb)

     
c Im part for the GPD H
c      

      IF(IKIND.EQ.1)THEN

         DO ISEGN = 1,2         
            IF(ISEGN.EQ.1)THEN
               
               XX = XI
            ELSE
               XX = -XI
            ENDIF
            
            CALL GPDH(XX,XI,d2g,dq2,1,IN,GPDU)
            CALL GPDH(XX,XI,d2g,dq2,2,IN,GPDD)
            CALL GPDH(XX,XI,d2g,dq2,3,IN,GPDS)
            
            gpdsu= dks*gpds

           
            GPDT= ((4.D0/9.D0*(gpdu+gpdsu)+1./9.D0*
     *           (gpdd+gpdsu)+1./9.*gpds))

            GPDIM (ISEGN)  = GPDT
         ENDDO
         DIMF = GPDIM(1) - GPDIM(2)
c         write(6,*) 'imh', dimf, dks
      ELSE
ccc
ccc Im part for the GPDE         
ccc     
         DO ISEGN = 1,2         
            IF(ISEGN.EQ.1)THEN
               
               XX = XI
            ELSE
               XX = -XI
            ENDIF
            
            CALL GPDE(XX,XI,d2g,dq2,1,IN,GPDU)
            CALL GPDE(XX,XI,d2g,dq2,2,IN,GPDD)
            CALL GPDE(XX,XI,d2g,dq2,3,IN,GPDS)
            
            gpdsu= gpds
            GPDT= ((4.D0/9.D0*(gpdu+gpdsu)+1./9.D0*
     *              (gpdd+gpdsu)+1./9.*gpds))
          

            GPDIM (ISEGN)  = GPDT
         ENDDO
         DIMF = GPDIM(1) - GPDIM(2)
c         write(6,*) 'ime',dimf
      ENDIF
     
      RETURN
      END 
      
      
      subroutine GPDH(X,XI,d2g,dq2,I,in,GPD)
      IMPLICIT REAl*8(A-H,O-Z)
      DIMENSION CONT(3),coe(3,4),HIJ(10),CONTAt(3),HIJI(10),hj(4)
      COMMON/PAR/D,A,B,DKS
      COMMON/MATR/C(4)
    
      CALL COEFFP(I,dq2)
      CALL PARAREGGEP(I,dq2)

      IF(I.EQ.1.or.i.eq.2) THEN
         
c         if(xi.lt.0.0001d0)then

         if(xi.le.0.000001d0)then
c            IF(X.LE.0.D0)THEN
c               WRITE(6,*)x, 'ATTENZIONE X NEGATIVE IN GPD(XI=0) H ', xi
c               STOP
c            ENDIF
            
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
                        DGRA = (XI**2-X)*(X1**DMIJ - X2**DMIJ)+DMIJ*XI*
     *                       (1.D0 -X)*(X1**DMIJ+x2**DMIJ)

                        HIJ(J)=(3.D0*GAMMA(DMIJ-1.D0)*DGRA*C(J))/
     *                      (2.D0*XI**3.D0* GAMMA(DMIJ+2.D0))
                     ELSE
CCCCCCCHVALJ ERBL

                        DGRA= XI**2-X + DMIJ*XI*(1.D0-X)
                        HIJ(J)=(3.D0*GAMMA(DMIJ-1.D0)*DGRA*C(J)*
     *                      X1**DMIJ)/(2.D0*XI**3.D0* GAMMA(DMIJ+2.D0))
                     ENDIF
                     
                     
                     cont(ii)=cont(ii)+hij(j)
                  enddo
                  
 12            ENDDO
               
               GPD=(CONT(1)+CONT(2))*DEXP(B*D2G)
            endif
         endif
      ELSE
         
         xt=dabs(x)
         if(xi.le.0.000001d0)then
            
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
                  DGRA =((DMIJ**2+2.D0)*(XI**2-XT)**2-(DMIJ**2-1.D0)*
     *                 (1.D0-XI**2)*
     -                 (XT**2-XI**2))* (X1**DMIJ -X2**DMIJ)+3.D0*DMIJ*
     *                 XI*(1.D0-XT)*(XI**2-XT)*(X1**DMIJ+X2**DMIJ)
                  
                  HIJI(J) = (15.D0*GAMMA(DMIJ-2.D0)*DGRA*C(J))/
     *                 (2.D0*XI**5*GAMMA(DMIJ+3.D0))
                  
               else
                  
CCCCCCCHSEAJ ERBL
                  
                  DGRA=(X1**DMIJ)*((DMIJ**2 +2.D0)*(XI**2-XT)**2.D0 +
     *                  3.D0*DMIJ*XI*(1.D0-XT)*(XI**2-XT)-
     *                 (DMIJ**2-1.D0)*(1.D0-XI**2)*(XT**2-XI**2))-
     -                 (X4**DMIJ)*((DMIJ**2+2.D0)*(XI**2+XT)**2.D0 +
     +                 3.D0*DMIJ*XI*(1.D0+XT)*(XI**2+XT)-(DMIJ**2-1.D0)*
     *                 (1.D0-XI**2)*(XT**2-XI**2))
                  
                  HIJI(J) = (15.D0*GAMMA(DMIJ-2.D0)*DGRA*C(J))/
     *                 (2.D0*XI**5.D0*GAMMA(DMIJ+3.D0))
                  
                  
               ENDIF
               
               
               CONTAT(III)= CONTAT(III)+HIJI(J)
               
            ENDDO

            
 66      ENDDO
         endif
         GPD = (CONTAT(1)+CONTAT(2))*DEXP(B*D2G)
         if (x.lt.0) GPD= -(CONTAT(1)+CONTAT(2))*DEXP(B*D2G)
         


      endif
      RETURN
      END

      SUBROUTINE COEFFP(I,dq2)
      IMPLICIT REAL*8(A-H,O-Z)
      common/Matr/c(4)
      
      DL = DLog(dq2/4.D0)
      IF(I.EQ.1)GO TO 30
      IF(I.EQ.2)GO TO 50
      IF(I.EQ.3)GO TO 20
 30   C(1) = 1.52D0 + 0.248d0*DL
      C(2) = 2.88D0 - 0.94d0*DL
      C(3) = -0.095D0*DL
      C(4) = 0.D0
      RETURN
 50   C(1) = 0.76D0 + 0.248*DL
      C(2) = 3.11D0 - 1.36*DL
      C(3) = -3.99D0 + 1.15*DL
      C(4) = 0.D0
      RETURN
 20   C(1) = 0.123D0 + 0.0003D0*DL
      C(2) = -0.327D0 - 0.004D0*DL
      C(3) = 0.692D0 - 0.068D0*DL
      C(4) = -0.486D0 + 0.038D0*DL
      END
      
      SUBROUTINE PARAREGGEP(I,dq2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PAR/D,A,B,DKS
      parameter (dm = 0.9383d0)
      
      DL = DLOG(dq2/4.D0)
      DKS = 1.d0 + 0.68d0/(1.d0 + 0.52d0*DL)  
      IF(I.EQ.1.OR.I.EQ.2)GO TO 1
      IF(I.EQ.3)GO TO 33
 1    D = 0.48D0
      A = 0.9D0
      B = 0.D0
      RETURN   
 33   D =  0.10D0 + 0.06*DL -0.027*DL*DL
      A = 0.15D0
      B = 2.58D0+0.25D0*DLOG(DM**2/(dq2+DM**2))
      RETURN
       END

      SUBROUTINE FORMF(NN,IFORM,DQM2,F1,F2,GE,GM)
CCCCCCC
CCCCCCCEVALUATE ELASTIC NUCLEON FORM FACTORS
CCCCCCCDIRAC F1 AND F2, SACHS GE AND GM
CCCCCCC
CCCCCCCDQM2 = VALUE OF THE FOUR-MOMENTUM TRANSFER (FM-2)
CCCCCCC
CCCCCCCNN  = 2 NEUTRON
CCCCCCCNN= 1  PROTON
CCCCCCC
CCCCCCCIFORM = 1 BLATNIK-ZOFKO
CCCCCCCIFORM = 2 GALSTER
CCCCCCCIFORM = 3 HOHLER
CCCCCCCIFORM = 4 GARI      ***********
CCCCCCCIFORM = 5 DIPOLE
CCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA DMOMG/0.06340592198D0/,DMPHI/0.03746990965D0/,
     1     DMOMGP/0.02780802576D0/,DMRHO/0.06654912153D0/,
     2     DMRHOP/0.02994710470D0/,DMRHOPP/0.01853868385D0/,
     3     AS/0.01837706810D0/,AV/0.0009255302930D0/,
     4     DMUS/0.440D0/,DMUV/2.353D0/,BS/0.03542742485D0/,
     5     BV/0.04282435971D0/,DM/4.75523795D0/
      DATA AMRO/0.776D0/,AMOM/0.784D0/,AKV/3.706D0/,AKS/-0.12D0/,
     1     GFRO/0.377D0/,GFOM/0.411D0/,AKRO/6.62D0/,AKOM/0.163D0/,
     2     ALA1/0.795D0/,ALA2/2.27D0/,ALQCD/0.29D0/,FI/1.02D0/
      TAU=0.25D0*DQM2/(DM*DM)
      GO TO (1010,1020,1030,1040,1050), IFORM
 1010 RS=1.D0/((1.D0+DMOMG*DQM2)*(1.D0+DMPHI*DQM2)*(1.D0+DMOMGP*DQM2))
      RV=1.D0/((1.D0+DMRHO*DQM2)*(1.D0+DMRHOP*DQM2)*(1.D0+DMRHOPP*DQM2))
      IF(NN.EQ.2) THEN
         GE=(0.5D0+AS*DQM2)*RS-(0.5D0+AV*DQM2)*RV
         GM=(DMUS+0.5D0*BS*DQM2)*RS-(DMUV+0.5D0*BV*DQM2)*RV
      ELSE
         GE=(0.5D0+AS*DQM2)*RS+(0.5D0+AV*DQM2)*RV
         GM=(DMUS+0.5D0*BS*DQM2)*RS+(DMUV+0.5D0*BV*DQM2)*RV
      ENDIF
      GO TO 1500
 1020 GE=1.D0/(1.D0+0.0560D0*DQM2)**2
      IF(NN.EQ.2) THEN
         GM=-1.913D0*GE
         GE=1.913D0*TAU*GE/(1.D0+5.6D0*TAU)
      ELSE
         GM=2.793D0*GE
      ENDIF
      GO TO 1500
 1030 Q2G=0.19731D0*0.19731D0*DQM2
      F1RO=0.5D0*(0.955D0+0.090D0/((1.D0+Q2G/0.355D0)**2))/
     /     (1.D0+Q2G/0.536D0)
      F2RO=0.5D0*(5.335D0+0.962D0/(1.D0+Q2G/0.268D0))/
     /     (1.D0+Q2G/0.603D0)
      TFO=AMOM*AMOM+Q2G
      F1O=0.71D0/TFO
      F2O=-0.11D0/TFO
      TFFI=1.02D0*1.02D0+Q2G
      F1FI=-0.64D0/TFFI
      F2FI=0.13D0/TFFI
      TFOI=1.8D0*1.8D0+Q2G
      F1OI=-0.13D0/TFOI
      F2OI=-0.02D0/TFOI
      TFR1=1.21D0*1.21D0+Q2G
      F1R1=0.05D0/TFR1
      F2R1=-1.99D0/TFR1
      TFR2=2.45D0*2.45D0+Q2G
      F1R2=-0.52D0/TFR2
      F2R2=0.20D0/TFR2
      TFR3=2.95D0*2.95D0+Q2G
      F1R3=0.28D0/TFR3
      F2R3=0.19D0/TFR3
      F1S=F1O+F1FI+F1OI
      F1V=F1RO+F1R1+F1R2+F1R3
      F2S=F2O+F2FI+F2OI
      F2V=F2RO+F2R1+F2R2+F2R3
      IF(NN.EQ.2) THEN
         F1=F1S-F1V
         F2=F2S-F2V
      ELSE
         F1=F1S+F1V
         F2=F2S+F2V
      ENDIF
      GE=F1-TAU*F2
      GM=F1+F2
      RETURN
 1040 Q2G=0.19731D0*0.19731D0*DQM2
      AMR2=AMRO*AMRO
      AMO2=AMOM*AMOM
      A1Q=ALA1*ALA1
      A2Q=ALA2*ALA2
      QCDQ=ALQCD*ALQCD
      Q2D=Q2G*DLOG((A2Q+Q2G)/QCDQ)/DLOG(A2Q/QCDQ)
      FL1=A1Q/(A1Q+Q2D)
      FL2=A2Q/(A2Q+Q2D)
      F1=FL1*FL2
      F2=FL1*FL2*FL2
      F1V=0.5D0*F1*(AMR2*GFRO/(AMR2+Q2G)+1.D0-GFRO)
      F2V=0.5D0*F2*(AMR2*GFRO*AKRO/(AMR2+Q2G)+AKV-GFRO*AKRO)
      F1S=0.5D0*F1*(AMO2*GFOM/(AMO2+Q2G)+1.D0-GFOM)
      F2S=0.5D0*F2*(AMO2*GFOM*AKOM/(AMO2+Q2G)+AKS-GFOM*AKOM)
      IF(NN.EQ.2) THEN
         F1=F1S-F1V
         F2=F2S-F2V
      ELSE
         F1=F1S+F1V
         F2=F2S+F2V
      ENDIF
      GE=F1-TAU*F2
      GM=F1+F2
      RETURN
 1050 GE=1.D0/(1.D0+0.0560D0*DQM2)**2
      IF(NN.EQ.2) THEN
         GM=-1.913D0*GE
         GE=0.D0
      ELSE
         GM=2.793D0*GE
      ENDIF
 1500 F1=(GE+TAU*GM)/(1.D0+TAU)
      F2=(GM-GE)/(1.D0+TAU)
      RETURN
      END


      subroutine REALCFFproton(xb,dq2,d2g,ikind,risretot)
      IMPLICIT REAl*8(A-H,O-Z)

    
      dimension xp(150),ww(150)
      external h,e
      dimension ax(2), nx(1)
      COMMON/PAR/D,A,B,DKS
c      common/skwe/xi
      common/kin/xbb,dq22,d2gg
      common/nn/in

      xbb = xb
      dq22 = dq2
      d2gg = d2g

c      write(*,*)'xi',xi, xb
      xi = xb/(2.d0-xb)
      aa = 0.d0
      bb = 1.d0
      dl = xi
      in = 1
      eps = 0.01d0
     
c
      if(ikind.eq.1)then
         ris = dcauch(h,aa,bb,dl,eps)
c         write(6,*)'h',ris, xi
      else
         ris = dcauch(e,aa,bb,dl,eps)
      endif
  
      mx = 1
      ax(1) = aa
      ax(2) = bb

      nx(1)=44
      
      
      CALL DMPPP(mx,NPOINTxpr,Nx,Ax,xp,ww)

      dintx=0.d0
      
ccccccc
cccccccCALCULATE THE PART OF THE INTEGRAL OVER X WITH KERNEL 1/(X+XI)
ccccccc
      
      do ipr = 1,npointxpr

         xpr = xp(ipr)
         dker = ( 1.d0/(xpr+xi))
                 
         if(ikind.eq.1)then
            CALL GPDH(xpr,XI,d2g,dq2,1,in,GPDU)
            CALL GPDH(xpr,XI,d2g,dq2,2,IN,GPDD)
            CALL GPDH(xpr,XI,d2g,dq2,3,IN,GPDS)
            CALL GPDH(-xpr,XI,d2g,dq2,1,in,GPDU1)
            CALL GPDH(-xpr,XI,d2g,dq2,2,IN,GPDD1)
            CALL GPDH(-xpr,XI,d2g,dq2,3,IN,GPDS1)
            gpdsu =DKS*gpds
            gpdsu1 = dks*gpds1

            
         else 
            CALL GPDE(xpr,XI,d2g,dq2,1,in,GPDU)
            CALL GPDE(xpr,XI,d2g,dq2,2,IN,GPDD)
            CALL GPDE(xpr,XI,d2g,dq2,3,IN,GPDS)
            CALL GPDE(-xpr,XI,d2g,dq2,1,in,GPDU1)
            CALL GPDE(-xpr,XI,d2g,dq2,2,IN,GPDD1)
            CALL GPDE(-xpr,XI,d2g,dq2,3,IN,GPDS1)
            gpdsu =gpds
            gpdsu1 = gpds1
            
         endif
         gpdut = gpdu-gpdu1
         gpddt = gpdd-gpdd1
         gpdst = gpds-gpds1
         gpdsut = gpdsu-gpdsu1

         
         GPDT=((4.D0/9.d0*(gpdut+gpdsut)+1.d0/9.d0*
     *        (gpddt+gpdsut)+1.d0/9.d0*gpdst))
         
      
      dintx = dintx +gpdt*dker*ww(ipr)
c      write(6,*)'gpdt', gpdt, dker, ww(ipr)
      enddo

      risretot =  dintx + ris 

      RETURN
      END

      double precision function h(xpr)

      implicit real*8(a-h,o-z)
      COMMON/PAR/D,A,B,DKS
      common/skwe/xi
      common/kin/xbb,dq22,d2gg
      common/nn/in

    
      xb = xbb
      dq2 = dq22
      d2g =  d2gg
      xi = xb/(2.d0-xb)
      
      CALL GPDH(xpr,XI,d2g,dq2,1,in,GPDU)
      CALL GPDH(xpr,XI,d2g,dq2,2,IN,GPDD)
      CALL GPDH(xpr,XI,d2g,dq2,3,IN,GPDS)

      CALL GPDH(-xpr,XI,d2g,dq2,1,in,GPDU1)
      CALL GPDH(-xpr,XI,d2g,dq2,2,IN,GPDD1)
      CALL GPDH(-xpr,XI,d2g,dq2,3,IN,GPDS1)
      gpdsu= dks*gpds

      gpdsu1 = dks*gpds1

      gpdut = gpdu-gpdu1
      gpddt = gpdd-gpdd1
      gpdst = gpds-gpds1
      gpdsut = gpdsu-gpdsu1

      GPDT=((4.D0/9.d0*(gpdut+gpdsut)+1.d0/9.d0*
     *     (gpddt+gpdsut)+1.d0/9.d0*gpdst))
      
      h = gpdt/(xpr-xi)

      return
      end

      double precision function e(xpr)

      implicit real*8(a-h,o-z)
      common/skwe/xi
      common/kin/xbb,dq22,d2gg
      common/nn/in

     
      xb = xbb
      dq2 = dq22
      d2g =  d2gg

            xi = xb/(2.d0-xb)

      CALL GPDE(xpr,XI,d2g,dq2,1,in,GPDU)
      CALL GPDE(xpr,XI,d2g,dq2,2,IN,GPDD)
      CALL GPDE(xpr,XI,d2g,dq2,3,IN,GPDS)
      CALL GPDE(-xpr,XI,d2g,dq2,1,in,GPDU1)
      CALL GPDE(-xpr,XI,d2g,dq2,2,IN,GPDD1)
      CALL GPDE(-xpr,XI,d2g,dq2,3,IN,GPDS1)
      
      gpdsu =gpds
      gpdsu1 = gpds1

      gpdut = gpdu-gpdu1
      gpddt = gpdd-gpdd1
      gpdst = gpds-gpds1
      gpdsut = gpdsu-gpdsu1

      GPDT=((4.D0/9.d0*(gpdut+gpdsut)+1.d0/9.d0*
     *     (gpddt+gpdsut)+1.d0/9.d0*gpdst))
      
      
      e = gpdt/(xpr-xi)

      return
      end


      SUBROUTINE DMPPP(MX,NPX,NX,AX,X,WX)                                 DGA00010
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
      CALL DVKNP(NPX)                                                    DGA00240
      CALL GPTP(POINT,NPX)                                               DGA00250
      CALL GWTP(WEIGHT,NPX)                                              DGA00260
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
      SUBROUTINE DVKNP(N)                                                DGA00450
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
      SUBROUTINE GPTP(P,N)                                               DGA01250
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
      SUBROUTINE GWTP(W,N)                                               DGA04340
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
 
      subroutine gpde(X,XI,D2G,dq2,I,nn,GPD)
       IMPLICIT REAl*8(A-H,O-Z)
       DIMENSION CONT(3),coe(3,4),HIJ(10),CONTAt(3),HIJI(10),hj(4)
       COMMON/PARE/DE,AE,BE,BETA,DKAP
       COMMON/MATRE/CE(8)
 
       CALL COEFFE(I,dq2)
       CALL PARAREGGEE(I,nn,dq2)
CCCCCCC II=1 DGLAP
CCCCCCC II=2 ERBL

       IF(I.EQ.1.or.i.eq.2) THEN


ccccc valence region (forward case)          
c             
          if(xi.le.0.000001d0)then
             
c             IF(X.LE.0.D0)THEN
c                WRITE(6,*)'ATTENZIONE X NEGATIVE IN GPD(XI=0)',x,xi
c                STOP
c             ENDIF
             
             GPD = DGAMMA(2.d0-dE+BETA)/(DGAMMA(1.d0-de)*DGAMMA(1.d0+
     +            BETA))*DKAP*x**(-De)*(1.d0-x)**BETA
             
cccc  valence region 
             
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
                   
                   
                   DO J=1,8
                      
                      DMIJ= 2.D0 + DFLOAT(j-1) - DE-(AE*D2G)
                      IF (II.EQ.1) THEN 
                         DGRA = (XI**2-X)*(X1**DMIJ - X2**DMIJ)+DMIJ*XI*
     *                        (1.D0 -X)*(X1**DMIJ+x2**DMIJ)

                         HIJ(J)=(3.D0*DGAMMA(DMIJ-1.D0)*DGRA*CE(J))/
     *                        (2.D0*XI**3.D0*DGAMMA(DMIJ+2.D0))
                      ELSE
CCCCCCC HVALJ ERBL

                         DGRA= XI**2-X + DMIJ*XI*(1.D0-X)
                         HIJ(J)=(3.D0*DGAMMA(DMIJ-1.D0)*DGRA*CE(J)*
     *                       X1**DMIJ)/(2.D0*XI**3.D0*DGAMMA(DMIJ+2.D0))
                      ENDIF
                      
                      
                      cont(ii)=cont(ii)+hij(j)
                   enddo
                   
 12             ENDDO
                
                GPD=(CONT(1)+CONT(2))*DEXP(BE*D2G)
             endif
          endif
       ELSE
          
          xt=dabs(x)
          if(xi.le.0.000001d0)then

c     N_s= -0.155 ( see Tab. 1 di https://arxiv.org/abs/0809.4126v1; checked with PARTONS)

          delta =-(1.1d0+0.06d0*dlog(dq2/4.d0)-0.0027*dlog(dq2/4.d0)**2)
          GPD = -0.155d0*xt**delta*(1.d0-x)**7
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
                   

                   DMIJ= 3.d0 + dfloat(j-1)-DE-(AE*D2G)
                   
                   IF(III.EQ.1) THEN        
                      DGRA =((DMIJ**2+2.D0)*(XI**2-XT)**2-(DMIJ**2-
     -                     1.D0)*(1.D0-XI**2)*(XT**2-XI**2))* 
     *                     (X1**DMIJ -X2**DMIJ)+3.D0*DMIJ*XI*(1.D0-XT)*
     *                     (XI**2-XT)*(X1**DMIJ+X2**DMIJ)
                      
                      HIJI(J) = (15.D0*DGAMMA(DMIJ-2.D0)*DGRA*CE(J))/
     /                      (2.D0*XI**5*DGAMMA(DMIJ+3.D0))
                      
                   else
                      
CCCCCCCHSEAJ ERBL
                      
                      DGRA=(X1**DMIJ)*((DMIJ**2 +2.D0)*(XI**2-XT)**2.D0 
     +                      + 3.D0*DMIJ*XI*(1.D0-XT)*(XI**2-XT)-
     -                     (DMIJ**2-1.D0)*(1.D0-XI**2)*(XT**2-XI**2))-
     -                     (X4**DMIJ)*((DMIJ**2+2.D0)*(XI**2+XT)**2.D0 +
     +                     3.D0*DMIJ*XI*(1.D0+XT)*(XI**2+XT)-
     -                     (DMIJ**2-1.D0)*(1.D0-XI**2)*(XT**2-XI**2))
                      
                      HIJI(J) = (15.D0*DGAMMA(DMIJ-2.D0)*DGRA*CE(J))/
     /                     (2.D0*XI**5.D0*DGAMMA(DMIJ+3.D0))
                      
                      
                   ENDIF
                   
                   
                   CONTAT(III)= CONTAT(III)+HIJI(J)
                   
                ENDDO

                
 66          ENDDO
             
             GPD = (CONTAT(1)+CONTAT(2))*DEXP(BE*D2G)
          endif
          if (x.lt.0) GPD= -(CONTAT(1)+CONTAT(2))*DEXP(BE*D2G)
          
       endif
       RETURN
       END


      SUBROUTINE COEFFE(I,dq2)
      IMPLICIT REAL*8(A-H,O-Z)
       common/Matre/cE(8)
      
c       DQ2 =4 .D0
       DL = DLOG(DQ2/4.D0)
       IF(I.EQ.1)GO TO 30
       IF(I.EQ.2)GO TO 50
       IF(I.EQ.3)GO TO 20
30      CE(1) = 2.2053d0
        CE(2) = -cE(1)
        CE(3) = 0.d0
        CE(4) = 0.D0
        CE(5) = 0.D0
        CE(6) = 0.D0
        CE(7) = 0.D0
        CE(8) = 0.D0
        
        RETURN
50       CE(1) = -3.114D0
         CE(2) = 8.096D0
         CE(3) = -6.477D0
         CE(4) = 1.295D0
         CE(5) = 0.1296D0
         CE(6) = 0.0362D0
         CE(7) = 0.014516D0
         CE(8) = 0.0070504D0
      RETURN
 20     CE(1) =-0.155d0
        CE(2) = -2.d0*cE(1)
        CE(3) = cE(1)
        CE(4) = 0.d0  
        CE(5) = 0.D0
        CE(6) = 0.D0
        CE(7) = 0.D0
        CE(8) = 0.D0
        
      return   
      END
     
      SUBROUTINE PARAREGGEE(I,nn,dq2)
      IMPLICIT REAL*8(A-H,O-Z)
   
       COMMON/PARE/DE,AE,BE,BETA,DKAP
    
       if(nn.eq.1)then
          DM = 0.9383d0
       else
          DM = 0.9396D0
       endif
      
       DL = DLOG(DQ2/4.D0)
      
       IF(I.EQ.1)GO TO 1
       IF(I.EQ.2)GO TO 13
       IF(I.EQ.3)GO TO 33
 1       DE = 0.48D0
         AE = 0.9D0
         BE = 0.D0
         BETA = 4.D0
         DKAP = 1.67D0
         RETURN
 13      DE = 0.48D0
         AE = 0.9D0
         BE = 0.D0
         BETA = 5.6D0
         DKAP = -2.03D0
         RETURN   
 33      DE = 1.1D0 + 0.06*DL -0.0027d0*DL*DL 
         AE = 0.15D0
         BE = 2.58D0 + 0.25D0*DLOG(dm**2/(DQ2+dm**2))
         RETURN
      
         END

      SUBROUTINE INTERFERENCEMULLER(xb,d2g,dq2,ebeam,ehad,PHI,dlambda,
     1    reh,ree,DINTREST)
      IMPLICIT REAl*8(A-H,O-Z)
      PARAMETER (pi=dacos(-1.d0),dm=0.9383d0)
      common/skwe/xi
      common/nn/in
     
      
      q = dsqrt(dq2)
      s = dm*dm +2.d0*ebeam*(ehad+dsqrt(ehad**2-dm**2))
      y = dq2/(s-dm**2)/xb
      EPS = 2.D0*DM*XB/Q
      Q2F = D2G*(1.D6/(197.31D0)**2)

      CALL IMCFFF1REST(xb,d2g,dq2,1,DIMF)
      CALL IMCFFF1REST(xb,d2g,dq2,2,DIMFe)
      CALL FORMF(IN,5,-Q2F,F1,F2,GEP,gmP)
c      call REALCFFproton(xb,dq2,d2g,1,reh)
c      call REALCFFproton(xb,dq2,d2g,2,ree)
c     write(6,*) 'reh', reh, 'ree',ree
c      ree =0.d0
c      reh=0.d0

      reheff =  (-2.d0*xi)/(1.d0+xi)*reh
      reeeff = (-2.d0*xi)/(1.d0+xi)*ree
      
      drec0 = f1*reh-d2g/(4.d0*dm*dm)*f2*ree
      drec0eff = f1*reheff-d2g/(4.d0*dm*dm)*f2*reeeff
      drec0delta = -xi**2*(f1+f2)*(reh+ree)
      dims2 = f1*(-2.d0*xi/(1.d0+xi))*dimf-
     -     d2g/(4.d0*dm*dm)*f2*dimfe*(-2.d0*xi/(1.d0+xi))
     
  
   
      DK = DSQRT((4.D0*(1.D0 - XB)*XB + EPS**2)*
     -     (-1.D0 + y + (y**2*EPS**2)/4.D0)*(D2G +
     -     (2.D0 * Q**2*(1.D0 - XB) + Q**2*EPS**2 -
     -     2.D0 * Q**2*(1.D0 - XB) * DSQRT(1.D0 + EPS**2))/
     -     (4.D0 * (1.D0 - XB)*XB + EPS**2))*
     -     (D2G + (2.D0*Q**2*(1.D0 - XB) +Q**2*EPS**2 +
     -     2.D0 *Q**2*(1.D0 - XB)*DSQRT(1.D0 +EPS**2))/
     -     (4.D0 * (1.D0 - XB)*XB + EPS**2)))/(2.D0*Q**2)


      DJEPS = (1.D0- Y -Y * EPS**2/2.D0)*(1.D0 + D2G/Q**2) - 
     1        (1.D0 - XB)*(2.D0 - Y)*D2G/Q**2
 
      pr1 = -1.d0/(y*(1.D0+EPS**2))*(DJEPS + 2.D0 * DK * DCOS(PHI))
      pr2 =  1.d0 + D2G/Q**2 + 1.d0/(y*(1.D0+EPS**2))*(DJEPS + 2.D0 *DK*
     1     DCOS(PHI))

      PROPA = PR1 * PR2

      dcom = 1.d0/(y**3*xb*d2g*propa)
  
      
CCCCCCC
CCCCCCCCfrom Muller et al. al.Nucl.Phys. B629 (2002) 323-392, DOI: 10.1016/S0550-3213(02)00144-X
CCCCCCC
    
 
      dc0 = ( -8.d0*(2.d0-y)*( (2.d0-y)**2/(1.d0-y)*dk*dk*drec0+
     +     d2g/dq2*(1.d0-y)*(2.d0-x)*(drec0 + drec0delta)) )*dcom
      
      dc1 = ( -2.d0*(2.d0-2.d0*y+y**2)*8.d0*dk*drec0 )*dcom*dcos(phi)
      
      dc2 = ( 16.d0*dk*dk/(2.d0-xb)*(y-2.d0)*drec0eff )*dcom*
     *     dcos(2.d0*phi)

      ds1 = 8.d0* dK *(2.D0- y)*y*(F1*dimf-F2*d2g/(4.d0*dm**2)*dimfe)*
     *     DSIN(PHI)*DLAMBDA*dcom
      
      ds2 = 16.d0*dk*dk/(2.d0-xb)*dlambda*y*dims2*dcom*dsin(2.d0*phi)

      dintrest =dc0+dc1+dc2+ds1+ds2

      RETURN
      END
    

        SUBROUTINE DVCSMULLER(xb,d2g,dq2,ebeam,ehad,PHI,dlambda,reh,ree,
     1    DVCSREST)
      IMPLICIT REAl*8(A-H,O-Z)
      PARAMETER (pi=dacos(-1.d0),dm=0.9383d0)
      common/nn/in
      common/skwe/xi
     
     
      q = dsqrt(dq2)
      s = dm*dm +2.d0*ebeam*(ehad+dsqrt(ehad**2-dm**2))
      y = dq2/(s-dm**2)/xb
      EPS = 2.D0*DM*XB/Q

c      call REALCFFproton(xb,dq2,d2g,1,reh)
      call IMCFFF1rest(xb,d2g,dq2,1,DIMFh)
c      call REALCFFproton(xb,dq2,d2g,2,ree)
c      ree = 0.d0
c      reh=0.d0
      call IMCFFF1rest(xb,d2g,dq2,2,DIMfe)
 
      dhh = dimfh**2+reh**2
      dhe = 2.d0*( ree*reh + DIMFe*DIMFh )
      dee = dimfe**2+ree**2
      
      cunp = 1.d0/(2.d0-xb)**2*(4.d0*(1.d0-xb)*dhh - xb**2*dhe -
     -     (xb**2+(2.d0-xb)**2*d2g/(4.d0*dm**2)*dee) )
      
      recunp =  -(2.d0*xi)/(1.d0+xi)*cunp
      
      imcunp =  0.d0

      dcom = 1.d0/(y**2*dq2)
      dc0dvcs = 2.d0*(2.d0-2.d0*y+y**2)*dcom*cunp
      dc1dvcs = (8.d0*dk*(2.d0-y)/(2.d0-xb)*recunp*dcos(phi))*dcom
      ds1dvcs = (8.d0*dk*(-dlambda*y)/(2.d0-xb)*imcunp*dsin(phi))*dcom

      dvcsrest = dc0dvcs+dc1dvcs+ds1dvcs

      RETURN
      END

