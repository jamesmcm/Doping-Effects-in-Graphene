c$$$Calculate xi for both lambda values (may be complex)
c$$$Calculate coefficients for lambda values and xi values - note lambda must match xi
c$$$Sub coefficients in to T formula
c$$$Take modulus and square
c$$$ REMEMBER THAT THE REAL N OF BLOCKS IS TWICE WHAT IS USED

      PROGRAM PLOTNOWRAP
      IMPLICIT NONE
      DOUBLE COMPLEX IM
      DOUBLE PRECISION E/-5.0/, LAMBDAP, LAMBDAN, SQRTARG
      INTEGER F/1/, N/1/
      CHARACTER*3 VALUE 
      DOUBLE COMPLEX XIP, XIN, AC1, AC2, BC1, BC2, CC1, CC2, DC1, DC2
      DOUBLE COMPLEX A, B, C, D, TP, TN

c$$$      set N to command line argument
      CALL GETARG(1, VALUE)
      READ(UNIT=VALUE, FMT=*) N
      
      IM=CMPLX(0,1)

      DO F = 1, 2002    
         LAMBDAP=E+1
         SQRTARG=((((E**2)*(LAMBDAP**2)) - 4*E*LAMBDAP))
         IF (SQRTARG .GE. 0) THEN
            XIP=DCMPLX(0.5*((E*LAMBDAP) -SQRT(SQRTARG) -2), 0)            
         ELSE
            XIP=DCMPLX(0.5*((E*LAMBDAP) -2), -0.5*SQRT(ABS(SQRTARG)))
         ENDIF
         AC2=(-1.0*(1+XIP))/((1.0/XIP) -XIP)
         AC1=1-AC2
         BC1=LAMBDAP/(XIP - (1.0/XIP))
         BC2=-1*BC1
         CC2=E/(XIP-(1.0/XIP))
         CC1=-1*CC2
         DC1=((E*LAMBDAP) - 1.0 - (1.0/XIP))/(XIP - (1.0/XIP))
         DC2 = 1-DC1
         A=(AC1*(XIP**N))+(AC2*(XIP**(-N)))
         B=(BC1*(XIP**N))+(BC2*(XIP**(-N)))
         C=(CC1*(XIP**N))+(CC2*(XIP**(-N)))
         D=(DC1*(XIP**N))+(DC2*(XIP**(-N)))
         TP=2.0/((A+D)+(IM*(C-B)))

         LAMBDAN=E-1
         SQRTARG=((((E**2)*(LAMBDAN**2)) - 4*E*LAMBDAN))
         IF (SQRTARG .GE. 0) THEN
            XIN=DCMPLX(0.5*((E*LAMBDAN) -SQRT(SQRTARG) -2), 0)            
         ELSE
            XIN=DCMPLX(0.5*((E*LAMBDAN) -2), -0.5*SQRT(ABS(SQRTARG)))
         ENDIF
c$$$         XIN=CMPLX(0.5*((E*LAMBDAN) -2), 0.5*SQRT(SQRTARG))
c$$$         XIN=CMPLX(0.5*((E*LAMBDAN) -2), 0.5*SQRT(COMPLEX
c$$$     +        ((((E**2)*(LAMBDAN**2)) - 4*E*LAMBDAN))))
         AC2=(-1.0*(1+XIN))/((1.0/XIN) -XIN)
         AC1=1-AC2
         BC1=LAMBDAN/(XIN - (1.0/XIN))
         BC2=-1*BC1
         CC2=E/(XIN-(1.0/XIN))
         CC1=-1*CC2
         DC1=((E*LAMBDAN) - 1.0 - (1.0/XIN))/(XIN - (1.0/XIN))
         DC2 = 1-DC1
         A=(AC1*(XIN**N))+(AC2*(XIN**(-N)))
         B=(BC1*(XIN**N))+(BC2*(XIN**(-N)))
         C=(CC1*(XIN**N))+(CC2*(XIN**(-N)))
         D=(DC1*(XIN**N))+(DC2*(XIN**(-N)))
         TN=2.0/((A+D)+(IM*(C-B)))

         WRITE (*,10) E,(ABS(TP))**2, (ABS(TN))**2



c$$$         SQRTARG=(((E**2)*(LAMBDAP**2)) - 4*E*LAMBDAP)
c$$$         IF (SQRTARG .GE. 0) THEN
c$$$            XIP=0.5*((E*LAMBDAP) - SQRT(SQRTARG) -2)
c$$$         ELSE
c$$$         XIP=(0.5*((E*LAMBDAP) -2), 0.5*SQRT(DABS(SQRTARG)))
         



         E=E+0.005
      END DO

 10   FORMAT (3ES15.5E4)
      END
