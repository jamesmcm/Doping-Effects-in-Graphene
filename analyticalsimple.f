c$$$Calculate xi for both lambda values (may be complex)
c$$$Calculate coefficients for lambda values and xi values - note lambda must match xi
c$$$Sub coefficients in to T formula
c$$$Take modulus and square
c$$$ REMEMBER THAT THE REAL N OF BLOCKS IS TWICE WHAT IS USED

      PROGRAM PLOTNOWRAP
      IMPLICIT NONE
      DOUBLE PRECISION E/-5.0/, LAMBDA, TP, TN, SINP, COSP, PHI, SQRTARG
     +     , REALXI
      INTEGER F/1/, N/1/
      CHARACTER*3 VALUE 


c$$$      set N to command line argument
      CALL GETARG(1, VALUE)
      READ(UNIT=VALUE, FMT=*) N
      
      DO F = 1, 2002    
c$$$  Do Lmabda = E+1
c$$$  Check if xi complex or real
c$$$  use appropriate case to calculate T
c$$$  repeat with lambda=E-1
c$$$  increment E
c$$$         LAMBDA=E+1.0
         LAMBDA=E
         SQRTARG=((((E**2)*(LAMBDA**2))/4.0) - (E*LAMBDA))
         IF (SQRTARG .GE. 0) THEN
            REALXI=((E*LAMBDA)/2.0) - SQRT(SQRTARG) - 1.0
            TP=4*((REALXI**N + REALXI**(-1.0*N))**2 + (((E+LAMBDA)/(
     +           REALXI-(REALXI**(-1.0))))*(REALXI**N - REALXI**(-1.0*N
     +           )))**2)/((REALXI**N + REALXI**(-1.0*N))**2 + (((E+
     +           LAMBDA)/(REALXI-(REALXI**(-1.0)))*(REALXI**N - REALXI**
     +           (-1.0*N))))**2)**2
         ELSE
            SINP=-1.0*SQRT((E*LAMBDA) - ((E**2 * LAMBDA**2)/4.0))
            COSP=((E*LAMBDA)/2.0) -1.0
            PHI=ATAN2(SINP, COSP)
            TP=4*(((2*COS(N*PHI))**2) + (((-E-LAMBDA)/(SIN(PHI))) 
     +           *SIN(N*PHI))**2) / (((4*(COS(N*PHI)**2)) + (((-E-LAMBDA
     +           )/SIN(PHI))*SIN(N*PHI))**2)**2)
         ENDIF
c$$$         LAMBDA=E-1.0
         LAMBDA=E
         SQRTARG=((((E**2)*(LAMBDA**2))/4.0) - (E*LAMBDA))
         IF (SQRTARG .GE. 0) THEN
            REALXI=((E*LAMBDA)/2.0) - SQRT(SQRTARG) - 1.0
            TN=4*((REALXI**N + REALXI**(-1.0*N))**2 + (((E+LAMBDA)/(
     +           REALXI-(REALXI**(-1.0))))*(REALXI**N - REALXI**(-1.0*N
     +           )))**2)/((REALXI**N + REALXI**(-1.0*N))**2 + (((E+
     +           LAMBDA)/(REALXI-(REALXI**(-1.0)))*(REALXI**N - REALXI**
     +           (-1.0*N))))**2)**2
         ELSE
            SINP=-1.0*SQRT((E*LAMBDA) - ((E**2 * LAMBDA**2)/4.0))
            COSP=((E*LAMBDA)/2.0) -1.0
            PHI=ATAN2(SINP, COSP)
            TN=4*(((2*COS(N*PHI))**2) + ((((-E-LAMBDA)/(SIN(PHI))) * 
     +           SIN(N*PHI))**2)) / ((4*((COS(N*PHI)**2)) + (((-E-LAMBDA
     +           )/SIN(PHI))*SIN(N*PHI))**2)**2)
         ENDIF

         WRITE (*,10) E,TP+TN



c$$$         SQRTARG=(((E**2)*(LAMBDAP**2)) - 4*E*LAMBDAP)
c$$$         IF (SQRTARG .GE. 0) THEN
c$$$            XIP=0.5*((E*LAMBDAP) - SQRT(SQRTARG) -2)
c$$$         ELSE
c$$$         XIP=(0.5*((E*LAMBDAP) -2), 0.5*SQRT(DABS(SQRTARG)))
         



         E=E+0.005
      END DO

 10   FORMAT (3ES15.5E4)
      END

      
c$$$      DOUBLE PRECISION FUNCTION CONDUCTANCE (TVALS, LIMX)
c$$$      INTEGER LIMX
c$$$      DOUBLE PRECISION TVALS(LIMX)
c$$$      DOUBLE PRECISION DDOT      
c$$$      CONDUCTANCE = DDOT (LIMX, TVALS, 1, TVALS, 1)
c$$$      RETURN
c$$$      END
