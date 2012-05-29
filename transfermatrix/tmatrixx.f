      DOUBLE PRECISION FUNCTION GETTRANSX(GAUGE, TVALS, LIMX, NSIZE,
     +   V, E, FLUX)
      IMPLICIT NONE
C     NSIZE = LIMY/2
      INTEGER I/1/, NSIZE, LIMX, LIMY
      DOUBLE PRECISION TVALS(NSIZE), E, FLUX
      DOUBLE PRECISION CHECKUNI
      EXTERNAL CHECKUNI
      DOUBLE PRECISION CHECKUNI2
      EXTERNAL CHECKUNI2
      DOUBLE COMPLEX MULT(2*NSIZE, 2*NSIZE)
      DOUBLE COMPLEX A(NSIZE, NSIZE), B(NSIZE, NSIZE),
     +               C(NSIZE, NSIZE), D(NSIZE, NSIZE)
      DOUBLE COMPLEX T(NSIZE, NSIZE),    TTILDE(NSIZE, NSIZE),
     +               R(NSIZE, NSIZE),    RTILDE(NSIZE, NSIZE),
     +               TINC(NSIZE, NSIZE), TTILDEINC(NSIZE, NSIZE),
     +               RINC(NSIZE, NSIZE), RTILDEINC(NSIZE, NSIZE)
      DOUBLE COMPLEX U(NSIZE,NSIZE)

      DOUBLE PRECISION V(LIMX,2*NSIZE)
      CHARACTER GAUGE

      CALL SQUNIT (MULT, 2*NSIZE)
      IF (GAUGE .EQ. 'X') THEN
          CALL FILLUXX(U, FLUX, NSIZE)
      ELSE
        CALL FILLUXY(U, FLUX, NSIZE)
      ENDIF

      CALL SQUNIT (T,      NSIZE)
      CALL SQUNIT (TTILDE, NSIZE)
      CALL SQZERO (R,      NSIZE)
      CALL SQZERO (RTILDE, NSIZE)

      DO I = 1, LIMX
            IF (GAUGE .EQ. 'X') THEN
              CALL CALCMULTXX(E, FLUX, I, MULT, LIMX, NSIZE, V)
            ELSE
              CALL CALCMULTXY(E, FLUX, I, MULT, LIMX, NSIZE, V)
            ENDIF
            CALL GENABCD(MULT, U, A,B,C,D, NSIZE)
            CALL GENTANDRINC(A, B, C, D,
     +                       TINC, RINC, TTILDEINC, RTILDEINC,
     +                       NSIZE)
            CALL UPDATETANDR(T,     R,    TTILDE,    RTILDE, 
     +                       TINC,  RINC, TTILDEINC, RTILDEINC, 
     +                       NSIZE)
      END DO
      CALL SQSVDVALS(T, TVALS, NSIZE)

      GETTRANSX = CHECKUNI2(T, R, TTILDE, RTILDE, NSIZE)
      RETURN
      END

      SUBROUTINE CALCMULTXY(E, FLUX, POS, MULT, LIMX, NSIZE, V)
C     NOTE LIMY MUST BE PASSED AS CONSTANT PARAMETER FOR THIS, LIMY MUST BE EVEN
      IMPLICIT NONE
      INTEGER LIMX, POS, IN, JN
C     WRAPY ignored
      INTEGER I/1/, NSIZE
      DOUBLE PRECISION E, FLUX
      DOUBLE COMPLEX CNUM, TAU1, TI, N3IJ, TJ

      DOUBLE COMPLEX MULT(NSIZE*2, NSIZE*2), 
     +     N3(NSIZE, NSIZE), N2(NSIZE, NSIZE)
      DOUBLE PRECISION V(LIMX,2*NSIZE)
      DOUBLE PRECISION EV2(NSIZE), EV3(NSIZE)
      DOUBLE PRECISION EV2JJ, EV3II
      INTEGER J, K

c$$$C     E-V values stored in column vector - in future
c$$$      DATA EV2/NSIZE*0.0/
c$$$      DATA EV3/NSIZE*0.0/


C     Odd is defined for odd leftmost column, even for even leftmost column
      CALL SQZERO (N3, NSIZE)
      CALL SQZERO (N2, NSIZE)

      CALL ZPOLAR (FLUX * (POS), TAU1)

      IF (MOD(POS,2) .EQ. 1) THEN
         CALL ZPOLAR (FLUX * (POS), TAU1)
      ENDIF
      IF (MOD(POS,2) .EQ. 0) THEN
         CALL ZPOLAR (FLUX * (-1) * (POS), TAU1)
      ENDIF

      DO I = 1, NSIZE

         N3(I, I)=1.0/TAU1
         N2(I, I)=1*TAU1

C     ((-2*MOD(I,2))+1)
c     The minus sign is due to the fact that tau1 is x->x+1 hopping,
c     while tau_2 is x+2->x+1 hopping

         IF (MOD(POS,2) .EQ. 0) THEN
            EV2(I)=E-V(POS, 2*I) 
            EV3(I)=E-V(POS, (2*I)-1)
         ELSE
            EV3(I)=E-V(POS, 2*I) 
            EV2(I)=E-V(POS, (2*I)-1) 
         ENDIF


C     E2 psi2 = N3 psi3 + psi1
C     E3 psi3 = N2 psi2 + psi4
         IF (I .NE. 1) THEN
            IF (MOD(POS,2) .EQ. 0) THEN
               N2(I,I-1)=1*TAU1
            ELSE
               N3(I,I-1)=1*TAU1
            ENDIF
         ENDIF
         IF (I .NE. NSIZE) THEN
            IF (MOD(POS,2) .EQ. 0) THEN
               N3(I,I+1)=1.0/TAU1
            ELSE
               N2(I,I+1)=1.0/TAU1
            ENDIF
         ENDIF
      ENDDO

c$$$      CALL ZPRINTM(N2,NSIZE,"N2 ")
c$$$      PRINT *, "---"
c$$$      CALL ZPRINTM(N3,NSIZE,"N3 ")
C     M={{-N3^-1, E*N3^-1},{-E*N3^-1, E^2 N3^-1 - N2^-1}}

c$$$      CALL INVERTMATRIX(N2, NSIZE)
C     I think that E will have to be changed to matrix of applicable E values to take V in to account
C     So we have {E-V_1, 0, 0...}{0, E-V_2, 0, ...} etc. 
C     N3 relates to E2, i.e. of second column
C     N2 relates to E3, i.e. of first column
C     Use looping


c$$$      DO I=1, NSIZE:
      CALL INVERTMATRIX(N3, NSIZE)
      
c$$$      MULT(1:NSIZE, 1:NSIZE)=-1*N3
c$$$      CALL SQDOT(N3EV2, N3, EV2, NSIZE)
c$$$      MULT(1:NSIZE, NSIZE+1:2*NSIZE)=N3EV2
c$$$      CALL SQDOT(EV3N3, EV3, N3, NSIZE)
c$$$      MULT(NSIZE+1:2*NSIZE, 1:NSIZE)=-1*EV3N3
c$$$      CALL SQDOT(EV3N3EV2, EV3N3, EV2, NSIZE)
c$$$
c$$$      MULT(NSIZE+1:2*NSIZE, NSIZE+1:2*NSIZE)=EV3N3EV2 - N2


      DO I = 1, NSIZE
         IN = I + NSIZE
         EV3II=EV3(I)

         DO J= 1, NSIZE
            JN = J + NSIZE

            N3IJ = N3(I,J)
            EV2JJ=EV2(J)
            MULT(I,   J)  =     - N3IJ 
            MULT(I,  JN)  =   EV2JJ * N3IJ  
            MULT(IN,  J)  = - EV3II * N3IJ 
            MULT(IN, JN)  = (EV3II * EV2JJ * N3IJ  - (N2(I, J)))
         ENDDO
      ENDDO
            

c$$$c$$$C     Old code
c$$$      MULT(1:NSIZE, 1:NSIZE)=-1*N3
c$$$      MULT(1:NSIZE, NSIZE+1:2*NSIZE)=E*N3
c$$$      MULT(NSIZE+1:2*NSIZE, 1:NSIZE)=-E*N3
c$$$      MULT(NSIZE+1:2*NSIZE, NSIZE+1:2*NSIZE)=(E*E*N3) - N2

c      PRINT*, "POS:", POS
c      CALL PRINTM(EV2, NSIZE, 'EV2: ')
c      PRINT *, "___"
c      CALL PRINTM(EV3, NSIZE, 'EV3: ')
c      PRINT *, "---"
      RETURN
      END
      
      SUBROUTINE CALCMULTXX(E, FLUX, POS, MULT, LIMX, NSIZE, V)
C     NOTE LIMY MUST BE PASSED AS CONSTANT PARAMETER FOR THIS, LIMY MUST BE EVEN
      IMPLICIT NONE
      INTEGER LIMX, POS
C     WRAPY ignored
      INTEGER I, J, IN, JN, NSIZE
      DOUBLE PRECISION E, FLUX
      DOUBLE COMPLEX CNUM, TAU1(NSIZE), TAU2(NSIZE), TI, N3IJ, TJ
      DOUBLE COMPLEX EV3II, EV2JJ

      DOUBLE COMPLEX MULT(NSIZE*2, NSIZE*2), 
     +     N3(NSIZE, NSIZE), N2(NSIZE, NSIZE)
      DOUBLE PRECISION V(LIMX,2*NSIZE), EV2(NSIZE), EV3(NSIZE)
C     Odd is defined for odd leftmost column, even for even leftmost column
      CALL SQZERO (N3, NSIZE)
      CALL SQZERO (N2, NSIZE)


      DO I = 1, NSIZE

         CALL ZPOLAR ( FLUX * 2 * I, TAU1(I))
c       The minus sign is due to the fact that tau1 is x->x+1 hopping,
c       while tau_2 is x+2->x+1 hopping
         CALL ZPOLAR (-FLUX * 2 * I, TAU2(I))

         IF (MOD(POS,2) .EQ. 0) THEN
            EV2(I)=E-V(POS, 2*I) 
            EV3(I)=E-V(POS, (2*I)-1)
         ELSE
            EV3(I)=E-V(POS, 2*I) 
            EV2(I)=E-V(POS, (2*I)-1) 
         ENDIF

         N3(I, I)=1
         N2(I, I)=1
         IF (I .NE. 1) THEN
            IF (MOD(POS,2) .EQ. 0) THEN
               N2(I,I-1)=1
            ELSE
               N3(I,I-1)=1
            ENDIF
         ENDIF
         IF (I .NE. NSIZE) THEN
            IF (MOD(POS,2) .EQ. 0) THEN
               N3(I,I+1)=1
            ELSE
               N2(I,I+1)=1
            ENDIF
         ENDIF
      ENDDO

C     M={{-N3^-1, E*N3^-1},{-E*N3^-1, E^2 N3^-1 - N2^-1}}

      CALL INVERTMATRIX(N3, NSIZE)
c
c     Note: strictly speaking, the phase should be different by one
c     FLUX between even and odd neighbouring slices. For an open boundary,
c     however, we can make a gauge transformation which eliminates the
c     extra phase. Indeed, after we do this, the fluxes through every loop
c     remain the same. If, for some reason, the offset is to be
c     reintroduced later, the U matrix has to be changed accordingly.
c 

c$$$      MULT(1:NSIZE, 1:NSIZE)=-1*N3
c$$$      CALL SQDOT(N3EV2, N3, EV2, NSIZE)
c$$$      MULT(1:NSIZE, NSIZE+1:2*NSIZE)=N3EV2
c$$$      CALL SQDOT(EV3N3, EV3, N3, NSIZE)
c$$$      MULT(NSIZE+1:2*NSIZE, 1:NSIZE)=-1*EV3N3
c$$$      CALL SQDOT(EV3N3EV2, EV3N3, EV2, NSIZE)
c$$$
c$$$      MULT(NSIZE+1:2*NSIZE, NSIZE+1:2*NSIZE)=EV3N3EV2 - N2


C     To introduce V here, must express E as matrix of (E-V) values
      DO I = 1, NSIZE
        IN = I + NSIZE
        TI = TAU2(I)
        EV3II=EV3(I)
        DO J = 1, NSIZE
          JN = J + NSIZE
          TJ = TAU1(J)
          N3IJ  = N3 (I, J)
          EV2JJ= EV2(J)
          MULT(I,   J)  =     - N3IJ * TJ
          MULT(I,  JN)  =   EV2JJ * N3IJ
          MULT(IN,  J)  = - EV3II * N3IJ * TJ / TI
          MULT(IN, JN)  = (EV3II * EV2JJ * N3IJ  - N2(I, J)) / TI
        ENDDO
      ENDDO
c      MULT(1:NSIZE, 1:NSIZE)=-1*N3
c      MULT(1:NSIZE, NSIZE+1:2*NSIZE)=E*N3
c      MULT(NSIZE+1:2*NSIZE, 1:NSIZE)=-E*N3
c      MULT(NSIZE+1:2*NSIZE, NSIZE+1:2*NSIZE)=(E*E*N3) - N2

      RETURN
      END

      SUBROUTINE FILLUXY(U, FLUX, LIMX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX U(LIMX, LIMX)
      DOUBLE PRECISION FLUX
      DOUBLE COMPLEX ZI/(0.0, 1.0)/
     
      CALL SQUNITZ (U, ZI, LIMX)

      RETURN
      END

      SUBROUTINE FILLUXX(U, FLUX, LIMX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX U(LIMX,LIMX)
      DOUBLE PRECISION FLUX
      DOUBLE COMPLEX   CNUM
      DOUBLE COMPLEX ZI/(0.0, 1.0)/      

      CALL SQZERO (U, LIMX)
      DO I = 1, LIMX
         CALL ZPOLAR(FLUX * 2 * I, CNUM)
         U(I, I) = ZI * CNUM
      ENDDO
     
      RETURN
      END
