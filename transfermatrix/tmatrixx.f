      DOUBLE PRECISION FUNCTION GETTRANSX(GAUGE, TVALS, LIMX, NSIZE,
     +   E, FLUX, V)
      IMPLICIT NONE
C     NSIZE = LIMY/2
      INTEGER I/1/, NSIZE, LIMX
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

      DOUBLE COMPLEX V(LIMX,2*NSIZE)
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
      INTEGER LIMX, POS
C     WRAPY ignored
      INTEGER I/1/, NSIZE
      DOUBLE PRECISION E, FLUX
      DOUBLE COMPLEX CNUM

      DOUBLE COMPLEX MULT(NSIZE*2, NSIZE*2), 
     +     N3(NSIZE, NSIZE), N2(NSIZE, NSIZE)
      DOUBLE COMPLEX V(LIMX,2*NSIZE)

C     Odd is defined for odd leftmost column, even for even leftmost column
      CALL SQZERO (N3, NSIZE)
      CALL SQZERO (N2, NSIZE)

      DO I = 1, NSIZE
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
c$$$      CALL INVERTMATRIX(N2, NSIZE)

C     I think that E will have to be changed to matrix of applicable E values to take V in to account
      MULT(1:NSIZE, 1:NSIZE)=-1*N3
      MULT(1:NSIZE, NSIZE+1:2*NSIZE)=E*N3
      MULT(NSIZE+1:2*NSIZE, 1:NSIZE)=-E*N3
      MULT(NSIZE+1:2*NSIZE, NSIZE+1:2*NSIZE)=(E*E*N3) - N2

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

      DOUBLE COMPLEX MULT(NSIZE*2, NSIZE*2), 
     +     N3(NSIZE, NSIZE), N2(NSIZE, NSIZE)
      DOUBLE COMPLEX V(LIMX,2*NSIZE)
C     Odd is defined for odd leftmost column, even for even leftmost column
      CALL SQZERO (N3, NSIZE)
      CALL SQZERO (N2, NSIZE)

      DO I = 1, NSIZE
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
      DO I = 1, NSIZE
        CALL ZPOLAR ( FLUX * 2 * I, TAU1(I))
c       The minus sign is due to the fact that tau1 is x->x+1 hopping,
c       while tau_2 is x+2->x+1 hopping
        CALL ZPOLAR (-FLUX * 2 * I, TAU2(I))
      ENDDO
      
C     To introduce V here, must express E as matrix of (E-V) values
      DO I = 1, NSIZE
        IN = I + NSIZE
        TI = TAU2(I)
        DO J = 1, NSIZE
          JN = J + NSIZE
          TJ = TAU1(J)
          N3IJ  = N3 (I, J)
          MULT(I,   J)  =     - N3IJ * TJ
          MULT(I,  JN)  =   E * N3IJ
          MULT(IN,  J)  = - E * N3IJ * TJ / TI
          MULT(IN, JN)  = (E * E * N3IJ  - N2(I, J)) / TI
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
