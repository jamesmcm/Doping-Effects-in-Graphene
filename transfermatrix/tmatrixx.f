      DOUBLE PRECISION FUNCTION GETTRANSX(GAUGE, TVALS, LIMX, NSIZE,
     +   E, FLUX)
      IMPLICIT NONE
C     NSIZE = LIMY/2
      INTEGER I/1/, NSIZE, LIMX
      DOUBLE PRECISION TVALS(NSIZE), E, FLUX
      DOUBLE PRECISION CHECKUNI
      EXTERNAL CHECKUNI
      DOUBLE PRECISION CHECKUNI2
      EXTERNAL CHECKUNI2
      DOUBLE COMPLEX MODD(2*NSIZE, 2*NSIZE), MEVEN(2*NSIZE, 2*NSIZE),
     +               MULT(2*NSIZE, 2*NSIZE)
      DOUBLE COMPLEX A(NSIZE, NSIZE), B(NSIZE, NSIZE),
     +               C(NSIZE, NSIZE), D(NSIZE, NSIZE)
      DOUBLE COMPLEX T(NSIZE, NSIZE),    TTILDE(NSIZE, NSIZE),
     +               R(NSIZE, NSIZE),    RTILDE(NSIZE, NSIZE),
     +               TINC(NSIZE, NSIZE), TTILDEINC(NSIZE, NSIZE),
     +               RINC(NSIZE, NSIZE), RTILDEINC(NSIZE, NSIZE)
      DOUBLE COMPLEX U(NSIZE,NSIZE)
      CHARACTER GAUGE

c      CALL CALCMULT(LIMX, WRAPX, MODD, MEVEN, E, FLUX)
c$$$      CALL CALCMULTNEW(LIMX, WRAPX, MODD, E, FLUX, 1)
c$$$      CALL CALCMULTNEW(LIMX, WRAPX, MEVEN, E, FLUX, 2)
      CALL CALCMULTXY(NSIZE, MODD, E, FLUX, 1)
      CALL CALCMULTXY(NSIZE, MEVEN, E, FLUX, 2)

c$$$  CALCMULT fills MODD, MEVEN - do multiplication in main loop
c$$$  Must decide whether we want zig-zag or armchair edges
C     For now I have left it as before so I can compare results

      CALL SQUNIT (MULT, 2*NSIZE)

      CALL FILLUYX(U, FLUX, NSIZE)


      CALL SQUNIT (T, NSIZE)
      CALL SQUNIT (TTILDE, NSIZE)
      CALL SQZERO (R, NSIZE)
      CALL SQZERO (RTILDE, NSIZE)

      DO I = 1, LIMX
            IF (MOD(I,2) .EQ. 1) THEN
               CALL SQCOPY (MODD, MULT, 2*NSIZE) 
            ELSE
               CALL SQCOPY(MEVEN, MULT, 2*NSIZE)
            END IF

        
            CALL GENABCD(MULT, U, A,B,C,D, NSIZE)
            CALL GENTANDRINC(A, B, C, D,
     +                       TINC, RINC, TTILDEINC, RTILDEINC,
     +                       NSIZE)
            CALL UPDATETANDR(T,     R,    TTILDE,    RTILDE, 
     +                       TINC,  RINC, TTILDEINC, RTILDEINC, 
     +                       NSIZE)
      END DO
      CALL SQSVDVALS(T, TVALS, NSIZE)

c$$$  CheckUni2 is slightly faster --- AVS
      GETTRANSX = CHECKUNI2(T, R, TTILDE, RTILDE, NSIZE)
      RETURN
      END

      SUBROUTINE CALCMULTXY(NSIZE, MULT, E, FLUX, POS)
C     NOTE LIMY MUST BE PASSED AS CONSTANT PARAMETER FOR THIS, LIMY MUST BE EVEN
      IMPLICIT NONE
      INTEGER LIMX, POS
C     WRAPY ignored
      INTEGER I/1/, NSIZE
      DOUBLE PRECISION E, FLUX
      DOUBLE COMPLEX CNUM

      DOUBLE COMPLEX MULT(NSIZE*2, NSIZE*2), 
     +     N3(NSIZE, NSIZE), N2(NSIZE, NSIZE)

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

      MULT(1:NSIZE, 1:NSIZE)=-1*N3
      MULT(1:NSIZE, NSIZE+1:2*NSIZE)=E*N3
      MULT(NSIZE+1:2*NSIZE, 1:NSIZE)=-E*N3
      MULT(NSIZE+1:2*NSIZE, NSIZE+1:2*NSIZE)=(E*E*N3) - N2

      RETURN
      END
