c$$$ 
      SUBROUTINE CALCMULT(LIMX, WRAPX, MODD, MEVEN, E, FLUX)
      IMPLICIT NONE
      INTEGER LIMX, WRAPX, SZ/1/
      INTEGER I/1/, NEIGH/1/
      DOUBLE PRECISION E, FLUX
      DOUBLE COMPLEX CNUM
c$$$  May need to move this
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX)

c$$$  HAMMERTIME! Program terminates here if LIMX is odd
      IF ((WRAPX .EQ. 1)) THEN
         IF ((MOD(LIMX,2) .NE. 0)) THEN
            WRITE (*,*) 'ERROR: LIMX must be even for physical results
     + for wrapped X'
            STOP
         ENDIF
      ENDIF

      SZ = 2 * LIMX
      CALL SQZERO (MODD,  2*LIMX)
      CALL SQZERO (MEVEN, 2*LIMX)


C$$$ FIRST ROW IS EVEN - WRAPX MAKES NO DIFF, SECOND ROW NOT, ETC.
C$$$ - WHAT MATTERS IS WHICH ROW IT IS CENTRED ON
C$$$ THERE ARE 2 TRANSFER MATRICES TO GENERATE
C$$$ THERE ARE 4 BLOCK SUBMATRICES TO FILL
C$$$ MODD DOESN'T DEPEND ON XWRAPPING, MEVEN DOES.   '
      DO I = 1, LIMX
C$$$ FILL TOP-RIGHT SUBMATRIX
         MODD(I, LIMX+I)=1
         MEVEN(I, LIMX+I)=1
C$$$ FILL BOTTOM-LEFT SUBMATRIX
         CALL ZPOLAR(2*FLUX*I, CNUM)
         MODD(I+LIMX, I)=-1*CNUM
         MEVEN(I+LIMX, I)=-1*CNUM
C$$$ FILL BOTTOM-RIGHT SUBMATRIX
         CALL ZPOLAR(FLUX*I, CNUM)
         MODD(LIMX+I,LIMX+I)=E*CNUM
         MEVEN(LIMX+I,LIMX+I)=E*CNUM

c$$$  Double-check this multiplication analytically at some stage

C$$$  THE FOLLOWING CODE WAS MODIFIED --- AVS
C$$$  NEIGHBOURING SITE FOR ODD ROW, ON THE LEFT/RIGHT, DEPENDING ON I
         NEIGH = I + (2*MOD(I,2)-1)
C$$$  NEIGHBOUR CAN BE < 0, OR > LIMX. IF WRAPX IS TRUE, THIS INDICATES
C$$$  A VALID SITE. THE FOLLOWING CODE IS A BIT UGLY, AS I AM NOT SURE
C$$$  WHAT IS MOD(-1, N) IN FORTRAN.
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
               NEIGH = NEIGH + LIMX
             ELSE
               NEIGH = MOD (NEIGH - 1, LIMX) + 1
             ENDIF
             MODD(LIMX + I, LIMX + NEIGH) = -1*CNUM
         END IF
C$$$     NOW REPEAT THE SAME FOR EVEN ROWS, SWAPPING LEFT AND RIGHT
C$$$     AGAIN, THE CODE IS NOW RATHER UGLY.
         NEIGH = I - (2 * MOD(I, 2)  - 1)
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
             NEIGH = LIMX
           ELSE
             NEIGH = MOD (NEIGH - 1, LIMX) + 1
           END IF
           MEVEN(LIMX + I, LIMX + NEIGH) = -1*CNUM
         END IF
      END DO
c$$$  Originally the first M matrix was set here
      RETURN
      END

      SUBROUTINE CALCMULTYX(LIMX, WRAPX, MULT, E, FLUX, POS)
      IMPLICIT NONE
      INTEGER LIMX, WRAPX
      DOUBLE COMPLEX MULT(2*LIMX, 2*LIMX)
      DOUBLE PRECISION E, FLUX
      INTEGER POS
      INTEGER I/1/, NEIGH/1/ 
      DOUBLE COMPLEX CNUM

c$$$  HAMMERTIME! Program terminates here if LIMX is odd
      IF ((WRAPX .EQ. 1)) THEN
         IF ((MOD(LIMX,2) .NE. 0)) THEN
            WRITE (*,*) 'ERROR: LIMX must be even for physical results
     + for wrapped X'
            STOP
         ENDIF
      ENDIF

      CALL SQZERO (MULT,  2*LIMX)


C$$$ FIRST ROW IS EVEN - WRAPX MAKES NO DIFF, SECOND ROW NOT, ETC.
C$$$ - WHAT MATTERS IS WHICH ROW IT IS CENTRED ON
C$$$ THERE ARE 2 TRANSFER MATRICES TO GENERATE
C$$$ THERE ARE 4 BLOCK SUBMATRICES TO FILL
C$$$ MODD DOESN'T DEPEND ON XWRAPPING, MEVEN DOES.   '
      DO I = 1, LIMX
C$$$ FILL TOP-RIGHT SUBMATRIX
         MULT(I, LIMX+I) = 1
C$$$ FILL BOTTOM-LEFT SUBMATRIX
         CALL ZPOLAR(2*FLUX*I, CNUM)
         MULT(I + LIMX, I) = -CNUM
C$$$ FILL BOTTOM-RIGHT SUBMATRIX
         CALL ZPOLAR(FLUX*I, CNUM)
         MULT(LIMX+I, LIMX+I) = E*CNUM

c$$$  Double-check this multiplication analytically at some stage

C$$$  THE FOLLOWING CODE WAS MODIFIED --- AVS
C$$$  NEIGHBOURING SITE FOR ODD ROW, ON THE LEFT/RIGHT, DEPENDING ON I
         NEIGH = I + (2*MOD(I,2)-1)
C$$$  NEIGHBOUR CAN BE < 0, OR > LIMX. IF WRAPX IS TRUE, THIS INDICATES
C$$$  A VALID SITE. THE FOLLOWING CODE IS A BIT UGLY, AS I AM NOT SURE
C$$$  WHAT IS MOD(-1, N) IN FORTRAN.
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
               NEIGH = NEIGH + LIMX
             ELSE
               NEIGH = MOD (NEIGH - 1, LIMX) + 1
             ENDIF
             IF (MOD(POS, 2).EQ.1) THEN
                 MULT(LIMX + I, LIMX + NEIGH) = - CNUM
             END IF
         END IF
C$$$     NOW REPEAT THE SAME FOR EVEN ROWS, SWAPPING LEFT AND RIGHT
C$$$     AGAIN, THE CODE IS NOW RATHER UGLY.
         NEIGH = I - (2 * MOD(I, 2)  - 1)
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
             NEIGH = LIMX
           ELSE
             NEIGH = MOD (NEIGH - 1, LIMX) + 1
           END IF
           IF (MOD (POS, 2).EQ.0) THEN
               MULT(LIMX + I, LIMX + NEIGH) = - CNUM
           END IF
         END IF
      END DO
      RETURN
      END


      DOUBLE PRECISION FUNCTION GETTRANS(TVALS, LIMX, LIMY,
     +   E, FLUX, WRAPX)
      IMPLICIT NONE
      INTEGER I/1/, LIMY, WRAPX, LIMX
      DOUBLE PRECISION TVALS(LIMX), E, FLUX
      DOUBLE PRECISION CHECKUNI
      EXTERNAL CHECKUNI
      DOUBLE PRECISION CHECKUNI2
      EXTERNAL CHECKUNI2
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX),
     +               MULT(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX)
      DOUBLE COMPLEX T(LIMX, LIMX),    TTILDE(LIMX, LIMX),
     +               R(LIMX, LIMX),    RTILDE(LIMX, LIMX),
     +               TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     +               RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)
      DOUBLE COMPLEX U(LIMX,LIMX)

c      CALL CALCMULT(LIMX, WRAPX, MODD, MEVEN, E, FLUX)
      CALL CALCMULTYX(LIMX, WRAPX, MODD, E, FLUX, 1)
      CALL CALCMULTYX(LIMX, WRAPX, MEVEN, E, FLUX, 2)
c$$$  CALCMULT fills MODD, MEVEN - do multiplication in main loop
c$$$  Must decide whether we want zig-zag or armchair edges
C     For now I have left it as before so I can compare results

      CALL SQUNIT (MULT, 2*LIMX)

      CALL FILLUYX(U, LIMX, FLUX)

      CALL SQUNIT (T, LIMX)
      CALL SQUNIT (TTILDE, LIMX)
      CALL SQZERO (R, LIMX)
      CALL SQZERO (RTILDE, LIMX)

      DO I = 1, LIMY
            IF (MOD(I,2) .EQ. 1) THEN
               CALL SQCOPY (MODD, MULT, 2*LIMX) 
            ELSE
               CALL SQCOPY(MEVEN, MULT, 2*LIMX)
            END IF

        
            CALL GENABCD(LIMX,MULT,A,B,C,D,U)
            CALL GENTANDRINC(LIMX, TINC, RINC, TTILDEINC, RTILDEINC,
     +           A, B, C,D)
            CALL UPDATETANDR(T,     R,    TTILDE,    RTILDE, 
     +                       TINC,  RINC, TTILDEINC, RTILDEINC, 
     +                       LIMX)
      END DO
      CALL SV_DECOMP(LIMX, T, TVALS)

c$$$  CheckUni2 is slightly faster --- AVS
      GETTRANS = CHECKUNI2(LIMX,T,R,TTILDE,RTILDE)
      RETURN
      END

      DOUBLE PRECISION FUNCTION GETTRANSXY(TVALS, LIMX, NSIZE,
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

c      CALL CALCMULT(LIMX, WRAPX, MODD, MEVEN, E, FLUX)
c$$$      CALL CALCMULTNEW(LIMX, WRAPX, MODD, E, FLUX, 1)
c$$$      CALL CALCMULTNEW(LIMX, WRAPX, MEVEN, E, FLUX, 2)
      CALL CALCMULTXY(NSIZE, MODD, E, FLUX, 1)
      CALL CALCMULTXY(NSIZE, MEVEN, E, FLUX, 2)

c$$$  CALCMULT fills MODD, MEVEN - do multiplication in main loop
c$$$  Must decide whether we want zig-zag or armchair edges
C     For now I have left it as before so I can compare results

      CALL SQUNIT (MULT, 2*NSIZE)

      CALL FILLUYX(U, NSIZE, FLUX)

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

        
            CALL GENABCD(NSIZE,MULT,A,B,C,D,U)
            CALL GENTANDRINC(NSIZE, TINC, RINC, TTILDEINC, RTILDEINC,
     +           A, B, C,D)
            CALL UPDATETANDR(T,     R,    TTILDE,    RTILDE, 
     +                       TINC,  RINC, TTILDEINC, RTILDEINC, 
     +                       NSIZE)
      END DO
      CALL SV_DECOMP(NSIZE, T, TVALS)

c$$$  CheckUni2 is slightly faster --- AVS
      GETTRANSXY = CHECKUNI2(NSIZE,T,R,TTILDE,RTILDE)
      RETURN
      END

c$$$  AIM TO CHANGE TO FILL(U,LIMX,FLUX)
      SUBROUTINE FILLUYX(U, LIMX, FLUX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX U(LIMX,LIMX)
      DOUBLE PRECISION FLUX
      DOUBLE COMPLEX   CNUM
      DOUBLE COMPLEX ZI/(0.0, 1.0)/      

      CALL SQZERO (U, LIMX)
      DO I = 1, LIMX
         CALL ZPOLAR(FLUX*I, CNUM)
         U(I,I)=ZI*CNUM
      ENDDO
     
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

      SUBROUTINE FILLUYY(U, FLUX, LIMX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX U(LIMX, LIMX)
      DOUBLE PRECISION FLUX
      DOUBLE COMPLEX ZI/(0.0, 1.0)/
     
      CALL SQUNITZ (U, ZI, LIMX)

      RETURN
      END
