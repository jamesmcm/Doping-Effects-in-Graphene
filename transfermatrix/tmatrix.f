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
     +               C(LIMX, LIMX), D(LIMX, LIMX),
     +               ABCD(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX T(LIMX, LIMX),    TTILDE(LIMX, LIMX),
     +               R(LIMX, LIMX),    RTILDE(LIMX, LIMX),
     +               TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     +               RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX),
     + AOLD(LIMX,LIMX),BOLD(LIMX,LIMX),COLD(LIMX,LIMX),DOLD(LIMX,LIMX),
     + U(LIMX,LIMX)

      CALL CALCMULT(LIMX, WRAPX, MODD, MEVEN, E, FLUX)
c$$$  CALCMULT fills MODD, MEVEN - do multiplication in main loop
c$$$  Must decide whether we want zig-zag or armchair edges
C     For now I have left it as before so I can compare results

      CALL SQUNIT (MULT, 2*LIMX)

      CALL FILLU(U, LIMX, FLUX)

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
            CALL UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     +           RTILDE, LIMX, RINC)
      END DO
      CALL SV_DECOMP(LIMX, T, TVALS)

c$$$  CheckUni2 is slightly faster --- AVS
      GETTRANS = CHECKUNI2(LIMX,T,R,TTILDE,RTILDE)
      RETURN
      END

c$$$  AIM TO CHANGE TO FILL(U,LIMX,FLUX)
      SUBROUTINE FILLU(U, LIMX, FLUX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX U(LIMX,LIMX)
      DOUBLE PRECISION SQRT05, FLUX
      DOUBLE COMPLEX ZISQRT05, CNUM
      DOUBLE COMPLEX ZI/(0.0, 1.0)/      
      CALL SQZERO (U, LIMX)
      DO I = 1, LIMX
         CALL ZPOLAR(FLUX*I, CNUM)
         U(I,I)=ZI*CNUM
      ENDDO
      CALL ZPRINTM(U, LIMX, "U :")

c$$$      CALL SQUNITZ(U, ZI, LIMX)
      


     
      RETURN
      END
