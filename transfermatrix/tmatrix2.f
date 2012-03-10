      SUBROUTINE CALCMULT2(LIMX, WRAPX, MULT, E, FLUX, POS)
      IMPLICIT NONE
      INTEGER LIMX, WRAPX, SZ/1/
      INTEGER I/1/, NEIGH/1/, POS
      DOUBLE PRECISION E, FLUX
      DOUBLE COMPLEX ZEROC / 0.0 /
      DOUBLE COMPLEX ZPLUS, ZMINUS, Z
c$$$  May need to move this
      DOUBLE COMPLEX MULT(2*LIMX, 2*LIMX)

c$$$  HAMMERTIME! Program terminates here if LIMX is odd
      IF ((WRAPX .EQ. 1)) THEN
         IF ((MOD(LIMX,2) .NE. 0)) THEN
            WRITE (*,*) 'ERROR: LIMX must be even for physical results
     + for wrapped X'
            STOP
         ENDIF
      ENDIF

      SZ = 2 * LIMX
      CALL ZLASET ('A', SZ, SZ, ZEROC, ZEROC, MULT, SZ)


C$$$ FIRST ROW IS EVEN - WRAPX MAKES NO DIFF, SECOND ROW NOT, ETC.
C$$$ - WHAT MATTERS IS WHICH ROW IT IS CENTRED ON
C$$$ THERE ARE 2 TRANSFER MATRICES TO GENERATE
C$$$ THERE ARE 4 BLOCK SUBMATRICES TO FILL
C$$$ MODD DOESN'T DEPEND ON XWRAPPING, MEVEN DOES.   '
      CALL ZPOLAR(FLUX*POS, ZPLUS)
      CALL ZPOLAR(-FLUX*POS, ZMINUS)
      DO I = 1, LIMX
C$$$ FILL TOP-RIGHT SUBMATRIX
         MULT(I, LIMX+I)=1
C$$$ FILL BOTTOM-LEFT SUBMATRIX
         MULT(I+LIMX, I)=-1
C$$$ FILL BOTTOM-RIGHT SUBMATRIX
         MULT(LIMX+I,LIMX+I)=E

c$$$  Double-check this multiplication analytically at some stage

C$$$  THE FOLLOWING CODE WAS MODIFIED --- AVS
C$$$  NEIGHBOURING SITE FOR ODD ROW, ON THE LEFT/RIGHT, DEPENDING ON I
         NEIGH = I + (2*MOD(I,2)-1)
         Z = ZPLUS**(2*MOD(I, 2)-1)
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
               MULT(LIMX + I, LIMX + NEIGH) = -1*Z
             ENDIF
         END IF
C$$$     NOW REPEAT THE SAME FOR EVEN ROWS, SWAPPING LEFT AND RIGHT
C$$$     AGAIN, THE CODE IS NOW RATHER UGLY.
         NEIGH = I - (2 * MOD(I, 2)  - 1)
         Z = ZPLUS**(2*MOD(I, 2)-1)
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
             NEIGH = LIMX
           ELSE
             NEIGH = MOD (NEIGH - 1, LIMX) + 1
           END IF
           IF (MOD(POS, 2).EQ.0) THEN
             MULT(LIMX + I, LIMX + NEIGH) = -1/Z
           ENDIF
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
      DOUBLE COMPLEX   ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX MULT(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX),
     +               ABCD(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX T(LIMX, LIMX),    TTILDE(LIMX, LIMX),
     +               R(LIMX, LIMX),    RTILDE(LIMX, LIMX),
     +               TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     +               RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX)

c$$$  CALCMULT fills MULT - do multiplication in main loop
c$$$  Must decide whether we want zig-zag or armchair edges
C     For now I have left it as before so I can compare results

      CALL ZLASET ('ALL', 2*LIMX, 2*LIMX, ZEROC, ONEC, MULT, 2*LIMX)

      CALL FILLOANDINVERT2(O, IO, LIMX, FLUX)

c      call zprintm (O, 2*LIMX, 'o: ')
c      call zprintm (IO, 2*LIMX, 'io: ')
      CALL ZLASET ('ALL', LIMX, LIMX, ZEROC, ONEC, T, LIMX)
      CALL ZLASET ('ALL', LIMX, LIMX, ZEROC, ZEROC, R, LIMX)
      CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
      CALL GENTANDRINC(LIMX, T, R, TTILDE, RTILDE, A, B, C, D)


      DO I = 1, LIMY
            CALL CALCMULT2(LIMX, WRAPX, MULT, E, FLUX, I)
c           CALL zprintm (MULT, 2*LIMX, 'm: ')
            CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
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

      SUBROUTINE FILLOANDINVERT2(O, IO, LIMX, FLUX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX)
      DOUBLE PRECISION SQRT05, FLUX
      DOUBLE COMPLEX ZISQRT05
      DOUBLE COMPLEX ZEROC/0.0/
      CALL ZLASET ('ALL', 2*LIMX, 2*LIMX, ZEROC, ZEROC, O, 2*LIMX)


c     It is slightly more efficient to calculate square root once
      SQRT05 = SQRT(0.5)
      ZISQRT05 = DCMPLX(0, SQRT05)
C$$$ GENERATE O-MATRIX
C$$$ O IS BLOCK MATRIX OF 1/SQRT(2) (1,1;I,-I)
         DO I = 1, LIMX
            O(I, I)=SQRT05
            O(I, LIMX+I)=SQRT05
            O(I+LIMX, I)=ZISQRT05
            O(I+LIMX, I+LIMX)=-ZISQRT05
c$$$  Hopefully this is correct - test analytically later
         ENDDO
         CALL ZCOPY(4*LIMX*LIMX, O, 1, IO, 1)
         CALL INVERTMATRIX(IO, 2*LIMX)
      RETURN
      END
