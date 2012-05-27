

      SUBROUTINE GENVAC(V, LIMX, LIMY, U, ALAT, BLAT)
      IMPLICIT NONE

      INTEGER LIMX, LIMY, I, J
      DOUBLE PRECISION V(LIMX,LIMY), U, ALAT, BLAT
      REAL RANDOM
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0

      CALL DLASET('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)

      DO I = 1, LIMX
         DO J = 1, LIMY
            CALL RANDOM_NUMBER(RANDOM)
            IF (MOD(I, 2) .EQ. MOD(J, 2)) THEN
               IF (RANDOM .LT. ALAT) THEN
                  V(I,J) = U
               END IF
            ELSE
               IF (RANDOM .LT. BLAT) THEN
                  V(I,J) = U
               END IF
            END IF
         END DO
      END DO

      RETURN
      END


      SUBROUTINE PRNTV(V, LIMX, LIMY)
      IMPLICIT NONE
      
      INTEGER LIMX, LIMY, I, J
      DOUBLE PRECISION V(LIMX,LIMY)

      DO J = 1, LIMY
         WRITE (*,10) (V(I,J), I = 1, LIMX)
      END DO

 10   FORMAT (10F8.1)

      RETURN
      END

      DOUBLE PRECISION FUNCTION ADDUPV(V, LIMX, LIMY)
      IMPLICIT NONE

      INTEGER I, J, LIMX, LIMY
      DOUBLE PRECISION V(LIMX,LIMY), T

      DO I = 1, LIMX
         DO J = 1, LIMY
            ADDUPV = ADDUPV + V(I,J)
         END DO
      END DO

      RETURN
      END

      SUBROUTINE MEANDEVV(TOTALV, TIMES, MEANV, DEVV)

      INTEGER TIMES, I
      DOUBLE PRECISION TOTALV(TIMES), MEANV, DEVV, STD(TIMES)

      DEVV = 0.0
      MEANV = 0.0
      DO I = 1, TIMES
         MEANV = MEANV + TOTALV(I)
      END DO
      MEANV = MEANV / TIMES
      DO I = 1, TIMES
         STD(I) = TOTALV(I) - MEANV
         STD(I) = STD(I) * STD(I)
         DEVV = DEVV + STD(I)
      END DO
      DEVV = SQRT(DEVV / TIMES)
      
      RETURN
      END
