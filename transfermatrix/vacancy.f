
c$$$  Subroutine to generate vacancies with differant P(x) for a,b
c$$$  sublattice
      SUBROUTINE GENVAC(V, LIMX, LIMY, U, ALAT, BLAT)
      IMPLICIT NONE

      INTEGER LIMX, LIMY, I, J
      INTEGER SEED/-4/
      DOUBLE PRECISION V(LIMX,LIMY), U, ALAT, BLAT
      REAL GETRAND
      REAL RANDOM
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0

      CALL DLASET('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)

      DO I = 1, LIMX
         DO J = 1, LIMY
c$$$            CALL RANDOM_NUMBER(RANDOM)
            CALL SAFESEED(SEED)
            RANDOM=GETRAND()
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

c     Function to print V (Potential lattice) 
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

c     
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
      DOUBLE PRECISION TOTALV(TIMES), MEANV, DEVV, STD(TIMES), RTIMES

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
      RTIMES = TIMES
      DEVV = DEVV / SQRT(RTIMES)
      
      RETURN
      END

      SUBROUTINE INCREMENTCOND(MEANCOND, ERRCOND, NEWCOND, RNUM, NE)
      IMPLICIT NONE

      INTEGER NE, I
      DOUBLE PRECISION MEANCOND(NE), ERRCOND(NE), NEWCOND(NE),
     +                 SIGMA, RNUM
c      INTEGER NUM

      DO I = 1, NE
         SIGMA = 0.0
c         RNUM = NUM
         SIGMA = ERRCOND(I) * SQRT(RNUM - 1.0)
         SIGMA = SIGMA * SIGMA
         SIGMA = SIGMA + (MEANCOND(I) * MEANCOND(I))
         SIGMA = SIGMA * (RNUM - 1.0)
         SIGMA = SIGMA + (NEWCOND(I) * NEWCOND(I))
         SIGMA = SIGMA / RNUM

         MEANCOND(I) = MEANCOND(I) * (RNUM - 1.0)
         MEANCOND(I) = (MEANCOND(I) + NEWCOND(I)) / RNUM

         SIGMA = SIGMA - (MEANCOND(I) * MEANCOND(I))
         SIGMA = SQRT(SIGMA)
         ERRCOND(I) = SIGMA / SQRT(RNUM)
      END DO

      RETURN
      END

      SUBROUTINE FAKEGETTRANS(TVALS, NTVALS)
      IMPLICIT NONE

      INTEGER NTVALS, I
      DOUBLE PRECISION TVALS(NTVALS), RANDOM

      DO I = 1, NTVALS
         CALL RANDOM_NUMBER(RANDOM)
         TVALS(I) = RANDOM
      END DO

      RETURN
      END
