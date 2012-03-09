      SUBROUTINE INVERTMATRIX(MATRIX, LIMX)
      IMPLICIT NONE
      INTEGER S
      INTEGER LIMX, PIVOT(LIMX, LIMX)
      DOUBLE COMPLEX MATRIX(LIMX, LIMX), WORK(LIMX*LIMX)
      CALL ZGETRF(LIMX, LIMX, MATRIX, LIMX, PIVOT, S)
      CALL ZGETRI(LIMX, MATRIX, LIMX, PIVOT, WORK, LIMX*LIMX, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'NON-INVERTABLE MATRIX WITH S=', S
         STOP
      END IF

      RETURN
      END

      SUBROUTINE SV_DECOMP(LIMX, MATRIX, OUTPUTS)

      INTEGER LIMX, MSIZE, S
      DOUBLE PRECISION SVALS(LIMX), OUTPUTS(LIMX), RWORK(5*LIMX)
      DOUBLE COMPLEX MATRIX(LIMX, LIMX), TEMP2(LIMX, LIMX),
     + SVCPY(LIMX, LIMX), WORK(4*LIMX*LIMX)

      MSIZE=4*LIMX*LIMX

C$$$ MAKE COPY OF MATRIX FOR SVD SINCE IT IS DESTROYED
C$$$      SVCPY=MATRIX
      CALL ZCOPY(LIMX*LIMX, MATRIX, 1, SVCPY, 1)
      CALL ZGESVD('N', 'N', LIMX, LIMX, SVCPY, LIMX, SVALS, TEMP2,
     + LIMX, TEMP2, LIMX , WORK, MSIZE, RWORK, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'SVD FAILED WITH S=', S
         STOP
      END IF

C$$$      OUTPUTS=SVALS
      CALL DCOPY(LIMX, SVALS, 1, OUTPUTS, 1)

      RETURN
      END

      DOUBLE PRECISION FUNCTION CHECKUNI(LIMX, T,R,TTILDE,RTILDE)
      IMPLICIT NONE

      INTEGER LIMX, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX T(LIMX,LIMX), BETA/-1/,ALPHA/1/,
     + R(LIMX,LIMX),TTILDE(LIMX,LIMX),RTILDE(LIMX,LIMX),
     + U(LIMX*2,LIMX*2), CHECK(LIMX*2,LIMX*2)


C$$$ TEST CASES OF T,R,T~ R~
C$$$ FILLS A UNIT MATRIX
      CALL ZLASET ('ALL', 2*LIMX, 2*LIMX, ZEROC, ONEC, CHECK, 2*LIMX)
      DO X=1, LIMX
         DO Y=1, LIMX
C$$$  TOP LEFT
            U(X,Y)=T(X,Y)
C$$$  BOTTOM LEFT
            U(X+LIMX,Y)=R(X,Y)
C$$$  TOP RIGHT
            U(X,Y+LIMX)=RTILDE(X,Y)
C$$$  BOTTOM RIGHT
            U(X+LIMX,Y+LIMX)=TTILDE(X,Y)
         END DO
      END DO

C$$$  ZGEMM HAS INBUILT FUNCTION TO FIND U**H
C$$   CALL PRINTT (U, 2*LIMX, 'U1 ')
      CALL ZGEMM('N', 'C', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, U,
     +     2*LIMX, U, 2*LIMX, BETA, CHECK, 2*LIMX)
C$$$  ZLANGE FINDS MATRIX NORM
C$$   CALL PRINTT (CHECK, 2 * LIMX, 'C2 ')
      CHECKUNI = ZLANGE('F', 2*LIMX, 2*LIMX, CHECK, 2*LIMX)
C$$   WRITE (*, *) 'CHECKUNI:', CHECKUNI
      RETURN
      END

C$$$  ROUTINE TO CHECK FOR UNITARY MATRICES

      DOUBLE PRECISION FUNCTION CHECKUNI2(LIMX, T,R,TTILDE,RTILDE)
      IMPLICIT NONE

      INTEGER LIMX
      DOUBLE PRECISION ZLANGE
      DOUBLE PRECISION DNORM1, DNORM2, DNORM3, DNORM
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX AU/0.0/, BU/1.0/, ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX T(LIMX,LIMX), R(LIMX,LIMX),
     +               TTILDE(LIMX,LIMX),RTILDE(LIMX,LIMX),
     +               CK(LIMX, LIMX)


C$$$ TEST CASES OF T,R,T~ R~

C$$$ FILLS A UNIT MATRIX

C      DO X = 1, 2*LIMX
C         DO Y = 1, 2*LIMX
C            CHECK(X, Y) = (0.0, 0.0)
C         END DO
C       CHECK(X, X) = 1.0
C      END DO

C$$$ ZGEMM HAS INBUILT FUNCTION TO FIND U**H
C$$   CALL PRINTT (U, 2*LIMX, 'U1 ')

c     Unitarity can be also checked using subblocks only
c     The definition of unitarity can be easily cast into the form
c     |T|^2 + |R|^2 = 1, T*R~ + R* T~ = 0

c     Check that T * T + R* R = 1   (* is the Hermitian conjugate)
c     Here I use zherk, which does T^+T
c     Beware: it computes only the upper-diagonal part!
      CALL ZLASET ('ALL', LIMX, LIMX, AU, BU, CK, LIMX)
      CALL ZHERK ('U', 'C', LIMX, LIMX, ONED, T, LIMX,
     +  -ONED, CK, LIMX)
      CALL ZHERK ('U', 'C', LIMX, LIMX, ONED, R, LIMX,
     +  ONED, CK, LIMX)
C$$$ ZLANGE FINDS MATRIX NORM
      DNORM1 = ZLANGE ('F', LIMX, LIMX, CK, LIMX)

c     Check that T~* T~ + R~* R~ = 1
      CALL ZLASET ('ALL', LIMX, LIMX, AU, BU, CK, LIMX)
      CALL ZHERK ('U', 'C', LIMX, LIMX, ONED, TTILDE, LIMX,
     +  -ONED, CK, LIMX)
      CALL ZHERK ('U', 'C', LIMX, LIMX, ONED, RTILDE, LIMX,
     +  ONED, CK, LIMX)
      DNORM2 = ZLANGE ('F', LIMX, LIMX, CK, LIMX)

c     Now check that T* R~ + R* T~ = 0
c     CALL ZLASET ('ALL', LIMX, LIMX, AU, AU, CK, LIMX)
      CALL ZGEMM ('C', 'N', LIMX, LIMX, LIMX, ONEC, T,
     +           LIMX, RTILDE, LIMX, ZEROC, CK, LIMX)
      CALL ZGEMM ('C', 'N', LIMX, LIMX, LIMX, ONEC, R,
     +           LIMX, TTILDE, LIMX, ONEC, CK, LIMX)
      DNORM3 = ZLANGE ('F', LIMX, LIMX, CK, LIMX)

C$$   CALL PRINTT (CHECK, 2 * LIMX, 'C2 ')
c     ZHERK does only the upper-triangular part. Therefore, the norm
c     DNORM1 should be roughly doubled, ditto DNORM2
c     There are two related off-diagonal blocks in the cross-product,
c     this doubles the norm DNORM3. (The relation is not exact,
c     since, the diagonal is not doubled)
      DNORM = 2.0* DNORM1 + 2.0 * DNORM2 + 2.0 * DNORM3
c     WRITE (*, *) 'CHECKUNI2:', DNORM1, DNORM2, DNORM3, DNORM
      CHECKUNI2 = DNORM
      RETURN
      END
C$$$  ROUTINE TO CHECK FOR UNITARY MATRICES

      DOUBLE PRECISION FUNCTION CHECKUNI3(LIMX, T,R,TTILDE,RTILDE)
      IMPLICIT NONE

      INTEGER LIMX, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX T(LIMX,LIMX),
     + R(LIMX,LIMX),TTILDE(LIMX,LIMX),RTILDE(LIMX,LIMX),
     + U(LIMX*2,LIMX*2), CHECK(LIMX*2,LIMX*2)


C$$$ TEST CASES OF T,R,T~ R~

C$$$ FILLS A UNIT MATRIX
      CALL ZLASET ('A', 2*LIMX, 2*LIMX, ZEROC, ONEC, CHECK, 2*LIMX)
      DO X=1, LIMX
         DO Y=1, LIMX
C$$$  TOP LEFT
            U(X,Y)=T(X,Y)
C$$$  BOTTOM LEFT
            U(X+LIMX,Y)=R(X,Y)
C$$$  TOP RIGHT
            U(X,Y+LIMX)=RTILDE(X,Y)
C$$$  BOTTOM RIGHT
            U(X+LIMX,Y+LIMX)=TTILDE(X,Y)
         END DO
      END DO

C$$$  ZHERK calculates U^H * U, and then subtracts 1.
C$$$  However, this is perfomed in the upper triangular part only
C$$   Therefore, we double the norm.
      CALL ZHERK('U', 'C', 2*LIMX, 2*LIMX, ONED, U, 2*LIMX,
     +     -ONED, CHECK, 2*LIMX)
C$$$  ZLANGE FINDS MATRIX NORM
      CHECKUNI3 = 2 * ZLANGE('F', 2*LIMX, 2*LIMX, CHECK, 2*LIMX)
      RETURN
      END
