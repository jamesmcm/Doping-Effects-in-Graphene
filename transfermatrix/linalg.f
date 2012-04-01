C
C     This function calculates dot product of two square matrices:
C
C         C = A . B
C      
C      N is the matrix size
C      
      SUBROUTINE SQDOT (C, A, B, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX A(N, N), B(N, N), C(N, N)
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
          
      CALL ZGEMM ('N', 'N',  N, N, N,
     &  ONEC,   A, N, B, N,
     &  ZEROC,  C, N)

      RETURN
      END
      
      SUBROUTINE SQDOTNH (C, A, B, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX A(N, N), B(N, N), C(N, N)
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
          
      CALL ZGEMM ('N', 'C',  N, N, N,
     &  ONEC,   A, N, B, N,
     &  ZEROC,  C, N)

      RETURN
      END
      
      SUBROUTINE SQDOTHN (C, A, B, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX A(N, N), B(N, N), C(N, N)
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
          
      CALL ZGEMM ('C', 'N',  N, N, N,
     &  ONEC,   A, N, B, N,
     &  ZEROC,  C, N)

      RETURN
      END
      
      SUBROUTINE SQUPDAXPY (Y, ALPHA, X, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX ALPHA
      DOUBLE COMPLEX X(N, N), Y(N, N)
          
      CALL ZAXPY (N*N, ALPHA, X, 1, Y, 1)
      
      RETURN
      END

      
C
C     This function calculates dot product of two square matrices,
C     and multiplies the result by a scalar alpha:     
C
C         C = alpha * A . B
C      
C      N is the matrix size
C      
      SUBROUTINE SQDOTAX (C, ALPHA, A, B, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX A(N, N), B(N, N), C(N, N)
      DOUBLE COMPLEX ALPHA, ZEROC/0.0/
     
      CALL ZGEMM ('N', 'N',  N, N, N,
     &  ALPHA,   A, N, B, N,
     &  ZEROC,   C, N)

      RETURN
      END

C
C     This function calculates dot product of two square matrices,
C     and updates C according to the formula
C
C         C := beta * C + alpha * A . B
C      
C      N is the matrix size
C            
      SUBROUTINE SQDOTUPD (BETA, C, ALPHA, A, B, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX A(N, N), B(N, N), C(N, N)
      DOUBLE COMPLEX BETA, ALPHA
     
      CALL ZGEMM ('N', 'N',  N, N, N,
     &  ALPHA,   A, N, B, N,
     &  BETA,  C, N)

      RETURN
      END
      
      SUBROUTINE INVERTMATRIX(MATRIX, N)
      IMPLICIT NONE
      INTEGER S
      INTEGER N, PIVOT(N, N)
      DOUBLE COMPLEX MATRIX(N, N), WORK(N*N)
      CALL ZGETRF(N, N, MATRIX, N, PIVOT, S)
      CALL ZGETRI(N, MATRIX, N, PIVOT, WORK, N*N, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'NON-INVERTABLE MATRIX WITH S=', S
         STOP
      END IF

      RETURN
      END

      DOUBLE PRECISION FUNCTION DNORMDIFF (A, B, N)
      DOUBLE COMPLEX A(N, N), B(N, N)
      DOUBLE COMPLEX D(N, N)
      DOUBLE PRECISION WORK (N)
      DOUBLE COMPLEX MINUSONE/-1.0/
      DOUBLE PRECISION ZLANGE
      
      CALL SQCOPY    (A, D, N)
      CALL SQUPDAXPY (D, MINUSONE, B, N)
      DNORMDIFF = ZLANGE ('F', N, N, D, N, WORK)

      RETURN
      END
      
      SUBROUTINE INVSOLVE(A, B, N)
c     Solve A := inv(A).B      
      IMPLICIT NONE
      INTEGER S
      INTEGER N, PIVOT(N, N)
      DOUBLE COMPLEX A(N, N), B(N, N), RHS(N, N)
c      DOUBLE COMPLEX SVA (N, N), PROD (N, N)
c      DOUBLE PRECISION dnormdiff
      
      CALL SQCOPY (B, RHS, N)
c     CALL SQCOPY (A, SVA, N)
      CALL ZGETRF(N, N, A, N, PIVOT, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'INVSOLVE: NON-INVERTABLE MATRIX WITH S=', S
         STOP
      END IF

      CALL ZGETRS('N', N, N, A, N, PIVOT, RHS, N, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'INVSOLVE: NON-SOLVABLE SYSTEM WITH S=', S
         STOP
      END IF
c      call sqdot (prod, SVA, rhs, N) 
c      write (*, *) 'quality: ', dnormdiff (B, PROD, N)
       
      CALL SQCOPY (RHS, A, N)
      
      RETURN
      END

      

      SUBROUTINE SQSVDVALS(N, MATRIX, OUTPUTS)
      INTEGER N
      DOUBLE COMPLEX MATRIX(N, N)
      DOUBLE PRECISION OUTPUTS(N)
      
      DOUBLE PRECISION RWORK(5*N)
      DOUBLE COMPLEX   TEMP2(N, N), SVCPY(N, N), WORK(4*N*N)
      INTEGER MSIZE, S

      MSIZE=4*N*N

C$$$ MAKE COPY OF MATRIX FOR SVD SINCE IT IS DESTROYED
C$$$      SVCPY=MATRIX
      CALL SQCOPY(MATRIX, SVCPY, N)
      CALL ZGESVD('N', 'N', N, N, SVCPY, N, OUTPUTS, TEMP2,
     + N, TEMP2, N , WORK, MSIZE, RWORK, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'SVD FAILED WITH S=', S
         STOP
      END IF
c      CALL DCOPY(N, SVALS, 1, OUTPUTS, 1)

      RETURN
      END
      
      SUBROUTINE SQSVDFULL(N, MATRIX, OUTPUTS, U, V)
      INTEGER N
      DOUBLE COMPLEX MATRIX(N, N),
     +               U(N, N), V(N, N)
      DOUBLE PRECISION OUTPUTS(N)
      
      INTEGER MSIZE, S
      DOUBLE PRECISION RWORK(5*N)
      DOUBLE COMPLEX   SVCPY(N, N),  WORK(4*N*N)

      MSIZE = 4*N*N

C$$$ MAKE COPY OF MATRIX FOR SVD SINCE IT IS DESTROYED
      CALL SQCOPY(MATRIX, SVCPY, N)
      
      CALL ZGESVD('N', 'N', N, N, SVCPY, N, OUTPUTS,
     +             U, N, V, N,
     +             WORK, MSIZE, RWORK, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'SVD FAILED WITH S=', S
         STOP

      END IF

c      CALL DCOPY(N, SVALS, 1, OUTPUTS, 1)

      RETURN
      END

      DOUBLE PRECISION FUNCTION CHECKUNI(N, T,R,TTILDE,RTILDE)
      IMPLICIT NONE

      INTEGER N, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
      DOUBLE COMPLEX T(N,N), BETA/-1/,ALPHA/1/,
     + R(N,N),TTILDE(N,N),RTILDE(N,N),
     + U(N*2,N*2), CHECK(N*2,N*2)


C$$$ TEST CASES OF T,R,T~ R~
C$$$ FILLS A UNIT MATRIX
      CALL SQUNIT (CHECK, 2*N)
      DO X=1, N
         DO Y=1, N
C$$$  TOP LEFT
            U(X,Y)=T(X,Y)
C$$$  BOTTOM LEFT
            U(X+N,Y)=R(X,Y)
C$$$  TOP RIGHT
            U(X,Y+N)=RTILDE(X,Y)
C$$$  BOTTOM RIGHT
            U(X+N,Y+N)=TTILDE(X,Y)
         END DO
      END DO

C$$$  ZGEMM HAS INBUILT FUNCTION TO FIND U**H
      CALL ZGEMM('N', 'C', 2*N, 2*N, 2*N, ALPHA, U,
     +     2*N, U, 2*N, BETA, CHECK, 2*N)
C$$$  ZLANGE FINDS MATRIX NORM
      CHECKUNI = ZLANGE('F', 2*N, 2*N, CHECK, 2*N)
      RETURN
      END

C$$$  ROUTINE TO CHECK FOR UNITARY MATRICES

      DOUBLE PRECISION FUNCTION CHECKUNI2(N, T,R,TTILDE,RTILDE)
      IMPLICIT NONE
 
      INTEGER N
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
      DOUBLE PRECISION DNORM1, DNORM2, DNORM3, DNORM
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX T(N,N), R(N,N),
     +               TTILDE(N,N),RTILDE(N,N),
     +               CK(N, N)


C$$$ TEST CASES OF T,R,T~ R~

C$$$ FILLS A UNIT MATRIX

C      DO X = 1, 2*N
C         DO Y = 1, 2*N
C            CHECK(X, Y) = (0.0, 0.0)
C         END DO
C       CHECK(X, X) = 1.0
C      END DO

C$$$ ZGEMM HAS INBUILT FUNCTION TO FIND U**H
C$$   CALL PRINTT (U, 2*N, 'U1 ')

c     Unitarity can be also checked using subblocks only
c     The definition of unitarity can be easily cast into the form
c     |T|^2 + |R|^2 = 1, T*R~ + R* T~ = 0

c     Check that T * T + R* R = 1   (* is the Hermitian conjugate)
c     Here I use zherk, which does T^+T
c     Beware: it computes only the upper-diagonal part!
      CALL SQUNIT (CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, T, N,
     +  -ONED, CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, R, N,
     +  ONED, CK, N)
C$$$ ZLANGE FINDS MATRIX NORM
      DNORM1 = ZLANGE ('F', N, N, CK, N)

c     Check that T~* T~ + R~* R~ = 1
      CALL SQUNIT (CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, TTILDE, N,
     +  -ONED, CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, RTILDE, N,
     +  ONED, CK, N)
      DNORM2 = ZLANGE ('F', N, N, CK, N)

c     Now check that T* R~ + R* T~ = 0
      CALL ZGEMM ('C', 'N', N, N, N, ONEC, T,
     +           N, RTILDE, N, ZEROC, CK, N)
      CALL ZGEMM ('C', 'N', N, N, N, ONEC, R,
     +           N, TTILDE, N, ONEC, CK, N)
      DNORM3 = ZLANGE ('F', N, N, CK, N)

C$$   CALL PRINTT (CHECK, 2 * N, 'C2 ')
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

      DOUBLE PRECISION FUNCTION CHECKUNI3(N, T,R,TTILDE,RTILDE)
      IMPLICIT NONE

      INTEGER N, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX T(N,N),
     + R(N,N),TTILDE(N,N),RTILDE(N,N),
     + U(N*2,N*2), CHECK(N*2,N*2)


C$$$ TEST CASES OF T,R,T~ R~

C$$$ FILLS A UNIT MATRIX
      CALL SQUNIT (CHECK, 2*N)
      DO X=1, N
         DO Y=1, N
C$$$  TOP LEFT
            U(X,Y)=T(X,Y)
C$$$  BOTTOM LEFT
            U(X+N,Y)=R(X,Y)
C$$$  TOP RIGHT
            U(X,Y+N)=RTILDE(X,Y)
C$$$  BOTTOM RIGHT
            U(X+N,Y+N)=TTILDE(X,Y)
         END DO
      END DO

C$$$  ZHERK calculates U^H * U, and then subtracts 1.
C$$$  However, this is perfomed in the upper triangular part only
C$$   Therefore, we double the norm.
      CALL ZHERK('U', 'C', 2*N, 2*N, ONED, U, 2*N,
     +     -ONED, CHECK, 2*N)
C$$$  ZLANGE FINDS MATRIX NORM
      CHECKUNI3 = 2 * ZLANGE('F', 2*N, 2*N, CHECK, 2*N)
      RETURN
      END
