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
      
      CALL ZGESVD('A', 'A', N, N, SVCPY, N, OUTPUTS,
     +             U, N, V, N,
     +             WORK, MSIZE, RWORK, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'SVD FAILED WITH S=', S
         STOP

      END IF

c      CALL DCOPY(N, SVALS, 1, OUTPUTS, 1)

      RETURN
      END

