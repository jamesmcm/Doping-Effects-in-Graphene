      PROGRAM TRANSFERM
      INTEGER, PARAMETER :: LIMX=2, LIMY=2, WRAPY=0, WRAPX=1,
     +     MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER PIVOT(2*LIMX, 2*LIMX), PIVOT2(LIMX, LIMX)
      INTEGER*4 I/1/, J/1/, E/2/, S/9/, K/1/
      DOUBLE PRECISION SVALS(LIMX), RWORK(5*LIMX)
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX), 
     +     MULT(2*LIMX, 2*LIMX), OUT(2*LIMX, 2*LIMX), ALPHA, BETA,
     +     O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX),
     +     TEMP(2*LIMX, 2*LIMX), ABCD(2*LIMX, 2*LIMX), A(LIMX, LIMX),
     +     B(LIMX, LIMX), C(LIMX, LIMX), D(LIMX, LIMX),
     +     WORK(MSIZE), TEMP2(LIMX, LIMX), T(LIMX, LIMX),
     +     TTILDE(LIMX, LIMX), R(LIMX, LIMX), RTILDE(LIMX, LIMX),
     +     WORK2(LIMX*LIMX), TEMP3(LIMX, LIMX), SVCPY(LIMX, LIMX)

      DATA MODD/MSIZE*0.0/, MEVEN/MSIZE*0.0/, O/MSIZE*0.0/,
     +     IO/MSIZE*0.0/, PIVOT/MSIZE*0/, ALPHA/1.0/,
     +     BETA/0.0/, TEMP/MSIZE*0.0/, PIVOT2/M2SIZE*0.0/


c$$$  First row is even - WRAPX makes no diff, second row not, etc.
c$$$  - what matters is which row it is centred on
c$$$  There are 2 transfer matrices to generate
c$$$  There are 4 block submatrices to fill
c$$$  MODD doesn't depend on xwrapping, MEVEN does.
      DO I = 1, LIMX
c$$$  Fill top-right submatrix
         MODD(I, LIMX+I)=1
         MEVEN(I, LIMX+I)=1
c$$$  Fill bottom-left submatrix
         MODD(I+LIMX, I)=-1
         MEVEN(I+LIMX, I)=-1
c$$$  Fill bottom-right  submatrix
         MODD(LIMX+I,LIMX+I)=E
         MEVEN(LIMX+I,LIMX+I)=E
c$$$  Smart trick to set alternating adjacent value to -1
         MODD(LIMX+I,LIMX+I+(2*MOD(I,2)-1))=-1
c$$$         Trick to set end values dependent on WRAPX
c$$$  Only the ends matter with regards to the  WRAPX effect
         MEVEN(LIMX+I, 2*LIMX+1-I)=-1*((-1*(I/(I*I)))+1)*((-1*(((LIMX-I)
     +        +1)/(((LIMX-I)+1)*((LIMX-I)+1))))+1) + (I/(I*I))*
     +        (-1*WRAPX) +(((LIMX-I)+1)/(((LIMX-I)+1)*((LIMX-I)+1)))
     +        *(-1*WRAPX)
      END DO	
      
      IF (MOD(LIMY,2) .EQ. 1) THEN
         MULT = MODD
      ELSE
         MULT = MEVEN
      END IF
         
      DO I = 1, LIMY-1
         IF (MOD(LIMY,2) .EQ. 1) THEN
            IF (MOD(I,2) .EQ. 1) THEN
c$$$               MULT = MATMUL(MULT,MEVEN)
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MEVEN, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            ELSE
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MODD, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            END IF
         ELSE
            IF (MOD(I,2) .EQ. 1) THEN
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MODD, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            ELSE
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MEVEN, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            END IF
         END IF
      END DO      

c$$$  MULT now holds the net transfer matrix for moving down vertically
c$$$      DO J = 1, 2*LIMX
c$$$         WRITE (*,20) (REAL(MULT(J,I)), I = 1, 2*LIMX)
c$$$      END DO

c$$$  Generate O-matrix	
c$$$  O is block matrix of 1/sqrt(2) (1,1;i,-i)
      DO I = 1, LIMX
         O(I, I)=1/SQRT(2.0)
         O(I, LIMX+I)=1/SQRT(2.0)
         O(I+LIMX, I)=DCMPLX(0, 1/SQRT(2.0))
         O(I+LIMX, I+LIMX)=DCMPLX(0 ,-1/SQRT(2.0))
      ENDDO

      IO=O
      CALL ZGETRF(2*LIMX, 2*LIMX, IO, 2*LIMX, PIVOT, S)
      CALL ZGETRI(2*LIMX, IO, 2*LIMX, PIVOT, WORK, 4*LIMX*LIMX, S)
      CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +     2*LIMX, O, 2*LIMX, BETA, TEMP, 2*LIMX)
      CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, IO, 
     +     2*LIMX, TEMP, 2*LIMX, BETA, ABCD, 2*LIMX)
c$$$
c$$$      PRINT *, 'ABCD matrix:'
c$$$
c$$$      DO J = 1, 2*LIMX
c$$$         DO I=1, 2*LIMX
c$$$            WRITE (*,30) REAL(ABCD(J,I)), ' + ', DIMAG(ABCD(J,I)), 'I'
c$$$         END DO
c$$$         PRINT *, '----'
c$$$      END DO
c$$$  This is Fortran 90 syntax, remove in future revision when BLAS/LAPACK subroutine is found
      A=ABCD(1:LIMX, 1:LIMX)
      B=ABCD(LIMX+1:2*LIMX, 1:LIMX)
      C=ABCD(1:LIMX, LIMX+1:2*LIMX)
      D=ABCD(LIMX+1:2*LIMX, LIMX+1:2*LIMX)
c$$$  I have verified that AD-BC=1 (identity matrix) as expected
c$$$  T~ = D^-1
      TTILDE=D
      CALL ZGETRF(LIMX, LIMX, TTILDE, LIMX, PIVOT2, S)
      CALL ZGETRI(LIMX, TTILDE, LIMX, PIVOT2, WORK2, LIMX*LIMX, S)
c$$$  R~ = BD^-1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, ALPHA, B, 
     +     LIMX, TTILDE, LIMX, BETA, RTILDE, LIMX)      
c$$$  R = -D^-1 C
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, ALPHA, TTILDE, 
     +     LIMX, C, LIMX, BETA, R, LIMX)      
c$$$  R=-R
c$$$  T=(A-)? BD^-1 C      
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, ALPHA, TTILDE, 
     +     LIMX, C, LIMX, BETA, TEMP2, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, ALPHA, B, 
     +     LIMX, TEMP2, LIMX, BETA, T, LIMX)      
      T=A-T
c$$$  make copy of matrix for SVD since it is destroyed
      SVCPY=T
      CALL ZGESVD('N', 'N', LIMX, LIMX, SVCPY, LIMX, SVALS, TEMP2, LIMX,
     +     TEMP2, LIMX , WORK, MSIZE, RWORK, S)

      PRINT *, 'T matrix SVD values squared:'
      DO J=1, LIMX
         WRITE (*,40) (SVALS(J)*SVALS(J))
      END DO

      SVCPY=R
      CALL ZGESVD('N', 'N', LIMX, LIMX, SVCPY, LIMX, SVALS, TEMP2, LIMX,
     +     TEMP2, LIMX , WORK, MSIZE, RWORK, S)
      print *,S
      PRINT *, 'R matrix SVD values squared:'
      DO J=1, LIMX
         WRITE (*,40) (SVALS(J)*SVALS(J))
      END DO
         
c$$$  so T^2 + R^2 =1 for SVD values, also verified with R~ and T~

 20   FORMAT (4F4.0)
 30   FORMAT (F8.4, A, F8.4, A)
 40   FORMAT (F8.4)

      END
