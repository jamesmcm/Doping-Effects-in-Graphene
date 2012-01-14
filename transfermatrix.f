      PROGRAM TRANSFERM
      INTEGER, PARAMETER :: LIMX=2, LIMY=2, WRAPY=0, WRAPX=1,
     +     MSIZE=4*LIMX*LIMX
      INTEGER PIVOT(2*LIMX, 2*LIMX)
      INTEGER*4 I/1/, J/1/, E/2/, S/9/, K/1/

      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX), 
     +     MULT(2*LIMX, 2*LIMX), OUT(2*LIMX, 2*LIMX), ALPHA, BETA
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX WORK(4*LIMX*LIMX)
      DATA MODD/MSIZE*0.0/, MEVEN/MSIZE*0.0/, O/MSIZE*0.0/,
     +     IO/MSIZE*0.0/, PIVOT/MSIZE*0/, ALPHA/1.0/,
     +     BETA/0.0/


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
c$$$               MULT = MATMUL(MULT,MODD)
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MODD, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            END IF
         ELSE
            IF (MOD(I,2) .EQ. 1) THEN
c$$$               MULT = MATMUL(MULT,MODD)
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MODD, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            ELSE
c$$$               MULT = MATMUL(MULT, MEVEN)
               CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +              2*LIMX, MEVEN, 2*LIMX, BETA, OUT, 2*LIMX)
               MULT=OUT
            END IF
         END IF
      END DO      

c$$$      CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MODD, 
c$$$     +     2*LIMX, MEVEN, 2*LIMX, BETA, OUT, 2*LIMX)

      DO J = 1, 2*LIMX
         WRITE (*,20) (REAL(OUT(J,I)), I = 1, 2*LIMX)
      END DO

c$$$  Generate O-matrix	
c$$$  O is block matrix of 1/sqrt(2) (1,1;i,-i)
c$$$      DO I = 1, LIMX
c$$$         O(I, I)=1/SQRT(2.0)
c$$$         O(I, LIMX+I)=1/SQRT(2.0)
c$$$         O(I+LIMX, I)=DCMPLX(0, 1/SQRT(2.0))
c$$$         O(I+LIMX, I+LIMX)=DCMPLX(0 ,-1/SQRT(2.0))
c$$$      ENDDO
c$$$
c$$$      IO=O
c$$$      CALL ZGETRF(2*LIMX, 2*LIMX, IO, 2*LIMX, PIVOT, S)
c$$$      CALL ZGETRI(2*LIMX, IO, 2*LIMX, PIVOT, WORK, 4*LIMX*LIMX, S)
c$$$      IF (S .EQ. 0) THEN
c$$$         PRINT *, 'GREAT SUCCESS'
c$$$      ELSE
c$$$         PRINT *, 'TERRIBLE FAILURE'
c$$$      END IF
c$$$
c$$$      PRINT *, 'O matrix:'
c$$$
c$$$      DO J = 1, 2*LIMX
c$$$         DO I=1, 2*LIMX
c$$$            WRITE (*,30) REAL(O(J,I)), ' + ', DIMAG(O(J,I)), 'I'
c$$$         END DO
c$$$         PRINT *, '----'
c$$$      END DO
c$$$
c$$$      PRINT *, 'IO matrix:'
c$$$
c$$$      DO J = 1, 2*LIMX
c$$$         DO I=1, 2*LIMX
c$$$            WRITE (*,30) REAL(IO(J,I)), ' + ', DIMAG(IO(J,I)), 'I'
c$$$         END DO
c$$$         PRINT *, '----'
c$$$      END DO

 20   FORMAT (4F4.0)
 30   FORMAT (F6.4, A, F6.4, A)


      END
