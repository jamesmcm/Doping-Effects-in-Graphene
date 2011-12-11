      PROGRAM TRANSFERM
      INTEGER, PARAMETER :: LIMX=2, LIMY=2, WRAPY=0, WRAPX=1,
     +     MSIZE=4*LIMX*LIMX
      INTEGER*4 I/1/, J/1/, E/2/

      REAL MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX), 
     +     MULT(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX)
      DATA MODD/MSIZE*0.0/, MEVEN/MSIZE*0.0/, O/MSIZE*0.0/


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
               MULT = MATMUL(MULT,MEVEN)
            ELSE
               MULT = MATMUL(MULT,MODD)
            END IF
         ELSE
            IF (MOD(I,2) .EQ. 1) THEN
               MULT = MATMUL(MULT,MODD)
            ELSE
               MULT = MATMUL(MULT, MEVEN)
            END IF
         END IF
      END DO      

c$$$      DO J = 1, 2*LIMX
c$$$         WRITE (*,20) (MULT(J,I), I = 1, 2*LIMX)
c$$$      END DO

c$$$  Generate O-matrix	
c$$$  O is block matrix of 1/sqrt(2) (1,1;i,-i)
      DO I = 1, LIMX
         O(I, I)=1/SQRT(2.0)
         O(I, LIMX+I)=1/SQRT(2.0)
         O(I+LIMX, I)=DCMPLX(0, 1/SQRT(2.0))
         O(I+LIMX, I+LIMX)=DCMPLX(0 ,-1/SQRT(2.0))
      ENDDO

      DO J = 1, 2*LIMX
         DO I=1, 2*LIMX
            WRITE (*,30) REAL(O(J,I)), ' + ', DIMAG(O(J,I)), 'I'
         END DO
      END DO

 20   FORMAT (8F3.0)
 30   FORMAT (F6.4, A, F6.4, A)

      END
