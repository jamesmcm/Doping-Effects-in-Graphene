c$$$  Want to loop over different energies and produce T^2 coefficients, check they match with analytical results

      PROGRAM TRANSFERM
      INTEGER, PARAMETER :: LIMX=2,WRAPY=0, WRAPX=1,
     +     MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER PIVOT(2*LIMX, 2*LIMX), PIVOT2(LIMX, LIMX)
      INTEGER*4 I/1/, J/1/, S/9/, K/1/, F/1/, LIMY/10/
      CHARACTER*3 VALUE

      DOUBLE PRECISION SVALS(LIMX), RWORK(5*LIMX), RVALS(LIMX),
     +     TVALS(LIMX), TTVALS(LIMX), RTVALS(LIMX), E/-5/
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


c$$$  reads command line argument as LIMY
      CALL GETARG(1, VALUE)
      READ(UNIT=VALUE, FMT=*) LIMY

      DO F = 1, 1001
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
c$$$  Trick to set end values dependent on WRAPX
c$$$  Only the ends matter with regards to the  WRAPX effect
            MEVEN(LIMX+I, 2*LIMX+1-I)=-1*((-1*(I/(I*I)))+1)*((-1*(((LIMX
     +           -I)+1)/(((LIMX-I)+1)*((LIMX-I)+1))))+1) + (I/(I*I))*
     +           (-1*WRAPX) +(((LIMX-I)+1)/(((LIMX-I)+1)*((LIMX-I)+1)))
     +           *(-1*WRAPX)
         END DO	
         
         IF (MOD(LIMY,2) .EQ. 1) THEN
            MULT = MODD
         ELSE
            MULT = MEVEN
         END IF
         
         DO I = 1, LIMY-1
            IF (MOD(LIMY,2) .EQ. 1) THEN
               IF (MOD(I,2) .EQ. 1) THEN
                  CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, 
     +                  MULT, 2*LIMX, MEVEN, 2*LIMX, BETA, OUT, 2*LIMX)
                  MULT=OUT
               ELSE
                  CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA,  
     +                 MULT,2*LIMX, MODD, 2*LIMX, BETA, OUT, 2*LIMX)
                  MULT=OUT
               END IF
            ELSE
               IF (MOD(I,2) .EQ. 1) THEN
                  CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA,  
     +                 MULT,2*LIMX, MODD, 2*LIMX, BETA, OUT, 2*LIMX)
                  MULT=OUT
               ELSE
                  CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA,  
     +                 MULT,2*LIMX, MEVEN, 2*LIMX, BETA, OUT, 2*LIMX)
                  MULT=OUT
               END IF
            END IF
         END DO      
         
c$$$  MULT now holds the net transfer matrix for moving down vertically
c$$$  Generate O-matrix	
c$$$  O is block matrix of 1/sqrt(2) (1,1;i,-i)
         DO I = 1, LIMX
            O(I, I)=1/SQRT(2.0)
            O(I, LIMX+I)=1/SQRT(2.0)
            O(I+LIMX, I)=DCMPLX(0, 1/SQRT(2.0))
            O(I+LIMX, I+LIMX)=DCMPLX(0 ,-1/SQRT(2.0))
         ENDDO
         
         IO=O
         CALL INVERSAS(IO, 2*LIMX)
         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, MULT, 
     +        2*LIMX, O, 2*LIMX, BETA, TEMP, 2*LIMX)
         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, IO, 
     +        2*LIMX, TEMP, 2*LIMX, BETA, ABCD, 2*LIMX)

         CALL AKIRAS(ABCD, LIMX, TVALS, RVALS, TTVALS, RTVALS)

         WRITE (*,10) E,(TVALS(I)*TVALS(I), I = 1, LIMX)

         E=E+0.01
      END DO

 10   FORMAT (3F10.6)
      END

      SUBROUTINE INVERSAS(X, LDM)
      DOUBLE COMPLEX X(LDM, LDM), WORK(4*LDM*LDM)
      INTEGER LDM, PIVOT(LDM, LDM), S

      CALL ZGETRF(LDM, LDM, X, LDM, PIVOT, S)
      CALL ZGETRI(LDM, X, LDM, PIVOT, WORK, LDM*LDM, S)
      END

      SUBROUTINE AKIRAS(ABCD, LIMX, TVALS, RVALS, TTVALS, RTVALS)
      DOUBLE PRECISION RWORK(5*LIMX), RVALS(LIMX), TVALS(LIMX),
     +     RTVALS(LIMX), TTVALS(LIMX)
      DOUBLE COMPLEX ALPHA, BETA, ABCD(2*LIMX, 2*LIMX), A(LIMX, LIMX),
     +     B(LIMX, LIMX), C(LIMX, LIMX), D(LIMX, LIMX),
     +     WORK(4*LIMX*LIMX), TEMP2(LIMX, LIMX), T(LIMX, LIMX),
     +     TTILDE(LIMX, LIMX), R(LIMX, LIMX), RTILDE(LIMX, LIMX),
     +     WORK2(LIMX*LIMX)
      DATA  ALPHA/1.0/, BETA/0.0/

c$$$  This is Fortran 90 syntax, remove in future revision when BLAS/LAPACK subroutine is found
      A=ABCD(1:LIMX, 1:LIMX)
      B=ABCD(LIMX+1:2*LIMX, 1:LIMX)
      C=ABCD(1:LIMX, LIMX+1:2*LIMX)
      D=ABCD(LIMX+1:2*LIMX, LIMX+1:2*LIMX)
c$$$  I have verified that AD-BC=1 (identity matrix) as expected
c$$$  T~ = D^-1
      TTILDE=D
      CALL INVERSAS(TTILDE, LIMX)
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
      CALL ZGESVD('N', 'N', LIMX, LIMX, T, LIMX, TVALS, TEMP2, 
     +     LIMX,TEMP2, LIMX , WORK, 4*LIMX*LIMX, RWORK, S)

      CALL ZGESVD('N', 'N', LIMX, LIMX, R, LIMX, RVALS, TEMP2, 
     +     LIMX,TEMP2, LIMX , WORK, 4*LIMX*LIMX, RWORK, S)
      CALL ZGESVD('N', 'N', LIMX, LIMX, RTILDE, LIMX, RTVALS, TEMP2, 
     +     LIMX,TEMP2, LIMX , WORK, 4*LIMX*LIMX, RWORK, S)
      CALL ZGESVD('N', 'N', LIMX, LIMX, TTILDE, LIMX, TTVALS, TEMP2, 
     +     LIMX,TEMP2, LIMX , WORK, 4*LIMX*LIMX, RWORK, S)
      END

      SUBROUTINE UPDVALS(NEWTS, CURTS, NEWRTS, CURRTS, NEWRS, 
     +     CURRS, CURTTS, NEWTTS, LIMX)
c$$$  Updates the values of the SVD values for the numerically stable algorithm
      INTEGER LIMX, G
      DOUBLE PRECISION NEWTS(LIMX), CURTS(LIMX), CURRTS(LIMX),
     +     NEWRS(LIMX), NEWRTS(LIMX), CURRS(LIMX), CURTTS(LIMX),
     +     NEWTTS(LIMX), TEMPT(LIMX), TEMPR(LIMX), TEMPTT(LIMX),
     +     TEMPRT(LIMX)

      DO G = 1, LIMX
         TEMPT(G)=NEWTS(G)*(1-(CURRTS(G)*NEWRS(G)))*CURTS(G)
         TEMPTT(G) = CURTTS(G) * (1 - (NEWRS(G) * CURRTS(G)))
     +        * NEWTTS(G)
         TEMPR(G)=CURRS(G) + ((CURTTS(G)*NEWRS(G)*CURTS(G))/
     +        (1-(NEWRS(G)*CURRTS(G))))
         TEMPRT(G)=NEWRTS(G) + ((NEWTS(G)*CURRTS(G)*NEWTTS(G))/
     +        (1-(CURRTS(G)*NEWRS(G))))
      END DO

      CURTS=TEMPT
      CURRS=TEMPR
      CURTTS=TEMPTT
      CURRTS=TEMPRT

      END

