c$$$ Want to loop over different energies and produce T^2 coefficients, check they match with analytical results

      PROGRAM TRANSFERMATIXTWO
      INTEGER, PARAMETER :: LIMX=2, WRAPY=0, WRAPX=1,
     + MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER PIVOT(2*LIMX, 2*LIMX), PIVOT2(LIMX, LIMX)
      INTEGER I/1/, J/1/, S/9/, K/1/, F/1/, X/1/, Y/1/,LIMY/10/       
      CHARACTER*3 VALUE
      DOUBLE PRECISION SVALS(LIMX), RWORK(5*LIMX), RVALS(LIMX),
     + TVALS(LIMX), E/-5/, TTVALS(LIMX), RTVALS(LIMX),SQUARE,COND/-1/, 
     + CHECKUNI
      
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX),
     + MULT(2*LIMX, 2*LIMX), OUTM(2*LIMX, 2*LIMX), ALPHA, BETA,
     + O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX),
     + TEMP(2*LIMX, 2*LIMX), ABCD(2*LIMX, 2*LIMX), A(LIMX, LIMX),
     + B(LIMX, LIMX), C(LIMX, LIMX), D(LIMX, LIMX),
     + WORK(MSIZE), TEMP2(LIMX, LIMX), T(LIMX, LIMX),
     + TTILDE(LIMX, LIMX), R(LIMX, LIMX), RTILDE(LIMX, LIMX),
     + TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX), RINC(LIMX, LIMX),
     + RTILDEINC(LIMX, LIMX), WORK2(LIMX*LIMX), TEMP3(LIMX, LIMX),
     + SVCPY(LIMX, LIMX)

      DATA MODD/MSIZE*0.0/, MEVEN/MSIZE*0.0/, O/MSIZE*0.0/,
     + IO/MSIZE*0.0/, PIVOT/MSIZE*0/, ALPHA/1.0/,
     + BETA/0.0/, TEMP/MSIZE*0.0/, PIVOT2/M2SIZE*0.0/

c$$$  reads command line argument as LIMY
      CALL GETARG(1, VALUE)
      READ(UNIT=VALUE, FMT=*) LIMY

      DO F = 1, 1001


         CALL CALCMULT(MULT, LIMX, LIMY, WRAPX, MODD, MEVEN, E)

         CALL FILLOANDINVERT(O, IO, LIMX)
         CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
         CALL GENTANDRINC(LIMX, T, R, TTILDE, RTILDE, A, B, C, D) 
	

		 
c$$$ ################################################################


         DO I = 1, LIMY-1
            IF (MOD(LIMY,2) .EQ. 1) THEN
               IF (MOD(I,2) .EQ. 1) THEN
                  CALL ZCOPY(4*LIMX*LIMX, MEVEN, 1, MULT, 1)		
               ELSE
                  CALL ZCOPY(4*LIMX*LIMX, MODD, 1, MULT, 1)		
               END IF
            ELSE
               IF (MOD(I,2) .EQ. 1) THEN
                  CALL ZCOPY(4*LIMX*LIMX, MODD, 1, MULT, 1)		
               ELSE
                  CALL ZCOPY(4*LIMX*LIMX, MEVEN, 1, MULT, 1)		
               END IF
            END IF


			
            CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
            CALL GENTANDRINC(LIMX, TINC, RINC, TTILDEINC, RTILDEINC, 
     +       A, B, C,D)
            CALL UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     +       RTILDE, LIMX, RINC)
	 
			
         END DO

c$$$ ################################################################



		 
         CALL SV_DECOMP(LIMX, T, TVALS)
         CALL SV_DECOMP(LIMX, R, RVALS)
         CALL SV_DECOMP(LIMX, TTILDE, TTVALS)
         CALL SV_DECOMP(LIMX, RTILDE, RTVALS)

c$$$ Output of main function-commented out while checking for unitarity
         		 

c$$$         CALL PRINTVECTOR(TVALS, LIMX, 'T ')
c$$$         CALL PRINTVECTOR(RVALS, LIMX, 'R ')
c$$$         CALL PRINTVECTOR(TTVALS, LIMX, 'T~')
c$$$         CALL PRINTVECTOR(RTVALS, LIMX, 'R~')

	 
          
c      DO X = 1, LIMX*LIMX
c         DO Y = 1, LIMX*LIMX
c            T(X, Y) = (0.0, 0.0)
c            R(X, Y) = (0.0, 0.0)
c           TTILDE(X, Y) = (0.0, 0.0)
c            RTILDE(X, Y) = (0.0, 0.0)
c         END DO
c        T(X, X) = 1.0
c        TTILDE(X, Y) = 1.0
c      END DO
    
c$$$  Error return type mismatch of function checkuni real(4)/real(8)      
      COND = CHECKUNI(LIMX,T,R,TTILDE,RTILDE)
     
      WRITE(*,50) E,COND,(TVALS(I)*TVALS(I), I = 1, LIMX)
         E=E+0.01
      END DO
c$$$ so T^2 + R^2 =1 for SVD values, also verified with R~ and T~

      
 50   FORMAT (F6.2,3ES15.5E2)      

      

      STOP
      END

c$$$ End of main program ###############################################
c$$$                     ###############################################
 
c$$$ Routine to generate a,b,c,d. 
      SUBROUTINE GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
	 
      INTEGER LIMX
      DOUBLE COMPLEX UNITY, ZERO, A(LIMX, LIMX), B(LIMX, LIMX),
     + C(LIMX, LIMX), D(LIMX, LIMX), MULT(2*LIMX, 2*LIMX),
     + O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX), TEMP(2*LIMX, 2*LIMX),
     + ABCD(2*LIMX, 2*LIMX)
	 
      UNITY = 1.0
      ZERO = 0.0	 
	 
c$$$ zgemm Multiplies matrices togther	  

         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, UNITY, MULT,
     +         2*LIMX, O, 2*LIMX, ZERO, TEMP, 2*LIMX)
         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, UNITY, IO,
     +         2*LIMX, TEMP, 2*LIMX, ZERO, ABCD, 2*LIMX)	  

c$$$ This is Fortran 90 syntax, remove in future revision when BLAS/LAPACK subroutine is found
         A=ABCD(1:LIMX, 1:LIMX)
         B=ABCD((LIMX+1):2*LIMX, 1:LIMX)
         C=ABCD(1:LIMX, LIMX+1:2*LIMX)
         D=ABCD((LIMX+1):2*LIMX, (LIMX+1):2*LIMX)

c$$$ I have verified that AD-BC=1 (identity matrix) as expected
      RETURN
      END
	  
c$$$ Routine to generate T and R f	  
      SUBROUTINE GENTANDRINC(LIMX, TINC,RINC,TTILDEINC,RTILDEINC,A,B,
     + C,D)
	 
      INTEGER LIMX, X, Y
      DOUBLE COMPLEX UNITY, ZERO, A(LIMX, LIMX),
     + C(LIMX, LIMX), D(LIMX, LIMX), TEMP2(LIMX, LIMX),
     + TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX), RINC(LIMX, LIMX),
     + RTILDEINC(LIMX, LIMX), B(LIMX, LIMX), UNITMATRIX(LIMX, LIMX)

      UNITY = 1.0
      ZERO = 0.0
      DO X = 1, LIMX
         DO Y=1, LIMX
            UNITMATRIX(X, Y) = (0.0, 0.0)
         END DO
         UNITMATRIX(X, X) = 1.0
      END DO
	 
c$$$ T~ = D^-1
      CALL ZCOPY(LIMX*LIMX, D, 1, TTILDEINC, 1)
      CALL INVERTMATRIX(TTILDEINC, LIMX)
c$$$ R~ = BD^-1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, B,
     +      LIMX, TTILDEINC, LIMX, ZERO, RTILDEINC, LIMX)
c$$$ R = -D^-1 C
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, -1*UNITY, TTILDEINC,
     +      LIMX, C, LIMX, ZERO, RINC, LIMX)
c$$$ R=-R
c$$$      RINC = -1*RINC   Produces non unitary matrix ? 
c$$$ T=(A-)? BD^-1 C
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDEINC,
     +      LIMX, C, LIMX, ZERO, TEMP2, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, B,
     +      LIMX, TEMP2, LIMX, ZERO, TINC, LIMX)
c$$$         TINC=A-TINC
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, A, LIMX,
     +      UNITMATRIX, LIMX, -1*UNITY, TINC, LIMX)

      RETURN
      END
	  
c$$$ With increasing energy updates T and R to new values.

      SUBROUTINE UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     + RTILDE, LIMX, RINC)
	  
      INTEGER LIMX
      DOUBLE COMPLEX T(LIMX, LIMX), TTILDE(LIMX, LIMX), R(LIMX, LIMX),
     + RTILDE(LIMX, LIMX), TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     + RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)

c$$$ Temporary local variables
      DOUBLE COMPLEX T1TEMP(LIMX, LIMX), TTILDE1TEMP(LIMX, LIMX),
     + R1TEMP(LIMX, LIMX), RTILDE1TEMP(LIMX, LIMX), T2TEMP(LIMX, LIMX),
     + TTILDE2TEMP(LIMX, LIMX),R2TEMP(LIMX, LIMX), UNITY, ZERO,
     + RTILDE2TEMP(LIMX, LIMX), BRACKET12(LIMX, LIMX),
     + BRACKET21(LIMX, LIMX), TRTEMP(LIMX, LIMX), TRTEMP2(LIMX, LIMX),
     + UNITMATRIX(LIMX, LIMX), ALLZERO(LIMX, LIMX)
      INTEGER X, Y
	  
      UNITY = 1.0
      ZERO = 0.0
	 
c$$$  Is it copy-by-reference, or by-value? 
c$$$  Use lapack functions, such as zcopy --- AVS
      CALL ZCOPY(LIMX*LIMX, T, 1, T1TEMP, 1)
      CALL ZCOPY(LIMX*LIMX, TTILDE, 1, TTILDE1TEMP, 1)
      CALL ZCOPY(LIMX*LIMX, R, 1, R1TEMP, 1)	  
      CALL ZCOPY(LIMX*LIMX, RTILDE, 1, RTILDE1TEMP, 1)	  
      CALL ZCOPY(LIMX*LIMX, TINC, 1, T2TEMP, 1)	  
      CALL ZCOPY(LIMX*LIMX, TTILDEINC, 1, TTILDE2TEMP, 1)	  
      CALL ZCOPY(LIMX*LIMX, RINC, 1, R2TEMP, 1)	  
      CALL ZCOPY(LIMX*LIMX, RTILDEINC, 1, RTILDE2TEMP, 1)
	  
      DO X = 1, LIMX
         DO Y = 1, LIMX
            ALLZERO(X, Y) = (0.0, 0.0)
            UNITMATRIX(X, Y) = (0.0, 0.0)
         END DO
         UNITMATRIX(X, X) = 1.0
      END DO
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP2, 1)

c$$$ BRACKET12 = (1 - RTILDE1.R2)^-1
c$$$      CALL PRINTT(RTILDE1TEMP, LIMX, 'Rt1')
c$$$      CALL PRINTT(R2TEMP, LIMX, 'R2 ')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, RTILDE1TEMP, LIMX,
     + R2TEMP, LIMX, ZERO, BRACKET12, LIMX)
c$$$      CALL PRINTT(BRACKET12, LIMX, 'b12')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + UNITMATRIX, LIMX, -1*UNITY, BRACKET12, LIMX)
c$$$      CALL PRINTT(BRACKET12, LIMX, 'b12')
      CALL INVERTMATRIX(BRACKET12, LIMX)
c$$$      CALL PRINTT(BRACKET12, LIMX, 'B12')
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)


c$$$ BRACKET21 = (1 - R2.RTILDE1)^-1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, R2TEMP, LIMX,
     + RTILDE1TEMP, LIMX, ZERO, BRACKET21, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + UNITMATRIX, LIMX, -1*UNITY, BRACKET21, LIMX)
      CALL INVERTMATRIX(BRACKET21, LIMX)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)


c$$$ T = T2.BRACKET12.T1
c$$$      CALL PRINTT(T1TEMP, LIMX, 'T1t')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET12, LIMX,
     + T1TEMP, LIMX, ZERO, TRTEMP, LIMX)
c$$$      CALL PRINTT(T2TEMP, LIMX, 'T2t')
c$$$      CALL PRINTT(TRTEMP, LIMX, 'TRt')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, T2TEMP, LIMX,
     + TRTEMP, LIMX, ZERO, T, LIMX)
c$$$      CALL PRINTT(T, LIMX, 'T  ')
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)


c$$$ TTILDE = TTILDE1.BRACKET21.TTILDE2
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET21, LIMX,
     + TTILDE2TEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDE1TEMP, LIMX,
     + TRTEMP, LIMX, ZERO, TTILDE, LIMX)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)


c$$$ R = R1 + TTILDE1.BRACKET21.R2.T1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, R2TEMP, LIMX,
     + T1TEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET21, LIMX,
     + TRTEMP, LIMX, ZERO, TRTEMP2, LIMX)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDE1TEMP, LIMX,
     + TRTEMP2, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + R1TEMP, LIMX, UNITY, TRTEMP, LIMX)
      CALL ZCOPY(LIMX*LIMX, TRTEMP, 1, R, 1)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP2, 1)


c$$$ RTILDE = RTILDE2 + T2.BRACKET12.RTILDE1.TTILDE2
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, RTILDE1TEMP, LIMX,
     + TTILDE2TEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET12, LIMX,
     + TRTEMP, LIMX, ZERO, TRTEMP2, LIMX)
      CALL ZCOPY(LIMX*LIMX, ALLZERO, 1, TRTEMP, 1)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, T2TEMP, LIMX,
     + TRTEMP2, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + RTILDE2TEMP, LIMX, UNITY, TRTEMP, LIMX)
      CALL ZCOPY(LIMX*LIMX, TRTEMP, 1, RTILDE, 1)
	  
      RETURN
      END
	  
c$$$ Create the 'o' matrix and invert it	  
      SUBROUTINE FILLOANDINVERT(O, IO, LIMX)
	  
      INTEGER LIMX, I
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX)
	  
c$$$ Generate O-matrix
c$$$ O is block matrix of 1/sqrt(2) (1,1;i,-i)
         DO I = 1, LIMX
            O(I, I)=1/SQRT(2.0)
            O(I, LIMX+I)=1/SQRT(2.0)
            O(I+LIMX, I)=DCMPLX(0, 1/SQRT(2.0))
            O(I+LIMX, I+LIMX)=DCMPLX(0 ,-1/SQRT(2.0))
         ENDDO
         CALL ZCOPY(4*LIMX*LIMX, O, 1, IO, 1)
         CALL INVERTMATRIX(IO, 2*LIMX)  
	  
      RETURN
      END
	  
c$$$ Routine to invert a matrix using LAPACK functions with a check	  
      SUBROUTINE INVERTMATRIX(MATRIX, LIMX)
	  
      INTEGER S
      INTEGER LIMX, PIVOT(LIMX, LIMX)
      DOUBLE COMPLEX MATRIX(LIMX, LIMX), WORK(LIMX*LIMX)
	  
      CALL ZGETRF(LIMX, LIMX, MATRIX, LIMX, PIVOT, S)
      CALL ZGETRI(LIMX, MATRIX, LIMX, PIVOT, WORK, LIMX*LIMX, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'Non-invertable matrix with S=', S
         STOP
      END IF
	  
      RETURN
      END
	  
c$$$ Function that creates a matrix 	  
      SUBROUTINE CALCMULT(MULT, LIMX, LIMY, WRAPX, MODD, MEVEN, E)
	  
      INTEGER LIMX, LIMY, WRAPX
      INTEGER I/1/
      DOUBLE PRECISION E
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX),
     + MULT(2*LIMX, 2*LIMX)


c$$$ First row is even - WRAPX makes no diff, second row not, etc.
c$$$ - what matters is which row it is centred on
c$$$ There are 2 transfer matrices to generate
c$$$ There are 4 block submatrices to fill
c$$$ MODD doesn't depend on xwrapping, MEVEN does.   '
      DO I = 1, LIMX
c$$$ Fill top-right submatrix
         MODD(I, LIMX+I)=1
         MEVEN(I, LIMX+I)=1
c$$$ Fill bottom-left submatrix
         MODD(I+LIMX, I)=-1
         MEVEN(I+LIMX, I)=-1
c$$$ Fill bottom-right submatrix
         MODD(LIMX+I,LIMX+I)=E
         MEVEN(LIMX+I,LIMX+I)=E
c$$$ Smart trick to set alternating adjacent value to -1
         MODD(LIMX+I,LIMX+I+(2*MOD(I,2)-1))=-1
c$$$ Trick to set end values dependent on WRAPX
c$$$ Only the ends matter with regards to the WRAPX effect
         MEVEN(LIMX+I, 2*LIMX+1-I)=-1*((-1*(I/(I*I)))+1)*((-1*
     +    (((LIMX-I)+1)/(((LIMX-I)+1)*((LIMX-I)+1))))+1) + (I/(I
     +    *I))*(-1*WRAPX) +(((LIMX-I)+1)/(((LIMX-I)+1)*((LIMX-I)
     +    +1)))*(-1*WRAPX)
      END DO
      IF (MOD(LIMY,2) .EQ. 1) THEN
c$$$         MULT=MODD
         CALL ZCOPY(4*LIMX*LIMX, MODD, 1, MULT, 1)
      ELSE
c$$$         MULT=MEVEN
         CALL ZCOPY(4*LIMX*LIMX, MEVEN, 1, MULT, 1)
      END IF
	  
      RETURN
      END
	  
c$$$ Single value decomposition of matrix

      SUBROUTINE SV_DECOMP(LIMX, MATRIX, OUTPUTS)
	  
      INTEGER LIMX, MSIZE, S
      DOUBLE PRECISION SVALS(LIMX), OUTPUTS(LIMX), RWORK(5*LIMX)
      DOUBLE COMPLEX MATRIX(LIMX, LIMX), TEMP2(LIMX, LIMX),
     + SVCPY(LIMX, LIMX), WORK(4*LIMX*LIMX)
	  
      MSIZE=4*LIMX*LIMX

c$$$ make copy of matrix for SVD since it is destroyed
c$$$      SVCPY=MATRIX
      CALL ZCOPY(LIMX*LIMX, MATRIX, 1, SVCPY, 1)
      CALL ZGESVD('N', 'N', LIMX, LIMX, SVCPY, LIMX, SVALS, TEMP2,
     + LIMX, TEMP2, LIMX , WORK, MSIZE, RWORK, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'SVD failed with S=', S
         STOP
      END IF

c$$$      OUTPUTS=SVALS
      CALL ZCOPY(LIMX, SVALS, 1, OUTPUTS, 1)
	  
      RETURN
      END

c$$$ Routine to print output to the screen	  
      SUBROUTINE PRINTVECTOR(INPUT, LIMX, MNAME)
	  
      INTEGER LIMX, I
      DOUBLE PRECISION INPUT(LIMX)
      CHARACTER*2 MNAME
	  
      WRITE (*,200) MNAME, (INPUT(I)*INPUT(I), I = 1, LIMX)
 200  FORMAT (A, ' =', ES15.5E2, ES15.5E2)
	  
      RETURN
      END
	  
      SUBROUTINE PRINTT(T, LIMX, MNAME)

      INTEGER LIMX
      DOUBLE COMPLEX T(LIMX, LIMX)
      CHARACTER*3 MNAME

      WRITE (*,400) MNAME, REAL(T(1, 1)), REAL(T(1, 2))
      WRITE (*,400) MNAME, REAL(T(2,1)), REAL(T(2, 2))
      WRITE (*,401) MNAME, 'i', AIMAG(T(1, 1)), AIMAG(T(1, 2))
      WRITE (*,401) MNAME, 'i', AIMAG(T(2,1)), AIMAG(T(2, 2))

 400  FORMAT (A, ES15.5E2, ES15.5E2)
 401  FORMAT (A, A, ES15.5E2, ES15.5E2)

      RETURN
      END

c$$$  Routine to check for unitary matrices
    
      DOUBLE PRECISION FUNCTION CHECKUNI(LIMX, T,R,TTILDE,RTILDE)
      IMPLICIT NONE

      INTEGER LIMX, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      DOUBLE COMPLEX T(LIMX,LIMX), BETA/-1/,ALPHA/1/,
     + R(LIMX,LIMX),TTILDE(LIMX,LIMX),RTILDE(LIMX,LIMX),
     + U(LIMX*LIMX,LIMX*LIMX), CHECK(LIMX*LIMX,LIMX*LIMX)
 

c$$$ test cases of t,r,t~ r~
    
c$$$ Fills a unit matrix

      DO X = 1, 2*LIMX
         DO Y = 1, 2*LIMX
            CHECK(X, Y) = (0.0, 0.0)
         END DO
       CHECK(X, X) = 1.0
      END DO
     
      DO X=1, LIMX
       DO Y=1, LIMX
c$$$ Top left 
       U(X,Y)=T(X,Y)      
c$$$  Bottom left
	 U(X+LIMX,Y)=R(X,Y)     
c$$$ Top right    
        U(X,Y+LIMX)=RTILDE(X,Y)
c$$$ Bottom right
	  U(X+LIMX,Y+LIMX)=TTILDE(X,Y)
        END DO
      END DO

c$$$ Zgemm has inbuilt function to find U**H

      CALL ZGEMM('N', 'C', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, U,
     +       2*LIMX, U, 2*LIMX, BETA, CHECK, 2*LIMX)
c$$$ Zlange finds matrix norm
      
      CHECKUNI = ZLANGE('F', 2*LIMX, 2*LIMX, CHECK, 2*LIMX)	
      

      RETURN       
      END
            













	      

