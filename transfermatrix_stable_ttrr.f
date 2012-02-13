c$$$ Want to loop over different energies and produce T^2 coefficients, check they match with analytical results

      PROGRAM TRANSFERMATIXTWO
      INTEGER, PARAMETER :: LIMX=2, LIMY=10, WRAPY=0, WRAPX=1,
     + MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER PIVOT(2*LIMX, 2*LIMX), PIVOT2(LIMX, LIMX)
      INTEGER*4 I/1/, J/1/, S/9/, K/1/, F/1/
      DOUBLE PRECISION SVALS(LIMX), RWORK(5*LIMX), RVALS(LIMX),
     + TVALS(LIMX), E/-5/, TTVALS(LIMX), RTVALS(LIMX), SQUARE
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

c$$$ WRITE (*,70) "LIMX:", LIMX, " LIMY:", LIMY
c$$$      DO F = 1, 1001
      E=0.0
c$$$ WRITE (*,50) "E value:", E

      CALL CALCMULT(MULT, LIMX, LIMY, WRAPX, MODD, MEVEN, E)

      CALL FILLOANDINVERT(O, IO, LIMX)
      CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
      CALL GENTANDRINC(LIMX, TTILDE, D, PIVOT2, B, RTILDE, C,
     + R, T, A)

      CALL SV_DECOMP(LIMX, T, TVALS)
      CALL SV_DECOMP(LIMX, R, RVALS)
      CALL SV_DECOMP(LIMX, TTILDE, TTVALS)
      CALL SV_DECOMP(LIMX, RTILDE, RTVALS)
	 
      CALL PRINTVECTOR(TVALS, LIMX, 'T ')
      CALL PRINTVECTOR(RVALS, LIMX, 'R ')
      CALL PRINTVECTOR(TTVALS, LIMX, 'T~')
      CALL PRINTVECTOR(RTVALS, LIMX, 'R~')
		 
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

c$$$         WRITE(*, 888) 'Mult', REAL(MULT(1, 1)), REAL(MULT(2, 1)),
c$$$     +    REAL(MULT(3,1)), REAL(MULT(4, 1))
c$$$         WRITE(*, 888) 'Mult', REAL(MULT(1, 2)), REAL(MULT(2, 2)),
c$$$     +    REAL(MULT(3,2)), REAL(MULT(4, 2))
c$$$         WRITE(*, 888) 'Mult', REAL(MULT(1, 3)), REAL(MULT(2, 3)),
c$$$     +    REAL(MULT(3,3)), REAL(MULT(4, 3))
c$$$         WRITE(*, 888) 'Mult', REAL(MULT(1, 4)), REAL(MULT(2, 4)),
c$$$     +    REAL(MULT(3,4)), REAL(MULT(4, 4))

 888  FORMAT (A, F8.4, F8.4, F8.4, F8.4)

			
         CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
         CALL GENTANDRINC(LIMX, TTILDEINC, D, PIVOT2, B, RTILDEINC, C,
     +    RINC, TINC, A)

         IF (I .EQ. LIMY) THEN
            CALL PRINTT(T, LIMX, 'T  ')
            CALL PRINTT(TTILDE, LIMX, 'Tt ')
            CALL PRINTT(R, LIMX, 'R  ')
            CALL PRINTT(RTILDE, LIMX, 'Rt ')
         END IF
	 
         CALL UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     +    RTILDE, LIMX, RINC)
	 
         CALL SV_DECOMP(LIMX, T, TVALS)
         CALL SV_DECOMP(LIMX, R, RVALS)
         CALL SV_DECOMP(LIMX, TTILDE, TTVALS)
         CALL SV_DECOMP(LIMX, RTILDE, RTVALS)
		 
         WRITE (*,*) ' '	 
         CALL PRINTVECTOR(TVALS, LIMX, 'T ')
         CALL PRINTVECTOR(RVALS, LIMX, 'R ')
         CALL PRINTVECTOR(TTVALS, LIMX, 'T~')
         CALL PRINTVECTOR(RTVALS, LIMX, 'R~')
			
      END DO

c$$$ ################################################################

      WRITE (*,*) ' '
      SQUARE=SQRT(0.5)
      DO I=1, LIMX
         DO J=1, LIMX
            T(I,J) = 0.0
            TINC(I,J) = 0.0
            TTILDE(I,J) = 0.0
            TTILDEINC(I,J) = 0.0
            R(I,J) = 0.0
            RINC(I,J) = 0.0
            RTILDE(I,J) = 0.0
            RTILDEINC(I,J) = 0.0
         END DO
         T(I,I) = SQUARE
         TINC(I,I) = SQUARE
         TTILDE(I,I) = SQUARE
         TTILDEINC(I,I) = SQUARE
         R(I,I) = SQUARE
         RINC(I,I) = SQUARE
         RTILDE(I,I) = SQUARE
         RTILDEINC(I,I) = SQUARE
      END DO
      CALL UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     + RTILDE, LIMX, RINC)
      CALL SV_DECOMP(LIMX, T, TVALS)
      CALL PRINTVECTOR(TVALS, LIMX, 'T ')
      CALL PRINTT(T, LIMX, 'T  ')


c$$$      DO I = 1, LIMY-1
c$$$         WRITE (*, 20) ALPHA
c$$$         WRITE (*, 20) BETA
c$$$         IF (I .GT. 1) THEN
c$$$            CALL UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
c$$$     +        RTILDE, LIMX)
c$$$         ELSE
c$$$		   T=TINC
c$$$            CALL ZCOPY(LIMX*LIMX, TINC, 1, T, 1)
c$$$           TTILDE=TTILDEINC
c$$$            CALL ZCOPY(LIMX*LIMX, TTILDEINC, 1, TTILDE, 1)
c$$$           R=RINC
c$$$            CALL ZCOPY(LIMX*LIMX, RINC, 1, R, 1)
c$$$           RTILDE=RTILDEINC
c$$$            CALL ZCOPY(LIMX*LIMX, RTILDEINC, 1, RTILDE, 1)
c$$$         END IF
c$$$      END DO

		 
c$$$         CALL SV_DECOMP(LIMX, TINC, TVALS)
c$$$         CALL SV_DECOMP(LIMX, RINC, RVALS)

c$$$ DO J=1, LIMX
c$$$ WRITE (*,60) "T^2 value: ", TVALS(J)*TVALS(J),
c$$$ + " R^2 value: ", RVALS(1+LIMX-J)*RVALS(1+LIMX-J)
c$$$ END DO


c$$$         WRITE (*,80) E, (TVALS(I)*TVALS(I), I = 1, LIMX)
c$$$ WRITE (*,80) "R^2 values: ",(RVALS(LIMX-I+1)*RVALS(LIMX-I+1),
c$$$ + I = 1, LIMX)



c$$$         E=E+0.01
c$$$      END DO
c$$$ so T^2 + R^2 =1 for SVD values, also verified with R~ and T~

 20   FORMAT (4F4.0)
 30   FORMAT (F8.4, A, F8.4, A)
 40   FORMAT (F8.4)
 50   FORMAT (A, F6.2)
 60   FORMAT (A, F8.4, A, F8.4)
 70   FORMAT (A, I6, A, I6)
 80   FORMAT (F8.4, ES15.5E2, ES15.5E2)
 
      STOP
      END
	  
 
      SUBROUTINE GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
	 
      INTEGER LIMX
      DOUBLE COMPLEX UNITY, ZERO, A(LIMX, LIMX), B(LIMX, LIMX),
     + C(LIMX, LIMX), D(LIMX, LIMX), MULT(2*LIMX, 2*LIMX),
     + O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX), TEMP(2*LIMX, 2*LIMX),
     + ABCD(2*LIMX, 2*LIMX)
	 
      UNITY = 1.0
      ZERO = 0.0	 
	 
c$$$         write (*, 990) REAL(MULT(1, 1)), REAL(MULT(1, 2)), 
c$$$     +                  REAL(MULT(2,1)), REAL(MULT(2, 2))
c$$$	     write (*, 990) REAL(MULT(3, 1)), REAL(MULT(3, 2)), 
c$$$     +                  REAL(MULT(4,1)), REAL(MULT(4, 2))
c$$$	     write (*, 990) REAL(MULT(1, 3)), REAL(MULT(1, 4)), 
c$$$     +                  REAL(MULT(2,3)), REAL(MULT(2, 4))
c$$$	     write (*, 990) REAL(MULT(3, 3)), REAL(MULT(3, 4)), 
c$$$     +                  REAL(MULT(4,3)), REAL(MULT(4, 4))
	  
         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, UNITY, MULT,
     +         2*LIMX, O, 2*LIMX, ZERO, TEMP, 2*LIMX)
         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, UNITY, IO,
     +         2*LIMX, TEMP, 2*LIMX, ZERO, ABCD, 2*LIMX)	  
c$$$
c$$$ PRINT *, 'ABCD matrix:'
c$$$
c$$$ DO J = 1, 2*LIMX
c$$$ DO I=1, 2*LIMX
c$$$ WRITE (*,30) REAL(ABCD(J,I)), ' + ', DIMAG(ABCD(J,I)), 'I'
c$$$ END DO
c$$$ PRINT *, '----'
c$$$ END DO
c$$$ This is Fortran 90 syntax, remove in future revision when BLAS/LAPACK subroutine is found
         A=ABCD(1:LIMX, 1:LIMX)
         B=ABCD((LIMX+1):2*LIMX, 1:LIMX)
         C=ABCD(1:LIMX, LIMX+1:2*LIMX)
         D=ABCD((LIMX+1):2*LIMX, (LIMX+1):2*LIMX)
c$$$         write (*, 991) REAL(A(1, 1)), REAL(A(1, 2)), REAL(A(2, 1)), REAL(A(2, 2))
c$$$         write (*, 991) AIMAG(A(1, 1)), AIMAG(A(1, 2)), AIMAG(A(2, 1)), AIMAG(A(2, 2))
c$$$         write (*, 992) REAL(B(1, 1)), REAL(B(1, 2)), REAL(B(2, 1)), REAL(B(2, 2))
c$$$         write (*, 992) AIMAG(B(1, 1)), AIMAG(B(1, 2)), AIMAG(B(2, 1)), AIMAG(B(2, 2))
c$$$         write (*, 993) REAL(C(1, 1)), REAL(C(1, 2)), REAL(C(2, 1)), REAL(C(2, 2))
c$$$         write (*, 993) AIMAG(C(1, 1)), AIMAG(C(1, 2)), AIMAG(C(2, 1)), AIMAG(C(2, 2))
c$$$         write (*, 994) REAL(D(1, 1)), REAL(D(1, 2)), REAL(D(2, 1)), REAL(D(2, 2))
c$$$         write (*, 994) AIMAG(D(1, 1)), AIMAG(D(1, 2)), AIMAG(D(2, 1)), AIMAG(D(2, 2))
c$$$         write (*, 995) 
 990     FORMAT ('M=',ES15.5E4,' ',ES15.5E4,' ',ES15.5E4,' ',ES15.5E4)
 991     FORMAT ('A=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 992     FORMAT ('B=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 993     FORMAT ('C=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 994     FORMAT ('D=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 995     FORMAT ('---------')
 996     FORMAT ('T=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)

c$$$ I have verified that AD-BC=1 (identity matrix) as expected
      RETURN
      END
	  
	  
      SUBROUTINE GENTANDRINC(LIMX, TTILDEINC, D, PIVOT2, B, RTILDEINC,
     + C, RINC, TINC, A)
	 
      INTEGER LIMX
      INTEGER PIVOT2(LIMX, LIMX)
      DOUBLE COMPLEX UNITY, ZERO, A(LIMX, LIMX),
     + C(LIMX, LIMX), D(LIMX, LIMX), TEMP2(LIMX, LIMX),
     + TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX), RINC(LIMX, LIMX),
     + RTILDEINC(LIMX, LIMX), WORK2(LIMX*LIMX), B(LIMX, LIMX)
      INTEGER*4 S/9/

      UNITY = 1.0
      ZERO = 0.0
	 
c$$$ T~ = D^-1
      CALL ZCOPY(LIMX*LIMX, D, 1, TTILDEINC, 1)
      CALL ZGETRF(LIMX, LIMX, TTILDEINC, LIMX, PIVOT2, S)
      CALL ZGETRI(LIMX, TTILDEINC, LIMX, PIVOT2, WORK2, LIMX*LIMX, S)
c$$$ R~ = BD^-1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, B,
     +      LIMX, TTILDEINC, LIMX, ZERO, RTILDEINC, LIMX)
c$$$ R = -D^-1 C
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDEINC,
     +      LIMX, C, LIMX, ZERO, RINC, LIMX)
c$$$ R=-R
c$$$ T=(A-)? BD^-1 C
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDEINC,
     +      LIMX, C, LIMX, ZERO, TEMP2, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, B,
     +      LIMX, TEMP2, LIMX, ZERO, TINC, LIMX)
         TINC=A-TINC
c$$$         write (*, 981) REAL(TINC(1, 1)), REAL(TINC(1, 2)), REAL(TINC(2, 1)), REAL(TINC(2, 2))
c$$$         write (*, 981) AIMAG(TINC(1, 1)), AIMAG(TINC(1, 2)), AIMAG(TINC(2, 1)), AIMAG(TINC(2, 2))
c$$$         write (*, 982) REAL(TTILDEINC(1, 1)), REAL(TTILDEINC(1, 2)), REAL(TTILDEINC(2, 1)), REAL(TTILDEINC(2, 2))
c$$$         write (*, 982) AIMAG(TTILDEINC(1, 1)), AIMAG(TTILDEINC(1, 2)), AIMAG(TTILDEINC(2, 1)), AIMAG(TTILDEINC(2, 2))
c$$$         write (*, 983) REAL(RINC(1, 1)), REAL(RINC(1, 2)), REAL(RINC(2, 1)), REAL(RINC(2, 2))
c$$$         write (*, 983) AIMAG(RINC(1, 1)), AIMAG(RINC(1, 2)), AIMAG(RINC(2, 1)), AIMAG(RINC(2, 2))
c$$$         write (*, 984) REAL(RTILDEINC(1, 1)), REAL(RTILDEINC(1, 2)), REAL(RTILDEINC(2, 1)), REAL(RTILDEINC(2, 2))
c$$$         write (*, 984) AIMAG(RTILDEINC(1, 1)), AIMAG(RTILDEINC(1, 2)), AIMAG(RTILDEINC(2, 1)), AIMAG(RTILDEINC(2, 2))
c$$$         write (*, 985) 
 980     FORMAT ('T=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 981     FORMAT ('T=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 982     FORMAT ('T~=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 983     FORMAT ('R=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 984     FORMAT ('T~=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
 985     FORMAT ('---------')
 986     FORMAT ('T=', f8.4, ' ', f8.4, ' ', f8.4, ' ', f8.4)
		 
      RETURN
      END
	  

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
     + BRACKET21(LIMX, LIMX), TRTEMP(LIMX, LIMX), ALLONE(LIMX, LIMX),
     + UNITMATRIX(LIMX, LIMX)
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
            ALLONE(X, Y) = 1.0
         END DO
         UNITMATRIX(X, X) = 1.0
      END DO

c$$$ BRACKET12 = (1 - RTILDE1.R2)^-1
c$$$      CALL PRINTT(RTILDE1TEMP, LIMX, 'Rt1')
c$$$      CALL PRINTT(R2TEMP, LIMX, 'R2 ')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, RTILDE1TEMP, LIMX,
     + R2TEMP, LIMX, ZERO, BRACKET12, LIMX)
c$$$      CALL PRINTT(BRACKET12, LIMX, 'b12')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + ALLONE, LIMX, -1*UNITY, BRACKET12, LIMX)
c$$$      CALL PRINTT(BRACKET12, LIMX, 'b12')
      CALL INVERTMATRIX(BRACKET12, LIMX)
c$$$      CALL PRINTT(BRACKET12, LIMX, 'B12')
c$$$ BRACKET21 = (1 - R2.RTILDE1)^-1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, R2TEMP, LIMX,
     + RTILDE1TEMP, LIMX, ZERO, BRACKET21, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + ALLONE, LIMX, -1*UNITY, BRACKET21, LIMX)
      CALL INVERTMATRIX(BRACKET21, LIMX)
c$$$ T = T2.BRACKET12.T1
c$$$      CALL PRINTT(T1TEMP, LIMX, 'T1t')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET12, LIMX,
     + T1TEMP, LIMX, ZERO, TRTEMP, LIMX)
c$$$      CALL PRINTT(T2TEMP, LIMX, 'T2t')
c$$$      CALL PRINTT(TRTEMP, LIMX, 'TRt')
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, T2TEMP, LIMX,
     + TRTEMP, LIMX, ZERO, T, LIMX)
c$$$      CALL PRINTT(T, LIMX, 'T  ')
c$$$ TTILDE = TTILDE1.BRACKET21.TTILDE2
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET21, LIMX,
     + TTILDE2TEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDE1TEMP, LIMX,
     + TRTEMP, LIMX, ZERO, TTILDE, LIMX)
c$$$ R = R1 + TTILDE1.BRACKET21.R2.T1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, R2TEMP, LIMX,
     + T1TEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET21, LIMX,
     + TRTEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDE1TEMP, LIMX,
     + TRTEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + R1TEMP, LIMX, UNITY, TRTEMP, LIMX)
      CALL ZCOPY(LIMX*LIMX, TRTEMP, 1, R, 1)
c$$$ RTILDE = RTILDE2 + T2.BRACKET12.RTILDE1.TTILDE2
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, RTILDE1TEMP, LIMX,
     + TTILDE2TEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, BRACKET12, LIMX,
     + TRTEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, T2TEMP, LIMX,
     + TRTEMP, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, UNITMATRIX, LIMX,
     + RTILDE2TEMP, LIMX, UNITY, TRTEMP, LIMX)
      CALL ZCOPY(LIMX*LIMX, TRTEMP, 1, RTILDE, 1)
	  
      RETURN
      END
	  
	  
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
	  
	  
      SUBROUTINE INVERTMATRIX(MATRIX, LIMX)
	  
      INTEGER*4 S
      INTEGER LIMX, PIVOT(LIMX, LIMX)
      DOUBLE COMPLEX MATRIX(LIMX, LIMX), WORK(LIMX*LIMX)
	  
      CALL ZGETRF(LIMX, LIMX, MATRIX, LIMX, PIVOT, S)
      CALL ZGETRI(LIMX, MATRIX, LIMX, PIVOT, WORK, LIMX*LIMX, S)
      IF (S .NE. 0) THEN
         WRITE (*,*) 'Non-invertable matrix'
         STOP
      END IF
	  
      RETURN
      END
	  
	  
      SUBROUTINE CALCMULT(MULT, LIMX, LIMY, WRAPX, MODD, MEVEN, E)
	  
      INTEGER LIMX, LIMY, WRAPX
      INTEGER*4 I/1/
      DOUBLE PRECISION E
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX),
     + MULT(2*LIMX, 2*LIMX)


c$$$ First row is even - WRAPX makes no diff, second row not, etc.
c$$$ - what matters is which row it is centred on
c$$$ There are 2 transfer matrices to generate
c$$$ There are 4 block submatrices to fill
c$$$ MODD doesn't depend on xwrapping, MEVEN does.
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
	  

      SUBROUTINE SV_DECOMP(LIMX, MATRIX, OUTPUTS)
	  
      INTEGER LIMX, MSIZE
      DOUBLE PRECISION SVALS(LIMX), OUTPUTS(LIMX), RWORK(5*LIMX)
      DOUBLE COMPLEX MATRIX(LIMX, LIMX), TEMP2(LIMX, LIMX),
     + SVCPY(LIMX, LIMX), WORK(4*LIMX*LIMX)
	  
      MSIZE=4*LIMX*LIMX

c$$$ make copy of matrix for SVD since it is destroyed
c$$$      SVCPY=MATRIX
      CALL ZCOPY(LIMX*LIMX, MATRIX, 1, SVCPY, 1)
      CALL ZGESVD('N', 'N', LIMX, LIMX, SVCPY, LIMX, SVALS, TEMP2,
     + LIMX, TEMP2, LIMX , WORK, MSIZE, RWORK, S)

c$$$      OUTPUTS=SVALS
      CALL ZCOPY(LIMX*LIMX, SVALS, 1, OUTPUTS, 1)
	  
      RETURN
      END

	  
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

      WRITE (*,*) MNAME, REAL(T(1, 1)), REAL(T(1, 2))
      WRITE (*,*) MNAME, REAL(T(2,1)), REAL(T(2, 2))
      WRITE (*,*) MNAME, 'i', AIMAG(T(1, 1)), AIMAG(T(1, 2))
      WRITE (*,*) MNAME, 'i', AIMAG(T(2,1)), AIMAG(T(2, 2))

      RETURN
      END
