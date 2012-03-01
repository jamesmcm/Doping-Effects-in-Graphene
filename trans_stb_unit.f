C$$$ WANT TO LOOP OVER DIFFERENT ENERGIES AND PRODUCE T^2 COEFFICIENTS, CHECK THEY MATCH WITH ANALYTICAL RESULTS

      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE
      
      DOUBLE PRECISION CHECKUNI
      DOUBLE PRECISION CHECKUNI2
      DOUBLE PRECISION CHECKUNI3
      DOUBLE PRECISION CONDUCTANCE
c$$   NB LIMX CHANGED TO 2       
      INTEGER, PARAMETER :: LIMX=100, WRAPY=0, WRAPX=1,
     + MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER I/1/, K/1/, F/1/,LIMY/100/
      CHARACTER*3 VALUE
      
      DOUBLE PRECISION RVALS(LIMX),  TVALS(LIMX),
     +                 TTVALS(LIMX), RTVALS(LIMX)

      DOUBLE PRECISION COND/-1.0/
      DOUBLE PRECISION E/-5/
      DOUBLE PRECISION G
      DOUBLE COMPLEX   ZEROC/0.0/, ONEC/1.0/
      DOUBLE PRECISION DLVAL
      
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX),
     +               MULT(2*LIMX, 2*LIMX), TEMP(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX),
     +               ABCD(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX T(LIMX, LIMX),    TTILDE(LIMX, LIMX),
     +               R(LIMX, LIMX),    RTILDE(LIMX, LIMX),
     +               TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     +               RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)
      DOUBLE COMPLEX O(2 *LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX)

      DATA MODD/MSIZE*0.0/, MEVEN/MSIZE*0.0/, O/MSIZE*0.0/,
     +     IO/MSIZE*0.0/,   TEMP/MSIZE*0.0/

C$$$  READS COMMAND LINE ARGUMENT AS LIMY

c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY
      DO F = 1, 100000
      

         CALL CALCMULT(MULT, LIMX, LIMY, WRAPX, MODD, MEVEN, E)

         CALL FILLOANDINVERT(O, IO, LIMX)
         CALL GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
         CALL GENTANDRINC(LIMX, T, R, TTILDE, RTILDE, A, B, C, D) 
C         CALL PRINTT (T, LIMX, 'T  ')
C         CALL PRINTT (TTILDE, LIMX, 'T~ ')
C         CALL PRINTT (R, LIMX, 'R  ')
C         CALL PRINTT (RTILDE, LIMX, 'R~ ')
C         COND = CHECKUNI (LIMX, T, R, TTILDE, RTILDE)      
C         WRITE (*, *) '1: I=', I, ' ', COND
	

		 
C$$$ ################################################################


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
C           COND = CHECKUNI (LIMX, T, R, TTILDE, RTILDE)      
C           WRITE (*, *) '2: I=', I, ' ', COND
			
         END DO

C$$$ ################################################################
            
         CALL SV_DECOMP(LIMX, T, TVALS)
c         CALL SV_DECOMP(LIMX, R, RVALS)
c         CALL SV_DECOMP(LIMX, TTILDE, TTVALS)
c         CALL SV_DECOMP(LIMX, RTILDE, RTVALS)

C$$$ OUTPUT OF MAIN FUNCTION-COMMENTED OUT WHILE CHECKING FOR UNITARITY
         		 

C$$$         CALL PRINTVECTOR(TVALS, LIMX, 'T ')
C$$$         CALL PRINTVECTOR(RVALS, LIMX, 'R ')
C$$$         CALL PRINTVECTOR(TTVALS, LIMX, 'T~')
C$$$         CALL PRINTVECTOR(RTVALS, LIMX, 'R~')

	 
C       CALL PRINTT (T, LIMX, 'T  ')
C       CALL PRINTT (TTILDE, LIMX, 'T~ ')
C       CALL PRINTT (R, LIMX, 'R  ')
C       CALL PRINTT (RTILDE, LIMX, 'R~ ')
       ZEROC = 0.0 
       ONEC = 1.0 
      CALL ZLASET ('ALL', LIMX, LIMX, ZEROC, ONEC, T, LIMX)
      CALL ZLASET ('ALL', LIMX, LIMX, ZEROC, ONEC, TTILDE, LIMX)
      CALL ZLASET ('ALL', LIMX, LIMX, ZEROC, ZEROC, R, LIMX)
      CALL ZLASET ('ALL', LIMX, LIMX, ZEROC, ZEROC, RTILDE, LIMX)
       
    
C$$$  ERROR RETURN TYPE MISMATCH OF FUNCTION CHECKUNI REAL(4)/REAL(8)      
      COND = CHECKUNI(LIMX,T,R,TTILDE,RTILDE)
      G    = CONDUCTANCE (TVALS, LIMX)

c$$$  WRITES THE VALUE OF THE CONDUCTANCE AND THE 
      WRITE(*,50) E, G, COND
  
c      WRITE(*,60) E,(TVALS(I)*TVALS(I), I = 1, LIMX)

c$$$ 'E' STEPS CONSISTANT WITH ANALYTICAL.C
      E=E+0.0001
      END DO
      
         


C$$$ SO T^2 + R^2 =1 FOR SVD VALUES, ALSO VERIFIED WITH R~ AND T~


      

 50   FORMAT (F8.5,15ES15.5E3)      
 60   FORMAT (15ES15.5E3)
     

      STOP
      END

C$$$ END OF MAIN PROGRAM ###############################################
C$$$                     ###############################################

      DOUBLE PRECISION FUNCTION CONDUCTANCE (TVALS, LIMX)
      INTEGER LIMX
      DOUBLE PRECISION TVALS(LIMX)
      DOUBLE PRECISION DDOT      
      CONDUCTANCE = DDOT (LIMX, TVALS, 1, TVALS, 1)
      RETURN
      END


 
C$$$ ROUTINE TO GENERATE A,B,C,D. 
      SUBROUTINE GENABCD(LIMX, MULT, O, IO, ABCD, A, B, C, D)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX UNITY, ZERO
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX)
      DOUBLE COMPLEX MULT(2*LIMX, 2*LIMX),
     +               O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX),
     +               TEMP(2*LIMX, 2*LIMX),
     + ABCD(2*LIMX, 2*LIMX)
	 
      UNITY = 1.0
      ZERO = 0.0	 
	 
C$$$ ZGEMM MULTIPLIES MATRICES TOGTHER	  

         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, UNITY, MULT,
     +         2*LIMX, O, 2*LIMX, ZERO, TEMP, 2*LIMX)
         CALL ZGEMM('N', 'N', 2*LIMX, 2*LIMX, 2*LIMX, UNITY, IO,
     +         2*LIMX, TEMP, 2*LIMX, ZERO, ABCD, 2*LIMX)	  

C$$$ THIS IS FORTRAN 90 SYNTAX, REMOVE IN FUTURE REVISION WHEN BLAS/LAPACK SUBROUTINE IS FOUND
         A=ABCD(1:LIMX, 1:LIMX)
         B=ABCD((LIMX+1):2*LIMX, 1:LIMX)
         C=ABCD(1:LIMX, LIMX+1:2*LIMX)
         D=ABCD((LIMX+1):2*LIMX, (LIMX+1):2*LIMX)

C$$$ I HAVE VERIFIED THAT AD-BC=1 (IDENTITY MATRIX) AS EXPECTED
      RETURN
      END
	  
C$$$ ROUTINE TO GENERATE T AND R F	  
      SUBROUTINE GENTANDRINC(LIMX, TINC,RINC,TTILDEINC,RTILDEINC,A,B,
     + C,D)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX UNITY, ZERO
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX),
     + TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     + RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)
c     UNITMATRIX is no longer needed     
c     UNITMATRIX(LIMX, LIMX)

      UNITY = 1.0
      ZERO = 0.0
c     Not needed
c     CALL ZLASET ('ALL', LIMX, LIMX, ZERO, ONE, UNITMATRIX, LIMX)
	 
C$$$ T~ = D^-1
      CALL ZCOPY(LIMX*LIMX, D, 1, TTILDEINC, 1)
      CALL INVERTMATRIX(TTILDEINC, LIMX)
C$$$ R~ = BD^-1
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, B,
     +      LIMX, TTILDEINC, LIMX, ZERO, RTILDEINC, LIMX)
C$$$ R = -D^-1 C
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, -1*UNITY, TTILDEINC,
     +      LIMX, C, LIMX, ZERO, RINC, LIMX)
C$$$ R=-R
C$$$ T=(A-)? BD^-1 C
C$$$ Thus, We can simply reuse R here: T = A + B * R
C$$$ Also, removed a slow call to zgemm, replacing it with zcopy
      CALL ZCOPY(LIMX*LIMX, A, 1, TINC, 1)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, UNITY, B,
     +      LIMX, RINC, LIMX, UNITY, TINC, LIMX)

      RETURN
      END
	  
C$$$ WITH INCREASING ENERGY UPDATES T AND R TO NEW VALUES.

      SUBROUTINE UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     + RTILDE, LIMX, RINC)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX T(LIMX, LIMX), TTILDE(LIMX, LIMX),
     +               R(LIMX, LIMX), RTILDE(LIMX, LIMX),
     +               TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     + RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)

C$$$ TEMPORARY LOCAL VARIABLES
      DOUBLE COMPLEX T1TEMP(LIMX, LIMX), TTILDE1TEMP(LIMX, LIMX),
     +               R1TEMP(LIMX, LIMX), RTILDE1TEMP(LIMX, LIMX)
      DOUBLE COMPLEX UNITY, ZERO
      DOUBLE COMPLEX BRACKET12(LIMX, LIMX), BRACKET21(LIMX, LIMX),
     +               TRTEMP(LIMX, LIMX), TRTEMP2(LIMX, LIMX),
     + UNITMATRIX(LIMX, LIMX), ALLZERO(LIMX, LIMX)
	  
      UNITY = 1.0
      ZERO = 0.0
	 
      CALL ZCOPY(LIMX*LIMX, T, 1, T1TEMP, 1)
      CALL ZCOPY(LIMX*LIMX, TTILDE, 1, TTILDE1TEMP, 1)
      CALL ZCOPY(LIMX*LIMX, R, 1, R1TEMP, 1)	  
      CALL ZCOPY(LIMX*LIMX, RTILDE, 1, RTILDE1TEMP, 1)	  
      CALL ZLASET ('ALL', LIMX, LIMX, ZERO, ZERO,  ALLZERO, LIMX)
      CALL ZLASET ('ALL', LIMX, LIMX, ZERO, UNITY, UNITMATRIX, LIMX)

C$$$ BRACKET12 = (1 - RTILDE1.R2)^-1
      CALL ZCOPY (LIMX*LIMX, UNITMATRIX, 1, BRACKET12, 1)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, -UNITY, RTILDE1TEMP, LIMX,
     + RINC, LIMX, UNITY, BRACKET12, LIMX)
      CALL INVERTMATRIX(BRACKET12, LIMX)

C$$$ BRACKET21 = (1 - R2.RTILDE1)^-1
      CALL ZCOPY (LIMX*LIMX, UNITMATRIX, 1, BRACKET21, 1)
      CALL ZGEMM('N', 'N', LIMX, LIMX, LIMX, -UNITY, RINC, LIMX,
     + RTILDE1TEMP, LIMX, UNITY, BRACKET21, LIMX)
      CALL INVERTMATRIX(BRACKET21, LIMX)


C$$$ T = T2.BRACKET12.T1
C$$$ Note: T2*Br12 can be reused later.
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TINC, LIMX,
     +  BRACKET12, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TRTEMP, LIMX,
     +  T1TEMP, LIMX, ZERO, T, LIMX)

C$$$ RTILDE = RTILDE2 + T2.BRACKET12.RTILDE1.TTILDE2
      CALL ZCOPY (LIMX*LIMX, RTILDEINC, 1, RTILDE, 1)
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TRTEMP, LIMX,
     +  RTILDE1TEMP, LIMX, ZERO, TRTEMP2, LIMX)

      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TRTEMP2, LIMX,
     +  TTILDEINC, LIMX, UNITY, RTILDE, LIMX)

C$$$ TTILDE = TTILDE1.BRACKET21.TTILDE2
C$$$ Note: T1~ * Br21 can be reused later. I reorganized the code
C$$$ to make it possible
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TTILDE1TEMP,
     +  LIMX, BRACKET21, LIMX, ZERO, TRTEMP, LIMX)
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TRTEMP,
     +  LIMX, TTILDEINC, LIMX, ZERO, TTILDE, LIMX)

C$$$ R = R1 + TTILDE1.BRACKET21.R2.T1. Note: TRTEMP from above is reused here
      CALL ZCOPY (LIMX*LIMX, R1TEMP, 1, R, 1)
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TRTEMP, LIMX,
     +  RINC, LIMX, ZERO, TRTEMP2, LIMX)
      CALL ZGEMM ('N', 'N', LIMX, LIMX, LIMX, UNITY, TRTEMP2, LIMX,
     +  T1TEMP, LIMX, UNITY, R, LIMX)
	  
      RETURN
      END
	  
C$$$ CREATE THE 'O' MATRIX AND INVERT IT	  
      SUBROUTINE FILLOANDINVERT(O, IO, LIMX)
      IMPLICIT NONE
      INTEGER LIMX, I
      DOUBLE COMPLEX O(2*LIMX, 2*LIMX), IO(2*LIMX, 2*LIMX)
      DOUBLE PRECISION SQRT05
      DOUBLE COMPLEX ZISQRT05
	  
c     It is slightly more efficient to calculat square root once 
      SQRT05 = SQRT(0.5)
      ZISQRT05 = DCMPLX(0, SQRT05)
C$$$ GENERATE O-MATRIX
C$$$ O IS BLOCK MATRIX OF 1/SQRT(2) (1,1;I,-I)
         DO I = 1, LIMX
           O(I, I)=SQRT05
           O(I, LIMX+I)=SQRT05
           O(I+LIMX, I)=ZISQRT05
           O(I+LIMX, I+LIMX)=-ZISQRT05
         ENDDO
         CALL ZCOPY(4*LIMX*LIMX, O, 1, IO, 1)
         CALL INVERTMATRIX(IO, 2*LIMX)  
	  
      RETURN
      END
	  
C$$$ ROUTINE TO INVERT A MATRIX USING LAPACK FUNCTIONS WITH A CHECK	  
      SUBROUTINE INVERTMATRIX(MATRIX, LIMX)
	  
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
	  
C$$$ FUNCTION THAT CREATES A MATRIX 	  
      SUBROUTINE CALCMULT(MULT, LIMX, LIMY, WRAPX, MODD, MEVEN, E)
	  
      INTEGER LIMX, LIMY, WRAPX, SZ/1/
      INTEGER I/1/, NEIGH/1/
      DOUBLE PRECISION E
      DOUBLE COMPLEX ZEROC / 0.0 / 
      DOUBLE COMPLEX MODD(2*LIMX, 2*LIMX), MEVEN(2*LIMX, 2*LIMX),
     + MULT(2*LIMX, 2*LIMX)

      IF ((MOD(LIMX,2) .NE. 0)) THEN
         WRITE (*,*) 'ERROR: LIMX must be even for physical results'
         STOP 
      ENDIF
c$$$  HAMMERTIME! Program terminates here if LIMX is odd

      SZ = 2 * LIMX
      
      CALL ZLASET ('A', SZ, SZ, ZEROC, ZEROC, MODD, SZ)
      CALL ZLASET ('A', SZ, SZ, ZEROC, ZEROC, MEVEN, SZ)
C$$$ FIRST ROW IS EVEN - WRAPX MAKES NO DIFF, SECOND ROW NOT, ETC.
C$$$ - WHAT MATTERS IS WHICH ROW IT IS CENTRED ON
C$$$ THERE ARE 2 TRANSFER MATRICES TO GENERATE
C$$$ THERE ARE 4 BLOCK SUBMATRICES TO FILL
C$$$ MODD DOESN'T DEPEND ON XWRAPPING, MEVEN DOES.   '
      DO I = 1, LIMX
C$$$ FILL TOP-RIGHT SUBMATRIX
         MODD(I, LIMX+I)=1
         MEVEN(I, LIMX+I)=1
C$$$ FILL BOTTOM-LEFT SUBMATRIX
         MODD(I+LIMX, I)=-1
         MEVEN(I+LIMX, I)=-1
C$$$ FILL BOTTOM-RIGHT SUBMATRIX
         MODD(LIMX+I,LIMX+I)=E
         MEVEN(LIMX+I,LIMX+I)=E

C$$$  THE FOLLOWING CODE WAS MODIFIED --- AVS
C$$$  NEIGHBOURING SITE FOR ODD ROW, ON THE LEFT/RIGHT, DEPENDING ON I
         NEIGH = I + (2*MOD(I,2)-1) 
C$$$     WRITE (*, *) '? I = ', I, ' NEIGH = ', NEIGH
C$$$  NEIGHBOUR CAN BE < 0, OR > LIMX. IF WRAPX IS TRUE, THIS INDICATES
C$$$  A VALID SITE. THE FOLLOWING CODE IS A BIT UGLY, AS I AM NOT SURE
C$$$  WHAT IS MOD(-1, N) IN FORTRAN.
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
               NEIGH = NEIGH + LIMX
             ELSE
               NEIGH = MOD (NEIGH - 1, LIMX) + 1
             ENDIF
C$$$         WRITE (*, *) 'PUT: I = ', I, ' NEIGH = ', NEIGH
             MODD(LIMX + I, LIMX + NEIGH) = -1
         END IF
C$$$     NOW REPEAT THE SAME FOR EVEN ROWS, SWAPPING LEFT AND RIGHT
C$$$     AGAIN, THE CODE IS NOW RATHER UGLY.
         NEIGH = I - (2 * MOD(I, 2)  - 1)
C$$$     WRITE (*, *) '? I = ', I, ' NEIGH = ', NEIGH
         IF (((NEIGH.LE.LIMX).AND.(NEIGH.GT.0)).OR.(WRAPX.GT.0)) THEN
           IF (NEIGH.LE.0) THEN
             NEIGH = LIMX
           ELSE
             NEIGH = MOD (NEIGH - 1, LIMX) + 1
           END IF
C$$$       WRITE (*, *) 'PUT: I = ', I, ' NEIGH = ', NEIGH
           MEVEN(LIMX + I, LIMX + NEIGH) = -1
         END IF
      END DO
      IF (MOD(LIMY,2) .EQ. 1) THEN
C$$$         MULT=MODD
         CALL ZCOPY(4*LIMX*LIMX, MODD, 1, MULT, 1)
      ELSE
C$$$         MULT=MEVEN
         CALL ZCOPY(4*LIMX*LIMX, MEVEN, 1, MULT, 1)
      END IF
c$$$      CALL PRINTM (MODD,  LIMX, 'MO ')	  
c$$$      CALL PRINTM (MEVEN, LIMX, 'ME ')	  
      RETURN
      END
	  
C$$$ SINGLE VALUE DECOMPOSITION OF MATRIX

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

C$$$ ROUTINE TO PRINT OUTPUT TO THE SCREEN	  
      SUBROUTINE PRINTVECTOR(INPUT, LIMX, MNAME)
	  
      INTEGER LIMX, I
      DOUBLE PRECISION INPUT(LIMX)
      CHARACTER*2 MNAME
	  
      WRITE (*,200) MNAME, (INPUT(I)*INPUT(I), I = 1, LIMX)
 200  FORMAT (A, ' =', ES15.5E2, ES15.5E2)
	  
      RETURN
      END
	  
      SUBROUTINE PRINTT(T, LIMX, MNAME)

      INTEGER LIMX, J, K
      DOUBLE COMPLEX T(LIMX, LIMX)
      CHARACTER*3 MNAME
      DO J=1, LIMX
        WRITE (*,400) MNAME, (REAL(T(J, K)), K=1, LIMX)
        WRITE (*,401) MNAME, 'I', (AIMAG(T(J, K)), K=1, LIMX)
      END DO
C      WRITE (*,400) MNAME, REAL(T(2,1)), REAL(T(2, 2))
C      WRITE (*,401) MNAME, 'I', AIMAG(T(1, 1)), AIMAG(T(1, 2))
C      WRITE (*,401) MNAME, 'I', AIMAG(T(2,1)), AIMAG(T(2, 2))

 400  FORMAT (A, 100ES15.5E3, ES15.5E3)
 401  FORMAT (A, A, 100ES15.5E3, ES15.5E3)

      RETURN
      END
      
      SUBROUTINE PRINTM(M, LIMX, MNAME)

      INTEGER LIMX, J, K
      DOUBLE COMPLEX M(2*LIMX, 2*LIMX)
      CHARACTER*3 MNAME
      DO J=1, 2*LIMX
        WRITE (*,500) MNAME, (REAL(M(J, K)), K=1, 2*LIMX)
      END DO
C      WRITE (*,400) MNAME, REAL(T(2,1)), REAL(T(2, 2))
C      WRITE (*,401) MNAME, 'I', AIMAG(T(1, 1)), AIMAG(T(1, 2))
C      WRITE (*,401) MNAME, 'I', AIMAG(T(2,1)), AIMAG(T(2, 2))

 500  FORMAT (A, 100F6.2)
 501  FORMAT (A, A, 100ES15.5E3, ES15.5E3)

      RETURN
      END

C$$$  ROUTINE TO CHECK FOR UNITARY MATRICES
    
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
C$$$ TOP LEFT 
       U(X,Y)=T(X,Y)      
C$$$  BOTTOM LEFT
	 U(X+LIMX,Y)=R(X,Y)     
C$$$ TOP RIGHT    
        U(X,Y+LIMX)=RTILDE(X,Y)
C$$$ BOTTOM RIGHT
	  U(X+LIMX,Y+LIMX)=TTILDE(X,Y)
        END DO
      END DO

C$$$ ZGEMM HAS INBUILT FUNCTION TO FIND U**H
C$$   CALL PRINTT (U, 2*LIMX, 'U1 ')
      CALL ZGEMM('N', 'C', 2*LIMX, 2*LIMX, 2*LIMX, ALPHA, U,
     +       2*LIMX, U, 2*LIMX, BETA, CHECK, 2*LIMX)
C$$$ ZLANGE FINDS MATRIX NORM
C$$   CALL PRINTT (CHECK, 2 * LIMX, 'C2 ')
      CHECKUNI = ZLANGE('F', 2*LIMX, 2*LIMX, CHECK, 2*LIMX)	
C$$      WRITE (*, *) 'CHECKUNI:', CHECKUNI   
   
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
c      
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
C$$$     TOP LEFT 
         U(X,Y)=T(X,Y)      
C$$$     BOTTOM LEFT
         U(X+LIMX,Y)=R(X,Y)
C$$$     TOP RIGHT    
         U(X,Y+LIMX)=RTILDE(X,Y)
C$$$     BOTTOM RIGHT
         U(X+LIMX,Y+LIMX)=TTILDE(X,Y)
        END DO
      END DO

C$$$  ZHERK calculates U^H * U, and then subtracts 1.
C$$$  However, this is perfomed in the upper triangular part only
C$$   Therefore, we double the norm.
      CALL ZHERK('U', 'C', 2*LIMX, 2*LIMX, ONED, U, 2*LIMX,
     +           -ONED, CHECK, 2*LIMX)
C$$$ ZLANGE FINDS MATRIX NORM
      CHECKUNI3 = 2 * ZLANGE('F', 2*LIMX, 2*LIMX, CHECK, 2*LIMX)
   
      RETURN       
      END
