      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE

      DOUBLE PRECISION CHECKUNI
      EXTERNAL CHECKUNI
      DOUBLE PRECISION CHECKUNI2
      EXTERNAL CHECKUNI2
      DOUBLE PRECISION CHECKUNI3
      EXTERNAL CHECKUNI3
      DOUBLE PRECISION CONDUCTANCE
      EXTERNAL CONDUCTANCE
c$$   NB LIMX CHANGED TO 2       
      INTEGER, PARAMETER :: LIMX=2, WRAPY=0, WRAPX=0,
     + MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER I/1/, K/1/, F/1/,LIMY/10/
      CHARACTER*3 VALUE
      
      DOUBLE PRECISION RVALS(LIMX),  TVALS(LIMX),
     +                 TTVALS(LIMX), RTVALS(LIMX)

      DOUBLE PRECISION COND/-1.0/
      DOUBLE PRECISION E/-3/
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

      DOUBLE COMPLEX CNUM
      DOUBLE PRECISION FLUX/0.0/
C$$$  READS COMMAND LINE ARGUMENT AS LIMY

c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY
      CALL FILLOANDINVERT(O, IO, LIMX, FLUX)
c$$$      CALL ZPRINTM (O,  LIMX, 'O ')	  
      DO F = 1, 1201
c$$$         ARG=5.0
c$$$         CNUM = DCMPLX(0,0)
c$$$         CALL ZPOLAR(ARG, CNUM)
c$$$         PRINT *, CNUM

         CALL CALCMULT(MULT, LIMX, WRAPX, MODD, MEVEN, E, FLUX)
c$$$  CALCMULT fills MODD, MEVEN - do multiplication in main loop
c$$$  Must decide whether we want zig-zag or armchair edges
C     For now I have left it as before so I can compare results
c$$$         PRINT *, '-----'
c$$$         CALL ZPRINTM (O,  LIMX, 'OO ')
c$$$         PRINT *, '-----'
c$$$         CALL ZPRINTM (IO,  LIMX, 'IO ')	  
c$$$
c$$$         STOP

         IF (MOD(LIMY,2) .EQ. 1) THEN
C$$$  MULT=MODD
            CALL ZCOPY(4*LIMX*LIMX, MODD, 1, MULT, 1)
         ELSE
C$$$  MULT=MEVEN
            CALL ZCOPY(4*LIMX*LIMX, MEVEN, 1, MULT, 1)
         END IF
c$$$         CALL PRINTM (MODD,  LIMX, 'MO ')	  
c$$$         CALL PRINTM (MEVEN, LIMX, 'ME ')

c$$$         CALL FILLOANDINVERT(O, IO, LIMX)
c$$$  This was moved outside the loop as it is unnecessary here, at the moment

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

c$$$  WRITES ENERGY, CONDUCTANCE, UNITARITY
      WRITE(*,50) E, G, COND
  
c      WRITE(*,60) E,(TVALS(I)*TVALS(I), I = 1, LIMX)

c$$$ 'E' STEPS CONSISTANT WITH ANALYTICAL.C
      E=E+0.005
      END DO
      
         


C$$$ SO T^2 + R^2 =1 FOR SVD VALUES, ALSO VERIFIED WITH R~ AND T~


      

 50   FORMAT (F8.5,15ES15.5E3)      
 60   FORMAT (15ES15.5E3)
     

      STOP
      END
