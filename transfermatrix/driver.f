      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE

      DOUBLE PRECISION CONDUCTANCE
      EXTERNAL CONDUCTANCE
      DOUBLE PRECISION GETTRANS
      EXTERNAL GETTRANS
c$$   NB LIMX CHANGED TO 2
      INTEGER, PARAMETER :: LIMX=2, WRAPY=0, WRAPX=0,
     + MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER F/1/, LIMY/10/
c$$$      CHARACTER*3 VALUE
      DOUBLE PRECISION TVALS(LIMX)

      DOUBLE PRECISION COND/-1.0/, FLUX/0.0/, E/-3/, G

C$$$  READS COMMAND LINE ARGUMENT AS LIMY
c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY


      DO F = 1, 1201

         COND = GETTRANS(TVALS, LIMX, LIMY, E, FLUX, WRAPX)
C$$$  ERROR RETURN TYPE MISMATCH OF FUNCTION CHECKUNI REAL(4)/REAL(8)
c$$$      COND = CHECKUNI(LIMX,T,R,TTILDE,RTILDE)
         G    = CONDUCTANCE (TVALS, LIMX)

c$$$  WRITES ENERGY, CONDUCTANCE, UNITARITY
         WRITE(*,50) E, G, COND
c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
         CALL FLUSH()

c      WRITE(*,60) E,(TVALS(I)*TVALS(I), I = 1, LIMX)

c$$$ 'E' STEPS CONSISTANT WITH ANALYTICAL.C
         E=-5+(F-1.0)*0.0001
      END DO

 50   FORMAT (F8.5,15ES15.5E3)

      STOP
      END
