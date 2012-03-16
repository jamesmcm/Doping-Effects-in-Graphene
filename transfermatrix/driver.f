      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE

      DOUBLE PRECISION CONDUCTANCE
      EXTERNAL CONDUCTANCE
      DOUBLE PRECISION GETTRANS
      EXTERNAL GETTRANS
c$$$ CHECKUNI2 IS SLIGHTLY FASTER
      DOUBLE PRECISION CHECKUNI2
      EXTERNAL CHECKUNI2
      CHARACTER*3 A,B,C,D,MA,MB,MC,MD
      

      INTEGER, PARAMETER :: LIMX=2, WRAPY=0, WRAPX=0,
     + MSIZE=4*LIMX*LIMX, M2SIZE=LIMX*LIMX
      INTEGER F/1/, LIMY/10/,I/1/
c$$$      CHARACTER*3 VALUE
      DOUBLE PRECISION TVALS(LIMX)
c$$$ E NEEDS TO BE THE SAME AS INCREMENT
      DOUBLE PRECISION CONDA/-1.0/,CONDB/-1.0/, FLUX/0.0/, E/-3/, G
      DOUBLE COMPLEX T(LIMX,LIMX),R(LIMX,LIMX),TTILDE(LIMX,LIMX),
     + RTILDE(LIMX,LIMX)

C$$$  READS COMMAND LINE ARGUMENT AS LIMY
c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY


      DO F = 1, 1201
        
c$$$  CHECKING THE NEW GENABCD *
         CONDA = GETTRANS(TVALS, LIMX, LIMY, E, FLUX, WRAPX)
C      CALL ZPRINTM(A,LIMX,'A')
C      CALL ZPRINTM(MA,LIMX,'MA')


C$$$  *END OF CHECKING NEW GENABCD

c         CONDB = CHECKUNI2(LIMX,T,R,TTILDE,RTILDE)
         G    = CONDUCTANCE (TVALS, LIMX)

c$$$  WRITES ENERGY, CONDUCTANCE, UNITARITY
         WRITE(*,50) E, G, CONDA

c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
         CALL FLUSH()

c      WRITE(*,60) E,(TVALS(I)*TVALS(I), I = 1, LIMX),CONDB

c$$$ 'E' STEPS CONSISTANT WITH ANALYTICAL.C
         E=-3+(F-1.0)*0.0001
      END DO

 50   FORMAT (F15.5,20ES20.5E3)
 60   FORMAT (F15.5,20ES20.5E3)

      STOP
      END
