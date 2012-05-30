      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE

      DOUBLE PRECISION CONDUCTANCE
      EXTERNAL CONDUCTANCE
      DOUBLE PRECISION GETTRANS
      EXTERNAL GETTRANS
      DOUBLE PRECISION CHECKUNI2
      DOUBLE PRECISION ZLANGE
      EXTERNAL CHECKUNI2

      DOUBLE PRECISION ADDUPV
      EXTERNAL ADDUPV

      
C     For X current,  LIMY MUST be even, LIMX MUST BE >=3
C     FOR Y current,  LIMX should be even if WRAPX = 1
      CHARACTER             CURRENT /'Y'/,
     +                      GAUGE   /'X'/
       
      INTEGER, PARAMETER :: LIMX  = 1000,
     +                      LIMY  = 1000,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      TIMES = 1000
      
      DOUBLE PRECISION      FLUX/0.0/
       
      DOUBLE PRECISION, PARAMETER :: EMIN = -3.0,
     +                               EMAX =  3.0
      INTEGER, PARAMETER ::          NE   =  600
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE PRECISION TVALS(MAXSIZE)
      INTEGER NTVALS
      DOUBLE PRECISION E, CONDA/-1.0/, G
      INTEGER IE/0/

      DOUBLE PRECISION V(LIMX,LIMY), TOTALV(TIMES), MEANV, DEVV
      DOUBLE PRECISION, PARAMETER :: U    = 1.0,
     +                               ALAT = 0.5,
     +                               BLAT = 0.5
    

C$$$  READS COMMAND LINE ARGUMENT AS LIMY
c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY
c$$$      CALL CALCMULTXY(LIMY/2, MODD, MEVEN, E, FLUX)

c$$$  Generates the matrix of vacancy locations
      DO IE = 1, TIMES
         CALL GENVAC(V, LIMX, LIMY, U, ALAT, BLAT)
         TOTALV(IE) = ADDUPV(V, LIMX, LIMY)
      END DO
      CALL MEANDEVV(TOTALV, TIMES, MEANV, DEVV)      

      WRITE (*,*) LIMX * LIMY, ',', ALAT, BLAT, ',', MEANV, '+/-', DEVV
c      WRITE (*,*) (TOTALV(IE), IE = 1, TIMES)

c      DO IE = 0, NE + 1
c         E = EMIN + ( (EMAX - EMIN) * IE) / NE
c         CONDA = GETTRANS(CURRENT, GAUGE,
c     +                   TVALS, NTVALS,
c     +                   LIMX, LIMY,
c     +                   E,    FLUX,
c     +                   WRAPX)
c         G    = CONDUCTANCE (TVALS, NTVALS)
c         WRITE(*,50) E, G, CONDA

c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
         CALL FLUSH()

c      END DO

 50   FORMAT (F15.5,20ES20.5E3)
      STOP
      END

