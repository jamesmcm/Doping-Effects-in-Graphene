      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE

      DOUBLE PRECISION CONDUCTANCE
      EXTERNAL CONDUCTANCE
      DOUBLE PRECISION GETTRANS
      EXTERNAL GETTRANS
      DOUBLE PRECISION CHECKUNI2
      DOUBLE PRECISION ZLANGE
      EXTERNAL CHECKUNI2


      
C     For X current,  LIMY MUST be even, LIMX MUST BE >=3
C     FOR Y current,  LIMX should be even if WRAPX = 1
      CHARACTER             CURRENT /'Y'/,
     +                      GAUGE   /'Y'/
       
      INTEGER, PARAMETER :: LIMX  = 3,
     +                      LIMY  = 10,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      VSIZE = LIMX*LIMY
      
      DOUBLE PRECISION      FLUX/0.0/
       
      DOUBLE PRECISION, PARAMETER :: EMIN = -3.0,
     +                               EMAX =  3.0
      INTEGER, PARAMETER ::          NE   =  600
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE  PRECISION TVALS(MAXSIZE)
      INTEGER NTVALS
      DOUBLE PRECISION E, CONDA/-1.0/, G
      INTEGER IE/0/
      DOUBLE COMPLEX V(LIMX,LIMY)
      
      DATA V/VSIZE*0.0/
C     Note that V is different size to every other matrix
c$$$  V - represents potentials of sites in real lattice
    

C$$$  READS COMMAND LINE ARGUMENT AS LIMY
c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY


      DO IE = 0, NE + 1
         E = EMIN + ( (EMAX - EMIN) * IE) / NE
C     Function to fill V here - for now just set to zeroes
         
         CONDA = GETTRANS(CURRENT, GAUGE,
     +                   TVALS, NTVALS,
     +                   LIMX, LIMY, V,
     +                   E,    FLUX,
     +                   WRAPX)
         G    = CONDUCTANCE (TVALS, NTVALS)
         WRITE(*,50) E, G, CONDA

c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
         CALL FLUSH()

      END DO

 50   FORMAT (F15.5,20ES20.5E3)
      STOP
      END
