      PROGRAM TRANSFERMATIXTWO
      IMPLICIT NONE

      DOUBLE PRECISION CONDUCTANCE
      EXTERNAL CONDUCTANCE
      DOUBLE PRECISION GETTRANS
      EXTERNAL GETTRANS
      DOUBLE PRECISION CHECKUNI2
      DOUBLE PRECISION ZLANGE
      EXTERNAL CHECKUNI2
      EXTERNAL TONE


      
C     For X current,  LIMY MUST be even, LIMX MUST BE >=3
C     FOR Y current,  LIMX should be even if WRAPX = 1
      CHARACTER             CURRENT /'X'/,
     +                      GAUGE   /'Y'/
       
      INTEGER, PARAMETER :: LIMX  = 11,
     +                      LIMY  = 10,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      VSIZE = LIMX*LIMY
      
      DOUBLE PRECISION      FLUX/0.0/
       
      DOUBLE PRECISION, PARAMETER :: EMIN = -3,
     +                               EMAX =  3
      INTEGER, PARAMETER ::          NE   =  2000
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE  PRECISION TVALS(MAXSIZE)
      INTEGER NTVALS
      DOUBLE PRECISION E, CONDA/-1.0/, G
      INTEGER IE/0/
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT /0.1/, POTA /0.2/, POTB /0.3/
      DOUBLE PRECISION HEIGHT /4.0/, WID /5.0/, GWID /0.01/
      
      DATA V/VSIZE*0.0/
C     Note that V is different size to every other matrix
c$$$  V - represents potentials of sites in real lattice
    

C$$$  READS COMMAND LINE ARGUMENT AS LIMY
c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY

      POT=2.0
c$$$      CALL TONE(V, POT, LIMX, LIMY)
      CALL TTWO(V, POT, (-1.0*POT), LIMX, LIMY)

c$$$      CALL TTHREE(V,POT,WID,HEIGHT,LIMX,LIMY)  
      DO IE = 0, NE + 1
         E = EMIN + ( (EMAX - EMIN) * IE) / NE

        
C      CALL TFOUR(V,POT,GWID,LIMX,LIMY)
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
