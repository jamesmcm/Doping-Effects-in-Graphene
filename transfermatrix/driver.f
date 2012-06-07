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


      DOUBLE PRECISION ADDUPV
      EXTERNAL ADDUPV

      
C     For X current,  LIMY MUST be even, LIMX MUST BE >=3
C     FOR Y current,  LIMX should be even if WRAPX = 1
      CHARACTER             CURRENT /'Y'/,
     +                      GAUGE   /'Y'/
       
<<<<<<< HEAD
      INTEGER, PARAMETER :: LIMX  = 1000,
     +                      LIMY  = 1000,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      TIMES = 1000
=======
      INTEGER, PARAMETER :: LIMX  = 10,
     +                      LIMY  = 10,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      VSIZE = LIMX*LIMY
>>>>>>> upstream/potential
      
      DOUBLE PRECISION      FLUX/0.0/
       
      DOUBLE PRECISION, PARAMETER :: EMIN = -3,
     +                               EMAX =  3
      INTEGER, PARAMETER ::          NE   =  2000
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE PRECISION TVALS(MAXSIZE)
      INTEGER NTVALS
      DOUBLE PRECISION E, CONDA/-1.0/, G
<<<<<<< HEAD
      INTEGER IE/0/

      DOUBLE PRECISION V(LIMX,LIMY), TOTALV(TIMES), MEANV, DEVV
      DOUBLE PRECISION, PARAMETER :: U    = 1.0,
     +                               ALAT = 0.5,
     +                               BLAT = 0.5
    

=======
      INTEGER IE/0/, J, K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT /0.5/, POTA /0.2/, POTB /0.3/
      DOUBLE PRECISION YDIM /0.0/, XDIM /5.0/, GWID /0.01/
      
      DATA V/VSIZE*0.0/
C     Note that V is different size to every other matrix
c$$$  V - represents potentials of sites in real lattice
c$$$
>>>>>>> upstream/potential
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

<<<<<<< HEAD
c      DO IE = 0, NE + 1
c         E = EMIN + ( (EMAX - EMIN) * IE) / NE
c         CONDA = GETTRANS(CURRENT, GAUGE,
c     +                   TVALS, NTVALS,
c     +                   LIMX, LIMY,
c     +                   E,    FLUX,
c     +                   WRAPX)
c         G    = CONDUCTANCE (TVALS, NTVALS)
c         WRITE(*,50) E, G, CONDA
=======
c$$$      CALL TONE(V, POT, LIMX, LIMY)
c$$$      CALL TTWO(V, POT, (-1.0*POT), LIMX, LIMY)

c      CALL TVERTADD(V,POT,XDIM,YDIM,LIMX,LIMY)  

      CALL THORADD(V,POT,XDIM,YDIM,LIMX,LIMY)

c$$$      DO J=1, LIMX
c$$$      WRITE(*,310) (V(J, K), K=1, LIMY)  
c$$$      ENDDO
c$$$ 310  FORMAT (100F10.6)
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
c$$$         WRITE(*,50) E, G, CONDA
>>>>>>> upstream/potential

c       WRITE(*,60) E, (TVALS(K), K=1, NTVALS)  
      
c$$$   
c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
         CALL FLUSH()

c      END DO

 50   FORMAT (F15.5,20ES20.5E3)
c$$$  Format 10F is dependant on NTVALS. 
 60   FORMAT (F15.5,16F10.6)

      STOP
      END

