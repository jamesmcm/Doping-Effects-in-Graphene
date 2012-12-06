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
       
      INTEGER, PARAMETER :: LIMX  = 50,
     +                      LIMY  = 50,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      VSIZE = LIMX*LIMY
      
      DOUBLE PRECISION      FLUX/0.000/
       
c      DOUBLE PRECISION, PARAMETER :: EMIN = -0.5,
c     +                               EMAX =  0.5
c$$$  NE Length of Conductance,mean,error
      INTEGER, PARAMETER ::          NE   =  35
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE PRECISION TVALS(MAXSIZE)
      INTEGER NTVALS
      DOUBLE PRECISION E, CONDA(NE), CONDA_TMP(NE), G(NE), GMEAN(NE),
     +                 GERROR(NE), LOGG(NE), LOGGERROR(NE),
     +                 LOGGMEAN(NE)
      INTEGER IE/0/, J, K, H
      DOUBLE PRECISION V(LIMX,LIMY), KREAL
c      DOUBLE PRECISION POT /2.0/, POTA /0.2/, POTB /0.3/
c      DOUBLE PRECISION HEIGHT /10.0/, WID /7.0/, GWID /0.01/

c$$$  NOVAC is number of renormalisations
      INTEGER, PARAMETER :: NOVAC = 5
      DOUBLE PRECISION, PARAMETER :: U = 1000.0,
     +                               ALAT = 0.015,
     +                               BLAT = 0.015,
     +                               ZERO = 0.0
      INTEGER SEED, INT_SEED/1/

      DOUBLE PRECISION ELIST(NE)

      ELIST = (/-0.5, -0.4, -0.3, -0.25, -0.225, -0.2, -0.175, -0.15,
     +          -0.125, -0.1, -0.075, -0.05, -0.025, -0.02, -0.015,
     +          -0.01, -0.005, 0.0, 0.005, 0.01, 0.015, 0.02, 0.025,
     +          0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25,
     +          0.3, 0.4, 0.5/)
      
      DATA V/VSIZE*0.0/
C     Note that V is different size to every other matrix
c$$$  V - represents potentials of sites in real lattice

C$$$  READS COMMAND LINE ARGUMENT AS LIMY
c      CALL GETARG(1, VALUE)
c      READ(UNIT=VALUE, FMT=*) LIMY

c      POT=100.0
c$$$      CALL TONE(V, POT, LIMX, LIMY)
c$$$      CALL TTWO(V, POT, (-1.0*POT), LIMX, LIMY)

c      CALL TTHREE(V,POT,WID,HEIGHT,LIMX,LIMY)  

c$$$      DO J=1, LIMX
c$$$      WRITE(*,310) (V(J, K), K=1, LIMY)  
c$$$      ENDDO
c$$$ 310  FORMAT (100F10.6)

      DO K = 1, NE
         GMEAN(K) = 0.0
         GERROR(K) = 0.0
         LOGGMEAN(K) = 0.0
         LOGGERROR(K) = 0.0
         CONDA(K) = 0.0
         CONDA_TMP(K) = 0.0
      END DO

      SEED = INT_SEED
      CALL SAFESEED(SEED)

      NTVALS = 1
      DO K = 1, NOVAC
         CALL GENVACFIX(V, LIMX, LIMY, U, ALAT, BLAT, SEED)
         DO IE = 0, NE - 1
c            E = EMIN + ( (EMAX - EMIN) * IE) / NE
             E = ELIST(IE + 1)
c            CALL FAKEGETTRANS(TVALS, NTVALS)
         
            CONDA_TMP(IE + 1) = GETTRANS(CURRENT, GAUGE,
     +                                   TVALS, NTVALS,
     +                                   LIMX, LIMY, V,
     +                                   E,    FLUX,
     +                                   WRAPX)
            IF (CONDA_TMP(IE + 1) .GT. CONDA(IE + 1)) THEN
               CONDA(IE + 1) = CONDA_TMP(IE + 1)
            END IF

            G(IE + 1) = CONDUCTANCE (TVALS, NTVALS)
            LOGG(IE + 1) = LOG(G(IE + 1))

c$$$         WRITE(*,50) E, G, CONDA

c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
            
         END DO
         KREAL = K
         CALL INCREMENTCOND(GMEAN, GERROR, G, KREAL, NE)
         CALL INCREMENTCOND(LOGGMEAN, LOGGERROR, LOGG, KREAL, NE)
c$$$      OPEN(UNIT=1, FILE="X_100x174_F=0-Y_a=b=0.015_wrap=00.dat")
      OPEN(UNIT=1, FILE="test2.dat")
c$$$ Writes parameters to the file with the intention for
c$$$ extending later.
c$$$  Commented out for convinced for plotting     
         WRITE (1,100) '#', ' Randomisations:', K, ';'
         WRITE (1,101) '#', ' X:', LIMX, ';  Y:', LIMY, ';  WrapX:',
     +                      WRAPX, ';  WrapY:', WRAPY, ';'
         WRITE (1,102) '#', ' nA:', ALAT, ';  nB:', BLAT, ';  U:', U,
     +                      ';'
         WRITE (1,103) '#', ' Current:', CURRENT, ';  Gauge:', GAUGE,
     +                      ';  Flux:', FLUX, ';'
         WRITE (1,100) '#', ' Data points:', NE, ';'
         WRITE (1,100) '#', ' Seed:', INT_SEED, ';'
         WRITE (1,104) '#'

         DO H = 0, NE - 1
            WRITE (1,60) ELIST(H+1), GMEAN(H+1),
     +                   GERROR(H+1), CONDA(H+1),
     +                   LOGGMEAN(H+1), LOGGERROR(H+1)
            CALL FLUSH()
         END DO
         CLOSE(1)
      END DO

 40   FORMAT (I15, 3F15.5, 4A, 5I15, 3F15.5)
 41   FORMAT (A, I15, 3F15.5, 4A, 5I15, 2A, F15.5, I15)
 50   FORMAT (F15.5,20ES20.5E3)
 60   FORMAT (F15.5, 6ES20.5E3)

 100  FORMAT (2A, I15, A)
 101  FORMAT (2A, I15, A, I15, A, I15, A, I15, A)
 102  FORMAT (2A, F15.5, A, F15.5, A, F15.5, A)
 103  FORMAT (6A, F15.5, A)
 104  FORMAT (A)

      STOP
      END


