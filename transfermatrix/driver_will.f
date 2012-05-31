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
       
      INTEGER, PARAMETER :: LIMX  = 100,
     +                      LIMY  = 173,
     +                      WRAPX = 0,
     +                      WRAPY = 0,
     +                      VSIZE = LIMX*LIMY
      
      DOUBLE PRECISION      FLUX/0.0/
       
      DOUBLE PRECISION, PARAMETER :: EMIN = -0.92,
     +                               EMAX =  1.02
      INTEGER, PARAMETER ::          NE   =  40
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE PRECISION TVALS(MAXSIZE)
      INTEGER NTVALS
      DOUBLE PRECISION E, CONDA/-1.0/, G(NE), GMEAN(NE), GERROR(NE)
      INTEGER IE/0/, J, K, H
      DOUBLE PRECISION V(LIMX,LIMY), KREAL
c      DOUBLE PRECISION POT /2.0/, POTA /0.2/, POTB /0.3/
c      DOUBLE PRECISION HEIGHT /10.0/, WID /7.0/, GWID /0.01/
      INTEGER, PARAMETER :: NOVAC = 1000000
      DOUBLE PRECISION, PARAMETER :: U = 100.0,
     +                               ALAT = 0.0,
     +                               BLAT = 0.0,
     +                               ZERO = 0.0
      
      DATA V/VSIZE*0.0/
C     Note that V is different size to every other matrix
c$$$  V - represents potentials of sites in real lattice
c$$$
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
      END DO

      NTVALS = 1
      DO K = 1, NOVAC
c         CALL GENVAC(V, LIMX, LIMY, U, ALAT, BLAT)
         DO IE = 0, NE - 1
            E = EMIN + ( (EMAX - EMIN) * IE) / NE
            CALL FAKEGETTRANS(TVALS, NTVALS)
C      CALL TFOUR(V,POT,GWID,LIMX,LIMY)
C     Function to fill V here - for now just set to zeroes
         
c         CONDA = GETTRANS(CURRENT, GAUGE,
c     +                   TVALS, NTVALS,
c     +                   LIMX, LIMY, V,
c     +                   E,    FLUX,
c     +                   WRAPX)
            G(IE)    = CONDUCTANCE(TVALS, NTVALS)
c$$$         WRITE(*,50) E, G, CONDA

c$$$     This is gfortran function to flush the output
c$$$     so that the data are written to the file immediately
c$$$     If you are not using gfortran, and cannot compile this,
c$$$     comment it out -- AVS
            CALL FLUSH()
         END DO
         KREAL = K
         CALL INCREMENTCOND(GMEAN, GERROR, G, KREAL, NE)

         OPEN(UNIT = 1, FILE = "conductance.dat")
c$$$ Writes parameters to the file with the intention for
c$$$ extending later.
            WRITE (1,40) K, U, ALAT, BLAT, ' ', CURRENT, ' ', GAUGE,
     +                  LIMX, LIMY, WRAPX, WRAPY, NE, EMIN, EMAX, FLUX
         DO H = 0, NE - 1
            WRITE (1,60) (EMIN+((EMAX-EMIN)*H)/NE), GMEAN(H+1),
     +                   GERROR(H+1)
         END DO
         CLOSE(1)
      END DO

 40   FORMAT (I15, 3F15.5, 4A, 5I15, 3F15.5)
 50   FORMAT (F15.5,20ES20.5E3)
 60   FORMAT (F15.5, 2ES20.5E3)
      STOP
      END


