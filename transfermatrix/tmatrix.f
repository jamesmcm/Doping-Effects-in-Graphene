      DOUBLE PRECISION FUNCTION GETTRANS(CURRENT, GAUGE,
     +                                   TVALS, NVALS,
     +                                   LIMX, LIMY, V,
     +                                   E, FLUX, WRAPX)
      IMPLICIT NONE
      CHARACTER CURRENT, GAUGE
      
      INTEGER, PARAMETER :: MAXSIZE = 10000
      DOUBLE  PRECISION TVALS(MAXSIZE)
      INTEGER NVALS

      DOUBLE PRECISION E, FLUX
      INTEGER LIMX, LIMY, WRAPX

      DOUBLE PRECISION GETTRANSX
      EXTERNAL GETTRANSX
      DOUBLE PRECISION GETTRANSY
      EXTERNAL GETTRANSY
      DOUBLE COMPLEX V(LIMX, LIMY)

      INTEGER I

      DO I = 1, MAXSIZE
         TVALS(I) = 0.0
      ENDDO

      GETTRANS = -1
      NVALS = 0
      


      IF (CURRENT .EQ. 'X') THEN
         IF (MOD(LIMY, 2) .NE. 0) THEN
            WRITE (*,*) 'ERROR, LIMY must be even for X current!'
            STOP
         END IF
         NVALS = LIMY / 2  
         GETTRANS = GETTRANSX(GAUGE, TVALS, LIMX, NVALS, E, FLUX)
      END IF
       
      IF (CURRENT .EQ. 'Y') THEN
            NVALS = LIMX 
C     Will focus on adding potential matrix to vertical consideration first
            GETTRANS = GETTRANSY(GAUGE, TVALS, NVALS, LIMY, V, E, FLUX,
     +                           WRAPX)
      END IF
      IF (NVALS .EQ. 0) THEN
            WRITE (*,*) 'Invalid current identifier (X and Y only)' 
            STOP
      END IF
          
      RETURN
      END
