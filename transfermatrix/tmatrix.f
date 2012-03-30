      DOUBLE PRECISION FUNCTION GETTRANS(CURRENT, GAUGE, TVALS, LIMX,
     +                                   LIMY, E, FLUX, WRAPX)
      IMPLICIT NONE
      INTEGER LIMY, WRAPX, LIMX
      DOUBLE PRECISION TVALS(LIMX), E, FLUX
      CHARACTER CURRENT, GAUGE

c      DOUBLE PRECISION GETTRANSX
c      EXTERNAL GETTRANSX
      DOUBLE PRECISION GETTRANSY
      EXTERNAL GETTRANSY

      IF (CURRENT .EQ. 'X') THEN
         WRITE (*,*) 'ERROR, X current unsuported!'
         STOP
c         GETTRANS = GETTRANSX(GAUGE, TVALS, LIMX, LIMY, E, FLUX, WRAPX)
      ELSE
         IF (CURRENT .EQ. 'Y') THEN
            GETTRANS = GETTRANSY(GAUGE, TVALS, LIMX, LIMY, E, FLUX,
     +                            WRAPX)
         ELSE
            WRITE (*,*) 'Invalid current identifier (X and Y only)' 
            STOP
         END IF
      END IF

      RETURN
      END