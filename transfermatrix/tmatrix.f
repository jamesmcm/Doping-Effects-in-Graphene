      DOUBLE PRECISION FUNCTION GETTRANS(CURRENT, GAUGE, TVALS, LIMX,
     +                                   LIMY, E, FLUX, WRAPX)
      IMPLICIT NONE
      INTEGER LIMY, WRAPX, LIMX
      DOUBLE PRECISION TVALS(LIMX), E, FLUX
      CHARACTER CURRENT, GAUGE

      DOUBLE PRECISION GETTRANSX
      EXTERNAL GETTRANSX
      DOUBLE PRECISION GETTRANSY
      EXTERNAL GETTRANSY

      IF (CURRENT .EQ. 'X') THEN
c$$$         WRITE (*,*) 'ERROR, X current unsuported!'
         IF (MOD(LIMY, 2) .NE. 0) THEN
            WRITE (*,*) 'ERROR, LIMY must be even for X current!'
            STOP
         ELSE
         GETTRANS = GETTRANSX(GAUGE, TVALS, LIMX, LIMY/2, E, FLUX)
      END IF
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
