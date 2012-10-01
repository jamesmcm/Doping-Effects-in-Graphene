c$$$  Tests for adding the potential.
C$$$  Whole matrix at one value

C$$$  Adds a constant value across whole of matrix      
      SUBROUTINE TONEADD(V,POT,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J, K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT
c      CHARACTER*3 T_1

      DO J=1, LIMX
         DO K =1,LIMY
            V(J,K)=V(J,K)+POT
         END DO
      END DO

c$$$      DO J=1, LIMX
c$$$      WRITE (*, 200) LIMX,LIMY, ' T_1', (V(J, K), K=1, LIMY)  
c$$$   
c$$$      END DO

      RETURN
      END


c$$$  Adds a constant value to A,B sublattices.      
      SUBROUTINE TTWOADD(V,POTA,POTB,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POTA,POTB
      
      DO J=1,LIMX
        DO K=1,LIMY
           IF((MOD(J,2) .EQ.0) .AND. (MOD(K,2) .EQ.0)) THEN 
              V(J,K)=V(J,K)+POTA
           ELSE IF((MOD(J,2).EQ.1) .AND. (MOD(K,2) .EQ.1)) THEN
              V(J,K) = V(J,K)+POTA
           ELSE
              V(J,K) = V(J,K) + POTB
           ENDIF
        END DO
      END DO

c$$$      DO J=1, LIMX
c$$$      WRITE (*, 300) LIMX,LIMY, ' T_2', (REAL(V(J, K)), K=1, LIMY)  
c$$$      END DO
 210  FORMAT (2I3,A)
 300  FORMAT (2I3,A, 100F6.2)

      RETURN
      END

c     Adds a vertical Sharp Barrier (In X) with a channel to the current matrix
      SUBROUTINE TVERTADD(V,POT,XDIM,YDIM,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT, ZERO/0.0/
      DOUBLE PRECISION XDIM, YDIM
    
c     Barrier 
C     N.B. YDIM corresponds to blocking across X (so blocking y-traversal)
      DO J=1,LIMX
           DO K= 1,LIMY    
            IF((K.GE.(LIMY/2-YDIM/2)).AND.(K .LT.(LIMY/2+YDIM/2)))THEN
                 V(J,K)= V(J,K)+POT
C     Also want checkerboard inside barrier
              ENDIF
      
c     Channel
              IF((K.GE.(LIMY/2-YDIM/2)).AND.(K .LT.(LIMY/2+YDIM/2)) 
     + .AND.(J.GT.(LIMX/2-XDIM/2)).AND.(J .LE.(LIMX/2+XDIM/2)))THEN 
            V(J,K)= ZERO
              ENDIF
            ENDDO
      ENDDO


      DO J=1, LIMX
c      WRITE(*,310) (V(J, K), K=1, LIMY)  
      ENDDO

 310  FORMAT (100F15.6)
      RETURN
      END

c$$$  Gaussian Barrier
      SUBROUTINE TFOUR(V,POT,GWID,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION DUM /0.0/
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT, ZERO/0.0/, GWID
      
      CALL DLASET ('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)
      
      DO J=1,LIMX
       DO K= 1, LIMY
            DUM=K*0.01
            V(J,K)=POT*
     +           DEXP(-1* (((DUM)-(LIMY/2)*0.01)**2/(2*GWID*GWID)))
C            WRITE(*,330) J,K, DUM, REAL(V(J,K))
       ENDDO
      ENDDO
     
      DO J=1, LIMX
c$$$c      WRITE (*, 320) LIMX,LIMY, ' T_4', (REAL(V(J, K)), K=1, LIMY)  
c      WRITE(*,340) (V(J, K), K=1, LIMY)
      ENDDO
            
 320  FORMAT (2I3,A, 100F10.6)
 330  FORMAT (2I3, 2F10.6)
 340  FORMAT (100F13.6)


      RETURN 
      END

c     Adds a horozontal (In Y) Sharp barrier with a channel
      SUBROUTINE THORADD(V,POT,XDIM,YDIM,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT, ZERO/0.0/
      DOUBLE PRECISION XDIM, YDIM
    
c     Barrier 
C     N.B. YDIM corresponds to blocking across X (so blocking y-traversal)
      
      DO J=1,LIMX
            DO K= 1,LIMY    
c     Horozontal barrier     
              IF( J.GE.(LIMX/2-XDIM/2).AND.(J.LT.(LIMX/2+XDIM/2 ))) THEN 
              V(J,K) = V(J,K)+POT
              ENDIF
       
c     Channel
              IF( J.GE.(LIMX/2-XDIM/2).AND.(J.LT.(LIMX/2+XDIM/2 )) 
     +       .AND.(K.GT.(LIMY/2-YDIM/2)).AND.(K.LE.(LIMY/2+YDIM/2)))THEN
             V(J,K)= ZERO
             ENDIF
           ENDDO
      ENDDO


      DO J=1, LIMX
c      WRITE(*,310) (V(J, K), K=1, LIMY)  
      ENDDO

 310  FORMAT (100F15.6)
      RETURN
      END

c$$$  Generates a random distribution of vacancies

      SUBROUTINE GENVAC(V, LIMX, LIMY, U, ALAT, BLAT)
      IMPLICIT NONE

      INTEGER LIMX, LIMY, I, J
      DOUBLE PRECISION V(LIMX,LIMY), U, ALAT, BLAT
      REAL RANDOM
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0

      CALL DLASET('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)

      DO I = 1, LIMX
         DO J = 1, LIMY
            CALL RANDOM_NUMBER(RANDOM)
            IF (MOD(I, 2) .EQ. MOD(J, 2)) THEN
               IF (RANDOM .LT. ALAT) THEN
                  V(I,J) = U
               END IF
            ELSE
               IF (RANDOM .LT. BLAT) THEN
                  V(I,J) = U
               END IF
            END IF
         END DO
      END DO

      RETURN
      END


c$$$  Calcualtes the running mean and error for the conductance
      SUBROUTINE INCREMENTCOND(MEANCOND, ERRCOND, NEWCOND, RNUM, NE)
      IMPLICIT NONE

      INTEGER NE, I
      DOUBLE PRECISION MEANCOND(NE), ERRCOND(NE), NEWCOND(NE),
     +                 SIGMA, RNUM
c      INTEGER NUM

      DO I = 1, NE
         SIGMA = 0.0
c         RNUM = NUM
         SIGMA = ERRCOND(I) * SQRT(RNUM - 1.0)
         SIGMA = SIGMA * SIGMA
         SIGMA = SIGMA + (MEANCOND(I) * MEANCOND(I))
         SIGMA = SIGMA * (RNUM - 1.0)
         SIGMA = SIGMA + (NEWCOND(I) * NEWCOND(I))
         SIGMA = SIGMA / RNUM

         MEANCOND(I) = MEANCOND(I) * (RNUM - 1.0)
         MEANCOND(I) = (MEANCOND(I) + NEWCOND(I)) / RNUM

         SIGMA = SIGMA - (MEANCOND(I) * MEANCOND(I))
         SIGMA = SQRT(SIGMA)
         ERRCOND(I) = SIGMA / SQRT(RNUM)
      END DO

      RETURN
      END

