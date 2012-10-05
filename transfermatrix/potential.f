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

c$$$c$$$  Generates a random distribution of vacancies
c$$$
c$$$      SUBROUTINE GENVAC(V, LIMX, LIMY, U, ALAT, BLAT)
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER LIMX, LIMY, I, J
c$$$      DOUBLE PRECISION V(LIMX,LIMY), U, ALAT, BLAT
c$$$      REAL RANDOM
c$$$      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0
c$$$
c$$$      CALL DLASET('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)
c$$$
c$$$      DO I = 1, LIMX
c$$$         DO J = 1, LIMY
c$$$            CALL RANDOM_NUMBER(RANDOM)
c$$$            IF (MOD(I, 2) .EQ. MOD(J, 2)) THEN
c$$$               IF (RANDOM .LT. ALAT) THEN
c$$$                  V(I,J) = U
c$$$               END IF
c$$$            ELSE
c$$$               IF (RANDOM .LT. BLAT) THEN
c$$$                  V(I,J) = U
c$$$               END IF
c$$$            END IF
c$$$         END DO
c$$$      END DO
c$$$
c$$$      RETURN
c$$$      END


c$$$c$$$  Calcualtes the running mean and error for the conductance
c$$$      SUBROUTINE INCREMENTCOND(MEANCOND, ERRCOND, NEWCOND, RNUM, NE)
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER NE, I
c$$$      DOUBLE PRECISION MEANCOND(NE), ERRCOND(NE), NEWCOND(NE),
c$$$     +                 SIGMA, RNUM
c$$$c      INTEGER NUM
c$$$
c$$$      DO I = 1, NE
c$$$         SIGMA = 0.0
c$$$c         RNUM = NUM
c$$$         SIGMA = ERRCOND(I) * SQRT(RNUM - 1.0)
c$$$         SIGMA = SIGMA * SIGMA
c$$$         SIGMA = SIGMA + (MEANCOND(I) * MEANCOND(I))
c$$$         SIGMA = SIGMA * (RNUM - 1.0)
c$$$         SIGMA = SIGMA + (NEWCOND(I) * NEWCOND(I))
c$$$         SIGMA = SIGMA / RNUM
c$$$
c$$$         MEANCOND(I) = MEANCOND(I) * (RNUM - 1.0)
c$$$         MEANCOND(I) = (MEANCOND(I) + NEWCOND(I)) / RNUM
c$$$
c$$$         SIGMA = SIGMA - (MEANCOND(I) * MEANCOND(I))
c$$$         SIGMA = SQRT(SIGMA)
c$$$         ERRCOND(I) = SIGMA / SQRT(RNUM)
c$$$      END DO
c$$$
c$$$      RETURN
c$$$      END

