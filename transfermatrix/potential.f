c$$$  Tests for adding the potential.
C$$$  Whole matrix at one value

      SUBROUTINE TONE(V,POT,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J, K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT
c      CHARACTER*3 T_1

      CALL DLASET ('ALL',LIMX,LIMY,POT,POT,V,LIMX)
 
c$$$      DO J=1, LIMX
c$$$      WRITE (*, 200) LIMX,LIMY, ' T_1', (V(J, K), K=1, LIMY)  
c$$$   
c$$$      END DO

 200  FORMAT (2I3,A, 100F6.2)

      RETURN
      END

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


c$$$  Differant for sublattice A and B
      SUBROUTINE TTWO(V,POTA,POTB,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POTA,POTB
      
      CALL DLASET ('ALL',LIMX,LIMY,POTB,POTB,V,LIMX)

      DO J=1,LIMX
        DO K=1,LIMY
      IF((MOD(J,2) .EQ.0) .AND. (MOD(K,2) .EQ.0)) V(J,K)=POTA
      IF((MOD(J,2).EQ.1) .AND. (MOD(K,2) .EQ.1)) V(J,K) = POTA
      END DO
      END DO

c$$$      DO J=1, LIMX
c$$$      WRITE (*, 300) LIMX,LIMY, ' T_2', (REAL(V(J, K)), K=1, LIMY)  
c$$$      END DO
 210  FORMAT (2I3,A)
 300  FORMAT (2I3,A, 100F6.2)

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


c     Sharp Barrier with a channel.
      SUBROUTINE TTHREE(V,POT,XDIM,YDIM,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION V(LIMX,LIMY)
      DOUBLE PRECISION POT, ZERO/0.0/
      DOUBLE PRECISION XDIM, YDIM

c      DO J=1, LIMX
c      WRITE(*,310) (V(J, K), K=1, LIMY)  
c      ENDDO
c      PRINT *, "---"
      CALL DLASET ('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)

c     Barrier 
      DO J=1,LIMX
            DO K= 1,LIMY    
             IF((K.GE.(LIMY/2-XDIM/2)).AND.(K .LT.(LIMY/2+XDIM/2)))THEN
                 V(J,K)=POT
              ENDIF
      
c     Channel
              IF((K.GE.(LIMY/2-XDIM/2)).AND.(K .LT.(LIMY/2+XDIM/2)) 
     + .AND.(J.GT.(LIMX/2.0-YDIM/2)).AND.(J .LE.(LIMX/2+YDIM/2)))THEN 
            V(J,K)= ZERO
              ENDIF
            ENDDO
      ENDDO
     
      DO J=1, LIMX
      WRITE(*,310) (V(J, K), K=1, LIMY)  
      ENDDO
 310  FORMAT (100F15.6)
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
            IF((K.GE.(LIMY/2-XDIM/2)).AND.(K .LT.(LIMY/2+XDIM/2)))THEN
                 V(J,K)= V(J,K)+POT
C     Also want checkerboard inside barrier
              ENDIF
      
c     Channel
              IF((K.GE.(LIMY/2-XDIM/2)).AND.(K .LT.(LIMY/2+XDIM/2)) 
     + .AND.(J.GT.(LIMX/2-YDIM/2)).AND.(J .LE.(LIMX/2+YDIM/2)))THEN 
            V(J,K)= ZERO
              ENDIF
            ENDDO
      ENDDO


      DO J=1, LIMX
      WRITE(*,310) (V(J, K), K=1, LIMY)  
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
              IF( J.GE.(LIMX/2-YDIM/2).AND.(J.LT.(LIMX/2+YDIM/2 ))) THEN 
              V(J,K) = V(J,K)+POT
              ENDIF
       
c     Channel
              IF( J.GE.(LIMX/2-YDIM/2).AND.(J.LT.(LIMX/2+YDIM/2 )) 
     +       .AND.(K.GT.(LIMY/2-XDIM/2)).AND.(K.LE.(LIMY/2+XDIM/2)))THEN
             V(J,K)= ZERO
             ENDIF
           ENDDO
      ENDDO


      DO J=1, LIMX
      WRITE(*,310) (V(J, K), K=1, LIMY)  
      ENDDO

 310  FORMAT (100F15.6)
      RETURN
      END














