c$$$  Tests for adding the potential.
C$$$  Whole matrix at one value

      SUBROUTINE TONE(V,POT,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J, K
      DOUBLE COMPLEX V(LIMX,LIMY)
      DOUBLE COMPLEX POT
c      CHARACTER*3 T_1

      CALL ZLASET ('ALL',LIMX,LIMY,POT,POT,V,LIMX)
 
      DO J=1, LIMX
      WRITE (*, 200) LIMX,LIMY, ' T_1', (REAL(V(J, K)), K=1, LIMY)  
   
      END DO

 200  FORMAT (2I3,A, 100F6.2)

      RETURN
      END


c$$$  Differant for sublattice A and B
      SUBROUTINE TTWO(V,POTA,POTB,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE COMPLEX V(LIMX,LIMY)
      DOUBLE COMPLEX POTA,POTB
      
      CALL ZLASET ('ALL',LIMX,LIMY,POTB,POTB,V,LIMX)

      DO J=1,LIMX
        DO K=1,LIMY
      IF((MOD(J,2) .EQ.0) .AND. (MOD(K,2) .EQ.0)) V(J,K)=POTA
      IF((MOD(J,2).EQ.1) .AND. (MOD(K,2) .EQ.1)) V(J,K) = POTA
      END DO
      END DO

      DO J=1, LIMX
      WRITE (*, 300) LIMX,LIMY, ' T_2', (REAL(V(J, K)), K=1, LIMY)  
      END DO
 210  FORMAT (2I3,A)
 300  FORMAT (2I3,A, 100F6.2)

      RETURN
      END

c     Sharp Barrier with a channel.
      SUBROUTINE TTHREE(V,POT,WID,HEIGHT,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K,WID,HEIGHT
      DOUBLE COMPLEX V(LIMX,LIMY)
      DOUBLE COMPLEX POT, ZERO
      
      CALL ZLASET ('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)
      
      DO J=1,LIMX
            DO K= 1,LIMY    
              IF((K.GE.(LIMY/2-WID/2)).AND.(K .LE.(LIMY/2+WID/2)))
     +  V(J,K)=POT 
              IF((K.GE.(LIMX/2-WID/2)).AND.(K .LE.(LIMX/2+WID/2))
     + .AND.(J.GE.(LIMY/2-WID/2)).AND.(J .LE.(LIMY/2+WID/2)))
     +  V(J,K)=ZERO 
            ENDDO
       
      ENDDO
 
      DO J=1, LIMX
      WRITE(*,310) (REAL(V(J, K)), K=1, LIMY)  
      ENDDO
 310  FORMAT (100F10.6)
      RETURN
      END

c$$$  Gaussian Barrier
      SUBROUTINE TFOUR(V,POT,GWID,LIMX,LIMY)
      IMPLICIT NONE 
      INTEGER LIMX,LIMY,J,K
      DOUBLE PRECISION DUM /0.0/
      DOUBLE COMPLEX V(LIMX,LIMY)
      DOUBLE COMPLEX POT, ZERO /0.0/, GWID
      
      CALL ZLASET ('ALL', LIMX, LIMY, ZERO, ZERO, V, LIMX)
      
      DO J=1,LIMX
       DO K= 1, LIMY
            DUM=K*0.01
            V(J,K)=ZEXP(-1* (((DUM)-(LIMY/2)*0.01)**2/(2*GWID*GWID)))
C            WRITE(*,330) J,K, DUM, REAL(V(J,K))
       ENDDO
      ENDDO
     
      DO J=1, LIMX
c      WRITE (*, 320) LIMX,LIMY, ' T_4', (REAL(V(J, K)), K=1, LIMY)  
      WRITE(*,340) (REAL(V(J, K)), K=1, LIMY)  
      ENDDO
      
 320  FORMAT (2I3,A, 100F10.6)
 330  FORMAT (2I3, 2F10.6)
 340  FORMAT (100F10.6)


      RETURN 
      END


















