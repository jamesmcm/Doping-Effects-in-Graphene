      DOUBLE PRECISION FUNCTION CHECKUNI(T, R, TTILDE, RTILDE, N)
      IMPLICIT NONE

      INTEGER N, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
      DOUBLE COMPLEX T(N,N), BETA/-1/,ALPHA/1/,
     + R(N,N),TTILDE(N,N),RTILDE(N,N),
     + U(N*2,N*2), CHECK(N*2,N*2)


C$$$ TEST CASES OF T,R,T~ R~
C$$$ FILLS A UNIT MATRIX
      CALL SQUNIT (CHECK, 2*N)
      DO X=1, N
         DO Y=1, N
C$$$  TOP LEFT
            U(X,Y)=T(X,Y)
C$$$  BOTTOM LEFT
            U(X+N,Y)=R(X,Y)
C$$$  TOP RIGHT
            U(X,Y+N)=RTILDE(X,Y)
C$$$  BOTTOM RIGHT
            U(X+N,Y+N)=TTILDE(X,Y)
         END DO
      END DO

C$$$  ZGEMM HAS INBUILT FUNCTION TO FIND U**H
      CALL ZGEMM('N', 'C', 2*N, 2*N, 2*N, ALPHA, U,
     +     2*N, U, 2*N, BETA, CHECK, 2*N)
C$$$  ZLANGE FINDS MATRIX NORM
      CHECKUNI = ZLANGE('F', 2*N, 2*N, CHECK, 2*N)
      RETURN
      END

C$$$  ROUTINE TO CHECK FOR UNITARY MATRICES

      DOUBLE PRECISION FUNCTION CHECKUNI2(T, R, TTILDE, RTILDE, N)
      IMPLICIT NONE
 
      INTEGER N
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
      DOUBLE PRECISION DNORM1, DNORM2, DNORM3, DNORM
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX T(N,N), R(N,N),
     +               TTILDE(N,N),RTILDE(N,N),
     +               CK(N, N)


C$$$ TEST CASES OF T,R,T~ R~

C$$$ FILLS A UNIT MATRIX

C      DO X = 1, 2*N
C         DO Y = 1, 2*N
C            CHECK(X, Y) = (0.0, 0.0)
C         END DO
C       CHECK(X, X) = 1.0
C      END DO

C$$$ ZGEMM HAS INBUILT FUNCTION TO FIND U**H
C$$   CALL PRINTT (U, 2*N, 'U1 ')

c     Unitarity can be also checked using subblocks only
c     The definition of unitarity can be easily cast into the form
c     |T|^2 + |R|^2 = 1, T*R~ + R* T~ = 0

c     Check that T * T + R* R = 1   (* is the Hermitian conjugate)
c     Here I use zherk, which does T^+T
c     Beware: it computes only the upper-diagonal part!
      CALL SQUNIT (CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, T, N,
     +  -ONED, CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, R, N,
     +  ONED, CK, N)
C$$$ ZLANGE FINDS MATRIX NORM
      DNORM1 = ZLANGE ('F', N, N, CK, N)

c     Check that T~* T~ + R~* R~ = 1
      CALL SQUNIT (CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, TTILDE, N,
     +  -ONED, CK, N)
      CALL ZHERK ('U', 'C', N, N, ONED, RTILDE, N,
     +  ONED, CK, N)
      DNORM2 = ZLANGE ('F', N, N, CK, N)

c     Now check that T* R~ + R* T~ = 0
      CALL ZGEMM ('C', 'N', N, N, N, ONEC, T,
     +           N, RTILDE, N, ZEROC, CK, N)
      CALL ZGEMM ('C', 'N', N, N, N, ONEC, R,
     +           N, TTILDE, N, ONEC, CK, N)
      DNORM3 = ZLANGE ('F', N, N, CK, N)

C$$   CALL PRINTT (CHECK, 2 * N, 'C2 ')
c     ZHERK does only the upper-triangular part. Therefore, the norm
c     DNORM1 should be roughly doubled, ditto DNORM2
c     There are two related off-diagonal blocks in the cross-product,
c     this doubles the norm DNORM3. (The relation is not exact,
c     since, the diagonal is not doubled)
      DNORM = 2.0* DNORM1 + 2.0 * DNORM2 + 2.0 * DNORM3
c     WRITE (*, *) 'CHECKUNI2:', DNORM1, DNORM2, DNORM3, DNORM
      CHECKUNI2 = DNORM
      RETURN
      END
C$$$  ROUTINE TO CHECK FOR UNITARY MATRICES

      DOUBLE PRECISION FUNCTION CHECKUNI3(T, R, TTILDE, RTILDE, N)
      IMPLICIT NONE

      INTEGER N, X/1/, Y/1/
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX T(N,      N), R(     N, N),
     +               TTILDE(N, N), RTILDE(N, N),
     +               U(N*2, N*2), CHECK(N*2, N*2)


C$$$ TEST CASES OF T,R,T~ R~

C$$$ FILLS A UNIT MATRIX
      CALL SQUNIT (CHECK, 2*N)
      DO X=1, N
         DO Y=1, N
C$$$  TOP LEFT
            U(X,Y)=T(X,Y)
C$$$  BOTTOM LEFT
            U(X+N,Y)=R(X,Y)
C$$$  TOP RIGHT
            U(X,Y+N)=RTILDE(X,Y)
C$$$  BOTTOM RIGHT
            U(X+N,Y+N)=TTILDE(X,Y)
         END DO
      END DO

C$$$  ZHERK calculates U^H * U, and then subtracts 1.
C$$$  However, this is perfomed in the upper triangular part only
C$$   Therefore, we double the norm.
      CALL ZHERK('U', 'C', 2*N, 2*N, ONED, U, 2*N,
     +     -ONED, CHECK, 2*N)
C$$$  ZLANGE FINDS MATRIX NORM
      CHECKUNI3 = 2 * ZLANGE('F', 2*N, 2*N, CHECK, 2*N)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CORRUNI(T, R, TTILDE, RTILDE,
     +                                  THRESHOLD, N)
      IMPLICIT NONE 
      INTEGER N
      DOUBLE COMPLEX T(N,N), R(N,N),
     +               TTILDE(N,N),RTILDE(N,N)
      DOUBLE PRECISION THRESHOLD 
      
      DOUBLE PRECISION DNORM1, DNORM2, DNORM3, DNORM
      DOUBLE PRECISION ONED/1.0/
      DOUBLE COMPLEX ZEROC/0.0/, ONEC/1.0/
      DOUBLE COMPLEX EPS1(N, N),  EPS2(N, N),
     +               EPS3 (N, N), EPS2H(N, N)
      INTEGER II, JJ 
      DOUBLE PRECISION ZLANGE
      EXTERNAL SQUNIT
 
c
c     By definition, unitary matrix obeys U^+ . U = 1
c     Numerically, the matrix can be different from a unitary, so that
c     the r.h.s.becomes 1 + epsilon. If this is the case, one can       
c     find the nearest unitary matrix. This amounts to changing U as
c
c          U := U - 1/2 * U . epsilon
c
c     Then, indeed, the first-order term in epsilon cancels epsilon
c     on the r.h.s. (It is also possible to implement a more symmetric
c     version of the algorithm: U: = U - 1/4 U . epsilon - 1/4 eta . U^+,
c     where the second term takes care of the defect in U.U^+ = 1 + eta.
c     However, I do not see a real benefit in doing so, while twice as many
c     flops are required to evaluate the two terms.)
      
c     Below, we correct the unitary matrix
c      
c              |  T    R~ |
c         U =  |          |
c              |  R    T~ |
c
c     according to this algorithm. We represent the r.h.s. as
c
c        T^+  . T   + R^+  . R  = 1 + eps_1
c        T~^+ . T~  + R~^+ . R~ = 1 + eps_3    
c        T^+  . R~  + R^+  . T~ = eps_2
c           
c     Then, eps_1 and eps_3 are Hermitian, and also      
c          
c        T~^+ . R  + R~^+ . T = eps_2^+      
c
c     Then, the update is performed as
c
c       T  := T  - 1/2 ( T . eps_1 + R~ . eps_2^+ )
c       R  := R  - 1/2 ( R . eps_1 + T~ . eps_2^+ )       
c       T~ := T~ - 1/2 ( R . eps_2 + T~ . eps_3   )       
c       R~ := R~ - 1/2 ( T . eps_2 + R~ . eps_3   )       
c

c     Calculate eps_1 =  T * T + R* R - 1   (* is the Hermitian conjugate)
c     Here I use zherk, which does T^+T
c     Beware: it computes only the upper-diagonal part!
      CALL SQUNIT (EPS1, N) 
c      CALL ZGEMM ('C', 'N', N, N, N,
c     &  ONEC, T, N, T, N,
c     &  -ONEC, EPS1, N)
c      CALL ZGEMM ('C', 'N', N, N, N,
c     &  ONEC, R, N, R, N,
c     &  ONEC, EPS1, N)
      CALL ZHERK ('U', 'C', N, N, ONED, T, N,
     +  -ONED, EPS1, N)
      CALL ZHERK ('U', 'C', N, N, ONED, R, N,
     +  ONED, EPS1, N)
      
      
      DNORM1 = ZLANGE ('F', N, N, EPS1, N)

c     Find eps_3 =  T~* T~ + R~* R~ - 1
      CALL SQUNIT (EPS3, N)
c      CALL ZGEMM ('C', 'N', N, N, N,
c     &  ONEC, TTILDE, N, TTILDE, N,
c     &  -ONEC, EPS3, N)
c      CALL ZGEMM ('C', 'N', N, N, N,
c     &  ONEC, RTILDE, N, RTILDE, N,
c     &  ONEC, EPS3, N)
      
      CALL ZHERK ('U', 'C', N, N, ONED, TTILDE, N,
     +  -ONED, EPS3, N)
      CALL ZHERK ('U', 'C', N, N, ONED, RTILDE, N,
     +  ONED, EPS3, N)
      DNORM3 = ZLANGE ('F', N, N, EPS3, N)

c     Now calculate eps_2 =  T* R~ + R* T~ 
      CALL SQZERO (EPS2, N)
      CALL ZGEMM ('C', 'N', N, N, N, ONEC, T,
     +           N, RTILDE, N, ZEROC, EPS2, N)
      CALL ZGEMM ('C', 'N', N, N, N, ONEC, R,
     +           N, TTILDE, N, ONEC, EPS2, N)
      DNORM2 = ZLANGE ('F', N, N, EPS2, N)

c     ZHERK did only the upper triangular part of EPS1 and EPS3
c     Here, we fill in the rest.      
      DO II = 1, N
        DO JJ = II + 1, N
          EPS1(JJ, II) = DCONJG (EPS1(II, JJ))
          EPS3(JJ, II) = DCONJG (EPS3(II, JJ))
        ENDDO
      ENDDO
       
      DO II=1, N
        DO JJ=1, N
          EPS2H(II, JJ) = DCONJG(EPS2(JJ, II))
         ENDDO
      ENDDO
      
      DNORM = DNORM1 + DNORM2 + DNORM3

      if (DNORM .GT. THRESHOLD) THEN
        write (*, *) '*** Unitarity correction: ', dnorm
c
c       Note that we do not use temporary variables to preserve
c       old values of T, R, etc. Doing so would only change
c       would eliminate terms of the order of  eps^2, which are
c       ignored by this procedure anyway.
c        
c       T  := T  - 1/2 ( T . eps_1 + R~ . eps_2^+ )
        CALL SQDOTUPD (ONEC, T, -0.5*ONEC, T,      EPS1,  N)
        CALL SQDOTUPD (ONEC, T, -0.5*ONEC, RTILDE, EPS2H, N)

c       R  := R  - 1/2 ( R . eps_1 + T~ . eps_2^+ )       
        CALL SQDOTUPD (ONEC, R, -0.5*ONEC, R,      EPS1,  N)
        CALL SQDOTUPD (ONEC, R, -0.5*ONEC, TTILDE, EPS2H, N)

c       R~ := R~ - 1/2 ( T . eps_2 + R~ . eps_3   )       
        CALL SQDOTUPD (ONEC, RTILDE, -0.5*ONEC, T,      EPS2, N)
        CALL SQDOTUPD (ONEC, RTILDE, -0.5*ONEC, RTILDE, EPS3, N)
        
c       T~ := T~ - 1/2 ( R . eps_2 + T~ . eps_3   )       
        CALL SQDOTUPD (ONEC, TTILDE, -0.5*ONEC, R,      EPS2, N)
        CALL SQDOTUPD (ONEC, TTILDE, -0.5*ONEC, TTILDE, EPS3, N)
      ENDIF

c     WRITE (*, *) 'CORRUNI:', DNORM1, DNORM2, DNORM3, DNORM
      CORRUNI = DNORM
      RETURN
      END

