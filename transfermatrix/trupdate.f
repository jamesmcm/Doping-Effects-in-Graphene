c$$   NEW SUBROUTINE TO GENERATE A,B,C,D
      SUBROUTINE GENABCD(MULT, U, A, B, C, D, N)
      IMPLICIT NONE
C    M(A,B,C,D) COMPONENTS OF MATRIX, U: PHASE MATRIX
      INTEGER N
      DOUBLE COMPLEX  MULT(2*N, 2*N)
      DOUBLE COMPLEX  U(N, N)
      DOUBLE COMPLEX  A(N, N), B(N, N),
     +                C(N, N), D(N, N)
      DOUBLE COMPLEX  MA(N, N)
c      , MB(N, N),
c     +                MC(N, N), MD(N, N)
      DOUBLE COMPLEX  BU (N, N),  UC (N, N),
     +                DU (N, N),  UDU(N,N),
     +                TEMPBC  (N, N)
      DOUBLE COMPLEX UNITY/1.0/, ZERO/0.0/,   NEG/-1.0/,
     +               HALF/0.5/,  NHALF/-0.5/
      DOUBLE COMPLEX UIC, UJ, MBIJ, MCIJ, MDIJ
      INTEGER I, J

c      MA = MULT(1:N, 1:N)
c      MB = MULT(1:N, (N+1):2*N)
c      MC = MULT((N+1):2*N, 1:N)
c      MD = MULT((N+1):2*N, (N+1):2*N)

c$$$ M(A,B,C,D) INITALISE TO ZERO'S
      CALL SQZERO(A, N)
      CALL SQZERO(B, N)
      CALL SQZERO(C, N)
      CALL SQZERO(D, N)
      CALL SQZERO(TEMPBC, N)

c     BU := Mb*U
c      CALL SQDOT(BU, MB, U, N)
c     DU := Md*U
c      CALL SQDOT(DU, MD, U, N)
c     UC := Udagger * Mc
c      CALL ZGEMM('C', 'N', N, N, N,
c     +     UNITY, U,  N,
c     +            MC, N,
c     +      ZERO, UC, N)
c     UDU := Udagger*DU = Udagger * Md * U     
c      CALL ZGEMM('C', 'N', N, N, N,
c     +     UNITY, U, N,
c     +           DU, N,
c     +     ZERO, UDU, N)

      DO I = 1, N
        UIC = DCONJG(U(I, I))
        DO J = 1, N
          MA(I, J) = MULT (I, J)
          MBIJ = MULT (I, N + J)
          MCIJ = MULT (I + N, J)
          MDIJ = MULT (I + N, J + N)
          
          UJ = U(J, J)
          BU(I, J)  =        MBIJ * UJ
          UC(I, J)  = UIC  * MCIJ
          UDU(I, J) = UIC  * MDIJ * UJ
        ENDDO
      ENDDO

c$$$ CREATED TEMP MATICIES -> GENERATE A
c     Micro-optimization: 16 calls to ZAXPY can be replaced by
c     10 calls to ZAXPY + 3 call to ZCOPY -- AVS
c             
c     A := 1/2 Ma
      CALL SQUPDAXPY (A, HALF, MA, N)
c     B := A = 1/2 Ma     
      CALL SQCOPY    (A, B, N)
c     A := 1/2 (Ma + U*Md*U)
      CALL SQUPDAXPY (A, HALF, UDU, N)
c     D := A = 1/2(Ma + U*Md*U)
      CALL SQCOPY    (A, D, N)
c     TEMPBC := 1/2 (Mb*u + u*Mc)      
      CALL SQUPDAXPY (TEMPBC, HALF, BU, N)      
      CALL SQUPDAXPY (TEMPBC, HALF, UC, N)
c     A := 1/2 (Ma + U*Md*U + Mb*u + u * Mc)
      CALL SQUPDAXPY (A, UNITY, TEMPBC, N)
c     D := 1/2(Ma + U*Md*U - Mb*u - u*Mc) 
      CALL SQUPDAXPY (D, NEG,   TEMPBC, N)
c     C:= B := 1/2 (Ma - U*Md*U)
      CALL SQUPDAXPY (B, NHALF, UDU, N)
      CALL SQCOPY    (B, C, N)
c     TEMPBC := TEMPBC - BU = 1/2(-Mb*U + U*Mc)      
      CALL SQUPDAXPY (TEMPBC, NEG, BU, N)
c     B  := 1/2(Ma - U*Md*U - Mb*U + U*Mc)
      CALL SQUPDAXPY (B, UNITY, TEMPBC, N)
c     C  := 1/2(Ma - U*Md*U + Mb*U - U*Mc)
      CALL SQUPDAXPY (C, NEG,   TEMPBC, N)
      
      RETURN
      END
 
C$$$ ROUTINE TO GENERATE T AND R
      SUBROUTINE GENTANDRINC(A, B, C, D,
     +                       TINC, RINC, TTILDEINC, RTILDEINC,
     +                       N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX UNITY/1.0/
      DOUBLE COMPLEX A(N, N), B(N, N),
     +               C(N, N), D(N, N),
     + TINC(N, N), TTILDEINC(N, N),
     + RINC(N, N), RTILDEINC(N, N)

C$$$ T~ = D^-1
      CALL SQCOPY (D, TTILDEINC, N)
      CALL INVERTMATRIX(TTILDEINC, N)
C$$$ R~ = BD^-1
      CALL SQDOT (RTILDEINC, B, TTILDEINC, N)
C$$$ R = -D^-1 C
      CALL SQDOTAX (RINC, -UNITY, TTILDEINC, C, N)
C$$$ T=(A-)? BD^-1 C
C$$$ Thus, We can simply reuse R here: T = A + B * R
      CALL SQCOPY(A, TINC, N)
      CALL SQDOTUPD(UNITY, TINC, UNITY, B, RINC, N)
      
      RETURN
      END

      SUBROUTINE UPDATETANDR(T,    R,    TTILDE, RTILDE, 
     +                       TINC, RINC, TTILDEINC, RTILDEINC, 
     +                       N )
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX T(N, N), TTILDE(N, N),
     +               R(N, N), RTILDE(N, N),
     +               TINC(N, N), TTILDEINC(N, N),
     +               RINC(N, N), RTILDEINC(N, N)

C$$$ TEMPORARY LOCAL VARIABLES
      DOUBLE COMPLEX T1TEMP(N, N), TTILDE1TEMP(N, N),
     +               R1TEMP(N, N), RTILDE1TEMP(N, N)
      DOUBLE COMPLEX BRACKET(N, N),
     +               TRTEMP(N, N)
      DOUBLE COMPLEX UNITY/1.0/
      
      CALL SQCOPY (T,      T1TEMP,      N)
      CALL SQCOPY (TTILDE, TTILDE1TEMP, N)
      CALL SQCOPY (R,      R1TEMP,      N)
      CALL SQCOPY (RTILDE, RTILDE1TEMP, N)
      
C$$$ BRACKET := (1 - RTILDE1.R2)^-1 . T1
      CALL SQUNIT   (BRACKET, N)
      CALL SQDOTUPD (UNITY, BRACKET, -UNITY, RTILDE1TEMP, RINC, N)
      CALL INVSOLVE (BRACKET, T1TEMP, N)
      
C$$$ T = T2.BRACKET
      CALL SQDOT  (T, TINC, BRACKET, N)
      
C$$$ R = R1 + TTILDE1.(1-R2.R1~)^-1.R2.T1
C$$$   = R1 + T1~.R_2.(1-R1~.R2)^-1.T1 = R1 + T1~.R2.BRACKET
      CALL SQDOT  (TRTEMP, RINC, BRACKET, N)
      CALL SQCOPY (R1TEMP, R, N)
      CALL SQDOTUPD (UNITY, R, UNITY, TTILDE1TEMP, TRTEMP, N)

C$$$ BRACKET := (1 - R2.RTILDE1)^-1.T2~
      CALL SQUNIT (BRACKET, N)
      CALL SQDOTUPD (UNITY, BRACKET, -UNITY, RINC, RTILDE1TEMP, N)
      CALL INVSOLVE (BRACKET, TTILDEINC, N)
      
C$$$  TTILDE = TTILDE1.BRACKET
      CALL SQDOT (TTILDE, TTILDE1TEMP, BRACKET, N)
      
C$$$ RTILDE = RTILDE2 + T2.inv(1 - R1~.R2).RTILDE1.TTILDE2
C$$$        = R2~ + T2.R1~.inv(1 - R2.R1~).T2~
C$$$        = R2~ + T2.R1~.Bracket      
      CALL SQDOT  (TRTEMP,   RTILDE1TEMP, BRACKET, N)
      CALL SQCOPY (RTILDEINC, RTILDE, N)
      CALL SQDOTUPD (UNITY, RTILDE, UNITY, TINC, TRTEMP, N)
      
      RETURN
      END
